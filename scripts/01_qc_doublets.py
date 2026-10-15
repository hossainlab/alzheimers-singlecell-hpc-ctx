"""Step 1 - per-sample QC, doublet removal, donor sex / 17q21.31 genotype inference.

Input : data/GSE163577_RAW/*.tar.gz (GEO filtered 10x matrices, 25 samples)
Output: data/processed/qc_counts.h5ad      (all QC-passing singlets, raw integer counts, NO down-sampling)
        tables/qc_cell_metrics.csv.gz       (every nucleus with its QC status)
        tables/qc_sample_summary.csv        (per-sample counts and thresholds)
        tables/sample_metadata.csv          (sample/donor/region/condition/batch/sex/genotype)
"""
import sys, tarfile, tempfile, shutil
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp

sys.path.insert(0, str(Path(__file__).parent))
from config import *  # noqa

sc.settings.verbosity = 0


def mad_outlier(x, nmads, tail="both"):
    med = np.median(x)
    mad = np.median(np.abs(x - med)) * 1.4826
    lo, hi = med - nmads * mad, med + nmads * mad
    if tail == "both":
        return (x < lo) | (x > hi), lo, hi
    return x > hi, -np.inf, hi


def read_sample(row) -> sc.AnnData:
    tmp = Path(tempfile.mkdtemp())
    with tarfile.open(row.file) as tf:
        tf.extractall(tmp, filter="data")
    mtx = next(tmp.rglob("matrix.mtx*")).parent
    ad = sc.read_10x_mtx(mtx, var_names="gene_symbols", make_unique=True)
    shutil.rmtree(tmp, ignore_errors=True)
    ad.obs_names = [f"{row.gsm}_{b}" for b in ad.obs_names]
    return ad


def main():
    st = sample_table()
    kept, cellrows, summ = [], [], []
    for row in st.itertuples():
        ad = read_sample(row)
        X = ad.X.tocsr()
        n_counts = np.asarray(X.sum(1)).ravel()
        n_genes = np.asarray((X > 0).sum(1)).ravel()
        mt = ad.var_names.str.startswith("MT-")
        pct_mt = np.asarray(X[:, mt].sum(1)).ravel() / np.maximum(n_counts, 1) * 100
        o_c, lo_c, hi_c = mad_outlier(np.log1p(n_counts), QC["nmads_counts"])
        o_g, lo_g, hi_g = mad_outlier(np.log1p(n_genes), QC["nmads_genes"])
        o_m, _, hi_m = mad_outlier(pct_mt, QC["nmads_mt"], tail="upper")
        thr_mt = float(np.clip(hi_m, QC["min_mt_ceiling"], QC["max_mt_pct"]))
        low_q = o_c | o_g | (n_genes < QC["min_genes"]) | (pct_mt > thr_mt)

        # doublets on QC-passing nuclei only (scrublet, per sample, fixed seed)
        sub = ad[~low_q].copy()
        sc.pp.scrublet(sub, random_state=SEED, verbose=False)
        dbl = pd.Series(False, index=ad.obs_names)
        pred = sub.obs["predicted_doublet"].astype(bool).values
        if pred.mean() < 0.01:
            # automatic threshold failed (no bimodal simulated-score distribution): fall back to the
            # fixed scrublet score cut-off of 0.25 so that every sample receives doublet filtering
            pred = sub.obs["doublet_score"].values > 0.25
            print(f"  {row.sample}: scrublet auto-threshold failed, fixed cut-off 0.25 -> {pred.sum()} doublets")
        dbl.loc[sub.obs_names] = pred
        dscore = pd.Series(np.nan, index=ad.obs_names)
        dscore.loc[sub.obs_names] = sub.obs["doublet_score"].values

        status = np.where(low_q, "Low quality", np.where(dbl.values, "Doublet", "Pass"))
        cellrows.append(pd.DataFrame(dict(cell=ad.obs_names, gsm=row.gsm, sample=row.sample,
                                          region=row.region, condition=row.condition,
                                          n_counts=n_counts, n_genes=n_genes, pct_mt=pct_mt,
                                          doublet_score=dscore.values, qc_status=status)))
        keep = status == "Pass"
        a = ad[keep].copy()
        a.X = sp.csr_matrix(a.X, dtype=np.float32)
        for k in ["gsm", "sample", "donor", "region", "condition", "batch"]:
            a.obs[k] = getattr(row, k)
        kept.append(a)
        summ.append(dict(gsm=row.gsm, sample=row.sample, donor=row.donor, region=row.region,
                         condition=row.condition, nuclei_input=ad.n_obs, low_quality=int(low_q.sum()),
                         doublets=int((status == "Doublet").sum()), nuclei_retained=int(keep.sum()),
                         retention_pct=round(100 * keep.mean(), 1),
                         umi_min=round(float(np.expm1(lo_c)), 1), umi_max=round(float(np.expm1(hi_c)), 1),
                         genes_min=max(QC["min_genes"], round(float(np.expm1(lo_g)), 1)),
                         genes_max=round(float(np.expm1(hi_g)), 1), mt_max_pct=round(thr_mt, 2),
                         median_umi_retained=float(np.median(n_counts[keep])),
                         median_genes_retained=float(np.median(n_genes[keep]))))
        print(f"{row.sample:22s} in={ad.n_obs:6d} lowQ={low_q.sum():5d} dbl={(status=='Doublet').sum():5d} kept={keep.sum():6d}", flush=True)

    summ = pd.DataFrame(summ)
    pd.concat(cellrows).to_csv(TABLES / "qc_cell_metrics.csv.gz", index=False)

    adata = sc.concat(kept, join="outer", fill_value=0)
    adata.X = sp.csr_matrix(adata.X)
    ncell = np.asarray((adata.X > 0).sum(0)).ravel()
    adata = adata[:, ncell >= QC["min_cells_per_gene"]].copy()

    # ---- donor sex and 17q21.31 genotype from pseudobulk CPM of each sample
    lib = np.asarray(adata.X.sum(1)).ravel()
    def sample_cpm(genes):
        genes = [g for g in genes if g in adata.var_names]
        M = adata[:, genes].X.toarray()
        df = pd.DataFrame(M, columns=genes); df["sample"] = adata.obs["sample"].values
        tot = pd.Series(lib).groupby(adata.obs["sample"].values).sum()
        return df.groupby("sample").sum().div(tot, axis=0) * 1e6
    ygenes = ["UTY", "USP9Y", "DDX3Y", "KDM5D", "NLGN4Y", "TTTY14"]
    cp = sample_cpm(["XIST"] + ygenes + ["ARL17B", "LRRC37A2", "KANSL1-AS1"])
    cp["chrY_sum"] = cp[[g for g in ygenes if g in cp]].sum(1)
    cp["sex"] = np.where(cp["XIST"] > cp["chrY_sum"], "F", "M")
    cp["ARL17B_cpm"] = cp["ARL17B"]
    cp["H2_17q21_carrier"] = cp["ARL17B"] > 100   # bimodal: ~10-60 vs 150-400 CPM
    meta = st.drop(columns="file").merge(
        cp[["XIST", "chrY_sum", "sex", "ARL17B_cpm", "H2_17q21_carrier"]].reset_index(), on="sample")
    # concordance of the two tissues from the same donor
    paired = meta.groupby("donor").filter(lambda d: len(d) == 2)
    conc = paired.groupby("donor").agg(sex_concordant=("sex", lambda s: s.nunique() == 1),
                                       genotype_concordant=("H2_17q21_carrier", lambda s: s.nunique() == 1))
    conc.to_csv(TABLES / "donor_pairing_concordance.csv")
    print("paired donors:", len(conc), "sex concordant:", int(conc.sex_concordant.sum()),
          "17q21 genotype concordant:", int(conc.genotype_concordant.sum()))
    meta.to_csv(TABLES / "sample_metadata.csv", index=False)
    summ = summ.merge(meta[["sample", "sex", "H2_17q21_carrier"]], on="sample")
    summ.to_csv(TABLES / "qc_sample_summary.csv", index=False)

    for k in ["sex", "H2_17q21_carrier"]:
        adata.obs[k] = adata.obs["sample"].map(meta.set_index("sample")[k]).values
    adata.obs["H2_17q21_carrier"] = adata.obs["H2_17q21_carrier"].astype(str)
    adata.write_h5ad(PROC / "qc_counts.h5ad", compression="gzip")
    print("final:", adata.shape)


if __name__ == "__main__":
    main()
