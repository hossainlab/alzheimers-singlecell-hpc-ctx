"""Step 4 - donor-level pseudobulk differential expression (AD vs Control) for the 5 NVU cell types.

Unit of replication = sample (one donor x one region). Counts are summed per sample x cell type
and tested with DESeq2 (pydeseq2), design  ~ sex + vascular_yield + condition.

* vascular_yield = logit of the sample's fraction of vascular nuclei (endothelial + pericyte +
  smooth muscle + fibroblast). VINE-seq enrichment efficiency differs systematically between AD
  and control samples, which changes both composition and the ambient-RNA background.
* Hippocampus (8 Control v 9 AD, independent donors) is the discovery set. Cortex (4 v 4) comes
  from 8 of the same donors and is reported as a second tissue, not as independent replication.

Per gene x cell type, two artefact flags are attached:
  flag_contamination : gene is >=5x more abundant in another cell type of the same region AND
                       its level in the target pseudobulk correlates (Spearman >= 0.5) with the
                       sample's fraction of that source cell type (ambient / doublet carry-over)
  flag_genotype_sex  : gene in the 17q21.31 copy-number region or a sex-chromosome gene
`robust` = padj < 0.05 and |log2FC| >= 0.5 and no flag.

Outputs: tables/pseudobulk_DE_<Region>_<CellType>.csv, tables/pseudobulk_DE_all.csv.gz,
         tables/pseudobulk_DE_summary.csv, tables/pseudobulk_DE_sensitivity_no_yield.csv,
         data/processed/pseudobulk/*.csv (count matrices)
"""
import sys, warnings, logging
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from config import *  # noqa

from pydeseq2.dds import DeseqDataSet
from pydeseq2.ds import DeseqStats
from pydeseq2.default_inference import DefaultInference

warnings.filterwarnings("ignore")
logging.disable(logging.WARNING)
INF = DefaultInference(n_cpus=8)
MIN_NUCLEI = 30
PB_DIR = PROC / "pseudobulk"; PB_DIR.mkdir(exist_ok=True)


def pseudobulk(ad, mask, key="sample"):
    X = ad.layers["counts"][mask]
    s = ad.obs.loc[mask, key].astype(str).values
    u, inv = np.unique(s, return_inverse=True)
    G = sp.csr_matrix((np.ones(len(s)), (inv, np.arange(len(s)))), shape=(len(u), len(s)))
    M = (G @ X).toarray()
    return pd.DataFrame(M, index=u, columns=ad.var_names), pd.Series(np.bincount(inv), index=u)


def run_deseq(C, meta, design):
    dds = DeseqDataSet(counts=C.astype(int), metadata=meta, design=design, refit_cooks=True,
                       inference=INF, quiet=True)
    dds.deseq2()
    st = DeseqStats(dds, contrast=["condition", "AD", "Control"], inference=INF, quiet=True)
    st.summary()
    res = st.results_df.copy()
    try:
        coef = [c for c in dds.varm["LFC"].columns if c.startswith("condition")][0]
        st.lfc_shrink(coeff=coef)
        res["log2FC_shrunk"] = st.results_df["log2FoldChange"].values
    except Exception:
        res["log2FC_shrunk"] = np.nan
    return res


def main():
    ad = sc.read_h5ad(PROC / "integrated.h5ad")
    ad = ad[ad.obs.cell_type != "Mixed"].copy()
    smeta = pd.read_csv(TABLES / "sample_metadata.csv").set_index("sample")
    comp = pd.crosstab(ad.obs["sample"], ad.obs["cell_type"])
    frac = comp.div(comp.sum(1), axis=0)
    assert "vascular_yield" in smeta, "run 02_integrate_annotate.py first"

    sexgenes = ad.var_names[ad.var_names.str.match(SEX_GENES_PATTERN)]
    allres, summ, sens = [], [], []
    for region in ["Hippocampus", "Cortex"]:
        rmask = (ad.obs.region == region).values
        # region-level pseudobulk CPM per cell type (for contamination source detection)
        cpm_ct, cpm_sample_ct = {}, {}
        for ct in comp.columns:
            m = rmask & (ad.obs.cell_type == ct).values
            if m.sum() < 200:
                continue
            C, n = pseudobulk(ad, m)
            C.to_csv(PB_DIR / f"{region}_{ct.replace(' ', '_')}_counts.csv.gz")
            n.to_csv(PB_DIR / f"{region}_{ct.replace(' ', '_')}_ncells.csv")
            cpm = C.div(C.sum(1), axis=0) * 1e6
            cpm_sample_ct[ct] = cpm
            cpm_ct[ct] = cpm.mean(0)
        cpm_ct = pd.DataFrame(cpm_ct)

        for ct in NVU5:
            C, n = pseudobulk(ad, rmask & (ad.obs.cell_type == ct).values)
            keep = n.index[n >= MIN_NUCLEI]
            C, n = C.loc[keep], n.loc[keep]
            meta = smeta.loc[keep, ["condition", "sex", "vascular_yield"]].copy()
            meta["condition"] = pd.Categorical(meta.condition, ["Control", "AD"])
            meta["vascular_yield"] = (meta.vascular_yield - meta.vascular_yield.mean()) / meta.vascular_yield.std()
            ng = meta.condition.value_counts().min()
            gk = (C >= 10).sum(0) >= ng
            gk &= (C.div(C.sum(1), axis=0) * 1e6 >= 5).sum(0) >= ng
            C = C.loc[:, gk]
            res = run_deseq(C, meta, "~sex + vascular_yield + condition")
            res0 = run_deseq(C, meta, "~sex + condition")

            # artefact flags
            tgt = cpm_ct[ct].reindex(res.index).fillna(0)
            others = cpm_ct.drop(columns=ct).reindex(res.index).fillna(0)
            src = others.idxmax(1)
            ratio = others.max(1) / (tgt + 0.5)
            rho = pd.Series(np.nan, index=res.index)
            cand = res.index[ratio >= 5]
            tcpm = cpm_sample_ct[ct].loc[keep]
            for g in cand:
                s = src[g]
                f = frac.loc[keep, s]
                rho[g] = stats.spearmanr(tcpm[g], f).statistic
            res["contamination_source"] = np.where(ratio >= 5, src, "")
            res["source_to_target_ratio"] = ratio.round(2)
            res["rho_with_source_fraction"] = rho.round(3)
            res["flag_contamination"] = (ratio >= 5) & (rho >= 0.5)
            res["flag_genotype_sex"] = res.index.isin(LOCUS_17Q21) | res.index.isin(sexgenes)
            res["robust"] = ((res.padj < 0.05) & (res.log2FoldChange.abs() >= 0.5)
                             & ~res.flag_contamination & ~res.flag_genotype_sex)
            res["gene"] = res.index
            res["region"], res["cell_type"] = region, ct
            res["n_control"], res["n_AD"] = int((meta.condition == "Control").sum()), int((meta.condition == "AD").sum())
            res = res.sort_values("pvalue")
            res.to_csv(TABLES / f"pseudobulk_DE_{region}_{ct}.csv", index=False)
            allres.append(res)
            sig = res.padj < 0.05
            summ.append(dict(region=region, cell_type=ct, n_control=res.n_control.iloc[0], n_AD=res.n_AD.iloc[0],
                             median_nuclei_per_sample=int(n.median()), genes_tested=len(res),
                             padj05=int(sig.sum()), padj05_flagged_contamination=int((sig & res.flag_contamination).sum()),
                             padj05_flagged_genotype_sex=int((sig & res.flag_genotype_sex).sum()),
                             robust_up=int((res.robust & (res.log2FoldChange > 0)).sum()),
                             robust_down=int((res.robust & (res.log2FoldChange < 0)).sum())))
            j = res0.reindex(res.index)
            sens.append(dict(region=region, cell_type=ct, padj05_with_yield=int(sig.sum()),
                             padj05_without_yield=int((j.padj < 0.05).sum()),
                             lfc_correlation=round(np.corrcoef(res.log2FoldChange, j.log2FoldChange)[0, 1], 3)))
            print(summ[-1], flush=True)

    pd.concat(allres).to_csv(TABLES / "pseudobulk_DE_all.csv.gz", index=False)
    pd.DataFrame(summ).to_csv(TABLES / "pseudobulk_DE_summary.csv", index=False)
    pd.DataFrame(sens).to_csv(TABLES / "pseudobulk_DE_sensitivity_no_yield.csv", index=False)


if __name__ == "__main__":
    main()
