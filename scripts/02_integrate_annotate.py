"""Step 2 - joint hippocampus + cortex integration, clustering and cell-type annotation.

Input : data/processed/qc_counts.h5ad
Output: data/processed/integrated.h5ad   (X = log-normalised, layers['counts'] = raw counts,
                                          obsm: X_pca, X_pca_harmony, X_umap; obs: leiden, cell_type)
        tables/cluster_annotation.csv     (per-cluster marker-module scores and assigned label)
        tables/celltype_counts_per_sample.csv
Annotation is cluster-level: each Leiden cluster gets the cell type whose canonical-marker module
score (config.MARKERS) is highest. Clusters whose best and second-best scores are within 25 % of
each other AND that co-express two lineages are labelled 'Mixed' (heterotypic doublets /
ambient-dominated) and removed before downstream analysis; this is reported in the table.
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc

sys.path.insert(0, str(Path(__file__).parent))
from config import *  # noqa

sc.settings.verbosity = 1
sc.settings.n_jobs = 8


def main(resolution=1.0):
    ad = sc.read_h5ad(PROC / "qc_counts.h5ad")
    ad.layers["counts"] = ad.X.copy()
    sc.pp.normalize_total(ad, target_sum=1e4)
    sc.pp.log1p(ad)

    # HVGs per sample (seurat_v3 on raw counts); exclude genes whose variance reflects donor
    # sex, 17q21.31 genotype, mitochondrial or ribosomal content rather than cell identity
    excl = (ad.var_names.str.match(SEX_GENES_PATTERN) | ad.var_names.isin(LOCUS_17Q21)
            | ad.var_names.str.startswith(("MT-", "RPL", "RPS")))
    hv = sc.pp.highly_variable_genes(ad[:, ~excl].copy(), n_top_genes=3000, flavor="seurat_v3",
                                     layer="counts", batch_key="sample", inplace=False)
    ad.var["highly_variable"] = False
    ad.var.loc[ad.var_names[~excl][hv["highly_variable"].values], "highly_variable"] = True

    sub = ad[:, ad.var.highly_variable].copy()
    sc.pp.scale(sub, max_value=10)
    sc.pp.pca(sub, n_comps=50, random_state=SEED)
    ad.obsm["X_pca"] = sub.obsm["X_pca"]
    ad.uns["pca"] = sub.uns["pca"]
    del sub
    harmony(ad, key="sample")
    sc.pp.neighbors(ad, n_neighbors=30, n_pcs=30, use_rep="X_pca_harmony", random_state=SEED)
    sc.tl.umap(ad, random_state=SEED, min_dist=0.3)
    sc.tl.leiden(ad, resolution=resolution, key_added="leiden", random_state=SEED,
                 flavor="igraph", n_iterations=2, directed=False)

    annotate(ad)
    ad.write_h5ad(PROC / "integrated.h5ad", compression="gzip")
    summarise(ad)


# Lineage families: labels inside a family are resolved on raw module scores (not 'Mixed')
FAMILY = {"Pericyte": "mural", "Smooth muscle": "mural", "Fibroblast": "mural",
          "Microglia": "myeloid", "Macrophage": "myeloid",
          "Excitatory neuron": "neuron", "Inhibitory neuron": "neuron"}
# One highly specific marker per lineage family, used to detect heterotypic-doublet clusters
KEY_MARKER = {"Endothelial": "CLDN5", "mural": "PDGFRB", "Astrocyte": "AQP4", "myeloid": "C3",
              "Oligodendrocyte": "MBP", "OPC": "PDGFRA", "neuron": "GAD1", "T cell": "CD3E", "Ependymal": "FOXJ1"}


def annotate(ad):
    """Cluster-level annotation.

    1. Module score per marker set (config.MARKERS); cluster label = highest z-scored mean score.
    2. If the best and second-best labels belong to the same lineage family (mural, myeloid,
       neuron) the label is chosen on the raw score within that family.
    3. Otherwise, a best-vs-second margin < 25 % marks the cluster 'Mixed'.
    4. Heterotypic-doublet check: a cluster is 'Mixed' if the key marker of a different family is
       detected in >= 30 % of its nuclei AND at >= 3x that marker's median detection in clusters
       of other families (MBP, the dominant ambient transcript, is thereby tolerated at background).
    'Mixed' nuclei are excluded from all downstream analyses.
    """
    fam = lambda c: FAMILY.get(c, c)
    for ct, genes in MARKERS.items():
        g = [x for x in genes if x in ad.var_names]
        sc.tl.score_genes(ad, g, score_name=f"score_{ct}", random_state=SEED, use_raw=False)
    sc_cols = [f"score_{c}" for c in MARKERS]
    cl = ad.obs.groupby("leiden", observed=True)[sc_cols].mean()
    cl.columns = list(MARKERS)
    z = (cl - cl.mean()) / cl.std()
    rank = z.apply(lambda r: r.sort_values(ascending=False).index.tolist(), axis=1)
    best = rank.str[0]; second = rank.str[1]
    srt = np.sort(z.values, 1)
    margin = pd.Series((srt[:, -1] - srt[:, -2]) / np.maximum(np.abs(srt[:, -1]), 1e-9), index=cl.index)
    label, reason = {}, {}
    for c in cl.index:
        b, s2 = best[c], second[c]
        if fam(b) == fam(s2) and b != s2:
            members = [m for m in MARKERS if fam(m) == fam(b)]
            # neurons are a minor, vascular-depleted population here: excitatory/inhibitory
            # scores are not separable at cluster level, so they are reported jointly
            label[c] = "Neuron" if fam(b) == "neuron" else cl.loc[c, members].idxmax()
            reason[c] = "family-resolved"
        elif margin[c] < 0.25:
            label[c] = "Mixed"; reason[c] = f"low margin ({b}/{s2})"
        else:
            label[c] = b; reason[c] = "marker score"
    # heterotypic doublet check on key-marker detection
    keys = {f: g for f, g in KEY_MARKER.items() if g in ad.var_names}
    det = pd.DataFrame((ad[:, list(keys.values())].X > 0).toarray(), columns=list(keys))
    det = det.groupby(ad.obs.leiden.values, observed=True).mean().reindex(cl.index)
    for c in cl.index:
        if label[c] == "Mixed":
            continue
        own = fam(label[c])
        for f in keys:
            if f == own:
                continue
            others = [k for k in cl.index if label[k] != "Mixed" and fam(label[k]) != f]
            bg = np.median(det.loc[others, f])
            if det.loc[c, f] >= max(0.30, 3 * bg):
                label[c] = "Mixed"; reason[c] = f"heterotypic: {own}+{f} ({keys[f]} {det.loc[c, f]:.2f} vs bg {bg:.2f})"
                break
    ann = cl.round(4).copy()
    ann["n_nuclei"] = ad.obs.leiden.value_counts().reindex(cl.index).values
    ann["best"], ann["second"], ann["margin"] = best, second, margin.round(3)
    ann = ann.join(det.add_prefix("detect_").round(3))
    ann["label"] = pd.Series(label); ann["reason"] = pd.Series(reason)
    ann.to_csv(TABLES / "cluster_annotation.csv")
    ad.obs["cell_type"] = ad.obs.leiden.map(ann["label"]).astype(str)
    print(ann[["n_nuclei", "best", "second", "margin", "label", "reason"]].to_string())


def summarise(ad):
    ct = pd.crosstab(ad.obs["sample"], ad.obs["cell_type"])
    ct.to_csv(TABLES / "celltype_counts_per_sample.csv")
    print(ct.sum())
    # per-sample vascular-enrichment yield (fraction of vascular nuclei among annotated nuclei);
    # used as a covariate downstream because VINE-seq yield differs between AD and control samples
    ct = ct.drop(columns=[c for c in ["Mixed"] if c in ct])
    vf = ct[[c for c in VASCULAR if c in ct]].sum(1) / ct.sum(1)
    smeta = pd.read_csv(TABLES / "sample_metadata.csv").set_index("sample")
    smeta["n_nuclei_annotated"] = ct.sum(1)
    smeta["vascular_fraction"] = vf.clip(1e-3, 1 - 1e-3)
    smeta["vascular_yield"] = np.log(smeta.vascular_fraction / (1 - smeta.vascular_fraction))
    smeta.to_csv(TABLES / "sample_metadata.csv")


def reannotate():
    """Re-run only the annotation on an existing integrated object (python 02_... --annotate-only)."""
    ad = sc.read_h5ad(PROC / "integrated.h5ad")
    annotate(ad)
    ad.write_h5ad(PROC / "integrated.h5ad", compression="gzip")
    summarise(ad)


if __name__ == "__main__":
    reannotate() if "--annotate-only" in sys.argv else main()
