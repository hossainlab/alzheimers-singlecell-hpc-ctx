"""Step 6 - Figure 3: cross-cell-type convergence and cross-region concordance (donor-level DE).

a  significant genes per cell type, split into robust vs artefact-flagged (ambient / 17q21.31 / sex)
b  overlap of robust DEGs across the 5 NVU cell types (UpSet-style), hippocampus
c  pairwise correlation of AD effects (log2FC) between cell types, hippocampus
d  genes robust in >= 2 hippocampal cell types: log2FC heat map across cell types and both regions
e  hippocampus vs cortex log2FC for hippocampal robust DEGs (same donors: concordance, not replication)
f  immediate-early genes claimed in the legacy README as 'universally convergent'

Outputs: figures/Fig3_Cross_Cell_Convergence/Fig3_convergence.{pdf,png}
         tables/convergent_genes.csv, tables/cross_region_concordance.csv, tables/celltype_lfc_correlation.csv
"""
import sys, itertools
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import TwoSlopeNorm
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from config import *  # noqa

OUT = FIGS / "Fig3_Cross_Cell_Convergence"
IEG = ["FOS", "FOSB", "JUN", "JUNB", "EGR1", "DUSP1", "GADD45B", "ATF3", "IER2", "NR4A1"]
SHORT = {"Microglia": "Micro", "Endothelial": "Endo", "Pericyte": "Peri", "Astrocyte": "Astro", "Oligodendrocyte": "Oligo"}


def main():
    set_style()
    de = pd.read_csv(TABLES / "pseudobulk_DE_all.csv.gz")
    H = de[de.region == "Hippocampus"]; C = de[de.region == "Cortex"]
    fig = plt.figure(figsize=(7.2, 7.6))
    gs = gridspec.GridSpec(3, 7, figure=fig, height_ratios=[1, 1.15, 1], hspace=0.75, wspace=1.1)

    # a: significant genes split by flag
    ax = fig.add_subplot(gs[0, 0:2])
    rows = []
    for r_ in ["Hippocampus", "Cortex"]:
        for ct in NVU5:
            d = de[(de.region == r_) & (de.cell_type == ct)]
            sig = (d.padj < 0.05) & (d.log2FoldChange.abs() >= 0.5)
            rows.append(dict(region=r_, ct=ct, robust=int((sig & d.robust).sum()),
                             ambient=int((sig & d.flag_contamination).sum()),
                             genotype=int((sig & d.flag_genotype_sex & ~d.flag_contamination).sum())))
    S = pd.DataFrame(rows)
    y = np.arange(len(S)); left = np.zeros(len(S))
    for col, colr, lab in [("robust", "#2A6F97", "Robust"), ("ambient", "#BDBDBD", "Ambient-flagged"),
                           ("genotype", "#7A7A7A", "17q21.31 / sex")]:
        ax.barh(y, S[col], left=left, color=colr, height=0.75, label=lab); left += S[col].values
    ax.set_yticks(y); ax.set_yticklabels([f"{'Hpc' if r=='Hippocampus' else 'Ctx'} {SHORT[c]}" for r, c in zip(S.region, S.ct)])
    ax.invert_yaxis(); ax.set_xlabel("DEGs (FDR < 0.05, |log$_2$FC| ≥ 0.5)")
    ax.legend(loc="lower right", fontsize=5, handlelength=0.8, bbox_to_anchor=(1.0, 0.0))
    ax.set_title("DEGs after artefact filtering")
    ax.axhline(4.5, color="0.6", lw=0.5)
    panel_letter(ax, "a", x=-0.5)

    # b: UpSet (hippocampus robust)
    sets = {ct: set(H[(H.cell_type == ct) & H.robust].gene) for ct in NVU5}
    allg = set().union(*sets.values())
    memb = pd.DataFrame({ct: [g in sets[ct] for g in sorted(allg)] for ct in NVU5}, index=sorted(allg))
    combos = memb.value_counts().reset_index(name="n").sort_values("n", ascending=False)
    combos = combos[combos[NVU5].sum(1) >= 1].head(14)
    gsb = gs[0, 3:7].subgridspec(2, 1, height_ratios=[1.3, 1], hspace=0.05)
    axb = fig.add_subplot(gsb[0]); axm = fig.add_subplot(gsb[1], sharex=axb)
    xs = np.arange(len(combos))
    k = combos[NVU5].sum(1).values
    axb.bar(xs, combos.n, color=["#2A6F97" if kk > 1 else "0.6" for kk in k], width=0.7)
    for x_, n_ in zip(xs, combos.n):
        axb.text(x_, n_, str(n_), ha="center", va="bottom", fontsize=5)
    axb.set_ylabel("Robust DEGs"); axb.tick_params(labelbottom=False)
    axb.set_title("Most hippocampal robust DEGs are cell-type specific")
    for i, ct in enumerate(NVU5):
        on = combos[ct].values
        axm.scatter(xs[~on], np.full((~on).sum(), i), s=10, color="0.85")
        axm.scatter(xs[on], np.full(on.sum(), i), s=10, color=CT_COLORS[ct])
    for x_, row in zip(xs, combos[NVU5].values):
        ii = np.where(row)[0]
        if len(ii) > 1:
            axm.plot([x_, x_], [ii.min(), ii.max()], color="0.3", lw=0.8, zorder=0)
    axm.set_yticks(range(5)); axm.set_yticklabels([f"{c} ({len(sets[c])})" for c in NVU5], fontsize=5.5)
    axm.set_ylim(4.6, -0.6); axm.set_xticks([]); axm.spines["bottom"].set_visible(False)
    panel_letter(axb, "b", x=-0.14)

    # c: LFC correlation between cell types (Hpc)
    ax = fig.add_subplot(gs[1, 0:2])
    W = H.pivot_table(index="gene", columns="cell_type", values="log2FoldChange")[NVU5].dropna()
    flagged = set(H[H.flag_contamination | H.flag_genotype_sex].gene)
    W = W[~W.index.isin(flagged)]
    Rm = W.corr(method="spearman")
    Rm.to_csv(TABLES / "celltype_lfc_correlation_hippocampus.csv")
    im = ax.imshow(Rm.values, cmap="RdBu_r", norm=TwoSlopeNorm(0, -1, 1))
    for i, j in itertools.product(range(5), range(5)):
        ax.text(j, i, f"{Rm.values[i, j]:.2f}", ha="center", va="center", fontsize=5,
                color="w" if abs(Rm.values[i, j]) > 0.6 else "k")
    ax.set_xticks(range(5)); ax.set_xticklabels([SHORT[c] for c in NVU5], rotation=45, ha="right")
    ax.set_yticks(range(5)); ax.set_yticklabels([SHORT[c] for c in NVU5])
    cax = ax.inset_axes([0.0, -0.42, 0.6, 0.04]); cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cb.set_label("Spearman ρ of log$_2$FC", fontsize=5, labelpad=1); cb.ax.tick_params(labelsize=4.5)
    ax.set_title(f"Shared AD effects, Hpc ({len(W):,} genes)")
    panel_letter(ax, "c", x=-0.45)

    # d: convergent genes heatmap
    ax = fig.add_subplot(gs[1, 3:7])
    cnt = memb.sum(1)
    conv = cnt[cnt >= 2].index.tolist()
    rows = []
    for g in conv:
        lf = H[(H.gene == g)].set_index("cell_type").reindex(NVU5)
        sgn = np.sign(lf.loc[[c for c in NVU5 if g in sets[c]], "log2FoldChange"])
        rows.append(dict(gene=g, n_celltypes_robust_hpc=int(cnt[g]), consistent_direction=bool(sgn.nunique() == 1),
                         celltypes=";".join([c for c in NVU5 if g in sets[c]]),
                         median_log2FC_hpc=float(lf.loc[[c for c in NVU5 if g in sets[c]], "log2FoldChange"].median())))
    CV = pd.DataFrame(rows)
    if len(CV):
        CV = CV.sort_values(["n_celltypes_robust_hpc", "median_log2FC_hpc"], ascending=[False, True])
        # cortex concordance per gene (direction in the same cell types)
        cc = []
        for r in CV.itertuples():
            cts_ = r.celltypes.split(";")
            lc = C[(C.gene == r.gene) & C.cell_type.isin(cts_)].log2FoldChange
            cc.append(np.nan if lc.empty else float(np.mean(np.sign(lc) == np.sign(r.median_log2FC_hpc))))
        CV["ctx_same_direction_frac"] = cc
    CV.to_csv(TABLES / "convergent_genes.csv", index=False)
    show = CV.head(30).gene.tolist() if len(CV) else []
    if show:
        cols = [(r_, ct) for r_ in ["Hippocampus", "Cortex"] for ct in NVU5]
        M = np.full((len(show), len(cols)), np.nan); Pm = np.full_like(M, np.nan)
        for j, (r_, ct) in enumerate(cols):
            d = de[(de.region == r_) & (de.cell_type == ct)].set_index("gene")
            M[:, j] = d.reindex(show).log2FoldChange.values; Pm[:, j] = d.reindex(show).padj.values
        im = ax.imshow(M.T, aspect="auto", cmap="RdBu_r", norm=TwoSlopeNorm(0, -2, 2), interpolation="nearest")
        for i, j in zip(*np.where(Pm.T < 0.05)):
            ax.text(j, i, "•", ha="center", va="center", fontsize=5, color="k")
        ax.set_yticks(range(len(cols))); ax.set_yticklabels([f"{'Hpc' if r=='Hippocampus' else 'Ctx'} {SHORT[c]}" for r, c in cols])
        ax.set_xticks(range(len(show))); ax.set_xticklabels(show, rotation=90, style="italic", fontsize=5)
        ax.axhline(4.5, color="k", lw=0.8)
        cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.01); cb.set_label("log$_2$FC (AD/Ctrl)", fontsize=5); cb.ax.tick_params(labelsize=5)
        ax.set_title(f"Robust in ≥ 2 hippocampal NVU types (n = {len(CV)}; • FDR < 0.05; white = n.t.)")
    else:
        ax.text(0.5, 0.5, "No gene robust in ≥ 2 cell types", ha="center", transform=ax.transAxes); ax.axis("off")
    panel_letter(ax, "d", x=-0.12)

    # e: cross-region concordance
    ax = fig.add_subplot(gs[2, 0:3])
    rows = []
    for ct in NVU5:
        h = H[(H.cell_type == ct) & H.robust].set_index("gene")
        c = C[C.cell_type == ct].set_index("gene").reindex(h.index).dropna(subset=["log2FoldChange"])
        h = h.loc[c.index]
        same = (np.sign(h.log2FoldChange) == np.sign(c.log2FoldChange))
        p = stats.binomtest(int(same.sum()), len(same), 0.5).pvalue if len(same) else np.nan
        rho = stats.spearmanr(h.log2FoldChange, c.log2FoldChange).statistic if len(same) > 3 else np.nan
        rows.append(dict(cell_type=ct, hpc_robust_tested_in_ctx=len(same), same_direction=int(same.sum()),
                         pct_same_direction=round(100 * same.mean(), 1) if len(same) else np.nan, binom_p=p,
                         spearman_rho=rho, also_ctx_padj05=int((c.padj < 0.05).sum())))
        ax.scatter(h.log2FoldChange, c.log2FoldChange, s=4, color=CT_COLORS[ct], lw=0, alpha=0.8,
                   label=f"{SHORT[ct]}: {rows[-1]['pct_same_direction']}% same sign (n={len(same)})")
    ax.axhline(0, color="0.7", lw=0.5); ax.axvline(0, color="0.7", lw=0.5)
    lim = 4; ax.plot([-lim, lim], [-lim, lim], color="0.6", ls=":", lw=0.6)
    ax.set_xlim(-lim, lim); ax.set_ylim(-lim, lim)
    ax.set_xlabel("Hippocampus log$_2$FC (robust DEGs)"); ax.set_ylabel("Cortex log$_2$FC, same donors")
    ax.legend(fontsize=4.8, loc="upper left", handletextpad=0.1, markerscale=1.5)
    ax.set_title("Same-donor cortex agrees in direction")
    pd.DataFrame(rows).to_csv(TABLES / "cross_region_concordance.csv", index=False)
    panel_letter(ax, "e", x=-0.2)

    # f: IEGs
    ax = fig.add_subplot(gs[2, 3:7])
    cols = [(r_, ct) for r_ in ["Hippocampus", "Cortex"] for ct in NVU5]
    M = np.full((len(IEG), len(cols)), np.nan); Pm = np.full_like(M, np.nan)
    for j, (r_, ct) in enumerate(cols):
        d = de[(de.region == r_) & (de.cell_type == ct)].set_index("gene")
        M[:, j] = d.reindex(IEG).log2FoldChange.values; Pm[:, j] = d.reindex(IEG).padj.values
    im = ax.imshow(M, aspect="auto", cmap="RdBu_r", norm=TwoSlopeNorm(0, -2, 2), interpolation="nearest")
    for i, j in zip(*np.where(np.isnan(M))):
        ax.text(j, i, "n.t.", ha="center", va="center", fontsize=4, color="0.5")
    for i, j in zip(*np.where(Pm < 0.05)):
        ax.text(j, i, "•", ha="center", va="center", fontsize=6)
    ax.set_yticks(range(len(IEG))); ax.set_yticklabels(IEG, style="italic")
    ax.set_xticks(range(len(cols))); ax.set_xticklabels([f"{'Hpc' if r=='Hippocampus' else 'Ctx'} {SHORT[c]}" for r, c in cols], rotation=90)
    ax.axvline(4.5, color="k", lw=0.8)
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02); cb.set_label("log$_2$FC", fontsize=5); cb.ax.tick_params(labelsize=5)
    ax.set_title("Immediate-early genes: no convergent AD effect (n.t. = not tested)")
    panel_letter(ax, "f", x=-0.2)

    save_fig(fig, OUT / "Fig3_convergence")
    print(pd.DataFrame(rows))


if __name__ == "__main__":
    main()
