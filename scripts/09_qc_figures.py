"""Step 9 - Supplementary quality-control figure (nuclei QC, doublets, donor sex and genotype).

Output: figures/Quality_Control/FigS1_quality_control.{pdf,png}
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from config import *  # noqa

OUT = FIGS / "Quality_Control"


def main():
    set_style()
    q = pd.read_csv(TABLES / "qc_cell_metrics.csv.gz")
    s = pd.read_csv(TABLES / "qc_sample_summary.csv")
    m = pd.read_csv(TABLES / "sample_metadata.csv").set_index("sample")
    order = (m.reset_index().assign(r=lambda d: d.region.map({"Hippocampus": 0, "Cortex": 1}),
                                    c=lambda d: d.condition.map({"Control": 0, "AD": 1}))
             .sort_values(["r", "c", "donor"])["sample"].tolist())
    lab = [f"{m.loc[x, 'donor']}" for x in order]
    passq = q[q.qc_status == "Pass"]

    fig = plt.figure(figsize=(7.2, 8.0))
    gs = gridspec.GridSpec(5, 4, figure=fig, height_ratios=[1, 1, 1, 1.1, 1.2], hspace=0.65, wspace=0.6)
    for k, (col, ylab, logy) in enumerate([("n_genes", "Genes per nucleus", True), ("n_counts", "UMIs per nucleus", True),
                                           ("pct_mt", "Mitochondrial reads (%)", False)]):
        ax = fig.add_subplot(gs[k, :])
        data = [passq.loc[passq["sample"] == x, col].values for x in order]
        if logy:
            data = [np.log10(np.maximum(d, 1)) for d in data]
        vp = ax.violinplot(data, showextrema=False, widths=0.8)
        for b, x in zip(vp["bodies"], order):
            b.set_facecolor(COND_COLORS[m.loc[x, "condition"]]); b.set_alpha(0.75); b.set_edgecolor("none")
        ax.scatter(range(1, len(order) + 1), [np.median(d) for d in data], s=4, color="k", zorder=3)
        ax.set_xlim(0.4, len(order) + 0.6)
        ax.set_ylabel(ylab + (" (log$_{10}$)" if logy else ""))
        nH = (m.loc[order].region == "Hippocampus").sum()
        ax.axvline(nH + 0.5, color="k", lw=0.6)
        ax.set_xticks(range(1, len(order) + 1)); ax.set_xticklabels(lab if k == 2 else [], rotation=90, fontsize=5)
        if k == 0:
            ax.text(nH / 2, 1.02, "Hippocampus", transform=ax.get_xaxis_transform(), ha="center", fontsize=6)
            ax.text(nH + (len(order) - nH) / 2, 1.02, "Cortex", transform=ax.get_xaxis_transform(), ha="center", fontsize=6)
            ax.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=c, alpha=0.75) for c in COND_COLORS.values()],
                      labels=list(COND_COLORS), loc="upper right", ncol=2)
            ax.set_title("Per-nucleus QC metrics after filtering (all retained nuclei, no down-sampling)", pad=12)
        panel_letter(ax, "abc"[k], x=-0.1)

    # d: retention
    ax = fig.add_subplot(gs[3, 0:3])
    S = s.set_index("sample").loc[order]
    b0 = np.zeros(len(order))
    for col, c, l in [("nuclei_retained", "#2A6F97", "Retained"), ("low_quality", "#E0A458", "Low quality"), ("doublets", "#7A7A7A", "Doublet")]:
        ax.bar(range(len(order)), S[col] / 1000, bottom=b0, color=c, width=0.8, label=l); b0 += S[col].values / 1000
    ax.set_xticks(range(len(order))); ax.set_xticklabels(lab, rotation=90, fontsize=5)
    ax.set_ylabel("Nuclei (thousands)"); ax.legend(ncol=3, loc="upper left"); ax.set_xlim(-0.6, len(order) - 0.4)
    ax.axvline(nH - 0.5, color="k", lw=0.6)
    ax.set_title(f"{int(S.nuclei_retained.sum()):,} of {int(S.nuclei_input.sum()):,} nuclei retained ({100*S.nuclei_retained.sum()/S.nuclei_input.sum():.1f}%)")
    panel_letter(ax, "d", x=-0.09)

    # e: median genes by condition (QC not confounded with condition?)
    ax = fig.add_subplot(gs[3, 3])
    for i, r_ in enumerate(["Hippocampus", "Cortex"]):
        for j, c_ in enumerate(["Control", "AD"]):
            v = S[(S.region == r_) & (S.condition == c_)].median_genes_retained.values
            x = i * 2.6 + j
            ax.scatter(x + np.random.default_rng(i + j).uniform(-0.15, 0.15, len(v)), v, s=8, color=COND_COLORS[c_], edgecolor="k", lw=0.3)
            ax.hlines(np.median(v), x - 0.3, x + 0.3, color="k", lw=0.8)
        a = S[(S.region == r_) & (S.condition == "AD")].median_genes_retained; b = S[(S.region == r_) & (S.condition == "Control")].median_genes_retained
        ax.text(i * 2.6 + 0.5, S.median_genes_retained.max() * 1.05, f"P = {stats.mannwhitneyu(a, b).pvalue:.2g}", ha="center", fontsize=5)
    ax.set_xticks([0.5, 3.1]); ax.set_xticklabels(["Hpc", "Ctx"]); ax.set_ylabel("Median genes / nucleus")
    ax.set_title("Library complexity", pad=10)
    panel_letter(ax, "e", x=-0.6)

    # f: sex
    ax = fig.add_subplot(gs[4, 0:2])
    for r_, mk in [("Hippocampus", "o"), ("Cortex", "^")]:
        d = m[m.region == r_]
        ax.scatter(np.log10(d.XIST + 1), np.log10(d.chrY_sum + 1), s=14, marker=mk, c=[COND_COLORS[c] for c in d.condition], edgecolor="k", lw=0.3, label=r_)
    ax.set_xlabel("log$_{10}$(XIST CPM + 1)"); ax.set_ylabel("log$_{10}$(chrY genes CPM + 1)")
    nf = (m.drop_duplicates("donor").sex == "F").sum()
    ax.set_title(f"Inferred sex: {nf} female, {m.donor.nunique() - nf} male donors")
    ax.legend(loc="center right", fontsize=5)
    panel_letter(ax, "f", x=-0.2)

    # g: 17q21.31 genotype, paired lines for same donor
    ax = fig.add_subplot(gs[4, 2:4])
    xpos = {"Hippocampus": 0, "Cortex": 1}
    for k, (dnr, d) in enumerate(m.groupby("donor")):
        d = d.sort_values("region", ascending=False)
        x = [xpos[r] for r in d.region]
        jitter = np.random.default_rng(k).uniform(-0.12, 0.12)
        xs = np.array(x) + jitter
        if len(d) == 2:
            ax.plot(xs, d.ARL17B_cpm, color="0.7", lw=0.5, zorder=1)
        ax.scatter(xs, d.ARL17B_cpm, s=12, color=[COND_COLORS[c] for c in d.condition], edgecolor="k", lw=0.3, zorder=3)
    ax.axhline(100, color="0.4", ls=":", lw=0.6)
    ax.set_xticks([0, 1]); ax.set_xticklabels(["Hippocampus", "Cortex"]); ax.set_xlim(-0.5, 1.5)
    ax.set_ylabel("ARL17B CPM (sample pseudobulk)")
    conc = pd.read_csv(TABLES / "donor_pairing_concordance.csv")
    ax.set_title(f"17q21.31 marker (ARL17B); lines join one donor's tissues ({int(conc.genotype_concordant.sum())}/{len(conc)} concordant)", fontsize=6)
    panel_letter(ax, "g", x=-0.2)

    save_fig(fig, OUT / "FigS1_quality_control")


if __name__ == "__main__":
    main()
