"""Step 3 - Figure 1: global single-nucleus landscape of hippocampus AND cortex (all cell types).

Output: figures/Fig1_Global_Landscape/Fig1_global_landscape.{pdf,png}
        tables/Fig1_composition_per_sample.csv, tables/Fig1_composition_tests.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from config import *  # noqa

OUT = FIGS / "Fig1_Global_Landscape"


def embed(ax, xy, colors, order=None, s=0.15, rng=np.random.default_rng(SEED), arrows=True):
    idx = rng.permutation(len(xy)) if order is None else order
    ax.scatter(xy[idx, 0], xy[idx, 1], c=np.asarray(colors)[idx], s=s, lw=0, rasterized=True)
    ax.set_xticks([]); ax.set_yticks([])
    for sp_ in ax.spines.values():
        sp_.set_visible(False)
    if not arrows:
        return
    kw = dict(xycoords="axes fraction", textcoords="axes fraction", annotation_clip=False,
              arrowprops=dict(arrowstyle="-|>", lw=0.6, color="k", shrinkA=0, shrinkB=0))
    ax.annotate("", xy=(0.13, 0.0), xytext=(0.0, 0.0), **kw)
    ax.annotate("", xy=(0.0, 0.13), xytext=(0.0, 0.0), **kw)
    ax.text(0.065, -0.015, "UMAP1", transform=ax.transAxes, fontsize=5, ha="center", va="top")
    ax.text(-0.015, 0.065, "UMAP2", transform=ax.transAxes, fontsize=5, ha="right", va="center", rotation=90)


def main():
    set_style()
    ad = sc.read_h5ad(PROC / "integrated.h5ad")
    ad = ad[ad.obs.cell_type != "Mixed"].copy()
    smeta = pd.read_csv(TABLES / "sample_metadata.csv").set_index("sample")
    cts = [c for c in CELLTYPE_ORDER if c in set(ad.obs.cell_type)]
    xy = ad.obsm["X_umap"]
    n_nuc = ad.n_obs

    fig = plt.figure(figsize=(7.2, 8.6))
    gs = gridspec.GridSpec(3, 6, figure=fig, height_ratios=[1.45, 0.95, 1.05], hspace=0.38, wspace=0.6)

    # a: cell types (labels de-overlapped with leader lines)
    ax = fig.add_subplot(gs[0, 0:4])
    embed(ax, xy, ad.obs.cell_type.map(CT_COLORS).values)
    ax.set_aspect("equal", adjustable="datalim")
    txt = []
    for ct in cts:
        m = (ad.obs.cell_type == ct).values
        cx, cy = np.median(xy[m, 0]), np.median(xy[m, 1])
        txt.append(ax.text(cx, cy, f"{ct} ({m.sum():,})", fontsize=5.5, ha="center", va="center", color="k",
                           bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.8)))
    from adjustText import adjust_text
    adjust_text(txt, ax=ax, arrowprops=dict(arrowstyle="-", lw=0.4, color="0.3"), expand=(1.15, 1.3),
                only_move={"text": "xy"})
    n_donor = smeta.donor.nunique()
    ax.set_title(f"{n_nuc:,} nuclei, {len(smeta)} samples, {n_donor} donors, {len(cts)} cell classes", pad=4)
    panel_letter(ax, "a", x=-0.03, y=1.03)

    # b, c: region and condition, stacked on the right
    sub = gs[0, 4:6].subgridspec(2, 1, hspace=0.35)
    small = np.random.default_rng(1).choice(n_nuc, min(n_nuc, 80000), replace=False)
    for k, (key, cmap, ttl) in enumerate([("region", REGION_COLORS, "Regions co-embed after integration"),
                                         ("condition", COND_COLORS, "Condition")]):
        ax = fig.add_subplot(sub[k])
        embed(ax, xy[small], ad.obs[key].map(cmap).values[small], s=0.25, arrows=False)
        ax.set_aspect("equal", adjustable="datalim")
        hs = [Line2D([], [], marker="o", ls="", color=c, ms=3, label=f"{l}\n({(ad.obs[key]==l).sum():,})")
              for l, c in cmap.items()]
        ax.legend(handles=hs, loc="center left", bbox_to_anchor=(-0.22, 0.5), ncol=1, handletextpad=0.1,
                  fontsize=5, borderaxespad=0, labelspacing=0.8)
        ax.set_title(ttl, pad=3)
        panel_letter(ax, "bc"[k], x=-0.12, y=1.03)

    # d: marker dot plot
    ax = fig.add_subplot(gs[1, :])
    genes = []
    for ct in cts:
        genes += [g for g in DISPLAY_MARKERS[ct] if g in ad.var_names and g not in genes]
    X = ad[:, genes].X
    expr = pd.DataFrame(X.toarray() if hasattr(X, "toarray") else X, columns=genes)
    expr["ct"] = ad.obs.cell_type.values
    mean = expr.groupby("ct").mean().loc[cts]
    frac = (expr.set_index("ct") > 0).groupby(level=0).mean().loc[cts]
    z = (mean - mean.min()) / (mean.max() - mean.min())
    yy, xx = np.meshgrid(np.arange(len(cts)), np.arange(len(genes)), indexing="ij")
    scat = ax.scatter(xx.ravel(), yy.ravel(), s=(frac.values.ravel() * 28) + 0.3, c=z.values.ravel(),
                      cmap="Reds", vmin=0, vmax=1, lw=0.2, edgecolors="0.6")
    ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes, rotation=90, style="italic")
    ax.set_yticks(range(len(cts))); ax.set_yticklabels(cts)
    for i, ct in enumerate(cts):
        ax.get_yticklabels()[i].set_color(CT_COLORS[ct])
    ax.set_xlim(-0.6, len(genes) - 0.4); ax.set_ylim(len(cts) - 0.4, -0.6)
    ax.set_title("Canonical markers confirm every annotated class")
    cax = ax.inset_axes([1.01, 0.55, 0.012, 0.4]); cb = fig.colorbar(scat, cax=cax)
    cb.set_label("Scaled mean\nexpression", fontsize=5); cb.ax.tick_params(labelsize=5)
    hl = [Line2D([], [], marker="o", ls="", mfc="0.5", mec="0.6", ms=np.sqrt(f_ * 28 + 0.3), label=lab)
          for f_, lab in [(0.25, "25%"), (0.5, "50%"), (1.0, "100%")]]
    ax.legend(handles=hl, title="Nuclei\nexpressing", loc="lower left", bbox_to_anchor=(1.005, -0.05), fontsize=5)
    panel_letter(ax, "d", x=-0.11)

    # e: per-sample composition
    comp = pd.crosstab(ad.obs["sample"], ad.obs.cell_type)[cts]
    frac_s = comp.div(comp.sum(1), axis=0) * 100
    order = (smeta.loc[frac_s.index].reset_index()
             .assign(r=lambda d: d.region.map({"Hippocampus": 0, "Cortex": 1}),
                     c=lambda d: d.condition.map({"Control": 0, "AD": 1}))
             .sort_values(["r", "c", "vascular_fraction"], ascending=[True, True, False])["sample"].tolist())
    ax = fig.add_subplot(gs[2, 0:4])
    bottom = np.zeros(len(order))
    for ct in cts:
        v = frac_s.loc[order, ct].values
        ax.bar(range(len(order)), v, bottom=bottom, color=CT_COLORS[ct], width=0.85, lw=0)
        bottom += v
    ax.set_xticks(range(len(order)))
    ax.set_xticklabels([smeta.loc[s, "donor"] for s in order], rotation=90, fontsize=5)
    ax.set_ylabel("Nuclei per sample (%)"); ax.set_ylim(0, 100); ax.set_xlim(-0.6, len(order) - 0.4)
    # group brackets
    grp = [(smeta.loc[s, "region"], smeta.loc[s, "condition"]) for s in order]
    start = 0
    for i in range(1, len(grp) + 1):
        if i == len(grp) or grp[i] != grp[start]:
            r_, c_ = grp[start]
            ax.plot([start - 0.35, i - 0.65], [103, 103], color=COND_COLORS[c_], lw=1.5, clip_on=False)
            ax.text((start + i - 1) / 2, 105, f"{'Hpc' if r_=='Hippocampus' else 'Ctx'} {c_}", ha="center",
                    va="bottom", fontsize=5.5, color=COND_COLORS[c_])
            start = i
    ax.set_xlabel("Donor", labelpad=1)
    ax.set_title("Composition varies with vascular-enrichment yield", pad=12)
    panel_letter(ax, "e", x=-0.09, y=1.1)

    # f: vascular fraction per sample (enrichment yield) by region x condition
    ax = fig.add_subplot(gs[2, 4:6])
    rows = []
    for i, (r_, lab) in enumerate([("Hippocampus", "Hpc"), ("Cortex", "Ctx")]):
        for j, c_ in enumerate(["Control", "AD"]):
            s = smeta[(smeta.region == r_) & (smeta.condition == c_)]
            v = s.vascular_fraction.values * 100
            x = i * 2.6 + j + np.random.default_rng(i * 3 + j).uniform(-0.15, 0.15, len(v))
            ax.scatter(x, v, s=12, color=COND_COLORS[c_], edgecolor="k", lw=0.3, zorder=3)
            ax.hlines(np.median(v), i * 2.6 + j - 0.3, i * 2.6 + j + 0.3, color="k", lw=0.9)
        a = smeta[(smeta.region == r_) & (smeta.condition == "AD")].vascular_fraction
        b = smeta[(smeta.region == r_) & (smeta.condition == "Control")].vascular_fraction
        p = stats.mannwhitneyu(a, b).pvalue
        rows.append(dict(region=r_, AD_median=a.median(), Control_median=b.median(), n_AD=len(a), n_Control=len(b), MWU_p=p))
        ax.text(i * 2.6 + 0.5, 97, f"P = {p:.2g}", ha="center", fontsize=5.5)
    ax.set_xticks([0, 1, 2.6, 3.6]); ax.set_xticklabels(["Ctrl", "AD", "Ctrl", "AD"])
    ax.text(0.5, -0.17, "Hippocampus", transform=ax.get_xaxis_transform(), ha="center", fontsize=6)
    ax.text(3.1, -0.17, "Cortex", transform=ax.get_xaxis_transform(), ha="center", fontsize=6)
    ax.set_ylabel("Vascular nuclei per sample (%)"); ax.set_ylim(0, 100)
    ax.set_title("AD samples yield fewer vascular nuclei", pad=12)
    panel_letter(ax, "f", x=-0.3, y=1.1)
    pd.DataFrame(rows).to_csv(TABLES / "Fig1_vascular_yield_tests.csv", index=False)

    # composition table + donor-level tests (descriptive; enrichment-dependent)
    out = frac_s.join(smeta[["donor", "region", "condition", "sex", "vascular_fraction"]])
    out.to_csv(TABLES / "Fig1_composition_per_sample.csv")
    tests = []
    for r_ in ["Hippocampus", "Cortex"]:
        d = out[out.region == r_]
        for ct in cts:
            a, b = d[d.condition == "AD"][ct], d[d.condition == "Control"][ct]
            tests.append(dict(region=r_, cell_type=ct, AD_median_pct=a.median(), Control_median_pct=b.median(),
                              MWU_p=stats.mannwhitneyu(a, b).pvalue))
    t = pd.DataFrame(tests); t["FDR"] = bh(t.MWU_p); t.to_csv(TABLES / "Fig1_composition_tests.csv", index=False)

    # legend for cell-type colours under panel e
    fig.legend(handles=[plt.Rectangle((0, 0), 1, 1, color=CT_COLORS[c]) for c in cts], labels=cts,
               loc="lower center", ncol=7, bbox_to_anchor=(0.5, 0.0), fontsize=5.5, handlelength=1, columnspacing=0.8)
    fig.subplots_adjust(bottom=0.1)
    save_fig(fig, OUT / "Fig1_global_landscape")
    print("saved", OUT)


if __name__ == "__main__":
    main()
