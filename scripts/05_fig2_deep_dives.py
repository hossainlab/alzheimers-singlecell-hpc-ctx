"""Step 5 - Figure 2 deep dives for the 5 neurovascular-unit (NVU) cell types.

For each of Microglia, Endothelial, Pericyte, Astrocyte, Oligodendrocyte (both regions jointly):
  a  sub-clustering UMAP (Harmony on sample)            b  same embedding by region / condition
  c  sub-cluster marker dot plot                        d  donor-level sub-cluster proportions
  e,f pseudobulk volcano, hippocampus / cortex          g  donor-level heatmap of robust DEGs
  h  pre-ranked GSEA (Hallmark, GO-BP, Reactome) on the DESeq2 Wald statistic
Artefact-flagged genes (ambient carry-over, 17q21.31, sex) are drawn as open grey circles and
excluded from GSEA ranking.

Outputs: figures/Fig2_Cell_Type_Deep_Dives/Fig2<letter>_<CellType>.{pdf,png}
         tables/Fig2_<CellType>_subcluster_markers.csv, _subcluster_proportions.csv, _subcluster_tests.csv
         tables/GSEA_<Region>_<CellType>.csv, tables/GSEA_all.csv
         data/processed/subclusters_<CellType>.h5ad
"""
import sys, warnings
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from matplotlib.colors import TwoSlopeNorm
from scipy import stats
import gseapy as gp

sys.path.insert(0, str(Path(__file__).parent))
from config import *  # noqa
from importlib import import_module
embed = import_module("03_fig1_global_landscape").embed

warnings.filterwarnings("ignore")
OUT = FIGS / "Fig2_Cell_Type_Deep_Dives"
LETTER = {"Microglia": "a", "Endothelial": "b", "Pericyte": "c", "Astrocyte": "d", "Oligodendrocyte": "e"}
SUB_PAL = ["#4477AA", "#EE6677", "#228833", "#CCBB44", "#66CCEE", "#AA3377", "#BBBBBB", "#8E6C8A",
           "#EE8866", "#44BB99", "#99DDFF", "#FFAABB"]
REPLOT = "--replot" in sys.argv
LIBS = {"Hallmark": "MSigDB_Hallmark_2020", "GO-BP": "GO_Biological_Process_2023", "Reactome": "Reactome_2022"}


def load_gmt(name):
    gs = {}
    for line in open(GENESETS / f"{name}.gmt", encoding="utf8"):
        p = line.rstrip("\n").split("\t")
        gs[p[0]] = [g for g in p[2:] if g]
    return gs


def subcluster(ad_ct, res=0.4):
    a = ad_ct.copy()
    excl = (a.var_names.str.match(SEX_GENES_PATTERN) | a.var_names.isin(LOCUS_17Q21)
            | a.var_names.str.startswith(("MT-", "RPL", "RPS")))
    expressed = np.asarray((a.layers["counts"] > 0).sum(0)).ravel() >= 20
    excl = excl | ~expressed
    try:   # per-sample seurat_v3 HVGs; LOESS can fail for samples with few nuclei of this type
        hv = sc.pp.highly_variable_genes(a[:, ~excl].copy(), n_top_genes=2000, flavor="seurat_v3",
                                         layer="counts", batch_key="sample", inplace=False)
    except ValueError:
        hv = sc.pp.highly_variable_genes(a[:, ~excl].copy(), n_top_genes=2000, flavor="seurat_v3",
                                         layer="counts", inplace=False)
    a.var["highly_variable"] = False
    a.var.loc[a.var_names[~excl][hv["highly_variable"].values], "highly_variable"] = True
    sub = a[:, a.var.highly_variable].copy()
    sc.pp.scale(sub, max_value=10)
    sc.pp.pca(sub, n_comps=30, random_state=SEED)
    a.obsm["X_pca"] = sub.obsm["X_pca"]
    harmony(a, key="sample")
    sc.pp.neighbors(a, n_neighbors=20, n_pcs=20, use_rep="X_pca_harmony", random_state=SEED)
    sc.tl.umap(a, random_state=SEED, min_dist=0.3)
    sc.tl.leiden(a, resolution=res, key_added="sub", random_state=SEED, flavor="igraph", n_iterations=2, directed=False)
    # merge tiny clusters (<1 % of nuclei) into nearest by PCA centroid
    vc = a.obs["sub"].value_counts()
    small = vc.index[vc < 0.01 * a.n_obs]
    if len(small) and len(small) < len(vc):
        cen = pd.DataFrame(a.obsm["X_pca_harmony"][:, :20]).groupby(a.obs["sub"].values).mean()
        big = [c for c in vc.index if c not in small]
        remap = {s: big[int(np.argmin(((cen.loc[big] - cen.loc[s]) ** 2).sum(1).values))] for s in small}
        a.obs["sub"] = a.obs["sub"].astype(str).replace(remap)
    order = a.obs["sub"].value_counts().index.tolist()
    a.obs["sub"] = pd.Categorical(a.obs["sub"].astype(str).map({c: f"S{i+1}" for i, c in enumerate(order)}),
                                  [f"S{i+1}" for i in range(len(order))])
    return a


def volcano(ax, de, title, n_label=8):
    d = de.dropna(subset=["padj"]).copy()
    d["nlp"] = -np.log10(d.pvalue.clip(lower=1e-300))
    flagged = d.flag_contamination | d.flag_genotype_sex
    ns = ~flagged & ~d.robust
    ax.scatter(d.log2FoldChange[ns], d.nlp[ns], s=2, color="0.75", lw=0, rasterized=True)
    ax.scatter(d.log2FoldChange[flagged], d.nlp[flagged], s=6, facecolor="none", edgecolor="0.45", lw=0.4,
               rasterized=True)
    up = d.robust & (d.log2FoldChange > 0); dn = d.robust & (d.log2FoldChange < 0)
    ax.scatter(d.log2FoldChange[up], d.nlp[up], s=5, color=COND_COLORS["AD"], lw=0)
    ax.scatter(d.log2FoldChange[dn], d.nlp[dn], s=5, color=COND_COLORS["Control"], lw=0)
    # FDR 0.05 line: smallest -log10 p among padj<0.05
    if (d.padj < 0.05).any():
        ax.axhline(d.loc[d.padj < 0.05, "nlp"].min(), color="0.4", ls=":", lw=0.6)
    lab = pd.concat([d[up].nsmallest(n_label // 2, "pvalue"), d[dn].nsmallest(n_label // 2, "pvalue")])
    try:
        from adjustText import adjust_text
        txt = [ax.text(r.log2FoldChange, r.nlp, r.gene, fontsize=5, style="italic") for r in lab.itertuples()]
        adjust_text(txt, ax=ax, arrowprops=dict(arrowstyle="-", lw=0.3, color="0.3"), expand=(1.2, 1.4))
    except Exception:
        pass
    lim = max(np.nanpercentile(np.abs(d.log2FoldChange), 99.8), d.loc[d.robust, "log2FoldChange"].abs().max() if d.robust.any() else 0) * 1.12
    ax.set_xlim(-lim, lim); ax.axvline(0, color="0.8", lw=0.5, zorder=0)
    ax.set_ylim(-0.3, d.nlp.max() * 1.15)
    ax.set_xlabel("log$_2$ fold change (AD / Control)"); ax.set_ylabel("$-$log$_{10}$ P")
    nU, nD = int(up.sum()), int(dn.sum())
    ax.set_title(title)
    ax.text(0.98, 0.03, f"{nU} up", transform=ax.transAxes, ha="right", va="bottom", color=COND_COLORS["AD"], fontsize=6)
    ax.text(0.02, 0.03, f"{nD} down", transform=ax.transAxes, ha="left", va="bottom", color=COND_COLORS["Control"], fontsize=6)
    ax.margins(y=0.05)


def run_gsea(de, gsets, region, ct):
    d = de[~(de.flag_contamination | de.flag_genotype_sex)].dropna(subset=["stat"])
    rnk = d.set_index("gene")["stat"].sort_values(ascending=False)
    rnk = rnk[~rnk.index.duplicated()]
    out = []
    for lab, gs in gsets.items():
        pre = gp.prerank(rnk=rnk, gene_sets=gs, min_size=15, max_size=500, permutation_num=1000,
                         seed=SEED, threads=4, outdir=None, verbose=False)
        r = pre.res2d.copy()
        r["library"] = lab
        out.append(r)
    r = pd.concat(out)
    r = r.rename(columns={"Term": "term", "NES": "NES", "FDR q-val": "FDR", "NOM p-val": "p"})
    r["NES"] = r.NES.astype(float); r["FDR"] = r.FDR.astype(float); r["p"] = r.p.astype(float)
    r["region"], r["cell_type"] = region, ct
    r = r.sort_values("FDR")
    r.to_csv(TABLES / f"GSEA_{region}_{ct}.csv", index=False)
    return r


def short(term, n=46):
    t = term.split(" (GO:")[0].split(" R-HSA")[0]
    return t if len(t) <= n else t[: n - 1] + "…"


def main():
    set_style()
    ad = None if REPLOT else sc.read_h5ad(PROC / "integrated.h5ad")
    smeta = pd.read_csv(TABLES / "sample_metadata.csv").set_index("sample")
    gsets = {k: load_gmt(v) for k, v in LIBS.items()}
    all_gsea = []
    for ct in NVU5:
        print("==", ct, flush=True)
        cached = REPLOT and (PROC / f"subclusters_{ct}.h5ad").exists()
        if cached:   # re-plot from saved subclustering / marker / GSEA results
            a = sc.read_h5ad(PROC / f"subclusters_{ct}.h5ad")
            mk = pd.read_csv(TABLES / f"Fig2_{ct}_subcluster_markers.csv")
        else:
            a = subcluster(ad[ad.obs.cell_type == ct])
            sc.tl.rank_genes_groups(a, "sub", method="wilcoxon", use_raw=False, pts=True)
            mk = sc.get.rank_genes_groups_df(a, None)
            mk = mk[(mk.logfoldchanges > 0.5) & (mk.pvals_adj < 0.01) & (mk.pct_nz_group > 0.25)]
            mk = mk[~mk.names.str.match(SEX_GENES_PATTERN) & ~mk.names.isin(LOCUS_17Q21) & ~mk.names.str.startswith(("MT-", "RPL", "RPS"))]
            mk.to_csv(TABLES / f"Fig2_{ct}_subcluster_markers.csv", index=False)
        subs = list(a.obs["sub"].cat.categories)
        pal = {s: SUB_PAL[i % len(SUB_PAL)] for i, s in enumerate(subs)}
        top = {s: mk[mk.group == s].sort_values("scores", ascending=False).names.head(3).tolist() for s in subs}
        a.obs["sub_label"] = a.obs["sub"].map({s: f"{s} {top[s][0] if top[s] else ''}" for s in subs})
        if not cached:
            a.write_h5ad(PROC / f"subclusters_{ct}.h5ad", compression="gzip")

        de = {r: pd.read_csv(TABLES / f"pseudobulk_DE_{r}_{ct}.csv") for r in ["Hippocampus", "Cortex"]}
        gsea = {r: (pd.read_csv(TABLES / f"GSEA_{r}_{ct}.csv") if cached else run_gsea(de[r], gsets, r, ct)) for r in de}
        all_gsea += list(gsea.values())

        # Layout: four rows read left-to-right, top-to-bottom (a | b b ; c | d ; e | f ; g | h).
        fig = plt.figure(figsize=(7.2, 8.6))
        # Two stacked grids so each inter-row gap is sized to what sits in it (UMAP legends vs rotated gene
        # names + volcano legend), leaving no empty bands.
        kw = dict(figure=fig, left=0.085, right=0.985)
        g_top = gridspec.GridSpec(2, 1, height_ratios=[0.92, 0.95], hspace=0.4, top=0.965, bottom=0.565, **kw)
        g_bot = gridspec.GridSpec(2, 1, height_ratios=[0.95, 1.05], hspace=0.72, top=0.475, bottom=0.06, **kw)
        r0 = g_top[0].subgridspec(1, 3, wspace=0.12)
        r1 = g_top[1].subgridspec(1, 2, width_ratios=[1.12, 1.0], wspace=0.22)
        r2 = g_bot[0].subgridspec(1, 2, wspace=0.22)
        r3 = g_bot[1].subgridspec(1, 3, width_ratios=[1.0, 0.2, 0.58], wspace=0.3)

        # a: subclusters
        ax = fig.add_subplot(r0[0, 0]); xy = a.obsm["X_umap"]
        embed(ax, xy, a.obs["sub"].map(pal).values, s=0.6)
        for s in subs:
            m = (a.obs["sub"] == s).values
            ax.text(np.median(xy[m, 0]), np.median(xy[m, 1]), s, fontsize=6, ha="center", va="center",
                    fontweight="bold", bbox=dict(boxstyle="round,pad=0.12", fc="white", ec="none", alpha=0.75))
        ax.set_aspect("equal", adjustable="datalim")
        ax.set_title(f"{ct} subclusters ({a.n_obs:,} nuclei)", pad=4)
        panel_letter(ax, "a", x=-0.1, y=1.02)
        # b: region / condition
        for k, (key, cmap, ttl) in enumerate([("region", REGION_COLORS, "Region"), ("condition", COND_COLORS, "Condition")]):
            ax = fig.add_subplot(r0[0, 1 + k])
            embed(ax, xy, a.obs[key].map(cmap).values, s=0.6, arrows=False)
            ax.set_aspect("equal", adjustable="datalim")
            ax.legend(handles=[Line2D([], [], marker="o", ls="", color=c, ms=3, label=l) for l, c in cmap.items()],
                      loc="upper left", bbox_to_anchor=(-0.02, 0.02), handletextpad=0.1, ncol=2, columnspacing=0.6)
            ax.set_title(ttl, pad=4, loc="center")
            if k == 0:
                panel_letter(ax, "b", x=-0.1, y=1.02)

        # c: marker dot plot
        ax = fig.add_subplot(r1[0, 0])
        genes = []
        for s in subs:
            genes += [g for g in top[s] if g not in genes]
        genes = genes[:24]
        E = pd.DataFrame(a[:, genes].X.toarray(), columns=genes); E["s"] = a.obs["sub"].values
        mean = E.groupby("s", observed=True).mean().loc[subs]; frac = (E.set_index("s") > 0).groupby(level=0, observed=True).mean().loc[subs]
        z = (mean - mean.mean()) / (mean.std() + 1e-9)
        yy, xx = np.meshgrid(np.arange(len(subs)), np.arange(len(genes)), indexing="ij")
        sca = ax.scatter(xx.ravel(), yy.ravel(), s=frac.values.ravel() * 22 + 0.3, c=z.values.ravel(), cmap="RdBu_r",
                         norm=TwoSlopeNorm(0, -2, 2), lw=0.2, edgecolors="0.6")
        ax.set_xticks(range(len(genes))); ax.set_xticklabels(genes, rotation=90, style="italic", fontsize=5)
        ax.set_yticks(range(len(subs))); ax.set_yticklabels(subs)
        for i, s in enumerate(subs):
            ax.get_yticklabels()[i].set_color(pal[s])
        ax.set_xlim(-0.6, len(genes) - 0.4); ax.set_ylim(len(subs) - 0.4, -0.6)
        cax = ax.inset_axes([0.72, 1.04, 0.28, 0.035]); cb = fig.colorbar(sca, cax=cax, orientation="horizontal")
        cb.ax.tick_params(labelsize=4.5, pad=1); cax.xaxis.set_ticks_position("top")
        cax.text(-0.04, 0.5, "z-score", transform=cax.transAxes, ha="right", va="center", fontsize=5)
        ax.set_title("Subcluster marker genes", y=1.02)
        panel_letter(ax, "c", x=-0.1, y=1.06)

        # d: donor-level subcluster proportions
        ax = fig.add_subplot(r1[0, 1])
        cnt = pd.crosstab(a.obs["sample"], a.obs["sub"]); prop = cnt.div(cnt.sum(1), axis=0) * 100
        prop = prop.join(smeta[["region", "condition", "donor"]])
        prop.to_csv(TABLES / f"Fig2_{ct}_subcluster_proportions.csv")
        tests = []
        w = 0.18
        for i, s in enumerate(subs):
            for j, (r_, c_) in enumerate([("Hippocampus", "Control"), ("Hippocampus", "AD"), ("Cortex", "Control"), ("Cortex", "AD")]):
                v = prop[(prop.region == r_) & (prop.condition == c_)][s].values
                x = i + (j - 1.5) * w
                ax.scatter(x + np.random.default_rng(i * 7 + j).uniform(-0.04, 0.04, len(v)), v, s=4,
                           color=COND_COLORS[c_], marker="o" if r_ == "Hippocampus" else "^", lw=0.2, edgecolor="k", zorder=3)
                ax.hlines(np.median(v), x - 0.07, x + 0.07, color="k", lw=0.7, zorder=4)
            for r_ in ["Hippocampus", "Cortex"]:
                q = prop[prop.region == r_]
                p = stats.mannwhitneyu(q[q.condition == "AD"][s], q[q.condition == "Control"][s]).pvalue
                tests.append(dict(cell_type=ct, subcluster=s, region=r_, marker_genes=",".join(top[s]),
                                  AD_median=q[q.condition == "AD"][s].median(), Control_median=q[q.condition == "Control"][s].median(), MWU_p=p))
        t = pd.DataFrame(tests); t["FDR"] = bh(t.MWU_p); t.to_csv(TABLES / f"Fig2_{ct}_subcluster_tests.csv", index=False)
        for i, s in enumerate(subs):
            for k, r_ in enumerate(["Hippocampus", "Cortex"]):
                f_ = t[(t.subcluster == s) & (t.region == r_)].FDR.iloc[0]
                if f_ < 0.05:
                    ax.text(i + (k - 0.5) * 2 * w, prop[subs].values.max() * 1.02, "*", ha="center", fontsize=7)
        ax.set_xticks(range(len(subs))); ax.set_xticklabels(subs)
        for i, s in enumerate(subs):
            ax.get_xticklabels()[i].set_color(pal[s])
        ax.set_ylabel("% of nuclei per sample")
        ax.legend(handles=[Line2D([], [], marker=m, ls="", mfc=COND_COLORS[c], mec="k", mew=0.2, ms=3.5, label=f"{r} {c}")
                           for r, m in [("Hpc", "o"), ("Ctx", "^")] for c in ["Control", "AD"]],
                  ncol=2, loc="upper right", fontsize=5, handletextpad=0.1, columnspacing=0.6)
        ax.set_title("Subcluster proportions per sample (* FDR < 0.05)", pad=6)
        panel_letter(ax, "d", x=-0.17, y=1.06)

        # e,f: volcanoes
        for k, r_ in enumerate(["Hippocampus", "Cortex"]):
            ax = fig.add_subplot(r2[0, k])
            n0, n1 = de[r_].n_control.iloc[0], de[r_].n_AD.iloc[0]
            volcano(ax, de[r_], f"{r_} ({n0} Control vs {n1} AD)")
            panel_letter(ax, "ef"[k], x=-0.15, y=1.06)
            if k == 0:
                vax = ax
        # shared legend for volcano glyphs, under panel e/f
        vax.legend(handles=[Line2D([], [], marker="o", ls="", color=COND_COLORS["AD"], ms=3, label="Robust up in AD"),
                            Line2D([], [], marker="o", ls="", color=COND_COLORS["Control"], ms=3, label="Robust down in AD"),
                            Line2D([], [], marker="o", ls="", mfc="none", mec="0.45", ms=3.5, label="Artefact-flagged (ambient / 17q21.31 / sex)"),
                            Line2D([], [], marker="o", ls="", color="0.75", ms=3, label="Not significant")],
                   loc="upper left", bbox_to_anchor=(0.0, -0.22), ncol=2, fontsize=5.5, handletextpad=0.1, columnspacing=1.2)

        # g: donor heatmap of robust genes (union, top by p)
        ax = fig.add_subplot(r3[0, 0])
        rob = pd.concat([de[r][de[r].robust] for r in de]).sort_values("pvalue").drop_duplicates("gene")
        hg = rob.gene.head(15).tolist()
        mats, cols = [], []
        for r_ in ["Hippocampus", "Cortex"]:
            C = pd.read_csv(PROC / "pseudobulk" / f"{r_}_{ct}_counts.csv.gz", index_col=0)
            n = pd.read_csv(PROC / "pseudobulk" / f"{r_}_{ct}_ncells.csv", index_col=0).iloc[:, 0]
            C = C.loc[n.index[n >= 30]]
            L = np.log2(C.div(C.sum(1), axis=0) * 1e6 + 1)
            L = L.reindex(columns=hg).fillna(0)
            o = smeta.loc[L.index].sort_values(["condition", "donor"], ascending=[False, True]).index
            mats.append(L.loc[o]); cols += [(r_, smeta.loc[s, "condition"]) for s in o]
        M = pd.concat(mats)
        Z = ((M - M.mean()) / (M.std() + 1e-9)).T
        if len(hg):
            im = ax.imshow(Z.values, aspect="auto", cmap="RdBu_r", norm=TwoSlopeNorm(0, -2.5, 2.5), interpolation="nearest")
            ax.set_yticks(range(len(hg))); ax.set_yticklabels(hg, style="italic", fontsize=5)
            ax.set_xticks([])
            for x, (r_, c_) in enumerate(cols):
                ax.add_patch(plt.Rectangle((x - 0.5, -1.9), 1, 0.8, color=COND_COLORS[c_], clip_on=False))
                ax.add_patch(plt.Rectangle((x - 0.5, -1.0), 1, 0.8, color=REGION_COLORS[r_], clip_on=False))
            nH = sum(1 for r_, _ in cols if r_ == "Hippocampus")
            ax.axvline(nH - 0.5, color="k", lw=0.8)
            ax.set_ylim(len(hg) - 0.5, -2.0)
            cax = ax.inset_axes([0.0, -0.2, 0.3, 0.045]); cb = fig.colorbar(im, cax=cax, orientation="horizontal")
            cb.set_label("z-scored log$_2$ CPM", fontsize=5, labelpad=1); cb.ax.tick_params(labelsize=4.5)
            ax.text(1.0, -0.12, f"Samples: Hpc n = {nH} | Ctx n = {len(cols) - nH}\nbars: condition (top), region (bottom)",
                    transform=ax.transAxes, ha="right", va="top", fontsize=5)
        else:
            ax.text(0.5, 0.5, "No robust DEGs", ha="center", va="center", transform=ax.transAxes); ax.axis("off")
        ax.set_title("Robust DEGs per sample", pad=14)
        panel_letter(ax, "g", x=-0.15, y=1.1)

        # h: GSEA dot plot (top 6 per region by FDR, union)
        ax = fig.add_subplot(r3[0, 1])
        G = pd.concat(gsea.values())
        sel = (G[G.FDR < 0.25].sort_values("FDR").groupby("region").head(6).term.unique().tolist())[:12]
        if sel:
            P = G[G.term.isin(sel)].pivot_table(index="term", columns="region", values="NES").reindex(columns=["Hippocampus", "Cortex"])
            Q = G[G.term.isin(sel)].pivot_table(index="term", columns="region", values="FDR").reindex(columns=["Hippocampus", "Cortex"])
            P = P.loc[P.mean(1).sort_values().index]; Q = Q.loc[P.index]
            yy, xx = np.meshgrid(np.arange(len(P)), np.arange(2), indexing="ij")
            sz = (-np.log10(Q.values.clip(1e-4, 1)) * 18 + 2).ravel()
            s_ = ax.scatter(xx.ravel(), yy.ravel(), s=sz, c=P.values.ravel(), cmap="RdBu_r", norm=TwoSlopeNorm(0, -3, 3),
                            edgecolors=np.where(Q.values.ravel() < 0.05, "k", "none"), lw=0.5)
            ax.set_yticks(range(len(P))); ax.set_yticklabels([short(t) for t in P.index], fontsize=5)
            ax.set_xticks([0, 1]); ax.set_xticklabels(["Hpc", "Ctx"]); ax.set_xlim(-0.55, 1.55); ax.set_ylim(len(P) - 0.4, -0.6)
            ax.yaxis.tick_right(); ax.spines["left"].set_visible(False); ax.spines["right"].set_visible(True)
            cax = ax.inset_axes([0.0, -0.2, 1.0, 0.045]); cb = fig.colorbar(s_, cax=cax, orientation="horizontal")
            cb.set_label("NES", fontsize=5, labelpad=1); cb.ax.tick_params(labelsize=4.5)
            ax.text(1.25, -0.2, "dot size: −log$_{10}$ FDR\nblack ring: FDR < 0.05", transform=ax.transAxes, fontsize=5, ha="left", va="top")
        else:
            ax.text(0.5, 0.5, "No gene set at FDR < 0.25", ha="center", va="center", transform=ax.transAxes); ax.axis("off")
        ax.set_title("Pathways (pre-ranked GSEA)", pad=14, loc="left")
        panel_letter(ax, "h", x=-0.55, y=1.1)

        save_fig(fig, OUT / f"Fig2{LETTER[ct]}_{ct}")
        plt.close(fig)
    pd.concat(all_gsea).to_csv(TABLES / "GSEA_all.csv", index=False)


if __name__ == "__main__":
    main()
