"""Step 7 - Figure 4: donor-level ligand-receptor (L-R) communication among the 5 NVU cell types.

Resource: LIANA 'consensus' L-R resource (complexes resolved by the minimum subunit).
Per sample (= donor x region) and per sender/receiver pair, an interaction is scored as
    score = sqrt( mean_expr(ligand | sender) * mean_expr(receptor | receiver) )
on log-normalised expression (LIANA 'lr_means'-style geometric mean), and is called 'expressed'
when ligand and receptor are each detected in >= 10 % of nuclei of their cell type.
Differential communication = per interaction, OLS of log2(score + 0.01) on condition + sex +
vascular_yield across samples (unit = sample), BH-FDR within region. Only interactions expressed
in >= 50 % of samples of either group are tested; interactions whose ligand or receptor is
artefact-flagged in the relevant cell type are reported but excluded from the summary panels.

Outputs: figures/Fig4_Intercellular_Communication/Fig4_communication.{pdf,png}
         tables/LR_scores_per_sample.csv.gz, tables/LR_differential_<Region>.csv
"""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import TwoSlopeNorm
import liana as li

sys.path.insert(0, str(Path(__file__).parent))
from config import *  # noqa

OUT = FIGS / "Fig4_Intercellular_Communication"
SHORT = {"Microglia": "Micro", "Endothelial": "Endo", "Pericyte": "Peri", "Astrocyte": "Astro", "Oligodendrocyte": "Oligo"}


def main():
    set_style()
    ad = sc.read_h5ad(PROC / "integrated.h5ad")
    ad = ad[ad.obs.cell_type.isin(NVU5)].copy()
    smeta = pd.read_csv(TABLES / "sample_metadata.csv").set_index("sample")
    res = li.resource.select_resource("consensus")
    # APP-SORL1 is a CellChatDB pair absent from the consensus resource; added so that the
    # legacy headline claim (APP->SORL1 hyper-signalling into microglia) can be tested directly
    res = pd.concat([res, pd.DataFrame({"ligand": ["APP"], "receptor": ["SORL1"]})], ignore_index=True).drop_duplicates()
    res["lig"] = res.ligand.str.split("_"); res["rec"] = res.receptor.str.split("_")
    genes = sorted({g for l in res.lig for g in l} | {g for r in res.rec for g in r})
    genes = [g for g in genes if g in ad.var_names]
    res = res[res.lig.apply(lambda l: all(g in ad.var_names for g in l)) & res.rec.apply(lambda r: all(g in ad.var_names for g in r))]
    X = ad[:, genes].X.tocsr() if sp.issparse(ad.X) else sp.csr_matrix(ad[:, genes].X)
    key = (ad.obs["sample"].astype(str) + "|" + ad.obs.cell_type.astype(str)).values
    u, inv = np.unique(key, return_inverse=True)
    G = sp.csr_matrix((np.ones(len(key)), (inv, np.arange(len(key)))), shape=(len(u), len(key)))
    n = np.asarray(G.sum(1)).ravel()
    mean = pd.DataFrame((G @ X).toarray() / n[:, None], index=u, columns=genes)
    frac = pd.DataFrame((G @ (X > 0).astype(float)).toarray() / n[:, None], index=u, columns=genes)
    ok_groups = set(u[n >= 30])
    # complex-level (min over subunits) mean and detection fraction, vectorised over groups
    cplx = sorted(set(res.ligand) | set(res.receptor))
    CM = pd.DataFrame({c: mean[c.split("_")].min(1) for c in cplx})
    CF = pd.DataFrame({c: frac[c.split("_")].min(1) for c in cplx})
    li_idx = CM.columns.get_indexer(res.ligand); re_idx = CM.columns.get_indexer(res.receptor)
    CMv, CFv, gidx = CM.values, CF.values, {g: i for i, g in enumerate(CM.index)}

    recs = []
    for s in smeta.index:
        for snd in NVU5:
            for rcv in NVU5:
                a, b = f"{s}|{snd}", f"{s}|{rcv}"
                if a not in ok_groups or b not in ok_groups:
                    continue
                ia, ib = gidx[a], gidx[b]
                lm, rm = CMv[ia, li_idx], CMv[ib, re_idx]
                lf, rf = CFv[ia, li_idx], CFv[ib, re_idx]
                recs.append(pd.DataFrame(dict(sample=s, sender=snd, receiver=rcv, ligand=res.ligand.values,
                                              receptor=res.receptor.values, score=np.sqrt(lm * rm),
                                              expressed=(lf >= 0.1) & (rf >= 0.1))))
    L = pd.concat(recs)
    L = L.join(smeta[["region", "condition", "sex", "vascular_yield"]], on="sample")
    L["interaction"] = L.sender + ">" + L.receiver + ":" + L.ligand + ">" + L.receptor
    L.to_csv(TABLES / "LR_scores_per_sample.csv.gz", index=False)

    de = pd.read_csv(TABLES / "pseudobulk_DE_all.csv.gz")
    flag = de[de.flag_contamination | de.flag_genotype_sex]
    flagset = set(zip(flag.region, flag.cell_type, flag.gene))

    diff = {}
    for r_ in ["Hippocampus", "Cortex"]:
        d = L[L.region == r_]
        grp = d.groupby(["interaction", "condition"]).expressed.mean().unstack()
        keep = grp.index[(grp.max(1) >= 0.5)]
        d = d[d.interaction.isin(keep)].copy()
        d["y"] = np.log2(d.score + 0.01)
        d["vy"] = (d.vascular_yield - d.vascular_yield.mean()) / d.vascular_yield.std()
        out = []
        for it, g in d.groupby("interaction"):
            if g.y.std() == 0:
                continue
            m = smf.ols("y ~ C(condition, Treatment('Control')) + sex + vy", data=g).fit()
            k = [c for c in m.params.index if c.startswith("C(condition")][0]
            snd, rest = it.split(">", 1); rcv, lr = rest.split(":", 1); lig, rec = lr.split(">")
            fl = any((r_, snd, x) in flagset for x in lig.split("_")) or any((r_, rcv, x) in flagset for x in rec.split("_"))
            out.append(dict(interaction=it, sender=snd, receiver=rcv, ligand=lig, receptor=rec, log2FC=m.params[k],
                            se=m.bse[k], p=m.pvalues[k], frac_expr_AD=grp.loc[it, "AD"], frac_expr_Ctrl=grp.loc[it, "Control"],
                            mean_score=g.score.mean(), flagged=fl))
        o = pd.DataFrame(out); o["FDR"] = bh(o.p); o = o.sort_values("p")
        o.to_csv(TABLES / f"LR_differential_{r_}.csv", index=False)
        diff[r_] = o
        print(r_, "tested", len(o), "FDR<0.05", int((o.FDR < 0.05).sum()), "FDR<0.05 unflagged", int(((o.FDR < 0.05) & ~o.flagged).sum()))

    fig = plt.figure(figsize=(7.2, 7.4))
    gs = gridspec.GridSpec(3, 4, figure=fig, height_ratios=[1, 1, 1.25], hspace=0.85, wspace=0.9)

    # a: number of expressed interactions per sample
    ax = fig.add_subplot(gs[0, 0:2])
    E = L[L.expressed].groupby(["sample"]).size().rename("n").to_frame().join(smeta[["region", "condition", "vascular_fraction"]])
    for i, r_ in enumerate(["Hippocampus", "Cortex"]):
        for j, c_ in enumerate(["Control", "AD"]):
            v = E[(E.region == r_) & (E.condition == c_)].n.values
            x = i * 2.6 + j
            ax.scatter(x + np.random.default_rng(i + j).uniform(-0.15, 0.15, len(v)), v, s=10, color=COND_COLORS[c_], edgecolor="k", lw=0.3, zorder=3)
            ax.hlines(np.median(v), x - 0.3, x + 0.3, color="k", lw=0.9)
        q = E[E.region == r_]
        from scipy.stats import mannwhitneyu
        p = mannwhitneyu(q[q.condition == "AD"].n, q[q.condition == "Control"].n).pvalue
        ax.text(i * 2.6 + 0.5, E.n.max() * 1.04, f"P = {p:.2g}", ha="center", fontsize=5.5)
    ax.set_xticks([0, 1, 2.6, 3.6]); ax.set_xticklabels(["Ctrl", "AD", "Ctrl", "AD"])
    ax.text(0.5, -0.2, "Hippocampus", transform=ax.get_xaxis_transform(), ha="center", fontsize=6)
    ax.text(3.1, -0.2, "Cortex", transform=ax.get_xaxis_transform(), ha="center", fontsize=6)
    ax.set_ylabel("Expressed L–R interactions\nper sample (5 NVU types)")
    ax.set_title("Expressed L–R repertoire per sample")
    panel_letter(ax, "a", x=-0.25)

    # b,c: sender x receiver median effect + count of nominally changed interactions
    for k, r_ in enumerate(["Hippocampus", "Cortex"]):
        ax = fig.add_subplot(gs[0, 2:4] if k == 0 else gs[1, 0:2])
        o = diff[r_][~diff[r_].flagged]
        med = o.groupby(["sender", "receiver"]).log2FC.median().unstack().reindex(index=NVU5, columns=NVU5)
        nom = o[o.p < 0.01].groupby(["sender", "receiver"]).size().unstack().reindex(index=NVU5, columns=NVU5).fillna(0)
        tot = o.groupby(["sender", "receiver"]).size().unstack().reindex(index=NVU5, columns=NVU5).fillna(0)
        im = ax.imshow(med.values, cmap="RdBu_r", norm=TwoSlopeNorm(0, -0.6, 0.6))
        for i in range(5):
            for j in range(5):
                ax.text(j, i, f"{int(nom.values[i, j])}/{int(tot.values[i, j])}", ha="center", va="center", fontsize=4.6)
        ax.set_xticks(range(5)); ax.set_xticklabels([SHORT[c] for c in NVU5], rotation=45, ha="right")
        ax.set_yticks(range(5)); ax.set_yticklabels([SHORT[c] for c in NVU5])
        ax.set_xlabel("Receiver"); ax.set_ylabel("Sender")
        cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.02); cb.set_label("Median log$_2$FC", fontsize=5, labelpad=1); cb.ax.tick_params(labelsize=5)
        nsig = int((o.FDR < 0.05).sum())
        ax.set_title(f"{r_}: {nsig} of {len(o):,} pairs at FDR < 0.05\n(cell text: nominal P < 0.01 / tested)", fontsize=6)
        panel_letter(ax, "bc"[k], x=-0.35)

    # d: per-donor APP->SORL1 into microglia (legacy headline claim), both regions
    ax = fig.add_subplot(gs[1, 2:4])
    d = L[(L.ligand == "APP") & (L.receptor == "SORL1") & (L.receiver == "Microglia")]
    d = d.groupby(["sample", "region", "condition"]).score.mean().reset_index()
    for i, r_ in enumerate(["Hippocampus", "Cortex"]):
        for j, c_ in enumerate(["Control", "AD"]):
            v = d[(d.region == r_) & (d.condition == c_)].score.values
            x = i * 2.6 + j
            ax.scatter(x + np.random.default_rng(5 + i + j).uniform(-0.15, 0.15, len(v)), v, s=10, color=COND_COLORS[c_], edgecolor="k", lw=0.3, zorder=3)
            if len(v):
                ax.hlines(np.median(v), x - 0.3, x + 0.3, color="k", lw=0.9)
        o = diff[r_]; sel = o[(o.ligand == "APP") & (o.receptor == "SORL1") & (o.receiver == "Microglia")]
        txt = f"min FDR = {sel.FDR.min():.2g}" if len(sel) else "not expressed"
        ax.text(i * 2.6 + 0.5, 1.02, txt, ha="center", fontsize=5.5, transform=ax.get_xaxis_transform())
    ax.set_xticks([0, 1, 2.6, 3.6]); ax.set_xticklabels(["Ctrl", "AD", "Ctrl", "AD"])
    ax.text(0.5, -0.2, "Hippocampus", transform=ax.get_xaxis_transform(), ha="center", fontsize=6)
    ax.text(3.1, -0.2, "Cortex", transform=ax.get_xaxis_transform(), ha="center", fontsize=6)
    ax.set_ylabel("APP→SORL1 score\n(into microglia)")
    ax.set_title("Legacy APP→SORL1 claim, per donor", pad=12)
    panel_letter(ax, "d", x=-0.25)

    # e: top differential interactions (hippocampus discovery) with cortex effect alongside
    ax = fig.add_subplot(gs[2, :])
    o = diff["Hippocampus"][~diff["Hippocampus"].flagged].sort_values("p").head(20)
    oc = diff["Cortex"].set_index("interaction")
    y = np.arange(len(o))
    ax.errorbar(o.log2FC, y - 0.15, xerr=1.96 * o.se, fmt="o", ms=3, color=REGION_COLORS["Hippocampus"], lw=0.6, label="Hippocampus (discovery)")
    c = oc.reindex(o.interaction)
    ax.errorbar(c.log2FC, y + 0.15, xerr=1.96 * c.se, fmt="^", ms=3, color=REGION_COLORS["Cortex"], lw=0.6, label="Cortex (same donors)")
    ax.axvline(0, color="0.6", lw=0.5)
    labs = [f"{SHORT[s]}→{SHORT[r]}  {l}→{rc}" for s, r, l, rc in zip(o.sender, o.receiver, o.ligand, o.receptor)]
    ax.set_yticks(y); ax.set_yticklabels(labs, fontsize=5); ax.invert_yaxis()
    ax.set_xlabel("log$_2$ fold change of interaction score (AD/Control) ± 95% CI")
    ax.legend(loc="lower left", fontsize=5)
    ax.set_title(f"Top 20 hippocampal L–R changes by nominal P (all FDR ≥ {o.FDR.min():.2f}); OLS with sex and yield")
    panel_letter(ax, "e", x=-0.2)

    save_fig(fig, OUT / "Fig4_communication")


if __name__ == "__main__":
    main()
