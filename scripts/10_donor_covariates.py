"""Step 10 - donor clinical covariates (Yang et al. 2022, Supplementary Table 1) and sensitivity analyses.

Input : data/41586_2021_4369_MOESM2_ESM.xlsx  (Nature 603:885, Supplementary Table 1; one row per sample)
        tables/sample_metadata.csv, tables/pseudobulk_DE_*.csv, data/processed/pseudobulk/*
Output: tables/donor_clinical_metadata.csv      donor-level age, sex, PMI, Braak, CERAD, APOE, pathology flags
        tables/sample_metadata.csv               (updated with the clinical columns)
        tables/covariate_balance.csv             AD vs control comparison of each covariate
        tables/hypoxia_module_score_hpc.csv      per-sample hypoxia-module score, 5 NVU cell types
        tables/hypoxia_module_covariate_models.csv  condition effect on the score under each adjustment
        tables/pseudobulk_DE_sensitivity_clinical.csv  hippocampal DE re-fitted with age + PMI and
                                                  without the donors with acute hypoxic neuropathology
        figures/Quality_Control/FigS2_donor_covariates.{pdf,png}

Mapping: the table's 'NCI k' / 'AD k' IDs equal the k in the GEO sample titles ('Control k Hpc' /
'AD k Hpc'); this is verified below by (i) sex inferred from XIST/chrY expression for all 17 donors and
(ii) the set of donors with a superior-frontal-cortex sample.

No agonal-factor score or brain pH is reported. Acute hypoxic neuropathology noted by the
neuropathologist ('acute hypoxic damage', 'acute hypoxic/ischemic changes', 'hypoxic-ischemic
encephalopathy') is used as the available proxy for terminal hypoxia.
"""
import re, sys, warnings, logging
from pathlib import Path

import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.lines import Line2D
from scipy import stats

sys.path.insert(0, str(Path(__file__).parent))
from config import *  # noqa

warnings.filterwarnings("ignore"); logging.disable(logging.WARNING)
SUPP = DATA / "41586_2021_4369_MOESM2_ESM.xlsx"
HYPOXIA_MODULE = ["ADAMTS9", "VEGFA", "HILPDA", "SLC2A3", "PDK1", "BNIP3L", "ELL2", "HIF1A"]
ACUTE_HYPOXIA = r"acute hypoxic|hypoxic[- /]ischemic|hypoxic damage"
INFARCT = r"infarct"
OUT = FIGS / "Quality_Control"


def load_clinical():
    raw = pd.read_excel(SUPP)
    raw.columns = [str(c).strip() for c in raw.columns]
    txt_cols = [c for c in raw.columns if c not in ["Patient ID", "Age", "Sex", "Brain region profiled",
                                                    "PMI (hours)", "Clinical AD or not", "AD on histopath",
                                                    "Braak", "CERAD", "ApoE genotype", "Overall brain weight"]]
    raw["notes"] = raw[txt_cols].apply(lambda r: " | ".join(str(v) for v in r if pd.notna(v)), axis=1)
    hpc = raw[raw["Brain region profiled"].str.strip() == "Hippocampus"].copy()
    ctx_ids = set(raw[raw["Brain region profiled"].str.contains("cortex", case=False)]["Patient ID"].str.strip())
    pid = hpc["Patient ID"].str.strip()
    out = pd.DataFrame({
        "patient_id": pid.values,
        "condition": np.where(pid.str.startswith("AD").values, "AD", "Control"),
        "title_index": pid.str.extract(r"(\d+)")[0].astype(int).values,
        "age": hpc["Age"].values,
        "sex_reported": hpc["Sex"].str.strip().values,
        "pmi_h": hpc["PMI (hours)"].astype(float).values,
        "braak": hpc["Braak"].astype(str).str.strip().replace({"nan": np.nan, "NR": np.nan}).values,
        "cerad": hpc["CERAD"].astype(str).str.strip().replace({"nan": np.nan, "NR": np.nan}).values,
        "apoe": hpc["ApoE genotype"].str.strip().values,
        "brain_weight_g": hpc["Overall brain weight"].values,
        "has_cortex_sample": pid.isin(ctx_ids).values,
        "acute_hypoxic_pathology": hpc["notes"].str.contains(ACUTE_HYPOXIA, case=False, regex=True).values,
        "infarct_noted": hpc["notes"].str.contains(INFARCT, case=False, regex=True).values,
    })
    out["apoe4_dose"] = out.apoe.str.count("4")
    return out


def module_scores(smeta):
    """Per-sample hypoxia-module score (mean z of log2 CPM across the module genes) per NVU cell type."""
    rows = []
    for region in ["Hippocampus", "Cortex"]:
        for ct in NVU5:
            C = pd.read_csv(PROC / "pseudobulk" / f"{region}_{ct}_counts.csv.gz", index_col=0)
            n = pd.read_csv(PROC / "pseudobulk" / f"{region}_{ct}_ncells.csv", index_col=0).iloc[:, 0]
            C = C.loc[n.index[n >= 30]]
            L = np.log2(C.div(C.sum(1), axis=0) * 1e6 + 1)
            g = [x for x in HYPOXIA_MODULE if x in L]
            z = ((L[g] - L[g].mean()) / L[g].std()).mean(1)
            rows.append(pd.DataFrame({"sample": z.index, "region": region, "cell_type": ct, "score": z.values,
                                      "n_genes": len(g)}))
    return pd.concat(rows)


def main():
    set_style()
    smeta = pd.read_csv(TABLES / "sample_metadata.csv")
    clin = load_clinical()

    # ---- map to GEO samples via title index and verify
    smeta["title_index"] = smeta.title.str.extract(r"(\d+)")[0].astype(int)
    hmeta = smeta[smeta.region == "Hippocampus"][["donor", "condition", "title_index", "sex"]]
    clin = clin.merge(hmeta, on=["condition", "title_index"], how="left", validate="one_to_one")
    assert clin.donor.notna().all(), "unmapped donors"
    clin["sex_concordant"] = clin.sex_reported == clin.sex
    ctx_donors = set(smeta[smeta.region == "Cortex"].donor)
    clin["cortex_concordant"] = clin.has_cortex_sample == clin.donor.isin(ctx_donors)
    print("sex concordant:", int(clin.sex_concordant.sum()), "/", len(clin),
          "| cortex subset concordant:", int(clin.cortex_concordant.sum()), "/", len(clin))
    assert clin.sex_concordant.all() and clin.cortex_concordant.all(), "donor mapping failed verification"
    clin = clin.drop(columns=["title_index", "sex"])
    clin.to_csv(TABLES / "donor_clinical_metadata.csv", index=False)

    keep = ["donor", "age", "pmi_h", "braak", "cerad", "apoe", "apoe4_dose", "acute_hypoxic_pathology",
            "infarct_noted", "brain_weight_g"]
    smeta = smeta.drop(columns=[c for c in keep[1:] + ["title_index"] if c in smeta]).merge(clin[keep], on="donor", how="left")
    smeta.to_csv(TABLES / "sample_metadata.csv", index=False)

    # ---- covariate balance (donor level)
    bal = []
    for v in ["age", "pmi_h", "apoe4_dose", "brain_weight_g"]:
        a, b = clin[clin.condition == "AD"][v], clin[clin.condition == "Control"][v]
        bal.append(dict(covariate=v, AD_median=a.median(), Control_median=b.median(), MWU_p=stats.mannwhitneyu(a, b).pvalue))
    for v in ["acute_hypoxic_pathology", "infarct_noted"]:
        t = pd.crosstab(clin.condition, clin[v]).reindex(columns=[False, True], fill_value=0)
        bal.append(dict(covariate=v, AD_median=int(t.loc["AD", True]), Control_median=int(t.loc["Control", True]),
                        MWU_p=stats.fisher_exact(t.values)[1]))
    bal = pd.DataFrame(bal); bal.to_csv(TABLES / "covariate_balance.csv", index=False)
    print(bal.round(3).to_string())

    # ---- hypoxia-module score vs covariates
    S = module_scores(smeta).merge(smeta, on=["sample", "region"], how="left")
    S["vy"] = S.vascular_yield
    S.to_csv(TABLES / "hypoxia_module_score_hpc.csv", index=False)
    models = {"condition + sex + yield (primary)": "score ~ C(condition, Treatment('Control')) + sex + vy",
              "+ age": "score ~ C(condition, Treatment('Control')) + sex + vy + age",
              "+ PMI": "score ~ C(condition, Treatment('Control')) + sex + vy + pmi_h",
              "+ APOE4 dose": "score ~ C(condition, Treatment('Control')) + sex + vy + apoe4_dose",
              "+ acute hypoxic pathology": "score ~ C(condition, Treatment('Control')) + sex + vy + acute_hypoxic_pathology",
              "+ age + PMI": "score ~ C(condition, Treatment('Control')) + sex + vy + age + pmi_h",
              "excluding acute-hypoxic donors": "score ~ C(condition, Treatment('Control')) + sex + vy"}
    res = []
    for ct in NVU5:
        d = S[(S.region == "Hippocampus") & (S.cell_type == ct)].copy()
        d["acute_hypoxic_pathology"] = d.acute_hypoxic_pathology.astype(int)
        for name, f in models.items():
            dd = d[~d.acute_hypoxic_pathology.astype(bool)] if name.startswith("excluding") else d
            m = smf.ols(f, data=dd).fit()
            k = [c for c in m.params.index if c.startswith("C(condition")][0]
            ci = m.conf_int().loc[k]
            res.append(dict(cell_type=ct, model=name, n=int(m.nobs), AD_effect=m.params[k], ci_low=ci[0], ci_high=ci[1], p=m.pvalues[k]))
        for v in ["age", "pmi_h", "apoe4_dose"]:
            r = stats.spearmanr(d[v], d.score)
            res.append(dict(cell_type=ct, model=f"spearman score~{v}", n=len(d), AD_effect=r.statistic, ci_low=np.nan, ci_high=np.nan, p=r.pvalue))
    R_ = pd.DataFrame(res); R_.to_csv(TABLES / "hypoxia_module_covariate_models.csv", index=False)
    print(R_[R_.model.str.startswith(("condition", "+ age + PMI", "excluding"))].round(4).to_string())

    # ---- effect sizes for validation planning (per-donor Cohen's d, hippocampus; inflated by selection)
    from statsmodels.stats.power import TTestIndPower
    pw = []
    for ct in ["Pericyte", "Endothelial", "Astrocyte", "Microglia"]:
        C = pd.read_csv(PROC / "pseudobulk" / f"Hippocampus_{ct}_counts.csv.gz", index_col=0)
        n = pd.read_csv(PROC / "pseudobulk" / f"Hippocampus_{ct}_ncells.csv", index_col=0).iloc[:, 0]
        C = C.loc[n.index[n >= 30]]
        L = np.log2(C.div(C.sum(1), axis=0) * 1e6 + 1)
        cond = smeta.set_index("sample").loc[L.index, "condition"]
        for g in [x for x in HYPOXIA_MODULE if x in L]:
            a, b = L.loc[cond == "AD", g], L.loc[cond == "Control", g]
            sp = np.sqrt(((len(a) - 1) * a.var() + (len(b) - 1) * b.var()) / (len(a) + len(b) - 2))
            dd = (b.mean() - a.mean()) / sp
            pw.append(dict(cell_type=ct, gene=g, cohen_d=round(dd, 2),
                           n_per_group_80pct=int(np.ceil(TTestIndPower().solve_power(effect_size=abs(dd), alpha=0.05, power=0.8)))))
    pd.DataFrame(pw).to_csv(TABLES / "validation_power_estimates.csv", index=False)

    # ---- DE sensitivity: hippocampus with age + PMI, and without acute-hypoxic donors
    from importlib import import_module
    de_mod = import_module("04_pseudobulk_de")
    sens = []
    hm = smeta[smeta.region == "Hippocampus"].set_index("sample")
    for ct in NVU5:
        prim = pd.read_csv(TABLES / f"pseudobulk_DE_Hippocampus_{ct}.csv").set_index("gene")
        C = pd.read_csv(PROC / "pseudobulk" / f"Hippocampus_{ct}_counts.csv.gz", index_col=0)
        n = pd.read_csv(PROC / "pseudobulk" / f"Hippocampus_{ct}_ncells.csv", index_col=0).iloc[:, 0]
        keep_s = n.index[n >= 30]
        C = C.loc[keep_s, prim.index]
        meta = hm.loc[keep_s, ["condition", "sex", "vascular_yield", "age", "pmi_h", "acute_hypoxic_pathology"]].copy()
        meta["condition"] = pd.Categorical(meta.condition, ["Control", "AD"])
        for v in ["vascular_yield", "age", "pmi_h"]:
            meta[v] = (meta[v] - meta[v].mean()) / meta[v].std()
        for name, design, mask in [("+ age + PMI", "~sex + vascular_yield + age + pmi_h + condition", np.ones(len(meta), bool)),
                                   ("excluding acute-hypoxic donors", "~sex + vascular_yield + condition", ~meta.acute_hypoxic_pathology.astype(bool).values)]:
            mm = meta[mask].copy()
            for v in ["vascular_yield", "age", "pmi_h"]:
                mm[v] = (mm[v] - mm[v].mean()) / mm[v].std()
            r = de_mod.run_deseq(C[mask], mm[["condition", "sex", "vascular_yield", "age", "pmi_h"]], design)
            rob0 = prim.index[prim.robust]
            sig = (r.padj < 0.05) & (r.log2FoldChange.abs() >= 0.5)
            same = np.sign(r.loc[rob0, "log2FoldChange"]) == np.sign(prim.loc[rob0, "log2FoldChange"])
            sens.append(dict(cell_type=ct, model=name, n_samples=int(mask.sum()), primary_robust=len(rob0),
                             retained_FDR05=int(sig.reindex(rob0).fillna(False).sum()),
                             same_sign_pct=round(100 * same.mean(), 1) if len(rob0) else np.nan,
                             lfc_correlation_all_genes=round(np.corrcoef(r.log2FoldChange, prim.loc[r.index, "log2FoldChange"])[0, 1], 3),
                             module_genes_FDR05=",".join([g for g in HYPOXIA_MODULE if g in r.index and r.loc[g, "padj"] < 0.05])))
            print(sens[-1], flush=True)
    pd.DataFrame(sens).to_csv(TABLES / "pseudobulk_DE_sensitivity_clinical.csv", index=False)

    # ---- Figure S2
    fig = plt.figure(figsize=(7.2, 6.2))
    gs = gridspec.GridSpec(2, 3, figure=fig, hspace=0.45, wspace=0.55, height_ratios=[1, 1.15])
    for k, (v, lab) in enumerate([("age", "Age at death (y)"), ("pmi_h", "Post-mortem interval (h)")]):
        ax = fig.add_subplot(gs[0, k])
        for j, c_ in enumerate(["Control", "AD"]):
            q = clin[clin.condition == c_]
            x = j + np.random.default_rng(j + k).uniform(-0.15, 0.15, len(q))
            ax.scatter(x, q[v], s=14, color=COND_COLORS[c_], edgecolor="k", lw=0.3, zorder=3)
            ax.scatter(x[q.acute_hypoxic_pathology.values], q[v][q.acute_hypoxic_pathology.values], s=60, facecolor="none",
                       edgecolor="k", lw=0.8, zorder=4)
            ax.hlines(q[v].median(), j - 0.3, j + 0.3, color="k", lw=0.9)
        p = bal.set_index("covariate").loc[v, "MWU_p"]
        ax.set_xticks([0, 1]); ax.set_xticklabels(["Control", "AD"]); ax.set_ylabel(lab)
        ax.set_title(f"{lab.split(' (')[0]}: P = {p:.2f}", pad=4)
        panel_letter(ax, "ab"[k], x=-0.3)
    ax = fig.add_subplot(gs[0, 2])
    d = S[(S.region == "Hippocampus") & (S.cell_type == "Pericyte")]
    for c_ in ["Control", "AD"]:
        q = d[d.condition == c_]
        ax.scatter(q.pmi_h, q.score, s=16, color=COND_COLORS[c_], edgecolor="k", lw=0.3, label=c_, zorder=3)
        hq = q[q.acute_hypoxic_pathology.astype(bool)]
        ax.scatter(hq.pmi_h, hq.score, s=70, facecolor="none", edgecolor="k", lw=0.8, zorder=4)
    rr = R_[(R_.cell_type == "Pericyte") & (R_.model == "spearman score~pmi_h")].iloc[0]
    ax.set_xlabel("Post-mortem interval (h)"); ax.set_ylabel("Pericyte hypoxia-module score")
    ax.set_title(f"Score vs PMI: ρ = {rr.AD_effect:.2f}, P = {rr.p:.2f}", pad=4)
    ax.legend(handles=[Line2D([], [], marker="o", ls="", color=COND_COLORS[c], mec="k", mew=0.3, ms=4, label=c) for c in COND_COLORS]
              + [Line2D([], [], marker="o", ls="", mfc="none", mec="k", ms=6, label="acute hypoxic pathology")], fontsize=5, loc="upper right")
    panel_letter(ax, "c", x=-0.3)

    ax = fig.add_subplot(gs[1, 0:2])
    mods = [m for m in models]
    cols = {ct: CT_COLORS[ct] for ct in NVU5}
    for i, ct in enumerate(NVU5):
        q = R_[(R_.cell_type == ct) & R_.model.isin(mods)].set_index("model").loc[mods]
        y = np.arange(len(mods)) + (i - 2) * 0.14
        ax.errorbar(q.AD_effect, y, xerr=[q.AD_effect - q.ci_low, q.ci_high - q.AD_effect], fmt="o", ms=3, lw=0.6, color=cols[ct], label=ct)
    ax.axvline(0, color="0.5", lw=0.6)
    ax.set_yticks(range(len(mods))); ax.set_yticklabels(mods); ax.invert_yaxis()
    ax.set_xlabel("AD − Control difference in hippocampal hypoxia-module score (± 95% CI)")
    ax.legend(fontsize=5, loc="upper center", bbox_to_anchor=(0.45, -0.17), ncol=5, columnspacing=0.8, handletextpad=0.2)
    ax.set_title("The AD reduction persists under every available adjustment", pad=4)
    panel_letter(ax, "d", x=-0.32)

    ax = fig.add_subplot(gs[1, 2])
    SS = pd.DataFrame(sens)
    x = np.arange(len(NVU5)); w = 0.38
    for k, (mname, colr) in enumerate([("+ age + PMI", "#2A6F97"), ("excluding acute-hypoxic donors", "#9CB9D1")]):
        q = SS[SS.model == mname].set_index("cell_type").loc[NVU5]
        ax.bar(x + (k - 0.5) * w, 100 * q.retained_FDR05 / q.primary_robust, w, color=colr, label=mname)
    ax.set_xticks(x); ax.set_xticklabels([{"Oligodendrocyte": "Oligodendro."}.get(c, c) for c in NVU5], rotation=40, ha="right")
    ax.set_ylim(0, 115); ax.set_ylabel("Primary robust DEGs still\nat FDR < 0.05 (%)")
    ax.legend(fontsize=5, loc="upper right", ncol=1)
    same = SS.same_sign_pct.min()
    ax.set_title(f"Hippocampal DE re-fitted\n(sign unchanged for {same:.0f}% of genes)", pad=4, fontsize=6)
    panel_letter(ax, "e", x=-0.4)
    save_fig(fig, OUT / "FigS2_donor_covariates")


if __name__ == "__main__":
    main()
