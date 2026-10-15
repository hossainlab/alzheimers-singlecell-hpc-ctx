"""Step 8 - Figure 5: transparent, computed target prioritisation.

Candidate genes = robust (padj < 0.05, |log2FC| >= 0.5, not artefact-flagged) in >= 1 hippocampal
NVU cell type (discovery set, 17 independent donors). Ranking is lexicographic, no hand-tuned weights:
  1. breadth   : number of hippocampal NVU cell types in which the gene is robust
  2. ctx_support: number of those cell types in which the cortex (same donors) effect has the same
                  sign and padj < 0.1
  3. evidence  : -log10 of the smallest hippocampal padj
Annotation (all retrieved programmatically, cached in data/annotation/):
  * UniProtKB/Swiss-Prot (human, reviewed): accession, protein name, subcellular location, TM/signal
  * mouse orthologue by gene symbol in reviewed mouse UniProtKB; human-mouse identity computed by
    global Needleman-Wunsch alignment (BLOSUM62, gap open -10, extend -0.5); identity = identical
    positions / human sequence length
  * STRING v12 functional network (combined score >= 0.7) among top candidates plus AD GWAS genes

Outputs: figures/Fig5_Target_Prioritization/Fig5_targets.{pdf,png}
         tables/target_prioritization.csv, tables/target_annotation_top.csv, tables/STRING_edges.csv
"""
import sys, json, time
from pathlib import Path

import numpy as np
import pandas as pd
import requests
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.colors import TwoSlopeNorm
from Bio import Align
from Bio.Align import substitution_matrices

sys.path.insert(0, str(Path(__file__).parent))
from config import *  # noqa

OUT = FIGS / "Fig5_Target_Prioritization"
ANN = DATA / "annotation"; ANN.mkdir(exist_ok=True)
AD_GENES = ["APOE", "APP", "PSEN1", "TREM2", "CLU", "BIN1", "PICALM", "SORL1", "ABCA7", "CD2AP", "PTK2B",
            "ADAM10", "MS4A6A", "INPP5D", "PLCG2", "CR1", "TYROBP", "MAPT"]
SHORT = {"Microglia": "Micro", "Endothelial": "Endo", "Pericyte": "Peri", "Astrocyte": "Astro", "Oligodendrocyte": "Oligo"}
TOPN = 20


def uniprot(gene, org):
    f = ANN / f"uniprot_{org}_{gene}.json"
    if f.exists():
        return json.loads(f.read_text())
    q = f"(gene_exact:{gene}) AND (organism_id:{org}) AND (reviewed:true)"
    fields = "accession,gene_primary,protein_name,cc_subcellular_location,ft_transmem,ft_signal,sequence,length"
    r = requests.get("https://rest.uniprot.org/uniprotkb/search", params=dict(query=q, fields=fields, format="json", size=5), timeout=60)
    r.raise_for_status()
    res = r.json().get("results", [])
    res = [x for x in res if x.get("genes") and x["genes"][0].get("geneName", {}).get("value", "").upper() == gene.upper()] or res
    out = res[0] if res else {}
    f.write_text(json.dumps(out)); time.sleep(0.2)
    return out


def parse_uniprot(e):
    if not e:
        return dict(accession=None, protein=None, location=None, n_tm=0, signal=False, seq=None)
    pn = e.get("proteinDescription", {}).get("recommendedName", {}).get("fullName", {}).get("value")
    locs = []
    for c in e.get("comments", []):
        if c.get("commentType") == "SUBCELLULAR LOCATION":
            for l in c.get("subcellularLocations", []):
                v = l.get("location", {}).get("value")
                if v and v not in locs:
                    locs.append(v)
    feats = e.get("features", [])
    return dict(accession=e.get("primaryAccession"), protein=pn, location="; ".join(locs[:4]) or None,
                n_tm=sum(f["type"] == "Transmembrane" for f in feats), signal=any(f["type"] == "Signal" for f in feats),
                seq=e.get("sequence", {}).get("value"))


def identity(a, b):
    al = Align.PairwiseAligner(mode="global", open_gap_score=-10, extend_gap_score=-0.5,
                               substitution_matrix=substitution_matrices.load("BLOSUM62"))
    aln = al.align(a, b)[0]
    ident = sum(1 for (s1, e1), (s2, e2) in zip(*aln.aligned) for i, j in zip(range(s1, e1), range(s2, e2)) if a[i] == b[j])
    return 100 * ident / len(a)


def access_class(loc, n_tm, signal):
    l = (loc or "").lower()
    if "secreted" in l or "extracellular" in l:
        return "Secreted / extracellular"
    if n_tm > 0 or "cell membrane" in l or "plasma membrane" in l:
        return "Cell surface / membrane"
    if l:
        return "Intracellular"
    return "Unknown"


def string_network(genes):
    """STRING v12: map symbols to STRING IDs (exact preferred-name matches only), then fetch edges."""
    base = "https://string-db.org/api/json"
    r = requests.post(f"{base}/get_string_ids", data=dict(identifiers="\r".join(genes), species=9606, limit=1,
                                                         echo_query=1, caller_identity="ad_nvu_pipeline"), timeout=120)
    r.raise_for_status()
    ids = pd.DataFrame(r.json())
    ids = ids[ids.queryItem.str.upper() == ids.preferredName.str.upper()]
    (ANN / "string_ids.json").write_text(ids.to_json(orient="records"))
    r = requests.post(f"{base}/network", data=dict(identifiers="\r".join(ids.stringId), species=9606,
                                                    required_score=700, caller_identity="ad_nvu_pipeline"), timeout=120)
    r.raise_for_status()
    (ANN / "string_network.json").write_text(r.text)
    E = pd.DataFrame(r.json())
    print("STRING mapped", len(ids), "of", len(genes), "| edges", len(E))
    return E


def main():
    set_style()
    de = pd.read_csv(TABLES / "pseudobulk_DE_all.csv.gz")
    H = de[de.region == "Hippocampus"]; C = de[de.region == "Cortex"].set_index(["cell_type", "gene"])
    rob = H[H.robust]
    rows = []
    for g, d in rob.groupby("gene"):
        cts = d.cell_type.tolist()
        sup = 0
        for r in d.itertuples():
            if (r.cell_type, g) in C.index:
                c = C.loc[(r.cell_type, g)]
                sup += int(np.sign(c.log2FoldChange) == np.sign(r.log2FoldChange) and c.padj < 0.1)
        rows.append(dict(gene=g, breadth_hpc=len(cts), celltypes=";".join(cts), ctx_support=sup,
                         direction="up" if d.log2FoldChange.median() > 0 else "down",
                         consistent_direction=bool(np.sign(d.log2FoldChange).nunique() == 1),
                         median_log2FC=d.log2FoldChange.median(), max_abs_log2FC=d.log2FoldChange.abs().max(),
                         min_padj=d.padj.min(), evidence=-np.log10(d.padj.min())))
    T = pd.DataFrame(rows)
    T = T[T.consistent_direction].sort_values(["breadth_hpc", "ctx_support", "evidence"], ascending=False).reset_index(drop=True)
    T["rank"] = np.arange(1, len(T) + 1)

    # annotation of the top candidates
    top = T.head(TOPN).copy()
    ann = []
    for g in top.gene:
        h = parse_uniprot(uniprot(g, 9606)); m = parse_uniprot(uniprot(g, 10090))
        idn = identity(h["seq"], m["seq"]) if h["seq"] and m["seq"] else np.nan
        ann.append(dict(gene=g, uniprot_human=h["accession"], protein=h["protein"], subcellular_location=h["location"],
                        n_transmembrane=h["n_tm"], signal_peptide=h["signal"], accessibility=access_class(h["location"], h["n_tm"], h["signal"]),
                        uniprot_mouse=m["accession"], human_mouse_identity_pct=round(idn, 1) if idn == idn else np.nan))
    A = pd.DataFrame(ann)
    top = top.merge(A, on="gene")
    T = T.merge(A, on="gene", how="left")
    T.to_csv(TABLES / "target_prioritization.csv", index=False)
    top.to_csv(TABLES / "target_annotation_top.csv", index=False)

    net_genes = list(dict.fromkeys(top.gene.tolist() + AD_GENES))
    E = string_network(net_genes)
    if len(E):
        E = E[["preferredName_A", "preferredName_B", "score"]].drop_duplicates()
        E = E[E.preferredName_A.isin(net_genes) | E.preferredName_B.isin(net_genes)]
    E.to_csv(TABLES / "STRING_edges.csv", index=False)

    # -------------------------------------------------- figure
    fig = plt.figure(figsize=(7.2, 7.8))
    gs = gridspec.GridSpec(2, 6, figure=fig, height_ratios=[1.35, 1], hspace=0.42, wspace=1.2)

    # a: evidence matrix
    ax = fig.add_subplot(gs[0, 0:4])
    cols = [(r_, ct) for r_ in ["Hippocampus", "Cortex"] for ct in NVU5]
    M = np.full((len(top), len(cols)), np.nan); P = np.full_like(M, np.nan)
    for j, (r_, ct) in enumerate(cols):
        d = de[(de.region == r_) & (de.cell_type == ct)].set_index("gene").reindex(top.gene)
        M[:, j] = d.log2FoldChange.values; P[:, j] = d.padj.values
    im = ax.imshow(M, aspect="auto", cmap="RdBu_r", norm=TwoSlopeNorm(0, -2.5, 2.5), interpolation="nearest")
    for i, j in zip(*np.where(P < 0.05)):
        ax.text(j, i, "•", ha="center", va="center", fontsize=6)
    for i, j in zip(*np.where(np.isnan(M))):
        ax.text(j, i, "n.t.", ha="center", va="center", fontsize=4, color="0.5")
    ax.set_yticks(range(len(top))); ax.set_yticklabels([f"{r} {g}" for r, g in zip(top["rank"], top.gene)], style="italic")
    ax.set_xticks(range(len(cols))); ax.set_xticklabels([f"{'Hpc' if r=='Hippocampus' else 'Ctx'} {SHORT[c]}" for r, c in cols], rotation=90)
    ax.axvline(4.5, color="k", lw=0.8)
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02); cb.set_label("log$_2$FC (AD/Control)", fontsize=5); cb.ax.tick_params(labelsize=5)
    ax.set_title(f"Top {TOPN} candidates ranked by hippocampal breadth, cortex support, evidence (• FDR < 0.05)")
    panel_letter(ax, "a", x=-0.18)

    # b: annotation columns
    ax = fig.add_subplot(gs[0, 4:6]); ax.axis("off")
    acc_col = {"Cell surface / membrane": "#1B9E77", "Secreted / extracellular": "#7570B3", "Intracellular": "#999999", "Unknown": "#DDDDDD"}
    ax.set_ylim(len(top) - 0.5, -1.5); ax.set_xlim(0, 3)
    ax.text(0.15, -1.1, "Breadth\n(Hpc)", fontsize=5, ha="center"); ax.text(1.0, -1.1, "Human–mouse\nidentity (%)", fontsize=5, ha="center")
    ax.text(2.2, -1.1, "Accessibility", fontsize=5, ha="center")
    for i, r in enumerate(top.itertuples()):
        ax.text(0.15, i, str(r.breadth_hpc), ha="center", va="center", fontsize=5.5)
        ax.text(1.0, i, "n/a" if pd.isna(r.human_mouse_identity_pct) else f"{r.human_mouse_identity_pct:.0f}", ha="center", va="center", fontsize=5.5)
        ax.add_patch(plt.Rectangle((1.6, i - 0.35), 0.25, 0.7, color=acc_col[r.accessibility]))
        ax.text(1.95, i, r.accessibility.split(" /")[0], va="center", fontsize=4.8)
    ax.set_title("Computed annotation", loc="left", y=1.04)
    panel_letter(ax, "b", x=-0.12, y=1.04)

    # c: per-donor expression of top 4 in their strongest hippocampal cell type
    smeta = pd.read_csv(TABLES / "sample_metadata.csv").set_index("sample")
    gsc = gs[1, 0:3].subgridspec(1, 4, wspace=0.6)
    for k, r in enumerate(top.head(4).itertuples()):
        ax = fig.add_subplot(gsc[k])
        d = H[(H.gene == r.gene) & H.robust].sort_values("padj").iloc[0]
        ct = d.cell_type
        for i, r_ in enumerate(["Hippocampus", "Cortex"]):
            Cn = pd.read_csv(PROC / "pseudobulk" / f"{r_}_{ct}_counts.csv.gz", index_col=0)
            nn = pd.read_csv(PROC / "pseudobulk" / f"{r_}_{ct}_ncells.csv", index_col=0).iloc[:, 0]
            Cn = Cn.loc[nn.index[nn >= 30]]
            v = np.log2(Cn[r.gene] / Cn.sum(1) * 1e6 + 1)
            for j, c_ in enumerate(["Control", "AD"]):
                vv = v[smeta.loc[v.index, "condition"] == c_].values
                x = i * 2.4 + j
                ax.scatter(x + np.random.default_rng(k + i + j).uniform(-0.15, 0.15, len(vv)), vv, s=6, color=COND_COLORS[c_], edgecolor="k", lw=0.25, zorder=3)
                ax.hlines(np.median(vv), x - 0.3, x + 0.3, color="k", lw=0.8)
        ax.set_xticks([0.5, 2.9]); ax.set_xticklabels(["Hpc", "Ctx"])
        ax.set_title(f"$\\it{{{r.gene}}}$\n{SHORT[ct]}", fontsize=6, loc="center")
        if k == 0:
            ax.set_ylabel("log$_2$(CPM + 1), per donor sample")
            panel_letter(ax, "c", x=-0.9, y=1.12)

    # d: STRING network
    ax = fig.add_subplot(gs[1, 3:6]); ax.axis("off")
    if len(E):
        import igraph as ig
        nodes = sorted(set(E.preferredName_A) | set(E.preferredName_B))
        g = ig.Graph(); g.add_vertices(nodes)
        g.add_edges([(a, b) for a, b in zip(E.preferredName_A, E.preferredName_B)])
        g.es["w"] = E.score.tolist()
        lay = np.array(g.layout_fruchterman_reingold(weights="w", niter=2000, seed=None).coords)
        for e_ in g.es:
            a, b = e_.tuple
            ax.plot(lay[[a, b], 0], lay[[a, b], 1], color="0.7", lw=0.4 + 1.2 * (e_["w"] - 0.7) / 0.3, zorder=1)
        cand = set(top.gene)
        for i, nm in enumerate(nodes):
            isc = nm in cand
            col = (COND_COLORS["AD"] if top.set_index("gene").loc[nm, "direction"] == "up" else COND_COLORS["Control"]) if isc else "#F2C14E"
            ax.scatter(*lay[i], s=40 if isc else 28, color=col, edgecolor="k", lw=0.4, zorder=3, marker="o" if isc else "s")
            ax.text(lay[i, 0], lay[i, 1] + 0.035 * np.ptp(lay[:, 1]), nm, fontsize=4.8, ha="center", style="italic", zorder=4)
        unm = sorted(set(top.gene) - set(json.loads((ANN / "string_ids.json").read_text()) and pd.read_json(ANN / "string_ids.json").queryItem))
        ax.set_title(f"STRING v12 (score ≥ 0.7): {len(cand & set(nodes))} of {TOPN} candidates connect"
                     + (f"\n(not in STRING: {', '.join(unm)})" if unm else ""), loc="left", fontsize=6)
        from matplotlib.lines import Line2D
        ax.legend(handles=[Line2D([], [], marker="o", ls="", mfc=COND_COLORS["AD"], mec="k", ms=4, label="Candidate, up in AD"),
                           Line2D([], [], marker="o", ls="", mfc=COND_COLORS["Control"], mec="k", ms=4, label="Candidate, down in AD"),
                           Line2D([], [], marker="s", ls="", mfc="#F2C14E", mec="k", ms=4, label="AD GWAS / reference gene")],
                  loc="upper left", fontsize=5, bbox_to_anchor=(0, -0.02), ncol=3, columnspacing=0.8)
    else:
        ax.text(0.5, 0.5, "No STRING edges at score ≥ 0.7", ha="center", transform=ax.transAxes)
    panel_letter(ax, "d", x=-0.06, y=1.04)

    save_fig(fig, OUT / "Fig5_targets")
    print(top[["rank", "gene", "breadth_hpc", "ctx_support", "direction", "min_padj", "accessibility", "human_mouse_identity_pct"]].to_string())


if __name__ == "__main__":
    main()
