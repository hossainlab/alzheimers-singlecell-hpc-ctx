"""Shared configuration for the GSE163577 (VINE-seq) re-analysis pipeline.

All scripts are run from the repository root, e.g. ``python scripts/02_qc_doublets.py``.
"""
from __future__ import annotations

import re
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RAW_DIR = ROOT / "data" / "GSE163577_RAW"
DATA = ROOT / "data"
PROC = DATA / "processed"
GENESETS = DATA / "genesets"
TABLES = ROOT / "tables"
FIGS = ROOT / "figures"
LOGS = ROOT / "logs"
for _p in (PROC, GENESETS, TABLES, FIGS, LOGS):
    _p.mkdir(parents=True, exist_ok=True)

SEED = 0

# ----------------------------------------------------------------------------- QC
QC = dict(
    min_genes=200,          # hard floor on detected genes per nucleus
    nmads_counts=5,         # per-sample MAD outlier rule on log total counts (both tails)
    nmads_genes=5,          # per-sample MAD outlier rule on log genes (both tails)
    nmads_mt=3,             # per-sample MAD rule on % mitochondrial (upper tail)
    max_mt_pct=5.0,         # absolute ceiling on % mitochondrial (nuclei)
    min_mt_ceiling=1.0,     # MAD-derived MT threshold is never set below this
    min_cells_per_gene=10,  # gene retained if detected in >= this many nuclei overall
)

# ------------------------------------------------------------------- cell types
NVU5 = ["Microglia", "Endothelial", "Pericyte", "Astrocyte", "Oligodendrocyte"]

MARKERS = {
    "Endothelial":       ["CLDN5", "FLT1", "VWF", "PECAM1", "ABCB1"],
    "Pericyte":          ["PDGFRB", "RGS5", "ABCC9", "KCNJ8"],
    "Smooth muscle":     ["ACTA2", "MYH11", "TAGLN"],
    "Fibroblast":        ["COL1A2", "DCN", "LUM", "CEMIP"],
    "Astrocyte":         ["AQP4", "GFAP", "SLC1A2", "GJA1"],
    "Microglia":         ["CSF1R", "P2RY12", "CX3CR1", "C3"],
    "Macrophage":        ["MRC1", "F13A1", "LYVE1", "CD163"],
    "T cell":            ["CD3E", "CD2", "IL7R", "SKAP1"],
    "Oligodendrocyte":   ["MBP", "MOG", "PLP1", "ST18"],
    "OPC":               ["PDGFRA", "VCAN", "OLIG2"],
    "Excitatory neuron": ["SLC17A7", "SATB2", "CAMK2A", "RBFOX3"],
    "Inhibitory neuron": ["GAD1", "GAD2", "SST", "VIP"],
    "Ependymal":         ["FOXJ1", "TTR", "CFAP299"],
}

DISPLAY_MARKERS = {**{k: v[:3] for k, v in MARKERS.items()}, "Neuron": ["RBFOX3", "GAD1", "SLC17A7"],
                   "Microglia": ["P2RY12", "CX3CR1", "C3"], "Oligodendrocyte": ["MOG", "ST18", "PLP1"]}

CELLTYPE_ORDER = ["Endothelial", "Pericyte", "Smooth muscle", "Fibroblast", "Astrocyte",
                  "Microglia", "Macrophage", "T cell", "Oligodendrocyte", "OPC",
                  "Neuron", "Excitatory neuron", "Inhibitory neuron", "Ependymal"]

# Okabe-Ito based, CVD-safe; vascular = blue/green family, glia = warm, neurons = purple/grey
CT_COLORS = {
    "Endothelial": "#0072B2", "Pericyte": "#56B4E9", "Smooth muscle": "#2A9D8F",
    "Fibroblast": "#8AB17D", "Astrocyte": "#E69F00", "Microglia": "#D55E00",
    "Macrophage": "#B5651D", "T cell": "#7F3C8D", "Oligodendrocyte": "#CC79A7",
    "OPC": "#E8A1C6", "Neuron": "#6A4C93", "Excitatory neuron": "#6A4C93", "Inhibitory neuron": "#9C89B8",
    "Ependymal": "#8C8C8C", "Unassigned": "#D0D0D0",
}
COND_COLORS = {"Control": "#4C72B0", "AD": "#C44E52"}
REGION_COLORS = {"Hippocampus": "#3B7D6E", "Cortex": "#B07AA1"}
VASCULAR = ["Endothelial", "Pericyte", "Smooth muscle", "Fibroblast"]

# Genes in the 17q21.31 inversion / copy-number region; expression is genotype-driven.
LOCUS_17Q21 = ["ARL17A", "ARL17B", "LRRC37A", "LRRC37A2", "LRRC37A4P", "KANSL1", "KANSL1-AS1",
               "MAPT-AS1", "CRHR1-IT1", "LINC02210", "NSF", "MAPK8IP1P2", "CRHR1", "SPPL2C", "MAPT"]
SEX_GENES_PATTERN = r"^(XIST|TSIX|UTY|DDX3Y|USP9Y|KDM5D|NLGN4Y|RPS4Y1|EIF1AY|ZFY|TTTY\d+|TXLNGY|PRKY)$"


# ------------------------------------------------------------- sample metadata
def sample_table() -> pd.DataFrame:
    """Parse GEO file names into sample metadata, harmonising donor IDs across regions.

    Cortex file tokens use the letter 'O' for zero (e.g. ``O4_11AD``); after normalisation each
    cortex sample maps onto exactly one hippocampal sample from the same donor (validated in
    02_qc_doublets.py by inferred sex and 17q21.31 genotype).
    """
    rows = []
    for f in sorted(RAW_DIR.glob("GSM*_filtered_feature_bc_matrix.tar.gz")):
        gsm, rest = f.name.split("_", 1)
        tok = rest.replace("_filtered_feature_bc_matrix.tar.gz", "")
        region = "Cortex" if tok.endswith("_Ctx") else "Hippocampus"
        tok = tok.replace("_Ctx", "")
        if region == "Cortex":
            tok = tok.replace("O", "0")             # 'O4_O9C' -> '04_09C'
        rows.append(dict(gsm=gsm, token=tok, region=region, file=str(f)))
    df = pd.DataFrame(rows)

    def parse(tok: str, region: str):
        if tok.startswith("H2004"):
            return "H2004", "H2004", "AD"
        if region == "Hippocampus":
            parts = tok.split("_")
            batch, num = parts[0], parts[1]
            cond = {"C": "Control", "AD": "AD"}.get(parts[2] if len(parts) > 2 else "", None)
            return batch, f"{batch}_{num}", cond
        # cortex, e.g. 01_10C / 03_03C / 04_11AD / 05_10AD after O->0
        m = re.match(r"^(\d\d)_(\d\d)(C|AD)$", tok)
        if not m:
            raise ValueError(f"unparsed cortex token {tok}")
        return m.group(1), f"{m.group(1)}_{m.group(2)}", {"C": "Control", "AD": "AD"}[m.group(3)]

    parsed = [parse(t, r) for t, r in zip(df.token, df.region)]
    df["batch"], df["donor"], df["cond_from_name"] = zip(*parsed)
    # Condition from the GEO sample title (authoritative; GSM4982089 '02_04' has no suffix)
    geo = pd.read_csv(ROOT / "legacy" / "tables" / "GSE163577_sample_metadata.csv")
    df = df.merge(geo[["gsm", "title", "condition"]], on="gsm", how="left")
    bad = df.cond_from_name.notna() & (df.cond_from_name != df.condition)
    assert not bad.any(), df[bad]
    df["sample"] = df.region.map({"Hippocampus": "Hpc", "Cortex": "Ctx"}) + "_" + df.donor + "_" + df.condition
    return df[["gsm", "sample", "donor", "region", "condition", "batch", "title", "file"]]


# ---------------------------------------------------------------- plotting
def set_style():
    import matplotlib as mpl
    mpl.rcParams.update({
        "font.family": "sans-serif", "font.sans-serif": ["Arial", "Helvetica", "DejaVu Sans"],
        "font.size": 7, "axes.titlesize": 7, "axes.labelsize": 7, "xtick.labelsize": 6,
        "ytick.labelsize": 6, "legend.fontsize": 6, "legend.title_fontsize": 6,
        "axes.linewidth": 0.6, "xtick.major.width": 0.6, "ytick.major.width": 0.6,
        "xtick.major.size": 2.5, "ytick.major.size": 2.5, "axes.spines.top": False,
        "axes.spines.right": False, "legend.frameon": False, "pdf.fonttype": 42, "ps.fonttype": 42,
        "savefig.dpi": 300, "figure.dpi": 150, "axes.titlelocation": "left",
        "axes.titleweight": "normal", "lines.linewidth": 0.8, "patch.linewidth": 0.5,
    })


def panel_letter(ax, letter, x=-0.12, y=1.06):
    t = ax.text(x, y, letter, transform=ax.transAxes, fontsize=9, fontweight="bold", va="bottom", ha="left")
    t.set_gid("panel_letter")
    ax.figure.__dict__.setdefault("_panel_axes", []).append((letter, ax))
    return t


SINGLE_DIR = ROOT / "single_figures"


def _single_name(stem: str, letter: str) -> str:
    """Fig1_global_landscape + a -> Fig1a ; Fig2c_Pericyte + b -> Fig2c_b ; FigS2_x + d -> FigS2d."""
    m = re.match(r"^(Fig(?:S)?\d+)([a-z]?)(?:_.*)?$", stem)
    base, sub = (m.group(1), m.group(2)) if m else (stem, "")
    return f"{base}{sub}_{letter}" if sub else f"{base}{letter}"


def export_single_panels(fig, stem: str, out_dir: Path = None):
    """Export every lettered panel of a composite figure as its own PNG + PDF.

    Panel letters and panel/figure titles are removed; axis labels, ticks, legends, colour bars and facet
    labels of multi-axes panels (e.g. one gene name per sub-axis) are kept.
    A panel = the lettered axes plus every axes created after it and before the next lettered axes
    (colour bars, marginal axes, sub-grids). Figure-level legends/texts go to the nearest panel.
    The composite figure is restored afterwards.
    """
    from matplotlib.transforms import Bbox
    panels = fig.__dict__.get("_panel_axes", [])
    if not panels:
        return []
    out_dir = Path(out_dir or SINGLE_DIR); out_dir.mkdir(parents=True, exist_ok=True)
    order = list(fig.axes)
    starts = sorted(((order.index(ax), L) for L, ax in panels if ax in order))
    groups = {}
    for k, (i0, L) in enumerate(starts):
        i1 = starts[k + 1][0] if k + 1 < len(starts) else len(order)
        groups[L] = [a for a in order[i0:i1] if a.get_visible()]
    fig.canvas.draw()
    r = fig.canvas.get_renderer()
    hidden = []                                      # (artist, previous visibility)

    def hide(a):
        hidden.append((a, a.get_visible())); a.set_visible(False)
    # Facet labels are kept: when >= 2 axes inside one panel carry titles (e.g. one gene per sub-axis),
    # those titles identify the sub-plots rather than title the panel.
    facet = set()
    for g in groups.values():
        titled = [a for a in g if any(t.get_text() for t in (a.title, a._left_title, a._right_title))]
        if len(titled) >= 2:
            facet.update(titled)
    for ax in order:
        for t in (ax.title, getattr(ax, "_left_title", None), getattr(ax, "_right_title", None)):
            if ax not in facet and t is not None and t.get_text():
                hide(t)
        for t in ax.texts:
            if t.get_gid() == "panel_letter":
                hide(t)
    if getattr(fig, "_suptitle", None) is not None:
        hide(fig._suptitle)
    for t in fig.texts:
        if t.get_gid() == "panel_letter":
            hide(t)
    # assign figure-level legends / texts to the nearest panel
    def centre(bb):
        return ((bb.x0 + bb.x1) / 2, (bb.y0 + bb.y1) / 2)
    gbb = {L: Bbox.union([a.get_tightbbox(r) for a in g]) for L, g in groups.items() if g}
    extras = {L: [] for L in gbb}
    for art in list(fig.legends) + [t for t in fig.texts if t.get_visible() and t.get_text()]:
        bb = art.get_window_extent(r); cx, cy = centre(bb)
        L = min(gbb, key=lambda k: (centre(gbb[k])[0] - cx) ** 2 + (centre(gbb[k])[1] - cy) ** 2)
        extras[L].append(art)
    written = []
    all_extra = [a for v in extras.values() for a in v]
    for L, g in groups.items():
        if not g:
            continue
        keep = set(g) | set(extras[L])
        vis = [(a, a.get_visible()) for a in order + all_extra]
        for a, _ in vis:
            a.set_visible(a in keep and _)
        fig.canvas.draw()
        bbs = [a.get_tightbbox(r) for a in g] + [e.get_window_extent(r) for e in extras[L]]
        bbs += [a.get_legend().get_window_extent(r) for a in g if a.get_legend() is not None]
        bb = Bbox.union(bbs).transformed(fig.dpi_scale_trans.inverted()).padded(0.04)
        name = _single_name(stem, L)
        fig.savefig(out_dir / f"{name}.pdf", bbox_inches=bb)
        fig.savefig(out_dir / f"{name}.png", dpi=300, bbox_inches=bb)
        written.append(name)
        for a, v in vis:
            a.set_visible(v)
    for a, v in reversed(hidden):
        a.set_visible(v)
    return written


def save_fig(fig, path_no_ext: Path, single_panels: bool = True):
    path_no_ext = Path(path_no_ext)
    path_no_ext.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(str(path_no_ext) + ".pdf", bbox_inches="tight")
    fig.savefig(str(path_no_ext) + ".png", dpi=300, bbox_inches="tight")
    if single_panels:
        export_single_panels(fig, path_no_ext.name)


def harmony(ad, key="sample", basis="X_pca", out="X_pca_harmony", max_iter=20):
    """Harmony integration (harmonypy) robust to the Z_corr orientation change in harmonypy >= 0.2."""
    import harmonypy
    ho = harmonypy.run_harmony(ad.obsm[basis], ad.obs, [key], max_iter_harmony=max_iter,
                               random_state=SEED, verbose=False)
    Z = np.asarray(ho.Z_corr)
    if Z.shape[0] != ad.n_obs:
        Z = Z.T
    ad.obsm[out] = Z


def bh(p):
    """Benjamini-Hochberg FDR (NaN-safe)."""
    p = np.asarray(p, float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    from statsmodels.stats.multitest import multipletests
    if ok.any():
        out[ok] = multipletests(p[ok], method="fdr_bh")[1]
    return out
