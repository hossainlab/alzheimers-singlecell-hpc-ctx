# Donor-level single-nucleus analysis of the neurovascular unit in Alzheimer's disease (GSE163577)

Re-analysis of the VINE-seq human brain vascular atlas (Yang *et al.*, *Nature* 603:885, 2022; GEO
**GSE163577**). The data are hippocampus from 9 AD and 8 control donors, plus superior frontal cortex from
**8 of the same donors** (4 AD, 4 control). The pipeline is written in Python (scanpy / Harmony / pydeseq2 /
LIANA resource / gseapy). It runs end-to-end from the 25 GEO filtered matrices, with no down-sampling and no
hand-entered results.

> The earlier R/Seurat/CellChat pipeline, together with its tables and figures, has been moved to
> [`legacy/`](legacy/README.md). It is kept only for provenance. That file lists each problem found in it and how
> it was fixed here: per-nucleus pseudoreplicated DE, genotype- and ambient-driven "targets", same-donor cortex
> presented as validation, hard-coded Fig. 5 validation values, and correlation networks labelled as STRING.

---

## Design principles

| Issue | How it is handled |
|---|---|
| Unit of replication | Each **sample** (one donor × one region). All AD-vs-control statistics are donor-level: pseudobulk DESeq2, Mann–Whitney tests on per-sample proportions, and per-sample OLS for ligand–receptor scores. |
| VINE-seq enrichment yield | The vascular fraction is far lower in AD samples (Fig. 1f). Its logit is a covariate in every DE model. Cell-type proportions are shown per donor and are described as yield-dependent. |
| Ambient RNA / heterotypic doublets | Doublets are called per sample with scrublet. Clusters co-expressing two lineages' key markers are removed (`tables/cluster_annotation.csv`). Each DE gene is flagged if it is ≥ 5× more abundant in another cell type **and** tracks that cell type's fraction across samples (ρ ≥ 0.5). |
| Donor genotype and sex | Sex is inferred from XIST/chrY. The 17q21.31 copy-number marker ARL17B is bimodal across donors (FigS1g). Genes in that locus and sex-chromosome genes are excluded from HVGs and flagged in DE. Sex is a covariate. |
| Same donors in both regions | Donor IDs are harmonised across regions; sex agrees for 8/8 pairs and genotype for 7/8. The hippocampus is the discovery set. Cortex is reported as **a second tissue from the same donors**, not as independent replication. |
| Annotation / validation | Every value is computed: UniProt localisation, human–mouse identity by Needleman–Wunsch alignment, STRING v12 edges (score ≥ 0.7) with exact-name ID mapping. Claims with no data behind them (5xFAD, ROSMAP, antibody catalogue) were removed. |

## Pipeline (`scripts/`)

| Step | Script | Main outputs |
|---|---|---|
| 1 | `01_qc_doublets.py` | Per-sample MAD QC (log UMI and log genes, 5 MAD; %MT 3 MAD capped at 5 %; ≥ 200 genes), scrublet (fixed 0.25 cut-off when the automatic threshold fails), sex and 17q21.31 inference → `data/processed/qc_counts.h5ad`, `tables/qc_*`, `tables/sample_metadata.csv`, `tables/donor_pairing_concordance.csv` |
| 2 | `02_integrate_annotate.py` | Joint Hpc + Ctx object, per-sample seurat_v3 HVGs, PCA, Harmony on sample, Leiden, rule-based cluster annotation, vascular yield → `data/processed/integrated.h5ad`, `tables/cluster_annotation.csv`, `tables/celltype_counts_per_sample.csv` |
| 3 | `03_fig1_global_landscape.py` | **Fig. 1**: all cells, both regions |
| 4 | `04_pseudobulk_de.py` | DESeq2 `~ sex + vascular_yield + condition` for 5 NVU cell types × 2 regions, artefact flags, sensitivity analysis without the yield term → `tables/pseudobulk_DE_*.csv` |
| 5 | `05_fig2_deep_dives.py` | **Fig. 2a–e**: one figure per NVU cell type (sub-clusters, donor-level proportions, volcanoes, donor heat map, pre-ranked GSEA) → `tables/Fig2_*`, `tables/GSEA_*` |
| 6 | `06_fig3_convergence.py` | **Fig. 3**: cross-cell-type convergence, Hpc vs Ctx concordance, test of the legacy IEG claim |
| 7 | `07_fig4_cell_communication.py` | **Fig. 4**: per-sample L–R scores (LIANA consensus resource + APP–SORL1), donor-level differential test |
| 8 | `08_fig5_target_prioritization.py` | **Fig. 5**: transparent ranking with UniProt/STRING annotation → `tables/target_prioritization.csv`, `tables/target_annotation_top.csv`, `tables/STRING_edges.csv` |
| 9 | `09_qc_figures.py` | **Fig. S1**: QC, retention, sex and genotype |
| 10 | `10_donor_covariates.py` | **Fig. S2**: donor clinical covariates from Yang *et al.* Supplementary Table 1 (age, PMI, Braak/CERAD, APOE, neuropathology notes), the hypoxia-module score, covariate-adjusted models and DE sensitivity → `tables/donor_clinical_metadata.csv`, `tables/covariate_balance.csv`, `tables/hypoxia_module_*.csv`, `tables/pseudobulk_DE_sensitivity_clinical.csv`, `tables/validation_power_estimates.csv` |

Run everything with `python scripts/run_all.py`, using the environment in `environment.yml`. Settings and the
sample/donor parser are in `scripts/config.py`. Gene-set libraries (MSigDB Hallmark 2020, GO-BP 2023, Reactome
2022 from Enrichr) are cached in `data/genesets/`. UniProt and STRING responses are cached in `data/annotation/`.

## Key results

**Data.** 319,449 nuclei were loaded and 269,956 passed QC and doublet removal (84.5 %). Of these, 12,931 were in
heterotypic or ambiguous clusters and were excluded, leaving **257,025 nuclei** in 11 classes
(`figures/Fig1_Global_Landscape`). Donors: 8 female, 9 male; sex is balanced between groups (Fisher P > 0.6).

**Composition.** AD samples contain far fewer vascular nuclei (hippocampus median 24 % vs 62 %, P = 0.021;
cortex 31 % vs 79 %, P = 0.029). Their parenchymal fraction (oligodendrocytes, astrocytes) is correspondingly
higher. This most likely reflects how well vessel isolation worked in each sample. It should not be read as
cell-type loss or gain in AD tissue.

**Differential expression (donor-level).**

| Region | Microglia | Endothelial | Pericyte | Astrocyte | Oligodendrocyte |
|---|---|---|---|---|---|
| Hippocampus (8 v 9), robust DEGs ↑/↓ | 4 / 25 | 29 / 89 | 150 / 210 | 9 / 33 | 15 / 50 |
| Cortex (4 v 4), robust DEGs ↑/↓ | 1 / 2 | 1 / 19 | 0 / 0 | 3 / 4 | 0 / 0 |

In cortex, vascular yield almost perfectly separates AD from control samples. Adjusting for it therefore
removes most of the power. Without the adjustment, 31–455 genes per cell type are called, but these cannot be
distinguished from composition and ambient effects (`tables/pseudobulk_DE_sensitivity_no_yield.csv`).

**Convergent program.** 47 genes are robust in ≥ 2 hippocampal NVU cell types, and almost all of them are
*lower* in AD. They form a hypoxia/HIF-1 and stress-response module: *ADAMTS9* (4 cell types), *VEGFA, ELL2,
BCL6, YBX3, HIF1A* (3 each), *HILPDA, SLC2A3, PDK1, BNIP3L, PDE10A*. Pre-ranked GSEA gives Hallmark *Hypoxia*
NES −2.1 to −2.7 in all five hippocampal cell types (FDR < 0.002), together with *TNF-α/NF-κB* and
*Inflammatory response*. The cortex of the same donors agrees in direction for 77–87 % of hippocampal robust
DEGs (binomial P from 8.8 × 10⁻³⁷ for pericytes to 5.5 × 10⁻⁴ for microglia; `tables/cross_region_concordance.csv`), although few cortical genes reach significance on their own. In pericytes, an
*ADAMTS9/ELL2/PCBP3*-high sub-state makes up a median 41 % of hippocampal pericytes in controls and 3 % in AD
(FDR = 0.044).

**Donor covariates (Fig. S2).** Yang *et al.* Supplementary Table 1 maps onto all 17 donors (sex agrees
17/17; the cortex subset agrees 17/17).
- *Balance.* AD and control donors are balanced for age (84 vs 78 y, P = 0.15), post-mortem interval
  (5.8 vs 5.0 h, P = 0.77), APOE ε4 dose and recorded acute hypoxic neuropathology (one donor per group).
- *Association with the score.* The hypoxia-module score is not associated with age, PMI or APOE ε4 in any
  cell type (|ρ| ≤ 0.23).
- *The AD effect persists.* In pericytes it is −1.48 module-score units (P = 0.0007) after adjusting for age and PMI, and
  −1.45 (P = 0.004) after excluding the hypoxic-pathology donors. It stays significant in all five cell
  types under every model.
- *Gene-level DE.* Re-fitting with age and PMI keeps 43–86 % of robust DEGs at FDR < 0.05, with no sign
  changes.

> **Interpretation caveat.** Hypoxia-inducible transcripts in post-mortem brain also depend on agonal state
> and tissue pH. Supplementary Table 1 does not report agonal score, brain pH or RIN. Recorded acute hypoxic
> neuropathology (2 donors) is only a coarse proxy. A terminal-hypoxia contribution is therefore **reduced
> but not excluded**. The remaining steps are obtaining those variables from the tissue source (the
> Stanford/VA/NIA Aging Clinical Research Center) and replicating in an independent cohort.

**Claims from the legacy analysis that are not supported:**
- ARL17B and the lncRNA "targets" follow donor 17q21.31 genotype.
- No convergent immediate-early-gene program (FOS/JUN/EGR1/DUSP1/GADD45B) is significant in more than one
  cell type (Fig. 3f).
- APP→SORL1 signalling into microglia does not differ (FDR 0.81 in hippocampus, 1.0 in cortex; Fig. 4d).
- No ligand–receptor interaction changes at donor-level FDR < 0.05 (Fig. 4).

**Candidates.** The ranking uses, in order, the number of hippocampal cell types where the gene is robust,
cortex support, and evidence (`tables/target_prioritization.csv`). The top hits are *ADAMTS9* (secreted;
human–mouse identity 89 %), *ELL2, VEGFA, BCL6, YBX3, HIF1A*. STRING links *HIF1A–STAT3–BCL6–SLC2A3* and
*COL4A2–LAMB1*; none of the top 20 connects directly to AD GWAS genes at score ≥ 0.7. VEGFA cannot be resolved
by symbol in STRING v12, and this is stated in Fig. 5d.

## Figures

All figures are vector PDF plus 300-dpi PNG, use CVD-safe colours, and have sentence-style panel titles. Every
summary is shown with per-donor points.

- `figures/Fig1_Global_Landscape/Fig1_global_landscape`: integrated UMAP of all 257,025 nuclei from both regions
  (cell type, region, condition), canonical-marker dot plot, per-sample composition, vascular yield.
- `figures/Fig2_Cell_Type_Deep_Dives/Fig2a_Microglia` … `Fig2e_Oligodendrocyte`: one deep-dive figure per NVU
  cell type, in a fixed four-row order. Each row reads left to right: (a, b) subcluster, region and
  condition UMAPs; (c, d) markers and per-sample proportions; (e, f) hippocampus and cortex volcano plots;
  (g, h) per-sample heatmap of robust DEGs and GSEA.
- `figures/Fig3_Cross_Cell_Convergence/Fig3_convergence`
- `figures/Fig4_Intercellular_Communication/Fig4_communication`
- `figures/Fig5_Target_Prioritization/Fig5_targets`
- `figures/Quality_Control/FigS1_quality_control`
- `figures/Quality_Control/FigS2_donor_covariates`

See [`figures/README.md`](figures/README.md) for panel-by-panel legends.

**Single panels.** `single_figures/` holds every lettered panel as a separate PNG (300 dpi) and PDF. Panel
letters and panel titles are removed; axes, legends, colour bars and facet labels are kept. They are written
automatically by `save_fig()` in `scripts/config.py` whenever a figure script runs. File names:
- `Fig1a` … `Fig1f` for Fig. 1, and likewise for Fig. 3, Fig. 4, Fig. 5, Fig. S1 and Fig. S2.
- `Fig2<cell type letter>_<panel>` for Fig. 2. For example, `Fig2c_e` is pericyte panel e, the hippocampal
  volcano plot.

## Limitations

1. Only filtered matrices are on GEO. SoupX/CellBender need raw droplets, so ambient RNA is handled through the
   yield covariate and the contamination flags, not by count correction.
2. Age, PMI, APOE, Braak/CERAD and neuropathology notes come from Yang *et al.* Supplementary Table 1 and
   are tested in `10_donor_covariates.py`. Agonal-factor score, brain pH and RIN are not reported anywhere and
   remain unmodelled. Braak stage is missing for 4 of 9 AD donors.
3. Cortex (4 v 4) is underpowered, and yield is nearly confounded with condition there. The cortex comes from
   the same donors as the hippocampus, so it cannot serve as independent replication.
4. VINE-seq enriches vessels. Parenchymal populations (neurons, OPCs) are small and are not analysed for DE.
5. Mouse orthologues are matched by gene symbol in reviewed UniProtKB. Primate-specific genes (e.g. *ADGRE2*)
   have no orthologue and are shown as n/a.

## License and attribution
Raw data: GEO GSE163577 (Yang, Vest, …, Wyss-Coray; *Nature* 2022). Code: see `LICENSE`.
