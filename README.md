# Single-Cell Transcriptomic Identification of Alzheimer’s Disease-Associated Molecular Targets and CNS–Peripheral Immune Signatures for Experimental Validation and Therapeutic Intervention

[![R](https://img.shields.io/badge/R-4.6.1-blue.svg)](https://www.r-project.org/)
[![Seurat](https://img.shields.io/badge/Seurat-v5.0-green.svg)](https://satijalab.org/seurat/)
[![CellChat](https://img.shields.io/badge/CellChat-v2.2-orange.svg)](https://github.com/jinworks/CellChat)
[![Dataset](https://img.shields.io/badge/GEO-GSE163577-purple.svg)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163577)

---

## 🔬 Study Overview

This repository houses the complete, reproducible end-to-end single-nucleus RNA sequencing (snRNA-seq) analytical workflow for the study **"Single-Cell Transcriptomic Identification of Alzheimer’s Disease-Associated Molecular Targets and CNS–Peripheral Immune Signatures for Experimental Validation and Therapeutic Intervention"**.

The analysis interrogates **GSE163577**, comprising high-throughput single-nucleus transcriptomes from the human **Prefrontal Cortex (Ctx)** and **Hippocampus (Hpc)** of clinically and neuropathologically verified **Alzheimer's Disease (AD)** donors alongside age-matched cognitively intact **Control** individuals. The analytical pipeline centers on resolving disease pathogenesis, neurovascular unit (NVU) disruption, cross-cell molecular convergence, and intercellular communication rewiring across five essential glial and vascular cell types:
1. **Microglia** (Innate neuroimmune surveillance and phagocytosis)
2. **Endothelial Cells** (Blood–brain barrier core endothelium)
3. **Pericytes** (Mural vascular integrity and basement membrane regulation)
4. **Astrocytes** (Metabolic coupling and neurovascular maintenance)
5. **Oligodendrocytes** (Myelination and trophic support)

---

## 📋 The 8-Step Analytical Framework

The study rigorously adheres to an 8-step single-cell systems biology framework:

| Step | Analysis Phase | Description | Key Scripts | Key Outputs & Visualizations |
| :---: | :--- | :--- | :--- | :--- |
| **1** | **QC & Preprocessing** | Nuclear filtering, ambient RNA assessment, mitochondrial cutoff, doublet removal, count normalization | `scripts/01_*.R`<br>`scripts/02*.R` | `tables/GSE163577_sample_metadata.csv`<br>`data/*_qc.rds` |
| **2** | **Integration & Clustering** | Batch correction, Harmony/Seurat PCA, UMAP dimensionality reduction, canonical marker annotation | `scripts/03a_integrate_cortex.R`<br>`scripts/03b_integrate_hippocampus.R` | `data/seurat_*_annotated.rds`<br>`Cortex/Fig1_*`, `Hippocampus/Fig1_*` |
| **3** | **Differential Expression (Core)** | Cell-type-specific Wilcoxon/MAST DEG testing across AD vs Control (FDR < 0.05, \|log2FC\| ≥ 0.25) | `scripts/04_*.R`<br>`scripts/process_hippocampus_*.R` | `tables/DEGs_*.csv`<br>`Cortex/Fig2_*`, `Hippocampus/Fig2_*` |
| **4** | **Functional Enrichment** | Hypergeometric overrepresentation testing across KEGG Pathways and Gene Ontology (BP, MF, CC) | `scripts/04_*.R`<br>`scripts/05_*.R` | `Cortex/Fig2f`, `Fig2g`<br>`Hippocampus/Fig2f`, `Fig2g` |
| **5** | **Cell–Cell Communication** | CellChat v2.2 mass-action modeling of ligand–receptor signaling networks, information flow, and circos paths | `scripts/06_*.R`<br>`scripts/finish_hippocampus_*.R` | `Cortex/Fig4_*`, `Hippocampus/Fig4_*`<br>`tables/CellChat_interactions_*.csv` |
| **6** | **PPI Networks** | STRING v12 / igraph protein–protein interaction network modeling, Louvain clustering, hub identification | `scripts/04_*.R`<br>`scripts/05_*.R` | `Cortex/Fig2h`, `Fig3d`<br>`Hippocampus/Fig2h`, `Fig3d` |
| **7** | **Target Prioritization & Validation** | Dual-region concordance (Cortex vs Hippocampus), CTPS ranking, and in silico cross-species orthology | `scripts/05_*.R`<br>`scripts/07_*.R` | `tables/consensus_target_prioritization_full.csv`<br>`Cross_Region_Comparison/` |
| **8** | **Single Panel Publication Exports** | High-resolution publication-ready vector PDFs and 300-DPI PNGs structured in clean figure folders | `scripts/reorganize_option1.ps1` | `Cortex/`, `Hippocampus/`, `Cross_Region_Comparison/` |

---

## 📁 Repository Structure

```
alzheimers-singlecell-hpc-ctx/
│
├── data/                            # Raw matrices and processed RDS objects
│   ├── GSE163577_RAW/               # Raw 10x Genomics matrices (.tar.gz) from GEO
│   ├── seurat_cortex_annotated.rds  # Full annotated Prefrontal Cortex snRNA-seq object
│   ├── seurat_hpc_annotated.rds     # Full annotated Hippocampus snRNA-seq object
│   ├── fig2_results_cortex.rds      # Target cell subclusters, DEGs, and enrichment results
│   ├── cellchat_test_merged.rds     # Comparative CellChat (Control vs AD) merged cortex object
│   └── cellchat_hpc_merged.rds      # Comparative CellChat (Control vs AD) merged hippocampus object
│
├── scripts/                         # Numbered, reproducible snRNA-seq analysis pipeline
│   ├── 01_data_exploration_and_metadata.R      # Step 1: GEO matrix discovery and metadata extraction
│   ├── 02a_process_cortex.R                    # Step 2a: Cortex snRNA-seq QC filtering and merging
│   ├── 02b_process_hippocampus.R               # Step 2b: Hippocampus snRNA-seq QC filtering and merging
│   ├── 02c_qc_and_filtering_figures.R          # Step 2c: Empirical QC figure suite & reviewer defense
│   ├── 03a_integrate_cortex.R                  # Step 3a: Cortex Harmony integration, annotation & Fig 1
│   ├── 03b_integrate_hippocampus.R             # Step 3b: Hippocampus Harmony integration, annotation & Fig 1
│   ├── 04a_cortex_subclustering_and_degs.R     # Step 4a: Cortex 6-lineage deep dives, Wilcoxon DEGs & Fig 2
│   ├── 04b_hippocampus_subclustering_and_degs.R # Step 4b: Hippocampus 6-lineage deep dives, DEGs & Fig 2
│   ├── 05a_cortex_cross_cell_convergence.R     # Step 5a: Cortex cross-cell convergence, CTPS ranking & Fig 3
│   ├── 05b_hippocampus_cross_cell_convergence.R # Step 5b: Hippocampus cross-cell convergence & Fig 3
│   ├── 06a_cortex_cell_cell_communication.R    # Step 6a: Cortex CellChat v2.2 crosstalk modeling & Fig 4
│   ├── 06b_hippocampus_cell_cell_communication.R # Step 6b: Hippocampus CellChat crosstalk modeling & Fig 4
│   ├── 07a_cortex_target_validation.R          # Step 7a: Cortex in silico target validation & Fig 5
│   ├── 07b_hippocampus_target_validation.R     # Step 7b: Hippocampus in silico target validation & Fig 5
│   ├── 08_cross_region_comparison.R            # Step 8: Head-to-head Cortex vs Hippocampus comparison
│   ├── process_cortex_all_analyses.R           # Master runner: Prefrontal Cortex end-to-end
│   ├── process_hippocampus_all_analyses.R      # Master runner: Hippocampus end-to-end
│   ├── explore_samples.py                      # Python utility: Sample discovery and statistics
│   ├── parse_metadata.py                       # Python utility: NCBI GEO metadata extractor
│   └── show_top_targets.py                     # Python utility: CTPS ranking inspection
│
├── tables/                          # Quantified differential expression, enrichment, and communication tables
│   ├── GSE163577_sample_metadata.csv           # Complete sample-level clinical and technical metadata
│   ├── sample_qc_retention_rates.csv           # Donor-level initial vs filtered nuclei retention rates
│   ├── DEGs_Cortex_*.csv                       # Cortex differential expression statistics per cell type
│   ├── DEGs_Hippocampus_*.csv                  # Hippocampus differential expression statistics per cell type
│   ├── consensus_target_prioritization_full.csv # Ranked cross-cell convergent therapeutic targets
│   ├── top5_prioritized_therapeutic_targets.csv # Top 5 multi-factorial therapeutic target candidates
│   ├── top5_canonical_protein_targets.csv      # Top 5 canonical protein-coding therapeutic targets
│   ├── target_in_silico_validation_metrics.csv # Cross-species orthology, identity, and suitability index
│   ├── recommended_antibodies_and_compounds.csv # Pre-validated antibody clones, vendors, and catalog numbers
│   ├── CellChat_interactions_*.csv             # Validated intercellular ligand–receptor pairs
│   └── CellChat_interactions_Hippocampus_*.csv # Validated Hippocampus ligand–receptor pairs
│
├── figures/                         # Publication figures suites (Option 1 Clean Hierarchy)
│   ├── Cortex/                      # Prefrontal Cortex figure suite
│   │   ├── Fig1_Global_Landscape/                  # Fig1a to Fig1e (UMAPs, proportions, markers)
│   │   ├── Fig2_Cell_Type_Deep_Dives/              # Subfolders (Astrocyte, Microglia, etc.) with panels + DEGs
│   │   ├── Fig3_Cross_Cell_Convergence/            # Fig3a to Fig3g (Venn, UpSet, KEGG, violins, ranking)
│   │   ├── Fig4_Intercellular_Communication/       # Fig4a to Fig4e (CellChat networks, circle plots, CSVs)
│   │   ├── Fig5_Target_Validation/                 # Fig5a to Fig5c (Conservation, druggability, network)
│   │   └── README.md                               # Complete analytical guide & CTPS mathematical justification
│   │
│   ├── Hippocampus/                 # Hippocampus figure suite
│   │   ├── Fig1_Global_Landscape/                  # Fig1a to Fig1e (UMAPs, proportions, markers)
│   │   ├── Fig2_Cell_Type_Deep_Dives/              # Subfolders (Astrocyte, Microglia, etc.) with panels + DEGs
│   │   ├── Fig3_Cross_Cell_Convergence/            # Fig3a to Fig3g (Venn, UpSet, KEGG, violins, ranking)
│   │   ├── Fig4_Intercellular_Communication/       # Fig4a to Fig4e (CellChat networks, circle plots, CSVs)
│   │   ├── Fig5_Target_Validation/                 # Fig5a to Fig5c (Conservation, druggability, network)
│   │   └── README.md                               # Complete analytical guide & CTPS mathematical justification
│   │
│   ├── Cross_Region_Comparison/     # Direct head-to-head comparative analysis (Cortex vs HPC)
│   │   ├── Fig_CrossRegion_01_celltype_proportions_comparison.*
│   │   ├── Fig_CrossRegion_02_cortex_expression_violins.*
│   │   ├── Fig_CrossRegion_03_hippocampus_expression_violins.*
│   │   ├── Fig_CrossRegion_04_deg_concordance_scatter.* (Pearson r per cell type)
│   │   ├── Fig5a–Fig5c.* (Cross-species homology, surfaceome druggability, mechanistic network)
│   │   └── README.md
│   │
│   └── Quality_Control/            # Reviewer Defense & Empirical Single-Nucleus QC Suite
│       ├── Fig_QC_01_metrics_violin_by_condition.* (Genes, UMIs, MT%, Ribo% by Condition)
│       ├── Fig_QC_02_metrics_bivariate_scatter_and_thresholds.* (Library complexity & MT cutoffs)
│       ├── Fig_QC_03_donor_batch_variability_violins.* (Donor stability across 25 samples)
│       ├── Fig_QC_04_celltype_quality_metrics.* (QC integrity across 6 NVU lineages)
│       ├── Fig_QC_05_cell_retention_and_filtering_summary.* (Post-QC yields & SOP thresholds)
│       └── README.md                               # Complete reviewer rebuttal defense guide
```

---

## 📊 Publication Figures Catalog

All analytical panels are strictly formatted for direct journal submission:
- **No chart-internal titles** (preventing redundant or cluttered text headers).
- **Legible, high-contrast axes, ticks, and labels** formatted at ≥10 pt.
- **Both vector graphic (`.pdf`) and 300-DPI raster (`.png`) formats**.
- **Organized symmetrically** across both [`figures/Cortex/`](figures/Cortex) and [`figures/Hippocampus/`](figures/Hippocampus), with multi-region comparisons in [`figures/Cross_Region_Comparison/`](figures/Cross_Region_Comparison).

### Figure 1: Single-Nucleus Landscape of Human Prefrontal Cortex and Hippocampus
- **Fig1a**: Prefrontal Cortex single-nucleus UMAP projection colored by clinical condition (Control vs AD).
- **Fig1b**: Prefrontal Cortex UMAP split side-by-side showing resolved cell-type distributions.
- **Fig1c**: Prefrontal Cortex stacked cell-type proportions with exact percentage labels for all populations.
- **Fig1d**: Prefrontal Cortex donor-level proportion boxplots showing cell frequency shifts.
- **Fig1e**: Prefrontal Cortex canonical cell-type marker gene dotplot validating cluster identities.
- **Fig1f**: Hippocampus single-nucleus UMAP projection colored by clinical condition.
- **Fig1g**: Hippocampus UMAP split side-by-side by cell type.
- **Fig1h**: Hippocampus stacked cell-type proportions with exact percentages.
- **Fig1i**: Hippocampus donor-level proportion boxplots across conditions.
- **Fig1j**: Hippocampus canonical marker gene dotplot.

### Figure 2: Neurovascular Unit (NVU) Target Cell Deep Dives (Repeated for 5 Cell Types)
For each of **Microglia**, **Endothelial Cells**, **Pericytes**, **Astrocytes**, and **Oligodendrocytes**:
- **Fig2a**: Fine-grained subclustering UMAP projection.
- **Fig2b**: Subcluster-defining top marker gene expression dotplot.
- **Fig2c**: Subcluster distribution and cell counts across Control vs AD.
- **Fig2d**: AD vs Control differential expression volcano plot highlighting top statistically significant DEGs.
- **Fig2e**: Barplot of the top 10 upregulated and downregulated DEGs ranked by log2 fold-change.
- **Fig2f**: Top statistically significant KEGG enriched signaling pathways.
- **Fig2g**: Top statistically significant Gene Ontology Biological Process (GO:BP) terms.
- **Fig2h**: STRING-derived Protein–Protein Interaction (PPI) network with Louvain community detection.

### Figure 3: Cross-Cell Molecular Convergence & Multi-Factorial Drug Target Prioritization
- **Fig3a**: UpSet plot (`Fig3a_upset_5celltypes`) and Venn diagram (`Fig3a_venn_5celltypes`) showing DEG overlap across all 5 NVU lineages.
- **Fig3b**: KEGG pathway enrichment for the shared convergent molecular targets.
- **Fig3c**: Gene Ontology Biological Process terms enriched in convergent genes.
- **Fig3d**: High-confidence PPI interaction network among convergent targets.
- **Fig3e**: Split violin plots of candidate target expression in Prefrontal Cortex (Control vs AD across 5 cell types).
- **Fig3f**: Split violin plots of candidate target expression in Hippocampus (validating dual-region conservation).
- **Fig3g**: Prioritization ranking barplot based on the quantitative Multi-Factorial Consensus Score ($S_{\text{consensus}}$).

### Figure 4: Comparative Intercellular Signaling & NVU Communication Rewiring
- **Fig4a**: Total intercellular interaction number (Control: 113, AD: 98) and interaction strength (Control: 7.072, AD: 8.183).
- **Fig4b**: Differential circular interaction network displaying net signaling shifts (red edges = increased in AD, blue edges = decreased in AD).
- **Fig4c**: Relative information flow ranking across 27 signaling pathways identifying AD-enriched vs Control-enriched cascades.
- **Fig4d**: Side-by-side circular network of **APP** signaling showing profound hyper-activation targeting Microglia in AD.
- **Fig4e**: Ligand–receptor bubble plot highlighting communication probability and significance of critical NVU interactions (APP–SORL1, CD46–JAG1, CNTN1–NOTCH2, NRG3–ERBB4).

### Figure 5: In Silico Target Validation, Cross-Species Conservation & Antibody Selection
- **Fig5a**: Human-to-Mouse sequence identity (%) and In-Vivo Suitability Index (0–100) across prioritized therapeutic targets.
- **Fig5b**: Subcellular surfaceome classification and extracellular epitope accessibility for antibody targeting.
- **Fig5c**: Mechanistic network linking prioritized candidates (*SORL1, CLEC5A, ADAMTS9, PCDH9, DUSP1, NOTCH2*) directly to canonical AD drivers (*APP, BACE1, TYROBP, TREM2, CLDN5, MAPK14*).

---

## 🧬 Key Biological Findings

1. **Neurovascular Unit Disruption**:
   - Endothelial cells and pericytes exhibit pronounced downregulation of basement membrane and structural extracellular matrix pathways (*Laminin*, *Collagen* signaling pathways show severe depletion in AD).
   - Pericytes undergo loss of contractile and barrier-supportive transcripts, matching documented blood–brain barrier breakdown in AD.

2. **Cross-Cell Lineage Convergence**:
   - Immediate early stress response and transcriptional regulator genes (**FOS**, **JUN**, **DUSP1**, **EGR1**) and damage-response mediators (**GADD45B**) converge across all five NVU cell types, indicating a universal cellular stress program.
   - Dual-region validation in both Prefrontal Cortex and Hippocampus confirms that this convergent stress signature is preserved across anatomically distinct, vulnerable brain regions.

3. **Intercellular Signaling Rewiring in Alzheimer's Disease**:
   - While total interaction numbers decline slightly in AD (113 in Control vs 98 in AD), total interaction strength increases markedly (7.072 to 8.183), driven by concentrated hyper-signaling into microglia.
   - **APP–SORL1 Signaling**: Amyloid Precursor Protein (*APP*) ligand–receptor signaling from astrocytes and oligodendrocytes onto microglial and neuronal receptors is heavily augmented in AD.
   - **Neuroimmune & Inflammatory Checkpoints**: Increased information flow is concentrated in *CD46*, *PTPR*, *TGFb*, *PSAP*, *RA*, *NRXN*, *APP*, and *NOTCH* cascades, establishing clear mechanistic targets for therapeutic intervention.

---

## 🚀 How to Run the Pipeline

### Prerequisites
- R (version ≥ 4.3.0) with packages: `Seurat` (v5), `CellChat` (v2.2), `igraph`, `clusterProfiler`, `org.Hs.eg.db`, `ggplot2`, `patchwork`, `VennDiagram`, `ggalluvial`, `readr`, `dplyr`.

### Execution Order
From the root directory of the project, execute the numbered scripts in sequence:

```powershell
# Step 1: Discover files and parse sample metadata
Rscript scripts/01_data_exploration_and_metadata.R

# Step 2: Quality control, filtering, and reviewer defense figures
Rscript scripts/02a_process_cortex.R
Rscript scripts/02b_process_hippocampus.R
Rscript scripts/02c_qc_and_filtering_figures.R

# Step 3: Harmony batch integration, clustering, and cell typing
Rscript scripts/03a_integrate_cortex.R
Rscript scripts/03b_integrate_hippocampus.R

# Step 4: NVU target cell subclustering, DEGs, KEGG/GO, and PPI networks
Rscript scripts/04a_cortex_subclustering_and_degs.R
Rscript scripts/04b_hippocampus_subclustering_and_degs.R

# Step 5: Cross-cell convergence, dual-region validation, and CTPS target prioritization
Rscript scripts/05a_cortex_cross_cell_convergence.R
Rscript scripts/05b_hippocampus_cross_cell_convergence.R

# Step 6: CellChat v2.2 cell–cell communication modeling
Rscript scripts/06a_cortex_cell_cell_communication.R
Rscript scripts/06b_hippocampus_cell_cell_communication.R

# Step 7: In silico target validation and antibody selection guide
Rscript scripts/07a_cortex_target_validation.R
Rscript scripts/07b_hippocampus_target_validation.R

# Step 8: Head-to-head cross-region comparative analysis
Rscript scripts/08_cross_region_comparison.R
```

---

## 📄 License & Attribution
The raw data originates from GEO accession [GSE163577](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163577). All code and analytical figures in this repository were developed for the research study *"Single-Cell Transcriptomic Identification of Alzheimer’s Disease-Associated Molecular Targets and CNS–Peripheral Immune Signatures for Experimental Validation and Therapeutic Intervention"*.
