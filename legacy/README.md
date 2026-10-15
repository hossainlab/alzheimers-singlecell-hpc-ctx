# Legacy (superseded) pipeline — kept for provenance only

This folder holds the original R/Seurat/CellChat pipeline, its tables, figures and log, moved here with `git mv`
so the history is preserved. **Do not use these outputs for publication.** A scientific audit found the following
problems, all of which are corrected by the current pipeline in `../scripts/`:

| # | Problem in the legacy pipeline | Where | Fix in current pipeline |
|---|---|---|---|
| 1 | AD vs Control DE tested per nucleus (Wilcoxon `FindMarkers`) → pseudoreplication, padj down to 1e-300 with 4 v 4 donors; only 5–22 % of reported DEGs survive donor-level testing | `04a/04b`, `process_hippocampus_all_analyses.R` | Donor-level pseudobulk DESeq2 with sex and vascular-yield covariates (`05_pseudobulk_de.py`) |
| 2 | Top "targets" driven by donor genotype (ARL17B: 17q21.31 copy-number polymorphism, 3/4 AD vs 0/4 control carriers in cortex) | CTPS tables | 17q21.31-locus genes flagged and excluded from prioritisation; per-donor genotype recorded |
| 3 | "Cross-cell convergent" genes track sample composition (vascular enrichment yield is much lower in AD samples) → ambient RNA carry-over | `05a/05b` | Per-gene contamination flag + yield covariate; convergence recomputed on clean genes |
| 4 | Cortex arm presented as independent validation, but it comes from 8 of the same 17 donors | README, `08_cross_region_comparison.R` | Donor IDs harmonised; cortex reported as a second tissue from the same donors |
| 5 | Fig. 5 "in-silico validation" (sequence identity, 5xFAD concordance, ROSMAP replication, suitability scores, antibody catalogue) hard-coded as literals | `07a/07b`, `execute_revision_checklist.R` | Identity computed from UniProt sequences by alignment; localisation from UniProt; unsupported claims removed |
| 6 | "STRING PPI" networks were thresholded per-nucleus expression correlations | `04*`, `05*` | Real STRING v12 API network (combined score ≥ 0.7) |
| 7 | Annotated objects silently down-sampled to 5,000 / 3,500 nuclei per sample; log and QC table disagree; cached DEG tables re-used | `02*`, `03*`, `04*` | No down-sampling; every output regenerated from raw GEO matrices by committed code |
| 8 | README claims (FOS/JUN/DUSP1/EGR1/GADD45B convergence; APP signalling onto neurons in cortex; "prefrontal cortex") not supported by own tables | `README.md` | README rewritten from regenerated results |
| 9 | Cell-type proportion shifts interpreted biologically in a vascular-enrichment (VINE-seq) protocol | Fig. 1d | Composition shown per donor, explicitly labelled as enrichment-yield dependent |
| 10 | CellChat on 500 pooled nuclei per group, no donor replication | `06a/06b` | LIANA per sample with donor-level differential testing |

The uncommitted edits to `scripts/02a_process_cortex.R` and `scripts/02b_process_hippocampus.R` that were present
when this folder was created have been preserved as-is in `legacy/scripts/`.
