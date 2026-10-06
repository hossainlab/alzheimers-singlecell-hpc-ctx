# Data Directory (GSE163577)

This directory is configured to store raw single-nucleus RNA sequencing matrices and processed Seurat objects generated during the analysis.

> [!NOTE]
> Due to GitHub file size limits (>100 MB per file, repository limit <2 GB), large raw `.tar.gz` feature-barcode matrices and processed `.rds` Seurat objects are excluded from git tracking via `.gitignore`. All scripts are designed to reproduce these objects deterministically from the raw public accessions.

---

## 📥 Downloading Raw Sequencing Matrices

Raw 10x Genomics feature-barcode matrices originate from NCBI GEO accession **[GSE163577](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE163577)**:

- **Study**: *Single-nucleus RNA-seq of human prefrontal cortex and hippocampus in Alzheimer's disease*
- **Samples**: 25 total samples across healthy controls and AD donors (GSM4982083 through GSM4982107).

Place the downloaded raw matrix archives directly in:
```
data/GSE163577_RAW/
```

Files should follow the format:
- `GSM4982083_01_10_C_filtered_feature_bc_matrix.tar.gz`
- `GSM4982084_01_05_AD_filtered_feature_bc_matrix.tar.gz`
- ...
- `GSM4982107_O5_1OAD_Ctx_filtered_feature_bc_matrix.tar.gz`

---

## 🔄 Generated Intermediate & Processed Objects

Executing the pipeline scripts will generate the following `.rds` files inside this directory:

| Filename | Description | Generating Script |
| :--- | :--- | :--- |
| `seurat_cortex_qc.rds` | Filtered Prefrontal Cortex nuclei post-QC | `scripts/02a_process_cortex.R` |
| `seurat_hpc_qc.rds` | Filtered Hippocampus nuclei post-QC | `scripts/02b_process_hippocampus.R` |
| `seurat_cortex_annotated.rds` | Harmony batch-corrected, annotated Cortex Seurat object | `scripts/03a_integrate_cortex.R` |
| `seurat_hpc_annotated.rds` | Harmony batch-corrected, annotated Hippocampus Seurat object | `scripts/03b_integrate_hippocampus.R` |
| `fig2_results_cortex.rds` | Cortex NVU subclustering, DEGs, and pathway results | `scripts/04a_cortex_subclustering_and_degs.R` |
| `subclusters_Hippocampus_*.rds` | Lineage-specific Hippocampus subclustering objects | `scripts/04b_hippocampus_subclustering_and_degs.R` |
| `cellchat_focused5_merged.rds` | Cortex CellChat comparative communication networks | `scripts/06a_cortex_cell_cell_communication.R` |
| `cellchat_hpc_merged.rds` | Hippocampus CellChat comparative communication networks | `scripts/06b_hippocampus_cell_cell_communication.R` |
