# ==============================================================================
# Script 02b: Process Hippocampus Samples (n = 17)
# GSM4982083 - GSM4982099 (8 Control, 9 AD)
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

cat("=== Starting Hippocampus Preprocessing ===\n")
meta_df <- read_csv("tables/GSE163577_sample_metadata.csv", show_col_types = FALSE)
hpc_meta <- meta_df %>% filter(region == "Hippocampus (Hpc)" | grepl("Hpc", region))

cat("Found", nrow(hpc_meta), "Hippocampus samples:\n")
print(hpc_meta[, c("gsm", "title", "condition", "cells")])

dir.create("data", showWarnings = FALSE)
sobj_list <- list()

for (i in 1:nrow(hpc_meta)) {
  row <- hpc_meta[i, ]
  gsm <- row$gsm
  raw_dir <- if (dir.exists("data/GSE163577_RAW")) "data/GSE163577_RAW" else "GSE163577_RAW"
  tar_path <- file.path(raw_dir, row$filename)
  
  cat(sprintf("\n[%d/%d] Reading %s (%s, %s)...\n", i, nrow(hpc_meta), gsm, row$title, row$condition))
  
  tmp_dir <- file.path(tempdir(), paste0("hpc_", gsm))
  if (dir.exists(tmp_dir)) unlink(tmp_dir, recursive = TRUE)
  dir.create(tmp_dir, recursive = TRUE)
  
  untar(tar_path, exdir = tmp_dir)
  sub_dirs <- list.dirs(tmp_dir, recursive = FALSE)
  data_dir <- if (length(sub_dirs) > 0) sub_dirs[1] else tmp_dir
  
  counts <- Read10X(data.dir = data_dir)
  unlink(tmp_dir, recursive = TRUE)
  
  sobj <- CreateSeuratObject(counts = counts, project = gsm, min.cells = 3, min.features = 200)
  sobj$GSM <- gsm
  sobj$Title <- row$title
  sobj$Region <- "Hippocampus"
  sobj$Condition <- row$condition
  sobj$Sample_ID <- paste0("Hpc_", row$condition, "_", gsm)
  
  sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")
  sobj[["percent.ribo"]] <- PercentageFeatureSet(sobj, pattern = "^RP[SL]")
  
  n_raw <- ncol(sobj)
  sobj <- subset(sobj, subset = nFeature_RNA >= 300 & nFeature_RNA <= 6500 & 
                                nCount_RNA >= 500 & nCount_RNA <= 25000 & 
                                percent.mt <= 5)
  n_clean <- ncol(sobj)
  cat(sprintf("  Raw: %d -> QC Passed: %d (%.1f%%)\n", n_raw, n_clean, 100 * n_clean / n_raw))
  
  # Sample up to 3,500 cells per sample to maintain balanced multi-donor representation
  if (ncol(sobj) > 3500) {
    set.seed(42)
    sobj <- subset(sobj, cells = sample(colnames(sobj), 3500))
    cat(sprintf("  Sampled to 3,500 representative cells.\n"))
  }
  
  sobj_list[[gsm]] <- sobj
  gc()
}

cat("\nMerging Hippocampus Seurat objects...\n")
hpc_merged <- merge(sobj_list[[1]], y = sobj_list[-1], add.cell.ids = names(sobj_list), project = "AD_Hippocampus")
rm(sobj_list)
gc()

cat("Hippocampus merged object dimensions:\n")
print(hpc_merged)
cat(sprintf("Total cells: %d across %d samples\n", ncol(hpc_merged), length(unique(hpc_merged$GSM))))
print(table(hpc_merged$Condition, hpc_merged$GSM))

saveRDS(hpc_merged, "data/seurat_hpc_qc.rds")
cat("Saved to data/seurat_hpc_qc.rds successfully!\n")
