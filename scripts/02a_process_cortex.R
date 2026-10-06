# ==============================================================================
# Script 02a: Process Cortex Samples (n = 8)
# GSM4982100 - GSM4982107 (4 Control, 4 AD)
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(readr)
  library(ggplot2)
})

cat("=== Starting Cortex Preprocessing ===\n")
meta_df <- read_csv("tables/GSE163577_sample_metadata.csv", show_col_types = FALSE)
ctx_meta <- meta_df %>% filter(region == "Cortex (Ctx)" | grepl("Ctx", region))

cat("Found", nrow(ctx_meta), "Cortex samples:\n")
print(ctx_meta[, c("gsm", "title", "condition", "cells")])

dir.create("data", showWarnings = FALSE)
sobj_list <- list()

# Process each Cortex sample
for (i in 1:nrow(ctx_meta)) {
  row <- ctx_meta[i, ]
  gsm <- row$gsm
  raw_dir <- if (dir.exists("data/GSE163577_RAW")) "data/GSE163577_RAW" else "GSE163577_RAW"
  tar_path <- file.path(raw_dir, row$filename)
  
  cat(sprintf("\n[%d/%d] Reading %s (%s, %s)...\n", i, nrow(ctx_meta), gsm, row$title, row$condition))
  
  tmp_dir <- file.path(tempdir(), paste0("ctx_", gsm))
  if (dir.exists(tmp_dir)) unlink(tmp_dir, recursive = TRUE)
  dir.create(tmp_dir, recursive = TRUE)
  
  untar(tar_path, exdir = tmp_dir)
  sub_dirs <- list.dirs(tmp_dir, recursive = FALSE)
  data_dir <- if (length(sub_dirs) > 0) sub_dirs[1] else tmp_dir
  
  counts <- Read10X(data.dir = data_dir)
  unlink(tmp_dir, recursive = TRUE)
  
  # Nuclei QC
  sobj <- CreateSeuratObject(counts = counts, project = gsm, min.cells = 3, min.features = 200)
  sobj$GSM <- gsm
  sobj$Title <- row$title
  sobj$Region <- "Cortex"
  sobj$Condition <- row$condition
  sobj$Sample_ID <- paste0("Ctx_", row$condition, "_", gsm)
  
  sobj[["percent.mt"]] <- PercentageFeatureSet(sobj, pattern = "^MT-")
  sobj[["percent.ribo"]] <- PercentageFeatureSet(sobj, pattern = "^RP[SL]")
  
  # Standard nuclear snRNA-seq filter
  n_raw <- ncol(sobj)
  sobj <- subset(sobj, subset = nFeature_RNA >= 300 & nFeature_RNA <= 6500 & 
                                nCount_RNA >= 500 & nCount_RNA <= 25000 & 
                                percent.mt <= 5)
  n_clean <- ncol(sobj)
  cat(sprintf("  Raw: %d -> QC Passed: %d (%.1f%%)\n", n_raw, n_clean, 100 * n_clean / n_raw))
  
  # Cap per sample to 5,000 cells max if needed for memory & balance, keeping high-quality cells
  if (ncol(sobj) > 5000) {
    set.seed(42)
    sobj <- subset(sobj, cells = sample(colnames(sobj), 5000))
    cat(sprintf("  Sampled to 5,000 representative cells for balanced representation.\n"))
  }
  
  sobj_list[[gsm]] <- sobj
  gc()
}

cat("\nMerging Cortex Seurat objects...\n")
cortex_merged <- merge(sobj_list[[1]], y = sobj_list[-1], add.cell.ids = names(sobj_list), project = "AD_Cortex")
rm(sobj_list)
gc()

cat("Cortex merged object dimensions:\n")
print(cortex_merged)
cat(sprintf("Total cells: %d across %d samples\n", ncol(cortex_merged), length(unique(cortex_merged$GSM))))
print(table(cortex_merged$Condition, cortex_merged$GSM))

saveRDS(cortex_merged, "data/seurat_cortex_qc.rds")
cat("Saved to data/seurat_cortex_qc.rds successfully!\n")
