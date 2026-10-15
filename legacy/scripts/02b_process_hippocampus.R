# ==============================================================================
# Script 02b: Process Hippocampus Samples (n = 17)
# GSM4982083 - GSM4982099 (8 Control, 9 AD)
#
# Quality control follows single-cell best practices (sc-best-practices / OSCA):
#   1. Per-sample QC metrics: library size, detected genes, MT%, ribosomal %, haemoglobin %
#   2. Per-sample, data-driven outlier detection with median absolute deviations (MADs)
#        - log10 library size and log10 detected genes: 5 MADs (two-sided, permissive)
#        - MT%: 3 MADs (upper), applied only above a 1% floor, plus a hard 5% ceiling.
#          Nuclear MT% is near zero, so the MAD can collapse to 0; the floor stops
#          nuclei with trace (<1%) MT reads being flagged as damaged.
#        - hard floor of 200 detected genes
#   3. Per-sample doublet detection with scDblFinder on QC-passing nuclei
#      Per-criterion failure counts (overlapping) are reported in the sample summary.
#   4. No downsampling: every nucleus passing QC and doublet removal is retained
#   5. Gene filtering after merge: genes detected in >= 10 nuclei
# Per-nucleus QC calls are written to tables/ for the QC figure suite (02c).
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(scuttle)
  library(scDblFinder)
  library(SingleCellExperiment)
  library(BiocParallel)
})

# QC parameters
nmads_counts   <- 5
nmads_genes    <- 5
nmads_mt       <- 3
max_mt_pct     <- 5
min_mt_pct     <- 1
min_genes      <- 200
min_gene_cells <- 10

cat("=== Starting Hippocampus Preprocessing ===\n")
meta_df <- read_csv("tables/GSE163577_sample_metadata.csv", show_col_types = FALSE)
hpc_meta <- meta_df %>% filter(region == "Hippocampus (Hpc)" | grepl("Hpc", region))

cat("Found", nrow(hpc_meta), "Hippocampus samples:\n")
print(hpc_meta[, c("gsm", "title", "condition", "cells")])

dir.create("data", showWarnings = FALSE)
dir.create("tables", showWarnings = FALSE)
sobj_list <- list()
qc_cells_list <- list()
qc_sample_list <- list()

# Process each Hippocampus sample
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

  sobj <- CreateSeuratObject(counts = counts, project = gsm, min.cells = 0, min.features = 0)
  sobj$GSM <- gsm
  sobj$Title <- row$title
  sobj$Region <- "Hippocampus"
  sobj$Condition <- row$condition
  sobj$Sample_ID <- paste0("Hpc_", row$condition, "_", gsm)

  # QC metrics (human gene symbols)
  sobj[["percent.mt"]]   <- PercentageFeatureSet(sobj, pattern = "^MT-")
  # Exact-match patterns: "^RP[SL][0-9]" also hits RPS6KA*/RPS6KB* kinases and "^HB[^(P)]"
  # hits HBEGF/HBS1L, so ribosomal and haemoglobin genes are matched explicitly.
  sobj[["percent.ribo"]] <- PercentageFeatureSet(sobj, pattern = "^RP[SL][0-9]+[A-Z]?[0-9]?$|^RPSA$|^RPLP[0-2]$")
  sobj[["percent.hb"]]   <- PercentageFeatureSet(sobj, pattern = "^HB(A1|A2|B|D|E1|G1|G2|M|Q1|Z)$")

  # Per-sample MAD-based outlier detection
  md <- sobj@meta.data
  out_counts <- isOutlier(md$nCount_RNA, nmads = nmads_counts, type = "both", log = TRUE)
  out_genes  <- isOutlier(md$nFeature_RNA, nmads = nmads_genes, type = "both", log = TRUE)
  out_mt     <- isOutlier(md$percent.mt, nmads = nmads_mt, type = "higher")
  thr_mt  <- min(max(attr(out_mt, "thresholds")[["higher"]], min_mt_pct), max_mt_pct)
  fail_umi   <- as.vector(out_counts)
  fail_genes <- as.vector(out_genes) | md$nFeature_RNA < min_genes
  fail_mt    <- md$percent.mt > thr_mt
  fail_qc    <- fail_umi | fail_genes | fail_mt

  # Doublet detection on QC-passing nuclei
  pass_cells <- colnames(sobj)[!fail_qc]
  sce <- SingleCellExperiment(list(counts = GetAssayData(sobj, layer = "counts")[, pass_cells]))
  set.seed(42)
  sce <- scDblFinder(sce, BPPARAM = SerialParam(RNGseed = 42))
  dbl_class <- setNames(as.character(sce$scDblFinder.class), colnames(sce))
  dbl_score <- setNames(sce$scDblFinder.score, colnames(sce))
  rm(sce)

  qc_status <- ifelse(fail_qc, "Low quality", ifelse(dbl_class[colnames(sobj)] == "doublet", "Doublet", "Pass"))
  sobj$qc_status <- qc_status
  sobj$scDblFinder_score <- unname(dbl_score[colnames(sobj)])

  qc_cells_list[[gsm]] <- sobj@meta.data %>%
    mutate(cell = paste0(gsm, "_", colnames(sobj))) %>%
    select(cell, GSM, Title, Region, Condition, nCount_RNA, nFeature_RNA,
           percent.mt, percent.ribo, percent.hb, scDblFinder_score, qc_status)

  thr_counts <- attr(out_counts, "thresholds")
  thr_genes  <- attr(out_genes, "thresholds")
  qc_sample_list[[gsm]] <- data.frame(
    GSM = gsm, Region = "Hippocampus", Condition = row$condition,
    Nuclei_Input = ncol(sobj),
    Low_Quality = sum(qc_status == "Low quality"),
    Doublets = sum(qc_status == "Doublet"),
    Nuclei_Retained = sum(qc_status == "Pass"),
    Fail_UMI = sum(fail_umi), Fail_Genes = sum(fail_genes), Fail_MT = sum(fail_mt),
    Min_UMI = thr_counts[["lower"]], Max_UMI = thr_counts[["higher"]],
    Min_Genes = max(thr_genes[["lower"]], min_genes), Max_Genes = thr_genes[["higher"]],
    Max_MT_Pct = thr_mt
  )

  sobj <- subset(sobj, cells = colnames(sobj)[qc_status == "Pass"])
  cat(sprintf("  Input: %d | Low quality: %d | Doublets: %d | Retained: %d (%.1f%%)\n",
              qc_sample_list[[gsm]]$Nuclei_Input, qc_sample_list[[gsm]]$Low_Quality,
              qc_sample_list[[gsm]]$Doublets, ncol(sobj),
              100 * ncol(sobj) / qc_sample_list[[gsm]]$Nuclei_Input))

  sobj_list[[gsm]] <- sobj
  gc()
}

qc_cells  <- bind_rows(qc_cells_list)
qc_sample <- bind_rows(qc_sample_list)
write_csv(qc_cells, "tables/qc_cell_metrics_hippocampus.csv.gz")
write_csv(qc_sample, "tables/qc_sample_summary_hippocampus.csv")
cat("\nPer-sample QC summary:\n")
print(qc_sample)

cat("\nMerging Hippocampus Seurat objects...\n")
hpc_merged <- merge(sobj_list[[1]], y = sobj_list[-1], add.cell.ids = names(sobj_list), project = "AD_Hippocampus")
rm(sobj_list)
gc()

# Gene filtering: keep genes detected in >= min_gene_cells nuclei
hpc_merged[["RNA"]] <- JoinLayers(hpc_merged[["RNA"]])
gene_cells <- Matrix::rowSums(GetAssayData(hpc_merged, layer = "counts") > 0)
hpc_merged <- subset(hpc_merged, features = names(gene_cells)[gene_cells >= min_gene_cells])

cat("Hippocampus merged object dimensions:\n")
print(hpc_merged)
cat(sprintf("Total cells: %d across %d samples\n", ncol(hpc_merged), length(unique(hpc_merged$GSM))))
print(table(hpc_merged$Condition, hpc_merged$GSM))

saveRDS(hpc_merged, "data/seurat_hpc_qc.rds")
cat("Saved to data/seurat_hpc_qc.rds successfully!\n")
