# ==============================================================================
# Script 01: GSE163577 Data Exploration and Metadata Generation
# Project: Single-Cell Transcriptomic Identification of AD Molecular Targets 
#          and CNS-Peripheral Immune Signatures
# ==============================================================================

suppressPackageStartupMessages({
  library(jsonlite)
  library(dplyr)
  library(readr)
})

cat("=== Step 1: Scanning GSE163577_RAW archive files ===\n")
raw_dir <- if (dir.exists("data/GSE163577_RAW")) "data/GSE163577_RAW" else "GSE163577_RAW"
tar_files <- list.files(raw_dir, pattern = "\\.tar\\.gz$", full.names = TRUE)

if (length(tar_files) == 0) {
  stop(sprintf("No .tar.gz files found in %s directory!", raw_dir))
}

cat("Found", length(tar_files), "raw sample archives in", raw_dir, "\n")

# Load existing parsed metadata if available, otherwise generate
meta_json <- "tables/sample_metadata.json"
cell_json <- "tables/sample_cell_counts.json"

if (file.exists(meta_json) && file.exists(cell_json)) {
  meta_df <- fromJSON(meta_json)
  cell_df <- fromJSON(cell_json)
  master_meta <- left_df <- merge(meta_df, cell_df[, c("gsm", "cells")], by = "gsm", all.x = TRUE)
} else {
  # Build metadata directly
  sample_info <- lapply(tar_files, function(f) {
    fname <- basename(f)
    gsm <- strsplit(fname, "_")[[1]][1]
    
    # Brain region
    region <- if (grepl("Ctx", fname, ignore.case = TRUE)) "Cortex" else "Hippocampus"
    
    # Condition
    cond <- if (grepl("_AD_|AD_Ctx|H2004_AD", fname, ignore.case = TRUE)) "AD" else "Control"
    
    # Donor ID estimation
    donor <- sub("_filtered_feature_bc_matrix.tar.gz", "", fname)
    donor <- sub(paste0(gsm, "_"), "", donor)
    
    data.frame(
      gsm = gsm,
      filename = fname,
      region = region,
      condition = cond,
      sample_id = donor,
      size_mb = round(file.info(f)$size / (1024 * 1024), 2),
      stringsAsFactors = FALSE
    )
  })
  master_meta <- do.call(rbind, sample_info)
}

# Create output directories if needed
dir.create("tables", showWarnings = FALSE)
dir.create("figures", showWarnings = FALSE)
dir.create("scripts", showWarnings = FALSE)

# Save master sample metadata table
write_csv(master_meta, "tables/GSE163577_sample_metadata.csv")
cat("\nSample Metadata Summary:\n")
print(table(master_meta$region, master_meta$condition))

cat("\nSummary saved to: tables/GSE163577_sample_metadata.csv\n")
cat("Step 1 complete!\n")
