# ==============================================================================
# Script 02c: Single-Nucleus Quality Control Figure Suite
# Run after 02a/02b (per-nucleus QC calls) and 03a/03b (cell type annotation).
# Every figure is a single panel (no multi-panel layouts, no panel tags, no titles,
# no background grid), exported as .pdf and .png, one set per brain region:
#   Fig_QC_01  UMI counts per sample before filtering, with per-sample MAD thresholds
#   Fig_QC_02  Detected genes per sample before filtering, with per-sample MAD thresholds
#   Fig_QC_03  Mitochondrial % per sample before filtering, with per-sample threshold
#   Fig_QC_04  UMI counts vs detected genes, coloured by QC outcome
#   Fig_QC_05  UMI counts vs mitochondrial %, coloured by QC outcome
#   Fig_QC_06  Nuclei retained / removed per sample (low quality, doublet)
#   Fig_QC_07  Donor-level median detected genes, Control vs AD (Wilcoxon on donors)
#   Fig_QC_08  Donor-level median mitochondrial %, Control vs AD (Wilcoxon on donors)
#   Fig_QC_09  Detected genes per nucleus by annotated cell type (post-QC)
#   Fig_QC_10  Mitochondrial % per nucleus by annotated cell type (post-QC)
# Outputs to: figures/Quality_Control/
# ==============================================================================
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(ggpubr)
  library(ggrastr)
  library(scales)
})

cat("======================================================================\n")
cat("Starting Step 02c: Generating Single-Nucleus QC Figure Suite\n")
cat("======================================================================\n")

out_dir <- "figures/Quality_Control"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Publication theme: white background, no grid, no titles
theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      plot.tag = element_blank(),
      axis.title = element_text(face = "bold", size = rel(1.0), color = "black"),
      axis.text = element_text(color = "black", size = rel(0.9)),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      legend.title = element_text(face = "bold", size = rel(0.95), color = "black"),
      legend.text = element_text(size = rel(0.85), color = "black"),
      legend.background = element_blank(),
      legend.key = element_blank(),
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

save_plot_pair <- function(plot_obj, out_path_no_ext, width = 7, height = 5.5) {
  ggsave(paste0(out_path_no_ext, ".pdf"), plot = plot_obj, width = width, height = height)
  ggsave(paste0(out_path_no_ext, ".png"), plot = plot_obj, width = width, height = height, dpi = 300)
  cat(sprintf("  Saved: %s (.pdf & .png)\n", basename(out_path_no_ext)))
}

cond_cols   <- c("Control" = "#4575b4", "AD" = "#d73027")
status_cols <- c("Pass" = "grey65", "Low quality" = "#E64B35", "Doublet" = "#3C5488")

# Remove figures from the previous multi-panel QC layout
old_files <- list.files(out_dir, pattern = "^Fig_QC_.*\\.(png|pdf)$", full.names = TRUE)
if (length(old_files) > 0) file.remove(old_files)

regions <- list(
  cortex      = list(label = "Cortex",      rds = "data/seurat_cortex_annotated.rds"),
  hippocampus = list(label = "Hippocampus", rds = "data/seurat_hpc_annotated.rds")
)

retention_all <- list()

for (key in names(regions)) {
  reg <- regions[[key]]
  cat(sprintf("\n>>> %s <<<\n", reg$label))

  qc_cells  <- read_csv(sprintf("tables/qc_cell_metrics_%s.csv.gz", key), show_col_types = FALSE)
  qc_sample <- read_csv(sprintf("tables/qc_sample_summary_%s.csv", key), show_col_types = FALSE)

  # Sample labels ordered by condition then GSM
  sample_order <- qc_sample %>% arrange(factor(Condition, levels = c("Control", "AD")), GSM) %>% pull(GSM)
  qc_cells <- qc_cells %>%
    mutate(GSM = factor(GSM, levels = sample_order),
           Condition = factor(Condition, levels = c("Control", "AD")),
           qc_status = factor(qc_status, levels = c("Pass", "Low quality", "Doublet")))
  qc_sample <- qc_sample %>%
    mutate(GSM = factor(GSM, levels = sample_order),
           Condition = factor(Condition, levels = c("Control", "AD")))

  x_sample <- theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

  # Fig_QC_01: UMI counts per sample with MAD thresholds
  p01 <- ggplot(qc_cells, aes(x = GSM, y = nCount_RNA, fill = Condition)) +
    geom_violin(scale = "width", alpha = 0.75, color = "black", linewidth = 0.25) +
    geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA, color = "black", linewidth = 0.3) +
    geom_errorbar(data = qc_sample, aes(x = GSM, ymin = Min_UMI, ymax = Min_UMI), inherit.aes = FALSE,
                  width = 0.7, color = "black", linetype = "dashed", linewidth = 0.5) +
    geom_errorbar(data = qc_sample, aes(x = GSM, ymin = Max_UMI, ymax = Max_UMI), inherit.aes = FALSE,
                  width = 0.7, color = "black", linetype = "dashed", linewidth = 0.5) +
    scale_y_log10(labels = comma) +
    scale_fill_manual(values = cond_cols) +
    labs(x = "Sample (GSM)", y = "UMI counts per nucleus", fill = "Condition") +
    theme_pub() + x_sample
  save_plot_pair(p01, file.path(out_dir, sprintf("Fig_QC_01_%s_umi_per_sample", key)),
                 width = max(6, 0.45 * length(sample_order) + 2.5))

  # Fig_QC_02: detected genes per sample with MAD thresholds
  p02 <- ggplot(qc_cells, aes(x = GSM, y = nFeature_RNA, fill = Condition)) +
    geom_violin(scale = "width", alpha = 0.75, color = "black", linewidth = 0.25) +
    geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA, color = "black", linewidth = 0.3) +
    geom_errorbar(data = qc_sample, aes(x = GSM, ymin = Min_Genes, ymax = Min_Genes), inherit.aes = FALSE,
                  width = 0.7, color = "black", linetype = "dashed", linewidth = 0.5) +
    geom_errorbar(data = qc_sample, aes(x = GSM, ymin = Max_Genes, ymax = Max_Genes), inherit.aes = FALSE,
                  width = 0.7, color = "black", linetype = "dashed", linewidth = 0.5) +
    scale_y_log10(labels = comma) +
    scale_fill_manual(values = cond_cols) +
    labs(x = "Sample (GSM)", y = "Detected genes per nucleus", fill = "Condition") +
    theme_pub() + x_sample
  save_plot_pair(p02, file.path(out_dir, sprintf("Fig_QC_02_%s_genes_per_sample", key)),
                 width = max(6, 0.45 * length(sample_order) + 2.5))

  # Fig_QC_03: mitochondrial % per sample with threshold
  p03 <- ggplot(qc_cells, aes(x = GSM, y = percent.mt, fill = Condition)) +
    geom_violin(scale = "width", alpha = 0.75, color = "black", linewidth = 0.25) +
    geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA, color = "black", linewidth = 0.3) +
    geom_errorbar(data = qc_sample, aes(x = GSM, ymin = Max_MT_Pct, ymax = Max_MT_Pct), inherit.aes = FALSE,
                  width = 0.7, color = "black", linetype = "dashed", linewidth = 0.5) +
    scale_fill_manual(values = cond_cols) +
    labs(x = "Sample (GSM)", y = "Mitochondrial reads (%)", fill = "Condition") +
    theme_pub() + x_sample
  save_plot_pair(p03, file.path(out_dir, sprintf("Fig_QC_03_%s_mt_pct_per_sample", key)),
                 width = max(6, 0.45 * length(sample_order) + 2.5))

  # Scatter plots: draw passing nuclei first so removed nuclei stay visible
  qc_scatter <- qc_cells %>% arrange(qc_status)
  status_n <- table(qc_cells$qc_status)
  status_labels <- sprintf("%s (n = %s)", names(status_n), comma(as.numeric(status_n)))
  names(status_labels) <- names(status_n)

  # Fig_QC_04: UMI vs genes coloured by QC outcome
  p04 <- ggplot(qc_scatter, aes(x = nCount_RNA, y = nFeature_RNA, color = qc_status)) +
    rasterise(geom_point(size = 0.3, alpha = 0.5), dpi = 300) +
    scale_x_log10(labels = comma) +
    scale_y_log10(labels = comma) +
    scale_color_manual(values = status_cols, labels = status_labels) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(x = "UMI counts per nucleus", y = "Detected genes per nucleus", color = "QC outcome") +
    theme_pub()
  save_plot_pair(p04, file.path(out_dir, sprintf("Fig_QC_04_%s_umi_vs_genes_qc_outcome", key)))

  # Fig_QC_05: UMI vs MT% coloured by QC outcome
  p05 <- ggplot(qc_scatter, aes(x = nCount_RNA, y = percent.mt, color = qc_status)) +
    rasterise(geom_point(size = 0.3, alpha = 0.5), dpi = 300) +
    scale_x_log10(labels = comma) +
    scale_color_manual(values = status_cols, labels = status_labels) +
    guides(color = guide_legend(override.aes = list(size = 3, alpha = 1))) +
    labs(x = "UMI counts per nucleus", y = "Mitochondrial reads (%)", color = "QC outcome") +
    theme_pub()
  save_plot_pair(p05, file.path(out_dir, sprintf("Fig_QC_05_%s_umi_vs_mt_pct_qc_outcome", key)))

  # Fig_QC_06: nuclei retained / removed per sample
  retention_long <- qc_sample %>%
    select(GSM, Condition, Pass = Nuclei_Retained, `Low quality` = Low_Quality, Doublet = Doublets) %>%
    pivot_longer(c(Pass, `Low quality`, Doublet), names_to = "qc_status", values_to = "Nuclei") %>%
    mutate(qc_status = factor(qc_status, levels = c("Doublet", "Low quality", "Pass")))
  p06 <- ggplot(retention_long, aes(x = GSM, y = Nuclei, fill = qc_status)) +
    geom_col(width = 0.75, color = "black", linewidth = 0.2) +
    scale_fill_manual(values = status_cols, breaks = c("Pass", "Low quality", "Doublet")) +
    scale_y_continuous(labels = comma, expand = expansion(mult = c(0, 0.03))) +
    labs(x = "Sample (GSM)", y = "Number of nuclei", fill = "QC outcome") +
    theme_pub() + x_sample
  save_plot_pair(p06, file.path(out_dir, sprintf("Fig_QC_06_%s_nuclei_retained_per_sample", key)),
                 width = max(6, 0.45 * length(sample_order) + 2.5))

  # Donor-level medians of retained nuclei (donor is the unit of replication)
  donor_qc <- qc_cells %>%
    filter(qc_status == "Pass") %>%
    group_by(GSM, Condition) %>%
    summarise(median_genes = median(nFeature_RNA), median_mt = median(percent.mt), .groups = "drop")

  donor_plot <- function(yvar, ylab) {
    ggplot(donor_qc, aes(x = Condition, y = .data[[yvar]], fill = Condition)) +
      geom_boxplot(width = 0.55, alpha = 0.75, outlier.shape = NA, color = "black", linewidth = 0.4) +
      geom_jitter(width = 0.12, size = 2.2, shape = 21, color = "black", stroke = 0.3) +
      stat_compare_means(method = "wilcox.test", label = "p.format", label.x.npc = "center", size = 3.8) +
      scale_fill_manual(values = cond_cols) +
      labs(x = "Condition", y = ylab) +
      theme_pub() + theme(legend.position = "none")
  }
  save_plot_pair(donor_plot("median_genes", "Donor median detected genes per nucleus"),
                 file.path(out_dir, sprintf("Fig_QC_07_%s_donor_median_genes", key)), width = 4, height = 5)
  save_plot_pair(donor_plot("median_mt", "Donor median mitochondrial reads (%)"),
                 file.path(out_dir, sprintf("Fig_QC_08_%s_donor_median_mt_pct", key)), width = 4, height = 5)

  # Post-QC metrics by annotated cell type (all annotated cell types)
  sobj <- readRDS(reg$rds)
  ct_meta <- sobj@meta.data %>% select(cell_type, nFeature_RNA, percent.mt)
  ct_order <- ct_meta %>% group_by(cell_type) %>% summarise(m = median(nFeature_RNA)) %>% arrange(desc(m)) %>% pull(cell_type)
  ct_meta$cell_type <- factor(ct_meta$cell_type, levels = ct_order)
  rm(sobj); gc()

  x_ct <- theme(axis.text.x = element_text(angle = 40, hjust = 1, vjust = 1))
  p09 <- ggplot(ct_meta, aes(x = cell_type, y = nFeature_RNA, fill = cell_type)) +
    geom_violin(scale = "width", alpha = 0.75, color = "black", linewidth = 0.25) +
    geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA, color = "black", linewidth = 0.3) +
    scale_y_log10(labels = comma) +
    labs(x = "Cell type", y = "Detected genes per nucleus") +
    theme_pub() + x_ct + theme(legend.position = "none")
  save_plot_pair(p09, file.path(out_dir, sprintf("Fig_QC_09_%s_celltype_genes", key)),
                 width = max(6, 0.6 * length(ct_order) + 2))

  p10 <- ggplot(ct_meta, aes(x = cell_type, y = percent.mt, fill = cell_type)) +
    geom_violin(scale = "width", alpha = 0.75, color = "black", linewidth = 0.25) +
    geom_boxplot(width = 0.12, fill = "white", outlier.shape = NA, color = "black", linewidth = 0.3) +
    labs(x = "Cell type", y = "Mitochondrial reads (%)") +
    theme_pub() + x_ct + theme(legend.position = "none")
  save_plot_pair(p10, file.path(out_dir, sprintf("Fig_QC_10_%s_celltype_mt_pct", key)),
                 width = max(6, 0.6 * length(ct_order) + 2))

  retention_all[[key]] <- qc_sample %>%
    mutate(Retention_Rate_Pct = round(Nuclei_Retained / Nuclei_Input * 100, 1))
}

retention_df <- bind_rows(retention_all)
write_csv(retention_df, "tables/sample_qc_retention_rates.csv")
cat("\nRetention summary saved to tables/sample_qc_retention_rates.csv\n")
print(retention_df %>% group_by(Region) %>%
        summarise(Input = sum(Nuclei_Input), Low_Quality = sum(Low_Quality),
                  Doublets = sum(Doublets), Retained = sum(Nuclei_Retained)))

cat("\n======================================================================\n")
cat("STEP 02c: QC FIGURE SUITE GENERATION COMPLETE!\n")
cat("======================================================================\n")
