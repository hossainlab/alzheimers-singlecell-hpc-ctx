# ==============================================================================
# Script 02c: Single-Nucleus Quality Control & Filtering Figure Suite
# Generates comprehensive publication QC panels to satisfy peer review standards:
#   - Fig_QC_01: Global QC metrics violins by condition & region (genes, UMIs, MT%, Ribo%)
#   - Fig_QC_02: Bivariate QC scatter plots (Complexity & MT exclusion thresholds)
#   - Fig_QC_03: Donor-by-donor batch variability across all 25 samples
#   - Fig_QC_04: Lineage-specific QC distributions across NVU cell types
#   - Fig_QC_05: Nuclei retention & filtering thresholds summary
# Outputs to: figures/Quality_Control/
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(patchwork)
  library(cowplot)
  library(ggpubr)
  library(scales)
})

cat("======================================================================\n")
cat("Starting Step 02c: Generating Single-Nucleus QC Figure Suite\n")
cat("======================================================================\n")

out_dir <- "figures/Quality_Control"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Publication styling theme
theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_blank(),
      plot.subtitle = element_blank(),
      axis.title = element_text(face = "bold", size = rel(1.0), color = "black"),
      axis.text = element_text(color = "black", size = rel(0.9)),
      axis.line = element_line(color = "black", linewidth = 0.5),
      axis.ticks = element_line(color = "black", linewidth = 0.5),
      legend.title = element_text(face = "bold", size = rel(0.95), color = "black"),
      legend.text = element_text(size = rel(0.85), color = "black"),
      legend.background = element_blank(),
      legend.key = element_blank(),
      panel.grid.major = element_line(color = "grey94", linewidth = 0.3),
      panel.grid.minor = element_blank(),
      strip.background = element_rect(fill = "grey92", color = NA),
      strip.text = element_text(face = "bold", size = rel(0.95), color = "black")
    )
}

save_plot_pair <- function(plot_obj, out_path_no_ext, width = 8, height = 6) {
  pdf_file <- paste0(out_path_no_ext, ".pdf")
  png_file <- paste0(out_path_no_ext, ".png")
  ggsave(pdf_file, plot = plot_obj, width = width, height = height)
  ggsave(png_file, plot = plot_obj, width = width, height = height, dpi = 300)
  cat(sprintf("  Saved: %s (.pdf & .png)\n", basename(out_path_no_ext)))
}

# 1. Load Seurat objects and sample metadata
cat("Loading annotated Seurat objects for QC profiling...\n")
cortex <- readRDS("data/seurat_cortex_annotated.rds")
hpc    <- readRDS("data/seurat_hpc_annotated.rds")

meta_ctx <- cortex@meta.data %>%
  select(GSM, Title, Region, Condition, nCount_RNA, nFeature_RNA, percent.mt, percent.ribo, cell_type)
meta_hpc <- hpc@meta.data %>%
  select(GSM, Title, Region, Condition, nCount_RNA, nFeature_RNA, percent.mt, percent.ribo, cell_type)

meta_combined <- bind_rows(meta_ctx, meta_hpc)
meta_combined$Condition <- factor(meta_combined$Condition, levels = c("Control", "AD"))
meta_combined$Region <- factor(meta_combined$Region, levels = c("Cortex", "Hippocampus"))

cat(sprintf("Total profiled nuclei: %d (Cortex: %d, HPC: %d)\n", 
            nrow(meta_combined), nrow(meta_ctx), nrow(meta_hpc)))

# ------------------------------------------------------------------------------
# PANEL 1: Fig_QC_01 - Global QC Violins by Condition & Region
# ------------------------------------------------------------------------------
cat("\nGenerating Fig_QC_01: Global QC metrics violins by condition & region...\n")

p1_feat <- ggplot(meta_combined, aes(x = Condition, y = nFeature_RNA, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.7, scale = "width", color = "black", linewidth = 0.3) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, color = "black", linewidth = 0.4) +
  facet_wrap(~Region) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x.npc = "center", vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("Control" = "#4575b4", "AD" = "#d73027")) +
  scale_y_continuous(labels = comma) +
  labs(x = "Condition", y = "Detected Genes per Nucleus (nFeature_RNA)") +
  theme_pub(base_size = 10) + theme(legend.position = "none")

p1_count <- ggplot(meta_combined, aes(x = Condition, y = nCount_RNA, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.7, scale = "width", color = "black", linewidth = 0.3) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, color = "black", linewidth = 0.4) +
  facet_wrap(~Region) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x.npc = "center", vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("Control" = "#4575b4", "AD" = "#d73027")) +
  scale_y_continuous(labels = comma) +
  labs(x = "Condition", y = "Total UMI Counts per Nucleus (nCount_RNA)") +
  theme_pub(base_size = 10) + theme(legend.position = "none")

p1_mt <- ggplot(meta_combined, aes(x = Condition, y = percent.mt, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.7, scale = "width", color = "black", linewidth = 0.3) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, color = "black", linewidth = 0.4) +
  facet_wrap(~Region) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x.npc = "center", vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("Control" = "#4575b4", "AD" = "#d73027")) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "red", linewidth = 0.5) +
  labs(x = "Condition", y = "Mitochondrial Reads (%)") +
  theme_pub(base_size = 10) + theme(legend.position = "none")

p1_ribo <- ggplot(meta_combined, aes(x = Condition, y = percent.ribo, fill = Condition)) +
  geom_violin(trim = FALSE, alpha = 0.7, scale = "width", color = "black", linewidth = 0.3) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, color = "black", linewidth = 0.4) +
  facet_wrap(~Region) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x.npc = "center", vjust = -0.5, size = 4) +
  scale_fill_manual(values = c("Control" = "#4575b4", "AD" = "#d73027")) +
  labs(x = "Condition", y = "Ribosomal Reads (%)") +
  theme_pub(base_size = 10) + theme(legend.position = "none")

p_qc_01 <- (p1_feat | p1_count) / (p1_mt | p1_ribo)
save_plot_pair(p_qc_01, file.path(out_dir, "Fig_QC_01_metrics_violin_by_condition"), width = 11, height = 9)

# ------------------------------------------------------------------------------
# PANEL 2: Fig_QC_02 - Bivariate Scatter & Library Complexity Correlation
# ------------------------------------------------------------------------------
cat("\nGenerating Fig_QC_02: Bivariate QC scatter & filtering thresholds...\n")

set.seed(42)
sampled_meta <- meta_combined %>%
  group_by(Region, Condition) %>%
  slice_sample(n = 3000) %>%
  ungroup()

p2_scatter_genes <- ggplot(sampled_meta, aes(x = nCount_RNA, y = nFeature_RNA, color = Region)) +
  geom_point(alpha = 0.35, size = 0.8) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.6) +
  stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 4) +
  scale_color_manual(values = c("Cortex" = "#3C5488", "Hippocampus" = "#00A087")) +
  scale_x_continuous(labels = comma) +
  scale_y_continuous(labels = comma) +
  geom_hline(yintercept = c(300, 6500), linetype = "dashed", color = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = c(500, 25000), linetype = "dashed", color = "grey40", linewidth = 0.4) +
  labs(x = "Total UMI Counts (nCount_RNA)", y = "Detected Genes (nFeature_RNA)", color = "Brain Region") +
  theme_pub(base_size = 10) + theme(legend.position = "bottom")

p2_scatter_mt <- ggplot(sampled_meta, aes(x = nCount_RNA, y = percent.mt, color = Region)) +
  geom_point(alpha = 0.35, size = 0.8) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "red", linewidth = 0.6) +
  annotate("text", x = 18000, y = 5.2, label = "Upper QC Cutoff: 5% MT", color = "red", fontface = "bold", size = 3.5) +
  scale_color_manual(values = c("Cortex" = "#3C5488", "Hippocampus" = "#00A087")) +
  scale_x_continuous(labels = comma) +
  labs(x = "Total UMI Counts (nCount_RNA)", y = "Mitochondrial Reads (%)", color = "Brain Region") +
  theme_pub(base_size = 10) + theme(legend.position = "bottom")

p_qc_02 <- (p2_scatter_genes | p2_scatter_mt)
save_plot_pair(p_qc_02, file.path(out_dir, "Fig_QC_02_metrics_bivariate_scatter_and_thresholds"), width = 12, height = 5.5)

# ------------------------------------------------------------------------------
# PANEL 3: Fig_QC_03 - Donor-by-Donor Batch Consistency (All 25 Samples)
# ------------------------------------------------------------------------------
cat("\nGenerating Fig_QC_03: Donor-by-donor batch variability...\n")

meta_combined$GSM_label <- paste0(meta_combined$GSM, " (", substr(meta_combined$Region, 1, 3), ")")

p3_donor_genes <- ggplot(meta_combined, aes(x = reorder(GSM_label, nFeature_RNA, FUN = median), y = nFeature_RNA, fill = Condition)) +
  geom_violin(scale = "width", alpha = 0.7, color = "black", linewidth = 0.25) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, color = "black", linewidth = 0.3) +
  scale_fill_manual(values = c("Control" = "#4575b4", "AD" = "#d73027")) +
  scale_y_continuous(labels = comma) +
  labs(x = "Individual Donor Sample (GSM)", y = "Detected Genes per Nucleus") +
  theme_pub(base_size = 9) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8, color = "black"), legend.position = "top")

p3_donor_mt <- ggplot(meta_combined, aes(x = reorder(GSM_label, percent.mt, FUN = median), y = percent.mt, fill = Condition)) +
  geom_violin(scale = "width", alpha = 0.7, color = "black", linewidth = 0.25) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, color = "black", linewidth = 0.3) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "red", linewidth = 0.5) +
  scale_fill_manual(values = c("Control" = "#4575b4", "AD" = "#d73027")) +
  labs(x = "Individual Donor Sample (GSM)", y = "Mitochondrial Reads (%)") +
  theme_pub(base_size = 9) +
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1, size = 8, color = "black"), legend.position = "none")

p_qc_03 <- p3_donor_genes / p3_donor_mt
save_plot_pair(p_qc_03, file.path(out_dir, "Fig_QC_03_donor_batch_variability_violins"), width = 14, height = 9.5)

# ------------------------------------------------------------------------------
# PANEL 4: Fig_QC_04 - Lineage-Specific QC Distributions (6 NVU Cell Types)
# ------------------------------------------------------------------------------
cat("\nGenerating Fig_QC_04: Lineage-specific QC metrics across NVU cell types...\n")

target_cts <- c("Astrocyte", "Endothelial", "Inhibitory neuron", "Microglia", "Oligodendrocyte", "Pericyte")
meta_target_cts <- meta_combined %>% filter(cell_type %in% target_cts)

palette_6 <- c(
  "Astrocyte"         = "#3C5488",
  "Endothelial"       = "#4DBBD5",
  "Inhibitory neuron" = "#8491B4",
  "Microglia"         = "#E64B35",
  "Oligodendrocyte"   = "#F39B7F",
  "Pericyte"          = "#00A087"
)

p4_ct_genes <- ggplot(meta_target_cts, aes(x = cell_type, y = nFeature_RNA, fill = cell_type)) +
  geom_violin(scale = "width", alpha = 0.75, color = "black", linewidth = 0.3) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, color = "black", linewidth = 0.4) +
  scale_fill_manual(values = palette_6) +
  scale_y_continuous(labels = comma) +
  labs(x = "Cell Lineage", y = "Detected Genes per Nucleus") +
  theme_pub(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, face = "bold", color = "black"), legend.position = "none")

p4_ct_umis <- ggplot(meta_target_cts, aes(x = cell_type, y = nCount_RNA, fill = cell_type)) +
  geom_violin(scale = "width", alpha = 0.75, color = "black", linewidth = 0.3) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, color = "black", linewidth = 0.4) +
  scale_fill_manual(values = palette_6) +
  scale_y_continuous(labels = comma) +
  labs(x = "Cell Lineage", y = "Total UMI Counts per Nucleus") +
  theme_pub(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, face = "bold", color = "black"), legend.position = "none")

p4_ct_mt <- ggplot(meta_target_cts, aes(x = cell_type, y = percent.mt, fill = cell_type)) +
  geom_violin(scale = "width", alpha = 0.75, color = "black", linewidth = 0.3) +
  geom_boxplot(width = 0.15, fill = "white", outlier.shape = NA, color = "black", linewidth = 0.4) +
  geom_hline(yintercept = 5, linetype = "dashed", color = "red", linewidth = 0.5) +
  scale_fill_manual(values = palette_6) +
  labs(x = "Cell Lineage", y = "Mitochondrial Reads (%)") +
  theme_pub(base_size = 10) +
  theme(axis.text.x = element_text(angle = 35, hjust = 1, face = "bold", color = "black"), legend.position = "none")

p_qc_04 <- (p4_ct_genes | p4_ct_umis | p4_ct_mt)
save_plot_pair(p_qc_04, file.path(out_dir, "Fig_QC_04_celltype_quality_metrics"), width = 14, height = 5.2)

# ------------------------------------------------------------------------------
# PANEL 5: Fig_QC_05 - Nuclei Retention & Filtering Thresholds Summary
# ------------------------------------------------------------------------------
cat("\nGenerating Fig_QC_05: Nuclei retention & filtering thresholds summary...\n")

sample_meta <- read_csv("tables/GSE163577_sample_metadata.csv", show_col_types = FALSE)

retention_df <- meta_combined %>%
  group_by(GSM) %>%
  summarise(Post_QC_Nuclei = n(), .groups = "drop") %>%
  left_join(sample_meta %>% select(gsm, region, condition, cells), by = c("GSM" = "gsm")) %>%
  rename(Initial_Nuclei = cells, Region = region, Condition = condition) %>%
  mutate(
    Retention_Rate_Pct = round(Post_QC_Nuclei / Initial_Nuclei * 100, 1),
    Region_Short = ifelse(grepl("Ctx", Region), "Cortex", "Hippocampus")
  )

p5_retention <- ggplot(retention_df, aes(x = reorder(GSM, Retention_Rate_Pct), y = Retention_Rate_Pct, fill = Condition)) +
  geom_col(width = 0.65, color = "black", linewidth = 0.3) +
  geom_text(aes(label = sprintf("%.1f%%", Retention_Rate_Pct)), hjust = -0.15, size = 2.9, fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values = c("Control" = "#4575b4", "AD" = "#d73027")) +
  scale_y_continuous(limits = c(0, 115), expand = c(0, 0)) +
  labs(x = "Sample ID (GSM)", y = "Nuclei Retention Rate Post-QC (%)") +
  theme_pub(base_size = 9) +
  theme(axis.text.y = element_text(size = 8, color = "black"), legend.position = "bottom")

# Threshold summary infocard table
threshold_table <- data.frame(
  Parameter = c("Brain Region", "Sequencing Platform", "Mitochondrial Cutoff", "Gene Detection (Min)", 
                "Gene Detection (Max)", "UMI Count (Min)", "UMI Count (Max)", "Mean Retention Rate"),
  Standard_Applied = c("Human Prefrontal Cortex & Hippocampus", "10x Genomics Chromium 3' snRNA-seq", 
                      "percent.mt <= 5.0% (snRNA-seq Standard)", "nFeature_RNA >= 300 genes", 
                      "nFeature_RNA <= 6,500 genes (Doublet removal)", "nCount_RNA >= 500 UMIs (Debris removal)", 
                      "nCount_RNA <= 25,000 UMIs", sprintf("%.1f%% across 25 donors", mean(retention_df$Retention_Rate_Pct))),
  stringsAsFactors = FALSE
)

p5_table <- ggtexttable(
  threshold_table,
  rows = NULL,
  theme = ttheme(
    base_style = "classic",
    base_size = 9,
    padding = unit(c(4, 4), "mm"),
    colnames.style = colnames_style(color = "white", fill = "#3C5488", face = "bold"),
    tbody.style = tbody_style(fill = c("grey96", "white"))
  )
)

p_qc_05 <- (p5_retention | p5_table) + plot_layout(widths = c(1.1, 1.2))
save_plot_pair(p_qc_05, file.path(out_dir, "Fig_QC_05_cell_retention_and_filtering_summary"), width = 14, height = 7.5)

write_csv(retention_df, "tables/sample_qc_retention_rates.csv")
cat("Retention rate metrics saved to tables/sample_qc_retention_rates.csv\n")

cat("\n======================================================================\n")
cat("STEP 02c: QC FIGURE SUITE GENERATION COMPLETE!\n")
cat("======================================================================\n")
