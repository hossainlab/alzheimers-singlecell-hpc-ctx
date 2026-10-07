# ==============================================================================
# Script 08: Cross-Region Comparison Suite (Cortex vs Hippocampus)
# Head-to-head comparison between human Prefrontal Cortex and Hippocampus
# - Cross-region cell-type composition comparison
# - Cross-region target expression violins
# - Cross-region DEG concordance scatter plot (Pearson correlation r)
# - In silico target validation panels
# Outputs to: figures/Cross_Region_Comparison/
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
  library(ggrepel)
})

cat("======================================================================\n")
cat("Starting Step 08: Cross-Region Comparison (Cortex vs Hippocampus)\n")
cat("======================================================================\n")

out_dir <- "figures/Cross_Region_Comparison"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

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

save_plot_pair <- function(plot_obj, out_path_no_ext, width = 7, height = 6) {
  pdf_file <- paste0(out_path_no_ext, ".pdf")
  png_file <- paste0(out_path_no_ext, ".png")
  ggsave(pdf_file, plot = plot_obj, width = width, height = height)
  ggsave(png_file, plot = plot_obj, width = width, height = height, dpi = 300)
  cat(sprintf("  Saved: %s (.pdf & .png)\n", basename(out_path_no_ext)))
}

cat("Loading annotated Seurat objects...\n")
cortex <- readRDS("data/seurat_cortex_annotated.rds")
hpc    <- readRDS("data/seurat_hpc_annotated.rds")

target_cts <- c("Astrocyte", "Endothelial", "Inhibitory neuron", "Microglia", "Oligodendrocyte", "Pericyte")
palette_6 <- c(
  "Astrocyte"         = "#3C5488",
  "Endothelial"       = "#4DBBD5",
  "Inhibitory neuron" = "#8491B4",
  "Microglia"         = "#E64B35",
  "Oligodendrocyte"   = "#F39B7F",
  "Pericyte"          = "#00A087"
)

# 1. Cross-Region Cell-Type Proportions Comparison
prop_ctx <- cortex@meta.data %>%
  filter(cell_type %in% target_cts) %>%
  group_by(Condition, cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Condition) %>%
  mutate(Proportion = n / sum(n) * 100, Region = "Cortex")

prop_hpc <- hpc@meta.data %>%
  filter(cell_type %in% target_cts) %>%
  group_by(Condition, cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Condition) %>%
  mutate(Proportion = n / sum(n) * 100, Region = "Hippocampus")

prop_all <- bind_rows(prop_ctx, prop_hpc)
prop_all$Region <- factor(prop_all$Region, levels = c("Cortex", "Hippocampus"))
prop_all$Condition <- factor(prop_all$Condition, levels = c("Control", "AD"))
prop_all$cell_type <- factor(prop_all$cell_type, levels = target_cts)

p_prop_comp <- ggplot(prop_all, aes(x = Condition, y = Proportion, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.72, color = "black", linewidth = 0.3) +
  facet_wrap(~Region) +
  scale_x_discrete(expand = expansion(mult = c(0.2, 0.2))) +
  scale_fill_manual(values = palette_6, name = "Cell Type") +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100.5)) +
  labs(x = "Condition", y = "Lineage Proportion (%)") +
  theme_pub() +
  theme(strip.text = element_text(face = "bold", size = 11, color = "black"), axis.text = element_text(face = "bold", color = "black"))
save_plot_pair(p_prop_comp, file.path(out_dir, "Fig_CrossRegion_01_celltype_proportions_comparison"), width = 6.8, height = 5.2)

# 2. Cross-Region Expression Violins of Top Targets
top_genes <- c("DUSP1", "FOS", "JUN", "EGR1", "SORL1", "ADAMTS9")
genes_ctx <- intersect(top_genes, rownames(cortex))
genes_hpc <- intersect(top_genes, rownames(hpc))

ctx_sub <- subset(cortex, subset = cell_type %in% target_cts)
hpc_sub <- subset(hpc, subset = cell_type %in% target_cts)

p_vln_ctx <- VlnPlot(ctx_sub, features = genes_ctx, group.by = "cell_type", split.by = "Condition",
                     cols = c("Control" = "#4575b4", "AD" = "#d73027"), pt.size = 0, combine = FALSE)
plots_cr_ctx <- lapply(seq_along(p_vln_ctx), function(i) {
  p_vln_ctx[[i]] + ggtitle(genes_ctx[i]) +
    theme_pub(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold.italic", size = 12, hjust = 0.5, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 8, color = "black"),
      legend.title = element_text(face = "bold", size = 9, color = "black")
    )
})
p_vln_ctx_comb <- wrap_plots(plots_cr_ctx, ncol = 3) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
save_plot_pair(p_vln_ctx_comb, file.path(out_dir, "Fig_CrossRegion_02_cortex_expression_violins"), width = 12, height = 7)

p_vln_hpc <- VlnPlot(hpc_sub, features = genes_hpc, group.by = "cell_type", split.by = "Condition",
                     cols = c("Control" = "#4575b4", "AD" = "#d73027"), pt.size = 0, combine = FALSE)
plots_cr_hpc <- lapply(seq_along(p_vln_hpc), function(i) {
  p_vln_hpc[[i]] + ggtitle(genes_hpc[i]) +
    theme_pub(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold.italic", size = 12, hjust = 0.5, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 8, color = "black"),
      legend.title = element_text(face = "bold", size = 9, color = "black")
    )
})
p_vln_hpc_comb <- wrap_plots(plots_cr_hpc, ncol = 3) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
save_plot_pair(p_vln_hpc_comb, file.path(out_dir, "Fig_CrossRegion_03_hippocampus_expression_violins"), width = 12, height = 7)

# 3. Cross-Region DEG Concordance (Correlation between Cortex log2FC and Hippocampus log2FC)
deg_pairs <- list(
  "Astrocyte"         = c("tables/DEGs_Cortex_Astrocyte.csv", "tables/DEGs_Hippocampus_Astrocyte.csv"),
  "Endothelial"       = c("tables/DEGs_Cortex_Endothelial.csv", "tables/DEGs_Hippocampus_Endothelial.csv"),
  "Inhibitory_neuron" = c("tables/DEGs_Cortex_Inhibitory_neuron.csv", "tables/DEGs_Hippocampus_Inhibitory_neuron.csv"),
  "Microglia"         = c("tables/DEGs_Cortex_Microglia.csv", "tables/DEGs_Hippocampus_Microglia.csv"),
  "Oligodendrocyte"   = c("tables/DEGs_Cortex_Oligodendrocyte.csv", "tables/DEGs_Hippocampus_Oligodendrocyte.csv"),
  "Pericyte"          = c("tables/DEGs_Cortex_Pericyte.csv", "tables/DEGs_Hippocampus_Pericyte.csv")
)

deg_concordance_list <- list()
for (ct_name in names(deg_pairs)) {
  f_ctx <- deg_pairs[[ct_name]][1]
  f_hpc <- deg_pairs[[ct_name]][2]
  if (file.exists(f_ctx) && file.exists(f_hpc)) {
    d_ctx <- read_csv(f_ctx, show_col_types = FALSE) %>% select(gene, avg_log2FC_ctx = avg_log2FC, padj_ctx = p_val_adj)
    d_hpc <- read_csv(f_hpc, show_col_types = FALSE) %>% select(gene, avg_log2FC_hpc = avg_log2FC, padj_hpc = p_val_adj)
    m <- inner_join(d_ctx, d_hpc, by = "gene")
    m$Cell_Type <- ct_name
    deg_concordance_list[[ct_name]] <- m
  }
}

if (length(deg_concordance_list) > 0) {
  df_concordance <- bind_rows(deg_concordance_list)
  p_concordance <- ggplot(df_concordance, aes(x = avg_log2FC_ctx, y = avg_log2FC_hpc, color = Cell_Type)) +
    geom_point(alpha = 0.5, size = 1.2) +
    geom_smooth(method = "lm", color = "black", linewidth = 0.6, se = FALSE) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey60") +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey60") +
    stat_cor(method = "pearson", label.x.npc = "left", label.y.npc = "top", size = 3.8, color = "black") +
    facet_wrap(~Cell_Type, scales = "free") +
    labs(x = "Cortex log2(Fold Change)", y = "Hippocampus log2(Fold Change)", color = "Cell Type") +
    theme_pub(base_size = 9) +
    theme(strip.text = element_text(face = "bold", size = 9, color = "black"), legend.position = "none")
  save_plot_pair(p_concordance, file.path(out_dir, "Fig_CrossRegion_04_deg_concordance_scatter"), width = 10, height = 7.5)
}

# 4. Synchronize In Silico Target Validation Figures into Cross_Region_Comparison
fig5_panels <- c(
  "Fig5a_target_cross_species_conservation",
  "Fig5b_target_druggability_and_surfaceome",
  "Fig5c_target_ppi_mechanistic_network"
)
for (p_base in fig5_panels) {
  for (ext in c(".pdf", ".png")) {
    src_f <- file.path("figures", "Cortex", "Fig5_Target_Validation", paste0(p_base, ext))
    if (!file.exists(src_f)) {
      src_f <- file.path("figures", "Hippocampus", "Fig5_Target_Validation", paste0(p_base, ext))
    }
    if (file.exists(src_f)) {
      file.copy(src_f, file.path(out_dir, paste0(p_base, ext)), overwrite = TRUE)
    }
  }
  cat(sprintf("  Ready: %s (.pdf & .png)\n", p_base))
}

cat("\n======================================================================\n")
cat("STEP 08: CROSS-REGION COMPARISON SUITE COMPLETE!\n")
cat("======================================================================\n")
