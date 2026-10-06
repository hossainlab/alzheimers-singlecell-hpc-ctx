# ==============================================================================
# Script 03b: Harmony Integration, Annotation, and Figure 1 for Hippocampus
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(cowplot)
  library(ggpubr)
  library(viridis)
  library(RColorBrewer)
})

cat("=== Step 3b: Loading Preprocessed Hippocampus Data ===\n")
hpc <- readRDS("data/seurat_hpc_qc.rds")
cat(sprintf("Loaded Hippocampus: %d genes across %d cells\n", nrow(hpc), ncol(hpc)))

cat("Joining layers for Seurat v5 compatibility...\n")
hpc[["RNA"]] <- JoinLayers(hpc[["RNA"]])

# Normalization & Variable Features
hpc <- NormalizeData(hpc, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
hpc <- FindVariableFeatures(hpc, selection.method = "vst", nfeatures = 2500, verbose = FALSE)
hpc <- ScaleData(hpc, features = VariableFeatures(hpc), verbose = FALSE)
hpc <- RunPCA(hpc, features = VariableFeatures(hpc), npcs = 30, verbose = FALSE)

cat("Running Harmony integration across donor/GSM batches...\n")
hpc <- RunHarmony(hpc, group.by.vars = "GSM", dims.use = 1:25, verbose = FALSE)

cat("Computing UMAP and graph-based clusters...\n")
hpc <- RunUMAP(hpc, reduction = "harmony", dims = 1:25, verbose = FALSE)
hpc <- FindNeighbors(hpc, reduction = "harmony", dims = 1:25, verbose = FALSE)
hpc <- FindClusters(hpc, resolution = 0.4, verbose = FALSE)

# Module scoring for canonical cell types
canonical_markers <- list(
  "Astrocyte"         = c("AQP4", "GFAP", "SLC1A2"),
  "Endothelial"       = c("CLDN5", "PECAM1", "VWF"),
  "Pericyte"          = c("PDGFRB", "RGS5", "ABCC9"),
  "Smooth Muscle"     = c("ACTA2", "TAGLN"),
  "Fibroblast"        = c("COL1A1", "COL1A2", "LUM"),
  "Microglia"         = c("PTPRC", "CX3CR1", "CSF1R", "CD68"),
  "T cells"           = c("CD3D", "TRAC", "CD4"),
  "Oligodendrocyte"   = c("MBP", "MOG", "PLP1"),
  "OPC"               = c("PDGFRA", "VCAN"),
  "Inhibitory neuron" = c("GAD1", "GAD2"),
  "Excitatory neuron" = c("SLC17A7", "CAMK2A"),
  "Ependymal"         = c("FOXJ1", "TTR")
)

# Score modules
cat("Assigning cell types based on canonical marker modules...\n")
for (ct in names(canonical_markers)) {
  genes <- intersect(canonical_markers[[ct]], rownames(hpc))
  if (length(genes) > 0) {
    hpc <- AddModuleScore(hpc, features = list(genes), name = paste0(gsub(" ", "_", ct), "_Score"), verbose = FALSE)
  }
}

# Determine top score per cluster
score_cols <- grep("_Score1$", colnames(hpc@meta.data), value = TRUE)
cluster_scores <- hpc@meta.data %>%
  group_by(seurat_clusters) %>%
  summarise(across(all_of(score_cols), mean))

cluster_annots <- character(nrow(cluster_scores))
names(cluster_annots) <- as.character(cluster_scores$seurat_clusters)

for (i in 1:nrow(cluster_scores)) {
  cl <- as.character(cluster_scores$seurat_clusters[i])
  scores <- as.numeric(cluster_scores[i, score_cols])
  best_idx <- which.max(scores)
  best_ct <- sub("_Score1$", "", score_cols[best_idx])
  best_ct <- gsub("_", " ", best_ct)
  cluster_annots[cl] <- best_ct
}

cat("Cluster assignments in Hippocampus:\n")
print(cluster_annots)

hpc$cell_type <- unname(cluster_annots[as.character(hpc$seurat_clusters)])
cat("Cell type counts in Hippocampus:\n")
print(table(hpc$cell_type, hpc$Condition))

# Save annotated object
saveRDS(hpc, "data/seurat_hpc_annotated.rds")
cat("Saved to data/seurat_hpc_annotated.rds\n")

# ------------------------------------------------------------------------------
# Generate Figure 1 Visualizations for Hippocampus
# ------------------------------------------------------------------------------
theme_pub <- function(base_size = 11) {
  theme_classic(base_size = base_size) +
    theme(
      plot.title = element_text(face = "bold", size = rel(1.1), hjust = 0),
      axis.title = element_text(face = "bold", size = rel(1.0)),
      axis.text = element_text(color = "black", size = rel(0.9)),
      legend.title = element_text(face = "bold", size = rel(0.9)),
      strip.background = element_rect(fill = "grey95", color = NA),
      strip.text = element_text(face = "bold")
    )
}

# 1. Panel A: UMAP colored by Condition
p1a <- DimPlot(
  hpc,
  group.by = "Condition",
  cols = c("Control" = "#4575b4", "AD" = "#d73027"),
  pt.size = 0.4,
  alpha = 0.8
) +
  labs(title = "Hippocampus: Cellular Landscape", x = "UMAP_1", y = "UMAP_2") +
  theme_pub()

# 2. Panel B: Split UMAP by Condition with Cell Types
p1b <- DimPlot(
  hpc,
  group.by = "cell_type",
  split.by = "Condition",
  label = TRUE,
  repel = TRUE,
  label.size = 3.2,
  pt.size = 0.4
) +
  labs(title = "Hippocampus: Control vs AD Cell Types", x = "UMAP_1", y = "UMAP_2") +
  theme_pub()

# 3. Panel I: Stacked Cell Type Proportions (%)
prop_df <- hpc@meta.data %>%
  group_by(Condition, cell_type) %>%
  summarise(count = n(), .groups = "drop") %>%
  group_by(Condition) %>%
  mutate(Proportion = count / sum(count) * 100)

p1i <- ggplot(prop_df, aes(x = Condition, y = Proportion, fill = cell_type)) +
  geom_bar(stat = "identity", width = 0.65, color = "black", size = 0.2) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 100.5)) +
  labs(title = "Cell Composition", x = "Condition", y = "Proportion (%)", fill = "Cell Type") +
  theme_pub()

# 4. Inset: Donor Proportion Boxplot with Wilcoxon test
donor_props <- hpc@meta.data %>%
  group_by(GSM, Condition, cell_type) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(GSM, Condition) %>%
  mutate(total = sum(n), prop = n / total * 100) %>%
  ungroup()

p_inset <- ggplot(donor_props, aes(x = Condition, y = prop, color = Condition)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.7, color = "black") +
  geom_jitter(width = 0.2, size = 1.8, alpha = 0.8) +
  facet_wrap(~cell_type, scales = "free_y", nrow = 2) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x.npc = "center", vjust = -0.5) +
  scale_color_manual(values = c("Control" = "#4575b4", "AD" = "#d73027")) +
  labs(title = "Donor-Level Proportions", x = "", y = "Donor proportion (%)") +
  theme_pub(base_size = 9) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

# 5. Marker DotPlot for Hippocampus
unique_markers <- unique(unlist(canonical_markers))
unique_markers <- intersect(unique_markers, rownames(hpc))

p_markers <- DotPlot(
  hpc,
  features = unique_markers,
  group.by = "cell_type",
  cols = c("lightgrey", "#b2182b"),
  dot.scale = 6
) +
  RotatedAxis() +
  labs(x = "Canonical Marker Genes", y = "Cell Type") +
  theme_pub(base_size = 11) +
  theme(axis.text.x = element_text(size = 9, face = "italic"))

# Export Figure 1 panel suite to figures/Hippocampus/Fig1_Global_Landscape/
out_fig1 <- "figures/Hippocampus/Fig1_Global_Landscape"
dir.create(out_fig1, showWarnings = FALSE, recursive = TRUE)

save_panel <- function(plt, base_name, w = 7, h = 6) {
  ggsave(file.path(out_fig1, paste0(base_name, ".png")), plot = plt, width = w, height = h, dpi = 300)
  ggsave(file.path(out_fig1, paste0(base_name, ".pdf")), plot = plt, width = w, height = h)
  cat(sprintf("  Saved %s (.pdf & .png)\n", base_name))
}

cat("\nSaving Hippocampus Figure 1 panels to:", out_fig1, "\n")
save_panel(p1a, "Fig1a_hippocampus_condition_umap", w = 7, h = 5.5)
save_panel(p1b, "Fig1b_hippocampus_celltype_split_umap", w = 10, h = 5.5)
save_panel(p1i, "Fig1c_hippocampus_celltype_stacked_bar", w = 6, h = 5.5)
save_panel(p_inset, "Fig1d_hippocampus_donor_proportion_boxplot", w = 8.5, h = 6.5)
save_panel(p_markers, "Fig1e_hippocampus_canonical_markers_dotplot", w = 12, h = 6)

cat("Step 3b complete!\n")
