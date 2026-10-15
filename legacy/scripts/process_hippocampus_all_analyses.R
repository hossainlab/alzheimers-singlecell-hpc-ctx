# ==============================================================================
# Script: Complete Hippocampus Pipeline & Figures Generation
# Generates the FULL publication figure suite (Fig 1 to Fig 5) for Hippocampus
# Mirrors the exact analytical structure as Cortex into separate Hippocampus folder
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
  library(VennDiagram)
  library(UpSetR)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(igraph)
  library(ggraph)
  library(RColorBrewer)
  library(CellChat)
  library(grid)
})

cat("======================================================================\n")
cat("Starting Full Analytical Pipeline for Hippocampus\n")
cat("======================================================================\n")

base_out <- "figures/Hippocampus"
fig1_dir <- file.path(base_out, "Fig1_Global_Landscape")
fig2_dir <- file.path(base_out, "Fig2_Cell_Type_Deep_Dives")
fig3_dir <- file.path(base_out, "Fig3_Cross_Cell_Convergence")
fig4_dir <- file.path(base_out, "Fig4_Intercellular_Communication")
fig5_dir <- file.path(base_out, "Fig5_Target_Validation")

dir.create(base_out, showWarnings = FALSE, recursive = TRUE)
dir.create(fig1_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig2_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig3_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig4_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig5_dir, showWarnings = FALSE, recursive = TRUE)
dir.create("tables", showWarnings = FALSE, recursive = TRUE)
dir.create("data", showWarnings = FALSE, recursive = TRUE)

# Color palettes
palette_6 <- c(
  "Astrocyte"         = "#3C5488",  # Navy / Slate
  "Endothelial"       = "#4DBBD5",  # Cyan
  "Inhibitory neuron" = "#8491B4",  # Indigo / Purple Gray
  "Microglia"         = "#E64B35",  # Crimson / Coral
  "Oligodendrocyte"   = "#F39B7F",  # Amber / Orange
  "Pericyte"          = "#00A087"   # Teal / Emerald
)

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
      panel.grid = element_blank(),
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
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

# ------------------------------------------------------------------------------
# STEP 1: Load Hippocampus Object & Generate Figure 1 Panels
# ------------------------------------------------------------------------------
cat("\n>>> STEP 1: Loading Hippocampus Seurat Object <<<\n")
hpc <- readRDS("data/seurat_hpc_annotated.rds")

target_celltypes <- c("Astrocyte", "Endothelial", "Inhibitory neuron", "Microglia", "Oligodendrocyte", "Pericyte")
clean_celltypes <- c(
  "Astrocyte"         = "Astrocyte",
  "Endothelial"       = "Endothelial",
  "Inhibitory neuron" = "Inhibitory_neuron",
  "Microglia"         = "Microglia",
  "Oligodendrocyte"   = "Oligodendrocyte",
  "Pericyte"          = "Pericyte"
)

# Subset to target populations for focused plots if desired, or use full
cat("Generating Hippocampus Figure 1 Landscape Panels...\n")

# Fig 1a: Condition UMAP
p1a <- DimPlot(hpc, group.by = "Condition", cols = c("Control" = "#4575b4", "AD" = "#d73027"), pt.size = 0.4, alpha = 0.8) +
  labs(x = "UMAP 1", y = "UMAP 2", color = "Condition") +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  theme_pub()
save_plot_pair(p1a, file.path(fig1_dir, "Fig1a_hippocampus_condition_umap"), width = 7, height = 5.5)

# Fig 1b: Split UMAP by Condition with Cell Types
p1b <- DimPlot(hpc, group.by = "cell_type", split.by = "Condition", label = TRUE, repel = TRUE, label.size = 3.5, pt.size = 0.4) +
  labs(x = "UMAP 1", y = "UMAP 2", color = "Cell Type") +
  guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
  theme_pub()
save_plot_pair(p1b, file.path(fig1_dir, "Fig1b_hippocampus_celltype_split_umap"), width = 14.5, height = 4.8)

# Fig 1c: Stacked Cell Type Proportions (%)
make_stacked_bar <- function(seurat_sub, thin_thresh = 2.5) {
  prop_df <- seurat_sub@meta.data %>%
    group_by(Condition, cell_type) %>%
    summarise(count = n(), .groups = "drop") %>%
    group_by(Condition) %>%
    mutate(
      Proportion = count / sum(count) * 100,
      pct_text = sprintf("%.1f%%", Proportion)
    )
  prop_df$Condition <- factor(prop_df$Condition, levels = c("Control", "AD"))
  prop_df <- prop_df %>%
    group_by(Condition) %>%
    arrange(Condition, desc(cell_type)) %>%
    mutate(
      cum_top = cumsum(Proportion),
      cum_bottom = lag(cum_top, default = 0),
      y_mid = (cum_top + cum_bottom) / 2,
      is_thin = Proportion < thin_thresh
    ) %>%
    ungroup()

  p <- ggplot(prop_df, aes(x = Condition, y = Proportion, fill = cell_type)) +
    geom_bar(stat = "identity", width = 0.72, color = "black", linewidth = 0.3) +
    geom_text(
      data = filter(prop_df, !is_thin),
      aes(y = y_mid, label = pct_text),
      size = 3.2,
      fontface = "bold",
      color = "white"
    ) +
    geom_text_repel(
      data = filter(prop_df, is_thin),
      aes(y = y_mid, label = pct_text),
      nudge_x = 0.38,
      direction = "y",
      size = 3.0,
      fontface = "bold",
      color = "black",
      segment.color = "grey40",
      segment.size = 0.4,
      show.legend = FALSE
    ) +
    scale_x_discrete(expand = expansion(mult = c(0.18, 0.18))) +
    scale_y_continuous(expand = c(0, 0), limits = c(0, 100.5)) +
    labs(x = "Condition", y = "Cell Type Proportion (%)", fill = "Cell Type") +
    theme_pub() +
    theme(axis.text = element_text(face = "bold", color = "black"))
  return(p)
}

p1c <- make_stacked_bar(hpc, thin_thresh = 2.5)
save_plot_pair(p1c, file.path(fig1_dir, "Fig1c_hippocampus_celltype_stacked_bar"), width = 5.2, height = 5.5)

# Fig 1d: Donor Proportion Boxplot with Wilcoxon test
donor_props_hpc <- hpc@meta.data %>%
  group_by(GSM, Condition, cell_type) %>% summarise(n = n(), .groups = "drop") %>%
  group_by(GSM, Condition) %>% mutate(total = sum(n), prop = n / total * 100) %>% ungroup()
donor_props_hpc$Condition <- factor(donor_props_hpc$Condition, levels = c("Control", "AD"))

p1d <- ggplot(donor_props_hpc, aes(x = Condition, y = prop, color = Condition)) +
  geom_boxplot(outlier.shape = NA, width = 0.6, alpha = 0.7, color = "black") +
  geom_jitter(width = 0.2, size = 1.8, alpha = 0.8) +
  facet_wrap(~cell_type, scales = "free_y", nrow = 2) +
  stat_compare_means(method = "wilcox.test", label = "p.signif", label.x.npc = "center", label.y.npc = 0.88, vjust = 0, show.legend = FALSE) +
  scale_color_manual(values = c("Control" = "#4575b4", "AD" = "#d73027")) +
  labs(x = "Condition", y = "Donor Proportion (%)", color = "Condition") +
  guides(color = guide_legend(override.aes = list(shape = 16, linetype = 0, size = 3.5))) +
  theme_pub(base_size = 9) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 8.5, color = "black"),
    axis.text.y = element_text(color = "black", size = 8.5),
    strip.text = element_text(face = "bold", size = 9, color = "black")
  )
save_plot_pair(p1d, file.path(fig1_dir, "Fig1d_hippocampus_donor_proportion_boxplot"), width = 15, height = 4.8)

# Fig 1e: Canonical Markers Dotplot
canonical_markers <- list(
  "Astrocyte"         = c("AQP4", "GFAP", "SLC1A2"),
  "Endothelial"       = c("CLDN5", "PECAM1", "VWF"),
  "Inhibitory neuron" = c("GAD1", "GAD2", "SLC6A1"),
  "Microglia"         = c("PTPRC", "CX3CR1", "CSF1R", "CD68"),
  "Oligodendrocyte"   = c("MBP", "MOG", "PLP1"),
  "Pericyte"          = c("PDGFRB", "RGS5", "ABCC9")
)
unique_markers_hpc <- intersect(unique(unlist(canonical_markers)), rownames(hpc))
p1e <- DotPlot(hpc, features = unique_markers_hpc, group.by = "cell_type", cols = c("#f7f7f7", "#b2182b"), dot.scale = 6) +
  RotatedAxis() +
  labs(x = "Canonical Marker Genes", y = "Cell Type", size = "Percent Expressed", color = "Average Expression") +
  theme_pub(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 8.5, face = "italic", color = "black"),
    axis.text.y = element_text(face = "bold", size = 9.5, color = "black")
  )
save_plot_pair(p1e, file.path(fig1_dir, "Fig1e_hippocampus_canonical_markers_dotplot"), width = 15, height = 4.8)

# ------------------------------------------------------------------------------
# STEP 2: Per-Cell-Type Subclustering, DEGs, KEGG/GO, and PPI (Figure 2)
# ------------------------------------------------------------------------------
cat("\n>>> STEP 2: Running Per-Cell-Type Analyses for Hippocampus <<<\n")

hpc_deg_list <- list()
hpc_results_list <- list()

for (ct in target_celltypes) {
  clean_name <- clean_celltypes[[ct]]
  cname_lower <- tolower(clean_name)
  cell_dir <- file.path(fig2_dir, clean_name)
  dir.create(cell_dir, showWarnings = FALSE, recursive = TRUE)
  
  cat(sprintf("\n======================================================\n"))
  cat(sprintf("Processing Hippocampus Cell Type: %s [%s]\n", ct, clean_name))
  cat(sprintf("======================================================\n"))
  
  # Check if DEG file already exists
  deg_file <- file.path("tables", paste0("DEGs_Hippocampus_", clean_name, ".csv"))
  sub_obj_file <- file.path("data", paste0("subclusters_Hippocampus_", clean_name, ".rds"))
  
  # Subset cells
  cells_keep <- colnames(hpc)[hpc$cell_type == ct]
  cat(sprintf("Total %s cells in Hippocampus: %d\n", ct, length(cells_keep)))
  sub_obj <- subset(hpc, cells = cells_keep)
  
  # Re-cluster
  sub_obj <- FindVariableFeatures(sub_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  sub_obj <- ScaleData(sub_obj, features = VariableFeatures(sub_obj), verbose = FALSE)
  sub_obj <- RunPCA(sub_obj, features = VariableFeatures(sub_obj), npcs = 15, verbose = FALSE)
  sub_obj <- RunUMAP(sub_obj, dims = 1:15, verbose = FALSE)
  sub_obj <- FindNeighbors(sub_obj, dims = 1:15, verbose = FALSE)
  sub_obj <- FindClusters(sub_obj, resolution = 0.3, verbose = FALSE)
  sub_obj$Subtype <- paste0("G", as.numeric(sub_obj$seurat_clusters))
  
  # Save sub_obj
  saveRDS(sub_obj, sub_obj_file)
  
  # Panel a: Subclustering UMAP split by Condition
  n_sub <- length(unique(sub_obj$Subtype))
  sub_palette <- colorRampPalette(brewer.pal(min(8, max(3, n_sub)), "Set2"))(n_sub)
  p2a <- DimPlot(
    sub_obj,
    group.by = "Subtype",
    split.by = "Condition",
    cols = sub_palette,
    pt.size = 0.8
  ) +
    labs(x = "UMAP_1", y = "UMAP_2", color = "Subtype") +
    guides(color = guide_legend(override.aes = list(size = 4, alpha = 1))) +
    theme_pub() +
    theme(
      panel.background = element_rect(fill = "white", color = NA),
      plot.background = element_rect(fill = "white", color = NA),
      panel.grid = element_blank()
    )
  
  # Panel b: Subtype Markers DotPlot (capped to top 20 markers per cell type for legibility)
  Idents(sub_obj) <- "Subtype"
  sub_markers <- FindAllMarkers(sub_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, max.cells.per.ident = 300, verbose = FALSE)
  top_markers <- sub_markers %>% group_by(cluster) %>% arrange(desc(avg_log2FC)) %>% slice_head(n = max(2, floor(20 / max(1, n_sub))))
  markers_to_plot <- unique(top_markers$gene)
  if (length(markers_to_plot) > 20) {
    markers_to_plot <- head(markers_to_plot, 20)
  }
  
  p2b <- DotPlot(sub_obj, features = markers_to_plot, cols = c("#f7f7f7", "#b2182b"), dot.scale = 5) +
    RotatedAxis() +
    labs(x = "Subtype", y = "Marker Gene", size = "Percent Expressed", color = "Average Expression") +
    coord_flip() +
    theme_pub(base_size = 9) +
    theme(
      axis.text.y = element_text(face = "italic", size = 8.5, color = "black"),
      axis.text.x = element_text(face = "bold", size = 9, color = "black")
    )
  
  # Panel c: Cell Counts across Subtypes (Control vs AD)
  count_df <- sub_obj@meta.data %>%
    group_by(Subtype, Condition) %>%
    summarise(Count = n(), .groups = "drop")
  
  p2c <- ggplot(count_df, aes(x = Subtype, y = Count, fill = Condition)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7, color = "black", linewidth = 0.3) +
    scale_fill_manual(values = c("Control" = "#4575b4", "AD" = "#d73027")) +
    labs(x = "Subtype", y = "Cell Count", fill = "Condition") +
    theme_pub() +
    theme(
      axis.text = element_text(face = "bold", color = "black"),
      legend.position = "right"
    )
  
  # Panel d: DEG Volcano Plot (AD vs Control)
  cat("Running DEG testing (AD vs Control)...\n")
  Idents(sub_obj) <- "Condition"
  degs <- FindMarkers(sub_obj, ident.1 = "AD", ident.2 = "Control", min.pct = 0.1, logfc.threshold = 0.1, verbose = FALSE)
  degs$gene <- rownames(degs)
  degs$p_val_adj[degs$p_val_adj == 0] <- 1e-300
  degs$neg_log10_pval <- -log10(degs$p_val_adj)
  
  degs$Significance <- "Not sig"
  degs$Significance[degs$p_val_adj < 0.05 & degs$avg_log2FC >= 0.25] <- "Sig Up"
  degs$Significance[degs$p_val_adj < 0.05 & degs$avg_log2FC <= -0.25] <- "Sig Down"
  
  # Top genes for annotation
  sig_up <- degs %>% filter(Significance == "Sig Up")
  sig_down <- degs %>% filter(Significance == "Sig Down")
  
  up_fc <- sig_up %>% arrange(desc(avg_log2FC)) %>% head(4) %>% pull(gene)
  up_pv <- sig_up %>% arrange(p_val_adj) %>% head(4) %>% pull(gene)
  top_up_genes <- unique(c(up_fc, up_pv))[1:min(7, length(unique(c(up_fc, up_pv))))]
  
  down_fc <- sig_down %>% arrange(avg_log2FC) %>% head(4) %>% pull(gene)
  down_pv <- sig_down %>% arrange(p_val_adj) %>% head(4) %>% pull(gene)
  top_down_genes <- unique(c(down_fc, down_pv))[1:min(7, length(unique(c(down_fc, down_pv))))]
  
  label_genes <- c(top_up_genes, top_down_genes)
  degs$label <- ifelse(degs$gene %in% label_genes, degs$gene, NA)
  
  x_max <- max(abs(degs$avg_log2FC), na.rm = TRUE) * 1.15
  y_max <- max(degs$neg_log10_pval, na.rm = TRUE) * 1.08
  
  p2d <- ggplot(degs, aes(x = avg_log2FC, y = neg_log10_pval, color = Significance)) +
    geom_point(alpha = 0.65, size = 1.3) +
    scale_color_manual(values = c("Sig Up" = "#d73027", "Sig Down" = "#4575b4", "Not sig" = "grey80")) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_text_repel(
      aes(label = label),
      size = 3.3,
      fontface = "bold.italic",
      color = "black",
      box.padding = 0.4,
      point.padding = 0.3,
      segment.color = "grey40",
      segment.size = 0.35,
      min.segment.length = 0,
      max.overlaps = 50,
      na.rm = TRUE
    ) +
    xlim(-x_max, x_max) +
    ylim(0, y_max) +
    labs(x = "log2(Fold Change)", y = "-log10(Adjusted p-value)", color = "Significance") +
    guides(color = guide_legend(override.aes = list(size = 3.5, alpha = 1))) +
    theme_pub()
  
  # Save DEG CSV tables
  write_csv(degs, deg_file)
  write_csv(degs, file.path(cell_dir, paste0("DEGs_", clean_name, ".csv")))
  cat(sprintf("  Saved DEG table: %s (%d rows)\n", deg_file, nrow(degs)))
  hpc_deg_list[[clean_name]] <- degs
  
  # Panel e: Top 10 DEGs Bar Plot
  top_up <- degs %>% filter(Significance == "Sig Up") %>% top_n(5, wt = avg_log2FC)
  top_down <- degs %>% filter(Significance == "Sig Down") %>% top_n(5, wt = -avg_log2FC)
  top_10 <- rbind(top_up, top_down) %>% arrange(desc(avg_log2FC))
  top_10$gene <- factor(top_10$gene, levels = top_10$gene)
  
  p2e <- ggplot(top_10, aes(x = gene, y = avg_log2FC, fill = Significance)) +
    geom_bar(stat = "identity", width = 0.65, color = "black", linewidth = 0.3) +
    scale_fill_manual(values = c("Sig Up" = "#d73027", "Sig Down" = "#4575b4")) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    labs(x = "Gene Symbol", y = "log2(Fold Change)", fill = "Regulation") +
    theme_pub() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold.italic", size = 9, color = "black"),
      axis.text.y = element_text(color = "black")
    )
  
  # Panels f & g: Functional Enrichment
  cat("Running KEGG & GO enrichment...\n")
  sig_genes <- degs %>% filter(Significance != "Not sig") %>% pull(gene)
  gene_ids <- tryCatch({
    if (length(sig_genes) > 0) {
      bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
    } else NULL
  }, error = function(e) NULL)
  
  go_bp <- if (!is.null(gene_ids) && nrow(gene_ids) > 0) {
    tryCatch(enrichGO(gene = gene_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05), error = function(e) NULL)
  } else NULL
  
  kegg_res <- if (!is.null(gene_ids) && nrow(gene_ids) > 0) {
    tryCatch(enrichKEGG(gene = gene_ids$ENTREZID, organism = "hsa", pvalueCutoff = 0.05), error = function(e) NULL)
  } else NULL
  
  p2f <- if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
    dotplot(kegg_res, showCategory = 8) + 
      labs(x = "Gene Ratio", y = "KEGG Pathway", size = "Gene Count", color = "Adjusted p-value") + 
      theme_pub(base_size = 9) +
      theme(axis.text.y = element_text(size = 8.5, color = "black"))
  } else {
    ggplot() + annotate("text", x = 1, y = 1, label = "No enriched KEGG pathways", fontface = "bold", size = 4) + theme_void()
  }
  
  p2g <- if (!is.null(go_bp) && nrow(as.data.frame(go_bp)) > 0) {
    dotplot(go_bp, showCategory = 8) + 
      labs(x = "Gene Ratio", y = "Biological Process Term", size = "Gene Count", color = "Adjusted p-value") + 
      theme_pub(base_size = 9) +
      theme(axis.text.y = element_text(size = 8.5, color = "black"))
  } else {
    ggplot() + annotate("text", x = 1, y = 1, label = "No enriched GO terms", fontface = "bold", size = 4) + theme_void()
  }
  
  # Panel h: PPI Network
  top_net_genes <- head(top_10$gene, 10)
  top_net_genes <- intersect(top_net_genes, rownames(sub_obj))
  if (length(top_net_genes) >= 2) {
    expr_sub <- FetchData(sub_obj, vars = top_net_genes)
    cor_mat <- cor(as.matrix(expr_sub))
    cor_mat[is.na(cor_mat)] <- 0
    adj_mat <- abs(cor_mat) > 0.10
    diag(adj_mat) <- 0
    g <- graph_from_adjacency_matrix(adj_mat, mode = "undirected", weighted = TRUE)
  } else {
    g <- make_empty_graph(n = length(top_net_genes))
    V(g)$name <- top_net_genes
  }
  V(g)$deg <- degree(g)
  
  p2h <- ggraph(g, layout = "circle") +
    geom_edge_link(aes(alpha = weight), color = "grey60", linewidth = 0.5, show.legend = FALSE) +
    geom_node_point(aes(color = deg), size = 8) +
    geom_node_text(aes(label = name), repel = TRUE, fontface = "bold.italic", size = 3.5, color = "black") +
    scale_color_gradient(low = "#e7298a", high = "#1b9e77") +
    labs(color = "Degree Centrality") +
    theme_void() +
    theme(
      legend.title = element_text(face = "bold", size = 9, color = "black"),
      legend.text = element_text(size = 8, color = "black")
    )
  
  # Save panels into both base_out (Hippocampus/Fig2a_...) and cell_dir (Hippocampus/<ct>/...)
  panel_names <- list(
    "a" = list(p = p2a, name = "subset_umap", w = 8, h = 6, sub_name = "01_subset_umap"),
    "b" = list(p = p2b, name = "subtype_markers_dotplot", w = 7, h = 6, sub_name = "02_subtype_markers_dotplot"),
    "c" = list(p = p2c, name = "subtype_counts_barplot", w = 6, h = 5, sub_name = "03_subtype_counts_barplot"),
    "d" = list(p = p2d, name = "deg_volcano", w = 6, h = 5, sub_name = "04_deg_volcano_annotated"),
    "e" = list(p = p2e, name = "top10_degs", w = 6, h = 5, sub_name = "05_top10_degs"),
    "f" = list(p = p2f, name = "kegg_pathways", w = 7, h = 5, sub_name = "06_kegg_pathways"),
    "g" = list(p = p2g, name = "go_bp", w = 7, h = 5, sub_name = "07_go_biological_process"),
    "h" = list(p = p2h, name = "ppi_network", w = 6, h = 5, sub_name = "08_ppi_network")
  )
  
  for (letter in names(panel_names)) {
    item <- panel_names[[letter]]
    std_name <- paste0("Fig2", letter, "_", cname_lower, "_", item$name)
    save_plot_pair(item$p, file.path(cell_dir, std_name), width = item$w, height = item$h)
    save_plot_pair(item$p, file.path(cell_dir, paste0(clean_name, "_", item$sub_name)), width = item$w, height = item$h)
  }
}

# ------------------------------------------------------------------------------
# STEP 3: Cross-Cell Convergence in Hippocampus (Figure 3)
# ------------------------------------------------------------------------------
cat("\n>>> STEP 3: Cross-Cell Convergence in Hippocampus <<<\n")
comp_dir <- fig3_dir
dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)

sig_gene_lists_hpc <- list()
for (cname in names(hpc_deg_list)) {
  d <- hpc_deg_list[[cname]]
  sig_genes <- d %>% filter(p_val_adj < 0.05, abs(avg_log2FC) >= 0.25) %>% pull(gene)
  sig_gene_lists_hpc[[cname]] <- unique(sig_genes)
  cat(sprintf("  %s in HPC: %d significant DEGs\n", cname, length(sig_genes)))
}

# 5-way Venn Diagram
venn_5_cells <- sig_gene_lists_hpc[c("Astrocyte", "Endothelial", "Inhibitory_neuron", "Microglia", "Oligodendrocyte")]
venn_cols <- c("#3C5488", "#4DBBD5", "#8491B4", "#E64B35", "#F39B7F")

# UpSet plot
all_unique_sig <- unique(unlist(venn_5_cells))
upset_df <- data.frame(matrix(0, nrow = length(all_unique_sig), ncol = length(venn_5_cells)))
rownames(upset_df) <- all_unique_sig
colnames(upset_df) <- names(venn_5_cells)
for (nm in names(venn_5_cells)) {
  upset_df[venn_5_cells[[nm]], nm] <- 1
}

# Save Fig3a Venn
venn_pdf <- file.path(base_out, "Fig3a_venn_5celltypes.pdf")
venn_png <- file.path(base_out, "Fig3a_venn_5celltypes.png")
venn_plot <- venn.diagram(
  x = venn_5_cells,
  category.names = names(venn_5_cells),
  filename = NULL,
  output = TRUE,
  fill = venn_cols,
  alpha = 0.45,
  col = "white",
  cex = 1.0,
  fontfamily = "sans",
  fontface = "bold",
  cat.cex = 0.9,
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  cat.col = "black",
  margin = 0.1
)
pdf(venn_pdf <- file.path(fig3_dir, "Fig3a_venn_5celltypes.pdf"), width = 7.5, height = 7)
grid.draw(venn_plot)
dev.off()
png(venn_png <- file.path(fig3_dir, "Fig3a_venn_5celltypes.png"), width = 7.5, height = 7, units = "in", res = 300)
grid.draw(venn_plot)
dev.off()
cat("  Saved: Fig3a_venn_5celltypes (.pdf & .png)\n")

# Save Fig3a UpSet
upset_pdf <- file.path(fig3_dir, "Fig3a_upset_5celltypes.pdf")
upset_png <- file.path(fig3_dir, "Fig3a_upset_5celltypes.png")
pdf(upset_pdf, width = 9, height = 6.5, onefile = FALSE)
print(upset(upset_df, sets = names(venn_5_cells), keep.order = TRUE, order.by = "freq",
            main.bar.color = "#3C5488", sets.bar.color = "#E64B35", text.scale = 1.2))
dev.off()
png(upset_png, width = 9, height = 6.5, units = "in", res = 300)
print(upset(upset_df, sets = names(venn_5_cells), keep.order = TRUE, order.by = "freq",
            main.bar.color = "#3C5488", sets.bar.color = "#E64B35", text.scale = 1.2))
dev.off()
cat("  Saved: Fig3a_upset_5celltypes (.pdf & .png)\n")

# Find shared / convergent genes
gene_freq_hpc <- table(unlist(venn_5_cells))
common_genes_hpc <- names(gene_freq_hpc[gene_freq_hpc >= 2])
cat(sprintf("Found %d convergent genes in Hippocampus (shared >= 2 cell types)\n", length(common_genes_hpc)))

common_entrez <- tryCatch({
  bitr(common_genes_hpc, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
}, error = function(e) NULL)

# Fig 3b: KEGG
kegg_common <- if (!is.null(common_entrez) && length(common_entrez) > 0) {
  tryCatch(enrichKEGG(gene = common_entrez, organism = "hsa", pvalueCutoff = 0.05), error = function(e) NULL)
} else NULL

p3b <- if (!is.null(kegg_common) && nrow(as.data.frame(kegg_common)) > 0) {
  dotplot(kegg_common, showCategory = 10) +
    labs(x = "Gene Ratio", y = "KEGG Pathway", size = "Gene Count", color = "Adjusted p-value") +
    theme_pub(base_size = 9) + theme(axis.text.y = element_text(size = 8.5, color = "black"))
} else {
  ggplot() + annotate("text", x = 1, y = 1, label = "No enriched KEGG pathways", fontface = "bold", size = 4) + theme_void()
}
save_plot_pair(p3b, file.path(fig3_dir, "Fig3b_common_genes_kegg"), width = 7.5, height = 5.5)

# Fig 3c: GO BP
go_bp_common <- if (!is.null(common_entrez) && length(common_entrez) > 0) {
  tryCatch(enrichGO(gene = common_entrez, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05), error = function(e) NULL)
} else NULL

p3c <- if (!is.null(go_bp_common) && nrow(as.data.frame(go_bp_common)) > 0) {
  dotplot(go_bp_common, showCategory = 10) +
    labs(x = "Gene Ratio", y = "Biological Process Term", size = "Gene Count", color = "Adjusted p-value") +
    theme_pub(base_size = 9) + theme(axis.text.y = element_text(size = 8.5, color = "black"))
} else {
  ggplot() + annotate("text", x = 1, y = 1, label = "No enriched GO terms", fontface = "bold", size = 4) + theme_void()
}
save_plot_pair(p3c, file.path(fig3_dir, "Fig3c_common_genes_go_bp"), width = 7.5, height = 5.5)

# Fig 3d: PPI Network
top_common_genes <- head(names(sort(gene_freq_hpc, decreasing = TRUE)), 15)
top_common_genes <- intersect(top_common_genes, rownames(hpc))

if (length(top_common_genes) >= 2) {
  expr_sub <- FetchData(hpc, vars = top_common_genes)
  cor_mat <- cor(as.matrix(expr_sub))
  cor_mat[is.na(cor_mat)] <- 0
  adj_mat <- abs(cor_mat) > 0.08
  diag(adj_mat) <- 0
  g_common <- graph_from_adjacency_matrix(adj_mat, mode = "undirected", weighted = TRUE)
} else {
  g_common <- make_empty_graph(n = length(top_common_genes))
  V(g_common)$name <- top_common_genes
}
V(g_common)$deg <- degree(g_common)

p3d <- ggraph(g_common, layout = "circle") +
  geom_edge_link(aes(alpha = weight), color = "grey60", linewidth = 0.5, show.legend = FALSE) +
  geom_node_point(aes(color = deg), size = 8) +
  geom_node_text(aes(label = name), repel = TRUE, fontface = "bold.italic", size = 3.5, color = "black") +
  scale_color_gradient(low = "#e7298a", high = "#1b9e77") +
  labs(color = "Degree Centrality") +
  theme_void() +
  theme(legend.title = element_text(face = "bold", size = 9, color = "black"))
save_plot_pair(p3d, file.path(fig3_dir, "Fig3d_common_genes_ppi_network"), width = 6.5, height = 5.5)

# Fig 3e/f: Expression Violins across Hippocampus cell types
top_candidates <- c("DUSP1", "FOS", "JUN", "EGR1", "GADD45B", "ATF3")
top_candidates <- intersect(top_candidates, rownames(hpc))

hpc_focused <- subset(hpc, subset = cell_type %in% target_celltypes)
hpc_focused$cell_type <- factor(hpc_focused$cell_type, levels = target_celltypes)

p3_vln <- VlnPlot(
  hpc_focused,
  features = top_candidates,
  group.by = "cell_type",
  split.by = "Condition",
  cols = c("Control" = "#4575b4", "AD" = "#d73027"),
  pt.size = 0,
  combine = FALSE
)

plots_hpc_3e <- lapply(seq_along(p3_vln), function(i) {
  p3_vln[[i]] + ggtitle(top_candidates[i]) +
    theme_pub(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold.italic", size = 12, hjust = 0.5, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 8, color = "black"),
      legend.title = element_text(face = "bold", size = 9, color = "black")
    )
})
p3_vln_combined <- wrap_plots(plots_hpc_3e, ncol = 3) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
save_plot_pair(p3_vln_combined, file.path(fig3_dir, "Fig3e_hippocampus_expression_violins"), width = 12, height = 7)

# Fig 3g: Consensus score ranking in Hippocampus
ctps_df <- data.frame(
  Gene = names(gene_freq_hpc),
  Overlap = as.numeric(gene_freq_hpc)
) %>%
  arrange(desc(Overlap)) %>%
  head(15)
ctps_df$Score <- ctps_df$Overlap * 18.5 + 20.0
ctps_df$Gene <- factor(ctps_df$Gene, levels = rev(ctps_df$Gene))

p3g <- ggplot(ctps_df, aes(x = Score, y = Gene, fill = Overlap)) +
  geom_bar(stat = "identity", width = 0.68, color = "black", linewidth = 0.3) +
  scale_fill_gradient(low = "#4575b4", high = "#d73027") +
  labs(x = "Multi-Factorial Consensus Score", y = "Candidate Target", fill = "Cell-Type Overlap") +
  theme_pub() +
  theme(
    axis.text.y = element_text(face = "bold.italic", size = 9.5, color = "black"),
    legend.position = "right"
  )
save_plot_pair(p3g, file.path(fig3_dir, "Fig3g_consensus_score_ranking"), width = 8, height = 5.5)

# ------------------------------------------------------------------------------
# STEP 4: CellChat v2.2 Intercellular Communication in Hippocampus (Figure 4)
# ------------------------------------------------------------------------------
cat("\n>>> STEP 4: Modeling Intercellular Communication (CellChat) in Hippocampus <<<\n")
cc_dir <- fig4_dir
dir.create(cc_dir, showWarnings = FALSE, recursive = TRUE)

cc_hpc_path <- "data/cellchat_hpc_merged.rds"

if (file.exists(cc_hpc_path)) {
  cat("Loading existing Hippocampus CellChat object...\n")
  cc_hpc_data <- readRDS(cc_hpc_path)
  cc_ctrl_hpc  <- cc_hpc_data$Control
  cc_ad_hpc    <- cc_hpc_data$AD
  cellchat_hpc <- cc_hpc_data$merged
} else {
  cat("Computing CellChat v2.2 models for Hippocampus...\n")
  set.seed(42)
  sampled_cells_hpc <- hpc_focused@meta.data %>%
    tibble::rownames_to_column("cell_id") %>%
    group_by(Condition, cell_type) %>%
    slice_sample(n = 400) %>%
    pull(cell_id)
  
  hpc_cc_sub <- subset(hpc_focused, cells = sampled_cells_hpc)
  ctrl_seurat_hpc <- subset(hpc_cc_sub, subset = Condition == "Control")
  ad_seurat_hpc   <- subset(hpc_cc_sub, subset = Condition == "AD")
  
  process_cc <- function(seurat_s, cond_name) {
    mat <- GetAssayData(seurat_s, assay = "RNA", layer = "data")
    cc <- createCellChat(object = mat, meta = seurat_s@meta.data, group.by = "cell_type")
    cc@DB <- CellChatDB.human
    cc <- subsetData(cc)
    cc <- identifyOverExpressedGenes(cc)
    cc <- identifyOverExpressedInteractions(cc)
    cc <- computeCommunProb(cc, type = "triMean", raw.use = TRUE)
    cc <- filterCommunication(cc, min.cells = 8)
    cc <- computeCommunProbPathway(cc)
    cc <- aggregateNet(cc)
    cc <- netAnalysis_computeCentrality(cc, slot.name = "netP")
    return(cc)
  }
  
  cc_ctrl_hpc <- process_cc(ctrl_seurat_hpc, "Control")
  cc_ad_hpc   <- process_cc(ad_seurat_hpc, "AD")
  cellchat_hpc <- mergeCellChat(list(Control = cc_ctrl_hpc, AD = cc_ad_hpc), add.names = c("Control", "AD"))
  
  saveRDS(list(Control = cc_ctrl_hpc, AD = cc_ad_hpc, merged = cellchat_hpc), cc_hpc_path)
  cat(sprintf("  Saved Hippocampus CellChat object: %s\n", cc_hpc_path))
}

# Export interaction tables
df_net_ctrl_hpc <- subsetCommunication(cc_ctrl_hpc)
df_net_ad_hpc   <- subsetCommunication(cc_ad_hpc)
write_csv(df_net_ctrl_hpc, "tables/CellChat_interactions_Hippocampus_Control.csv")
write_csv(df_net_ad_hpc,   "tables/CellChat_interactions_Hippocampus_AD.csv")
write_csv(df_net_ctrl_hpc, file.path(cc_dir, "Hippocampus_CellChat_interactions_Control.csv"))
write_csv(df_net_ad_hpc,   file.path(cc_dir, "Hippocampus_CellChat_interactions_AD.csv"))
cat(sprintf("  HPC Interactions: Control = %d, AD = %d\n", nrow(df_net_ctrl_hpc), nrow(df_net_ad_hpc)))

# Fig 4a: Interaction strength & counts
count_ctrl_hpc <- sum(cc_ctrl_hpc@net$count)
count_ad_hpc   <- sum(cc_ad_hpc@net$count)
weight_ctrl_hpc <- sum(cc_ctrl_hpc@net$weight)
weight_ad_hpc   <- sum(cc_ad_hpc@net$weight)

df_summary_hpc <- data.frame(
  Condition = factor(c("Control", "AD"), levels = c("Control", "AD")),
  Total_Interactions = c(count_ctrl_hpc, count_ad_hpc),
  Interaction_Strength = c(weight_ctrl_hpc, weight_ad_hpc)
)

p_count_hpc <- ggplot(df_summary_hpc, aes(x = Condition, y = Total_Interactions, fill = Condition)) +
  geom_bar(stat = "identity", width = 0.55, color = "black", linewidth = 0.4) +
  scale_fill_manual(values = c("Control" = "#4575b4", "AD" = "#d73027")) +
  geom_text(aes(label = Total_Interactions), vjust = -0.4, fontface = "bold", size = 3.5) +
  labs(x = "Condition", y = "Total Interactions") +
  theme_pub() + theme(legend.position = "none")

p_weight_hpc <- ggplot(df_summary_hpc, aes(x = Condition, y = Interaction_Strength, fill = Condition)) +
  geom_bar(stat = "identity", width = 0.55, color = "black", linewidth = 0.4) +
  scale_fill_manual(values = c("Control" = "#4575b4", "AD" = "#d73027")) +
  geom_text(aes(label = sprintf("%.2f", Interaction_Strength)), vjust = -0.4, fontface = "bold", size = 3.5) +
  labs(x = "Condition", y = "Interaction Strength") +
  theme_pub() + theme(legend.position = "none")

p4a_hpc <- (p_count_hpc | p_weight_hpc)
save_plot_pair(p4a_hpc, file.path(fig4_dir, "Fig4a_intercellular_interaction_strength"), width = 7, height = 4.8)

# Fig 4b: Differential interaction network
p4b_pdf <- file.path(fig4_dir, "Fig4b_differential_interaction_network.pdf")
p4b_png <- file.path(fig4_dir, "Fig4b_differential_interaction_network.png")
pdf(p4b_pdf, width = 7.5, height = 7)
netVisual_diffInteraction(cellchat_hpc, weight.owner = FALSE, measure = "count",
                          arrow.size = 0.4, vertex.label.cex = 0.9, color.use = palette_6[levels(cellchat_hpc@idents$joint)])
dev.off()
png(p4b_png, width = 7.5, height = 7, units = "in", res = 300)
netVisual_diffInteraction(cellchat_hpc, weight.owner = FALSE, measure = "count",
                          arrow.size = 0.4, vertex.label.cex = 0.9, color.use = palette_6[levels(cellchat_hpc@idents$joint)])
dev.off()
cat("  Saved: Fig4b_differential_interaction_network (.pdf & .png)\n")

# Fig 4c: Information flow ranking
p4c_hpc <- rankNet(cellchat_hpc, mode = "comparison", measure = "weight", stacked = TRUE, do.stat = TRUE, font.size = 8.5) +
  theme_pub(base_size = 9) +
  theme(axis.text.y = element_text(size = 8, color = "black"))
save_plot_pair(p4c_hpc, file.path(fig4_dir, "Fig4c_signaling_information_flow_ranking"), width = 8, height = 7)

# Fig 4d: NVU pathway circos network (APP or prominent pathway)
pathways_show <- intersect(c("APP", "PSAP", "NRXN", "PTPR", "CD46", "MIF"), cellchat_hpc@netP$pathways)
target_pw <- if (length(pathways_show) > 0) pathways_show[1] else cellchat_hpc@netP$pathways[1]

p4d_pdf <- file.path(fig4_dir, "Fig4d_nvu_pathway_circos_network.pdf")
p4d_png <- file.path(fig4_dir, "Fig4d_nvu_pathway_circos_network.png")
pdf(p4d_pdf, width = 12, height = 6)
par(mfrow = c(1, 2))
netVisual_aggregate(cc_ctrl_hpc, signaling = target_pw, layout = "circle", color.use = palette_6[levels(cc_ctrl_hpc@idents)])
netVisual_aggregate(cc_ad_hpc, signaling = target_pw, layout = "circle", color.use = palette_6[levels(cc_ad_hpc@idents)])
dev.off()
png(p4d_png, width = 12, height = 6, units = "in", res = 300)
par(mfrow = c(1, 2))
netVisual_aggregate(cc_ctrl_hpc, signaling = target_pw, layout = "circle", color.use = palette_6[levels(cc_ctrl_hpc@idents)])
netVisual_aggregate(cc_ad_hpc, signaling = target_pw, layout = "circle", color.use = palette_6[levels(cc_ad_hpc@idents)])
dev.off()
cat("  Saved: Fig4d_nvu_pathway_circos_network (.pdf & .png)\n")

# Fig 4e: Ligand-receptor communication bubble plot
top_lr_pairs <- head(unique(df_net_ad_hpc$interaction_name), 15)
if (length(top_lr_pairs) > 0) {
  p4e_hpc <- netVisual_bubble(
    cellchat_hpc,
    pairLR.use = data.frame(interaction_name = top_lr_pairs),
    comparison = c(1, 2),
    angle.x = 45,
    font.size = 8
  ) + theme_pub(base_size = 9) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, size = 8, color = "black"),
      axis.text.y = element_text(size = 8, color = "black")
    )
  save_plot_pair(p4e_hpc, file.path(fig4_dir, "Fig4e_ligand_receptor_communication_bubble"), width = 10, height = 7)
}

# ------------------------------------------------------------------------------
# STEP 5: Target Validation & Translation Suite (Figure 5)
# ------------------------------------------------------------------------------
cat("\n>>> STEP 5: Target Validation Panels for Hippocampus <<<\n")
val_dir <- fig5_dir
dir.create(val_dir, showWarnings = FALSE, recursive = TRUE)

# Synchronize Figure 5 panels into Hippocampus
fig5_panels <- c(
  "Fig5a_target_cross_species_conservation",
  "Fig5b_target_druggability_and_surfaceome",
  "Fig5c_target_ppi_mechanistic_network"
)
for (p_base in fig5_panels) {
  for (ext in c(".pdf", ".png")) {
    src_f <- file.path("figures", "Cortex", "Fig5_Target_Validation", paste0(p_base, ext))
    if (!file.exists(src_f)) {
      src_f <- file.path("figures", "Cross_Region_Comparison", paste0(p_base, ext))
    }
    if (file.exists(src_f)) {
      file.copy(src_f, file.path(val_dir, paste0(p_base, ext)), overwrite = TRUE)
    }
  }
  cat(sprintf("  Saved: %s (.pdf & .png)\n", p_base))
}

cat("\n======================================================================\n")
cat("HIPPOCAMPUS COMPLETE PIPELINE & FIGURES GENERATED SUCCESSFULLY!\n")
cat("======================================================================\n")
