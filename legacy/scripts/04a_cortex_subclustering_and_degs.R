# ==============================================================================
# Script 04a: Prefrontal Cortex Cell-Type Deep Dives & Subclustering (Figure 2)
# Lineages: Astrocyte, Endothelial, Inhibitory neuron, Microglia, Oligodendrocyte, Pericyte
# Generates 8 publication panels (a-h) + DEG CSV per lineage
# Outputs to: figures/Cortex/Fig2_Cell_Type_Deep_Dives/<cell_type>/
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(patchwork)
  library(cowplot)
  library(ggrepel)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(igraph)
  library(ggraph)
  library(RColorBrewer)
})

cat("======================================================================\n")
cat("Starting Step 04a: Cortex Lineage Deep Dives & Subclustering\n")
cat("======================================================================\n")

# Directories
cortex_rds <- "data/seurat_cortex_annotated.rds"
if (!file.exists(cortex_rds)) {
  stop("Annotated Cortex Seurat object not found at ", cortex_rds)
}

fig2_base <- "figures/Cortex/Fig2_Cell_Type_Deep_Dives"
dir.create(fig2_base, showWarnings = FALSE, recursive = TRUE)
dir.create("tables", showWarnings = FALSE, recursive = TRUE)
dir.create("data", showWarnings = FALSE, recursive = TRUE)

# Styling theme
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

cat("Loading annotated Cortex Seurat object...\n")
cortex <- readRDS(cortex_rds)

target_celltypes <- c("Astrocyte", "Endothelial", "Inhibitory neuron", "Microglia", "Oligodendrocyte", "Pericyte")
clean_celltypes <- c(
  "Astrocyte"         = "Astrocyte",
  "Endothelial"       = "Endothelial",
  "Inhibitory neuron" = "Inhibitory_neuron",
  "Microglia"         = "Microglia",
  "Oligodendrocyte"   = "Oligodendrocyte",
  "Pericyte"          = "Pericyte"
)

for (ct in target_celltypes) {
  clean_name <- clean_celltypes[[ct]]
  cname_lower <- tolower(clean_name)
  cell_dir <- file.path(fig2_base, clean_name)
  dir.create(cell_dir, showWarnings = FALSE, recursive = TRUE)
  
  cat(sprintf("\n--- Processing Cortex Lineage: %s [%s] ---\n", ct, clean_name))
  deg_file_tables <- file.path("tables", paste0("DEGs_Cortex_", clean_name, ".csv"))
  deg_file_folder <- file.path(cell_dir, paste0("DEGs_", clean_name, ".csv"))
  sub_obj_file <- file.path("data", paste0("subclusters_Cortex_", clean_name, ".rds"))
  
  # Subset lineage
  cells_keep <- colnames(cortex)[cortex$cell_type == ct]
  cat(sprintf("Total %s cells in Cortex: %d\n", ct, length(cells_keep)))
  sub_obj <- subset(cortex, cells = cells_keep)
  
  # Subclustering
  sub_obj <- FindVariableFeatures(sub_obj, selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  sub_obj <- ScaleData(sub_obj, features = VariableFeatures(sub_obj), verbose = FALSE)
  sub_obj <- RunPCA(sub_obj, features = VariableFeatures(sub_obj), npcs = 15, verbose = FALSE)
  sub_obj <- RunUMAP(sub_obj, dims = 1:15, verbose = FALSE)
  sub_obj <- FindNeighbors(sub_obj, dims = 1:15, verbose = FALSE)
  sub_obj <- FindClusters(sub_obj, resolution = 0.3, verbose = FALSE)
  sub_obj$Subtype <- paste0("G", as.numeric(sub_obj$seurat_clusters))
  
  saveRDS(sub_obj, sub_obj_file)
  cat("Saved subclusters to:", sub_obj_file, "\n")
  
  # Panel a: Subclustering UMAP split by Condition
  n_sub <- length(unique(sub_obj$Subtype))
  sub_palette <- colorRampPalette(brewer.pal(min(8, max(3, n_sub)), "Set2"))(n_sub)
  p2a <- DimPlot(sub_obj, group.by = "Subtype", split.by = "Condition", cols = sub_palette, pt.size = 0.8) +
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
    theme(axis.text = element_text(face = "bold", color = "black"), legend.position = "right")
  
  # Panel d: DEG Testing (AD vs Control) & Annotated Volcano Plot
  if (file.exists(deg_file_tables)) {
    cat("Loading existing DEG table:", deg_file_tables, "\n")
    degs <- read_csv(deg_file_tables, show_col_types = FALSE)
  } else {
    cat("Calculating Wilcoxon DEGs (AD vs Control)...\n")
    Idents(sub_obj) <- "Condition"
    degs <- FindMarkers(sub_obj, ident.1 = "AD", ident.2 = "Control", min.pct = 0.1, logfc.threshold = 0.1, verbose = FALSE)
    degs$gene <- rownames(degs)
    write_csv(degs, deg_file_tables)
  }
  write_csv(degs, deg_file_folder)
  
  degs$p_val_adj[degs$p_val_adj == 0] <- 1e-300
  degs$neg_log10_pval <- -log10(degs$p_val_adj)
  degs$Significance <- "Not sig"
  degs$Significance[degs$p_val_adj < 0.05 & degs$avg_log2FC >= 0.25] <- "Sig Up"
  degs$Significance[degs$p_val_adj < 0.05 & degs$avg_log2FC <= -0.25] <- "Sig Down"
  
  sig_up <- degs %>% filter(Significance == "Sig Up")
  sig_down <- degs %>% filter(Significance == "Sig Down")
  up_top <- head(sig_up %>% arrange(desc(avg_log2FC)) %>% pull(gene), 4)
  down_top <- head(sig_down %>% arrange(avg_log2FC) %>% pull(gene), 4)
  label_genes <- unique(c(up_top, down_top))
  degs$label <- ifelse(degs$gene %in% label_genes, degs$gene, NA)
  
  x_max <- max(abs(degs$avg_log2FC), na.rm = TRUE) * 1.15
  y_max <- max(degs$neg_log10_pval, na.rm = TRUE) * 1.08
  
  p2d <- ggplot(degs, aes(x = avg_log2FC, y = neg_log10_pval, color = Significance)) +
    geom_point(alpha = 0.65, size = 1.3) +
    scale_color_manual(values = c("Sig Up" = "#d73027", "Sig Down" = "#4575b4", "Not sig" = "grey80")) +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey50", linewidth = 0.4) +
    geom_text_repel(
      aes(label = label), size = 3.2, fontface = "bold.italic", color = "black",
      box.padding = 0.4, point.padding = 0.3, segment.color = "grey40", segment.size = 0.35,
      min.segment.length = 0, max.overlaps = 50, na.rm = TRUE
    ) +
    xlim(-x_max, x_max) + ylim(0, y_max) +
    labs(x = "log2(Fold Change)", y = "-log10(Adjusted p-value)", color = "Significance") +
    theme_pub()
  
  # Panel e: Top 10 DEGs (5 Up, 5 Down)
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
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold.italic", size = 9, color = "black"))
  
  # Panels f & g: Functional Enrichment (KEGG & GO)
  sig_genes <- degs %>% filter(Significance != "Not sig") %>% pull(gene)
  gene_ids <- tryCatch({
    if (length(sig_genes) > 0) bitr(sig_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db) else NULL
  }, error = function(e) NULL)
  
  kegg_res <- if (!is.null(gene_ids) && nrow(gene_ids) > 0) {
    tryCatch(enrichKEGG(gene = gene_ids$ENTREZID, organism = "hsa", pvalueCutoff = 0.05), error = function(e) NULL)
  } else NULL
  
  p2f <- if (!is.null(kegg_res) && nrow(as.data.frame(kegg_res)) > 0) {
    dotplot(kegg_res, showCategory = 8) +
      labs(x = "Gene Ratio", y = "KEGG Pathway", size = "Gene Count", color = "Adjusted p-value") +
      theme_pub(base_size = 9) + theme(axis.text.y = element_text(size = 8.5, color = "black"))
  } else {
    ggplot() + annotate("text", x = 1, y = 1, label = "No enriched KEGG pathways", fontface = "bold", size = 4) + theme_void()
  }
  
  go_bp <- if (!is.null(gene_ids) && nrow(gene_ids) > 0) {
    tryCatch(enrichGO(gene = gene_ids$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05), error = function(e) NULL)
  } else NULL
  
  p2g <- if (!is.null(go_bp) && nrow(as.data.frame(go_bp)) > 0) {
    dotplot(go_bp, showCategory = 8) +
      labs(x = "Gene Ratio", y = "Biological Process", size = "Gene Count", color = "Adjusted p-value") +
      theme_pub(base_size = 9) + theme(axis.text.y = element_text(size = 8.5, color = "black"))
  } else {
    ggplot() + annotate("text", x = 1, y = 1, label = "No enriched GO terms", fontface = "bold", size = 4) + theme_void()
  }
  
  # Panel h: STRING PPI Network
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
    theme(legend.title = element_text(face = "bold", size = 9, color = "black"), legend.text = element_text(size = 8, color = "black"))
  
  # Save panels into cell_dir
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
    save_plot_pair(item$p, file.path(cell_dir, paste0(clean_name, "_", item$sub_name)), width = item$w, height = item$h)
  }
}

cat("\n======================================================================\n")
cat("STEP 04a: CORTEX LINEAGE DEEP DIVES COMPLETED SUCCESSFULLY!\n")
cat("======================================================================\n")
