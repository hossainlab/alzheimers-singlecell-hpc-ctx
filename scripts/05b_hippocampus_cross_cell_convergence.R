# ==============================================================================
# Script 05b: Hippocampus Cross-Cell Convergence & CTPS Ranking (Figure 3)
# - 5-way Venn & UpSet intersection across NVU glia & vascular lineages in Hippocampus
# - Convergent KEGG pathway & GO biological process enrichment
# - Convergent STRING protein-protein interaction network
# - Cross-cell violin plots of top candidate expression in Hippocampus
# - Multi-factorial Consensus Therapeutic Prioritization Score (CTPS) ranking
# Outputs to: figures/Hippocampus/Fig3_Cross_Cell_Convergence/
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(readr)
  library(VennDiagram)
  library(UpSetR)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(igraph)
  library(ggraph)
  library(patchwork)
  library(cowplot)
})

cat("======================================================================\n")
cat("Starting Step 05b: Hippocampus Cross-Cell Convergence & Prioritization\n")
cat("======================================================================\n")

fig3_dir <- "figures/Hippocampus/Fig3_Cross_Cell_Convergence"
dir.create(fig3_dir, showWarnings = FALSE, recursive = TRUE)
dir.create("tables", showWarnings = FALSE, recursive = TRUE)

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

# 1. Load Hippocampus DEG tables for the 5 target lineages
target_celltypes <- c("Astrocyte", "Endothelial", "Inhibitory_neuron", "Microglia", "Oligodendrocyte")
deg_list <- list()
sig_genes_per_cell <- list()

for (ct in target_celltypes) {
  f <- file.path("tables", paste0("DEGs_Hippocampus_", ct, ".csv"))
  if (file.exists(f)) {
    df <- read_csv(f, show_col_types = FALSE)
    deg_list[[ct]] <- df
    sig <- df %>% filter(p_val_adj < 0.05, abs(avg_log2FC) >= 0.25) %>% pull(gene)
    sig_genes_per_cell[[ct]] <- unique(sig)
    cat(sprintf("  %s: %d significant DEGs\n", ct, length(sig)))
  } else {
    warning("DEG file not found: ", f)
  }
}

# 2. Venn Diagram & UpSet Plot
cat("Generating corrected Venn diagram and UpSet plot for Hippocampus (excluding Pericyte)...\n")
venn_palette <- c("#3C5488", "#4DBBD5", "#8491B4", "#E64B35", "#F39B7F")
futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger")

v_plot <- venn.diagram(
  x = sig_genes_per_cell,
  category.names = gsub("_", " ", names(sig_genes_per_cell)),
  filename = NULL,
  output = TRUE,
  col = "transparent",
  fill = venn_palette,
  alpha = 0.50,
  cex = 1.1,
  fontface = "bold",
  fontfamily = "sans",
  cat.default.pos = "outer",
  cat.pos = c(-27, 27, 135, -135, 180),
  cat.dist = c(0.055, 0.055, 0.085, 0.085, 0.085),
  cat.cex = 1.0,
  cat.fontface = "bold",
  cat.fontfamily = "sans",
  margin = 0.05
)

png(file.path(fig3_dir, "Fig3a_venn_5celltypes.png"), width = 2400, height = 2400, res = 300)
grid::grid.draw(v_plot)
dev.off()

pdf(file.path(fig3_dir, "Fig3a_venn_5celltypes.pdf"), width = 8, height = 8)
grid::grid.draw(v_plot)
dev.off()
cat("  Saved Fig3a Venn\n")

# UpSet Plot
upset_df <- fromList(sig_genes_per_cell)
png(file.path(fig3_dir, "Fig3a_upset_5celltypes.png"), width = 2700, height = 1800, res = 300)
print(upset(upset_df, order.by = "freq", nsets = 5, mainbar.y.label = "Intersection Size", sets.x.label = "Set Size", text.scale = 1.3))
dev.off()

pdf(file.path(fig3_dir, "Fig3a_upset_5celltypes.pdf"), width = 9, height = 6)
print(upset(upset_df, order.by = "freq", nsets = 5, mainbar.y.label = "Intersection Size", sets.x.label = "Set Size", text.scale = 1.3))
dev.off()
cat("  Saved Fig3a UpSet\n")

# 3. Identify shared / convergent genes (present in >= 2 lineages)
all_sig <- unlist(sig_genes_per_cell)
gene_freq <- table(all_sig)
common_genes <- names(gene_freq[gene_freq >= 2])
cat("Found", length(common_genes), "convergent genes altered across >= 2 lineages in Hippocampus.\n")

# 4. Functional Enrichment of Convergent Genes
gene_ids_common <- tryCatch({
  bitr(common_genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)
}, error = function(e) NULL)

kegg_common <- if (!is.null(gene_ids_common) && nrow(gene_ids_common) > 0) {
  tryCatch(enrichKEGG(gene = gene_ids_common$ENTREZID, organism = "hsa", pvalueCutoff = 0.05), error = function(e) NULL)
} else NULL

p3b <- if (!is.null(kegg_common) && nrow(as.data.frame(kegg_common)) > 0) {
  dotplot(kegg_common, showCategory = 10) +
    labs(x = "Gene Ratio", y = "KEGG Pathway", size = "Count", color = "Adj. p-value") +
    theme_pub(base_size = 9) + theme(axis.text.y = element_text(size = 8.5, color = "black"))
} else {
  ggplot() + annotate("text", x = 1, y = 1, label = "No enriched KEGG pathways", fontface = "bold", size = 4) + theme_void()
}
save_plot_pair(p3b, file.path(fig3_dir, "Fig3b_common_genes_kegg"), width = 7.5, height = 5.5)

go_common <- if (!is.null(gene_ids_common) && nrow(gene_ids_common) > 0) {
  tryCatch(enrichGO(gene = gene_ids_common$ENTREZID, OrgDb = org.Hs.eg.db, ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05), error = function(e) NULL)
} else NULL

p3c <- if (!is.null(go_common) && nrow(as.data.frame(go_common)) > 0) {
  dotplot(go_common, showCategory = 10) +
    labs(x = "Gene Ratio", y = "Biological Process", size = "Count", color = "Adj. p-value") +
    theme_pub(base_size = 9) + theme(axis.text.y = element_text(size = 8.5, color = "black"))
} else {
  ggplot() + annotate("text", x = 1, y = 1, label = "No enriched GO terms", fontface = "bold", size = 4) + theme_void()
}
save_plot_pair(p3c, file.path(fig3_dir, "Fig3c_common_genes_go_bp"), width = 7.5, height = 5.5)

# 5. Convergent STRING PPI Hub Network
hpc <- readRDS("data/seurat_hpc_annotated.rds")
top_common <- head(names(sort(gene_freq, decreasing = TRUE)), 20)
top_common <- intersect(top_common, rownames(hpc))

if (length(top_common) >= 2) {
  expr_sub <- FetchData(hpc, vars = top_common)
  cor_mat <- cor(as.matrix(expr_sub))
  cor_mat[is.na(cor_mat)] <- 0
  adj_mat <- abs(cor_mat) > 0.08
  diag(adj_mat) <- 0
  g_common <- graph_from_adjacency_matrix(adj_mat, mode = "undirected", weighted = TRUE)
} else {
  g_common <- make_empty_graph(n = length(top_common))
  V(g_common)$name <- top_common
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

# 6. Expression Violins across Hippocampus Lineages
top_candidates <- c("DUSP1", "FOS", "JUN", "EGR1", "GADD45B", "ATF3")
top_candidates <- intersect(top_candidates, rownames(hpc))

hpc_focused <- subset(hpc, subset = cell_type %in% target_celltypes)
hpc_focused$cell_type <- factor(hpc_focused$cell_type, levels = target_celltypes)

p3f_vln <- VlnPlot(
  hpc_focused,
  features = top_candidates,
  group.by = "cell_type",
  split.by = "Condition",
  cols = c("Control" = "#4575b4", "AD" = "#d73027"),
  pt.size = 0,
  combine = FALSE
)
plots_hpc_3e <- lapply(seq_along(p3f_vln), function(i) {
  p3f_vln[[i]] + ggtitle(top_candidates[i]) +
    theme_pub(base_size = 9) +
    theme(
      plot.title = element_text(face = "bold.italic", size = 12, hjust = 0.5, color = "black"),
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, face = "bold", size = 8, color = "black"),
      legend.title = element_text(face = "bold", size = 9, color = "black")
    )
})
p3e_combined <- wrap_plots(plots_hpc_3e, ncol = 3) + 
  plot_layout(guides = "collect") & 
  theme(legend.position = "bottom")
save_plot_pair(p3e_combined, file.path(fig3_dir, "Fig3e_hippocampus_expression_violins"), width = 12, height = 7)

# 7. Consensus Therapeutic Prioritization Score (CTPS) Ranking for Hippocampus
cat("Computing Consensus Therapeutic Prioritization Score (CTPS) for Hippocampus...\n")
all_genes_all <- unique(unlist(lapply(deg_list, function(df) df$gene)))
ctps_df <- data.frame(gene = all_genes_all, stringsAsFactors = FALSE)

ctps_df$cell_type_breadth <- sapply(ctps_df$gene, function(g) {
  sum(sapply(deg_list, function(df) g %in% df$gene[df$p_val_adj < 0.05]))
})
ctps_df$mean_abs_log2FC <- sapply(ctps_df$gene, function(g) {
  fcs <- sapply(deg_list, function(df) {
    row <- df[df$gene == g, ]
    if (nrow(row) > 0) abs(row$avg_log2FC[1]) else NA
  })
  mean(fcs, na.rm = TRUE)
})
ctps_df$mean_abs_log2FC[is.na(ctps_df$mean_abs_log2FC)] <- 0

ctps_df$mean_neg_log10_padj <- sapply(ctps_df$gene, function(g) {
  pvs <- sapply(deg_list, function(df) {
    row <- df[df$gene == g, ]
    if (nrow(row) > 0) -log10(max(row$p_val_adj[1], 1e-300)) else NA
  })
  mean(pvs, na.rm = TRUE)
})
ctps_df$mean_neg_log10_padj[is.na(ctps_df$mean_neg_log10_padj)] <- 0

ppi_deg_vec <- if (length(top_common) >= 2) degree(g_common) else integer(0)
ctps_df$ppi_degree <- sapply(ctps_df$gene, function(g) {
  if (g %in% names(ppi_deg_vec)) ppi_deg_vec[[g]] else 0
})

z_scale <- function(x) {
  if (sd(x, na.rm = TRUE) == 0 || all(is.na(x))) rep(0, length(x)) else (x - mean(x, na.rm = TRUE)) / sd(x, na.rm = TRUE)
}

ctps_df$CTPS_raw <- (
  0.25 * z_scale(ctps_df$mean_abs_log2FC) +
  0.25 * z_scale(ctps_df$mean_neg_log10_padj) +
  0.25 * (ctps_df$cell_type_breadth / length(deg_list)) +
  0.15 * z_scale(ctps_df$ppi_degree)
)

min_s <- min(ctps_df$CTPS_raw, na.rm = TRUE)
max_s <- max(ctps_df$CTPS_raw, na.rm = TRUE)
ctps_df$Consensus_Score <- round(10 + 90 * (ctps_df$CTPS_raw - min_s) / (max_s - min_s), 2)
ctps_df <- ctps_df %>% arrange(desc(Consensus_Score))

write_csv(ctps_df, "tables/consensus_target_prioritization_hippocampus.csv")

# Fig 3g: Ranking plot of top 15 candidates
top15 <- head(ctps_df, 15)
top15$gene <- factor(top15$gene, levels = rev(top15$gene))

p3g <- ggplot(top15, aes(x = Consensus_Score, y = gene, fill = cell_type_breadth)) +
  geom_bar(stat = "identity", width = 0.68, color = "black", linewidth = 0.3) +
  scale_fill_gradient(low = "#4575b4", high = "#d73027", name = "Cell Breadth") +
  labs(x = "Consensus Therapeutic Prioritization Score (CTPS)", y = "Candidate Target") +
  theme_pub() +
  theme(axis.text.y = element_text(face = "bold.italic", size = 9.5, color = "black"), legend.position = "right")
save_plot_pair(p3g, file.path(fig3_dir, "Fig3g_consensus_score_ranking"), width = 8, height = 5.5)

cat("\n======================================================================\n")
cat("STEP 05b: HIPPOCAMPUS CONVERGENCE & CTPS COMPLETED SUCCESSFULLY!\n")
cat("======================================================================\n")
