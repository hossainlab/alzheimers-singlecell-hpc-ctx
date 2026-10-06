# ==============================================================================
# Script 06b: Hippocampus Intercellular Crosstalk (CellChat v2.2 - Figure 4)
# - Total interaction count & interaction strength in Hippocampus (Control vs AD)
# - Clean circular differential interaction network
# - Information flow ranking of conserved vs altered signaling pathways
# - APP signaling circular pathway network across NVU lineages
# - Increased ligand-receptor signaling pairs bubble plot
# Outputs to: figures/Hippocampus/Fig4_Intercellular_Communication/
# ==============================================================================

suppressPackageStartupMessages({
  library(CellChat)
  library(Seurat)
  library(igraph)
  library(dplyr)
  library(readr)
  library(ggplot2)
  library(patchwork)
  library(grid)
})

cat("======================================================================\n")
cat("Starting Step 06b: Hippocampus Intercellular Communication Modeling\n")
cat("======================================================================\n")

fig4_dir <- "figures/Hippocampus/Fig4_Intercellular_Communication"
dir.create(fig4_dir, showWarnings = FALSE, recursive = TRUE)
dir.create("tables", showWarnings = FALSE, recursive = TRUE)
dir.create("data", showWarnings = FALSE, recursive = TRUE)

palette_6 <- c(
  "Astrocyte"         = "#3C5488",
  "Endothelial"       = "#4DBBD5",
  "Inhibitory neuron" = "#8491B4",
  "Microglia"         = "#E64B35",
  "Oligodendrocyte"   = "#F39B7F",
  "Pericyte"          = "#00A087"
)

# Load existing precomputed Hippocampus CellChat or compute
cc_hpc_path <- "data/cellchat_hpc_merged.rds"

if (file.exists(cc_hpc_path)) {
  cat("Loading existing Hippocampus CellChat object from:", cc_hpc_path, "\n")
  cc_hpc_data <- readRDS(cc_hpc_path)
  cc_ctrl_hpc  <- cc_hpc_data$Control
  cc_ad_hpc    <- cc_hpc_data$AD
  cellchat_hpc <- cc_hpc_data$merged
} else {
  cat("Computing CellChat objects from annotated Hippocampus Seurat object...\n")
  hpc <- readRDS("data/seurat_hpc_annotated.rds")
  target_cts <- c("Astrocyte", "Endothelial", "Inhibitory neuron", "Microglia", "Oligodendrocyte", "Pericyte")
  sub_hpc <- subset(hpc, subset = cell_type %in% target_cts)
  
  set.seed(42)
  sampled_cells <- sub_hpc@meta.data %>%
    tibble::rownames_to_column("cell_id") %>%
    group_by(Condition, cell_type) %>%
    slice_sample(n = 500) %>%
    pull(cell_id)
  sub_hpc <- subset(sub_hpc, cells = sampled_cells)
  
  ctrl_seurat <- subset(sub_hpc, subset = Condition == "Control")
  ad_seurat   <- subset(sub_hpc, subset = Condition == "AD")
  
  process_cc <- function(seurat_sub, cond_name) {
    mat <- GetAssayData(seurat_sub, assay = "RNA", layer = "data")
    cc <- createCellChat(object = mat, meta = seurat_sub@meta.data, group.by = "cell_type")
    cc@DB <- CellChatDB.human
    cc <- subsetData(cc)
    cc <- identifyOverExpressedGenes(cc)
    cc <- identifyOverExpressedInteractions(cc)
    cc <- computeCommunProb(cc, type = "triMean", raw.use = TRUE)
    cc <- filterCommunication(cc, min.cells = 10)
    cc <- computeCommunProbPathway(cc)
    cc <- aggregateNet(cc)
    return(cc)
  }
  
  cc_ctrl_hpc <- process_cc(ctrl_seurat, "Control")
  cc_ad_hpc   <- process_cc(ad_seurat, "AD")
  cellchat_hpc <- mergeCellChat(list(Control = cc_ctrl_hpc, AD = cc_ad_hpc), add.names = c("Control", "AD"))
  saveRDS(list(Control = cc_ctrl_hpc, AD = cc_ad_hpc, merged = cellchat_hpc), cc_hpc_path)
}

# 1. Export Interaction Tables
df_net_ctrl_hpc <- subsetCommunication(cc_ctrl_hpc)
df_net_ad_hpc   <- subsetCommunication(cc_ad_hpc)
write_csv(df_net_ctrl_hpc, "tables/CellChat_interactions_Hippocampus_Control.csv")
write_csv(df_net_ad_hpc,   "tables/CellChat_interactions_Hippocampus_AD.csv")
write_csv(df_net_ctrl_hpc, file.path(fig4_dir, "Hippocampus_CellChat_interactions_Control.csv"))
write_csv(df_net_ad_hpc,   file.path(fig4_dir, "Hippocampus_CellChat_interactions_AD.csv"))

# 2. Fig 4a: Number of Interactions & Interaction Strength
p4a_count_hpc  <- compareInteractions(cellchat_hpc, show.legend = FALSE, group = c(1, 2), measure = "count")
p4a_weight_hpc <- compareInteractions(cellchat_hpc, show.legend = FALSE, group = c(1, 2), measure = "weight")

pdf(file.path(fig4_dir, "Fig4a_intercellular_interaction_strength.pdf"), width = 7.0, height = 4.2)
print(p4a_count_hpc + p4a_weight_hpc)
dev.off()

png(file.path(fig4_dir, "Fig4a_intercellular_interaction_strength.png"), width = 2100, height = 1260, res = 300)
print(p4a_count_hpc + p4a_weight_hpc)
dev.off()
cat("  Saved Fig4a\n")

# Clean circos helper
plot_clean_circos <- function(mat, color_vec = palette_6, is_diff = FALSE, main_title = "") {
  cell_names <- rownames(mat)
  node_cols <- color_vec[cell_names]
  node_cols[is.na(node_cols)] <- "#999999"
  
  diag(mat) <- 0
  g <- graph_from_adjacency_matrix(mat, mode = "directed", weighted = TRUE)
  coords <- layout_in_circle(g)
  
  if (is_diff) {
    edge_weights <- E(g)$weight
    edge_cols <- ifelse(edge_weights > 0, "#d73027", "#4575b4")
    edge_widths <- 1 + 5 * (abs(edge_weights) / max(abs(edge_weights), 1e-6))
  } else {
    edge_weights <- E(g)$weight
    edge_cols <- rgb(0.2, 0.2, 0.2, 0.35)
    edge_widths <- 1 + 4 * (edge_weights / max(edge_weights, 1e-6))
  }
  
  plot(
    g,
    layout = coords,
    vertex.color = node_cols,
    vertex.size = 28,
    vertex.frame.color = "white",
    vertex.frame.width = 1.5,
    vertex.label = cell_names,
    vertex.label.family = "sans",
    vertex.label.font = 2,
    vertex.label.cex = 0.85,
    vertex.label.dist = 2.4,
    vertex.label.color = "black",
    edge.color = edge_cols,
    edge.width = edge_widths,
    edge.arrow.size = 0.45,
    edge.curved = 0.25,
    margin = c(0.2, 0.2, 0.2, 0.2),
    main = main_title
  )
}

# 3. Fig 4b: Differential Interaction Network
net_diff_hpc <- cc_ad_hpc@net$count - cc_ctrl_hpc@net$count
pdf(file.path(fig4_dir, "Fig4b_differential_interaction_network.pdf"), width = 7.5, height = 7.5)
plot_clean_circos(net_diff_hpc, is_diff = TRUE)
dev.off()

png(file.path(fig4_dir, "Fig4b_differential_interaction_network.png"), width = 2250, height = 2250, res = 300)
plot_clean_circos(net_diff_hpc, is_diff = TRUE)
dev.off()
cat("  Saved Fig4b\n")

# 4. Fig 4c: Signaling Information Flow Ranking
p4c_hpc <- rankNet(cellchat_hpc, mode = "comparison", stacked = TRUE, do.stat = TRUE) +
  theme(axis.text = element_text(size = 9, color = "black"), axis.title = element_text(face = "bold", size = 10))

pdf(file.path(fig4_dir, "Fig4c_signaling_information_flow_ranking.pdf"), width = 7.5, height = 6.5)
print(p4c_hpc)
dev.off()

png(file.path(fig4_dir, "Fig4c_signaling_information_flow_ranking.png"), width = 2250, height = 1950, res = 300)
print(p4c_hpc)
dev.off()
cat("  Saved Fig4c\n")

# 5. Fig 4d: APP Signaling Pathway Network
pdf(file.path(fig4_dir, "Fig4d_nvu_pathway_circos_network.pdf"), width = 11, height = 5.5)
par(mfrow = c(1, 2), mar = c(2, 2, 3, 2))
if ("APP" %in% names(cc_ctrl_hpc@netP$pathways)) {
  plot_clean_circos(cc_ctrl_hpc@netP$prob[, , "APP"], is_diff = FALSE, main_title = "Hippocampus Control (APP)")
} else plot.new()
if ("APP" %in% names(cc_ad_hpc@netP$pathways)) {
  plot_clean_circos(cc_ad_hpc@netP$prob[, , "APP"], is_diff = FALSE, main_title = "Hippocampus AD (APP)")
} else plot.new()
dev.off()

png(file.path(fig4_dir, "Fig4d_nvu_pathway_circos_network.png"), width = 3300, height = 1650, res = 300)
par(mfrow = c(1, 2), mar = c(2, 2, 3, 2))
if ("APP" %in% names(cc_ctrl_hpc@netP$pathways)) {
  plot_clean_circos(cc_ctrl_hpc@netP$prob[, , "APP"], is_diff = FALSE, main_title = "Hippocampus Control (APP)")
} else plot.new()
if ("APP" %in% names(cc_ad_hpc@netP$pathways)) {
  plot_clean_circos(cc_ad_hpc@netP$prob[, , "APP"], is_diff = FALSE, main_title = "Hippocampus AD (APP)")
} else plot.new()
dev.off()
cat("  Saved Fig4d\n")

# 6. Fig 4e: Ligand-Receptor Communication Bubble Plot
p4e_hpc <- netVisual_bubble(
  cellchat_hpc,
  sources.use = c("Microglia", "Endothelial", "Astrocyte"),
  targets.use = c("Astrocyte", "Endothelial", "Pericyte"),
  comparison = c(1, 2),
  max.dataset = 2,
  angle.x = 45,
  remove.isolate = TRUE
) +
  theme(axis.text = element_text(size = 8.5, color = "black"), axis.title = element_text(face = "bold", size = 9.5))

pdf(file.path(fig4_dir, "Fig4e_ligand_receptor_communication_bubble.pdf"), width = 10, height = 7)
print(p4e_hpc)
dev.off()

png(file.path(fig4_dir, "Fig4e_ligand_receptor_communication_bubble.png"), width = 3000, height = 2100, res = 300)
print(p4e_hpc)
dev.off()
cat("  Saved Fig4e\n")

cat("\n======================================================================\n")
cat("STEP 06b: HIPPOCAMPUS CELLCHAT MODELING COMPLETED SUCCESSFULLY!\n")
cat("======================================================================\n")
