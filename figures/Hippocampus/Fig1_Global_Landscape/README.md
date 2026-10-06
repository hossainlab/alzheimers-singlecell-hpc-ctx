# Figure 1: Global Tissue Landscape (Hippocampus)

> **Plain-Language Summary for Collaborators**  
> This folder contains the unbiased cellular census of the human Hippocampus — the brain's primary memory hub. It shows the distribution of all cell types living in the hippocampal tissue and quantifies whether their numbers changed in Alzheimer's Disease.

---

## 🖼️ Panel Walkthrough (What You Are Looking At)

### 1. `Fig1a_hippocampus_condition_umap.*` — Clinical Condition Map
- **What it shows**: A 2D map (UMAP) where every dot represents an individual hippocampal cell nucleus. Blue dots are healthy Control donors; red dots are Alzheimer's Disease donors.
- **In plain English**: The thorough mixing of red and blue cells across all clusters verifies that our single-nucleus profiling is free of batch distortion and technical artifacts.

### 2. `Fig1b_hippocampus_celltype_split_umap.*` — Cell Type Neighborhoods
- **What it shows**: The hippocampal cell map split side-by-side (Control on the left, Alzheimer's on the right) and colored by the 6 major cell lineages:
  - Astrocytes (`#3C5488`, Navy)
  - Endothelial cells (`#4DBBD5`, Cyan)
  - Inhibitory neurons (`#8491B4`, Slate)
  - Microglia (`#E64B35`, Red)
  - Oligodendrocytes (`#F39B7F`, Peach)
  - Pericytes (`#00A087`, Teal)

### 3. `Fig1c_hippocampus_celltype_stacked_bar.*` — Total Cell Population Census
- **What it shows**: A 100% stacked bar chart displaying the overall proportion of each cell type in Control vs Alzheimer's donors, with exact percentages printed on each block.

### 4. `Fig1d_hippocampus_donor_proportion_boxplot.*` — Patient-by-Patient Statistical Test
- **What it shows**: Donor-level frequency boxplots across all patients, testing statistical significance with the Wilcoxon rank-sum test.
- **Key finding**: Confirms a concordant loss of microvascular endothelial cells and expansion of inflammatory microglia in the hippocampal memory formation center.

### 5. `Fig1e_hippocampus_canonical_markers_dotplot.*` — Quality-Control Identity Badges
- **What it shows**: Confirms that each hippocampal cell type expresses its canonical identity markers (*AQP4/GFAP* for astrocytes, *CLDN5/PECAM1* for blood vessels, *PTPRC/CX3CR1* for microglia, *MBP/MOG* for oligodendrocytes, *PDGFRB/RGS5* for pericytes, *GAD1/GAD2* for inhibitory neurons).

---

## 📄 File Formats
Each panel is provided in two publication-ready formats:
- `.png`: High-resolution 300 DPI raster image for presentations and reports.
- `.pdf`: Vector graphic without internal title clutter, optimized for manuscript typesetting.
