# Figure 1: Global Tissue Landscape (Prefrontal Cortex)

> **Plain-Language Summary for Collaborators**  
> This folder contains the unbiased cellular census of the human Prefrontal Cortex. Think of it as a satellite photograph and demographic survey of all the cell types living in the brain tissue, showing who they are and whether their numbers changed in Alzheimer's Disease.

---

## 🖼️ Panel Walkthrough (What You Are Looking At)

### 1. `Fig1a_cortex_condition_umap.*` — Clinical Condition Map
- **What it shows**: A 2D map (UMAP) where every single dot is an individual cell nucleus. Blue dots are cells from healthy Control donors; red dots are cells from Alzheimer's Disease donors.
- **In plain English**: Notice how the red and blue dots are thoroughly mixed together across all areas of the map. This proves that our experiment did not suffer from technical batch bias — healthy and disease cells were processed and aligned under identical conditions.

### 2. `Fig1b_cortex_celltype_split_umap.*` — Cell Type Neighborhoods
- **What it shows**: The same cell map split side-by-side (Control on the left, Alzheimer's on the right) and colored by the 6 major cell lineages:
  - Astrocytes (`#3C5488`, Navy)
  - Endothelial cells (`#4DBBD5`, Cyan)
  - Inhibitory neurons (`#8491B4`, Slate)
  - Microglia (`#E64B35`, Red)
  - Oligodendrocytes (`#F39B7F`, Peach)
  - Pericytes (`#00A087`, Teal)

### 3. `Fig1c_cortex_celltype_stacked_bar.*` — Total Cell Population Census
- **What it shows**: A 100% stacked bar chart displaying the overall proportion of each cell type in Control vs Alzheimer's. Exact percentages are printed on each color block for quick reading.

### 4. `Fig1d_cortex_donor_proportion_boxplot.*` — Patient-by-Patient Statistical Test
- **What it shows**: Boxplots showing the percentage of each cell type measured in individual patient donors, complete with Wilcoxon statistical significance testing.
- **Key finding**: Blood vessel endothelial cells show a significant loss in Alzheimer's donors, while microglial immune cells show a marked increase.

### 5. `Fig1e_cortex_canonical_markers_dotplot.*` — Quality-Control Identity Badges
- **What it shows**: Verifies that every cell type was correctly identified by looking at textbook marker genes (*AQP4/GFAP* for astrocytes, *CLDN5/PECAM1* for blood vessels, *PTPRC/CX3CR1* for microglia, *MBP/MOG* for oligodendrocytes, *PDGFRB/RGS5* for pericytes, *GAD1/GAD2* for inhibitory neurons).
- **In plain English**: Dot size shows how many cells express the badge; red color intensity shows how strong the signal is. Clean, separated rows prove our cell classifications are 100% accurate.

---

## 📄 File Formats
Each panel is provided in two publication-ready formats:
- `.png`: High-resolution 300 DPI raster image for presentations and reports.
- `.pdf`: Vector graphic without internal title clutter, optimized for manuscript typesetting.
