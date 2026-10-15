# Figure 2: Lineage-Specific Deep Dives (Prefrontal Cortex)

> **Plain-Language Summary for Collaborators**  
> In Step 1, we took a bird’s-eye census of the whole tissue. In Step 2, we zoom directly into each of the 6 major cell lineages to open the hood and see what specific machinery caught fire or failed inside that cell type during Alzheimer's Disease.

---

## 📁 Subfolder Organization
This directory contains 6 dedicated folders, one for each investigated cell lineage:
- `Astrocyte/` — Brain support and metabolic cells
- `Endothelial/` — Capillary blood–brain barrier walls
- `Inhibitory_neuron/` — Cortical inhibitory interneurons
- `Microglia/` — Resident immune firefighters
- `Oligodendrocyte/` — Myelin electrical insulators
- `Pericyte/` — Vascular mural stabilizing cells

---

## 🖼️ Standard 8-Panel Suite in Each Cell-Type Folder
Inside every cell-type folder, you will find the exact same 8-panel suite and a full numerical DEG spreadsheet:

### 1. `Fig2a_<celltype>_subset_umap.*` — Subcluster Specialization Map
- **What it shows**: High-resolution re-clustering of only this cell lineage, split side-by-side by Control vs Alzheimer's.
- **In plain English**: Cells are not clones of each other; they specialize into distinct functional sub-teams. This map shows which sub-teams emerge or expand under disease stress.

### 2. `Fig2b_<celltype>_subtype_markers_dotplot.*` — Subtype Genetic Signatures
- **What it shows**: The distinctive marker genes defining each sub-team. Larger dots mean more cells have the gene; darker red means higher activity.

### 3. `Fig2c_<celltype>_subtype_counts_barplot.*` — Sub-Team Balance Shifts
- **What it shows**: Direct cell counts for each sub-team in Control (blue) vs Alzheimer's (red). Highlights disease-associated cellular shifts.

### 4. `Fig2d_<celltype>_deg_volcano.*` — The Disease Volcano Plot
- **What it shows**: The core statistical test for what is altered in Alzheimer's.
- **How to read it in plain English**:
  - **Horizontal axis (x-axis)**: The Fold Change. Points to the right mean the gene is **increased** in Alzheimer's; points to the left mean the gene is **decreased**.
  - **Vertical axis (y-axis)**: Statistical confidence ($-\log_{10} p$-value). The higher up a dot is, the more certain we are it is not a statistical accident.
  - **Colors**: **Red dots** are significantly Upregulated ($p_{\text{adj}} < 0.05, \log_2\text{FC} \ge 0.25$); **blue dots** are significantly Downregulated; **grey dots** are unchanged. Key disease genes are directly labeled.

### 5. `Fig2e_<celltype>_top10_degs.*` — Top 10 Broken Tools
- **What it shows**: A straightforward bar chart displaying the 5 largest increases (red) and 5 largest decreases (blue) in this cell type.

### 6. `Fig2f_<celltype>_kegg_pathways.*` — KEGG Biological Assembly Lines
- **What it shows**: Instead of looking at individual genes, this tests whether entire biochemical assembly lines are broken (e.g., oxidative phosphorylation, MAPK signaling, cytokine cascades).

### 7. `Fig2g_<celltype>_go_bp.*` — Gene Ontology Biological Processes
- **What it shows**: The broader biological duties affected in this cell type (e.g., blood vessel maintenance, immune activation, synaptic transmission).

### 8. `Fig2h_<celltype>_ppi_network.*` — Protein Physical Teamwork Web
- **What it shows**: An interactome web based on STRING database protein interactions. It reveals which of the altered proteins physically touch and work together inside the cell.

### 9. `DEGs_<celltype>.csv` — Full Differential Expression Table
- **What it contains**: Complete, un-truncated numerical results for every detected gene:
  - `gene`: Official HGNC gene symbol
  - `p_val`: Unadjusted Wilcoxon p-value
  - `avg_log2FC`: Log2-transformed fold change (positive = elevated in AD, negative = depleted in AD)
  - `pct.1`: Fraction of AD cells expressing the gene
  - `pct.2`: Fraction of Control cells expressing the gene
  - `p_val_adj`: Benjamini–Hochberg false-discovery-rate adjusted p-value
