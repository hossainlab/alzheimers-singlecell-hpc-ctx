# Figure 2: Lineage-Specific Deep Dives (Hippocampus)

> **Plain-Language Summary for Collaborators**  
> In Step 1, we surveyed the overall cellular composition of the Hippocampus. In Step 2, we zoom into each of the 6 hippocampal cell lineages to uncover the exact molecular pathways disrupted inside each cell type during Alzheimer's Disease.

---

## 📁 Subfolder Organization
This directory contains 6 dedicated folders, one for each investigated hippocampal lineage:
- `Astrocyte/` — Hippocampal support and synaptic metabolic cells
- `Endothelial/` — Capillary blood–brain barrier walls
- `Inhibitory_neuron/` — Hippocampal GABAergic interneurons regulating memory rhythms
- `Microglia/` — Resident immune sentinels responding to tau tangles and amyloid
- `Oligodendrocyte/` — Myelin electrical insulation for memory circuits
- `Pericyte/` — Vascular mural stabilizing cells

---

## 🖼️ Standard 8-Panel Suite in Each Cell-Type Folder
Inside every cell-type folder, you will find the exact same 8-panel suite and a full numerical DEG spreadsheet:

### 1. `Fig2a_<celltype>_subset_umap.*` — Subcluster Specialization Map
- **What it shows**: High-resolution re-clustering of only this cell lineage, split side-by-side by Control vs Alzheimer's.
- **In plain English**: Reveals which functional sub-teams emerge or expand under disease stress in the memory center.

### 2. `Fig2b_<celltype>_subtype_markers_dotplot.*` — Subtype Genetic Signatures
- **What it shows**: The distinctive marker genes defining each sub-team. Larger dots mean more cells have the gene; darker red means higher activity.

### 3. `Fig2c_<celltype>_subtype_counts_barplot.*` — Sub-Team Balance Shifts
- **What it shows**: Direct cell counts for each sub-team in Control (blue) vs Alzheimer's (red).

### 4. `Fig2d_<celltype>_deg_volcano.*` — The Disease Volcano Plot
- **What it shows**: The statistical test for what is altered in Alzheimer's.
- **How to read it in plain English**:
  - **Horizontal axis (x-axis)**: The Fold Change. Points to the right mean the gene is **increased** in Alzheimer's; points to the left mean the gene is **decreased**.
  - **Vertical axis (y-axis)**: Statistical confidence ($-\log_{10} p$-value). The higher up a dot is, the more certain we are it is not a statistical accident.
  - **Colors**: **Red dots** are significantly Upregulated ($p_{\text{adj}} < 0.05, \log_2\text{FC} \ge 0.25$); **blue dots** are significantly Downregulated; **grey dots** are unchanged.

### 5. `Fig2e_<celltype>_top10_degs.*` — Top 10 Broken Tools
- **What it shows**: A bar chart displaying the 5 largest increases (red) and 5 largest decreases (blue) in this hippocampal cell type.

### 6. `Fig2f_<celltype>_kegg_pathways.*` — KEGG Biological Assembly Lines
- **What it shows**: Tests whether entire biochemical pathways (e.g., oxidative phosphorylation, MAPK signaling, synaptic transmission) are malfunctioning.

### 7. `Fig2g_<celltype>_go_bp.*` — Gene Ontology Biological Processes
- **What it shows**: Broader biological duties affected in this cell type (e.g., blood–brain barrier maintenance, immune activation, synaptic plasticity).

### 8. `Fig2h_<celltype>_ppi_network.*` — Protein Physical Teamwork Web
- **What it shows**: An interactome web based on STRING database interactions showing which altered proteins physically bind and cooperate with each other.

### 9. `DEGs_<celltype>.csv` — Full Differential Expression Table
- **What it contains**: Complete, un-truncated numerical results for every detected gene (gene, p_val, avg_log2FC, pct.1, pct.2, p_val_adj).
