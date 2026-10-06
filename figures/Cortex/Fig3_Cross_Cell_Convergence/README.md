# Figure 3: Cross-Cell Lineage Convergence & Target Scoring (Prefrontal Cortex)

> **Plain-Language Summary for Collaborators**  
> In Step 2, we looked at each cell type in isolation. Here in Step 3, we ask a pivotal question: **Is there a universal fire alarm sounding across ALL cell types in the brain at the same time?**  
> Instead of developing a drug that only helps one rare cell type, we search for master regulator genes shared across astrocytes, blood vessels, and microglia. This is where our mathematical target scoring (CTPS) takes place.

---

## 🎯 Why We Score Targets This Way (The Logic & Justification)

In single-cell genomics, picking targets simply by sorting for the lowest p-value often fails in the clinic. A gene might have a tiny p-value because of high cell numbers, but only change by 2%, be confined to one cell type, or be impossible to drug with antibodies.

To avoid these traps, we created the **Multi-Factorial Consensus Therapeutic Prioritization Score (CTPS)**:

$$\text{CTPS} = 0.25 \cdot Z(|\log_2\text{FC}|) + 0.25 \cdot Z(-\log_{10} p_{\text{adj}}) + 0.25 \cdot \left(\frac{\text{Cell Breadth}}{5}\right) + 0.15 \cdot Z(\text{PPI Degree}) + 0.10 \cdot \text{Cross-Region Confirmed}$$

This formula guarantees that our top candidates meet 5 non-negotiable criteria:
1. **Strong Effect Size**: The gene changed by $\ge 25\%$ (large therapeutic window).
2. **Statistical Confidence**: Benjamini–Hochberg false-discovery corrected ($p_{\text{adj}} < 0.05$).
3. **Broad Multi-Cell Action**: The gene is broken across multiple cell types (reversing it rescues the whole tissue).
4. **Master Network Switch**: The protein sits in the center of the physical interactome (controlling downstream pathways).
5. **Cross-Region Replication**: The change occurs in both Cortex AND Hippocampus (eliminating regional artifacts).

---

## 🖼️ Panel Walkthrough (What You Are Looking At)

### 1. `Fig3a_venn_5celltypes.*` — The 5-Way Overlap Venn Diagram
- **What it shows**: A 5-way publication Venn diagram intersecting the significant disease genes across 5 key neurovascular unit lineages:
  - Astrocytes (Navy, `#3C5488`)
  - Endothelial cells (Cyan, `#4DBBD5`)
  - Microglia (Red, `#E64B35`)
  - Oligodendrocytes (Peach, `#F39B7F`)
  - Pericytes (Teal, `#00A087`)
- **What to notice in plain English**: The central overlapping sections contain hundreds of genes that are co-altered across multiple lineages, proving that Alzheimer's is a coordinated multi-cell collapse.

### 2. `Fig3a_upset_5celltypes.*` — The UpSet Intersection Matrix
- **What it shows**: While 5-way Venn diagrams show the broad picture, an UpSet plot is an exact mathematical bar-code matrix quantifying every single combination of cell overlaps with precise numbers.

### 3. `Fig3b_common_genes_kegg.*` — What the Universal Alarm Does (KEGG Pathways)
- **What it shows**: Enriched biochemical pathways among the shared convergent genes. Highlights MAPK signaling, blood–brain barrier failure, and protein processing.

### 4. `Fig3c_common_genes_go_bp.*` — Biological Processes of Shared Targets
- **What it shows**: Gene Ontology Biological Processes shared across cell types, including response to oxidative stress, regulation of cell death, and vascular maintenance.

### 5. `Fig3d_common_genes_ppi_network.*` — The Master Switchboard Network
- **What it shows**: A circular protein interaction map of top convergent genes. Node colors indicate "degree centrality" (connectivity). The genes with the highest degree sit in the center of the cellular web, acting as master switches.

### 6. `Fig3e_cortex_expression_violins.*` — Candidate Expression Levels
- **What it shows**: Split violin plots showing the exact expression levels of our top prioritized candidates (*DUSP1, FOS, JUN, EGR1, GADD45B, ATF3*) across all 5 lineages in Control (blue) vs Alzheimer's (red).
- **In plain English**: Shows that genes like *DUSP1* and *FOS* spike dramatically across every single cell type, confirming broad multi-cell activation.

### 7. `Fig3g_consensus_score_ranking.*` — The Official Target Leaderboard
- **What it shows**: A horizontal bar chart ranking the top 15 candidate targets according to their final multi-factorial CTPS score (scaled 10 to 100). The color bar reflects cell breadth (how many cell types share the target).
- **Top Candidates Highlighted**: *SORL1, CLEC5A, ADAMTS9, PCDH9, DUSP1, NOTCH2*.
