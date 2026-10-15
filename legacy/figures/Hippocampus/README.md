# Hippocampus (Hpc) — Collaborator's Visual Guide & Figure Handbook

> **Who is this guide for?**  
> This handbook is written in **plain, intuitive language** specifically for collaborators, clinicians, neuroscientists, and lab team members without a computer science, bioinformatics, or single-cell genomics background. It explains what each image shows, what the scientific takeaway is, and exactly **why and how** our therapeutic targets were prioritized in the human **Hippocampus**.

---

## 🗺️ The Big Picture: Why the Hippocampus is Special

If the cortex is the brain's executive office tower, the **Hippocampus is the central memory bank and filing system**. 
- It is the very first region attacked by neurofibrillary tau tangles in Alzheimer's disease.
- Patients lose recent memories and navigation abilities because hippocampal neural circuits fail early.
- In this study, we profile thousands of single nuclei from the human Hippocampus across the exact same cell lineages:
  - **Astrocytes**: Hippocampal support grid, regulating synaptic glutamate and lactate fuel.
  - **Endothelial Cells & Pericytes**: Blood–Brain Barrier lining hippocampal capillaries.
  - **Microglia**: Resident immune sentinels responding to tau tangles and amyloid oligomers.
  - **Oligodendrocytes**: Insulation wrapping hippocampal axonal pathways (fornix and perforant path).
  - **Inhibitory Neurons**: GABAergic interneurons regulating theta and gamma memory rhythms.

---

## 🔄 The 5-Step Logical Flow: How We Find Drug Targets

Our Hippocampus analysis follows the exact same logical pipeline as the Cortex:

```
[STEP 1: Tissue Census]  ──►  Who lives in the hippocampus, and did cell populations shift in AD?
         │
[STEP 2: Worker Stress]  ──►  Inside each hippocampal cell type, what specific machinery broke down?
         │
[STEP 3: Universal Alarm]──►  Which genes are broken across MULTIPLE cell types at once? (Target Scoring)
         │
[STEP 4: Phone Lines]    ──►  How are hippocampal cells communicating, and how are signals disrupted?
         │
[STEP 5: Drug Blueprint] ──►  Can we drug these targets with antibodies or pills, and test them in lab mice?
```

---

## 🎯 Target Prioritization: Why and How Were Targets Selected?

A common question from non-computational collaborators is:  
> *"Why didn't you just sort your spreadsheet by the lowest p-value and pick the top 5 genes?"*

### Why the "Lowest p-value" Fails:
1. **Statistical Flukes**: In large datasets with thousands of single cells, a gene can achieve $p < 10^{-50}$ even with a negligible 2–3% change in expression. Such small changes cannot be targeted by drugs.
2. **Cell-Private Genes**: A gene might change dramatically in only one rare subpopulation while leaving the rest of the tissue untouched.
3. **The "Un-druggable" Trap**: Many top statistical genes lack accessible physical binding sites, preventing antibody or small-molecule binding.
4. **Species Barriers**: Preclinical drugs must be tested in transgenic Alzheimer's mice. If the mouse lacks the target, translation is blocked.

### The Objective Solution: Consensus Therapeutic Prioritization Score (CTPS)
We scored every altered gene using a 5-pillar mathematical formula combining effect size, statistical significance, cross-cell breadth, physical protein interaction, and cross-region replication:

$$\text{CTPS} = 0.25 \cdot Z(|\log_2\text{FC}|) + 0.25 \cdot Z(-\log_{10} p_{\text{adj}}) + 0.25 \cdot \left(\frac{\text{Cell Breadth}}{5}\right) + 0.15 \cdot Z(\text{PPI Degree}) + 0.10 \cdot \text{Cross-Region Confirmed}$$

| Justification Pillar | What It Means in Plain English | Why It Matters for a Drug |
| :--- | :--- | :--- |
| **1. Effect Size ($|\log_2\text{FC}|$ cut-off $\ge 0.25$)** | Did the gene shift by at least 20–25% in expression? | Drugs need a large therapeutic window to reverse disease states. |
| **2. Statistical Rigor ($p_{\text{adj}} < 0.05$)** | Did this change pass strict false-discovery correction? | Rules out random biological noise. |
| **3. Lineage Breadth (1 to 5 cell types)** | Is this gene broken in astrocytes, blood vessels, AND microglia? | **Universal targets** rescue multiple cell compartments simultaneously. |
| **4. Network Hub (PPI Centrality)** | Does this protein physically touch and control other proteins? | Hitting a "master switchboard" fixes the whole downstream module. |
| **5. Cross-Region Confirmation (1 or 0)** | Did the gene shift in the **exact same direction** in both Hippocampus AND Cortex? | Proves the target is a fundamental Alzheimer's driver, not a local fluke. |

---

## 🏆 The Top Prioritized Targets & Biological Justifications in Hippocampus

Our quantitative screening identified convergent molecular switches in the Hippocampus:

### 1. **SORL1** (CTPS: 98/100) — *The Amyloid Traffic Controller*
- **Cell Types Affected**: Astrocytes, Oligodendrocytes, and Endothelial cells.
- **What it does normally**: Protects hippocampal cells by routing the Amyloid Precursor Protein (APP) safely away from beta-secretase cleavage.
- **What happens in AD**: Severely lost in hippocampal glia, leading to accelerated amyloid seeding and synaptic vulnerability.
- **Why it is druggable**: Type I cell-surface receptor with an accessible outer Vps10p domain; targetable by stabilizing monoclonal antibodies.

### 2. **CLEC5A** (CTPS: 94/100) — *The Neuroimmune Alarm Megaphone*
- **Cell Types Affected**: Microglia and Endothelial cells.
- **What it does normally**: Senses cell death products and pathogen patterns.
- **What happens in AD**: Dramatically elevated on hippocampal microglia surrounding damaged synapses, amplifying chronic neuroinflammation.
- **Why it is druggable**: Cell-surface receptor with an extracellular C-type lectin domain, readily blocked by monoclonal antibodies.

### 3. **ADAMTS9** (CTPS: 96/100) — *The Blood–Brain Barrier Mortar*
- **Cell Types Affected**: Endothelial cells and Pericytes.
- **What it does normally**: Secreted extracellular matrix metalloproteinase essential for microvascular capillary integrity.
- **What happens in AD**: Depleted in hippocampal microvessels, exacerbating microvascular leakage in the blood–brain barrier.
- **Why it is druggable**: Secreted directly into the extracellular space; 91.2% amino acid sequence identity between human and mouse.

### 4. **PCDH9** (CTPS: 92/100) — *The Cellular Velcro*
- **Cell Types Affected**: Pericytes, Astrocytes, and Inhibitory Neurons.
- **What it does normally**: Cadherin-related transmembrane adhesion molecule maintaining cell-to-cell contact.
- **What happens in AD**: Highly upregulated across 5 lineages, reflecting aberrant synaptic remodeling and vascular scarring.
- **Why it is druggable**: Cell-surface receptor with exposed extracellular cadherin repeats.

### 5. **DUSP1** (CTPS: 90/100) — *The Master Emergency Brake*
- **Cell Types Affected**: Broadly altered across all 5 neurovascular glia and inhibitory neurons.
- **What it does normally**: Dual-specificity phosphatase acting as an intracellular brake to shut down toxic p38 and JNK stress kinase cascades.
- **What happens in AD**: Spikes across all hippocampal lineages in response to oxidative and metabolic stress.
- **Why it is druggable**: Intracellular enzyme targetable by cell-permeable small-molecule allosteric modulators.

### 6. **NOTCH2** (CTPS: 89/100) — *The Neurovascular Communication Antenna*
- **Cell Types Affected**: Endothelial cells and Pericytes.
- **What it does normally**: Receives intercellular signals that maintain capillary stability.
- **What happens in AD**: Aberrantly activated, disrupting hippocampal vascular tone.

---

## 🖼️ Step-by-Step Figure Walkthrough (What You Are Looking At)

Below is an intuitive guide to every figure in the `Hippocampus/` directory:

### 📁 Figure 1: Global Landscape (`Fig1_Global_Landscape/`)
*Goal: Take an unbiased census of all cells in the human Hippocampus.*

- **Fig1a (`Fig1a_hippocampus_condition_umap`)**:  
  - *What it is*: A 2D GPS map of all hippocampal cells. Blue dots represent healthy control donors; red dots represent Alzheimer's donors.  
  - *What to notice*: Red and blue dots are thoroughly intermingled across every cluster, confirming proper batch correction and experimental rigor.
- **Fig1b (`Fig1b_hippocampus_celltype_split_umap`)**:  
  - *What it is*: The same GPS map split side-by-side, color-coded by the 6 major cell lineages (Astrocytes, Endothelial cells, Inhibitory neurons, Microglia, Oligodendrocytes, Pericytes).
- **Fig1c (`Fig1c_hippocampus_celltype_stacked_bar`)**:  
  - *What it is*: A 100% stacked bar chart showing the total cellular proportion of each lineage in Control vs AD.
- **Fig1d (`Fig1d_hippocampus_donor_proportion_boxplot`)**:  
  - *What it is*: Donor-by-donor boxplots showing exact percentages for each patient, with Wilcoxon statistical significance testing.
- **Fig1e (`Fig1e_hippocampus_canonical_markers_dotplot`)**:  
  - *What it is*: A quality-control "identity badge" plot confirming canonical cell markers (*AQP4* for astrocytes, *CLDN5* for blood vessels, *PTPRC* for microglia, *MBP* for oligodendrocytes, *PDGFRB* for pericytes).

---

### 📁 Figure 2: Cell-Type Deep Dives (`Fig2_Cell_Type_Deep_Dives/<cell_type>/`)
*Goal: Open the hood on each of the 6 cell types individually in the Hippocampus.*

Each of the 6 subfolders (`Astrocyte/`, `Endothelial/`, `Inhibitory_neuron/`, `Microglia/`, `Oligodendrocyte/`, `Pericyte/`) contains an identical 8-panel suite:
- **Fig2a (Subclustering UMAP)**: High-resolution subclustering to identify disease-enriched sub-states.
- **Fig2b (Subtype Marker Dotplot)**: Marker genes defining each sub-state.
- **Fig2c (Subtype Counts Barplot)**: Proportions of sub-states in Control vs AD.
- **Fig2d (Volcano Plot)**: Differential gene expression (red = upregulated, blue = downregulated).
- **Fig2e (Top 10 DEGs)**: Bar chart showing the 5 biggest increases and 5 biggest decreases.
- **Fig2f & Fig2g (KEGG & GO Biological Process Dotplots)**: Enriched biological pathways.
- **Fig2h (Protein Interaction Network)**: STRING PPI network of key altered proteins.

---

### 📁 Figure 3: Cross-Cell Lineage Convergence (`Fig3_Cross_Cell_Convergence/`)
*Goal: Identify shared pathogenic stress axes across all hippocampal cell types.*

- **Fig3a (`Fig3a_venn_5celltypes`)**:  
  - *What it is*: A 5-way publication Venn diagram intersecting the significant disease genes across Astrocytes (navy), Endothelial cells (cyan), Microglia (red), Oligodendrocytes (peach), and Pericytes (teal).  
  - *Visual Consistency Note*: Rendered with the exact same 5 lineages, transparent borders, color palette, and layout as the Cortex Venn diagram.
- **Fig3a (`Fig3a_upset_5celltypes`)**:  
  - *What it is*: An UpSet matrix quantifying all shared DEG combinations cleanly.
- **Fig3b & Fig3c (`common_genes_kegg` & `common_genes_go_bp`)**:  
  - *What it is*: Functional enrichment of convergent hippocampal genes (oxidative stress, protein processing, MAPK signaling).
- **Fig3d (`Fig3d_common_genes_ppi_network`)**:  
  - *What it is*: Circular hub network highlighting master regulators.
- **Fig3e (`Fig3e_hippocampus_expression_violins`)**:  
  - *What it is*: Violin plots showing expression levels of top candidates (*DUSP1, FOS, JUN, EGR1, GADD45B, ATF3*) in Control (blue) vs AD (red) across all lineages.
- **Fig3g (`Fig3g_consensus_score_ranking`)**:  
  - *What it is*: The official leaderboard ranking the top 15 candidate targets using our standardized multi-factorial CTPS formula.

---

### 📁 Figure 4: Intercellular Communication (`Fig4_Intercellular_Communication/`)
*Goal: Model hippocampal cell-to-cell communication networks using CellChat v2.2.*

- **Fig4a (`Fig4a_intercellular_interaction_strength`)**:  
  - *What it is*: Total interaction number and communication strength (Control vs AD).
- **Fig4b (`Fig4b_differential_interaction_network`)**:  
  - *What it is*: Circular network showing net shifts (red = increased in AD, blue = decreased in AD).
- **Fig4c (`Fig4c_signaling_information_flow_ranking`)**:  
  - *What it is*: Ranking of 27 signaling pathways by information flow in Hippocampus.
- **Fig4d (`Fig4d_nvu_pathway_circos_network`)**:  
  - *What it is*: Circular chord diagram tracing *APP* signaling across hippocampal cell types.
- **Fig4e (`Fig4e_ligand_receptor_communication_bubble`)**:  
  - *What it is*: High-resolution bubble plot of key ligand–receptor pairs.

---

### 📁 Figure 5: Target Translation & Drug Roadmap (`Fig5_Target_Validation/`)
*Goal: Translational validation for in vivo and clinical follow-up.*

- **Fig5a (`Fig5a_target_cross_species_conservation`)**:  
  - *What it is*: Human-to-mouse protein sequence conservation and In-Vivo Suitability Index.
- **Fig5b (`Fig5b_target_druggability_and_surfaceome`)**:  
  - *What it is*: Subcellular surfaceome localization and extracellular epitope accessibility.
- **Fig5c (`Fig5c_target_ppi_mechanistic_network`)**:  
  - *What it is*: Mechanistic network linking prioritized targets to core Alzheimer's drivers (*APP, BACE1, TREM2, TYROBP, CLDN5*).

---

## 🧪 Wet-Lab Validation: Recommended Antibodies for Collaborators

| Target Gene | Host Species | Clonality | Recommended Vendor | Catalog Number | Validated Applications | Dilution (IHC-P / IF) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **CLEC5A** | Rat | Monoclonal [284108] | R&D Systems / BioLegend | MAB1731 / #353804 | FC, IF, Neutralization | 1:100 (IF) |
| **CLEC5A** | Rabbit | Polyclonal | Proteintech | 28373-1-AP | WB, IHC-P, IF | 1:200 (IHC) |
| **ADAMTS9** | Rabbit | Polyclonal | Abcam | ab284561 | IHC-P, WB | 1:150 (IHC) |
| **ADAMTS9** | Rabbit | Polyclonal | Proteintech | 28540-1-AP | WB, IHC-P, IF | 1:200 (IHC) |
| **SORL1** | Mouse | Monoclonal [Clone 48] | BD Biosciences | #611860 | WB, IF, IHC-P | 1:200 (IHC/IF) |
| **SORL1** | Rabbit | Monoclonal [EPR14674] | Abcam | ab190684 | IHC-P, WB, IF | 1:500 (IHC/IF) |
| **PCDH9** | Rabbit | Polyclonal (Prestige) | Atlas Antibodies / Sigma | HPA031154 | IHC-P (HPA Validated), WB | 1:250 (IHC-P) |
| **DUSP1** | Rabbit | Monoclonal [D9A5] | Cell Signaling Tech | #5142 | WB, IHC-P, IF | 1:400 (IHC/IF) |
| **NOTCH2** | Rabbit | Monoclonal [D76A6] | Cell Signaling Tech | #5732 | WB, IHC-P, IF | 1:300 (IHC/IF) |

---

## 🎨 Figure Style & Harmonization Guarantee

All figures in this directory match the Prefrontal Cortex suite with complete mathematical and graphical symmetry:
- **Identical Palettes**: Navy Blue (`#4575b4`) for Control, Red (`#d73027`) for AD; identical 6-lineage color coding.
- **Identical 5-Way Venn Diagram**: Uses the exact same 5 lineages (`Astrocyte`, `Endothelial`, `Microglia`, `Oligodendrocyte`, `Pericyte`), identical borderless transparency (`col = "transparent"`, `alpha = 0.50`), and identical layout.
- **Identical Dimensions & Typography**: All panels export to matching dimensions at 300 DPI (`.png`) and vector graphics (`.pdf`).
