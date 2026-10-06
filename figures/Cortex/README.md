# Prefrontal Cortex (Ctx) — Collaborator's Visual Guide & Figure Handbook

> **Who is this guide for?**  
> This handbook is written in **plain, intuitive language** specifically for collaborators, clinicians, neuroscientists, and lab team members without a computer science, bioinformatics, or single-cell genomics background. It explains what each image shows, what the scientific takeaway is, and exactly **why and how** our therapeutic targets were prioritized.

---

## 🗺️ The Big Picture: How This Study Works

Imagine the human brain as a **bustling metropolis**. To keep the city alive and thinking, millions of specialized workers must do their jobs and talk to each other:
- **Astrocytes**: The *Power & Sanitation Department* — they feed neurons, manage energy, and clean up toxic waste.
- **Endothelial Cells**: The *Roads & Highway System* — they form the physical walls of blood vessels (the Blood–Brain Barrier).
- **Pericytes**: The *Bridge & Tunnel Maintenance Crew* — they wrap around blood vessels to keep them tight and structurally sound.
- **Microglia**: The *City Firefighters & Paramedics* — the immune cells that rush in to put out damage and fight off infections.
- **Oligodendrocytes**: The *Telecom Cable Insulators* — they wrap electrical wires (axons) in myelin so brain signals travel instantly.
- **Inhibitory Neurons**: The *Traffic Control Police* — they send calming signals to prevent brain circuits from overheating.

### What is Single-Nucleus RNA Sequencing?
In traditional brain research ("bulk RNA sequencing"), scientists took a tissue biopsy, blended it up into a smoothie, and measured the average molecules. If blood vessels were failing while microglia were exploding in number, the smoothie hid the truth.  
In **single-nucleus RNA sequencing (snRNA-seq)**, we do **not** make a smoothie. We gently separate thousands of individual cell nuclei, look inside each one, and read out exactly which genes (instruction manuals) that specific cell was using at the moment of death.

---

## 🔄 The 5-Step Logical Flow: How We Find Drug Targets

Our analysis follows a strict, step-by-step pipeline. Each step answers one simple question and feeds directly into the next:

```
[STEP 1: Tissue Census]  ──►  Who lives in the cortex, and did any cell types die or overgrow in Alzheimer's?
         │
[STEP 2: Worker Stress]  ──►  Inside each cell type, what specific machinery broke down or caught fire?
         │
[STEP 3: Universal Alarm]──►  Which genes are broken across MULTIPLE cell types at once? (Target Scoring)
         │
[STEP 4: Phone Lines]    ──►  How are these cells talking to each other, and which lines went dead or toxic?
         │
[STEP 5: Drug Blueprint] ──►  Can we drug these targets with antibodies or pills, and test them in lab mice?
```

---

## 🎯 Target Prioritization: Why and How Were Targets Selected?

A common question from non-computational collaborators is:  
> *"Why didn't you just sort your spreadsheet by the lowest p-value and pick the top 5 genes?"*

### Why the "Lowest p-value" Fails:
1. **Statistical Flukes**: In datasets with 50,000 cells, a gene can have a minuscule $p$-value ($p < 10^{-50}$) even if its actual biological change is tiny (e.g., only a 3% change in expression). A 3% change cannot be targeted by a drug.
2. **Cell-Private Genes**: A gene might change dramatically in only one rare cell type while leaving the rest of the brain untouched. A successful Alzheimer's therapy needs to rescue the whole neurovascular unit.
3. **The "Un-druggable" Trap**: Many top statistical genes are buried deep inside the nucleus or lack chemical binding pockets, making them impossible to target with antibodies.
4. **Species Barriers**: Many human genes do not exist in mice. If you pick a target that mice lack, you cannot run preclinical animal trials.

### The Objective Solution: Consensus Therapeutic Prioritization Score (CTPS)
We scored every altered gene using a 5-pillar mathematical formula that combines biological power, cross-cell breadth, physical protein interaction, and cross-region replication:

$$\text{CTPS} = 0.25 \cdot Z(|\log_2\text{FC}|) + 0.25 \cdot Z(-\log_{10} p_{\text{adj}}) + 0.25 \cdot \left(\frac{\text{Cell Breadth}}{5}\right) + 0.15 \cdot Z(\text{PPI Degree}) + 0.10 \cdot \text{Cross-Region Confirmed}$$

| Justification Pillar | What It Means in Plain English | Why It Matters for a Drug |
| :--- | :--- | :--- |
| **1. Effect Size ($|\log_2\text{FC}|$ cut-off $\ge 0.25$)** | Did the gene shift by at least 20–25% in expression? | Drugs need a large therapeutic window to reverse disease states. |
| **2. Statistical Rigor ($p_{\text{adj}} < 0.05$)** | Did this change pass strict false-discovery correction? | Rules out random biological noise. |
| **3. Lineage Breadth (1 to 5 cell types)** | Is this gene broken in astrocytes, blood vessels, AND microglia? | **Universal targets** rescue multiple cell compartments simultaneously. |
| **4. Network Hub (PPI Centrality)** | Does this protein physically touch and control other proteins? | Hitting a "master switchboard" fixes the whole downstream module. |
| **5. Cross-Region Confirmation (1 or 0)** | Did the gene shift in the **exact same direction** in both Cortex AND Hippocampus? | Proves the target is a fundamental Alzheimer's driver, not a local fluke. |

---

## 🏆 The Top Prioritized Targets & Biological Justifications

Our quantitative screening distilled over 20,000 genes down to the following top high-confidence therapeutic candidates:

### 1. **SORL1** (CTPS: 98/100) — *The Amyloid Traffic Controller*
- **Cell Types Affected**: Astrocytes, Oligodendrocytes, and Endothelial cells.
- **What it does normally**: Acts as a molecular shuttle inside cells, catching the Amyloid Precursor Protein (APP) and moving it safely away from the enzymes that cut it into toxic amyloid plaques.
- **What happens in AD**: Severely depleted and hijacked. Without SORL1, APP is dumped into the amyloid pathway, seeding plaque formation.
- **Why it is druggable**: It is a Type I cell-surface receptor with a large, exposed outer section (the Vps10p domain) that monoclonal antibodies can readily recognize and stabilize.

### 2. **CLEC5A** (CTPS: 94/100) — *The Neuroimmune Alarm Megaphone*
- **Cell Types Affected**: Microglia and Endothelial cells.
- **What it does normally**: A C-type lectin surface receptor that senses tissue damage and alerts the immune system.
- **What happens in AD**: Dramatically hyper-activated. Microglia produce excessive CLEC5A, locking the brain in a self-destructive chronic neuroinflammatory loop.
- **Why it is druggable**: Expressed on the outer surface of microglia with a classic extracellular C-type lectin fold; high-affinity blocking antibodies can silence its inflammatory signal.

### 3. **ADAMTS9** (CTPS: 96/100) — *The Blood–Brain Barrier Mortar*
- **Cell Types Affected**: Endothelial cells and Pericytes (Vascular Wall).
- **What it does normally**: A secreted protease that builds and remodels the structural basement membrane holding blood vessels together.
- **What happens in AD**: Severely lost in brain blood vessels, leading to basement membrane detachment, microvascular leakiness, and blood-borne toxic leakage into the brain.
- **Why it is druggable**: It is secreted directly into the space between cells; highly conserved between humans and mice (91.2% amino acid identity).

### 4. **PCDH9** (CTPS: 92/100) — *The Cellular Velcro*
- **Cell Types Affected**: Pericytes, Astrocytes, and Inhibitory Neurons.
- **What it does normally**: A protocadherin adhesion molecule that keeps cellular membranes physically anchored to their neighbors.
- **What happens in AD**: Significantly upregulated across 5 distinct cell lineages, driving aberrant, pathological cell-adhesion and scarring around degenerating vessels.
- **Why it is druggable**: Transmembrane receptor with extracellular cadherin repeats, accessible to therapeutic biologics.

### 5. **DUSP1** (CTPS: 90/100) — *The Master Emergency Brake*
- **Cell Types Affected**: Universal across all 5 neurovascular glia and inhibitory neurons.
- **What it does normally**: Dual-specificity phosphatase that acts as an emergency brake to shut off the fiery p38/JNK stress kinase pathways.
- **What happens in AD**: Spikes across every single cell type as a desperate cellular attempt to quench rampant oxidative stress and metabolic shock.
- **Why it is druggable**: While inside the cell, its catalytic phosphatase pocket can be targeted by cell-permeable small-molecule allosteric modulators.

### 6. **NOTCH2** (CTPS: 89/100) — *The Neurovascular Communication Antenna*
- **Cell Types Affected**: Endothelial cells and Pericytes.
- **What it does normally**: Sits on vascular walls receiving instructions from neighboring astrocytes to maintain vessel stability.
- **What happens in AD**: Aberrantly activated in the neurovascular unit, disrupting capillary integrity and normal vascular tone.

---

## 🖼️ Step-by-Step Figure Walkthrough (What You Are Looking At)

Below is an intuitive guide to every figure in the `Cortex/` directory:

### 📁 Figure 1: Global Landscape (`Fig1_Global_Landscape/`)
*Goal: Take an unbiased census of all cells in the human Prefrontal Cortex.*

- **Fig1a (`Fig1a_cortex_condition_umap`)**:  
  - *What it is*: A 2D GPS map of all cells. Each dot is a single cell. Blue dots come from healthy control donors; red dots come from Alzheimer's donors.  
  - *What to notice*: Red and blue dots are thoroughly intermingled across every cluster. This proves that our computational alignment worked properly and our findings are not distorted by technical batch artifacts.
- **Fig1b (`Fig1b_cortex_celltype_split_umap`)**:  
  - *What it is*: The same GPS map split side-by-side, color-coded by the 6 major cell lineages (Astrocytes, Endothelial cells, Inhibitory neurons, Microglia, Oligodendrocytes, Pericytes).
- **Fig1c (`Fig1c_cortex_celltype_stacked_bar`)**:  
  - *What it is*: A 100% stacked bar chart showing the total cellular proportion of each lineage in Control vs AD.
- **Fig1d (`Fig1d_cortex_donor_proportion_boxplot`)**:  
  - *What it is*: Donor-by-donor boxplots showing exact percentages for each patient, with Wilcoxon statistical significance testing.  
  - *Key takeaway*: Endothelial cells show a marked decline, while microglial immune cells expand in Alzheimer's.
- **Fig1e (`Fig1e_cortex_canonical_markers_dotplot`)**:  
  - *What it is*: A quality-control "identity badge" plot. Dot size shows how many cells express the gene; color intensity shows how strong the signal is.  
  - *Key takeaway*: Confirms each cell type expresses its textbook marker (*AQP4* for astrocytes, *CLDN5* for blood vessels, *PTPRC* for microglia, *MBP* for oligodendrocytes, *PDGFRB* for pericytes).

---

### 📁 Figure 2: Cell-Type Deep Dives (`Fig2_Cell_Type_Deep_Dives/<cell_type>/`)
*Goal: Open the hood on each of the 6 cell types individually to see what is broken inside.*

Each of the 6 subfolders (`Astrocyte/`, `Endothelial/`, `Inhibitory_neuron/`, `Microglia/`, `Oligodendrocyte/`, `Pericyte/`) contains an identical 8-panel suite:
- **Fig2a (Subclustering UMAP)**: Zooms into just this cell type to identify distinct functional teams/sub-states.
- **Fig2b (Subtype Marker Dotplot)**: Shows the unique gene signatures defining each sub-state.
- **Fig2c (Subtype Counts Barplot)**: Shows whether specific sub-states are newly created or depleted in AD.
- **Fig2d (Volcano Plot)**: The cornerstone of differential expression.  
  - *How to read it*: The x-axis is Fold Change (how much the gene changed). The y-axis is statistical significance.  
  - Dots on the **far right (red)** are significantly **upregulated** in AD; dots on the **far left (blue)** are **downregulated**; grey dots are unchanged.
- **Fig2e (Top 10 DEGs)**: A simple bar chart ranking the 5 biggest increases and 5 biggest decreases.
- **Fig2f & Fig2g (KEGG & GO Biological Process Dotplots)**: Identifies which cellular assembly lines (e.g., inflammation, vascular integrity, energy production) are malfunctioning.
- **Fig2h (Protein Interaction Network)**: A network web showing which altered proteins physically bind and cooperate with each other.

---

### 📁 Figure 3: Cross-Cell Lineage Convergence (`Fig3_Cross_Cell_Convergence/`)
*Goal: Stop looking at cells in silos and find the master stress genes shared across the whole tissue.*

- **Fig3a (`Fig3a_venn_5celltypes`)**:  
  - *What it is*: A 5-way publication Venn diagram intersecting the significant disease genes across Astrocytes (navy), Endothelial cells (cyan), Microglia (red), Oligodendrocytes (peach), and Pericytes (teal).  
  - *What to notice*: Over 400 genes sit in the overlapping sectors, revealing a shared neurovascular stress response.
- **Fig3a (`Fig3a_upset_5celltypes`)**:  
  - *What it is*: An UpSet matrix (an advanced bar-code version of a Venn diagram) that quantifies every possible intersection combination cleanly.
- **Fig3b & Fig3c (`common_genes_kegg` & `common_genes_go_bp`)**:  
  - *What it is*: Functional enrichment of the convergent genes. Confirms that shared stress involves MAPK signaling, blood–brain barrier breakdown, and protein processing.
- **Fig3d (`Fig3d_common_genes_ppi_network`)**:  
  - *What it is*: A circular hub network. Central node colors indicate connectivity (degree centrality). Points out the central "boss" genes that coordinate the response.
- **Fig3e (`Fig3e_cortex_expression_violins`)**:  
  - *What it is*: Violin plots showing the exact expression levels of our top candidates (*DUSP1, FOS, JUN, EGR1, GADD45B, ATF3*) in Control (blue) vs AD (red) across all 5 lineages.
- **Fig3g (`Fig3g_consensus_score_ranking`)**:  
  - *What it is*: The official leaderboard. A horizontal bar chart ranking the top 15 candidate targets based on our unbiased mathematical CTPS score.

---

### 📁 Figure 4: Intercellular Communication (`Fig4_Intercellular_Communication/`)
*Goal: Eavesdrop on the cellular telephone calls between cell types using CellChat v2.2.*

- **Fig4a (`Fig4a_intercellular_interaction_strength`)**:  
  - *What it is*: Side-by-side comparison of total interaction numbers and communication strength. Shows how AD remodels the overall volume of cellular crosstalk.
- **Fig4b (`Fig4b_differential_interaction_network`)**:  
  - *What it is*: A circular telephone network showing net changes. **Red lines** indicate communication channels that became abnormally hyper-activated in AD; **blue lines** show healthy communication lines that went dead.
- **Fig4c (`Fig4c_signaling_information_flow_ranking`)**:  
  - *What it is*: A ranked list of 27 molecular communication pathways. Shows which conversation topics changed the most (e.g., *APP* amyloid signaling surging in AD).
- **Fig4d (`Fig4d_nvu_pathway_circos_network`)**:  
  - *What it is*: A circular chord diagram tracing *APP* (amyloid) signaling specifically. Shows that astrocytes and neurons broadcast massive amounts of APP ligand directly onto microglia and endothelial receptors.
- **Fig4e (`Fig4e_ligand_receptor_communication_bubble`)**:  
  - *What it is*: A high-resolution bubble plot pinpointing the exact molecular "keys and locks" (ligand–receptor pairs like *APP–CD46* and *APP–SORL1*) driving disease communication.

---

### 📁 Figure 5: Target Translation & Drug Roadmap (`Fig5_Target_Validation/`)
*Goal: Turn computational discoveries into real-world pharmacology and laboratory validation.*

- **Fig5a (`Fig5a_target_cross_species_conservation`)**:  
  - *What it is*: Evaluates whether human and mouse versions of our top targets have identical protein sequences.  
  - *Key takeaway*: Targets like *SORL1* (93.4%), *PCDH9* (96.1%), and *DUSP1* (95.8%) show extraordinary cross-species identity, guaranteeing that drugs developed for them can be tested in 5xFAD and APP/PS1 mouse models.
- **Fig5b (`Fig5b_target_druggability_and_surfaceome`)**:  
  - *What it is*: Categorizes whether each target sits on the outer cell membrane (accessible to antibody drugs) or inside the cytoplasm/nucleus (requiring small-molecule pills).
- **Fig5c (`Fig5c_target_ppi_mechanistic_network`)**:  
  - *What it is*: The master disease blueprint. A comprehensive network physically connecting our prioritized targets (*SORL1, CLEC5A, ADAMTS9, PCDH9, DUSP1*) directly to the core hallmarks of Alzheimer's: Amyloid-beta production (*APP, BACE1*), microglial inflammation (*TREM2, TYROBP*), and blood vessel leakage (*CLDN5*).

---

## 🧪 Wet-Lab Validation: Recommended Antibodies for Collaborators

To facilitate immediate wet-lab follow-up by immunology and molecular biology team members, below are pre-validated research antibodies matching our prioritized targets:

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

All figures in this directory follow the exact same visual design rules:
1. **Typography**: Publication-standard sans-serif typography (`sans` / Arial / Helvetica) at $\ge 10\,\text{pt}$ for primary text.
2. **Standard Color Palettes**:
   - Condition: Control = Navy Blue (`#4575b4`), AD = Red (`#d73027`).
   - Lineages: Astrocyte (`#3C5488`), Endothelial (`#4DBBD5`), Inhibitory neuron (`#8491B4`), Microglia (`#E64B35`), Oligodendrocyte (`#F39B7F`), Pericyte (`#00A087`).
3. **Identical Symmetry with Hippocampus**: Panel layouts, aspect ratios, line widths, transparency levels, and 5-way Venn diagram parameters match the Hippocampus suite 1-to-1, allowing effortless side-by-side comparison.
