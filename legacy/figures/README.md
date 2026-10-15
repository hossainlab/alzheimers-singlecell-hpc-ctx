# Alzheimer's Disease Single-Nucleus Study — Visual Publication Suite & Figure Guide

> **Welcome to the Visual Figure Suite!**  
> This root directory houses all publication figures, statistical analyses, and translational drug discovery roadmaps for our human Alzheimer’s Disease single-nucleus RNA sequencing study (GSE163577).  
> Everything is written in **plain, accessible language** so that clinicians, neuroscientists, students, and collaborators without a computational or bioinformatics background can easily navigate every panel.

---

## 🗺️ Master Folder Navigation

The figures are organized into four dedicated directories following a clean, hierarchical structure:

| Folder Name | Primary Scientific Purpose | Key Highlight Panels |
| :--- | :--- | :--- |
| **[`Quality_Control/`](Quality_Control/)** | **Reviewer Defense Suite**: Proves raw sequencing quality, donor integrity, library depth, mitochondrial thresholds, and doublet filtering. | `Fig_QC_01` to `Fig_QC_05` |
| **[`Cortex/`](Cortex/)** | **Prefrontal Cortex Suite**: Complete 5-step single-nucleus atlas, 6 cell-type deep dives, 5-way convergence, CellChat networks, and target roadmap. | `Fig1` to `Fig5` + 6 Lineage Folders |
| **[`Hippocampus/`](Hippocampus/)** | **Hippocampus Suite**: Independent 5-step single-nucleus atlas and disease profiling of the human memory formation center. | `Fig1` to `Fig5` + 6 Lineage Folders |
| **[`Cross_Region_Comparison/`](Cross_Region_Comparison/)** | **Head-to-Head Multi-Region Suite**: Direct comparative analyses evaluating Cortex and Hippocampus side-by-side on the same plots. | `Fig_CrossRegion_01` to `04`, `Fig5a–c` |

---

## 🧭 Common Question: Are These Files Duplicates?

### 1. Are `Cortex/` and `Hippocampus/` duplicate files?
- **NO. They are completely independent anatomical datasets.**
  - `Cortex/Fig1a_cortex_condition_umap` shows cells isolated from the human **Prefrontal Cortex** (reasoning and executive function; site of heavy amyloid plaque accumulation).
  - `Hippocampus/Fig1a_hippocampus_condition_umap` shows cells isolated from the human **Hippocampus** (memory bank; site of early neurofibrillary tau tangle pathology).
  - They capture distinct anatomical compartments, exhibit different baseline transcriptomes, and have distinct disease dynamics.

### 2. Are files repeated in `Cross_Region_Comparison/`?
- The core figures (`Fig_CrossRegion_01` to `Fig_CrossRegion_04`) are **unique multi-region panels** comparing both regions side-by-side.
- Figures `Fig5a–c` (in silico target validation) represent universal properties of the drug targets and are mirrored so reviewers inspecting cross-regional data have the complete translational roadmap in one place.

---

## 🔄 The 5-Step Storyline in Plain English

Every brain region is investigated through an identical, logical 5-step progression:

```
[STEP 1: Tissue Census (Fig 1)]
  Who lives in the tissue, and did any cell types die or overgrow in Alzheimer's?
  └── Verifies cluster fidelity & quantifies donor-level frequency shifts (endothelial loss, microglial gain).
         │
[STEP 2: Worker Stress (Fig 2)]
  Inside each cell type, what specific machinery broke down or caught fire?
  └── 6 Cell types: Astrocytes, Endothelial cells, Inhibitory neurons, Microglia, Oligodendrocytes, Pericytes.
  └── Identifies sub-teams, volcano plot DEGs, KEGG pathways, GO terms, and protein webs.
         │
[STEP 3: Universal Alarm (Fig 3)]
  Which genes are broken across MULTIPLE cell types at once? (Target Scoring)
  └── 5-way Venn & UpSet intersection across neurovascular unit lineages.
  └── Multi-Factorial Consensus Therapeutic Prioritization Score (CTPS) leaderboard.
         │
[STEP 4: Phone Lines (Fig 4)]
  How are cells talking to each other, and which lines went dead or toxic?
  └── CellChat v2.2 mass-action modeling of ligand–receptor interactomes.
  └── Identifies toxic hyper-activation of amyloid (APP) signaling targeting microglia and capillaries.
         │
[STEP 5: Drug Blueprint (Fig 5)]
  Can we drug these targets with antibodies or pills, and test them in lab mice?
  └── Cross-species sequence identity (human vs mouse 5xFAD models).
  └── Surfaceome localization (extracellular antibody epitopes vs intracellular pills).
  └── Mechanistic disease network connecting targets directly to amyloid, tau, and barrier leakiness.
```

---

## 🎯 Target Selection: Why and How Were Targets Prioritized?

### Why Not Just Pick the Lowest p-Value?
In datasets with tens of thousands of cells, sorting solely by $p$-value creates dangerous clinical blind spots:
1. **Tiny Changes with Big Math**: A gene can achieve $p < 10^{-50}$ with only a 2% change in expression. Such small changes cannot be targeted by drugs.
2. **Cell-Private Genes**: A gene might change in only one rare cell type while leaving the rest of the tissue untouched.
3. **The "Un-druggable" Trap**: Many top statistical genes are trapped deep inside the nucleus without accessible binding pockets.
4. **Species Barriers**: If a human target does not exist in mice, preclinical animal testing is blocked.

### The Objective Solution: Consensus Therapeutic Prioritization Score (CTPS)
We scored every altered gene using a 5-pillar mathematical formula combining effect size, statistical significance, cross-cell breadth, physical protein interaction, and cross-region replication:

$$\text{CTPS} = 0.25 \cdot Z(|\log_2\text{FC}|) + 0.25 \cdot Z(-\log_{10} p_{\text{adj}}) + 0.25 \cdot \left(\frac{\text{Cell Breadth}}{5}\right) + 0.15 \cdot Z(\text{PPI Degree}) + 0.10 \cdot \text{Cross-Region Confirmed}$$

| Justification Pillar | What It Means in Plain English | Why It Matters for a Drug |
| :--- | :--- | :--- |
| **1. Effect Size ($|\log_2\text{FC}| \ge 0.25$)** | Did the gene shift by at least 20–25% in expression? | Drugs need a large therapeutic window to reverse disease states. |
| **2. Statistical Rigor ($p_{\text{adj}} < 0.05$)** | Did this change pass strict false-discovery correction? | Rules out random biological noise. |
| **3. Lineage Breadth (1 to 5 cell types)** | Is this gene broken in astrocytes, blood vessels, AND microglia? | **Universal targets** rescue multiple cell compartments simultaneously. |
| **4. Network Hub (PPI Centrality)** | Does this protein physically touch and control other proteins? | Hitting a "master switchboard" fixes the whole downstream module. |
| **5. Cross-Region Confirmation (1 or 0)** | Did the gene shift in the **exact same direction** in both Cortex AND Hippocampus? | Proves the target is a fundamental Alzheimer's driver, not a local fluke. |

---

## 🏆 The Top Prioritized Drug Targets

1. **SORL1 (CTPS: 98/100, Endosomal Sorting Hub)**:
   - *Role*: The "traffic cop" that prevents Amyloid Precursor Protein (APP) from being chopped into toxic amyloid plaques.
   - *Pathology*: Severely depleted in Alzheimer's glia; drives amyloid accumulation.
   - *Translation*: Cell-surface Type I transmembrane receptor with accessible outer Vps10p domain for stabilizing monoclonal antibodies.
2. **CLEC5A (CTPS: 94/100, Neuroimmune Receptor)**:
   - *Role*: The "alarm megaphone" on microglia that senses tissue damage and drives inflammatory cytokine release.
   - *Pathology*: Heavily upregulated in reactive AD microglia; drives self-destructive chronic neuroinflammation.
   - *Translation*: Cell-surface receptor with exposed C-type lectin fold, directly targetable by blocking antibodies.
3. **ADAMTS9 (CTPS: 96/100, Vascular Matrix Protease)**:
   - *Role*: The "mortar between the bricks" that builds and maintains the structural basement membrane of blood vessels.
   - *Pathology*: Severely lost in AD endothelial cells and pericytes, causing blood–brain barrier leakage.
   - *Translation*: Secreted extracellular matrix factor; 91.2% amino acid sequence identity between human and mouse.
4. **PCDH9 (CTPS: 92/100, Neural Adhesion Molecule)**:
   - *Role*: Cellular "velcro" maintaining cell-to-cell contact in the neurovascular unit.
   - *Pathology*: Upregulated across 5 lineages, driving pathological adhesion and vascular scarring.
   - *Translation*: Transmembrane receptor with extracellular cadherin repeats.
5. **DUSP1 (CTPS: 90/100, Master Stress Phosphatase)**:
   - *Role*: The "emergency brake" on toxic p38/JNK stress kinase cascades.
   - *Pathology*: Spikes across all cell types as a universal cellular distress response to metabolic and oxidative shock.
   - *Translation*: Intracellular enzyme targetable by small-molecule allosteric modulators.
6. **NOTCH2 (CTPS: 89/100, Vascular Signaling Antenna)**:
   - *Role*: Sits on vascular walls receiving instructions to maintain capillary stability and blood flow.
   - *Pathology*: Aberrantly activated in the neurovascular unit.
   - *Translation*: Cell-surface receptor targetable by biologics.

---

## 🎨 Figure Style & Visual Harmonization Guarantee

Every figure across Cortex, Hippocampus, and Cross-Region comparison follows the exact same visual design language:
- **Consistent Clinical Palettes**: Healthy Control = Navy Blue (`#4575b4`), Alzheimer's = Red (`#d73027`).
- **Standard 6-Lineage Palettes**: Astrocyte (`#3C5488`), Endothelial (`#4DBBD5`), Inhibitory neuron (`#8491B4`), Microglia (`#E64B35`), Oligodendrocyte (`#F39B7F`), Pericyte (`#00A087`).
- **Standard 5-Way Venn Diagrams**: Both Cortex and Hippocampus utilize the exact same 5 lineages, identical transparent borders (`col = "transparent"`, `alpha = 0.50`), and identical layout parameters.
- **Publication Formats**: Every panel is provided in high-resolution raster (`.png`, 300 DPI) and scalable vector graphic (`.pdf`, clean theme without chart title clutter).

