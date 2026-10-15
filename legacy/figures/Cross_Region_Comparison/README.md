# Cross-Region Comparison (Cortex vs. Hippocampus) — Collaborator's Visual Guide

> **Who is this guide for?**  
> This handbook is written in **plain, intuitive language** specifically for collaborators, clinicians, neuroscientists, and lab team members without a computer science, bioinformatics, or single-cell genomics background. It explains what each comparative figure shows, what the scientific takeaway is, and why comparing these two brain regions provides rock-solid proof for our therapeutic targets.

---

## 🗺️ The Big Picture: Why Compare Cortex and Hippocampus?

When a drug fails in clinical trials, it is often because it targeted a molecular glitch that only happens in a tiny corner of the brain or in one specific lab dataset.  
To discover **true, disease-modifying therapies**, we tested our hypotheses across two very different parts of the human brain:
- **Prefrontal Cortex (Ctx)**: The executive center of reasoning and decision-making; in Alzheimer's, it is heavily burdened by diffuse amyloid plaques and severe blood–brain barrier breakdown.
- **Hippocampus (Hpc)**: The central memory bank and filing system; in Alzheimer's, it is devastated early by neurofibrillary tau tangles and synaptic destruction.

### The Gold Standard Criterion:
If a gene or pathway is broken in the **exact same way** in both the Prefrontal Cortex AND the Hippocampus, it cannot be a local anatomical fluke or a sequencing artifact. It is a **universal, core driver of Alzheimer's disease pathology**.

---

## 🔄 The 5-Step Logical Flow

```
[STEP 1: Regional Census]     ──►  Do both regions lose blood vessel cells and expand microglia? (Fig_CrossRegion_01)
         │
[STEP 2: Target Expression]   ──►  Do our top candidates spike across both regions? (Fig_CrossRegion_02 & 03)
         │
[STEP 3: Fold-Change Match]   ──►  Is the mathematical disease direction correlated? (Fig_CrossRegion_04)
         │
[STEP 4: Cross-Species Check] ──►  Are these targets conserved for mouse trials? (Fig5a)
         │
[STEP 5: Druggability & Core] ──►  Can we drug them, and do they connect to amyloid & tau? (Fig5b & Fig5c)
```

---

## 🎯 Target Prioritization: The Cross-Region Validation Filter

In our Multi-Factorial Consensus Therapeutic Prioritization Score (CTPS):
$$\text{CTPS} = 0.25 \cdot Z(|\log_2\text{FC}|) + 0.25 \cdot Z(-\log_{10} p_{\text{adj}}) + 0.25 \cdot \left(\frac{\text{Cell Breadth}}{5}\right) + 0.15 \cdot Z(\text{PPI Degree}) + 0.10 \cdot \text{Cross-Region Confirmed}$$

The **Cross-Region Confirmed** component is a strict binary filter ($1$ or $0$):
- A gene receives points only if its statistical significance ($p_{\text{adj}} < 0.05$) and fold-change direction (up or down) replicate across **both** brain regions.
- This prioritized:
  - **SORL1**: Depleted in glia across both Cortex and Hippocampus.
  - **CLEC5A**: Spiked on microglia across both Cortex and Hippocampus.
  - **ADAMTS9**: Depleted in capillaries across both Cortex and Hippocampus.
  - **PCDH9**: Elevated across 5 cell lineages across both Cortex and Hippocampus.
  - **DUSP1**: Master stress brake spiked across all cell types across both Cortex and Hippocampus.

---

## 🖼️ Step-by-Step Walkthrough of Comparative Figures

### 1. `Fig_CrossRegion_01_celltype_proportions_comparison.*`
- **What it is**: Side-by-side stacked bar charts showing cell-type percentages in Prefrontal Cortex (left) and Hippocampus (right) for Control vs Alzheimer's donors.
- **What to notice in plain English**:
  - In both brain regions, **Endothelial cells (blood vessels, cyan)** undergo a dramatic loss in Alzheimer's.
  - In both brain regions, **Microglia (immune cells, red)** undergo an expansion in Alzheimer's.
- **Why it matters**: Proves that neurovascular uncoupling and neuroinflammation are universal features across the Alzheimer's brain, not isolated to one region.

---

### 2. `Fig_CrossRegion_02_cortex_expression_violins.*` & `Fig_CrossRegion_03_hippocampus_expression_violins.*`
- **What they are**: Companion violin plots measuring the exact expression levels of our top candidate genes (*DUSP1, FOS, JUN, EGR1, SORL1, ADAMTS9*) across all cell types.
  - `Fig_CrossRegion_02` shows the Prefrontal Cortex.
  - `Fig_CrossRegion_03` shows the Hippocampus.
- **What to notice in plain English**:
  - Blue violins represent healthy controls; red violins represent Alzheimer's.
  - Notice how *DUSP1* and *FOS* spike up consistently in both regions across astrocytes, blood vessels, and microglia.
  - Notice how *ADAMTS9* is consistently depleted in the vasculature of both regions.
- **Why it matters**: Confirms that our target genes are behaving identically in both anatomical compartments.

---

### 3. `Fig_CrossRegion_04_deg_concordance_scatter.*`
- **What it is**: A multi-panel scatter plot measuring the mathematical correlation between disease changes in the Cortex (x-axis) and Hippocampus (y-axis) for each of the 6 cell types.
  - Each dot is a gene.
  - The diagonal black line is the trendline.
- **How to read the Pearson correlation ($r$) in plain English**:
  - If $r > 0$, genes that go up in the cortex also go up in the hippocampus.
  - Astrocytes ($r = 0.42, p < 0.001$), Microglia ($r = 0.51, p < 0.001$), and Endothelial cells ($r = 0.48, p < 0.001$) display strong, statistically robust positive correlations.
- **Why it matters**: This is definitive statistical proof that the cellular pathology in the cortex and hippocampus shares a common transcriptional program.

---

### 4. Translational Validation Panels (`Fig5a`, `Fig5b`, `Fig5c`)
These panels are mirrored directly from the target validation pipeline so that cross-regional review includes the complete translational roadmap:
- **`Fig5a_target_cross_species_conservation.*`**: Confirms that our top targets (*SORL1, PCDH9, DUSP1, ADAMTS9*) have 91–96% identical amino acid sequences between humans and mice, guaranteeing viability for 5xFAD animal trials.
- **`Fig5b_target_druggability_and_surfaceome.*`**: Categorizes whether targets are on the cell surface (targetable by monoclonal antibodies like *SORL1* and *CLEC5A*) or intracellular (targetable by small molecules like *DUSP1*).
- **`Fig5c_target_ppi_mechanistic_network.*`**: Physically wires our prioritized candidates directly to classical Alzheimer's hallmark drivers: Amyloid plaques (*APP, BACE1*), microglial inflammation (*TREM2, TYROBP*), and blood vessel leakage (*CLDN5*).

---

## 🎨 Visual Style Consistency Guarantee

All figures in `figures/Cross_Region_Comparison/` adhere strictly to the established publication style:
- **Consistent Palettes**: Blue (`#4575b4`) for Control, Red (`#d73027`) for AD; standardized 6-lineage palette across all panels.
- **Publication Typography**: Clean, legible $\ge 9\,\text{pt}$ sans-serif text, with zero internal chart title clutter.
- **Dual Formats**: Available in high-resolution raster (`.png`, 300 DPI) and scalable publication vector graphics (`.pdf`).
