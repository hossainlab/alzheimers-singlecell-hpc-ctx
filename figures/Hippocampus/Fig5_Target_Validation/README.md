# Figure 5: In Silico Target Validation & Drug Translation (Hippocampus)

> **Plain-Language Summary for Collaborators**  
> Finding a target on a computer screen is only the first step. To make a clinical difference in Alzheimer's Disease, a target must be **testable in animal models**, **druggable with real medicines (antibodies or pills)**, and **mechanistically tied to the core causes of hippocampal neurodegeneration**.  
> Figure 5 bridges computational discovery to wet-lab pharmacology and preclinical validation.

---

## 🖼️ Panel Walkthrough (What You Are Looking At)

### 1. `Fig5a_target_cross_species_conservation.*` — Can We Test This in Lab Mice?
- **What it shows**: Evaluates whether human and mouse versions of our top prioritized target proteins have identical amino acid sequences.
- **In plain English**: Before testing a drug in human patients, it must be validated in preclinical Alzheimer's mouse models (such as 5xFAD or APP/PS1). High sequence identity guarantees that preclinical trials in mice will faithfully model human pharmacology.
- **Key findings**:
  - **PCDH9**: 96.1% sequence identity (98.0% functional domain match)
  - **DUSP1**: 95.8% sequence identity (97.4% functional domain match)
  - **SORL1**: 93.4% sequence identity (96.2% functional domain match)
  - **ADAMTS9**: 91.2% sequence identity (94.8% functional domain match)
  - All four achieve high **In-Vivo Suitability Scores (>90/100)**.

### 2. `Fig5b_target_druggability_and_surfaceome.*` — How Can We Deliver the Drug?
- **What it shows**: Classifies where each protein lives in the cell (the "Surfaceome").
- **In plain English**:
  - **Cell Surface / Secreted (Green)**: Proteins on the cell membrane (*SORL1, CLEC5A, PCDH9*) or secreted into the extracellular matrix (*ADAMTS9*). Directly targetable by **monoclonal antibodies**.
  - **Intracellular (Blue)**: Proteins located inside the cell (*DUSP1*). Targetable by **cell-permeable small-molecule pills**.

### 3. `Fig5c_target_ppi_mechanistic_network.*` — The Master Disease Blueprint
- **What it shows**: A comprehensive molecular circuit diagram physically connecting our prioritized targets (*SORL1, CLEC5A, ADAMTS9, PCDH9, DUSP1*) directly to classical Alzheimer's drivers:
  - **Amyloid-Beta Pathway**: Links *SORL1* to *APP* and *BACE1*.
  - **Neuroinflammation Pathway**: Links *CLEC5A* to *TREM2* and *TYROBP*.
  - **Blood–Brain Barrier Pathway**: Links *ADAMTS9* to *CLDN5*.
  - **Stress Kinase Pathway**: Links *DUSP1* to *MAPK14 (p38)* and *JUN*.

---

## 🧪 Wet-Lab Reagents for Immediate Testing
Validated research antibodies, catalog numbers, host species, and recommended dilutions for each prioritized target are cataloged in `tables/recommended_antibodies_and_compounds.csv` and detailed in the main `figures/Hippocampus/README.md`.
