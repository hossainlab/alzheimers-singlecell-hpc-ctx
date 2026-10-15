# Interpretation of the figures

This document explains what each figure shows, why each analysis was chosen, and how far each conclusion
can be pushed. Every number comes from a table in `tables/`; the source is given in brackets. Panel legends
are in [`figures/README.md`](../figures/README.md).

**Reading conventions**

- **Unit of inference.** The statistical unit is the *sample* (one donor × one brain region). Nuclei are never
  treated as replicates for AD-vs-control questions.
- **Hippocampus is the discovery set.** It has 8 control and 9 AD donors, all independent.
- **Cortex is not independent.** It has 4 control and 4 AD samples, all from donors who also gave hippocampus.
  It is a second tissue from the same people and cannot replicate the hippocampal results.
- **"Robust" DEG** means FDR < 0.05, |log₂FC| ≥ 0.5, and no artefact flag (ambient carry-over, 17q21.31 locus,
  or sex chromosome). The model is pseudobulk DESeq2, `~ sex + vascular_yield + condition`.

---

## Figure S1. Quality control (`figures/Quality_Control/FigS1_quality_control`)

**What it shows.**
- 319,449 nuclei were loaded and 269,956 kept (84.5 %). Of those removed, 18,863 were low quality and 30,630
  were doublets [`qc_sample_summary.csv`].
- Retention per sample ranged from 56 % (Hpc 01_06 control, 3,859 low-quality nuclei) to 97 %.
- Median genes per nucleus did not differ by condition (Mann–Whitney P = 0.37 in hippocampus, 0.89 in cortex).
- Inferred sex: 8 female and 9 male donors (AD 5 F / 4 M; control 3 F / 5 M), so sex is not confounded with
  condition.
- In every one of the 8 donors who gave both tissues, the sex inferred from hippocampus and cortex agrees. The
  ARL17B (17q21.31) genotype class agrees in 7 of 8 [`donor_pairing_concordance.csv`].

**Why these analyses.**
- Thresholds are MAD-based within each sample, not fixed. GEO libraries differ more than twofold in depth, and
  a fixed cut-off would remove whole samples unevenly.
- Doublets are called per sample, because 10x doublet rates depend on how many nuclei were loaded in each run.
- Sex and genotype checks test the assumption that cortex and hippocampus files with the same donor code come
  from the same person. This decides whether the cortex can be used as replication, and the answer is no.

**What can be concluded.**
- Library quality is comparable between AD and control samples, so it cannot explain the DE results.
- ARL17B expression is either ~10–60 or ~150–400 CPM, depending on donor. That bimodality reflects a germline
  copy-number polymorphism, not a disease state. It explains why ARL17B looked like the "top target" in the
  legacy analysis: in cortex, 3 of 4 AD donors fall in the high class and 0 of 4 controls do.
- One pair (04_09) disagrees on genotype class. That could be a mislabel or a region-specific expression
  difference, and it should be checked against the original donor key.


---

## Figure S2. Donor clinical covariates (`figures/Quality_Control/FigS2_donor_covariates`)

**What it shows.**
- **a, b.** Age and post-mortem interval of the 17 donors, from Yang *et al.* Supplementary Table 1. Neither
  differs between AD and control. Circles mark the one donor per group with acute hypoxic neuropathology.
- **c.** The pericyte hypoxia-module score does not track post-mortem interval (ρ = −0.12, P = 0.63).
- **d.** The AD − control difference in the module score for each NVU cell type, under seven models: the
  primary model plus age, PMI, APOE ε4 dose, hypoxic pathology, age + PMI, and exclusion of the two
  hypoxic-pathology donors. Every estimate is negative and every 95 % CI excludes zero.
- **e.** Re-fitting the hippocampal pseudobulk DESeq2 with age + PMI keeps most robust DEGs significant.
  Excluding the two hypoxic donors loses power, but the sign of every primary robust DEG is unchanged in
  both re-fits.

**Why these analyses.** The main finding is a lower hypoxia program in AD. The leading alternative
explanation is that the groups differed in how they died. The supplementary table is the only donor-level
clinical information available, so every usable variable was tested two ways:
1. one covariate at a time against a single summary score, which keeps power at n = 17;
2. in the full gene-level model, which checks that individual genes hold up.

**What can be concluded.** No recorded covariate explains or weakens the AD effect, and the covariates are
balanced between groups. The remaining uncertainty is limited to variables that were never recorded: agonal
score, brain pH and RIN. These require the brain-bank records or a matched cohort
([`02_wet_lab_validation_plan.md`](02_wet_lab_validation_plan.md), Stage 1).
---

## Figure 1. Landscape of both regions (`figures/Fig1_Global_Landscape/Fig1_global_landscape`)

**What it shows.**
- **a–c.** 257,025 nuclei in 11 classes from both regions on one Harmony-integrated embedding. The largest
  classes are endothelial (73,202), oligodendrocyte (65,786), pericyte (38,166), astrocyte (37,740), smooth
  muscle (15,097) and microglia (9,164). Hippocampus and cortex nuclei mix within every cluster (b), so region
  is not driving cluster identity.
- **d.** Each class expresses its canonical markers. Examples: *CLDN5/FLT1/VWF* in endothelium;
  *PDGFRB/RGS5/ABCC9* in pericytes; *P2RY12/CX3CR1/C3* in microglia; *MOG/ST18/PLP1* in oligodendrocytes.
- **Removed clusters.** 12,931 nuclei in seven clusters were excluded [`cluster_annotation.csv`]. In each, a
  marker of a second lineage was detected in ≥ 30 % of nuclei and at ≥ 3× background. For example, MBP was
  detected in 80–86 % of nuclei in endothelial-, myeloid- and mural-labelled clusters. That pattern indicates
  heterotypic doublets or nuclei dominated by ambient RNA.
- **e, f.** Per-sample composition. AD samples contain far fewer vascular nuclei: hippocampus median 24 % vs
  62 % (P = 0.021); cortex 31 % vs 79 % (P = 0.029). The parenchymal share rises correspondingly; for example,
  oligodendrocytes are 44 % vs 15 % in hippocampus [`Fig1_composition_tests.csv`].

**Why these analyses.**
- **Joint embedding.** Integrating both regions in one embedding means every cell class is defined the same
  way in hippocampus and cortex. This is required before comparing regions.
- **Cluster-level annotation with an explicit doublet rule.** Marker-score annotation alone labelled
  heterotypic clusters as clean cell types. That was one source of the legacy "convergent" genes, such as
  *PLP1* and *MBP* appearing in every cell type.
- **Composition per donor.** Showing composition per sample, not pooled, lets the reader see that the spread
  between donors is as large as the difference between groups.

**What can be concluded.**
- VINE-seq enriches vessels by physical isolation. How well that works depends on the tissue, and AD tissue
  gave systematically less enrichment.
- Fig. 1e–f therefore describes **technical yield**, not loss of vascular cells in AD. No proportion test
  survives FDR < 0.05 after correction across classes (lowest FDR 0.079).
- The yield difference has two consequences that the rest of the pipeline handles:
  1. *Ambient background.* AD libraries contain more parenchymal ambient RNA. Each DE gene is therefore tested
     for whether it tracks the source cell type's abundance, and flagged if it does.
  2. *Confounding.* Yield is partly confounded with condition, so the logit of each sample's vascular fraction
     is a covariate in every DE model.

---

## Figure 2a–e. NVU cell-type deep dives (`figures/Fig2_Cell_Type_Deep_Dives/`)

**Common logic, all five figures.**
- **Sub-clusters (panels a–c)** are descriptive. Their markers come from a within-cell-type Wilcoxon test,
  which is acceptable for labelling states but not for AD inference.
- **Sub-cluster proportions (d)** are tested per sample, and only shifts with BH-FDR < 0.05 are marked.
- **Volcanoes and heat maps (e–g)** are donor-level.
- **Pathway analysis (h)** is pre-ranked GSEA on the DESeq2 Wald statistic with flagged genes removed. Ranking
  keeps information from genes below the significance cut-off, which matters most for the underpowered cortex.

### 2a Microglia (9,164 nuclei, 8 sub-clusters)
- **DE.** Hippocampus has 29 robust DEGs (4 up, 25 down). The down-regulated genes include *SLC2A3*
  (GLUT3), *SLC11A1*, *FOSL2*, *ADAMTS9*, *HIF1A-AS2* and *ADGRE2*; *HDAC9* and *PDGFC* are up.
- **Pathways.** Hallmark TNF-α/NF-κB (NES −2.37), Hypoxia (−2.11) and Inflammatory response (−1.93) are lower
  in AD.
- **Sub-clusters.** No sub-cluster shift reaches FDR < 0.05. Two trends sit just above it (FDR 0.09): a
  *PCDH9/IL1RAPL1* group (S7) and an *FLT1/MECOM* group (S4). Both carry oligodendrocyte or endothelial marker
  transcripts, so they may be partly ambient.
- **Cortex.** Only *F13A1* and *CD163* are down. Both are border-associated macrophage markers, which more
  likely reflects how many perivascular macrophages were captured than a change in microglial state.
- **Interpretation.** There is no evidence for a disease-associated microglia (DAM) activation signature in
  this vessel-enriched material. The dominant change is *lower* expression of inducible stress and
  inflammatory genes.

### 2b Endothelial cells (73,202 nuclei, 8 sub-clusters)
- **DE.** Hippocampus has 118 robust DEGs (29 up, 89 down). Down: *MYOF, HS3ST2, PYGL, ADAMTS9, ELL2, PCK1,
  ATP1A1*. Up: *SLC39A10, TNFSF10, GPR85, KLF12, SLCO1A2*.
- **Pathways.** Hypoxia (NES −2.40), TNF-α/NF-κB (−2.11), IL-4/IL-13 and IL-10 signalling are lower in AD.
- **Sub-clusters.** One shift passes FDR: the *TSHZ2/PDE4D* sub-cluster (S2) is slightly smaller in AD
  (22.8 % vs 26.0 %; FDR = 0.039).
- **Cortex.** 20 robust DEGs, nearly all down: *SELE, NAMPT, ACKR1* (venular), *ALOX5, LIPG*. Cortical GSEA is
  led by lower ribosome biogenesis and Myc targets.
- **Interpretation.** Endothelial cells show the same lower hypoxia and inflammatory tone. Lower *SELE* and
  *ACKR1* in cortex points to less venular activation in AD. Neither change fits a simple "BBB breakdown"
  account, and neither can be separated from terminal state.

### 2c Pericytes (38,166 nuclei, 4 sub-clusters): the largest effect
- **DE.** Hippocampus has 360 robust DEGs (150 up, 210 down), the largest of any cell type. Down: *SLC6A17,
  PCBP3, STC1, CEBPD, HIPK2, ADAMTS9, IGFBP4*. Up: *GRID2, SESN3, PDP1, SLC20A2, CYP7B1*.
- **Pathways.** TNF-α/NF-κB (NES −2.83), IL-6/JAK/STAT3 (−2.30), Inflammatory response (−2.20), Hypoxia
  (−2.18) and Angiogenesis (−2.16) are lower in AD.
- **Sub-clusters.** The *ADAMTS9/ELL2/PCBP3*-high state (S2) is a median 41.4 % of pericytes in control
  hippocampus and 2.8 % in AD (FDR = 0.044). The *ATP1A2/SLC20A2/LAMA2* homeostatic state (S1) rises
  reciprocally, from 41.7 % to 78.3 %. Cortex shows the same direction (35.8 % → 4.4 %), but it is not
  significant with 4 vs 4 samples.
- **Interpretation.** The pericyte DE is largely this state switch. S2 is defined by the same inducible genes
  that make up the convergent program. Yang *et al.* reported vulnerability of ECM-maintaining pericytes in AD;
  that result and this one may be two views of the same shift. *SLC20A2* (higher in AD) is a
  primary-familial-brain-calcification gene, so pericyte phosphate handling is a secondary lead.

### 2d Astrocytes (37,740 nuclei, 8 sub-clusters)
- **DE.** Hippocampus has 42 robust DEGs (9 up, 33 down). Down: *GLUL, ALDH1L1, VEGFA, HILPDA, F3, ENO1*. Up:
  *SLC39A10, LRP6, SEMA3A*.
- **Pathways.** Hypoxia (NES −2.72, the strongest of any contrast) and glycolytic process (−2.36) are lower in
  AD. Cilium movement and axoneme assembly are higher; these probably come from ependymal-adjacent astrocytes.
- **Sub-clusters.** No shift reaches FDR < 0.05.
- **Cortex.** Seven robust DEGs (up: *DACH1, HES5, CAMK4*; down: *ADAMTS9, CPAMD8, SLC38A1*). Cortical GSEA
  shows lower HSF1/heat-shock and interferon-γ responses.
- **Interpretation.** Lower *GLUL* (glutamine synthetase) and *ALDH1L1* would, if confirmed at protein level,
  indicate reduced astrocytic metabolic support. The co-regulated loss of *VEGFA, HILPDA, ENO1* again matches
  a lower HIF-1 state.

### 2e Oligodendrocytes (65,786 nuclei, 4 sub-clusters)
- **DE.** Hippocampus has 65 robust DEGs (15 up, 50 down). Down: *ADAMTS9* (the most significant gene), *SETD7,
  SH3PXD2B, KLF9, MGAT5*. Up: *PCDH10, DCC, SCN8A, SESN3*.
- **Pathways.** Hypoxia (NES −2.13) is lower. Ribosomal and translation-elongation sets are higher (NES
  +2.5 to +2.6).
- **Cortex.** No robust DEGs, but GSEA shows the opposite ribosomal direction. A regional difference that
  reverses sign is a warning sign of technical, library-level effects.
- **Interpretation.** Oligodendrocytes are the most common nucleus type in AD samples because of low vessel
  yield. Their AD-specific signal is small apart from the shared hypoxia module.

---

## Figure 3. Convergence and concordance (`figures/Fig3_Cross_Cell_Convergence/Fig3_convergence`)

**What it shows.**
- **a.** Few significant genes are removed by the artefact flags once the yield covariate is in the model:
  4 of 124 in hippocampal endothelium and 7 of 381 in pericytes. Without the covariate, the legacy-style result
  of hundreds of "shared" genes would reappear (sensitivity table).
- **b.** Most robust DEGs are specific to one cell type: 322 to pericytes only, 93 to endothelium only, 57 to
  oligodendrocytes only.
- **c.** AD effects are correlated between cell types (Spearman ρ 0.22–0.50). The strongest correlation is
  endothelium–pericyte (0.50) [`celltype_lfc_correlation_hippocampus.csv`].
- **d.** 47 genes are robust in ≥ 2 hippocampal cell types, and 39 of those are lower in AD
  [`convergent_genes.csv`]. The leading genes are *ADAMTS9* (4 cell types; median log₂FC −3.3), *VEGFA, ELL2,
  YBX3, BCL6, HIF1A* (3 each), then *HILPDA, PDK1, SLC2A3, BNIP3L, LAMB1, COL4A2, STAT3, SMAD7*.
- **e.** For hippocampal robust DEGs, the same donors' cortex has the same sign in 76.9–87.0 % of genes per cell
  type (binomial P from 8.8 × 10⁻³⁷ in pericytes to 5.5 × 10⁻⁴ in microglia)
  [`cross_region_concordance.csv`].
- **f.** No immediate-early gene (*FOS, FOSB, JUN, JUNB, EGR1, DUSP1, GADD45B, ATF3, IER2, NR4A1*) is
  significant in more than one contrast.

**Why these analyses.**
- **Panel a** was added to show how much of the legacy "convergence" came from composition rather than biology.
- **The UpSet plot (b)** shows that specificity, not convergence, is the rule.
- **Correlation (c)** measures sharing without depending on a significance threshold.
- **Panel f** directly tests the legacy README claim of a universal immediate-early-gene program. That claim is
  not supported.

**What can be concluded.**
- **The shared module is coherent.** *VEGFA, HILPDA, PDK1, SLC2A3, BNIP3L, ADAMTS9, HIF1A* and *ELL2* are all
  hypoxia-inducible. In the cell types where it is detected, this module is **lower in AD** across the
  vascular, glial and myeloid compartments.
- **It is not a yield artefact.** A per-sample hypoxia score in hippocampal pericytes (mean z of *ADAMTS9,
  VEGFA, HILPDA, SLC2A3, PDK1, BNIP3L, ELL2, HIF1A*) separates the groups (Mann–Whitney P = 0.001). The
  difference holds after adjusting for yield and sex (OLS P = 0.0004 for condition; P = 0.26 for yield)
  [`hypoxia_module_score_hpc.csv`]. Low-yield controls (01_06, 04_09) score high and high-yield AD
  samples (04_11, 01_05) score low.
- **Donor covariates do not explain it (Fig. S2).** Yang *et al.* Supplementary Table 1
  (`data/41586_2021_4369_MOESM2_ESM.xlsx`) maps onto the 17 donors with 17/17 sex and 17/17 cortex-subset
  agreement [`donor_clinical_metadata.csv`].
  - *Balance.* AD and control donors do not differ in age (median 84 vs 78 y, P = 0.15), post-mortem
    interval (5.8 vs 5.0 h, P = 0.77), APOE ε4 dose (median 1 vs 1), acute hypoxic neuropathology (1 donor
    each) or noted infarcts (3 vs 1, P = 0.58) [`covariate_balance.csv`].
  - *No association with the score.* The module score is not associated with age, PMI or APOE ε4 dose in
    any cell type (|ρ| ≤ 0.23, all P ≥ 0.38).
  - *The AD effect survives every adjustment.* In pericytes the AD effect is −1.44 (95 % CI −2.10 to −0.78,
    P = 0.0004). It is −1.48 (P = 0.0007) after adding age and PMI, and −1.45 (P = 0.004) after removing the
    two donors with acute hypoxic neuropathology. Every cell type stays significant
    [`hypoxia_module_covariate_models.csv`].
  - *Gene-level DE is stable.* Re-fitting with age and PMI keeps 43–86 % of hippocampal robust DEGs at
    FDR < 0.05, with gene-wide log₂FC correlation 0.95–0.99 and no sign changes. Removing the two
    hypoxic-pathology donors (n = 15) costs power (2–48 % retained) but changes no sign
    [`pseudobulk_DE_sensitivity_clinical.csv`].
  - *Informative single cases.* The control donor with hypoxic-ischaemic encephalopathy (NCI 7) has a high
    score (+0.90), as agonal hypoxia would predict. The AD donor with acute hypoxic changes (AD 1) is still
    low (−0.88). This fits a blunted response in AD, but it rests on a single donor per group.
- **What remains open.** The table reports no agonal-factor score, brain pH or RIN. A neuropathological note of
  acute hypoxic damage is an insensitive proxy for terminal hypoxia, so an agonal confound is *reduced* but
  not *excluded*. Braak stage is reported for only 5 of the 9 AD donors (III, IV, VI, VI, VI), which is too few
  for a dose–response test.
- **Cortical agreement is consistency, not replication.** It shows the module is a property of the donor, not of
  one tissue block, which fits *both* an AD effect and an agonal effect.

---

## Figure 4. Cell–cell communication (`figures/Fig4_Intercellular_Communication/Fig4_communication`)

**What it shows.**
- **a.** The number of expressed ligand–receptor pairs per sample does not differ by condition. Hippocampus
  medians are 1,826 (control) vs 686 (AD), P = 0.48, because the spread between donors is very wide. Cortex is
  1,668 vs 1,821, P = 0.89.
- **b, c.** Of 1,638 hippocampal and 2,476 cortical interactions tested, none reaches FDR < 0.05; the lowest
  FDR values are 0.68 and 0.997. Nominal P < 0.01 occurs for 15 and 7 interactions, about what chance predicts
  for this many tests.
- **d.** The legacy headline, APP→SORL1 hyper-signalling into microglia, is not supported (FDR 0.81 in
  hippocampus, 1.0 in cortex).
- **e.** The top nominal hits are mostly pericyte/endothelial basement-membrane–integrin pairs (*COL4A1/2,
  LAMC1, LAMA2/3 → ITGA7/ITGB1, ITGAV/ITGB8*). Their cortical estimates have wide intervals.

**Why this design.**
- CellChat, as used before, compared one pooled object per condition. Its permutation P values test labels
  across nuclei, not across donors, which gave apparent significance from a single comparison.
- Here each interaction is scored per sample and tested across samples with the same covariates as DE.
- APP–SORL1 is a CellChatDB pair missing from the LIANA consensus resource. It was added explicitly so the
  legacy claim could be tested on its own terms.

**What can be concluded.**
- With 17 hippocampal donors there is no detectable rewiring of NVU ligand–receptor communication.
- Any intercellular mechanism proposed for the hypoxia module (for example, lower VEGFA → KDR/FLT1 signalling)
  should be tested directly, not inferred from these scores. Expression-product scores cannot detect changes
  in protein secretion or binding.

---

## Figure 5. Target prioritisation (`figures/Fig5_Target_Prioritization/Fig5_targets`)

**What it shows.**
- **a.** The top 20 candidates, ranked first by the number of hippocampal cell types in which the gene is
  robust, then by cortex support, then by the smallest FDR.
- **b.** Computed annotation for each candidate. Examples:
  - *ADAMTS9* (UniProt Q9P2N4; secreted extracellular-matrix protease; human–mouse identity 88.9 %)
  - *VEGFA* (P15692; secreted; 83.8 %)
  - *HIF1A* (Q16665; nuclear; 90.7 %)
  - *SLC2A3* (P11169; cell membrane; 82.9 %)
  - *ADGRE2*: primate-specific, with no mouse orthologue.
- **c.** Per-donor expression of the top four candidates. The two groups barely overlap; *ADAMTS9* in
  oligodendrocytes is the clearest case.
- **d.** STRING v12 links *HIF1A–STAT3–BCL6–SLC2A3* (HIF1A–STAT3 score 0.998) and *COL4A2–LAMB1* (0.862).
  None of the 20 candidates connects to the 18 AD GWAS/reference genes at score ≥ 0.7. *VEGFA* cannot be
  resolved by symbol in STRING v12; querying the name returns COL18A1, so it was excluded rather than
  mis-mapped.

**Why this design.**
- The legacy score added z-scores of log-fold-change, −log₁₀ P (pseudoreplicated), breadth and a
  "PPI degree" that came from an expression-correlation network, with weights chosen by hand.
- A lexicographic ranking has no free weights and is easy to audit.
- Annotation values are retrieved and computed, not typed in, and the API responses are cached in
  `data/annotation/`.

**What can be concluded.**
- The candidates are really one biological module, not independent targets. Their ranking tells us which
  genes report the module most strongly and most widely, not which should be drugged.
- *ADAMTS9* is the most reproducible single gene: down in 4 of 5 hippocampal cell types and in cortical
  astrocytes. It is secreted, so it can be measured as protein in tissue or CSF.
- Because the module may reflect terminal state, none of these genes should be described as a therapeutic
  target until the validation in [`02_wet_lab_validation_plan.md`](02_wet_lab_validation_plan.md) separates
  disease from agonal effects.

---

## Overall conclusions

1. **Composition.** VINE-seq captures fewer vascular nuclei from AD tissue. This must be modelled; it should
   not be interpreted as cell loss.
2. **The main finding.** After donor-level testing and artefact control, the dominant transcriptional
   difference across the hippocampal NVU is a **lower hypoxia/HIF-1 and inflammatory-stress response in AD**.
   It is strongest in pericytes, where an *ADAMTS9/ELL2/PCBP3*-high state largely disappears.
3. **Robustness.** The difference is not explained by yield, sex, library quality, age, post-mortem
   interval, APOE ε4 dose or recorded acute hypoxic neuropathology. It is consistent in the same donors'
   cortex. Agonal state and brain pH are still unmeasured, so a terminal-hypoxia contribution cannot be fully
   excluded.
4. **Legacy claims not supported.** ARL17B and the lncRNA "targets", the universal immediate-early-gene
   program, APP→SORL1 hyper-signalling and NVU communication rewiring do not hold at donor level.
5. **Possible meaning if confirmed.** If the module survives agonal matching, it would suggest that
   AD-affected neurovascular cells respond less to hypoxic stress. That would be a testable mechanism linking
   AD to vascular vulnerability.
