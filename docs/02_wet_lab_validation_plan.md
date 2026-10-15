# Wet-lab validation plan

## Questions to answer

The computational analysis ([`01_figure_interpretations.md`](01_figure_interpretations.md)) supports one main
finding. **In AD hippocampus, the hypoxia/HIF-1 and inflammatory-stress transcriptional module is lower across
neurovascular-unit (NVU) cell types, and the effect is strongest in pericytes.** The module's lead genes are
*ADAMTS9, VEGFA, HILPDA, SLC2A3, PDK1, BNIP3L, ELL2* and *HIF1A*. Part of this is a near-complete loss of an
*ADAMTS9/ELL2/PCBP3*-high pericyte state.

The validation programme must answer three questions, in this order:

| # | Question | Why it comes first |
|---|---|---|
| Q1 | Does the difference hold when AD and control donors are **matched for agonal state, tissue pH and post-mortem interval (PMI)**? | If not, it is a terminal-state artefact and nothing downstream should go ahead. |
| Q2 | Is it present **in situ**, in the expected cell types and at protein level, independent of vessel-isolation yield? | snRNA-seq of isolated vessels cannot separate cell-intrinsic change from capture bias. |
| Q3 | Is the HIF response of AD-relevant NVU cells **intrinsically blunted**, and does that change vascular function? | This turns a correlation into a testable mechanism and decides whether any gene is a target. |

The plan is organised as gated stages. Each stage has a pre-specified decision rule, so a negative result
stops further spending.

---

## Stage 0: metadata and in-silico replication (no bench work; ~4 weeks)

1. **Donor covariates.**

   *Done.* Yang *et al.* Supplementary Table 1 (`data/41586_2021_4369_MOESM2_ESM.xlsx`) was parsed by
   `scripts/10_donor_covariates.py`. Its donor IDs were verified against our samples: sex agrees for 17/17
   donors and the cortex subset for 17/17. Results (Fig. S2):
   - AD and control donors are balanced for age, post-mortem interval, APOE ε4 dose and recorded acute hypoxic
     neuropathology.
   - The hypoxia-module score is not associated with any of these.
   - The AD reduction persists under every adjustment (pericytes −1.48, P = 0.0007 with age + PMI) and when
     the two hypoxic-pathology donors are excluded.

   *Still missing.* Agonal-factor score, cause and manner of death, terminal ventilation and illness
   duration, brain pH and RIN. The table does not report these; request them from the Stanford/VA/NIA Aging
   Clinical Research Center (ACRC), which Yang *et al.* name as the tissue source.
   - Treat agonal score, pH and RIN as confounders, using the same two analyses as in `10_donor_covariates.py`.
   - Use Braak/CERAD/APOE for a dose–response test, not as adjustment, because they are part of the disease.
     Braak is currently reported for only 5 of 9 AD donors.

2. **Replicate in independent snRNA-seq cohorts.** Prefer cohorts with vascular coverage and agonal/pH data:
   SEA-AD (MTG; detailed donor metadata); the ROSMAP multiregion cerebrovascular snRNA-seq of Sun *et al.*,
   *Nat. Neurosci.* 26:970–982 (2023), with 22,514 vascular nuclei from six regions including hippocampus;
   and the vascular dissection of Garcia *et al.*, *Nature* 603:893–899 (2022). Apply the same pseudobulk
   model and score the same 8-gene module.

3. **Decision gate G0.**

   | Result | Action |
   |---|---|
   | Module effect survives adjustment for the remaining agonal/pH variables (or these cannot be obtained) **and** replicates in ≥ 1 independent cohort | Proceed to Stage 1 with full scope. Matching on agonal score and pH in Stage 1 then becomes the decisive test. |
   | Effect disappears after adjustment | Report it as a terminal-state confound. Re-scope Stage 1 to the pericyte-state and ADAMTS9 questions only, using an agonal-matched cohort. |

---

## Stage 1: in-situ validation in human hippocampus (~4–6 months)

### Cohort
- **Groups.** Hippocampus (CA1, subiculum, entorhinal border) from neuropathologically confirmed AD
  (Braak V–VI, CERAD moderate–frequent) and from controls without cognitive impairment (Braak 0–II).
- **Matching.** Groups are frequency-matched for age (± 5 y), sex, PMI (< 12 h preferred), brain pH (± 0.2)
  and agonal score.
- **Calibration group.** Include donors with documented prolonged agonal hypoxia (for example, extended
  ventilation or hypoxic-ischaemic injury) regardless of diagnosis. This measures how strongly agonal hypoxia
  alone moves the module.
- **Sample size.** Per-donor effect sizes in this dataset are Cohen's d ≈ 1.1–2.0 (median by gene:
  *ADAMTS9* 1.4, *VEGFA* 1.8, *HILPDA* 2.0, *SLC2A3* 1.7, *PDK1* 1.1;
  [`tables/validation_power_estimates.csv`](../tables/validation_power_estimates.csv)). These estimates come
  from the discovery data and are inflated by selection (the "winner's curse"), so assume d ≈ 0.8–1.0. That
  requires **n = 17–26 donors per group** for 80 % power at α = 0.05 (two-sided). Plan **n = 20 AD + 20
  control + 8 calibration**.
- **Blinding and pre-registration.** Regions of interest and quantification are done blind to diagnosis. The
  primary endpoint and analysis are pre-registered.

### 1A. Multiplex RNA in situ hybridisation (RNAscope HiPlex or equivalent smFISH)
- **Probe panel 1 (module × cell type):**
  - *ADAMTS9, VEGFA, HILPDA, SLC2A3*
  - cell identity: *PDGFRB* (pericyte), *CLDN5* (endothelium), *AQP4* (astrocyte), *P2RY12* (microglia)
- **Probe panel 2 (pericyte state):** *ADAMTS9, ELL2, PCBP3* (state S2); *ATP1A2, SLC20A2* (state S1);
  *PDGFRB, KCNJ8* (pericyte identity).
- **Controls.** *PPIB/POLR2A* (positive control, also used as an RNA-integrity covariate) and *DapB*
  (negative control) on every slide.
- **Readout.** Puncta per cell within each identity-defined cell type, from ≥ 300 cells per cell type per
  donor. Also the fraction of *PDGFRB⁺* cells that are S2-high, using a threshold defined on control tissue
  before unblinding.
- **Primary endpoint.** Donor-mean module score (z-averaged puncta for *ADAMTS9, VEGFA, HILPDA, SLC2A3*) in
  *PDGFRB⁺* pericytes, AD vs control. Analysis is a linear mixed model with donor as a random effect and
  agonal score, pH, PMI, sex and positive-control puncta as covariates.
- **Why ISH first.** It counts transcripts in cells identified in place, so it is unaffected by vessel-isolation
  yield and ambient RNA. Those are the two artefacts that dominated the sequencing data.

### 1B. Protein-level confirmation (multiplex immunofluorescence)
- **Targets:** ADAMTS9; VEGF-A; GLUT3 (SLC2A3); HIF-1α; and the ADAMTS cleavage product of versican.
  - HIF-1α degrades within minutes after death, so it is interpretable only with short PMI. Score nuclear
    HIF-1α in pericytes as exploratory.
  - The versican fragment is detected with a neoepitope antibody to the C-terminal DPEAAE sequence. It reads
    out ADAMTS9 *activity* in the perivascular matrix.
- **Vessel masks:** collagen IV or *Ulex europaeus* lectin (vessels), CD13/PDGFRβ (pericytes), AQP4
  (astrocyte endfeet).
- **Readouts:**
  - perivascular intensity of each target, normalised to vessel area
  - pericyte coverage (CD13⁺ area / lectin⁺ area)
  - BBB leakage (extravascular fibrinogen and IgG)
- **Antibody qualification.** Use only clones with independent validation: knockout/knockdown evidence,
  Human Protein Atlas *enhanced* validation, or two antibodies to non-overlapping epitopes that agree. Confirm
  each antibody on siRNA-treated human pericytes (Stage 3) before tissue use. Catalogue numbers are not
  pre-specified here.

### 1C. Isolated-microvessel bulk assays (fresh-frozen tissue from the same donors)
- **Isolation.** Microvessels are isolated by density-gradient/filtration (the VINE-seq first step). Record
  yield as vessel protein per mg tissue.
- **Readouts:**
  - RT-qPCR of the 8-gene module, normalised to *CLDN5* and *PDGFRB*, which corrects for the
    endothelial/mural share
  - western blot for ADAMTS9, VEGF-A and GLUT3
- **Purpose.** Tests whether yield differs by diagnosis in an independent isolation, and whether the module
  difference holds *within* the vessel fraction.

### Decision gate G1

| Result | Action |
|---|---|
| Primary ISH endpoint is significant after covariate adjustment, in the same direction as the snRNA-seq (lower in AD), **and** the calibration group shows the module rises with agonal hypoxia | The effect is AD-associated beyond terminal state. Proceed to Stage 2. |
| The difference is explained by agonal or pH covariates | Stop the mechanism work. Report the module as a terminal-state confound in post-mortem vascular transcriptomics, which is itself a useful methodological result. |

---

## Stage 2: does an AD milieu blunt the NVU hypoxic response? (~4 months, in parallel with late Stage 1)

**Cells.** Primary human brain vascular pericytes, human brain microvascular endothelial cells and human
astrocytes, each from ≥ 3 donors/lots. iPSC-derived pericytes and BMEC-like cells from isogenic *APOE3/E4* or
*APP*-mutant lines are a genetic arm.

**Design (2 × 4 factorial per cell type):**
- **Oxygen:** normoxia (21 % O₂) vs hypoxia (1 % O₂, 6 h and 24 h). A pharmacological arm uses the
  prolyl-hydroxylase inhibitor DMOG to separate HIF stabilisation from O₂ effects.
- **AD-relevant exposure:** vehicle; Aβ42 oligomers (0.5–1 µM, characterised by SEC/AFM); TNF-α + IL-1β
  (10 ng/mL, inflammatory milieu); conditioned medium from AD-patient iPSC astrocytes or microglia.

**Readouts:**
- **RNA:** RT-qPCR of the module (*ADAMTS9, VEGFA, HILPDA, SLC2A3, PDK1, BNIP3L, ELL2*) and canonical HIF
  targets (*CA9, EGLN3*). Bulk RNA-seq on the 24 h arm.
- **HIF activity:** HIF-1α/HIF-2α protein (western, nuclear fraction) and an HRE-luciferase reporter.
- **Secretion:** VEGF-A (ELISA) and ADAMTS9 in the medium.
- **Metabolism:** glycolytic response (Seahorse ECAR/OCR) and glucose uptake (2-NBDG).

**Prediction if the finding is real.** AD-relevant exposures reduce hypoxic induction of the module by
≥ 30 %, measured as the interaction term in a two-way ANOVA with donor/lot as a blocking factor.

---

## Stage 3: function of the lead genes (~4–6 months; only if G1 passes)

1. **ADAMTS9 in pericytes** (siRNA and CRISPRi; rescue with a catalytically dead E→A mutant vs wild type):
   - versican and aggrecan cleavage (DPEAAE neoepitope)
   - deposition of basement-membrane proteins (collagen IV, laminin β1; both lower in AD here)
   - pericyte migration and contractility (collagen-gel contraction)
   - in pericyte–endothelial co-culture: TEER, 4/70 kDa dextran permeability and tube stability
2. **Restoring the HIF response.** In AD-milieu-treated pericytes and endothelial cells from Stage 2, test
   whether pharmacological HIF stabilisation (DMOG or roxadustat) or *VEGFA* re-expression rescues barrier
   and coverage readouts. This tests whether the blunted response matters functionally.
3. **Pericyte state S2.** Sort *PDGFRB⁺* pericytes by ADAMTS9 reporter level (or index-sort with ISH) and
   profile them. The question is whether S2 is a hypoxia-inducible *state* that resting pericytes enter,
   rather than a fixed subtype. That distinction decides whether its loss in AD reflects reduced
   inducibility or cell loss.

---

## Stage 4: in vivo confirmation in mouse (~6 months; optional, after Stage 3)

- **Orthology.** The lead genes are well conserved: ADAMTS9 88.9 %, VEGFA 83.8 %, HIF1A 90.7 %, SLC2A3 82.9 %
  human–mouse protein identity, computed in `tables/target_annotation_top.csv`. *ADGRE2* has no mouse
  orthologue and is excluded from mouse work.
- **Models.** An amyloid model (*APP^NL-G-F* knock-in, chosen over transgenic overexpression to avoid
  artefacts) and wild-type littermates. Ages are 6 and 12 months; n = 10 per group per age, both sexes.
- **Challenge.** Acute hypoxia (8 % O₂, 6 h) or bilateral carotid stenosis (chronic hypoperfusion). Use
  pimonidazole to mark hypoxic tissue.
- **Readouts:**
  - hippocampal microvessel qPCR/RNA-seq and RNAscope of the module in *Pdgfrb⁺* pericytes
  - pericyte coverage, BBB leakage (Evans blue, fibrinogen) and cerebral blood flow (laser speckle)
- **Prediction.** *APP^NL-G-F* mice show a smaller induction of the module after hypoxia than wild-type mice,
  and this coincides with reduced pericyte coverage or increased leakage.
- **Optional causal test.** Pericyte-specific *Adamts9* deletion (*Pdgfrb-CreER*) crossed onto
  *APP^NL-G-F*.

---

## Statistics and reporting (all stages)

- **Replication unit.** The donor (human), the independent cell lot (in vitro) or the animal (in vivo).
  Technical replicates are averaged first.
- **Pre-registration.** Primary endpoints, covariates and exclusion criteria are registered before
  unblinding. Multiplicity across secondary endpoints is controlled by BH-FDR.
- **Analysis.** Mixed models with donor/lot random effects; effect sizes are reported with 95 % CIs.
  Negative results are reported, as are the agonal-calibration results.
- **Data release.** Imaging data and per-cell quantification are deposited, together with the analysis code,
  in this repository.

## Timeline and go/no-go summary

| Months | Stage | Go criterion |
|---|---|---|
| 0–1 | 0: metadata + independent cohorts | Effect survives agonal/pH adjustment and replicates in ≥ 1 cohort |
| 1–7 | 1: ISH + IF + microvessel assays (n = 20 + 20 + 8) | Primary ISH endpoint significant, lower in AD, after covariates |
| 4–8 | 2: in vitro AD-milieu × hypoxia | Significant AD-milieu × hypoxia interaction on module induction |
| 7–13 | 3: ADAMTS9 / HIF functional tests | Knockdown or rescue changes barrier, ECM or coverage readouts |
| 10–16 | 4: mouse model (optional) | Blunted hypoxic induction plus a vascular phenotype in *APP^NL-G-F* |

## Main risks and mitigations

| Risk | Mitigation |
|---|---|
| The module is an agonal artefact | Stage 0 and the Stage 1 calibration group are designed to detect this early and cheaply. |
| HIF-1α protein does not survive post-mortem delay | Use transcript (ISH) and stable downstream proteins (VEGF-A, GLUT3, versican neoepitope) as primary readouts. |
| Antibody non-specificity (ADAMTS family members are similar) | Knockdown-validated antibodies; ISH as an orthogonal readout. |
| In vitro pericytes lose their in vivo identity | Use low passage (≤ 5); confirm *PDGFRB/KCNJ8/ABCC9* expression; add an iPSC-pericyte arm. |
| Effect specific to the hippocampus | Stage 1 includes an entorhinal/subiculum ROI. The cortex is analysed only as a secondary endpoint. |
