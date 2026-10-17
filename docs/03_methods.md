# Methods

This section describes the re-analysis pipeline in `scripts/` (steps 01–10, run in order by `scripts/run_all.py`).
Every parameter given here is the value used in the code. Shared settings (QC thresholds, marker sets, excluded
gene lists, random seed) are in `scripts/config.py`. All stochastic steps use seed 0.

## Data and study design

**Source data.** Filtered single-nucleus count matrices (10x Genomics, CellRanger output) for 25 samples
were downloaded from GEO series GSE163577. These are the VINE-seq vascular-enriched nuclei of Yang *et al.*,
*Nature* 603:885 (2022).
- **Hippocampus:** 17 donors, 8 cognitively normal controls and 9 with Alzheimer's disease (AD).
- **Superior frontal cortex:** 8 samples from 8 of the same donors (4 control, 4 AD).

Raw reads were not re-processed, and no ambient-RNA correction (SoupX, CellBender) was possible, because
GEO provides only filtered matrices.

**Sample and donor identifiers.** Sample metadata were parsed from the GEO file names (`config.sample_table()`).
- Cortex tokens use the letter "O" in place of zero (e.g. `O4_11AD`); after normalising this, each cortex sample
  maps to exactly one hippocampal sample from the same donor.
- Condition labels were taken from the GEO sample titles and checked against the file-name suffix wherever one
  exists.
- The donor pairing was verified independently by inferred sex and 17q21.31 genotype (see QC below).

**Unit of replication.** The unit of replication is the sample (one donor × one region); nuclei are never
treated as replicates.
- The hippocampus (17 independent donors) is the discovery set.
- The cortex samples come from 8 of the same donors. They are analysed as a second tissue from the same
  individuals and are reported as cross-region concordance, not as independent replication.

**Donor clinical data.** Age, sex, post-mortem interval (PMI), Braak stage, CERAD score, APOE genotype, brain
weight and neuropathology notes came from Yang *et al.* Supplementary Table 1
(`data/41586_2021_4369_MOESM2_ESM.xlsx`), as described in "Donor covariates" below.

## Quality control and doublet removal (`01_qc_doublets.py`)

**Per-nucleus metrics.** Each sample was processed separately, computing total UMI counts, number of detected
genes and the percentage of mitochondrial (`MT-`) reads.

**Low-quality nuclei.** A nucleus was classed as low quality if any of the following held:
1. log(1 + total counts) lay more than 5 median absolute deviations (MADs) from the sample median, in either
   direction;
2. log(1 + detected genes) lay more than 5 MADs from the sample median, in either direction;
3. fewer than 200 genes were detected;
4. the mitochondrial percentage exceeded a per-sample upper limit. This limit was 3 MADs above the sample
   median, bounded to the range 1–5 %.

**Doublets.** Scrublet (`scanpy.pp.scrublet`, seed 0) was run per sample on the nuclei that passed QC.
- Nuclei called doublets by Scrublet's automatic threshold were removed.
- Where the automatic threshold could not be determined, a fixed doublet-score cut-off of 0.25 was applied, so
  that every sample was filtered.

**Gene filter.** Genes detected in fewer than 10 retained nuclei across the whole dataset were removed.

**Result.** Of 319,449 input nuclei, 18,863 were low quality and 30,630 were doublets; 269,956 nuclei were
retained. No down-sampling was applied. Per-sample thresholds and counts are in `tables/qc_sample_summary.csv`,
and every nucleus with its QC status is in `tables/qc_cell_metrics.csv.gz`.

**Donor sex and 17q21.31 genotype.** Both were inferred from sample-level pseudobulk counts per million (CPM).
- **Sex:** female if XIST CPM exceeded the summed CPM of six Y-chromosome genes (*UTY, USP9Y, DDX3Y, KDM5D,
  NLGN4Y, TTTY14*), otherwise male.
- **17q21.31 genotype:** *ARL17B* CPM is bimodal across samples (about 10–60 vs 150–400). Samples above
  100 CPM were classed as carriers of the copy-number-variable 17q21.31 haplotype.
- **Concordance:** sex and genotype calls were compared between the two tissues of each paired donor
  (`tables/donor_pairing_concordance.csv`).

## Normalisation, integration and clustering (`02_integrate_annotate.py`)

**Normalisation.** Counts were normalised to 10,000 per nucleus and log(1 + x)-transformed; raw counts were
kept in a separate layer.

**Highly variable genes (HVGs).** 3,000 HVGs were selected with the Seurat v3 method on raw counts, with
`batch_key = sample`. Genes excluded before selection:
- sex-chromosome genes (XIST, TSIX and the Y-linked genes in `config.SEX_GENES_PATTERN`);
- the 17q21.31 locus genes (`config.LOCUS_17Q21`);
- mitochondrial (`MT-`) and ribosomal (`RPL`/`RPS`) genes.

These genes vary with donor sex, genotype or technical quality rather than with cell identity.

**Integration.** Scaled HVG expression (clipped at 10) was reduced to 50 principal components. Harmony
(harmonypy, at most 20 iterations, batch variable = sample) was then applied, so hippocampus and cortex are
integrated in one shared embedding.

**Graph, clustering and UMAP.**
- k-nearest-neighbour graph: k = 30, built on the first 30 Harmony components.
- Leiden clustering: igraph implementation, resolution 1.0, 2 iterations.
- UMAP: min_dist 0.3.

**Annotation.** Each Leiden cluster was labelled as a whole, not nucleus by nucleus.
1. A module score (`scanpy.tl.score_genes`) was computed for each canonical marker set in `config.MARKERS`
   (13 sets, e.g. endothelial *CLDN5, FLT1, VWF, PECAM1, ABCB1*; pericyte *PDGFRB, RGS5, ABCC9, KCNJ8*).
   The scores were averaged per cluster and z-scored across clusters.
2. Each cluster took the label with the highest z-score.
   - If the best and second-best labels belonged to the same lineage family (mural = pericyte, smooth muscle,
     fibroblast; myeloid = microglia, macrophage; neuron = excitatory, inhibitory), the label was chosen on
     the raw score within that family.
   - Excitatory and inhibitory neurons could not be separated at cluster level and were labelled jointly as
     "Neuron".
3. A cluster was labelled **Mixed** in either of two cases:
   - the relative margin between its best and second-best z-scores was below 25 %;
   - **heterotypic doublets:** the key marker of another lineage family (*CLDN5, PDGFRB, AQP4, C3, MBP, PDGFRA,
     GAD1, CD3E, FOXJ1*) was detected in at least 30 % of its nuclei *and* in at least three times that marker's
     median detection rate in clusters of other families. The threefold background rule tolerates the
     ambient *MBP* signal present in all clusters.

**Result.** 7 of 26 clusters (12,931 nuclei) were labelled Mixed and excluded from all downstream analyses,
leaving 257,025 nuclei in 11 classes. Per-cluster scores, margins, detection rates and the reason for each label
are in `tables/cluster_annotation.csv`.

**Vascular yield.** VINE-seq enriches vascular nuclei, and the enrichment is less efficient in AD samples. For
each sample, the fraction of annotated nuclei that were endothelial, pericyte, smooth muscle or fibroblast was
computed. This fraction was bounded to [0.001, 0.999] and logit-transformed (`vascular_yield`). Because yield
changes both composition and the ambient-RNA background, it is a covariate in all differential analyses.

## Composition (`03_fig1_global_landscape.py`)

Per-sample cell-type proportions and vascular fractions were compared between AD and control within each region
by two-sided Mann–Whitney U tests, with Benjamini–Hochberg (BH) FDR correction across cell types and regions.
Because composition depends on enrichment efficiency, these comparisons are descriptive. Outputs:
`tables/Fig1_composition_per_sample.csv`, `Fig1_composition_tests.csv` and `Fig1_vascular_yield_tests.csv`.

## Pseudobulk differential expression (`04_pseudobulk_de.py`)

**Cell types.** Five neurovascular-unit (NVU) cell types were analysed in each region: microglia,
endothelial cells, pericytes, astrocytes and oligodendrocytes.

**Pseudobulk construction.**
- Raw counts were summed per sample within each cell type.
- Samples with fewer than 30 nuclei of that cell type were dropped.
- A gene was kept if it had at least 10 counts and at least 5 CPM in at least *n* samples, where *n* is the size
  of the smaller condition group.

**Model.** DESeq2 as implemented in pydeseq2 was fitted with design `~ sex + vascular_yield + condition`.
- Vascular yield was z-scored within each test.
- Outliers flagged by Cook's distance were refitted.
- The AD-vs-control Wald test was used, with BH-adjusted *P* values.
- Shrunken log2 fold changes (pydeseq2 `lfc_shrink`) are reported alongside the maximum-likelihood estimates.
  Thresholds use the maximum-likelihood log2 fold change.

**Artefact flags.** Two flags are attached to each gene in each cell type.
- **`flag_contamination` (ambient RNA or doublet carry-over).** Both conditions must hold:
  1. the gene's mean CPM in some other cell type of the same region is at least 5× its CPM in the target cell
     type (+0.5 pseudocount);
  2. across samples, its CPM in the target pseudobulk correlates with that sample's fraction of the source cell
     type (Spearman ρ ≥ 0.5).
- **`flag_genotype_sex`.** The gene lies in the 17q21.31 locus or is a sex-chromosome gene. Expression of these
  genes is set by genotype or sex.

**Robust DEG.** A gene is called a robust differentially expressed gene if BH-adjusted *P* < 0.05, |log2 fold
change| ≥ 0.5 and neither flag is set. Every downstream analysis uses only robust DEGs.

**Sensitivity analysis.** Each model was refitted without the yield term (`~ sex + condition`). The number of
genes at FDR < 0.05 and the correlation of log2 fold changes between the two models are reported in
`tables/pseudobulk_DE_sensitivity_no_yield.csv`.

## Cell-type deep dives (`05_fig2_deep_dives.py`)

**Sub-clustering.** Each NVU cell type was sub-clustered with both regions together.
- Genes detected in fewer than 20 nuclei, and the sex, 17q21.31, mitochondrial and ribosomal genes, were
  excluded.
- 2,000 Seurat v3 HVGs were selected per sample. Where the LOESS fit failed for samples with few nuclei,
  selection was done without the batch key.
- The data were reduced to 30 principal components, Harmony was applied (batch variable = sample), and a k = 20
  neighbour graph was built on 20 components.
- Leiden clustering was run at resolution 0.4.
- Sub-clusters with fewer than 1 % of nuclei were merged into the nearest sub-cluster by Harmony-space centroid
  distance. Sub-clusters were then renamed S1…Sn by size.

**Marker genes.** Markers per sub-cluster came from Wilcoxon rank-sum tests: log2 fold change > 0.5, adjusted
*P* < 0.01, detection in more than 25 % of sub-cluster nuclei, with sex, 17q21.31, mitochondrial and ribosomal
genes excluded. Markers describe the sub-clusters only and are not used for AD-vs-control inference.

**Sub-cluster proportions.** Per-sample sub-cluster proportions were compared between AD and control within each
region by two-sided Mann–Whitney U tests, with BH FDR correction across all sub-clusters and both regions of
that cell type (`tables/Fig2_<cell type>_subcluster_tests.csv`).

**Donor-level heatmaps.** These show the 15 robust DEGs with the smallest *P* values (union of both regions), as
z-scored log2(CPM + 1) per sample.

**Gene-set enrichment.** Pre-ranked GSEA (gseapy `prerank`) was run on the DESeq2 Wald statistic, with
artefact-flagged genes removed before ranking.
- Gene-set sizes: 15–500 genes.
- Permutations: 1,000 gene-set permutations.
- Libraries: three Enrichr libraries, MSigDB_Hallmark_2020, GO_Biological_Process_2023 and Reactome_2022,
  downloaded from maayanlab.cloud and stored in `data/genesets/`.

Each library was tested separately. The figures show up to 12 terms with FDR < 0.25 (the 6 smallest FDRs per
region, union); full results are in `tables/GSEA_<Region>_<CellType>.csv`.

## Cross-cell-type convergence and cross-region concordance (`06_fig3_convergence.py`)

**Convergent genes.** These are genes robust in at least two of the five hippocampal NVU cell types
(`tables/convergent_genes.csv`).

**Similarity of AD effects between cell types.** This was measured as the Spearman correlation of hippocampal
log2 fold changes. The correlation used unflagged genes tested in all five cell types.

**Cross-region concordance.** For each cell type, the hippocampal robust DEGs that were also tested in cortex
were compared between the two regions.
- **Sign agreement:** the share of these genes changing in the same direction in both regions, tested against
  50 % with an exact two-sided binomial test.
- **Effect agreement:** the Spearman correlation of hippocampal and cortical log2 fold changes.

Because the cortex donors are a subset of the hippocampal donors, concordance indicates within-individual
consistency across tissues. It is not independent replication.

**Immediate-early genes.** The immediate-early genes (*FOS, FOSB, JUN, JUNB, EGR1, DUSP1, GADD45B, ATF3, IER2,
NR4A1*) were shown across all ten region × cell-type contrasts. This tests a claim made in the superseded
analysis that these genes were "universally convergent".

## Ligand–receptor communication (`07_fig4_cell_communication.py`)

**Resource.** Ligand–receptor pairs came from the LIANA "consensus" resource (liana 0.1.9).
- Each multi-subunit complex was represented by its least-expressed subunit (minimum mean expression and minimum
  detection fraction).
- APP–SORL1 is not in the consensus resource. It was added so that the earlier claim about this pair could be
  tested directly.

**Interaction score.** For every sample and every sender–receiver pair of NVU cell types (each with at least 30
nuclei in that sample), the score was the geometric mean of the ligand's mean log-normalised expression in the
sender and the receptor's mean in the receiver. This follows LIANA's `lr_means`. An interaction was called
expressed when the ligand and the receptor were each detected in at least 10 % of nuclei of their cell type.

**Differential communication.**
- Interactions expressed in at least 50 % of samples of either condition were tested within each region.
- For each interaction, an ordinary least squares (OLS) regression was fitted to log2(score + 0.01) on condition,
  sex and z-scored vascular yield, with one observation per sample.
- BH FDR was computed within each region.
- Interactions whose ligand or receptor was artefact-flagged in the relevant cell type are reported, but excluded
  from the summary panels.
- The number of expressed interactions per sample was compared between conditions by Mann–Whitney U test.

## Target prioritisation (`08_fig5_target_prioritization.py`)

**Candidates.** Candidates were genes robust in at least one hippocampal NVU cell type.

**Ranking.** Candidates were ranked lexicographically, with no hand-set weights, on three criteria in order:
1. **Breadth:** the number of hippocampal NVU cell types in which the gene is robust.
2. **Cortical support:** the number of those cell types in which the cortex effect has the same sign and
   BH-adjusted *P* < 0.1.
3. **Evidence:** −log10 of the smallest hippocampal adjusted *P*.

**Annotation of the top 20 candidates.** All annotations were retrieved programmatically and cached in
`data/annotation/`.
- **Protein data:** accession, protein name, subcellular location, transmembrane and signal-peptide features from
  reviewed human UniProtKB/Swiss-Prot entries (UniProt REST API).
- **Mouse orthologue:** found by gene symbol among reviewed mouse entries.
- **Human–mouse sequence identity:** from a global Needleman–Wunsch alignment (Biopython `PairwiseAligner`,
  BLOSUM62, gap open −10, gap extend −0.5), computed as identical aligned positions divided by human sequence
  length.
- **Accessibility class** (secreted/extracellular, cell surface/membrane, intracellular, unknown): assigned from
  the UniProt location and features.

**Functional network.** A STRING v12 network (human, combined score ≥ 0.7) was retrieved for the top candidates
together with a fixed panel of AD GWAS and reference genes (`tables/STRING_edges.csv`).
- Gene symbols were mapped to STRING identifiers by exact preferred-name match only.
- Symbols without an exact match (VEGFA) are reported as unmapped rather than replaced by a fuzzy match.

## Donor covariates and sensitivity analyses (`10_donor_covariates.py`)

**Mapping the clinical table.** In Supplementary Table 1 of Yang *et al.*, the identifiers "NCI *k*" and
"AD *k*" correspond to the index *k* in the GEO sample titles ("Control *k* Hpc", "AD *k* Hpc"). The mapping was
accepted only after two independent checks, both of which passed for 17/17 donors:
1. reported sex agreed with the sex inferred from XIST and Y-chromosome expression;
2. the donors listed with a cortex sample matched the donors with a GEO cortex sample.

**Variables extracted.**
- age, PMI, Braak stage, CERAD score, APOE genotype (converted to ε4 allele dose) and brain weight;
- two binary flags parsed from the neuropathology notes by regular expression:
  - acute hypoxic neuropathology ("acute hypoxic", "hypoxic-ischemic", "hypoxic damage");
  - recorded infarcts.

The table reports no agonal-factor score, brain pH or RIN. Acute hypoxic neuropathology is therefore used as the
only available proxy for terminal hypoxia.

**Covariate balance.** Continuous covariates were compared between AD and control donors by Mann–Whitney U test,
and binary covariates by Fisher's exact test (`tables/covariate_balance.csv`).

**Hypoxia-module score.** The module contains eight genes: *ADAMTS9, VEGFA, HILPDA, SLC2A3, PDK1, BNIP3L, ELL2,
HIF1A*. These are the canonical hypoxia-inducible (HIF-target) genes among the hippocampal convergent DEGs.
- Per sample and cell type, each module gene's log2(CPM + 1) was z-scored across samples, and the eight z-scores
  were averaged.
- Samples with fewer than 30 nuclei were excluded.
- Because the genes were chosen from the same data, the module's AD effect size is inflated by selection. The
  score is used to test whether covariates explain the effect, not to test the effect itself.

**Score models.** For each hippocampal NVU cell type, the score was regressed by OLS on condition + sex +
vascular yield (the primary model). The primary model was then extended by adding each covariate in turn:
- age;
- PMI;
- APOE ε4 dose;
- acute hypoxic pathology;
- age and PMI together.

The primary model was also refitted after excluding the two donors with acute hypoxic pathology. Spearman
correlations of the score with age, PMI and APOE ε4 dose were computed separately.

**DE sensitivity.** The hippocampal DESeq2 models were refitted in two ways:
1. with design `~ sex + vascular_yield + age + pmi_h + condition`;
2. with the primary design after excluding the two donors with acute hypoxic pathology (n = 15).

For each refit, the following were recorded: the share of the primary robust DEGs still at FDR < 0.05, the
share keeping the same sign, and the correlation of log2 fold changes with the primary fit across all genes
(`tables/pseudobulk_DE_sensitivity_clinical.csv`).

**Effect sizes for validation planning.** For each hypoxia-module gene in hippocampal pericytes, endothelial cells,
astrocytes and microglia, per-donor Cohen's *d* between control and AD was computed on log2(CPM + 1). The number of donors per group needed for 80 % power at two-sided
α = 0.05 was derived with a two-sample *t*-test power calculation (statsmodels `TTestIndPower`), reported in
`tables/validation_power_estimates.csv`. These *d* values come from genes selected for significance in the same
data, so they overestimate the true effect; the validation plan uses smaller, assumed effect sizes.

## Statistics

- **Tests:** all tests are two-sided. Multiple testing is controlled with the Benjamini–Hochberg procedure within
  the family stated for each analysis.
- **Replication unit:** the sample (donor × region) in every between-condition test.
- **No pooled-nuclei tests:** no test treats nuclei as independent observations, except the within-cell-type
  Wilcoxon tests used to describe sub-cluster markers.
- **Groups:** the condition contrast is always AD versus control, with control as the reference level.
- **Displayed summaries:** these are medians with per-sample points.

## Figures

**Generation.** Figures were produced with matplotlib using shared settings (`config.set_style`):
- colour-blind-safe Okabe–Ito-based palettes;
- Arial/Helvetica text at 5–9 pt, with TrueType fonts embedded in the PDFs;
- vector PDF plus 300-dpi PNG output.

**Single panels.** Each lettered panel is also exported on its own to `single_figures/`, without its panel
letter or title, by `config.export_single_panels`. Panel legends are in `figures/README.md`.

## Software

| Package | Version |
|---|---|
| Python | 3.12.15 |
| scanpy | 1.12.4 |
| anndata | 0.13.4 |
| harmonypy | 0.2.0 |
| python-igraph / leidenalg | 1.0.0 / 0.12.0 |
| pydeseq2 | 0.5.4 |
| gseapy | 1.3.1 |
| liana | 0.1.9 |
| statsmodels | 0.15.0 |
| SciPy | 1.18.1 |
| Biopython | 1.88 |
| NumPy / pandas | 2.5.4 / 3.0.6 |
| matplotlib | 3.11.2 |

The conda environment is specified in `environment.yml`. The pipeline is run with `python scripts/run_all.py`.
External resources were cached on first use: Enrichr gene-set libraries, UniProt entries and STRING responses.

## Differences from the original analysis of Yang *et al.* (2022)

- **Unit of inference.** DE is done on donor-level pseudobulk rather than on pooled nuclei.
- **Yield covariate.** Differential models include vascular-enrichment yield as a covariate.
- **Artefact flags.** Genes are flagged for ambient or doublet carry-over and for 17q21.31 or sex dependence,
  and flagged genes are excluded from the robust set.
- **Regions.** Hippocampus and cortex are integrated in one embedding, but cortex is treated as a second tissue
  from the same donors rather than as replication.
- **Ambient RNA.** Filtered GEO matrices were used, so ambient RNA could only be flagged after the fact, not
  removed.
