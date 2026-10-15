# Figure legends

Every figure is produced by the script named in brackets and is saved as a vector PDF (fonts embedded as
Type 42) and as a 300-dpi PNG. Hpc = hippocampus (8 control, 9 AD donors); Ctx = superior frontal cortex
(4 control, 4 AD, all from donors who also contributed hippocampus). Points are samples (donor × region),
horizontal bars are medians, and all tests are donor-level.

## Figure 1 | Single-nucleus landscape of both regions (`03_fig1_global_landscape.py`)
**a**, Integrated UMAP of all 257,025 retained nuclei from 25 samples (Harmony on sample), coloured by cell
class; n per class in brackets. **b**, The same embedding coloured by region. **c**, The same embedding coloured
by condition. **d**, Canonical markers per class (dot size = fraction of nuclei expressing; colour = min–max
scaled mean). **e**, Cell-class composition of each sample, grouped by region and condition and ordered by
vascular fraction. **f**, Vascular nuclei (endothelial + pericyte + smooth muscle + fibroblast) per sample;
two-sided Mann–Whitney P. Composition reflects VINE-seq vessel-enrichment yield, not tissue abundance
(`tables/Fig1_composition_tests.csv`).

## Figure 2a–e | NVU cell-type deep dives (`05_fig2_deep_dives.py`)
Every figure uses the same four-row layout, read left to right: a, b / c, d / e, f / g, h.
One figure each for microglia (2a), endothelial cells (2b), pericytes (2c), astrocytes (2d) and
oligodendrocytes (2e). **a**, Sub-clusters (Harmony-integrated, Leiden 0.4; clusters < 1 % merged).
**b**, Region and condition on the same embedding. **c**, Top Wilcoxon markers per sub-cluster (descriptive only;
not used for AD-vs-control inference). **d**, Sub-cluster proportions per sample (* Mann–Whitney BH-FDR < 0.05).
**e, f**, Pseudobulk DESeq2 volcanoes, `~ sex + vascular_yield + condition`. Filled colour = robust
(FDR < 0.05, |log₂FC| ≥ 0.5, not flagged); open grey = ambient / 17q21.31 / sex-flagged; dotted line = FDR 0.05.
**g**, Per-sample z-scored log₂CPM of the top robust DEGs (both regions). **h**, Pre-ranked GSEA on the
DESeq2 Wald statistic with flagged genes removed (MSigDB Hallmark 2020, GO-BP 2023, Reactome 2022; 1,000
permutations). Up to 12 terms with FDR < 0.25 are shown.

## Figure 3 | Convergence and cross-region concordance (`06_fig3_convergence.py`)
**a**, DEGs per region × cell type, split into robust vs artefact-flagged. **b**, UpSet plot of hippocampal robust
DEGs across the five cell types. **c**, Spearman correlation of hippocampal AD log₂FC between cell types
(unflagged genes tested in all five). **d**, The 47 genes robust in ≥ 2 hippocampal cell types, with log₂FC in all
10 contrasts (• FDR < 0.05). **e**, Hippocampal robust DEGs compared with the cortex effect in the same
donors; per-cell-type sign agreement (binomial test in `tables/cross_region_concordance.csv`).
**f**, Immediate-early genes that the legacy analysis described as universally convergent.

## Figure 4 | Ligand–receptor communication at donor level (`07_fig4_cell_communication.py`)
**a**, Expressed L–R interactions per sample (ligand and receptor each detected in ≥ 10 % of nuclei).
**b, c**, Median AD log₂FC of interaction scores per sender → receiver pair. Cell text gives the number with
nominal P < 0.01 over the number tested. No pair reaches FDR < 0.05. **d**, APP→SORL1 into microglia per sample
(the legacy headline claim). **e**, The 20 hippocampal interactions with the smallest nominal P, with 95 % CI,
and the cortex estimate for the same interactions.

## Figure 5 | Computed target prioritisation (`08_fig5_target_prioritization.py`)
**a**, The 20 top-ranked candidates. Ranking is lexicographic: hippocampal breadth, then cortex support, then
evidence. log₂FC is shown in all 10 contrasts (• FDR < 0.05; n.t. = not tested). **b**, Annotation retrieved
programmatically: breadth, human–mouse protein identity (global alignment, BLOSUM62) and accessibility
(UniProt subcellular location / TM / signal peptide). **c**, Per-sample expression of the top four candidates
in the cell type with the strongest effect. **d**, STRING v12 network (combined score ≥ 0.7) among the
candidates and 18 AD GWAS / reference genes.

## Figure S1 | Quality control (`09_qc_figures.py`)
**a–c**, Genes, UMIs and % mitochondrial reads per retained nucleus by sample. **d**, Nuclei retained or removed
per sample. **e**, Median genes per nucleus by condition (Mann–Whitney). **f**, Sex inferred from XIST vs chrY
genes. **g**, ARL17B (17q21.31 copy-number marker) per sample; lines join the two tissues of one donor.

## Figure S2 | Donor clinical covariates (`10_donor_covariates.py`)
Covariates come from Yang *et al.* 2022, Supplementary Table 1, mapped to samples by GEO title index and
verified by sex (17/17) and cortex subset (17/17).
**a, b**, Age and post-mortem interval per donor (Mann–Whitney P). Black circles mark acute hypoxic
neuropathology. **c**, Pericyte hypoxia-module score (mean z of *ADAMTS9, VEGFA, HILPDA, SLC2A3, PDK1, BNIP3L,
ELL2, HIF1A*) against PMI (Spearman). **d**, AD − control difference in the module score (OLS ± 95 % CI) for
each NVU cell type under the primary model, each added covariate, and exclusion of the two hypoxic-pathology
donors. **e**, Share of primary robust DEGs still at FDR < 0.05 when the hippocampal DESeq2 model is re-fitted
with age + PMI, or without the two hypoxic-pathology donors. No gene changes sign in either re-fit.

## Single panels (`single_figures/`)
Each lettered panel above is also exported on its own as `<figure><panel>.png/.pdf` (for example `Fig1a` or `FigS2d`); for Fig. 2 the name is `Fig2<cell type>_<panel>` (for example `Fig2c_e`). Panel letters and titles are removed; the legends above apply unchanged.
