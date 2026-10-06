# Single-Nucleus Quality Control (QC) & Peer-Review Defense Guide

## 🔬 Peer-Review Overview & Objective

In single-cell and single-nucleus RNA sequencing (snRNA-seq) studies of neurodegenerative disease, peer reviewers routinely scrutinize quality control protocols for potential technical artifacts, including:
1. **Confounding by Disease State**: Did neurodegenerative tissue from Alzheimer's Disease (AD) patients suffer greater cellular breakdown, introducing false-positive differential expression?
2. **Doublet Contamination**: Were nuclear multiplets / doublets removed to avoid spurious hybrid cell states?
3. **Ambient RNA & Empty Droplets**: Were low-complexity droplets stripped of ambient RNA?
4. **Mitochondrial Fraction in Nuclei**: Why was a 5% mitochondrial threshold selected for single-nucleus sequencing?

This dedicated quality control suite ([`figures/Quality_Control/`](.)) provides rigorous empirical evidence, statistical tests, and visual documentation to completely rebut reviewer critique.

---

## 📊 Figure Suite Catalog

### [Fig_QC_01: Global Quality Control Metrics by Condition & Region](Fig_QC_01_metrics_violin_by_condition.png)
- **Panels**:
  - `nFeature_RNA` (Detected genes per nucleus)
  - `nCount_RNA` (Total unique molecular identifiers [UMIs] per nucleus)
  - `percent.mt` (Mitochondrial transcript percentage)
  - `percent.ribo` (Ribosomal transcript percentage)
- **Reviewer Defense Takeaway**:
  - Evaluated across **95,051 total nuclei** (**40,000 Cortex**, **55,051 Hippocampus**).
  - Wilcoxon rank-sum statistical annotations confirm that library complexity and quality metrics show no systematic bias between Control and AD donors, proving downstream DEGs reflect true biology rather than technical degradation.

### [Fig_QC_02: Bivariate Quality Metrics Correlation & Thresholds](Fig_QC_02_metrics_bivariate_scatter_and_thresholds.png)
- **Panels**:
  - **Left**: `nCount_RNA` vs `nFeature_RNA` bivariate distribution with linear regression and Pearson correlation coefficient ($R = 0.90$, $p < 2.2 \times 10^{-16}$). Demonstrates high library complexity and clear separation from empty droplets.
  - **Right**: `nCount_RNA` vs `percent.mt` with red dashed line at the 5.0% threshold.
- **Reviewer Defense Takeaway**:
  - Confirms exclusion of apoptotic/damaged nuclei and defines the bivariate boundary eliminating doublets ($> 6,500$ genes / $> 25,000$ UMIs) and debris ($< 300$ genes / $< 500$ UMIs).

### [Fig_QC_03: Donor-by-Donor Batch Variability (All 25 Samples)](Fig_QC_03_donor_batch_variability_violins.png)
- **Panels**:
  - Sample-resolved distribution of `nFeature_RNA` and `percent.mt` across all 25 donors (GSM4982083–GSM4982107).
- **Reviewer Defense Takeaway**:
  - Demonstrates inter-donor uniformity across frozen brain batches without single-sample dropout or outlier sequencing failure.

### [Fig_QC_04: Lineage-Specific Quality Distributions across 6 NVU Cell Types](Fig_QC_04_celltype_quality_metrics.png)
- **Panels**:
  - Gene count, UMI count, and MT% across Astrocytes, Endothelial cells, Inhibitory neurons, Microglia, Oligodendrocytes, and Pericytes.
- **Reviewer Defense Takeaway**:
  - Small vascular and immune cell types (Endothelial cells, Pericytes, Microglia) have intrinsically lower cytoplasmic volume than large neurons, yet exhibit tight, high-quality gene detection distributions well above noise thresholds.

### [Fig_QC_05: Nuclei Retention & Filtering Thresholds Summary](Fig_QC_05_cell_retention_and_filtering_summary.png)
- **Panels**:
  - Left: Retention rate (%) post-QC per individual donor sample.
  - Right: Standard Operating Procedure (SOP) threshold summary table for peer review.
- **Reviewer Defense Takeaway**:
  - High mean post-QC retention rate demonstrates minimal data loss while filtering out all technical artifacts.

---

## 📋 Quantitative Quality Control Thresholds Applied

| Parameter | Threshold | Rationale for Peer Reviewers |
|---|---|---|
| **Detected Genes (`nFeature_RNA`)** | $300 \le \text{Genes} \le 6,500$ | Lower cutoff eliminates empty droplets and ambient RNA; upper cutoff removes nuclear doublets/multiplets. |
| **Total UMIs (`nCount_RNA`)** | $500 \le \text{UMIs} \le 25,000$ | Eliminates low-yield cellular debris while censoring PCR hyper-amplification artifacts. |
| **Mitochondrial Fraction (`percent.mt`)** | $\le 5.0\%$ | **snRNA-seq Standard**: Single-cell suspensions from fresh tissue allow 10–20% MT due to intact cytoplasm; single-nucleus suspensions from frozen postmortem tissue should have $\le 5\%$ MT because mitochondria are cytoplasmic organelles excluded during nuclear isolation. |
| **Doublet & Outlier Handling** | Bivariate filtering + Harmony | Graph-based clustering and Harmony batch alignment prevent batch-specific artificial subclusters. |

---

## 💬 Sample Responses to Common Reviewer Questions

> **Reviewer Question 1**: *"Did postmortem tissue degradation in Alzheimer's disease donors cause reduced library quality or artificially elevated mitochondrial read counts compared to age-matched controls?"*
> 
> **Author Response**: *"We quantitatively addressed this concern by profiling all quality control metrics across clinical conditions (Fig_QC_01). In both the Prefrontal Cortex and Hippocampus, detected gene counts (nFeature_RNA), total UMI counts (nCount_RNA), and mitochondrial read percentages showed equivalent distributions between AD and Control donors (Wilcoxon rank-sum test, $p > 0.05$). Furthermore, donor-by-donor violin profiling across all 25 individual samples (Fig_QC_03) demonstrated consistent sequencing depth and library complexity across batches."*

---

> **Reviewer Question 2**: *"How were doublets and empty droplets separated from genuine nuclei?"*
> 
> **Author Response**: *"Nuclei were filtered using stringent bivariate boundaries (Fig_QC_02). Droplets with fewer than 300 detected genes or fewer than 500 total UMIs were removed to eliminate ambient RNA contamination. Potential homotypic and heterotypic nuclear doublets were eliminated by applying an upper threshold of 6,500 detected genes and 25,000 total UMIs, followed by Harmony integration to verify that cluster topology was driven by cell lineage rather than library size."*
