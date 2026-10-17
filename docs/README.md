# Documentation

| Document | Content |
|---|---|
| [`01_figure_interpretations.md`](01_figure_interpretations.md) | For each figure (S1, 1, 2a–e, 3, 4, 5): what it shows, why the analysis was chosen, and what can be concluded |
| [`02_wet_lab_validation_plan.md`](02_wet_lab_validation_plan.md) | Staged, gated plan: metadata and replication → in-situ human validation → in vitro mechanism → functional tests → mouse |
| [`03_methods.md`](03_methods.md) | Methods for every pipeline step (01–10): parameters, models, statistics, software versions |
| [`../figures/README.md`](../figures/README.md) | Panel-by-panel figure legends |
| [`../legacy/README.md`](../legacy/README.md) | Defects of the superseded R pipeline and how each was fixed |

Supporting tables are all produced by `scripts/10_donor_covariates.py`:
- `tables/donor_clinical_metadata.csv`: donor clinical data (Yang *et al.* Supplementary Table 1), mapped and verified
- `tables/covariate_balance.csv`: AD vs control comparison of each covariate
- `tables/hypoxia_module_score_hpc.csv`: per-sample hypoxia-module score, all 5 NVU cell types
- `tables/hypoxia_module_covariate_models.csv`: AD effect on the score under each adjustment
- `tables/pseudobulk_DE_sensitivity_clinical.csv`: hippocampal DE re-fitted with age + PMI, or without hypoxic-pathology donors
- `tables/validation_power_estimates.csv`: per-donor effect sizes (Cohen's d) and sample-size estimates
