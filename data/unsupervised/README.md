# Archived plotting data for Figures 3–5

`plot_data.RData` contains the saved values used to draw Figures 3–5. Loading the file creates one list, `plot_data`, with data frames, matrices and AUC vectors. The plotting script reads these values directly.

## Load and inspect

From the repository root:

```r
saved <- new.env()
load("data/unsupervised/plot_data.RData", envir = saved)
plot_data <- saved$plot_data
names(plot_data)
head(plot_data$variance_data)
dim(plot_data$adjusted_matrix)
```

The RData file uses R serialization version 2 and xz compression. Its SHA-256 is:

```text
f8c4f1a818f238c3bffcebb6ba736ea1e110239884ac2f01df1e61cfbea54ddb
```

## Object dictionary

| Element of `plot_data` | Size | Stored content | Used by |
|------------------------|------|----------------|---------|
| `schema_version` | Scalar integer `1L` | Plot-data format version | Input check |
| `variance_data` | 18 × 6 | `factor`, `view`, `value`, `group`, `total_variance`, `percentage`; saved relative variance contributions | Figure 3a |
| `covariate_data` | 117 × 12 | Covariate/interaction (`var`), condition (`Status`), omics (`Modality`), `factor`, archived model/test fields, `p.value`, `logp` and significance `label` | Figure 3b |
| `loading_top12` | 108 × 4 | `feature`, `factor`, signed loading `value`, `view`; fixed top 12 absolute loadings for each of 3 omics × 3 factors | Figure 3c, including separate factor displays |
| `rf_importance` | Named list of 3 tables | `proteomics`: 54 rows; `lipidomics`: 44; `metabolomics`: 4. Eight columns: `Control`, `Cushings`, `MeanDecreaseAccuracy`, `MeanDecreaseGini`, `Feature`, `Factor`, `Group`, `Omics` | Figure 4a–c |
| `auc_summary` | 4 × 5 | `Factor`, `Mean_AUC`, `SD_AUC`, `Upper`, `Lower` | Figure 4d |
| `auc_cv_results` | 4 named vectors, 10 values each | Saved AUC values for Factor1, Factor2, Factor3 and Combined; retained as provenance, not used to recalculate Figure 4d | Reference only |
| `pca_data` | 30 × 8 | Saved `PC1`, `PC2`, `sample`, `Status`, `Sex`, `Age`, `Batch`, `Type`; 15 horses before and the same 15 after adjustment | Figure 4e–f |
| `adjusted_matrix` | 102 × 15 numeric matrix | Saved adjusted heatmap values, with feature and sample IDs | Figure 5 |
| `raw_matrix` | 102 × 15 numeric matrix | Saved, already scaled heatmap values before adjustment; this is not raw instrument data | Raw heatmap companion |
| `annotation_row` | 102 × 3 | `Polarity`, `Factor`, `Omics`; row names match matrix feature IDs | Both heatmaps |
| `annotation_col` | 15 × 3 | `Sex`, `Age`, `Batch`; row names match matrix sample IDs | Both heatmaps |
| `model_variance` | 108 × 4 | `Model`, `Factor`, `Omics`, `AverageVariance` | Archived variance heatmap |
| `model_redundancy` | 8 × 2 | `Model`, `NumRedundant` from the saved redundant-factor plot | Archived redundancy bars |
| `model_scores` | 8 × 5 | `Model`, `TotalVariance`, `VarianceSD`, `NumRedundant`, `Score` from the saved score plot | Archived model-score bars |
| `model_total_variance` | 16 × 4 | `Model`, `Condition`, `Omics`, `TotalVariance`; the archived `Omics` column is missing (`NA`) throughout | Archived total-variance bars |

The heatmap has 102 rows and 15 columns, representing eight control horses and seven horses with PPID. PEF1 appears in two factor groups, stored as `PEF1` and `PEF1.1`; both rows have identical values and are displayed as PEF1. There are 101 distinct display names. `Polarity` denotes loading direction, not ionization mode.

Matrix values, row order and sample order are preserved. The script sets `cluster_rows = FALSE`, `cluster_cols = FALSE` and leaves `scale = "none"` in `pheatmap`. The adjusted matrix is already standardized by row.

Figure 4d uses the stored summary columns exactly. `Lower = max(Mean_AUC - SD_AUC, 0)` and `Upper = min(Mean_AUC + SD_AUC, 1)`. These are clipped SD bounds, not 95% confidence intervals. Mean values in Factor1, Factor2, Factor3, Combined order are 0.675, 0.83125, 0.8625, 0.896875. The plotting script neither pools nor recomputes the ten stored AUC values.

The model-comparison plot tables are separate archived display objects. In particular, the redundancy counts in `model_redundancy` and `model_scores` are not interchangeable and are not recalculated here. Missing omics labels in `model_total_variance` are retained; its plot uses only the existing model and condition fields.

## Provenance and extraction

The extraction was prepared on 2026-10-28 from the following existing project results:

| Current element | Original saved source |
|-----------------|-----------------------|
| `variance_data` | `adjust_heatmap.RData`: `plotgg$data` |
| `covariate_data` | `adjust_heatmap.RData`: `p3$data` |
| `pca_data` | `adjust_heatmap.RData`: `p_all$data`, restricted to the eight plotting/ID columns listed above |
| `auc_summary`, `auc_cv_results` | Same-named objects in `adjust_heatmap.RData` |
| `adjusted_matrix` | `adjust_heatmap.RData`: `input_matrix` |
| `raw_matrix` | `adjust_heatmap.RData`: `sample_data_scaled_matrix_raw` |
| `annotation_row`, `annotation_col` | Same-named objects in `adjust_heatmap.RData` |
| `model_variance` | `adjust_heatmap.RData`: `heatmap_plot$data` |
| `model_redundancy` | `adjust_heatmap.RData`: `redundant_factors_plot$data` |
| `model_scores` | `adjust_heatmap.RData`: `model_scores_plot$data` |
| `model_total_variance` | `adjust_heatmap.RData`: `total_variance_plot$data` |
| `loading_top12` | Archived loading values from full-feature `0120/model.hdf5`; the existing top-12 absolute-loading display rule was applied once during extraction, separately for each view/factor |
| `rf_importance` | Existing `0211/<omics>/importance_top_factor*_<positive-or-negative>.csv` files; the existing top-10 importance display rule was applied once per factor/polarity, retaining ties |

The original `adjust_heatmap.RData` archive has MD5 `fe58f0084a535a0d7f86f94fb36511a3`. Only the listed plotting data were extracted. Values were not refitted or rescaled, and the original workspace and HDF5 files are not needed to draw the figures.

Figure 3 traces to the full-feature model (834 proteins, 3,906 metabolomic features and 1,619 lipidomic features). The saved RF selection traces to a reduced-feature model (834 proteins, 115 metabolomic features and 625 lipidomic features). Both archived MOFA models have 45 modality-specific IDs for 15 horses, without shared observed IDs across views. The compact plotting data preserve the historical outputs; they do not repair this sample-pairing issue or reconcile the two model lineages. See the [notes in the main README](../../README.md#notes-on-the-archived-results).

## Output files

The main script is [`equine_multiomics_unsupervised_analysis.R`](../../equine_multiomics_unsupervised_analysis.R). Each filename below is written with `.pdf`, `.svg` and `.png` extensions under `outputs/unsupervised/`. `all` generates all 19 sets; the corresponding individual command generates only its group.

| Run argument | Relative output stem | Content |
|--------------|----------------------|---------|
| `figure3` | `figures/Figure3_replay` | Assembled Figure 3 |
| `figure3` | `panels/Figure3A_variance` | Variance bars |
| `figure3` | `panels/Figure3B_covariates` | Covariate associations |
| `figure3` | `panels/Figure3C_loadings` | All nine loading panels |
| `figure3` | `panels/Figure3C_Factor1` | Factor 1 loadings |
| `figure3` | `panels/Figure3C_Factor2` | Factor 2 loadings |
| `figure3` | `panels/Figure3C_Factor3` | Factor 3 loadings |
| `figure4` | `figures/Figure4_replay` | Assembled Figure 4 |
| `figure4` | `panels/Figure4_RF_proteomics` | Proteomics importance |
| `figure4` | `panels/Figure4_RF_lipidomics` | Lipidomics importance |
| `figure4` | `panels/Figure4_RF_metabolomics` | Metabolomics importance |
| `figure4` | `panels/Figure4D_AUC` | Stored AUC summary |
| `figure4` | `panels/Figure4EF_PCA` | Before and adjusted PCA facets |
| `figure5` | `figures/Figure5_adjusted_replay` | Adjusted heatmap |
| `supplement` | `supplement/Supplement_raw_heatmap` | Before-adjustment heatmap |
| `supplement` | `supplement/Supplement_heatmap_plot` | Archived model variance |
| `supplement` | `supplement/Supplement_redundant_factors_plot` | Archived redundancy counts |
| `supplement` | `supplement/Supplement_model_scores_plot` | Archived model scores |
| `supplement` | `supplement/Supplement_total_variance_plot` | Archived total variance |

Generated PDFs, SVGs, PNGs and `sessionInfo.txt` are local outputs and are ignored by the supplied `outputs/unsupervised/.gitignore`. The data are covered by the repository's existing [CC BY 4.0 data license](../LICENSE).
