# Equine Multi-omics Biomarker Discovery

## 1. Introduction

This repository contains code and saved results from a study of plasma proteomics, lipidomics and metabolomics in horses with pituitary pars intermedia dysfunction (PPID). The cohort included eight control horses and seven horses with PPID.

The supervised analysis compares molecular features between groups. The unsupervised plotting script redraws Figures 3–5 from saved results. Known issues in the archived MOFA analysis are described in Section 3.

## 2. Supervised Biomarker Discovery

**Script:** [`equine_multiomics_supervised_analysis.R`](equine_multiomics_supervised_analysis.R)

Supervised biomarker discovery is done by differential expression analysis using [limma](https://bioconductor.org/packages/limma/) package.
Other required packages are outlined in script. 
Preprocessing steps are detailed in line 19-105 of the script.

| Figure | Description | Section of script |
|--------|-------------|-------------------|
| Figure 2a-d | lipidomics, metabolomics, proteomics, and multi-omics combined PCA | Line 108-193 |
| Figure 2e | multi-omics heatmap | Line 195-352 |
| Figure S1 | multi-omics volcano plot | Line 354-379 |
| Figure S3 | multi-omics volcano plot | Line 381-503 |


## 3. Unsupervised Biomarker Discovery

**Script:** [`equine_multiomics_unsupervised_analysis.R`](equine_multiomics_unsupervised_analysis.R)  
**Saved plotting data:** [`data/unsupervised/plot_data.RData`](data/unsupervised/plot_data.RData)  
**Data dictionary and provenance:** [`data/unsupervised/README.md`](data/unsupervised/README.md)

The script loads one RData file and draws Figures 3–5 from saved tables and matrices. Loadings, RF importance scores, AUC summaries, PCA coordinates and heatmap values are already stored. The script handles plotting and file export only.

### Requirements and quick start

Tested with R 4.4.1 on macOS and ggplot2 4.0.2, tidytext 0.4.3, patchwork 1.3.2, pheatmap 1.0.13, scico 1.5.0, svglite 2.2.2 and ragg 1.5.1. Install the plotting packages in R:

```r
install.packages(c(
  "ggplot2", "tidytext", "patchwork", "pheatmap",
  "scico", "svglite", "ragg"
))
```

From the repository root, run:

```bash
Rscript --vanilla equine_multiomics_unsupervised_analysis.R all
```

The script writes PDF, editable SVG and 300 dpi PNG files to `outputs/unsupervised/`, along with `sessionInfo.txt`. This folder is ignored by Git. All plotting inputs are included in the RData file; the original workspaces and model files are not needed.

Individual figure groups can be generated with:

```bash
Rscript --vanilla equine_multiomics_unsupervised_analysis.R figure3
Rscript --vanilla equine_multiomics_unsupervised_analysis.R figure4
Rscript --vanilla equine_multiomics_unsupervised_analysis.R figure5
Rscript --vanilla equine_multiomics_unsupervised_analysis.R supplement
```

An optional second argument sets the output directory:

```bash
Rscript --vanilla equine_multiomics_unsupervised_analysis.R figure3 my_figures
```

The script locates its RData file relative to the script, so an absolute script path also works from another working directory. In RStudio, open the repository as the working directory and run `source("equine_multiomics_unsupervised_analysis.R")` to generate all figures.

### Figures, data and code

All data objects in the table belong to `plot_data`, loaded from `data/unsupervised/plot_data.RData`. The last column gives the argument to pass after the script name.

| Figure | Description | Saved data object | Section of script | Run argument |
|--------|-------------|-------------------|-------------------|--------------|
| Figure 3a | Variance contribution by factor, omics and condition; 18 bars | `variance_data` | [Figure 3A-B: variance](equine_multiomics_unsupervised_analysis.R#L93) | `figure3` |
| Figure 3b | Archived covariate associations and significance labels; 117 rows | `covariate_data` | [Figure 3A-B: covariates](equine_multiomics_unsupervised_analysis.R#L105) | `figure3` |
| Figure 3c | Archived top 12 absolute loadings per omics/factor; 108 rows | `loading_top12` | [Figure 3C](equine_multiomics_unsupervised_analysis.R#L120) | `figure3` |
| Figure 3, assembled | Combine panels a–c | Objects for Figure 3a–c | [Figure 3 assembly](equine_multiomics_unsupervised_analysis.R#L141) | `figure3` |
| Figure 4a | Proteomics RF importance; 54 displayed rows | `rf_importance$proteomics` | [Figure 4A-C](equine_multiomics_unsupervised_analysis.R#L146) | `figure4` |
| Figure 4b | Lipidomics RF importance; 44 displayed rows | `rf_importance$lipidomics` | [Figure 4A-C](equine_multiomics_unsupervised_analysis.R#L146) | `figure4` |
| Figure 4c | Metabolomics RF importance; 4 displayed rows | `rf_importance$metabolomics` | [Figure 4A-C](equine_multiomics_unsupervised_analysis.R#L146) | `figure4` |
| Figure 4d | Saved mean AUC with clipped mean ± SD error bars | `auc_summary` | [Figure 4D](equine_multiomics_unsupervised_analysis.R#L164) | `figure4` |
| Figure 4e–f | Saved PCA coordinates before and after adjustment; two facets, 15 horses each | `pca_data` | [Figure 4E-F](equine_multiomics_unsupervised_analysis.R#L178) | `figure4` |
| Figure 4, assembled | Combine RF, AUC and before/after PCA panels | Objects for Figure 4a–f | [Figure 4 assembly](equine_multiomics_unsupervised_analysis.R#L191) | `figure4` |
| Figure 5 | Saved adjusted heatmap; 102 rows × 15 horses | `adjusted_matrix`, `annotation_row`, `annotation_col` | [Figure 5 / heatmap helper](equine_multiomics_unsupervised_analysis.R#L197) | `figure5` |
| Raw heatmap companion | Saved matrix before adjustment; 102 × 15 | `raw_matrix`, `annotation_row`, `annotation_col` | [Raw heatmap](equine_multiomics_unsupervised_analysis.R#L226) | `supplement` |
| Archived model comparison: variance | Average variance by model, factor and omics | `model_variance` | [Model variance](equine_multiomics_unsupervised_analysis.R#L235) | `supplement` |
| Archived model comparison: redundancy | Stored redundant-factor counts | `model_redundancy` | [Model redundancy](equine_multiomics_unsupervised_analysis.R#L240) | `supplement` |
| Archived model comparison: score | Stored model scores | `model_scores` | [Model scores](equine_multiomics_unsupervised_analysis.R#L245) | `supplement` |
| Archived model comparison: total variance | Stored condition-level summed variance values | `model_total_variance` | [Total variance](equine_multiomics_unsupervised_analysis.R#L249) | `supplement` |

The additional archived plots have descriptive filenames; no final supplementary figure numbers are assigned here. The full list of 19 exported plot sets is in the [data and output documentation](data/unsupervised/README.md#output-files).

### Notes on the archived results

The archived MOFA analysis encoded 15 horses as 45 modality-specific sample IDs, with no observed IDs shared across omics views. Figure 3 comes from the full-feature model, `0120/model.hdf5`. The RF feature groups come from a reduced-feature model, `mofa_testing/mofa3/model.hdf5`. Their factor numbers should not be assumed to identify the same factors. The downstream RF inputs, PCA and heatmap still represent 15 horses. Redrawing these results does not resolve the model issues.

Several details also differ from the manuscript text or layout:

- Figure 3a shows relative contributions, normalized to sum to 100% within each condition. The six totals quoted in the Results match the saved values.
- Figure 3b shows association P values, not percentages of variance explained. The saved stars use unadjusted P values.
- Figure 3c contains 12 features per omics/factor panel, as in the archived figure. The manuscript caption says 10. Positive and negative loadings describe factor direction, not disease-related increases or decreases.
- Figure 4d preserves AUC means of 0.675, 0.83125, 0.8625 and 0.896875. Its bounds are mean ± SD, clipped to [0, 1], rather than 95% confidence intervals. Ten AUC values per group are saved; the manuscript image says 100 runs.
- Figure 4e–f uses saved PCA coordinates from 834 proteins, with 15 horses before and after adjustment. In the current combined replay, both facets share panel label e; the manuscript labels the adjusted facet f separately.
- Figure 5 uses the saved adjusted matrix. All 1,530 heatmap cell colors in the manuscript image match this matrix under the original color scale; the image title incorrectly says “Raw Data Version”. Rows are grouped by omics, factor and loading direction, with clustering disabled. The 102 displayed rows include PEF1 twice in different factor groups.
- Source IDs are retained in the RData file. Display labels use PPID for `Cushings`, Adjusted for `Adjusted1`, and PEF1 for `PEF1.1`. The saved metabolite name `TETRAHYDROURIDINE` differs from the manuscript labels TRU and TRU.isomer; the annotation should be checked before renaming it.
- The archived total-variance table has no usable omics labels, so its existing condition-level values are plotted together.

The script preserves numerical results while rebuilding the layout. These outputs are for checking the archived figures, not for confirming the manuscript's biological or predictive claims.

The source code (R scripts) in this repository is released under the
[MIT License](LICENSE), an OSI-approved open source license.
The data and metadata are released under
[CC BY 4.0](data/LICENSE).

