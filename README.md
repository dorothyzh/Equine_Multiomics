# Equine Multi-omics Biomarker Discovery

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.23754598.svg)](https://doi.org/10.5281/zenodo.23754598)

## 1. Introduction

This repository contains code and saved results from a study of plasma proteomics, lipidomics and metabolomics in horses with pituitary pars intermedia dysfunction (PPID). The cohort included eight control horses and seven horses with PPID.

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

### Requirements and quick start

Tested with R 4.4.1 on macOS and ggplot2 4.0.2, tidytext 0.4.3, patchwork 1.3.2, pheatmap 1.0.13, scico 1.5.0, svglite 2.2.2, ragg 1.5.1 and scales 1.4.0. Install the plotting packages in R:

```r
install.packages(c(
  "ggplot2", "tidytext", "patchwork", "pheatmap",
  "scico", "svglite", "ragg", "scales"
))
```

From the repository root, run:

```bash
Rscript --vanilla equine_multiomics_unsupervised_analysis.R all
```

The script writes PDF, editable SVG and 300 dpi PNG files to `outputs/unsupervised/`, along with `sessionInfo.txt`. This folder is ignored by Git.

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
| Figure 3a | Variance contribution by factor, omics and condition; 18 bars | `variance_data` | [Figure 3A-B: variance](equine_multiomics_unsupervised_analysis.R#L100) | `figure3` |
| Figure 3b | Archived covariate associations and significance labels; 117 rows | `covariate_data` | [Figure 3A-B: covariates](equine_multiomics_unsupervised_analysis.R#L112) | `figure3` |
| Figure 3c | Archived top 12 absolute loadings per omics/factor; 108 rows | `loading_top12` | [Figure 3C](equine_multiomics_unsupervised_analysis.R#L127) | `figure3` |
| Figure 3, assembled | Combine panels a–c | Objects for Figure 3a–c | [Figure 3 assembly](equine_multiomics_unsupervised_analysis.R#L148) | `figure3` |
| Figure 4a | Proteomics RF importance; 54 displayed rows | `rf_importance$proteomics` | [Figure 4A-C](equine_multiomics_unsupervised_analysis.R#L153) | `figure4` |
| Figure 4b | Lipidomics RF importance; 44 displayed rows | `rf_importance$lipidomics` | [Figure 4A-C](equine_multiomics_unsupervised_analysis.R#L153) | `figure4` |
| Figure 4c | Metabolomics RF importance; 4 displayed rows | `rf_importance$metabolomics` | [Figure 4A-C](equine_multiomics_unsupervised_analysis.R#L153) | `figure4` |
| Figure 4d | Saved mean AUC with clipped mean ± SD error bars | `auc_summary` | [Figure 4D](equine_multiomics_unsupervised_analysis.R#L171) | `figure4` |
| Figure 4e–f | Saved PCA coordinates before and after adjustment; two facets, 15 horses each | `pca_data` | [Figure 4E-F](equine_multiomics_unsupervised_analysis.R#L185) | `figure4` |
| Figure 4, assembled | Combine RF, AUC and before/after PCA panels | Objects for Figure 4a–f | [Figure 4 assembly](equine_multiomics_unsupervised_analysis.R#L198) | `figure4` |
| Figure 5 | Saved adjusted heatmap; 102 rows × 15 horses | `adjusted_matrix`, `annotation_row`, `annotation_col` | [Figure 5 / heatmap helper](equine_multiomics_unsupervised_analysis.R#L204) | `figure5` |
| Raw heatmap companion | Saved matrix before adjustment; 102 × 15 | `raw_matrix`, `annotation_row`, `annotation_col` | [Raw heatmap](equine_multiomics_unsupervised_analysis.R#L233) | `supplement` |
| Archived model comparison: variance | Average variance by model, factor and omics | `model_variance` | [Model variance](equine_multiomics_unsupervised_analysis.R#L242) | `supplement` |
| Archived model comparison: redundancy | Stored redundant-factor counts | `model_redundancy` | [Model redundancy](equine_multiomics_unsupervised_analysis.R#L247) | `supplement` |
| Archived model comparison: score | Stored model scores | `model_scores` | [Model scores](equine_multiomics_unsupervised_analysis.R#L252) | `supplement` |
| Archived model comparison: total variance | Stored condition-level summed variance values | `model_total_variance` | [Total variance](equine_multiomics_unsupervised_analysis.R#L256) | `supplement` |

The full list of 19 exported plot sets is in the [data and output documentation](data/unsupervised/README.md#output-files).

### Input checks and CSV export

Validate the saved-data structure, sample identifiers, significance labels and AUC summary consistency:

```bash
Rscript --vanilla equine_multiomics_unsupervised_analysis.R --check
```

Export the numerical figure inputs to 17 CSV tables:

```bash
Rscript --vanilla scripts/export_unsupervised_tables.R
# Or choose a destination:
Rscript --vanilla scripts/export_unsupervised_tables.R my_tables
```

The exporter preserves source feature/sample identifiers and stored row order. Its default destination is `outputs/unsupervised/tables/`.

## 4. Versioned archive and citation

Version: [`v1.0.0`](https://github.com/dorothyzh/Equine_Multiomics/releases/tag/v1.0.0). Zenodo DOI: [10.5281/zenodo.23754598](https://doi.org/10.5281/zenodo.23754598).

[Validation results](docs/validation-results.json) · [R environment](docs/unsupervised-sessionInfo.txt)

Citation metadata, including all authors in manuscript order, are provided in [CITATION.cff](CITATION.cff) and [.zenodo.json](.zenodo.json). The repository is connected to Zenodo for GitHub Release archiving.

The source code (R scripts) in this repository is released under the
[MIT License](LICENSE), an OSI-approved open source license.
The data and metadata are released under
[CC BY 4.0](data/LICENSE).

