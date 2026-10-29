# Validation of the archived figure-replay software

Version: `v1.0.0`  
Checked: 2026-10-28

## Commands

Run from the repository root:

```bash
Rscript --vanilla equine_multiomics_unsupervised_analysis.R --check
Rscript --vanilla scripts/test_validate_unsupervised.R
Rscript --vanilla equine_multiomics_unsupervised_analysis.R all
Rscript --vanilla scripts/export_unsupervised_tables.R
```

The first two commands use base R. The plotting command uses the packages listed in the main README. The CSV exporter uses base R and performs the same data checks before writing its tables.

## Checks performed

- The archived object schema and required table columns and dimensions.
- Relative variance contributions summing to 100 per condition.
- Significance labels matching their stored unadjusted P values.
- Heatmap feature and sample IDs matching their annotations and matrix order.
- The same 15 horse identities in each saved PCA group and the heatmaps, accounting for the documented `Cushings_`/`PPID_` naming convention.
- Four sets of ten archived AUC values, their stored means and SDs, and clipped mean ± SD bounds.
- Intact inputs accepted and six deliberately altered inputs rejected (annotation order, sample ID, AUC mean, significance label, matrix dimensions and missing object).
- All 19 plot sets exported to PDF, SVG and PNG (57 files).
- All 19 rendered PNGs identical by SHA-256 to the corresponding outputs before the validation and metadata changes.
- All 17 CSV exports written, with explicit matrix/annotation identifiers and stored-run indices for AUC values.

The source `data/unsupervised/plot_data.RData` has SHA-256:

```text
f8c4f1a818f238c3bffcebb6ba736ea1e110239884ac2f01df1e61cfbea54ddb
```

Its bytes and stored numerical results are preserved in this release. The validator inspects the input data and consistency of stored summaries.

[validation-results.json](validation-results.json) records output counts and PNG checksums. [unsupervised-sessionInfo.txt](unsupervised-sessionInfo.txt) records the R environment used for these commands.
