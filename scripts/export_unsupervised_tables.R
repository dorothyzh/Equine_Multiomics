# Export the archived figure inputs as CSV; no fitted values are recalculated.
# Usage: Rscript --vanilla scripts/export_unsupervised_tables.R [output_directory]
options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (!length(file_arg)) stop("Run this exporter with Rscript.")
root <- dirname(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
out <- if (length(args)) args[1] else file.path(root, "outputs", "unsupervised", "tables")
saved <- new.env(parent = emptyenv())
loaded <- load(file.path(root, "data", "unsupervised", "plot_data.RData"), envir = saved)
if (!identical(loaded, "plot_data")) stop("Expected exactly one archived object named plot_data.")
source(file.path(root, "scripts", "validate_unsupervised_data.R"))
b <- saved$plot_data
validate_unsupervised_data(b)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
write_table <- function(x, name, row_id = NULL) {
  if (!is.null(row_id)) {
    ids <- rownames(x)
    x <- data.frame(ids, x, check.names = FALSE)
    names(x)[1] <- row_id
  }
  write.csv(x, file.path(out, paste0(name, ".csv")), row.names = FALSE, na = "NA")
}
for (name in c("variance_data", "covariate_data", "loading_top12", "auc_summary", "pca_data",
  "model_variance", "model_redundancy", "model_scores", "model_total_variance")) write_table(b[[name]], name)
for (name in names(b$rf_importance)) write_table(b$rf_importance[[name]], paste0("rf_importance_", name))
for (name in c("adjusted_matrix", "raw_matrix", "annotation_row")) write_table(b[[name]], name, "source_feature_id")
write_table(b$annotation_col, "annotation_col", "sample_id")
runs <- do.call(rbind, lapply(names(b$auc_cv_results), function(name) data.frame(
  Factor = name, stored_run_index = seq_along(b$auc_cv_results[[name]]), AUC = b$auc_cv_results[[name]])))
write_table(runs, "auc_cv_results")
cat("Exported 17 CSV tables to ", normalizePath(out), "\n", sep = "")
