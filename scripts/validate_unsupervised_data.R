# Validate the archived figure inputs without fitting any model or changing values.
# Sourced by the plotting entry point; uses base R only.
validate_unsupervised_data <- function(b) {
  check <- function(ok, message) {
    if (!isTRUE(ok)) stop("Archived plotting data: ", message, call. = FALSE)
  }
  table_check <- function(name, columns, rows) {
    x <- b[[name]]
    check(is.data.frame(x), paste(name, "must be a data frame"))
    check(all(columns %in% names(x)), paste(name, "has missing columns"))
    check(nrow(x) == rows, paste(name, "has an unexpected row count"))
  }
  finite <- function(x) is.numeric(x) && all(is.finite(x))
  check(is.list(b) && identical(b$schema_version, 1L), "unsupported schema_version")
  required <- c("variance_data", "covariate_data", "loading_top12", "rf_importance",
    "auc_summary", "auc_cv_results", "pca_data", "adjusted_matrix", "raw_matrix",
    "annotation_row", "annotation_col", "model_variance", "model_redundancy",
    "model_scores", "model_total_variance")
  check(all(required %in% names(b)), "required objects are missing")
  table_check("variance_data", c("factor", "view", "percentage", "group"), 18L)
  v <- b$variance_data
  check(finite(v$percentage) && all(v$percentage >= 0), "invalid variance percentages")
  check(all(abs(tapply(v$percentage, v$group, sum) - 100) < 1e-8),
    "relative variance contributions must sum to 100 per condition")
  table_check("covariate_data", c("var", "Status", "Modality", "factor", "p.value", "logp", "label"), 117L)
  c <- b$covariate_data
  check(finite(c$p.value) && all(c$p.value >= 0 & c$p.value <= 1), "invalid association P values")
  check(finite(c$logp), "invalid association colour values")
  stars <- ifelse(c$p.value < .001, "***", ifelse(c$p.value < .01, "**", ifelse(c$p.value < .05, "*", "")))
  check(identical(stars, as.character(c$label)), "significance labels differ from stored unadjusted P values")
  table_check("loading_top12", c("feature", "factor", "value", "view"), 108L)
  check(finite(b$loading_top12$value), "invalid loadings")
  check(all(table(b$loading_top12$view, b$loading_top12$factor) == 12L),
    "expected 12 saved loadings per view and factor")
  views <- c("proteomics", "lipidomics", "metabolomics")
  check(identical(names(b$rf_importance), views), "RF view names/order differ")
  for (i in seq_along(views)) {
    x <- b$rf_importance[[i]]
    check(is.data.frame(x) && nrow(x) == c(54L, 44L, 4L)[i] &&
      all(c("Feature", "MeanDecreaseAccuracy", "Factor", "Group") %in% names(x)),
      paste("unexpected RF table:", views[i]))
    check(finite(x$MeanDecreaseAccuracy), paste("invalid RF importance:", views[i]))
  }
  table_check("auc_summary", c("Factor", "Mean_AUC", "SD_AUC", "Upper", "Lower"), 4L)
  a <- b$auc_summary
  factors <- c("Factor1", "Factor2", "Factor3", "Combined")
  check(identical(as.character(a$Factor), factors) && identical(names(b$auc_cv_results), factors),
    "AUC factor identities/order differ")
  check(all(lengths(b$auc_cv_results) == 10L), "expected 10 archived AUC values per group")
  check(all(vapply(b$auc_cv_results, function(x) finite(x) && all(x >= 0 & x <= 1), logical(1))),
    "invalid archived AUC values")
  near <- function(x, y) finite(x) && all(abs(x - y) < 1e-12)
  check(near(a$Mean_AUC, vapply(b$auc_cv_results, mean, numeric(1))) &&
    near(a$SD_AUC, vapply(b$auc_cv_results, sd, numeric(1))), "AUC summaries differ from saved runs")
  check(near(a$Upper, pmin(a$Mean_AUC + a$SD_AUC, 1)) &&
    near(a$Lower, pmax(a$Mean_AUC - a$SD_AUC, 0)), "AUC bounds differ from clipped mean +/- SD")
  for (name in c("adjusted_matrix", "raw_matrix")) {
    m <- b[[name]]
    check(is.matrix(m) && identical(dim(m), c(102L, 15L)) && finite(m), paste("invalid", name))
    check(!is.null(rownames(m)) && !anyDuplicated(rownames(m)) &&
      !is.null(colnames(m)) && !anyDuplicated(colnames(m)), paste(name, "needs unique row and sample IDs"))
    check(identical(dimnames(m), dimnames(b$adjusted_matrix)), "heatmap IDs/order differ")
  }
  table_check("annotation_row", c("Polarity", "Factor", "Omics"), 102L)
  table_check("annotation_col", c("Sex", "Age", "Batch"), 15L)
  check(identical(rownames(b$annotation_row), rownames(b$adjusted_matrix)) &&
    identical(rownames(b$annotation_col), colnames(b$adjusted_matrix)), "heatmap annotations are misaligned")
  check(!anyNA(b$annotation_row) && !anyNA(b$annotation_col), "heatmap annotations contain missing values")
  table_check("pca_data", c("PC1", "PC2", "sample", "Status", "Sex", "Age", "Batch", "Type"), 30L)
  p <- b$pca_data
  check(finite(p$PC1) && finite(p$PC2) && setequal(unique(p$Type), c("Before", "Adjusted1")), "invalid PCA coordinates/groups")
  for (type in c("Before", "Adjusted1")) {
    # The PCA archive retains Cushings_ IDs; the heatmaps use PPID_ for the same horses.
    ids <- sub("^Cushings_", "PPID_", p$sample[p$Type == type])
    check(length(ids) == 15L && !anyDuplicated(ids) && setequal(ids, colnames(b$adjusted_matrix)),
      paste("PCA sample identities differ for", type))
  }
  table_check("model_variance", c("Model", "Factor", "Omics", "AverageVariance"), 108L)
  table_check("model_redundancy", c("Model", "NumRedundant"), 8L)
  table_check("model_scores", c("Model", "Score"), 8L)
  table_check("model_total_variance", c("Model", "Condition", "Omics", "TotalVariance"), 16L)
  # Missing Omics in this historical table is documented and deliberately retained.
  check(all(is.na(b$model_total_variance$Omics)), "archived total-variance Omics labels changed")
  invisible(TRUE)
}
