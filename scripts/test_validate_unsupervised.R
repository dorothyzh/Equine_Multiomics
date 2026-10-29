# Regression checks for detection of misaligned or inconsistent saved figure inputs.
# Usage: Rscript --vanilla scripts/test_validate_unsupervised.R
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
if (!length(file_arg)) stop("Run these checks with Rscript.")
root <- dirname(dirname(normalizePath(sub("^--file=", "", file_arg[1]))))
source(file.path(root, "scripts", "validate_unsupervised_data.R"))
saved <- new.env(parent = emptyenv())
load(file.path(root, "data", "unsupervised", "plot_data.RData"), envir = saved)
b <- saved$plot_data
stopifnot(isTRUE(validate_unsupervised_data(b)))
reject <- function(x, text) {
  error <- tryCatch({ validate_unsupervised_data(x); NULL }, error = identity)
  stopifnot(inherits(error, "error"), grepl(text, conditionMessage(error), fixed = TRUE))
}
x <- b
x$annotation_col <- x$annotation_col[15:1, , drop = FALSE]
reject(x, "heatmap annotations are misaligned")
x <- b
x$pca_data$sample[1] <- "unknown_horse"
reject(x, "PCA sample identities differ")
x <- b
x$auc_summary$Mean_AUC[1] <- x$auc_summary$Mean_AUC[1] + 0.01
reject(x, "AUC summaries differ")
x <- b
x$covariate_data$label[1] <- "***"
reject(x, "significance labels differ")
x <- b
x$adjusted_matrix <- x$adjusted_matrix[-1, , drop = FALSE]
reject(x, "invalid adjusted_matrix")
x <- b
x$auc_cv_results <- NULL
reject(x, "required objects are missing")
cat("PASS: intact archive accepted; six deliberately altered inputs rejected.\n")
