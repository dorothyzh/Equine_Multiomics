# Equine multi-omics: plot archived unsupervised results from RData.
# Usage: Rscript
# equine_multiomics_unsupervised_analysis.R [all|figure3|figure4|figure5|supplement]
# [output_directory]
options(stringsAsFactors = FALSE)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) && args[1] %in% c("-h", "--help")) {
  cat("Usage: Rscript equine_multiomics_unsupervised_analysis.R [all|figure3|figure4|figure5|supplement|--check] [output_directory]\n")
  quit(status = 0)
}
target <- if (length(args)) args[1] else "all"
if (!target %in% c("all", "figure3", "figure4", "figure5", "supplement", "--check")) {
  stop("Unknown figure: ", target, ". Use --help for available commands.")
}
file_arg <- grep("^--file=", commandArgs(), value = TRUE)
root <- if (length(file_arg)) dirname(normalizePath(sub("^--file=", "", file_arg[1]))) else getwd()
saved <- new.env(parent = emptyenv())
loaded <- load(file.path(root, "data", "unsupervised", "plot_data.RData"), envir = saved)
if (!identical(loaded, "plot_data")) stop("Expected exactly one archived object named plot_data.")
b <- saved$plot_data
source(file.path(root, "scripts", "validate_unsupervised_data.R"), local = TRUE)
validate_unsupervised_data(b)
if (target == "--check") {
  cat("PASS: archived figure-input structure, sample alignment and stored AUC summaries.\n")
  cat("Checked: saved-data schema, heatmap/PCA sample identities, significance labels and mean +/- SD bounds.\n")
  quit(status = 0)
}
packages <- c("ggplot2", "tidytext", "patchwork", "pheatmap", "scico", "svglite", "ragg", "scales")
missing <- packages[!vapply(packages, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) stop("Install these R packages first: ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({
  library(ggplot2)
  library(tidytext)
  library(patchwork)
  library(pheatmap)
  library(grid)
  library(scico)
  library(svglite)
  library(ragg)
})
output_root <- if (length(args) > 1L) args[2] else file.path(root, "outputs", "unsupervised")
for (d in c("figures", "panels", "supplement")) {
  dir.create(file.path(output_root, d), showWarnings = FALSE, recursive = TRUE)
}
font_family <- "Helvetica"
open_pdf <- function(path, w, h) {
  grDevices::pdf(path, width = w, height = h, family = font_family, useDingbats = FALSE)
}
theme_set(theme_gray(base_size = 11, base_family = font_family))

export_plot <- function(plot, stem, w, h, folder = "panels") {
  probe <- tempfile(fileext = ".pdf")
  open_pdf(probe, w, h)
  g <- if (inherits(plot, "patchwork"))
    patchworkGrob(plot) else if (inherits(plot, "ggplot"))
    ggplotGrob(plot) else plot
  # Free-scale facets have different feature-label widths. Give repeated left-axis lanes the same
  # width within each plot so plot-area gutters stay uniform.
  balance_axes <- function(gt) {
    if (!inherits(gt, "gtable"))
      return(gt)
    for (i in seq_along(gt$grobs)) if (inherits(gt$grobs[[i]], "gtable"))
      gt$grobs[[i]] <- balance_axes(gt$grobs[[i]])
    ax <- gt$layout[grepl("^axis-l-[0-9]+-[0-9]+$", gt$layout$name), , drop = FALSE]
    cols <- unique(ax$l)
    if (length(cols) > 1) {
      maxw <- do.call(unit.pmax, lapply(cols, function(i) gt$widths[i]))
      for (i in cols) gt$widths[i] <- maxw
    }
    gt
  }
  g <- balance_axes(g)
  dev.off()
  unlink(probe)
  out <- file.path(output_root, folder, stem)
  open_pdf(paste0(out, ".pdf"), w, h)
  grid.newpage()
  grid.draw(g)
  dev.off()
  svglite(paste0(out, ".svg"), width = w, height = h)
  grid.newpage()
  grid.draw(g)
  dev.off()
  agg_png(paste0(out, ".png"), width = w, height = h, units = "in", res = 300, background = "white")
  grid.newpage()
  grid.draw(g)
  dev.off()
  cat("EXPORTED", folder, stem, "\n")
  invisible(g)
}

replay_theme <- function(size = 11) {
  theme(text = element_text(family = font_family, size = size), plot.title = element_text(size = size +
    1), strip.text = element_text(size = size), axis.text = element_text(size = size - 1), legend.text = element_text(size = size -
    1), legend.title = element_text(size = size), plot.margin = margin(6, 7, 6, 7))
}

# Figure 3A-B: saved variance percentages and covariate statistics.
if (target %in% c("all", "figure3")) {
  d3a <- b$variance_data
  d3a$view <- factor(d3a$view, levels = c("proteomics", "lipidomics", "metabolomics"))
  d3a$group <- factor(ifelse(as.character(d3a$group) == "Cushings", "PPID", as.character(d3a$group)),
    levels = c("Control", "PPID"))
  d3a$label_y <- d3a$percentage + ifelse(d3a$percentage < 5, ifelse(d3a$view == "lipidomics", 2.3,
    1.1), 0.4)
  p3a <- ggplot(d3a, aes(factor, percentage, fill = view)) + geom_col(position = "dodge") + geom_text(aes(y = label_y,
    label = sprintf("%.1f%%", percentage)), position = position_dodge(0.9), vjust = -0.2, size = 2.8) +
    facet_wrap(~group) + scale_fill_manual(values = c(proteomics = "#239021", lipidomics = "#EB9A1B",
    metabolomics = "#1D90FF")) + labs(title = "Percentage Contribution of Variance by Factor and View",
    x = "Factor", y = "Percentage (%)", fill = "View") + theme_minimal(base_family = font_family) +
    replay_theme(11) + theme(legend.position = "bottom", panel.grid = element_blank())
  d3b <- b$covariate_data
  d3b$Status <- factor(ifelse(as.character(d3b$Status) == "Cushings", "PPID", as.character(d3b$Status)),
    levels = c("Control", "PPID"))
  d3b$Modality <- factor(as.character(d3b$Modality), levels = c("lipidomics", "Metabolomics", "Proteomics"),
    labels = c("Lipidomics", "Metabolomics", "Proteomics"))
  d3b$var <- factor(d3b$var, levels = c("Age", "Age_Batch", "Age_Batch_Sex", "Age_Sex", "Batch", "Batch_Sex",
    "Sex"))
  p3b <- ggplot(d3b, aes(factor, var, fill = logp)) + geom_tile(color = "white") + geom_text(aes(label = label),
    size = 2.8) + facet_grid(Status ~ Modality) + scale_x_continuous("Factor", breaks = 1:3, expand = c(0,
    0)) + scale_y_discrete(expand = c(0, 0)) + scale_fill_scico(palette = "devon", direction = -1,
    limits = c(0, 50), breaks = c(0, 50), labels = c("0", ">50"), name = expression(-log[10](P))) +
    coord_fixed() + labs(y = NULL) + theme_gray(base_family = font_family) + replay_theme(10) + theme(panel.grid = element_blank())
  export_plot(p3a, "Figure3A_variance", 8, 7)
  export_plot(p3b, "Figure3B_covariates", 7, 6)

  # Figure 3C: fixed archived top 12 loadings per view/factor (108 rows).
  top_w <- b$loading_top12
  top_w$view <- factor(top_w$view, levels = c("proteomics", "metabolomics", "lipidomics"), labels = c("Proteomics",
    "Metabolomics", "Lipidomics"))
  top_w$factor <- factor(top_w$factor, levels = c("Factor1", "Factor2", "Factor3"))
  top_w$Direction <- factor(ifelse(top_w$value < 0, "Negative", "Positive"), levels = c("Negative",
    "Positive"))
  top_w$panel <- interaction(top_w$view, top_w$factor, lex.order = TRUE)
  make_loading <- function(d, ncol = 3, size = 10) {
    ggplot(d, aes(x = reorder_within(feature, abs(value), panel), y = value, color = Direction)) +
      geom_hline(yintercept = 0, color = "grey88", linewidth = 0.4) + geom_linerange(aes(ymin = 0,
      ymax = value), linewidth = 0.35) + geom_point(size = 1.1) + scale_x_reordered() + coord_flip() +
      facet_wrap(~view + factor, scales = "free", ncol = ncol) + scale_color_manual(values = c(Negative = "#F8766D",
      Positive = "#619CFF"), drop = FALSE) + labs(x = NULL, y = "Factor loading", color = "Direction") +
      replay_theme(size) + theme(legend.position = "bottom")
  }
  p3c <- make_loading(top_w, 3, 10)
  export_plot(p3c, "Figure3C_loadings", 13, 12)
  for (f in levels(top_w$factor)) export_plot(make_loading(top_w[top_w$factor == f, ], 1, 11), paste0("Figure3C_",
    f), 6, 11)
  # Keep comparable facets at consistent sizes within each plotting role.
  fig3 <- (p3a + p3b + plot_layout(widths = c(1, 1)))/p3c + plot_layout(heights = c(1, 1.55)) + plot_annotation(tag_levels = "a")
  export_plot(fig3, "Figure3_replay", 14, 15, "figures")

}

# Figure 4A-C: fixed archived RF importance values.
if (target %in% c("all", "figure4")) {
  rfplots <- list()
  for (v in names(b$rf_importance)) {
    d <- b$rf_importance[[v]]
    d$Factor <- factor(d$Factor, levels = c("Factor1", "Factor2", "Factor3"))
    d$Group <- factor(d$Group, levels = c("negative", "positive"))
    d$FeatureKey <- reorder_within(d$Feature, -d$MeanDecreaseAccuracy, d$Factor)
    p <- ggplot(d, aes(x = FeatureKey, y = MeanDecreaseAccuracy, fill = Group)) + geom_col(position = "dodge",
      width = 0.8) + facet_wrap(~Factor, scales = "free") + coord_flip() + scale_x_reordered() +
      scale_y_continuous(breaks = scales::breaks_pretty(n = 3)) + scale_fill_manual(values = c(negative = "red",
      positive = "blue")) + labs(title = paste0("Feature Importance (", tools::toTitleCase(v),
      ")"), x = "Feature", y = "Mean decrease accuracy", fill = "Group") + theme_minimal(base_size = 11,
      base_family = font_family) + replay_theme(10) + theme(legend.position = "bottom")
    rfplots[[v]] <- p
    export_plot(p, paste0("Figure4_RF_", v), 10, 6)
  }

  # Figure 4D: saved summaries, with their original mean +/- SD error bars.
  auc <- b$auc_summary
  auc$Factor <- factor(auc$Factor, levels = c("Factor1", "Factor2", "Factor3", "Combined"))
  auc_colors <- c(Factor1 = "#FFFF00", Factor2 = "#CF8F78", Factor3 = "#A020F0", Combined = "#1F78B4")
  stopifnot(all(abs(auc$Upper - pmin(auc$Mean_AUC + auc$SD_AUC, 1)) < 1e-12), all(abs(auc$Lower - pmax(auc$Mean_AUC -
    auc$SD_AUC, 0)) < 1e-12))
  p4d <- ggplot(auc, aes(Factor, Mean_AUC, fill = Factor)) + geom_col(width = 0.6, alpha = 0.9, color = "black") +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.2, linewidth = 0.65) + geom_text(aes(label = sprintf("%.3f",
    Mean_AUC), y = Upper + 0.04), size = 3.5) + scale_fill_manual(values = auc_colors) + scale_y_continuous(limits = c(0,
    1.15), breaks = c(0, 0.3, 0.6, 0.9, 1.2)) + labs(title = "Mean AUC across stored validation runs",
    x = "Feature Group (Factor)", y = "Mean AUC (+/- SD)") + theme_classic(base_size = 12, base_family = font_family) +
    replay_theme(11) + theme(legend.position = "none", plot.title = element_text(hjust = 0.5))
  export_plot(p4d, "Figure4D_AUC", 7, 5)

  # Figure 4E-F: print the saved PCA coordinates and original scale mappings.
  d4ef <- b$pca_data
  d4ef$Status <- factor(ifelse(as.character(d4ef$Status) == "Cushings", "PPID", as.character(d4ef$Status)),
    levels = c("Control", "PPID"))
  d4ef$Type <- factor(ifelse(d4ef$Type == "Adjusted1", "Adjusted", as.character(d4ef$Type)), levels = c("Before",
    "Adjusted"))
  d4ef$Batch <- factor(d4ef$Batch, levels = c(1, 2, 3))
  p4ef <- ggplot(d4ef, aes(PC1, PC2, color = Batch, shape = Sex, size = Age, alpha = Status)) + geom_point() +
    scale_size(range = c(2, 6)) + scale_alpha_manual(values = c(Control = 0.5, PPID = 1)) + scale_color_manual(values = c(`1` = "#8CD3C7",
    `2` = "#FDFEB3", `3` = "#BEBADA")) + scale_shape_manual(values = c(Female = 16, Male = 17)) +
    facet_wrap(~Type) + labs(title = "PCA Comparison: Before vs Adjusted", color = "Batch", shape = "Sex",
    size = "Age", alpha = "Status") + theme_minimal(base_family = font_family) + replay_theme(11)
  export_plot(p4ef, "Figure4EF_PCA", 12, 6)
  fig4 <- (rfplots$proteomics + rfplots$lipidomics)/(rfplots$metabolomics + p4d)/p4ef + plot_layout(heights = c(1.1,
    1, 1.1)) + plot_annotation(tag_levels = "a")
  export_plot(fig4, "Figure4_replay", 14, 14, "figures")

}

# Figure 5 / raw heatmap companion: saved matrices, fixed row and column order.
if (target %in% c("all", "figure5", "supplement")) {
  annrow <- b$annotation_row
  annrow <- annrow[rownames(b$adjusted_matrix), , drop = FALSE]
  annrow <- annrow[, c("Polarity", "Factor", "Omics")]
  annrow$Omics <- factor(annrow$Omics, levels = c("Proteomics", "Lipidomics", "Metabolomics"))
  annrow$Factor <- factor(annrow$Factor, levels = c("Factor1", "Factor2", "Factor3"))
  annrow$Polarity <- factor(annrow$Polarity, levels = c("Positive", "Negative"))
  anncol <- b$annotation_col
  anncol$Batch <- factor(anncol$Batch, levels = c(1, 2, 3))
  ac <- list(Omics = c(Proteomics = "#239021", Lipidomics = "#EB9A1B", Metabolomics = "#1D90FF"), Factor = c(Factor1 = "#FFFF00",
    Factor2 = "#CF8F78", Factor3 = "#A020F0"), Polarity = c(Positive = "firebrick", Negative = "steelblue"),
    Sex = c(Male = "#1f77b4", Female = "#e377c2"), Age = colorRampPalette(c("#f7fcf5", "#00441b"))(100),
    Batch = c(`1` = "#8CD3C7", `2` = "#FDFEB3", `3` = "#BEBADA"))
  gaps <- which(vapply(seq_len(nrow(annrow) - 1), function(i) any(as.character(annrow[i, ]) != as.character(annrow[i +
    1, ])), logical(1)))
  heatmap_plot <- function(mat, title) {
    stopifnot(nrow(mat) == nrow(annrow), ncol(mat) == 15, all(is.finite(mat)))
    mat <- mat[rownames(b$adjusted_matrix), , drop = FALSE]
    stopifnot(identical(rownames(mat), rownames(annrow)), identical(colnames(mat), rownames(anncol)))
    pheatmap(mat, labels_row = sub("^PEF1[.]1$", "PEF1", rownames(mat)), annotation_row = annrow,
      annotation_col = anncol, annotation_colors = ac, gaps_row = gaps, cluster_rows = FALSE, cluster_cols = FALSE,
      color = colorRampPalette(c("blue", "white", "red"))(100), main = title, angle_col = 90, fontsize = 10,
      fontsize_row = 9, fontsize_col = 10, border_color = "grey60", silent = TRUE)$gtable
  }
  if (target %in% c("all", "figure5")) {
    ph5 <- heatmap_plot(b$adjusted_matrix, "Omics Heatmap - Top Features (Adjusted Data)")
    export_plot(ph5, "Figure5_adjusted_replay", 10, 16, "figures")
  }
  if (target %in% c("all", "supplement")) {
    phraw <- heatmap_plot(b$raw_matrix, "Omics Heatmap - Top Features (Raw Data)")
    export_plot(phraw, "Supplement_raw_heatmap", 10, 16, "supplement")
  }
}

# Additional archived model-comparison plots (no manuscript numbering assigned).
if (target %in% c("all", "supplement")) {
  supp_plots <- list()
  dh <- b$model_variance
  supp_plots$heatmap_plot <- ggplot(dh, aes(Factor, Model, fill = AverageVariance)) + geom_tile(color = "white") +
    facet_wrap(~Omics) + scale_fill_gradient(low = "white", high = "blue") + scale_x_discrete(labels = function(x) sub("Factor",
    "", x)) + labs(title = "Average Variance Explained", x = "Factor", y = "Model", fill = "Variance (%)") +
    theme_minimal(base_family = font_family) + replay_theme(11)
  dr <- b$model_redundancy
  supp_plots$redundant_factors_plot <- ggplot(dr, aes(Model, NumRedundant, fill = Model == "MOFA3")) +
    geom_col() + scale_fill_manual(values = c(`FALSE` = "grey70", `TRUE` = "#4DAF4A")) + labs(title = "Redundant Factors across Archived Models",
    x = "Model", y = "Redundant factors") + theme_minimal(base_family = font_family) + replay_theme(11) +
    theme(legend.position = "none")
  ds <- b$model_scores
  supp_plots$model_scores_plot <- ggplot(ds, aes(reorder(Model, -Score), Score, fill = Model == "MOFA3")) +
    geom_col() + scale_fill_manual(values = c(`FALSE` = "grey70", `TRUE` = "#4DAF4A")) + labs(title = "Archived Model Scores",
    x = "Model", y = "Score") + theme_minimal(base_family = font_family) + replay_theme(11) + theme(legend.position = "none")
  dt <- b$model_total_variance
  supp_plots$total_variance_plot <- ggplot(dt, aes(Model, TotalVariance, fill = Condition)) + geom_col(position = "dodge") +
    scale_fill_manual(values = c(Control = "#A6CEE3", PPID = "#FB9A99", Cushings = "#FB9A99")) +
    labs(title = "Total Explained Variance across Archived Models", x = "Model", y = "Summed variance contribution (%)") +
    theme_minimal(base_family = font_family) + replay_theme(11)
  for (n in names(supp_plots)) export_plot(supp_plots[[n]], paste0("Supplement_", n), if (n == "total_variance_plot")
    12 else 9, 6, "supplement")
}
writeLines(capture.output(sessionInfo()), file.path(output_root, "sessionInfo.txt"))
cat("Finished: ", target, "\nOutputs: ", normalizePath(output_root), "\n", sep = "")
