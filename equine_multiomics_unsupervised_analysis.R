# Load required libraries
library(MOFA2)

# Load data (assumed to be preprocessed and in list form)
data <- list(
  proteomics = read.csv("proteomics_data.csv"),
  metabolomics = read.csv("metabolomics_data.csv"),
  lipidomics = read.csv("lipidomics_data.csv")
)

# Metadata
meta <- read.csv("metadata.csv")

# Normalize and prepare data for MOFA
data_norm <- lapply(data, function(x) normalize.quantiles(as.matrix(x)))
data_list <- lapply(data_norm, function(x) list(data = x, meta = meta))

# Create MOFA object
mofa_obj <- create_mofa(data_list)

# Set options
model_options <- get_default_model_options(mofa_obj)
model_options$num_factors <- 3

training_options <- get_default_training_options(mofa_obj)

# Train the model
mofa_obj <- prepare_mofa(mofa_obj, model_options, training_options)
model <- run_mofa(mofa_obj, "mofa_model.hdf5", use_basilisk = TRUE)


###

# Load necessary libraries
library(ggplot2)

# Define color mappings
batch_colors <- c("#8CD3C7", "#FDFEB3", "#BEBADA")
age_colors <- colorRampPalette(c("#FFFFFF", "#00441B"))(100)
omic_colors <- c("Lipidomics" = "#EB9A1B", "Metabolomics" = "#1D90FF", "Proteomics" = "#239021")
status_colors <- c("Control" = "#D91E36", "PPID" = "#372B8D")

# Set plotting order
p[[1]]$data$view <- factor(p[[1]]$data$view, levels = c("proteomics", "lipidomics", "metabolomics"))

# Replace group names
p[[1]]$data$group <- factor(p[[1]]$data$group, levels = c("Control", "Cushings"), labels = c("Control", "PPID"))

# Custom color mapping
custom_colors <- c("proteomics" = "#239021", "lipidomics" = "#EB9A1B", "metabolomics" = "#1D90FF")

# Plotting
plotgg <- ggplot(p[[1]]$data, aes(x = factor, y = percentage, fill = view)) +
  geom_bar(stat = "identity", position = "dodge") +
  geom_text(aes(label = sprintf("%.1f%%", percentage)), position = position_dodge(0.9), vjust = -0.5, size = 3) +
  facet_wrap(~ group) +
  labs(title = "Percentage Contribution of Variance by Factor and View", x = "Factor", y = "Percentage (%)") +
  scale_fill_manual(values = custom_colors) +
  theme_minimal() +
  theme(legend.position = "bottom", text = element_text(size = 12))

# Save as SVG
ggsave("factor_all_variance_plot.svg", plot = plotgg, width = 8, height = 8, dpi = 300)


library(ggplot2)

# Transpose data for PCA calculation
pp_t <- t(pp)
pp_adjusted_t <- t(pp_adjusted)

# Perform PCA
pca_before <- prcomp(pp_t, scale. = TRUE)
pca_after1 <- prcomp(pp_adjusted_t, scale. = TRUE)

# Combine PCA data with metadata
pca_before_df <- data.frame(pca_before$x[, 1:2], meta, Type = "Before")
pca_after1_df <- data.frame(pca_after1$x[, 1:2], meta, Type = "Adjusted1")
pca_all_df <- rbind(pca_before_df, pca_after1_df)

# Force Batch to be a factor with ordered levels
pca_all_df$Batch <- factor(pca_all_df$Batch, levels = c(1, 2, 3))

# Custom batch colors
batch_colors <- c("1" = "#8CD3C7", "2" = "#FDFEB3", "3" = "#BEBADA")

# Plot PCA comparison
p_all <- ggplot(pca_all_df, aes(x = PC1, y = PC2, color = Batch, shape = as.factor(Sex), size = Age, alpha = as.factor(Status))) +
  geom_point() +
  scale_size(range = c(2, 6)) +
  scale_alpha_manual(values = c(0.5, 1)) +
  scale_color_manual(values = batch_colors) +
  facet_wrap(~Type) +
  labs(title = "PCA Comparison: Before vs Adjusted1", color = "Batch", shape = "Sex", size = "Age", alpha = "Status") +
  theme_minimal()

# Show plot
print(p_all)

# Save as SVG
ggsave("PCA_Before_vs_Adjusted1.svg", plot = p_all, width = 12, height = 6, dpi = 300)


library(pheatmap)
library(RColorBrewer)
library(svglite)

# Define color mappings for heatmap
omics_colors <- c("lipidomics" = "#EB9A1B", "metabolomics" = "#1D90FF", "proteomics" = "#239021")
polarity_colors <- c("positive" = "firebrick", "negative" = "steelblue")
factor_levels <- unique(row_annotation$Factor)
factor_colors <- setNames(colorRampPalette(c("yellow", "purple"))(length(factor_levels)), factor_levels)
sex_colors <- c("Male" = "#1f77b4", "Female" = "#e377c2")
age_unique <- sort(unique(meta$Age))
age_colors <- setNames(colorRampPalette(c("#f7fcf5", "#00441b"))(length(age_unique)), as.character(age_unique))
batch_colors <- c("1" = "#8CD3C7", "2" = "#FDFEB3", "3" = "#BEBADA")

# Row annotation data
annotation_row_colors <- data.frame(Polarity = annotation_row00$Polarity, Factor = annotation_row00$Factor, Omics = annotation_row00$Omics)
rownames(annotation_row_colors) <- rownames(sample_data_scaled_matrix2)

# Column annotation data
annotation_col <- meta[match(colnames(sample_data_scaled_matrix2), meta$sample), c("Sex", "Age", "Batch")]
rownames(annotation_col) <- colnames(sample_data_scaled_matrix2)

# Force Batch to be a factor
annotation_col$Batch <- factor(annotation_col$Batch, levels = c(1, 2, 3))

# Combine all color mappings
annotation_colors <- list(
  Omics = omics_colors,
  Polarity = polarity_colors,
  Factor = factor_colors,
  Sex = sex_colors,
  Age = age_colors,
  Batch = batch_colors
)

# Calculate row gaps for split annotations
omics_positions <- which(diff(as.numeric(factor(annotation_row$Omics))) != 0)
factor_positions <- which(diff(as.numeric(factor(annotation_row$Factor))) != 0)
polarity_positions <- which(diff(as.numeric(factor(annotation_row$Polarity))) != 0)
gaps_row <- sort(unique(c(omics_positions, factor_positions, polarity_positions)))
sample_data_scaled_matrix <- as.matrix(sample_data_scaled)

# Start creating SVG heatmap
svg("omics_heatmap_with_col_annotations.svg", width = 10, height = 16)

pheatmap(
  sample_data_scaled_matrix2, 
  labels_row = annotation_row$Gene_ID,  
  annotation_row = annotation_row_colors,  
  annotation_col = annotation_col,  
  annotation_colors = annotation_colors,    
  gaps_row = gaps_row,  
  cluster_rows = FALSE, 
  cluster_cols = FALSE,  
  color = colorRampPalette(c("blue", "white", "red"))(100),
  main = "Omics Heatmap - Top Features with Sample Annotations",
  angle_col = 90
)

dev.off()
