setwd("~/Cushings")

#package
library(readr)
library(readxl)
library(missForest)
library(doParallel)
library(factoextra) #pca
library(tibble)
library(limma)
library(ggplot2)
library(VennDiagram)
library("RColorBrewer")
library(ComplexHeatmap)
library(IMIFA)
library(dplyr)
library(tidyverse)
library(ggrepel)
library(matrixStats)

#load data
metadata <- read_excel("Horse_Metadata.xlsx")
Combined_with_Proteins <- read_csv("Combined_with_Proteins.csv")
prot_insolution_comblibrary <- read.csv("InSolution_Combined_Library.csv")
Total_Proteins_BarChart <- read_excel("Total_Proteins_BarChart.xlsx")

#rename metadata
metadata$Status[metadata$Status == "Cushings"] <- "PPID"


#proteomics data reshape
prot_insolution_comblibrary <- as.data.frame(prot_insolution_comblibrary[-1,])
prot_insolution_comblibrary[,-1] <- sapply(prot_insolution_comblibrary[,-1], as.numeric)
rownames(prot_insolution_comblibrary) <- prot_insolution_comblibrary[,1]
prot_insolution_comblibrary <-  prot_insolution_comblibrary[,-1]

#lipidomics and metabolomics data
lipo <- Combined_with_Proteins[875:2493,]
met <- Combined_with_Proteins[2494:6400,]

#log2 transform
prot_insolution_comblibrary <- log2(prot_insolution_comblibrary)
met <- log2(met)
lipo <- log2(lipo)

#replace INF by 0
met[met == -Inf] <- 0
lipo[lipo == -Inf] <- 0

#find features with at least 3 non-NA in proteomics
prot_insolution_comblibrary <- prot_insolution_comblibrary[which(rowSums(is.na(prot_insolution_comblibrary)) <= ncol(prot_insolution_comblibrary)-3),]


#imputation
#proteomics
#random forest imputation by group
registerDoParallel(cores = 6)
prot_insolution_comblibrary_rf_c <- missForest(prot_insolution_comblibrary[,1:8], maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
prot_insolution_comblibrary_rf_t <- missForest(prot_insolution_comblibrary[,9:15], maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
prot_insolution_comblibrary_rf_c <- prot_insolution_comblibrary_rf_c$ximp
prot_insolution_comblibrary_rf_t <- prot_insolution_comblibrary_rf_t$ximp
prot_insolution_comblibrary_rf_ct <- cbind(prot_insolution_comblibrary_rf_c, prot_insolution_comblibrary_rf_t)
write.csv(prot_insolution_comblibrary_rf_ct, "cushings_proteomics_comblibrary_imputed.csv")

#metabolomics lipidomomics imputation
#impute 
for (i in 1:nrow(met)) {
  if (length(which(met[i,] <= 0)) > 0) {
    print(paste(i, rownames(met)[i], colnames(met)[which(met[i,] <= 0)], min(met[i, which(met[i,] > 0)])/5, sep = " "))
    met[i, met[i,] <= 0] <- min(met[i, which(met[i,] > 0)])/5
  } else {
    next
  }
}

for (i in 1:nrow(lipo)) {
  if (length(which(lipo[i,] <= 0)) > 0) {
    print(paste(i, rownames(lipo)[i], colnames(lipo)[which(lipo[i,] <= 0)], min(lipo[i, which(lipo[i,] > 0)])/5, sep = " "))
    lipo[i, lipo[i,] <= 0] <- min(lipo[i, which(lipo[i,] > 0)])/5
  } else {
    next
  }
}

#pareto scaling
met <- pareto_scale(met)
lipo <- pareto_scale(lipo)

write.csv(met, "cushings_metabolomics_imputed.csv")
write.csv(lipo, "cushings_lipodomics_imputed.csv")


#pca
#proteomics
PCA_insolution_comblibrary_rf <- prcomp(t(prot_insolution_comblibrary_rf_ct))
fviz_pca_ind(PCA_insolution_comblibrary_rf,
             geom = "point", #"point" or "text"
             habillage = metadata$Status,
             repel = TRUE,
             addEllipses = TRUE,
             ellipse.level = 0.95,
             palette  = c("Control" = "orange","PPID" = "purple")) +
  labs(title = "Proteomics") +
  theme(plot.title = element_text(size = 25, hjust = 0.5),
        axis.title = element_text(size = 15))
ggsave("pca_prot_final.svg",
       plot = last_plot(),
       path = "figures/",
       device = "svg",
       width = 1500,
       height = 1500,
       units = "px",
       dpi = 300)

#metabolomics
PCA_met <- prcomp(t(met))
fviz_pca_ind(PCA_met,
             geom = "point", #"point" or "text"
             habillage = metadata$Status,
             #repel = TRUE,
             addEllipses = TRUE,
             ellipse.level = 0.95,
             palette  = c("Control" = "orange","PPID" = "purple")) +
  labs(title = "Metabolomics") +
  theme(plot.title = element_text(size = 25, hjust = 0.5),
        axis.title = element_text(size = 15))
ggsave("pca_met_final.svg",
       plot = last_plot(),
       path = "figures/",
       device = "svg",
       width = 1500,
       height = 1500,
       units = "px",
       dpi = 300)

#lipodomics
PCA_lipo <- prcomp(t(lipo))
fviz_pca_ind(PCA_lipo,
             geom = "point", #"point" or "text"
             habillage = metadata$Status,
             repel = TRUE,
             addEllipses = TRUE,
             ellipse.level = 0.95,
             palette  = c("Control" = "orange","PPID" = "purple")) +
  labs(title = "lipidomics") +
  theme(plot.title = element_text(size = 25, hjust = 0.5),
        axis.title = element_text(size = 15))
ggsave("pca_lipo_final.svg",
       plot = last_plot(),
       path = "figures/",
       device = "svg",
       width = 1500,
       height = 1500,
       units = "px",
       dpi = 300)

#combined omics
sum <- rbind(prot_insolution_comblibrary_rf_ct, met, lipo)
PCA_sum <- prcomp(t(sum))
fviz_pca_ind(PCA_sum,
             geom = "point", #"point" or "text"
             habillage = metadata$Status,
             repel = TRUE,
             addEllipses = TRUE,
             ellipse.level = 0.95,
             palette  = c("Control" = "orange","PPID" = "purple")) +
  labs(title = "Combined") +
  theme(plot.title = element_text(size = 25, hjust = 0.5),
        axis.title = element_text(size = 15))
ggsave("pca_combine_final.svg",
       plot = last_plot(),
       path = "figures/",
       device = "svg",
       width = 1500,
       height = 1500,
       units = "px",
       dpi = 300)
write.csv(sum, "multiomics_data.csv")


#differential expression analysis
#limma
design <- model.matrix(~ 0 + Status + Sex + Age + Batch, data = metadata)
#proteomics
fit.prot_insolution_comblibrary <- lmFit(prot_insolution_comblibrary_rf_ct, design)
con.prot_insolution_comblibrary <- makeContrasts(StatusCushings - StatusControl, levels = design)
fit.con.prot_insolution_comblibrary <- contrasts.fit(fit.prot_insolution_comblibrary, con.prot_insolution_comblibrary)
bayes.prot_insolution_comblibrary <- eBayes(fit.con.prot_insolution_comblibrary)
result.prot_insolution_comblibrary <- topTable(bayes.prot_insolution_comblibrary, sort.by="none", n=Inf, adjust = "fdr")
sig.prot_insolution_comblibrary <- result.prot_insolution_comblibrary[result.prot_insolution_comblibrary$adj.P.Val < 0.05 & abs(result.prot_insolution_comblibrary$logFC) > 1,]
write.csv(result.prot_insolution_comblibrary, "result_prot_insolution_comblibrary.csv")
write.csv(sig.prot_insolution_comblibrary, "sig_prot_insolution_comblibrary.csv")

#metabolomics
fit.met <- lmFit(met, design)
con.met <- makeContrasts(StatusCushings - StatusControl, levels = design)
fit.con.met <- contrasts.fit(fit.met, con.met)
bayes.met <- eBayes(fit.con.met)
result.met <- topTable(bayes.met, sort.by="none", n=Inf, adjust = "fdr")
sig.met <- result.met[result.met$adj.P.Val < 0.05 & abs(result.met$logFC) > 1,]
write.csv(result.met, "result_metabolomics.csv")
write.csv(sig.met, "sig_metabolomics.csv")

#lipo
fit.lipo <- lmFit(lipo, design)
con.lipo <- makeContrasts(StatusCushings - StatusControl, levels = design)
fit.con.lipo <- contrasts.fit(fit.lipo, con.lipo)
bayes.lipo <- eBayes(fit.con.lipo)
result.lipo <- topTable(bayes.lipo, sort.by="none", n=Inf, adjust = "fdr")
sig.lipo <- result.lipo[result.lipo$adj.P.Val < 0.05 & abs(result.lipo$logFC) > 1,]
write.csv(result.lipo, "result_lipodomics.csv")
write.csv(sig.lipo, "sig_lipodomics.csv")

#combined
fit.sum <- lmFit(sum, design)
con.sum <- makeContrasts(StatusCushings - StatusControl, levels = design)
fit.con.sum <- contrasts.fit(fit.sum, con.sum)
bayes.sum <- eBayes(fit.con.sum)
result.sum <- topTable(bayes.sum, sort.by="none", n=Inf, adjust = "fdr")
sig.sum <- result.sum[result.sum$adj.P.Val < 0.05 & abs(result.sum$logFC) > 1,]
write.csv(result.lipo, "result_combined.csv")
write.csv(sig.sum, "sig_combined.csv")


#heatmap
# Define colors for each levels of qualitative variables
col = list(Condition = c("Control" = "orange","PPID" = "purple"),
           Batch = c("1" = "#8CD3C7", "2" = "#FDFEB3", "3" = "#BEBADA"),
           Age = colorRamp2(c(0, 32), c("#FFFFFF", "#00441B")))
# Create the heatmap annotation
ha <- HeatmapAnnotation(Condition = metadata$Status,
                        Batch = metadata$Batch,
                        Age = metadata$Age,
                        col = col)

#proteomics
#top features from metaboanalyst t-test
top_prot_ttest <- c("HSD17B4", "PHB2","CCT7","RPS16","DDR1",
                    "SACM1L","ENO3","THBS4","BCAP31","COL15A1",
                    "COMT","IL6ST","NPEPPS","TBCB","HSPA4",
                    "CDON","CTSA","ARPC5","COMP","MRC2")
#get limma regressed data
prot_limma_all <- fitted(fit.prot_insolution_comblibrary)
prot_limma_ttest <- as.data.frame(prot_limma_all[rownames(prot_limma_all) %in% top_prot_ttest,])
colnames(prot_limma_ttest) <- colnames(met_limma)

prot_limma_ttest$feature <- rownames(prot_limma_ttest)
prot_limma_ttest <- prot_limma_ttest %>% 
  pivot_longer(cols = 1:(ncol(.)-1)) %>%
  group_by(feature) %>%
  mutate(value = as.numeric(scale(value))) %>% #scale within each feature
  ungroup() %>% 
  pivot_wider(names_from = name, values_from = value) %>%
  as.data.frame()
rownames(prot_limma_ttest) <- prot_limma_ttest$feature
prot_limma_ttest <- prot_limma_ttest[,-1]

svg(filename = "figures/heatmap_proteomics_relativeintensity_ttest.svg",
    width = 12, height = 6.67, #unit = inches
    family = "arial")
Heatmap(prot_limma_ttest,
        name = paste("Relative", "Intensity", sep = "\n"),
        column_title = "Sample",
        row_title = "Proteomics",
        clustering_method_columns = "ward.D2",
        top_annotation = ha,
        width = unit(20, "cm"), 
        height = unit(10, "cm"),
        cluster_columns = FALSE)
dev.off()

#metabolites
#top features from metaboanalyst t-test
top_met_ttest <- c("5.16_173.0599m/z", "1.23_858.1001m/z", "1.29_311.9875m/z","11.71_327.2898m/z",
                   "11.71_543.3143m/z","11.92_262.6368m/z","1.40_188.0347m/z","1.83_651.0505m/z",
                   "11.74_498.3432n","11.71_259.1554m/z","11.71_541.3011m/z","1.13_292.0785m/z",
                   "1.73_821.0088m/z","1.76_776.9483m/z","L-allo-Threonine","10.78_398.2539m/z",
                   "1.42_164.0386m/z","11.62_424.3419m/z","11.71_605.2856m/z","1.13_276.1045m/z")
#get limma regressed data
met_limma_all <- fitted(fit.met)
met_limma_ttest <- as.data.frame(met_limma_all[rownames(met_limma_all) %in% top_met_ttest,])
colnames(met_limma_ttest) <- colnames(met_limma)

met_limma_ttest$feature <- rownames(met_limma_ttest)
met_limma_ttest <- met_limma_ttest %>% 
  pivot_longer(cols = 1:(ncol(.)-1)) %>%
  group_by(feature) %>%
  mutate(value = as.numeric(scale(value))) %>% #scale within each feature
  ungroup() %>% 
  pivot_wider(names_from = name, values_from = value) %>%
  as.data.frame()
rownames(met_limma_ttest) <- met_limma_ttest$feature
met_limma_ttest <- met_limma_ttest[,-1]

svg(filename = "figures/heatmap_metabolomics_relativeintensity_ttest.svg",
    width = 12, height = 6.67, #unit = inches
    family = "arial")
Heatmap(met_limma_ttest,
        name = paste("Relative", "Intensity", sep = "\n"),
        column_title = "Sample",
        row_title = "Metabolomics",
        clustering_method_columns = "ward.D2",
        top_annotation = ha,
        width = unit(20, "cm"), 
        height = unit(10, "cm"),
        cluster_columns = FALSE)
dev.off()

#lipidomics
top_lipo_ttest <- c("11.845_673.5946","11.833_669.6359","TG O-30:0","11.831_770.7064",
                    "11.837_739.7128","CE 21:1","LPC 14:0","TG 48:3",
                    "TG O-51:9","TG O-42:6","11.831_371.3664","DG 48:9",
                    "11.843_1058.913","DG 51:8","3.511_724.4121",
                    "TG 50:5","7.664_1129.62341","11.832_281.2454","Cer 43:2;O2","1.489_511.2237")
#get limma regressed data
lipo_limma_all <- fitted(fit.lipo)
lipo_limma_ttest <- as.data.frame(lipo_limma_all[rownames(lipo_limma_all) %in% top_lipo_ttest,])
colnames(lipo_limma_ttest) <- colnames(met_limma)

lipo_limma_ttest$feature <- rownames(lipo_limma_ttest)
lipo_limma_ttest <- lipo_limma_ttest %>% 
  pivot_longer(cols = 1:(ncol(.)-1)) %>%
  group_by(feature) %>%
  mutate(value = as.numeric(scale(value))) %>% #scale within each feature
  ungroup() %>% 
  pivot_wider(names_from = name, values_from = value) %>%
  as.data.frame()
rownames(lipo_limma_ttest) <- lipo_limma_ttest$feature
lipo_limma_ttest <- lipo_limma_ttest[,-1]

svg(filename = "figures/heatmap_lipodomics_relativeintensity_ttest.svg",
    width = 12, height = 6.67, #unit = inches
    family = "arial")
Heatmap(lipo_limma_ttest,
        name = paste("Relative", "Intensity", sep = "\n"),
        column_title = "Sample",
        row_title = "Lipodomics",
        clustering_method_rows = "ward.D2",
        top_annotation = ha,
        width = unit(20, "cm"), 
        height = unit(10, "cm"),
        cluster_columns = FALSE)
dev.off()


#volcano plot
#limma result of omics combined analysis
result.sum$category <- NA
result.sum$category[rownames(result.sum) %in% rownames(prot) & result.sum$adj.P.Val < 0.05 & abs(result.sum$logFC) > 1] <- "protein"
result.sum$category[rownames(result.sum) %in% rownames(met) & result.sum$adj.P.Val < 0.05 & abs(result.sum$logFC) > 1] <- "metabolite"
result.sum$category[rownames(result.sum) %in% rownames(lipo) & result.sum$adj.P.Val < 0.05 & abs(result.sum$logFC) > 1] <- "lipid"
result.sum$category[is.na(result.sum$category)] <- "ns"
result.sum$category <- factor(result.sum$category, levels = c("protein", "lipid", "metabolite", "ns"))

ggplot(result.sum, aes(x = logFC, y = -log10(adj.P.Val), color = category)) + 
  geom_point(alpha = 0.5) + 
  scale_color_manual(values = c("lipid" = "#EB9A1B", "metabolite" = "#1D90FF", "protein" = "#239021", "ns" = "grey"),
                     labels = c("lipid" = "Lipid", "metabolite" = "Metabolite", "protein" = "Protein", "ns" = "Not significant")) +
  theme_bw() +
  labs(x = bquote(paste("log"["2"]*"(Fold Change), PPID vs. Control")),
       y = bquote(paste("-log"["10"]*"(adjusted p-value)")),
       color = "Category")
ggsave("volcano_sum.svg",
       plot = last_plot(),
       path = "figures/",
       device = "svg",
       width = 3000,
       height = 2000,
       units = "px",
       dpi = 300)


#lipid class 
lipo_class <- data.frame(lipid = rownames(lipo)[4:625],
                         class = NA,
                         length = NA,
                         saturation = NA)

lipo_class$class <- str_split_i(lipo_class$lipid, " ", 1)
#get lipid length and saturation from lipid name
for (i in 1:length(lipo_class$lipid)) {
  lipo_class$length[i] <- substr(lipo_class$lipid[i], 
                                 str_locate(lipo_class$lipid[i], ":")[1]-2,
                                 str_locate(lipo_class$lipid[i], ":")[1]-1)
  lipo_class$saturation[i] <- substr(lipo_class$lipid[i], 
                                     str_locate(lipo_class$lipid[i], ":")[1]+1,
                                     str_locate(lipo_class$lipid[i], ":")[1]+1)
}
#combine lipid class with lipodomics abundance
lipo_class <- cbind(lipo_class, lipo[4:625,])

#length effect
lipo_length_sum <- data.frame(length = unique(lipo_class$length),
                              mean_cont = NA,
                              mean_trea = NA,
                              se_cont = NA,
                              se_trea = NA)

#get average intensity for features in each length by conditions
for (i in 1:length(lipo_length_sum$length)) {
  if (nrow(lipo_class[lipo_class$length == lipo_length_sum$length[i],]) == 1) {
    lipo_length_sum$mean_cont[i] <- mean(colMeans(lipo_class[lipo_class$length == lipo_length_sum$length[i],5:12]))
    lipo_length_sum$mean_trea[i] <- mean(colMeans(lipo_class[lipo_class$length == lipo_length_sum$length[i],13:19]))
    lipo_length_sum$se_cont[i] <- 0
    lipo_length_sum$se_trea[i] <- 0
  } else if (nrow(lipo_class[lipo_class$length == lipo_length_sum$length[i],]) > 1) {
    lipo_length_sum$mean_cont[i] <- mean(colMeans(lipo_class[lipo_class$length == lipo_length_sum$length[i],5:12]))
    lipo_length_sum$mean_trea[i] <- mean(colMeans(lipo_class[lipo_class$length == lipo_length_sum$length[i],13:19]))
    lipo_length_sum$se_cont[i] <- sd(colMeans(lipo_class[lipo_class$length == lipo_length_sum$length[i],5:12]))/sqrt(7)
    lipo_length_sum$se_trea[i] <- sd(colMeans(lipo_class[lipo_class$length == lipo_length_sum$length[i],5:12]))/sqrt(6)
  }
}

lipo_length_sum <- lipo_length_sum %>%
  pivot_longer(cols = 2:5,
               names_to = c(".value", "condition"),
               names_pattern = "(.*)_(cont|trea)") 
#barplot
ggplot(data = lipo_length_sum, aes(x = length, y = mean, fill = condition)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean-1.96*se, ymax = mean+1.96*se), width=.2,
                position=position_dodge(.9)) +
  scale_fill_manual(values = c("cont" = "orange", "trea" = "purple"),
                    labels = c("cont" = "Control", "trea" = "PPID")) +
  labs(x = "Lipid Chain Length",
       y = bquote(paste("Mean log"["2"]*"Intensity")),
       fill = "Condition") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12))
ggsave("barplot_lipidlength_intensity.svg",
       plot = last_plot(),
       path = "figures/",
       device = "svg",
       width = 3300,
       height = 1496,
       units = "px",
       dpi = 300)

#saturation effect
lipo_saturation_sum <- data.frame(saturation = unique(lipo_class$saturation),
                                  mean_cont = NA,
                                  mean_trea = NA,
                                  se_cont = NA,
                                  se_trea = NA)

#get average intensity for features in each saturation by conditions
for (i in 1:length(lipo_saturation_sum$saturation)) {
  if (nrow(lipo_class[lipo_class$saturation == lipo_saturation_sum$saturation[i],]) == 1) {
    lipo_saturation_sum$mean_cont[i] <- mean(colMeans(lipo_class[lipo_class$saturation == lipo_saturation_sum$saturation[i],5:12]))
    lipo_saturation_sum$mean_trea[i] <- mean(colMeans(lipo_class[lipo_class$saturation == lipo_saturation_sum$saturation[i],13:19]))
    lipo_saturation_sum$se_cont[i] <- 0
    lipo_saturation_sum$se_trea[i] <- 0
  } else if (nrow(lipo_class[lipo_class$saturation == lipo_saturation_sum$saturation[i],]) > 1) {
    lipo_saturation_sum$mean_cont[i] <- mean(colMeans(lipo_class[lipo_class$saturation == lipo_saturation_sum$saturation[i],5:12]))
    lipo_saturation_sum$mean_trea[i] <- mean(colMeans(lipo_class[lipo_class$saturation == lipo_saturation_sum$saturation[i],13:19]))
    lipo_saturation_sum$se_cont[i] <- sd(colMeans(lipo_class[lipo_class$saturation == lipo_saturation_sum$saturation[i],5:12]))/sqrt(7)
    lipo_saturation_sum$se_trea[i] <- sd(colMeans(lipo_class[lipo_class$saturation == lipo_saturation_sum$saturation[i],5:12]))/sqrt(6)
  }
}

lipo_saturation_sum <- lipo_saturation_sum %>%
  pivot_longer(cols = 2:5,
               names_to = c(".value", "condition"),
               names_pattern = "(.*)_(cont|trea)") 

#barplot
ggplot(data = lipo_saturation_sum, aes(x = saturation, y = mean, fill = condition)) +
  geom_bar(stat = "identity", position = position_dodge()) +
  geom_errorbar(aes(ymin = mean-1.96*se, ymax = mean+1.96*se), width=.2,
                position=position_dodge(.9)) +
  scale_fill_manual(values = c("cont" = "orange", "trea" = "purple"),
                    labels = c("cont" = "Control", "trea" = "PPID")) +
  labs(x = "Lipid Unsaturation Number",
       y = bquote(paste("Mean log"["2"]*"Intensity")),
       fill = "Condition") +
  theme_bw() +
  theme(axis.text.x = element_text(size = 10),
        axis.text.y = element_text(size = 10),
        axis.title = element_text(size = 12),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 12),
        legend.position = "none")
ggsave("barplot_lipidsaturation_intensity.svg",
       plot = last_plot(),
       path = "figures/",
       device = "svg",
       width = 1496,
       height = 1496,
       units = "px",
       dpi = 300)
