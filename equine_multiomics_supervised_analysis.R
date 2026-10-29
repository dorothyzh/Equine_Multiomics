setwd(dir = "~")

#packages
library(readr)
library(readxl)
library(dplyr)
library(tidyverse)
library(matrixStats) #rowSds
library(ggplot2)
library(missForest) #random forest imputation
library(doParallel) #parallel computing in random forest
library(IMIFA) #pareto scaling
library(factoextra) #pca
library(limma) #limma for differential expression analysis
library("RColorBrewer")
library(ComplexHeatmap) #heatmap
library(circlize) #colorRamp2 palette

#data import
#proteomics
prot <- read.csv("InSolution_Combined_Library.csv")
#metabolomics
met <- read.csv("data/Metabolites.csv")
#lipidomics
lipo <- read.csv("data/Lipids.csv")
#meta data
metadata <- read.csv("data/Metadata.csv")
#rename metadata
metadata$Status[metadata$Status == "Cushings"] <- "PPID"

#remove first row
prot <- as.data.frame(prot[-1,])
met <- as.data.frame(met[-1,])
lipo <- as.data.frame(lipo[-1,])

#store first column
feature_prot <- prot[,1]
feature_met <- met[,1]
feature_lipo <- lipo[,1]

#convert columns from character to numeric
prot <- as.data.frame(sapply(prot[,-1], as.numeric))
met <- as.data.frame(sapply(met[,-1], as.numeric))
lipo <- as.data.frame(sapply(lipo[,-1], as.numeric))

#feature as row names
rownames(prot) <- feature_prot
rownames(met) <- feature_met
rownames(lipo) <- feature_lipo

#colnames Cushings -> PPID
colnames(prot) <- gsub("Cushings", "PPID", colnames(prot))
colnames(met) <- gsub("Cushings", "PPID", colnames(met))
colnames(lipo) <- gsub("Cushings", "PPID", colnames(lipo))



#log2 transformation
prot <- log2(prot)
met <- log2(met)
lipo <- log2(lipo)

#convert -Inf to 0
met[met == -Inf] <- 0
lipo[lipo == -Inf] <- 0

#find features with at least 3 non-NA
prot <- prot[which(rowSums(is.na(prot)) <= ncol(prot)-3),]

#proteomics random forest imputation by group
registerDoParallel(cores = 6)
prot_rf_c <- missForest(prot[,1:8], maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
prot_rf_t <- missForest(prot[,9:15], maxiter = 10, ntree = 100, parallelize = "variable", verbose = TRUE)
prot_rf_c <- prot_rf_c$ximp
prot_rf_t <- prot_rf_t$ximp
prot_rf_ct <- cbind(prot_rf_c, prot_rf_t)

#metabolomics lipidomics imputation by 1/5 of the minimal value per feature
for (i in 1:nrow(met)) {
  if (length(which(met[i,] <= 0)) > 0) {
    print(paste(i, rownames(met)[i], colnames(met)[which(met[i,] <= 0)], #which sample is missing
                min(met[i, which(met[i,] > 0)])/5, sep = " ")) #what is the min of the feature
    met[i, met[i,] <= 0] <- min(met[i, which(met[i,] > 0)])/5
  } else {
    next
  }
}

#lipidomics
for (i in 1:nrow(lipo)) {
  if (length(which(lipo[i,] <= 0)) > 0) {
    print(paste(i, rownames(lipo)[i], colnames(lipo)[which(lipo[i,] <= 0)], #which sample is missing
                min(lipo[i, which(lipo[i,] > 0)])/5, sep = " ")) #what is the min of the feature
    lipo[i, lipo[i,] <= 0] <- min(lipo[i, which(lipo[i,] > 0)])/5
  } else {
    next
  }
}

met_pre <- met
lipo_pre <- lipo

#pareto scaling
met <- as.data.frame(pareto_scale(met))
lipo <- as.data.frame(pareto_scale(lipo))


#Figure 2 PCA + heatmap
#a ipodomics pca
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

#b metabolomics pca
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

# c proteomics pca
PCA_prot <- prcomp(t(prot_rf_ct))
fviz_pca_ind(PCA_prot,
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

#d combined omics pca
sum <- rbind(prot_rf_ct, met, lipo)
write.csv(sum, "multiomics_data.csv")
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

#e combined omics heatmap 
#differential expression analysis
#limma design matrix
design <- model.matrix(~ 0 + Status + Sex + Age + Batch, data = metadata)
#proteomics
fit.prot <- lmFit(prot_rf_ct, design)
con.prot <- makeContrasts(StatusPPID - StatusControl, levels = design)
fit.con.prot <- contrasts.fit(fit.prot, con.prot)
bayes.prot <- eBayes(fit.con.prot)
result.prot <- topTable(bayes.prot, sort.by="none", n=Inf, adjust = "fdr")
sig.prot <- result.prot[result.prot$adj.P.Val < 0.05 & abs(result.prot$logFC) > 1,]
write.csv(sig.prot, "sig_prot_insolution_comblibrary.csv")

#metabolomics
fit.met <- lmFit(met, design)
con.met <- makeContrasts(StatusPPID - StatusControl, levels = design)
fit.con.met <- contrasts.fit(fit.met, con.met)
bayes.met <- eBayes(fit.con.met)
result.met <- topTable(bayes.met, sort.by="none", n=Inf, adjust = "fdr")
sig.met <- result.met[result.met$adj.P.Val < 0.05 & abs(result.met$logFC) > 1,]
write.csv(sig.met, "sig_metabolomics.csv")

#lipo
fit.lipo <- lmFit(lipo, design)
con.lipo <- makeContrasts(StatusPPID - StatusControl, levels = design)
fit.con.lipo <- contrasts.fit(fit.lipo, con.lipo)
bayes.lipo <- eBayes(fit.con.lipo)
result.lipo <- topTable(bayes.lipo, sort.by="none", n=Inf, adjust = "fdr")
sig.lipo <- result.lipo[result.lipo$adj.P.Val < 0.05 & abs(result.lipo$logFC) > 1,]
write.csv(sig.lipo, "sig_lipodomics.csv")

#combined
fit.sum <- lmFit(sum, design)
con.sum <- makeContrasts(StatusPPID - StatusControl, levels = design)
fit.con.sum <- contrasts.fit(fit.sum, con.sum)
bayes.sum <- eBayes(fit.con.sum)
result.sum <- topTable(bayes.sum, sort.by="none", n=Inf, adjust = "fdr")
sig.sum <- result.sum[result.sum$adj.P.Val < 0.05 & abs(result.sum$logFC) > 1,]
write.csv(sig.sum, "sig_combined.csv")

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
prot_limma_all <- fitted(fit.prot)
prot_limma_ttest <- as.data.frame(prot_limma_all[rownames(prot_limma_all) %in% top_prot_ttest,])
colnames(prot_limma_ttest) <- colnames(prot)

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
colnames(met_limma_ttest) <- colnames(met)

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
colnames(lipo_limma_ttest) <- colnames(lipo)

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

#Figure S1
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

#Figure S3
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
lipo_class <- cbind(lipo_class, lipo_pre[4:625,])

#Figure S3a saturation effect
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

#Figure S3b length effect
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












