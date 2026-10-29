# Equine Multi-omics Biomarker Discovery

## 1. Introduction

**Background**\
Pituitary pars intermedia dysfunction (PPID) is a progressive neuroendocrine disorder that primarily affects aged horses, caused by the degeneration of hypothalamic dopaminergic neurons, leading to abnormal pituitary function. Diagnosing PPID is currently hampered by assay variability and the complexity of early-stage presentations.

**Results**\
This study presents an unsupervised method of integrating proteomics, lipidomics, and metabolomics to reveal cross-omic trends in PPID, overlooked by traditional analysis methods. Using multi-omics factor analysis, we find trends not seen in the conventional approach, such as divergence of C4BPA with triglycerides in the factor explaining the most lipidomic variance between the control and PPID-diagnosed horses, contradictory to current literary knowledge. Using a random forest classifier, we highlight potential biomarkers for PPID, which can be further explored as potential diagnostic and therapeutic targets.

**Conclusions**\
This study demonstrates an integrated and unsupervised multi-omics data analysis of equine PPID, highlighting the importance of complementary sample preparation procedures and statistical analysis methodologies to capture a broader image of the underlying biology. More extensive studies are required to verify these methods as diagnostic tools and validate potential biomarkers observed in the proteomic and lipidomic alterations. 


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

<!-- Describe the unsupervised analysis here. -->

| Figure | Description | Section of script |
|--------|-------------|-------------------|
| Figure X | <!-- description --> | <!-- e.g. line numbers or section header --> |
| Figure X | <!-- description --> | <!-- e.g. line numbers or section header --> |

The source code (R scripts) in this repository is released under the
[MIT License](LICENSE), an OSI-approved open source license.
The data and metadata are released under
[CC BY 4.0](data/LICENSE).
