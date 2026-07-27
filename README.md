# LSK-Analysis
<img width="1868" height="899" alt="image" src="https://github.com/user-attachments/assets/ad3020ef-e80e-4366-a187-e7e88038e844" />

# LSK DATASET
LSK scRNA .h5 file is available at: 
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM8960398
# Requirements
**Software**
R (version 4.2 or newer recommended)
RStudio (recommended)
Required Packages

The analysis uses both CRAN and Bioconductor packages.

CRAN Packages

Seurat

Signac

tidyverse

hdf5r

ggpubr

rstatix

stringr

UpSetR

**Bioconductor Packages**

clusterProfiler

org.Mm.eg.db

JASPAR2020

TFBSTools

BSgenome.Mmusculus.UCSC.mm10

# Input Files

The workflow expects a standard 10x Genomics Multiome output directory.

outs/
├── filtered_feature_bc_matrix.h5
└── atac_fragments.tsv.gz

The data directory is specified in the script as

data_dir <- "D:/LSK_DB/outs"

# LSK scRNA-seq Analysis Summary
This script analyzes scRNA-seq data from LSK cells, focusing on Ly6a (Sca-1) expression.

Setup & QC: Loads data, then creates a Seurat object, and performs QC filtering (nFeature, nCount, %mt). Filters cells based on established quality metrics:nFeature_RNA= 200 to 4000, nCount_RNA= 200 to  10000 and Max mito content: < 10%.

Core Analysis: Normalizes, scales, runs PCA, UMAP, and clustering.

Ly6a Grouping: Subsets Ly6a-positive cells and divides them into Low, Medium, and High expression groups (tertiles).

Cell Cycle: Scores cell cycle (G1, S, G2/M) and plots its distribution across the Ly6a groups.

GSEA: Finds Differentially Expressed Genes (DEGs) between Ly6a High and Low cells. 
Runs Gene Set Enrichment Analysis (GSEA) on the ranked DEGs (using GO:BP) and generates a divergent bar plot of the top enriched pathways.

ORA: Also performs standard GO Over-Representation Analysis (ORA) on significant up/downregulated genes.

Pearson Correlation: Association of Ly6a expression with cell cycle activity in Ly6a-positive cells. 

