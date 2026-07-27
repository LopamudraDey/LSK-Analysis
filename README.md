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

outs
filtered_feature_bc_matrix.h5

atac_fragments.tsv.gz

The data directory is specified in the script as

data_dir <- "D:/LSK_DB/outs"


# scRNA-seq and scATAC-seq Analysis of LSK Cells with Ly6a (Sca-1) Expression Profiling

This script analyzes **10X Genomics Multiome RNA + ATAC sequencing data from LSK cells**, focusing on **Ly6a (Sca-1) expression-associated transcriptional and chromatin states**.

**Setup & QC:** Loads multiome data (RNA + ATAC), creates a Seurat object with integrated RNA and ATAC assays, and performs quality control filtering. RNA QC includes filtering based on gene counts, UMI counts, and mitochondrial content (**nFeature_RNA = 200–4000, nCount_RNA = 200–10000, mitochondrial content <10%**). ATAC quality metrics including fragment counts, nucleosome signal, and TSS enrichment are also evaluated and filtered.

**Core RNA Analysis:** Normalizes RNA expression data, identifies variable features, scales expression values, performs PCA, generates UMAP embeddings, and clusters cells using Seurat.

**Ly6a Expression Grouping:** Identifies Ly6a-positive cells and divides them into three expression groups based on tertiles: **Ly6a Low, Medium, and High** expression populations. Ly6a groups are visualized in RNA UMAP space.

**Cell Cycle Analysis:** Calculates cell-cycle scores and assigns cells into **G1, S, and G2/M phases**. Cell-cycle distributions are compared across Ly6a expression groups using UMAP, heatmaps, and proportion plots.

**Differential Expression Analysis:** Identifies Differentially Expressed Genes (DEGs) between **Ly6a High and Ly6a Low cells** using Seurat FindMarkers analysis.

**GSEA:** Performs Gene Set Enrichment Analysis (GSEA) using ranked DEGs and Gene Ontology Biological Process (**GO:BP**) gene sets. Generates enrichment plots and divergent bar plots showing positively and negatively enriched pathways.

**ORA:** Performs Gene Ontology Over-Representation Analysis (ORA) separately on Ly6a High and Ly6a Low enriched genes to identify significantly associated biological processes.

**Pearson Correlation:** Evaluates the relationship between **Ly6a expression and cell-cycle activity scores (S-phase and G2/M scores)** in Ly6a-positive cells using Pearson correlation analysis.

**ATAC Chromatin Analysis:** Processes chromatin accessibility data from Ly6a-positive cells using Signac, including TF-IDF normalization, LSI dimensional reduction, ATAC UMAP visualization, and comparison of chromatin landscapes across Ly6a expression groups.

**Differential Accessibility Analysis:** Identifies Ly6a-associated differential chromatin accessibility peaks, including High vs Low, Medium vs Low, and High vs Medium comparisons, with G1-phase-specific analysis to minimize cell-cycle effects.

**Motif Enrichment Analysis:** Performs transcription factor motif enrichment analysis on differential ATAC peaks using JASPAR motifs. Generates motif enrichment tables, motif preference comparisons, volcano plots, and sequence logos for enriched transcription factors.

**TF Expression Validation:** Validates candidate transcription factors at the RNA level by comparing expression across Ly6a Low, Medium, and High groups using statistical testing and visualization. 

