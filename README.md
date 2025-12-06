# LSK-Analysis

LSK scRNA-seq Analysis SummaryThis script analyzes scRNA-seq data from LSK cells, focusing on Ly6a (Sca-1) expression.

Setup & QC: Loads data, then creates a Seurat object, and performs QC filtering (nFeature, nCount, %mt).

Core Analysis: Normalizes, scales, runs PCA, UMAP, and clustering.

Ly6a Grouping: Subsets Ly6a-positive cells and divides them into Low, Medium, and High expression groups (tertiles).

Cell Cycle: Scores cell cycle (G1, S, G2/M) and plots its distribution across the Ly6a groups.

GSEA: Finds Differentially Expressed Genes (DEGs) between Ly6a High and Low cells. 
Runs Gene Set Enrichment Analysis (GSEA) on the ranked DEGs (using GO:BP) and generates a divergent bar plot of the top enriched pathways.

ORA: Also performs standard GO Over-Representation Analysis (ORA) on significant up/downregulated genes.
