###############################
# Install Required Packages
###############################
###############################
# Install Packages
###############################

# Install BiocManager if needed
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# CRAN packages
cran_pkgs <- c(
  "Seurat",
  "Signac",
  "tidyverse",
  "hdf5r",
  "ggpubr",
  "rstatix",
  "stringr",
  "UpSetR"
)

for (pkg in cran_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    install.packages(pkg)
}

# Bioconductor packages
bioc_pkgs <- c(
  "clusterProfiler",
  "org.Mm.eg.db",
  "JASPAR2020",
  "TFBSTools",
  "BSgenome.Mmusculus.UCSC.mm10"
)

for (pkg in bioc_pkgs) {
  if (!requireNamespace(pkg, quietly = TRUE))
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
}

###############################
# Load Libraries
###############################

library(Seurat)
library(Signac)
library(tidyverse)
library(hdf5r)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggpubr)
library(rstatix)
library(stringr)
library(ggplot2)
library(dplyr)
library(tidyr)
library(JASPAR2020)
library(TFBSTools)
library(BSgenome.Mmusculus.UCSC.mm10)
library(UpSetR)
###############################
# 2. Load 10X Data
###############################
data_dir  <- "D:/LSK_DB/outs"
h5_path   <- file.path(data_dir, "filtered_feature_bc_matrix.h5")
frag_path <- file.path(data_dir, "atac_fragments.tsv.gz")

###################################################################################################
# 3. Load multiome data (RNA + ATAC)
###################################################################################################
counts_multi <- Read10X_h5(h5_path)
str(counts_multi)   

rna_counts  <- counts_multi$`Gene Expression`
peak_counts <- counts_multi$Peaks

# Keep only standard chromosomes in the peak set 
grange_all  <- StringToGRanges(rownames(peak_counts), sep = c(":", "-"))
std_chroms  <- as.logical(seqnames(grange_all) %in% standardChromosomes(grange_all))
peak_counts <- peak_counts[std_chroms, ]

###################################################################################################
# 4. Gene annotation for the ATAC assay
###################################################################################################
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
seqlevelsStyle(annotations) <- "UCSC"
genome(annotations) <- "mm10"

###############################
# 3. Create Seurat Object & QC Parameters
###############################

seurat_obj <- CreateSeuratObject(counts = rna_counts, project = "LSK_Multiome", assay = "RNA")

chrom_assay <- CreateChromatinAssay(
  counts     = peak_counts,
  sep        = c(":", "-"),
  fragments  = frag_path,
  genome     = "mm10",
  annotation = annotations,
  min.cells  = 10
)

# Restrict to barcodes present in both modalities before attaching
common_cells <- intersect(colnames(seurat_obj), colnames(chrom_assay))
seurat_obj   <- subset(seurat_obj, cells = common_cells)
seurat_obj[["ATAC"]] <- subset(chrom_assay, cells = common_cells)
seurat_obj
###################################################################################################
# 4. Joint QC (RNA + ATAC) and filtering
###################################################################################################
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")

DefaultAssay(seurat_obj) <- "ATAC"
seurat_obj <- NucleosomeSignal(seurat_obj)
seurat_obj <- TSSEnrichment(seurat_obj, fast = TRUE)

# QC thresholds - tune to your data (check VlnPlot(seurat_obj, features=...) first if unsure)
min_nFeature   <- 200
max_nFeature   <- 4000
min_nCount     <- 200
max_nCount     <- 10000
max_percent_mt <- 10
min_nCount_ATAC <- 1000
max_nCount_ATAC <- 100000
max_nucleosome  <- 2
min_TSS         <- 2

VlnPlot(seurat_obj,
        features = c("nFeature_RNA", "nCount_RNA", "percent.mt",
                     "nCount_ATAC", "TSS.enrichment", "nucleosome_signal"),
        ncol = 3, pt.size = 0)

seurat_obj <- subset(
  seurat_obj,
  subset = nFeature_RNA > min_nFeature & nFeature_RNA < max_nFeature &
    nCount_RNA > min_nCount & nCount_RNA < max_nCount &
    percent.mt < max_percent_mt &
    nCount_ATAC > min_nCount_ATAC & nCount_ATAC < max_nCount_ATAC &
    nucleosome_signal < max_nucleosome & TSS.enrichment > min_TSS
)

###############################
# 5. Downstream Analysis
###############################
DefaultAssay(seurat_obj) <- "RNA"
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst")
seurat_obj <- RunPCA(seurat_obj, npcs = 30, verbose = FALSE)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30, reduction.name = "umap.rna")
DimPlot(seurat_obj, reduction = "umap.rna")
###############################
# 6. Ly6a & Kit Expression
###############################
ly6a_expr <- FetchData(seurat_obj, "Ly6a")
kit_expr  <- FetchData(seurat_obj, "Kit")
sum(ly6a_expr$Ly6a > 0)
sum(kit_expr$Kit > 0)
cor(ly6a_expr$Ly6a, kit_expr$Kit)
VlnPlot(seurat_obj, features = c("Ly6a", "Kit"), pt.size = 0.1)

###################################################################################################
# 7. Cell cycle scoring
###################################################################################################
cc.genes <- Seurat::cc.genes.updated.2019
cc.genes.mouse <- lapply(cc.genes, function(genes) str_to_title(genes))

seurat_obj <- CellCycleScoring(seurat_obj,
                               s.features = cc.genes.mouse$s.genes,
                               g2m.features = cc.genes.mouse$g2m.genes,
                               assay = "RNA", layer = "data")
seurat_obj$Phase <- factor(seurat_obj$Phase, levels = c("G1", "S", "G2M"))
DimPlot(seurat_obj, reduction = "umap.rna", group.by = "Phase",
        cols = c("G1" = "red", "S" = "green", "G2M" = "blue"), pt.size = 1)

###################################################################################################
# 8. Ly6a Low / Medium / High grouping (Ly6a-positive cells only)
###################################################################################################
ly6a_present_cell <- subset(seurat_obj, subset = Ly6a > 0)
ly6a_vals <- FetchData(ly6a_present_cell, "Ly6a")$Ly6a

cuts <- quantile(ly6a_vals, probs = c(0, 1/3, 2/3, 1), na.rm = TRUE)
cuts <- unique(cuts)

ly6a_present_cell$expression_group <- cut(
  ly6a_vals, breaks = cuts,
  labels = c("Low", "Medium", "High")[seq_len(length(cuts) - 1)],
  include.lowest = TRUE
)
table(ly6a_present_cell$expression_group,ly6a_present_cell$Phase)

# Check distribution
table(ly6a_present_cell$expression_group)

# Visualize Ly6a groups
DimPlot(ly6a_present_cell, group.by = "expression_group",
        cols = c("Low"="green","Medium"="blue","High"="red"), pt.size=1)

###############################
# 7. Cell Cycle Distribution by Ly6a Group
###############################
counts <- table(ly6a_present_cell$expression_group, ly6a_present_cell$Phase)

# 2. Convert to data frame, calculate row percentages, and format factors
df <- as.data.frame(counts) %>%
  rename(
    expression_group = Var1,
    Phase = Var2,
    count = Freq
  ) %>%
  group_by(expression_group) %>%
  mutate(
    percent = (count / sum(count)) * 100,
    # Label formatting: strip trailing .0 for whole percentages (like 73%)
    label = paste0(round(percent, 1), "%"),
    label = gsub("\\.0%", "%", label) 
  ) %>%
  ungroup() %>%
  mutate(
    # Set factor levels to ensure exact axis ordering
    expression_group = factor(expression_group, levels = c("Low", "Medium", "High")),
    Phase = factor(Phase, levels = c("G1", "S", "G2M"), labels = c("G1", "S", "G2/M"))
  )

# 3. Generate Heatmap
ggplot(df, aes(x = Phase, y = expression_group, fill = percent)) +
  geom_tile(color = "white", linewidth = 1) +
  geom_text(aes(label = label), color = "black", size = 5.5) +
  scale_fill_gradient(
    low = "#edf3f8", 
    high = "#4682b4", 
    name = "percent",
    limits = c(0, 80),
    breaks = c(20, 40, 60)
  ) +
  labs(
    title = "Heatmap of Cell Cycle Phase Distribution by Ly6a Expression Group",
    x = "Cell Cycle Phase",
    y = "Ly6a Expression Group"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(hjust = 0.5, face = "plain", size = 16),
    panel.grid = element_blank(),
    
    
    axis.text.y = element_text(color = "black", size = 16, face = "bold"), # High/Medium/Low size
    axis.text.x = element_text(color = "black", size = 14),                # G1/S/G2/M size
    axis.title = element_text(size = 14),                                 # Axis titles size
    
    legend.title = element_text(size = 11),
    legend.position = "right"
  )
DimPlot(ly6a_present_cell, reduction = "umap.rna", group.by = "expression_group",
        cols = c("Low" = "green", "Medium" = "blue", "High" = "red"), pt.size = 2)+ggtitle("Ly6a Expression group")

# Bar plot
ggplot(cellcycle_summary, aes(x = Phase, y = percent, fill = expression_group)) +
  geom_bar(stat="identity", position=position_dodge()) +
  scale_fill_manual(values=c("Low"="#228B22","Medium"="blue","High"="red")) +
  ylab("Percentage") + xlab("Cell Cycle Phase") +
  theme_minimal()

 #  Ly6a  Cell Cycle Distribution                      
ggplot(ly6a_present_cell@meta.data, aes(x = expression_group, fill = Phase)) +
  geom_bar(position = "fill") +
  theme_minimal() +
  labs(y = "Proportion", title = "Cell Cycle Distribution in Ly6a High, Medium, and Low")

#################################################################################################################################################################################################
library(tidyverse)
library(hdf5r)
library(clusterProfiler)
library(org.Mm.eg.db)
library(enrichplot)

###############################
# 8. Differential Expression
###############################
hl_cells <- subset(ly6a_present_cell, subset = expression_group %in% c("Low", "High"))
Idents(hl_cells) <- hl_cells$expression_group

markers <- FindMarkers(
  hl_cells,
  ident.1 = "High",
  ident.2 = "Low",
  logfc.threshold = 0.25,
  min.pct = 0.1
)

###############################
# 9. GSEA / GO Analysis
###############################
head(markers)

markers <- as.data.frame(markers) %>%
  rownames_to_column("gene")

# Rank genes by avg_log2FC (descending)
gsea_ranking <- markers %>%
  dplyr::filter(!is.na(avg_log2FC)) %>%
  dplyr::arrange(desc(avg_log2FC)) %>%
  dplyr::select(gene, avg_log2FC)

# Remove duplicate genes (keeps highest logFC copy because of prior sort)
gsea_ranking <- gsea_ranking[!duplicated(gsea_ranking$gene), ]

# Convert to named numeric vector, sorted strictly decreasing
ranked_genes <- gsea_ranking$avg_log2FC
names(ranked_genes) <- gsea_ranking$gene
ranked_genes <- sort(ranked_genes, decreasing = TRUE)

head(ranked_genes)

# -----------------------------
#  Run GSEA using GO:BP
# -----------------------------
gsea_go <- gseGO(
  geneList = ranked_genes,
  OrgDb = org.Mm.eg.db,
  ont = "BP",
  keyType = "SYMBOL",
  minGSSize = 10,
  maxGSSize = 500,
  pvalueCutoff = 0.05,
  verbose = TRUE
)

# -----------------------------
#  Plot top results (dotplot)
# -----------------------------
dotplot(gsea_go, showCategory = 20) +
  ggtitle("GSEA (GO Biological Process): Ly6a High vs Low") +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.text.x = element_text(size = 13),
    axis.text.y = element_text(size = 13),
    axis.title = element_text(size = 15),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 12)
  )

# Save results to .csv
write.csv(gsea_go@result, "Ly6a_gsea_go_results.csv", row.names = FALSE)

# ------------------------------------------
# Select top 15 Positive and Negative NES
# ------------------------------------------
# gsea_go is an S4 gseaResult object -- convert to data frame before dplyr verbs
gsea_df <- as.data.frame(gsea_go)

top_positive <- gsea_df %>%
  arrange(desc(NES)) %>%
  head(15)

top_negative <- gsea_df %>%
  arrange(NES) %>%
  head(15)

# --------------------------------------------
# Merge positive + negative (dedupe just in case of overlap)
# --------------------------------------------
merged_gseaall <- rbind(top_positive, top_negative) %>%
  distinct(ID, .keep_all = TRUE)

# Format p.adjust for labels
merged_gseaall$p_label <- format(
  merged_gseaall$p.adjust,
  scientific = TRUE,
  digits = 2
)

# ----------------------------------------------------
# Plot in a single combined panel with p-values
# ----------------------------------------------------
ggplot(merged_gseaall,
       aes(x = reorder(Description, NES),
           y = NES,
           fill = NES)) +

  geom_bar(stat = "identity") +

  # Add p-values on plot
  geom_text(aes(label = p_label),
            hjust = ifelse(merged_gseaall$NES > 0, -0.15, 1.15),
            size = 4.5) +

  coord_flip() +

  scale_fill_gradient2(
    low = "blue",
    mid = "white",
    high = "red",
    midpoint = 0
  ) +

  labs(
    title = "Top 15 GO Terms Enriched in Ly6a High vs Ly6a Low LSK Cells",
    x = "GO Terms",
    y = "Normalized Enrichment Score (NES)",
    fill = "NES"
  ) +

  theme_minimal(base_size = 16) +
  theme(
    axis.text.y = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    axis.title = element_text(size = 15, face = "bold"),
    plot.title = element_text(size = 18, face = "bold"),
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12)
  ) +

  # Make room for labels outside bars (both directions)
  expand_limits(y = c(min(merged_gseaall$NES) * 1.2,
                       max(merged_gseaall$NES) * 1.2))

 ##############GO Analysis #####################################
                         
# Convert gene symbols to Entrez IDs
Idents(ly6a_present_cell) <- "expression_group"  # Set group identity
markers <- FindMarkers(ly6a_present_cell, ident.1 = "High", ident.2 = "Low", logfc.threshold = 0.25)
head(markers)
gene_list <- markers %>% rownames_to_column("gene") %>%
  dplyr::filter(p_val_adj < 0.05) %>%
  dplyr::arrange(desc(avg_log2FC)) %>%
  dplyr::pull(gene)

entrez_ids <- mapIds(org.Mm.eg.db, keys = gene_list, keytype = "SYMBOL", column = "ENTREZID")

# Perform GO enrichment analysis
go_results <- enrichGO(gene = entrez_ids, OrgDb = org.Mm.eg.db, keyType = "ENTREZID", ont = "BP")

# Plot results
dotplot(go_results)

table(ly6a_present_cell$expression_group, ly6a_present_cell$Phase)

summary(FetchData(ly6a_present_cell, vars = "Cd19"))
sum(FetchData(ly6a_present_cell, vars = "Cd19") > 0)

rownames(ly6a_present_cell)


# Upregulated in Sca-1 High
#genes_high <- rownames(markers[markers$avg_log2FC > 0 & markers$p_val_adj < 0.05, ])
upregulated_genes <- rownames(markers[markers$avg_log2FC > 0, ])

# Upregulated in Sca-1 Low
genes_low <- rownames(markers[markers$avg_log2FC < 0 & markers$p_val_adj < 0.05, ])
downregulated_genes <- rownames(markers[markers$avg_log2FC < 0, ])
library(clusterProfiler)
library(org.Mm.eg.db)  # Mouse genome annotation

go_low1 <- enrichGO(gene = downregulated_genes,
                   OrgDb = org.Mm.eg.db,
                   keyType = "SYMBOL",
                   ont = "BP",  # Biological Process
                   pAdjustMethod = "BH",
                   pvalueCutoff = 0.05)


go_high <- enrichGO(gene = upregulated_genes, 
                    OrgDb = org.Mm.eg.db, 
                    keyType = "SYMBOL", 
                    ont = "BP",   # "BP" for Biological Process
                    pAdjustMethod = "BH", 
                    pvalueCutoff = 0.05)

# View top enriched terms
head(go_high@result)
write.csv(go_high,"go_highv2.csv")
# Visualization
dotplot(go_high, showCategory = 10)  # Show top 10 GO terms
go_low <- enrichGO(gene = genes_low, 
                   OrgDb = org.Mm.eg.db, 
                   keyType = "SYMBOL", 
                   ont = "BP",   
                   pAdjustMethod = "BH", 
                   pvalueCutoff = 0.05)

# View results
head(go_low@result)
write.csv(go_low1,"go_lowv2.csv")
# Visualization
dotplot(go_low, showCategory = 10)

#############################################################################################################################
# 10. Ly6a expression vs cell cycle activity scores
#############################################################################################################################
ly6a_pos <- subset(seurat_obj, subset = Ly6a > 0)
ly6a_pos
df <- FetchData(ly6a_pos, vars = c("Ly6a", "S.Score", "G2M.Score"))
df$Ly6a_log <- log1p(df$Ly6a)
cor.test(df$Ly6a_log, df$S.Score, method = "pearson")

cor.test(df$Ly6a_log, df$G2M.Score, method = "pearson")



p_s <- ggplot(df, aes(x = log1p(Ly6a), y = S.Score)) +
  geom_point(alpha = 0.25, size = 0.7) +
  geom_smooth(method = "lm", se = TRUE) +
  annotate("text", x = max(log1p(df$Ly6a)), y = max(df$S.Score),
           label = "r = -0.296, p < 2.2e-16", hjust = 1, size = 5) +
  labs(
    x = bquote("log1p("*italic("Ly6a")*")"), 
    y = "S-phase score",
    title = bquote("Association between" ~ italic("Ly6a") ~ "and S activity")
  ) + # Fixed: Added closing parenthesis here
  theme_minimal(base_size = 16)

p_g2m <- ggplot(df, aes(x = log1p(Ly6a), y = G2M.Score)) +
  geom_point(alpha = 0.25, size = 0.7) +
  geom_smooth(method = "lm", se = TRUE) +
  annotate("text", x = max(log1p(df$Ly6a)), y = max(df$G2M.Score),
           label = "r = -0.301, p < 2.2e-16", hjust = 1, size = 5) +
  labs(
    # Use bquote to italicize specific text
    x = bquote("log1p("*italic("Ly6a")*")"), 
    y = "G2/M score",
    title = bquote("Association between" ~italic("Ly6a") ~ "and G2/M activity")
  ) +
  theme_minimal(base_size = 16)

patchwork::wrap_plots(p_s, p_g2m)

patchwork::wrap_plots(pear_s, pear_g2m)

##############################################################################################################################
11. Process ATAC Space Specifically for Ly6a-Present Cells
##############################################################################################################################
DefaultAssay(ly6a_present_cell) <- "ATAC"

# Run the standard Signac ATAC workflow on this subset
ly6a_present_cell <- RunTFIDF(ly6a_present_cell)
ly6a_present_cell <- FindTopFeatures(ly6a_present_cell, min.cutoff = "q0")
ly6a_present_cell <- RunSVD(ly6a_present_cell, reduction.name = "lsi", reduction.key = "LSI_")

# Run UMAP on the ATAC assay for these specific cells (dropping LSI component 1)
ly6a_present_cell <- RunUMAP(
  ly6a_present_cell, 
  reduction = "lsi", 
  dims = 2:30, 
  reduction.name = "umap.atac", 
  reduction.key = "ATACUMAP_"
)

# Harmonize the cell cycle factor labels
ly6a_present_cell$Phase <- factor(ly6a_present_cell$Phase, levels = c("G1", "S", "G2M"))

p1 <- DimPlot(
  ly6a_present_cell, 
  reduction = "umap.atac", 
  group.by = "expression_group", 
  pt.size = 1.2
) +
  scale_color_manual(values = c("Low" = "green", "Medium" = "blue", "High" = "red")) +
  theme_minimal() +
  ggtitle("ATAC Space: Ly6a Expression Groups") +
  theme(plot.title = element_text(face = "bold", size = 12))

# Plot B: Co-mapped cell cycle phases for the exact same cells
p2 <- DimPlot(
  ly6a_present_cell, 
  reduction = "umap.atac", 
  group.by = "Phase", 
  pt.size = 1.2
) +
  scale_color_manual(values = c("G1" = "red", "S" = "green", "G2M" = "blue")) +
  theme_minimal() +
  ggtitle("ATAC Space: RNA Cell Cycle") +
  theme(plot.title = element_text(face = "bold", size = 12))

# Display side-by-side
print(p1 + p2)
# For RNA Assay Comparison

DefaultAssay(ly6a_present_cell) <- "RNA"
p1 <- DimPlot(
  ly6a_present_cell, 
  reduction = "umap.rna", 
  group.by = "expression_group", 
  pt.size = 1.2
) +
  scale_color_manual(values = c("Low" = "green", "Medium" = "blue", "High" = "red")) +
  theme_minimal() +
  ggtitle("RNA Space: Ly6a Expression Groups") +
  theme(plot.title = element_text(face = "bold", size = 12))
p2<-DimPlot(
  ly6a_present_cell, 
  reduction = "umap.rna", 
  group.by = "Phase", 
  pt.size = 1.2
) +
  scale_color_manual(values = c("G1" = "red", "S" = "green", "G2M" = "blue")) +
  theme_minimal() +
  ggtitle("RNA Space: Cell Cycle Phase Distribution") +
  theme(plot.title = element_text(face = "bold", size = 12))

print (p1+p2)

############## Cell-Cycle Split ATAC Assay###############################################
                         
DefaultAssay(ly6a_present_cell) <- "ATAC"
p_split <- DimPlot(
  ly6a_present_cell, 
  reduction = "umap.atac", 
  group.by = "expression_group", 
  split.by = "Phase",
  pt.size = 1
) +
  scale_color_manual(values = c("Low" = "green", "Medium" = "blue", "High" = "red")) +
  theme_bw() +
  ggtitle("ATAC Landscape of Ly6a Groups Split by Phase")
p_split

Idents(ly6a_present_cell) <- ly6a_present_cell$expression_group
#####################################################################################
12. Peak Difference Plots
#####################################################################################                        
library(Signac)
library(Seurat)
library(VennDiagram)
library(grid)

DefaultAssay(ly6a_present_cell) <- "ATAC"

# Get peak accessibility matrix
peak_matrix <- GetAssayData(ly6a_present_cell, slot = "counts")

# Cells in each group
low_cells  <- WhichCells(ly6a_present_cell, idents = "Low")
med_cells  <- WhichCells(ly6a_present_cell, idents = "Medium")
high_cells <- WhichCells(ly6a_present_cell, idents = "High")

# Peaks present in at least 5% of cells
min_pct <- 0.05

low_peaks <- rownames(peak_matrix)[
  Matrix::rowMeans(peak_matrix[, low_cells] > 0) >= min_pct
]

med_peaks <- rownames(peak_matrix)[
  Matrix::rowMeans(peak_matrix[, med_cells] > 0) >= min_pct
]

high_peaks <- rownames(peak_matrix)[
  Matrix::rowMeans(peak_matrix[, high_cells] > 0) >= min_pct
]

cat("Low peaks:", length(low_peaks), "\n")
cat("Medium peaks:", length(med_peaks), "\n")
cat("High peaks:", length(high_peaks), "\n")

venn.plot <- venn.diagram(
  x = list(
    Low = low_peaks,
    Medium = med_peaks,
    High = high_peaks
  ),
  filename = NULL,
  fill = c("forestgreen", "royalblue", "firebrick"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.4,
  cat.fontface = "bold"
)

grid.newpage()
grid.draw(venn.plot)
                         
#####################################################################################
   13. G1 Phase Only                      
####################################################################################                       
g1_cells <- subset(
  ly6a_present_cell,
  subset = Phase == "G1"
)

Idents(g1_cells) <- g1_cells$expression_group
DefaultAssay(g1_cells) <- "ATAC"

peak_matrix <- GetAssayData(
  g1_cells,
  slot = "counts"
)

low_cells <- WhichCells(g1_cells, idents = "Low")
med_cells <- WhichCells(g1_cells, idents = "Medium")
high_cells <- WhichCells(g1_cells, idents = "High")


min_pct <- 0.05

low_peaks <- rownames(peak_matrix)[
  Matrix::rowMeans(peak_matrix[,low_cells] > 0) >= min_pct
]

med_peaks <- rownames(peak_matrix)[
  Matrix::rowMeans(peak_matrix[,med_cells] > 0) >= min_pct
]

high_peaks <- rownames(peak_matrix)[
  Matrix::rowMeans(peak_matrix[,high_cells] > 0) >= min_pct
]


cat("G1 Low peaks:", length(low_peaks), "\n")
cat("G1 Medium peaks:", length(med_peaks), "\n")
cat("G1 High peaks:", length(high_peaks), "\n")

venn.plot <- venn.diagram(
  x = list(
    Low = low_peaks,
    Medium = med_peaks,
    High = high_peaks
  ),
  filename = NULL,
  fill = c("forestgreen", "royalblue", "firebrick"),
  alpha = 0.5,
  cex = 1.5,
  cat.cex = 1.4,
  cat.fontface = "bold",
  main = "G1 Ly6a ATAC Accessible Peak Overlap"
)

grid.newpage()
grid.draw(venn.plot)


######################UpSet Plot#############################################

peak_list <- fromList(
  list(
    Low = low_peaks,
    Medium = med_peaks,
    High = high_peaks
  )
)

upset(
  peak_list,
  
  nsets = 3,
  nintersects = NA,
  keep.order = TRUE,
  order.by = "freq",
  
  # Larger dots and connecting lines
  point.size = 11,
  line.size = 4,
  
  # Bar colors
  main.bar.color = "black",
  sets.bar.color = "steelblue4",
  
  # Axis labels
  mainbar.y.label = "G1 Phase: Number of Accessible ATAC Peaks",
  sets.x.label = "Accessible ATAC peaks",
  
  # Show counts above bars
  show.numbers = "yes",
  
  # Increase all text sizes
  text.scale = c(
    2.5,  # Main bar title
    2.2,  # Main bar tick labels
    2.2,  # Set size title
    1.7,  # Set names
    2.8,  # Numbers above bars
    2.8   # Matrix labels
  )
)
################################################################################
# 14. Generate the G1-Specific Differential Peaks (DA Peaks) Object
################################################################################

# Isolate only the resting G1 population to eliminate proliferation bias
g1_cells <- subset(ly6a_present_cell, subset = Phase == "G1")

# Set the active identity to your Ly6a expression groups
Idents(g1_cells) <- g1_cells$expression_group
g1_cells$expression_group
# We include 'nCount_ATAC' as a latent variable to correct for sequencing depth differences.
 # High vs Low
                         
da_peaks_g1_high_low <- FindMarkers(
  object = g1_cells,
  ident.1 = "High",
  ident.2 = "Low",
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)


write.csv(da_peaks_g1_high_low,"Ly6a_da_peaks_g1_High_Low.csv")
 # Medium vs Low
                         
  da_peaks_g1_medium_low <- FindMarkers(
  object = g1_cells,
  ident.1 = "Medium",
  ident.2 = "Low",
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)

write.csv(da_peaks_g1_medium_low,"Ly6a_da_peaks_g1_Medium_low.csv")
                         
 # High vs Medium
                         
da_peaks_g1_high_medium <- FindMarkers(
  object = g1_cells,
  ident.1 = "High",
  ident.2 = "Medium",
  test.use = 'LR',
  latent.vars = 'nCount_ATAC'
)

write.csv(da_peaks_g1_high_medium,"Ly6a_da_peaks_g1_High_Medium.csv")

################################################################################
# 1. Select candidate DA peaks
################################################################################                         
da_peaks_g1<-da_peaks_g1_high_low     
                         
 top_high_peaks <- rownames(
  da_peaks_g1[
    da_peaks_g1$avg_log2FC > 0.25 &
      da_peaks_g1$p_val < 0.01 &
      da_peaks_g1$pct.1 > 0.05,
  ]
)

# Ly6a-Low accessible peaks
top_low_peaks <- rownames(
  da_peaks_g1[
    da_peaks_g1$avg_log2FC < -0.25 &
      da_peaks_g1$p_val < 0.01 &
      da_peaks_g1$pct.2 > 0.05,
  ]
)

print(paste("Ly6a-High candidate peaks:", length(top_high_peaks)))
print(paste("Ly6a-Low candidate peaks:", length(top_low_peaks)))


################################################################################
# 2. Run motif enrichment
################################################################################

DefaultAssay(g1_cells) <- "ATAC"

g1_cells <- AddMotifs(
  object = g1_cells,
  genome = BSgenome.Mmusculus.UCSC.mm10,
  pfm = pfm
)

motif_matrix <- GetMotifData(
  object = g1_cells,
  assay = "ATAC",
  slot = "data"
)

dim(motif_matrix)

head(rownames(motif_matrix))
sum(top_high_peaks %in% rownames(g1_cells[["ATAC"]]))
length(top_high_peaks)
# Background = all tested peaks
background_peaks <- rownames(da_peaks_g1)

# High motif enrichment
enriched_motifs_high <- FindMotifs(
  object = g1_cells,
  features = top_high_peaks,
  background = background_peaks
)


# Low motif enrichment
enriched_motifs_low <- FindMotifs(
  object = g1_cells,
  features = top_low_peaks,
  background = background_peaks
)


################################################################################
# 3. Save motif tables
################################################################################

write.csv(
  enriched_motifs_high,
  "Ly6a_High_ATAC_motif_enrichment(G1_cells).csv",
  row.names = FALSE
)

write.csv(
  enriched_motifs_low,
  "Ly6a_Low_ATAC_motif_enrichment_(G1 cells).csv",
  row.names = FALSE
)


################################################################################
# 4. Merge High and Low motif results
################################################################################

high_data <- enriched_motifs_high
high_data$motif_id <- rownames(high_data)

low_data <- enriched_motifs_low
low_data$motif_id <- rownames(low_data)


merged_data <- merge(
  high_data,
  low_data,
  by = "motif_id",
  suffixes = c("_high", "_low")
)


################################################################################
# 5. Calculate High vs Low motif preference
################################################################################

merged_data$log2_enrichment_ratio <- log2(
  merged_data$fold.enrichment_high /
    merged_data$fold.enrichment_low
)


# significance score
merged_data$log_p_min <- -log10(
  pmin(
    merged_data$pvalue_high,
    merged_data$pvalue_low
  )
)


# Significant if enriched on either side
merged_data$is_sig <- (
  merged_data$p.adjust_high < 0.05 |
    merged_data$p.adjust_low < 0.05
)


# Ranking score for labeling
merged_data$ranking_score <-
  abs(merged_data$log2_enrichment_ratio) *
  merged_data$log_p_min


################################################################################
# 6. Select top TFs to label
################################################################################

top_TFs <- merged_data %>%
  arrange(desc(ranking_score)) %>%
  slice_head(n = 25)


merged_data$label <- NA

merged_data$label[
  merged_data$motif_id %in% top_TFs$motif_id
] <- merged_data$motif.name_high[
  match(
    merged_data$motif_id[
      merged_data$motif_id %in% top_TFs$motif_id
    ],
    merged_data$motif_id
  )
]

# Keep only significant motifs
sig_motifs <- subset(
  merged_data,
  is_sig == TRUE
)
################################################################################
# 7. Plot two-sided motif volcano
################################################################################
y = max(merged_data$log_p_min) - 0.5
p_motif_volcano <- ggplot(
  sig_motifs,
  aes(
    x = log2_enrichment_ratio,
    y = log_p_min
  )
) +
  
  geom_point(
    aes(color = is_sig),
    size = 3,
    alpha = 0.8
  ) +
  
  scale_color_manual(
    values = c(
      "TRUE" = "firebrick",
      "FALSE" = "grey70"
    )
  ) +
  
  geom_text_repel(
    aes(label = label),
    size = 6,
    max.overlaps = Inf,
    box.padding = 0.5,
    segment.color = "grey50"
  ) +
  
  geom_vline(
    xintercept = 0,
    linetype = "dashed"
  ) +
  
  annotate(
    "text",
    x = min(merged_data$log2_enrichment_ratio),
    y = max(merged_data$log_p_min),
    label = "",
    hjust = 0,
    size = 6,
    color = "#667543",
    fontface = "bold"
  ) +
  
  annotate(
    "text",
    x = max(merged_data$log2_enrichment_ratio),
    y = max(merged_data$log_p_min),
    label = "",
    hjust = 1,
    size = 6,
    color = "black",
    fontface = "bold"
  ) +
  
  labs(
    title = "Ly6a-High vs Ly6a-Low G1 ATAC Motif Enrichment",
    subtitle = "DA peak cutoff: |log2FC| > 0.25; p < 0.01; accessibility >5% cells",
    x = "log2(Fold enrichment High/ Low)",
    y = "-log10(minimum motif p-value)",
    color = "Adjusted significant"
  ) +
  theme_classic() +
  theme(
    plot.title = element_text(
      size = 18,
      face = "bold",
      hjust = 0.5
    ),
    axis.title = element_text(
      size = 18,
      face = "bold"
    ),
    plot.subtitle = element_text(
      size = 14,
      face = "italic",
      hjust = 0.5
    ),
    axis.text = element_text(
      size = 14
    ),
    legend.title = element_text(
      size = 14,
      face = "bold"
    ),
    legend.text = element_text(
      size = 14
    ),
    axis.text.x = element_text(
      size = 18,
      face = "bold"
    ),
    axis.text.y = element_text(
      size = 18,
      face = "bold"
    )
  )


print(p_motif_volcano)
                        
#############################################################
 Print Sequence Logos of Top TFs                        
##############################################################
top_motif_ids <- top_TFs$motif_id[1:4] # Select top 4 motif IDs

# Generate and customize the MotifPlot
p_motifs <- MotifPlot(
  object = g1_cells,
  motifs = top_motif_ids,
  assay = "ATAC"
) + 
  theme(
    # Motif title names (e.g., TF name / Motif ID above each logo)
    strip.text = element_text(size = 18, face = "bold"),
    
    # Axis labels ("Position", "Bits", etc.)
    axis.title = element_text(size = 16, face = "bold"),
    
    # Axis tick labels (position numbers: 1, 2, 3...)
    axis.text = element_text(size = 14, face = "bold")
  )

print(p_motifs)       
###############################################################
Zfx in RNA Assay
###############################################################  


                         library(Seurat)
library(ggplot2)

DefaultAssay(g1_cells) <- "RNA"
Idents(g1_cells) <- g1_cells$expression_group
 zfx_data <- FetchData(
   g1_cells,
   vars = c("Zfx", "expression_group")
 )
 p_low_high <- wilcox.test(
       Zfx ~ expression_group,
      data = subset(zfx_data, expression_group %in% c("Low","High"))
    )$p.value
 p_low_med <- wilcox.test(
     Zfx ~ expression_group,
     data = subset(zfx_data, expression_group %in% c("Low","Medium"))
   )$p.value
 p_med_high <- wilcox.test(
     Zfx ~ expression_group,
     data = subset(zfx_data, expression_group %in% c("Medium","High"))
   )$p.value
 
 # Print p-values
 p_low_high

p_low_med

p_med_high

# VlnPlot
p_zfx <- VlnPlot(
  g1_cells,
  features = "Zfx",
  group.by = "expression_group",
  pt.size = 0
) +
  theme_classic() +
  ggtitle("Zfx RNA Expression in G1 Ly6a Groups") +
  theme(
    plot.title = element_text(size = 18, face = "bold"),
    axis.text = element_text(size = 12)
  )


                         
