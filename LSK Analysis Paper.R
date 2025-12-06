###############################
# Install Required Packages
###############################
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")

cran_pkgs <- c("tidyverse", "hdf5r", "ggpubr", "rstatix", "stringr")
for (p in cran_pkgs) if (!requireNamespace(p, quietly = TRUE)) install.packages(p)

if (!requireNamespace("Seurat", quietly = TRUE)) install.packages("Seurat")

if (!requireNamespace("clusterProfiler", quietly = TRUE)) BiocManager::install("clusterProfiler")
if (!requireNamespace("org.Mm.eg.db", quietly = TRUE)) BiocManager::install("org.Mm.eg.db")

###############################
# 1. Load Libraries
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

###############################
# 2. Load 10X Data
###############################
counts <- Read10X_h5("C:/Users/lopde33/RNADB/filtered_feature_bc_matrix.h5")
str(counts)

###############################
# 3. Create Seurat Object & QC Parameters
###############################

seurat_obj <- CreateSeuratObject(counts = counts$`Gene Expression`, project = "LSK_Analysis", assay = "RNA")
seurat_obj[["percent.mt"]] <- PercentageFeatureSet(seurat_obj, pattern = "^mt-")
min_nFeature <- 200
max_nFeature <- 4000
min_nCount <- 200
max_nCount <- 10000
max_percent_mt <- 10



# QC 
seurat_obj <- subset(seurat_obj, 
                     subset = nFeature_RNA > min_nFeature & nFeature_RNA < max_nFeature &
                       nCount_RNA > min_nCount & nCount_RNA < max_nCount &
                       percent.mt < max_percent_mt)

# Normalize and scale
seurat_obj <- NormalizeData(seurat_obj)
seurat_obj <- ScaleData(seurat_obj)

###############################
# 3. Downstream Analysis
###############################
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst")
seurat_obj <- RunPCA(seurat_obj, npcs = 30, verbose = FALSE)
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:30)
seurat_obj <- FindClusters(seurat_obj, resolution = 0.5)
seurat_obj <- RunUMAP(seurat_obj, dims = 1:30)
DimPlot(seurat_obj)

###############################
# 4. Ly6a & Kit Expression
###############################
VlnPlot(seurat_obj, features = c("Ly6a", "Kit"), pt.size = 0.1)
ly6a_expr <- FetchData(seurat_obj, "Ly6a")
kit_expr <- FetchData(seurat_obj, "Kit")
cor(ly6a_expr$Ly6a, kit_expr$Kit)

###############################
# 5. Cell Cycle Scoring
###############################
cc.genes <- Seurat::cc.genes.updated.2019
cc.genes.mouse <- lapply(cc.genes, function(genes) str_to_title(genes))

seurat_obj <- CellCycleScoring(seurat_obj, 
                               s.features = cc.genes.mouse$s.genes,
                               g2m.features = cc.genes.mouse$g2m.genes,
                               assay = "RNA", layer = "data")

seurat_obj$Phase <- factor(seurat_obj$Phase, levels = c("G1", "S", "G2/M"))
DimPlot(seurat_obj, group.by = "Phase", cols = c("G1"="red","S"="green","G2/M"="blue"), pt.size=1)

###############################
# 6. Ly6a Positive Cells & Grouping
###############################
ly6a_present_cell <- subset(seurat_obj, subset = Ly6a > 0)
ly6a_present_cell$Ly6a_expression <- FetchData(ly6a_present_cell, "Ly6a")$Ly6a

ly6a_present_cell$expression_group <- cut(ly6a_present_cell$Ly6a_expression,
                                          breaks = quantile(ly6a_present_cell$Ly6a_expression, probs = 0:3/3),
                                          labels = c("Low", "Medium", "High"),
                                          include.lowest = TRUE)

table(ly6a_present_cell$expression_group)

# Visualize Ly6a groups
DimPlot(ly6a_present_cell, group.by = "expression_group",
        cols = c("Low"="green","Medium"="blue","High"="red"), pt.size=1)

###############################
# 7. Cell Cycle Distribution by Ly6a Group
###############################
ly6a_present_cell$Phase <- seurat_obj$Phase[colnames(ly6a_present_cell)]
cellcycle_summary <- ly6a_present_cell@meta.data %>%
  group_by(expression_group, Phase) %>%
  summarise(count=n()) %>%
  group_by(expression_group) %>%
  mutate(percent = count / sum(count) * 100)

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

###############################
# 8. Differential Expression: Ly6a High vs Low
###############################
hl_cells <- subset(ly6a_present_cell, subset = expression_group %in% c("Low","High"))
Idents(hl_cells) <- hl_cells$expression_group

markers_high_vs_low <- FindMarkers(hl_cells, ident.1="High", ident.2="Low", logfc.threshold=0.25, min.pct=0.1)
head(markers_high_vs_low)

###############################
# 9. GSEA / GO Analysis
###############################
gene_list <- markers_high_vs_low$avg_log2FC
names(gene_list) <- rownames(markers_high_vs_low)
gene_list <- sort(gene_list, decreasing = TRUE)

gene.df <- bitr(names(gene_list), fromType="SYMBOL", toType="ENTREZID", OrgDb=org.Mm.eg.db)
gene_list <- gene_list[gene.df$SYMBOL]
names(gene_list) <- gene.df$ENTREZID

gsea_go <- gseGO(geneList=gene_list, OrgDb=org.Mm.eg.db, ont="BP", keyType="ENTREZID",
                 minGSSize=10, maxGSSize=500, pvalueCutoff=0.05, verbose=TRUE)

# Top 15 GO terms
top_terms <- as.data.frame(gsea_go) %>%
  arrange(p.adjust) %>% head(15)

# Bar plot: NES and -log10(p.adjust)
ggplot(top_terms, aes(x=reorder(Description, -p.adjust), y=-log10(p.adjust), fill=NES)) +
  geom_bar(stat="identity") +
  coord_flip() +
  scale_fill_gradient2(low="blue", mid="white", high="red", midpoint=0) +
  labs(title="Top 15 Enriched GO Terms (Ly6a High vs Low)",
       x="GO Term", y="-log10(adjusted p-value)") +
  theme_minimal()



