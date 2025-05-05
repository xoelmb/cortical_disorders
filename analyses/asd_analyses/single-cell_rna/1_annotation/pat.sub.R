## Subclustering of cluster 10, a putative patterning center cluster
library(Seurat)
library(dplyr)
library(ggplot2)
source("./sc.fun.R")

file_prefix <- "PAT_sub"

### subset to cluster 10 & reprocess
ipsc <- readRDS(file = "./data/Inte_CellCycle_SeuScore_Seurat_t2_v1.seurat.rds")
ipsc$RNA_snn_res.1.2 <- as.character(ipsc$RNA_snn_res.1.2)
pat <- subset(ipsc, RNA_snn_res.1.2 %in% "10")
pat <- FindNeighbors(pat, dims = 1:30, reduction = "pca") %>%
                    FindClusters(., resolution = 1, algorithm = 3)
set.seed(0)
pat <- RunUMAP(pat, dims = 1:30, reduction = "pca", umap.method = "umap-learn", metric = "correlation")


DimFig(pat, group.by = c("RNA_snn_res.1", "condition", "cell_origin"), raster = FALSE, file_name = paste0(file_prefix, "_meta"), label = TRUE, pt.size = 0.5, plot.scale = 1, cols = "glasbey", na.value = NA)
FeatureFig(pat, features = c("FGF8", "ZIC1", "IRX3", "WNT4", "PAX3"), raster = FALSE, file_name = paste0(file_prefix, "_SEL"))

ipsc$RNA_snn_res.1.2[colnames(ipsc) %in% colnames(pat)[pat$RNA_snn_res.1 %in% "3"]] <- "PAT"
DimFig(ipsc, group.by = c("RNA_snn_res.1.2", "condition", "cell_origin"), raster = FALSE, file_name = paste0(file_prefix, "_full"), label = TRUE, pt.size = 0.1, plot.scale = 2, cols = "glasbey", na.value = NA)

pat$identity <- ifelse(pat$RNA_snn_res.1 %in% "3", "FGF17-like", "10")
saveRDS(pat, file = "./data/Inte_CellCycle_SeuScore_Seurat_t2_v1.organizer.seurat.rds")

