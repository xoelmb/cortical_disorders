library(dplyr)
library(Seurat)
library(Matrix)
library(ggplot2)
source("./sc.fun.R")
source("preprocess.fun.R")


seu <- readRDS(file = paste0("./data/Inte_CellCycle_SeuScore_Seurat_t2_v", 1, ".slim.rds"))

oldseu <- readRDS("../preprocess2/data/Inte_CellCycle_SeuScore_Seurat_v1.slim.rds")


seu$oldcluster <- NA
sh_cells <- intersect(colnames(oldseu), colnames(seu))
seu@meta.data[sh_cells, "oldcluster"] <- as.character(oldseu@meta.data[sh_cells, "seurat_clusters"])

file_prefix <- paste0("Inte_CellCycle_SeuScore_Seurat_t2_v", 1)
DimFig(seu, group.by = "oldcluster", raster = FALSE, file_name = paste0(file_prefix, "_old_cluster"), label = TRUE, pt.size = 0.1, plot.scale = 1.2)


