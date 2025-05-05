args <- commandArgs(trailingOnly = TRUE)
library(dplyr)
library(Seurat)
library(Matrix)
library(ggplot2)
source("./sc.fun.R")
source("preprocess.fun.R")


library(future)
plan("multisession", workers = 6)
options(future.globals.maxSize = 80*1000*1024^2)


seu <- readRDS(file = "./data/ASD_iPSC_object_rmdoublets.rds")

hvgversion <- args[1]

## Convert object to list
seu_list <- SplitObject(seu, split.by = "cell_origin")
hvg <- readRDS(file = paste0("./data/ASD_iPSC_HVG_set", hvgversion, ".rds"))


## Integration
seu <- Integratelist.seurat.cellcycle(obj.list = seu_list, hvg = hvg, file_name = paste0("Inte_CellCycle_SeuScore_Seurat_t2_v", hvgversion), input_dir = "./data/", inte.dims = 1:30, cluster.dims = 1:30, reference = NULL, do.cluster = TRUE)


## Visualization
seu <- readRDS(file = paste0("./data/Inte_CellCycle_SeuScore_Seurat_t2_v", hvgversion, ".seurat.rds"))



file_prefix <- paste0("Inte_CellCycle_SeuScore_Seurat_t2_v", hvgversion)
DimFig(seu, group.by = c("Phase", "seurat_clusters", "cell_origin", "condition"), raster = FALSE, file_name = paste0(file_prefix, "_meta"), label = TRUE, pt.size = 0.1, plot.scale = 1.2)
features <- c("IRX3", "WNT4", "WNT1", "PAX3", "PCDH8", "RSPO3", "FGF8", "SP8", "ZIC1", "ZIC3", "ZIC4", "SOX2", "NES", "MKI67", "PCNA", "FGF17", "RSPO3", "TRIM71", "EOMES", "NEUROD2", "DCX", "DLX1", "OLIG2", "PLP1", "FOXD3", "EGFR") %>%
		intersect(., rownames(seu))
FeatureFig(seu, features = features, file_name = paste0(file_prefix, "_expr"), plot.scale = 1, pt.size = 0.1)
FeatureFig(seu, features = c("S.Score", "G2M.Score"), file_name = paste0(file_prefix, "_phase"), plot.scale = 1)
FeatureFig(seu, features = c("nCount_RNA", "nFeature_RNA", "percent.mt", "percent.ribo", "scrublet_score", "scrublet_assign"), file_name = paste0(file_prefix, "_quality"), plot.scale = 1)



## Dotplot
features <- c("IRX3", "WNT4", "WNT1", "PAX3", "PCDH8", "RSPO3", "FGF8", "SP8", "ZIC1", "ZIC3", "ZIC4")
p <- DotPlot(seu, features = features, cols = c("lightgrey", "red"), dot.scale = 5, scale.by = "radius", dot.min = 0.025) +
			theme_classic()+
			RotatedAxis() +
			coord_flip() + 
			theme(axis.ticks = element_line(size = 0.2), axis.line = element_line(size = 0.2),
				axis.text.y = element_text(size = rel(0.8)), axis.text.x = element_text(size = rel(0.7)), axis.title = element_blank(), legend.position = "bottom")
pdf(paste0("./report/", file_prefix, ".selected.dotplot.pdf"), width = 6, height = 4, useDingbats = FALSE)
print(p)
dev.off()


















