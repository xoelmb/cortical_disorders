## Combine the annotations
library(Seurat)
library(dplyr)
source("./sc.fun.R")

## Manually annotate the cell clusters
ipsc <- readRDS(file = "./data/Inte_CellCycle_SeuScore_Seurat_t2_v1.seurat.rds")
ipsc$RNA_snn_res.1.2 <- as.character(ipsc$RNA_snn_res.1.2)

ipsc$anno_cell <- NA
ipsc$anno_cell[ipsc$RNA_snn_res.1.2 %in% as.character(c(6,7,12,14,15,16,17))] <- "High-mito"
ipsc$anno_cell[ipsc$RNA_snn_res.1.2 %in% as.character(c(10,13,19,8))] <- "Mes prog"
ipsc$anno_cell[ipsc$RNA_snn_res.1.2 %in% as.character(c(20))] <- "Neral crest/Mes"
ipsc$anno_cell[ipsc$RNA_snn_res.1.2 %in% as.character(c(18))] <- "Neuron"
ipsc$anno_cell[ipsc$RNA_snn_res.1.2 %in% as.character(c(0,1,3,4,11,5,2,9))] <- "NSC"


## Add organizer identities
pat <- readRDS(file = "./data/Inte_CellCycle_SeuScore_Seurat_t2_v1.organizer.seurat.rds")
ipsc$anno_cell[colnames(ipsc) %in% colnames(pat)[pat$identity %in% "FGF17-like"]] <- "FGF17-like"


## Add NSC identities (cell-level)
pred.sum <- readRDS(file = paste0("./data/MI_NSC_CorrCyc_v1_label-transfer_", "NSCdor", ".pred.sum.rds"))
pred.sum <- filter(pred.sum, cell %in% colnames(ipsc)[ipsc$anno_cell %in% "NSC"])
ipsc@meta.data[pred.sum$cell, "anno_cell"] <- pred.sum$subtype2_new


nsc <- readRDS(file = paste0("./data/Inte_CellCycle_SeuScore_Seurat_t2_v1.NSC.seurat.rds"))
ipsc$anno_cluster <- ipsc$anno_cell
ipsc@meta.data[colnames(nsc), "anno_cluster"] <- nsc$anno


file_prefix <- "Anno_Final"
DimFig(ipsc, group.by = c("anno_cell", "anno_cluster", "condition", "cell_origin"), raster = FALSE, file_name = paste0(file_prefix, "_full"), label = TRUE, pt.size = 0.1, plot.scale = 2, cols = "glasbey", na.value = NA)

saveRDS(ipsc, file = "./data/Inte_CellCycle_SeuScore_Seurat_annotated.rds")



##---------------------------------------------------------------------------
## Sankey plot showing the correspondance between clusters & macaque cortical NSCs
library(Seurat)
library(tidyverse)
source("./sc.fun.R")
source("./preprocess.fun.R")

ipsc <- readRDS(file = "./data/Inte_CellCycle_SeuScore_Seurat_t2_v1.seurat.rds")
ipsc$RNA_snn_res.1.2 <- as.character(ipsc$RNA_snn_res.1.2)

## UMAP showing the raw clusters
p <- DimPlot(ipsc, group.by = "RNA_snn_res.1.2", label.size = 4, raster = TRUE, label = TRUE, pt.size = 1.5, raster.dpi = c(512, 512)) +
		theme(legend.position = "right")
pdf(paste0("./report/Anno_Final_seurat_clusters.pdf"), width = 6.5, height = 5)
print(p)
dev.off()



pred.sum <- readRDS(file = paste0("./data/MI_NSC_CorrCyc_v1_label-transfer_", "NSCdor", ".pred.sum.rds"))
nsc_cells <- colnames(ipsc)[ipsc$RNA_snn_res.1.2 %in% as.character(c(0,1,3,4,11,5,2,9))]

pmeta <- pred.sum %>%
			filter(cell %in% nsc_cells) %>%
			mutate(cluster = ipsc$RNA_snn_res.1.2[match(cell, colnames(ipsc))]) %>%
			select(-nhits) %>%
			as.data.frame()


p <- plot_ggsankey(meta = pmeta, 
                    x_col = "cluster", 
                    x_ord = as.character(c(0,1,3,4,11,5,2,9)), 
                    next_col = "subtype2_new", 
                    next_ord = c("RG early", "RG late"), 
                    type = "sankey")

pdf(paste0("./report/Anno_Final_NSC_macaque_v1_sankey.pdf"), width = 6, height = 6)
print(p)
dev.off()




##-----------------------------------------------------------------------------
## Plot some common marker genes
library(Seurat)
library(tidyverse)
source("./sc.fun.R")
source("./preprocess.fun.R")

ipsc <- readRDS(file = "./data/Inte_CellCycle_SeuScore_Seurat_t2_v1.seurat.rds")
ipsc$RNA_snn_res.1.2 <- as.character(ipsc$RNA_snn_res.1.2)
Idents(ipsc) <- "RNA_snn_res.1.2"


asd <- subset(ipsc, condition %in% "ASD")
cls_ord <- as.character(c(0,1,3,4,11,5,2,9,18,10,13,19,8,20,6,7,12,14,15,16,17))
custom_genes <- c("FGF17", "FGF8", "FGF18", "SP8", "ZIC1", "ZIC3", "ZIC4", "RSPO3", "RSPO2", "ARX", "WNT2B", "WNT8B", "LMX1A", "TTR", "FOXJ1", "BMP4", "SHH", "NKX2-1", "SOX2", "NES", "VIM", "PAX6", "HMGA2", "CCND1", "STMN2", "SLC1A3", "HOPX", "FABP7", "PTN", "NRG1", "CRYAB", "CXCL12", "MKI67", "TOP2A", "HMGB2", "CDK1", "NUSAP1", "PCNA", "MCM5", "TYMS", "FEN1", "MCM2", "FN1", "LUM", "FOXD3", "PLP1", "EOMES", "NEUROG1", "ASCL1", "NHLH1", "PPP1R17", "DCX", "SLA", "NEUROD4", "NEUROD2", "NEUROD6", "SATB2", "TLE4", "SOX5", "FEZF2", "FOXP2", "DLX1", "DLX2", "GAD1", "GAD2", "MEIS2", "FOXP1", "EGFR", "OLIG2", "PDGFRA", "MOBP", "MBP", "AQP4", "GFAP", "PTPRC", "C1QB", "C1QC", "CLDN5", "FLT1", "ACTA2", "CEMIP", "MT-ND1", "MT-ND2", "MT-CO1", "MT-CO2")
p2 <- DotPlot(ipsc, features = custom_genes, cols = c("lightgrey", "red"), dot.scale = 5, scale.by = "size", dot.min = 0.025) +
			theme_classic()+
			RotatedAxis() + 
			scale_y_discrete(limits = rev(cls_ord)) +
			labs(title = paste0("Custom markers")) +
			theme(axis.ticks = element_line(size = 0.2), axis.line = element_line(size = 0.2),
				axis.text.y = element_text(size = rel(0.8)), axis.text.x = element_text(size = rel(0.8)), axis.title = element_blank(), legend.position = "bottom")
pdf(paste0("./report/Anno_Final_custom_markers.dotplot.pdf"), width = 12, height = 5.5, useDingbats = FALSE)
print(p2)
dev.off()














