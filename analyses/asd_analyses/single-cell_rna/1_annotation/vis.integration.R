## Visualize the integration results
library(Seurat)
library(dplyr)
library(ggplot2)
source("./sc.fun.R")


##------------------------------------------------------------------------------------
## Full integration with all the macaque cell types
mac_reg <- "full-FC-OC"
seuo <- readRDS(file = paste0("./data/MI_Full_CorrCyc_", mac_reg, ".seurat.rds"))
seuo$dataset <- ifelse(grepl("^D8", seuo$inte.batch), "ASD_iPSC", "Macaque")


seuo$macaque_anno <- seuo$cell_subclass
seuo$macaque_anno[seuo$macaque_anno %in% c("enIPC", "Excitatory neurons")] <- "ExNs"
seuo$macaque_anno[seuo$macaque_anno %in% c("inIPC", "Inhibitory neurons")] <- "InNs"
seuo$macaque_anno[seuo$subtype2 %in% "AntVen NKX2-1"] <- "Patterning centers"
seuo$macaque_anno[seuo$subtype2 %in% c("Astro", "OPC-Oligo", "gIPC")] <- "Gliogenesis"
seuo$macaque_anno[seuo$macaque_anno %in% "Early subtypes"] <- NA

ipsc <- readRDS(file = "./data/Inte_CellCycle_SeuScore_Seurat_annotated.rds")
seuo@meta.data[colnames(ipsc), "anno_cell"] <- ipsc$anno_cell
seuo@meta.data[colnames(ipsc), "anno_cluster"] <- ipsc$anno_cluster

##seuo$ipsc_anno <- seuo$RNA_snn_res.1.2

file_prefix <- paste0("Figure_MCInte_Full_", mac_reg)
DimFig(seuo, group.by = c("dataset", "anno_cell", "anno_cluster", "macaque_anno"), raster = FALSE, file_name = paste0(file_prefix, "_meta"), label = TRUE, pt.size = 0.1, plot.scale = 2, cols = "glasbey", na.value = NA)


##------------------------------------------------------------------------------------
## Integration with only macaque NSCs
for (mac_reg in c("NSCdor", "NSCdor-PAT-Mes", "NSC-GE", "NSC-PAT")){
    seuo <- readRDS(file = paste0("./data/MI_NSC_CorrCyc_v1_", mac_reg, ".seurat.rds"))
    seuo$dataset <- ifelse(grepl("^D8", seuo$inte.batch), "ASD_iPSC", "Macaque")
    seuo$macaque_anno <- seuo$subtype2
    seuo$ipsc_anno <- seuo$RNA_snn_res.1.2
    ipsc <- readRDS(file = "./data/Inte_CellCycle_SeuScore_Seurat_annotated.rds")
    seuo@meta.data[colnames(ipsc), "anno_cell"] <- ipsc$anno_cell
    seuo@meta.data[colnames(ipsc), "anno_cluster"] <- ipsc$anno_cluster


    file_prefix <- paste0("Figure_MCInte_NSC_", mac_reg)
    DimFig(seuo, group.by = c("dataset", "anno_cell", "anno_cluster", "macaque_anno"), raster = FALSE, file_name = paste0(file_prefix, "_meta"), label = TRUE, pt.size = 0.1, plot.scale = 2, cols = "glasbey", na.value = NA)
}













