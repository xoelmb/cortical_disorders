args <- commandArgs(trailingOnly = TRUE)
library(Seurat)
library(dplyr)
library(ggplot2)
 
mac_reg <- c("NSCdor", "NSCdor-PAT-Mes", "NSC-GE", "NSC-PAT")[as.numeric(args[1])]
umapFile <- paste0("./data/MI_NSC_CorrCyc_v1_", mac_reg, ".NSCslim.rds")
if (!file.exists(umapFile)){
    seuo <- readRDS(file = paste0("./data/MI_NSC_CorrCyc_v1_", mac_reg, ".seurat.rds"))
    seuo$dataset <- ifelse(grepl("^D8", seuo$inte.batch), "ASD_iPSC", "Macaque")

    ## Subset to ASD_iPSC data & specific macaque subtypes
    sel_stps <- c("GE RG NKX2-1", "GE RG HMGA2", "GE RG early", "GE RG late", "Mesenchymal", "RG early", "RG late")
    sel_cls <- as.character(c(0, 1, 2, 3, 4, 5, 8, 9, 10, 11, 13, 19, 20))
    cells <- c(colnames(seuo)[seuo$dataset %in% "ASD_iPSC" & seuo$RNA_snn_res.1.2 %in% sel_cls],
                colnames(seuo)[seuo$dataset %in% "Macaque" & seuo$subtype2 %in% sel_stps])


    subseu <- seuo[, cells]
    subseu <- RunUMAP(subseu, dims = 1:30, umap.method = "umap-learn", metric = "correlation")
    saveRDS(subseu, file = umapFile)
}

subseu <- readRDS(file = umapFile)

source("./sc.fun.R")
file_prefix <- paste0("MC_Inte_NSC_CorrectCycle_v1_", mac_reg, "_NSC")
DimFig(subseu, group.by = c("RNA_snn_res.1.2", "subtype2", "subtype", "cell_subclass", "dataset", "age_group", "broad_region"), raster = FALSE, file_name = paste0(file_prefix, "_meta"), label = TRUE, pt.size = 0.1, plot.scale = 2, cols = "glasbey", na.value = NA)


##dsq --job-file nsc.inte.sub.task.txt --batch-file nsc.inte.sub.sh -c 4 -p day --mem-per-cpu=20G -J umap --max-jobs 4 -N 1 -t 1-00:00:00



