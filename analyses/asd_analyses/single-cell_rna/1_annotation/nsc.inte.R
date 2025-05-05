## Pairwise integration with the macaque data from different regions
args <- commandArgs(trailingOnly = TRUE)
library(Seurat)
library(dplyr)
library(ggplot2)
source("./sc.fun.R")
source("./preprocess.fun.R")


library(future)
plan("multicore", workers = 8)
options(future.globals.maxSize = 40*1000*1024^2)


## Subset the data of interest for the data integration
mac_reg <- c("NSCdor", "NSCdor-PAT-Mes", "NSC-GE", "NSC-PAT")[as.numeric(args[1])]
mac <- readRDS(file = paste0("./data/Macaque_", mac_reg, "_orthogene.rds"))
hvg1 <- VariableFeatures(mac)


## Load the ASD data
ipsc <- readRDS(file = "./data/Inte_CellCycle_SeuScore_Seurat_t2_v1.seurat.rds")
ipsc$inte.batch <- ipsc$cell_origin
ipsc$inte.batch2 <- ipsc$cell_origin
hvg2 <- readRDS(file = "./data/ASD_iPSC_HVG_set1.rds")


## Combined the two datasets using the shared genes
sh_genes <- intersect(rownames(mac), rownames(ipsc))
cbn <- merge(x = mac[sh_genes, ], y = ipsc[sh_genes, ])
hvg <- union(hvg1, hvg2) %>%
			intersect(., sh_genes)

seu_list <- SplitObject(cbn, split.by = "inte.batch2")

rm(cbn, mac, ipsc)
seu <- Integratelist.seurat.cellcycle(obj.list = seu_list, hvg = hvg, file_name = paste0("MI_NSC_CorrCyc_v1_", mac_reg), input_dir = "./data/", inte.dims = 1:30, cluster.dims = 1:30, reference = NULL, do.cluster = FALSE)


###-----------------------------------------------------------------
## Visualize the integration results
seuo <- readRDS(file = paste0("./data/MI_NSC_CorrCyc_v1_", mac_reg, ".seurat.rds"))
seuo$dataset <- ifelse(grepl("^D8", seuo$inte.batch), "ASD_iPSC", "Macaque")


file_prefix <- paste0("MI_NSC_CorrCyc_v1_", mac_reg)
DimFig(seuo, group.by = c("RNA_snn_res.1.2", "subtype", "subtype2", "cell_subclass", "dataset", "age_group", "broad_region"), raster = FALSE, file_name = paste0(file_prefix, "_meta"), label = TRUE, pt.size = 0.1, plot.scale = 2, cols = "glasbey", na.value = NA)

##dsq --job-file nsc.inte.task.txt --batch-file nsc.inte.sh -c 8 -p day --mem-per-cpu=40G -J nsc --max-jobs 4 -N 1 -t 1-00:00:00



