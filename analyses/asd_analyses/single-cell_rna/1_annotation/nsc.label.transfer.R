## Do label transfer in the integration between macaque and 
library(Seurat)
library(dplyr)
library(ggplot2)
library(parallel)
library(foreach)
source("./sc.fun.R")
source("./preprocess.fun.R")

mac_reg <- "NSCdor"

###-----------------------------------------------------------------
## Perform the custom label transfer (transfer early and late NSCs)
## The default seurat label transfer doesn't work well
seu <- readRDS(file = paste0("./data/MI_NSC_CorrCyc_v1_", mac_reg, ".seurat.rds"))


## subset the RG early to have the same size as RG-late
all_mac_ctps <- levels(as.factor(seu$subtype2))
mac_size <- table(seu$subtype2) %>% sort()
set.seed(0)
rm_cells <- lapply(names(mac_size)[2:length(mac_size)], function(cls){
        cells <- colnames(seu)[seu$subtype2 %in% cls]
        rm_cells <- sample(cells, length(cells)-mac_size[1], replace = FALSE)
    }) %>%
        unlist()
seu <- seu[, setdiff(colnames(seu), rm_cells)]


## For each iteration, subset to 90% ref cells and do label transfer
set.seed(42)
cell.list <- list()
for (ii in 1:100){
    sub_rm <- lapply(names(mac_size), function(cls){
        cells <- colnames(seu)[seu$subtype2 %in% cls]
        rm_cells <- sample(cells, ceiling(length(cells) * 0.1), replace = FALSE)
    }) %>%
        unlist()
    cell.list[[ii]] <- setdiff(colnames(seu), sub_rm)
}

predictions <- LabelTransfer.parallel(object = seu, subset.list = cell.list, reduction = "pca", dims.use = 30, transfer_cols = "subtype2", k = 20, nreps = 100, nCores = 4)
saveRDS(predictions, file = paste0("./data/MI_NSC_CorrCyc_v1_label-transfer_", mac_reg, ".predictions.rds"))


###-----------------------------------------------------------------
## Summarize the predictions
predictions <- readRDS(file = paste0("./data/MI_NSC_CorrCyc_v1_label-transfer_", mac_reg, ".predictions.rds"))
predictions <- predictions %>%
                filter(grepl("^D8", cell))

pred.sum <- predictions %>%
                filter(replicate != 100) %>% ## to have odd number, so the most abundant is either early or late
                group_by(cell, subtype2_new) %>%
                summarize(nhits = n()) %>%
                ungroup() %>%
                group_by(cell) %>%
                filter(nhits == max(nhits)) %>%
                ungroup()
saveRDS(pred.sum, file = paste0("./data/MI_NSC_CorrCyc_v1_label-transfer_", mac_reg, ".pred.sum.rds"))


###-----------------------------------------------------------------
## Obtain the cluster-level predictions
ipsc <- readRDS(file = "./data/Inte_CellCycle_SeuScore_Seurat_t2_v1.seurat.rds")
ipsc$RNA_snn_res.1.2 <- as.character(ipsc$RNA_snn_res.1.2)
ipsc <- subset(ipsc, RNA_snn_res.1.2 %in% as.character(c(0,1,3,4,11,5,2,9)))
ipsc <- FindNeighbors(ipsc, dims = 1:30, reduction = "pca") %>%
                    FindClusters(., resolution = 2, algorithm = 3)
ipsc$identity <- paste0("nsc", as.character(ipsc$seurat_clusters))



pred.sum <- readRDS(file = paste0("./data/MI_NSC_CorrCyc_v1_label-transfer_", mac_reg, ".pred.sum.rds")) %>%
            filter(cell %in% colnames(ipsc)) %>%
            mutate(identity = ipsc$identity[match(cell, colnames(ipsc))]) %>%
            group_by(identity, subtype2_new) %>%
            summarize(nhits = n()) %>%
            ungroup() %>%
            group_by(identity) %>%
            mutate(prop = nhits/sum(nhits)) %>%
            filter(prop == max(prop)) %>%
            ungroup()
cls2label <- setNames(pred.sum$subtype2_new, pred.sum$identity)
ipsc$anno <- cls2label[ipsc$identity]
saveRDS(ipsc, file = paste0("./data/Inte_CellCycle_SeuScore_Seurat_t2_v1.NSC.seurat.rds"))



###-----------------------------------------------------------------
## Visualization of the labels
mac_reg <- "NSCdor"
pred.sum <- readRDS(file = paste0("./data/MI_NSC_CorrCyc_v1_label-transfer_", mac_reg, ".pred.sum.rds"))
seu <- readRDS(file = paste0("./data/MI_NSC_CorrCyc_v1_", mac_reg, ".seurat.rds"))
seu$anno_cell <- seu$subtype2
seu@meta.data[pred.sum$cell, "anno_cell"] <- pred.sum$subtype2_new
seu$dataset <- ifelse(grepl("^D8", seu$inte.batch), "ASD_iPSC", "Macaque")
seu$anno_cell[(!as.character(seu$RNA_snn_res.1.2) %in% as.character(c(0,1,3,4,11,5,2,9))) &
                (seu$dataset %in% "ASD_iPSC")] <- NA


ipsc <- readRDS(file = paste0("./data/Inte_CellCycle_SeuScore_Seurat_t2_v1.NSC.seurat.rds"))
seu$anno_cluster <- seu$subtype2
seu@meta.data[colnames(ipsc), "anno_cluster"] <- ipsc$anno



file_prefix <- paste0("MI_NSC_CorrCyc_v1_label-transfer_", mac_reg)
p_dim <- DimPlot(object = seu, group.by = c("anno_cluster", "anno_cell"), pt.size = 0.1, raster = FALSE, split.by = "dataset", cols = "glasbey") &
			theme(plot.title = element_text(size = 8, face = "bold"), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())
jpeg(paste0("./report/", file_prefix, ".split.jpeg"), width = 6.5 * 2, height = 6 * 2, unit = "in", res = 300)
print(p_dim)
dev.off()


DimFig(seu, group.by = c("RNA_snn_res.1.2", "anno_cell", "anno_cluster", "subtype2"), raster = FALSE, file_name = paste0(file_prefix, "_meta"), label = TRUE, pt.size = 0.1, plot.scale = 2, cols = "glasbey", na.value = NA)




