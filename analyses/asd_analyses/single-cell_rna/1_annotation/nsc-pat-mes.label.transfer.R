## Do label transfer in the integration between macaque and 
library(Seurat)
library(dplyr)
library(ggplot2)
library(parallel)
library(foreach)
source("./sc.fun.R")
source("./preprocess.fun.R")

mac_reg <- "NSCdor-PAT-Mes"

###-----------------------------------------------------------------
## Perform the custom label transfer (transfer early and late NSCs)
## The default seurat label transfer doesn't work well
seu <- readRDS(file = paste0("./data/MI_NSC_CorrCyc_v1_", mac_reg, ".seurat.rds"))


## subset the RG early to have the same size as RG-late
all_mac_ctps <- levels(as.factor(seu$subtype2))
mac_size <- table(seu$subtype2) %>% sort()
set.seed(0)
rm_cells <- lapply(all_mac_ctps, function(cls){
        cells <- colnames(seu)[seu$subtype2 %in% cls]
        if (length(cells) > 8000){
            rm_cells <- sample(cells, length(cells)-8000, replace = FALSE)
        } else {
            rm_cells <- c()
        }
        return(rm_cells)
    }) %>%
        unlist()
seu2 <- seu[, setdiff(colnames(seu), rm_cells)]


## For each iteration, subset to 90% ref cells and do label transfer
set.seed(42)
cell.list <- list()
for (ii in 1:100){
    sub_rm <- lapply(names(mac_size), function(cls){
        cells <- colnames(seu2)[seu2$subtype2 %in% cls]
        rm_cells <- sample(cells, ceiling(length(cells) * 0.1), replace = FALSE)
    }) %>%
        unlist()
    cell.list[[ii]] <- setdiff(colnames(seu2), sub_rm)
}

predictions <- LabelTransfer.parallel(object = seu2, subset.list = cell.list, reduction = "pca", dims.use = 30, transfer_cols = "subtype2", k = 20, nreps = 100, nCores = 20)
saveRDS(predictions, file = paste0("./data/MI_NSC_CorrCyc_v1_label-transfer_", mac_reg, ".predictions.rds"))


###-----------------------------------------------------------------
## Summarize the predictions
predictions <- readRDS(file = paste0("./data/MI_NSC_CorrCyc_v1_label-transfer_", mac_reg, ".predictions.rds"))
predictions <- predictions %>%
                filter(grepl("^D8", cell))

pred.sum <- predictions %>%
                group_by(cell, subtype2_new) %>%
                summarize(nhits = n()) %>%
                ungroup() %>%
                group_by(cell) %>%
                filter(nhits == max(nhits)) %>%
                ungroup()
saveRDS(pred.sum, file = paste0("./data/MI_NSC_CorrCyc_v1_label-transfer_", mac_reg, ".pred.sum.rds"))



##---------------------------------------------------------------------------
## Sankey plot showing the correspondance between clusters & macaque cortical NSCs
library(Seurat)
library(tidyverse)
source("./sc.fun.R")
source("./preprocess.fun.R")

ipsc <- readRDS(file = "./data/Inte_CellCycle_SeuScore_Seurat_t2_v1.seurat.rds")
ipsc$RNA_snn_res.1.2 <- as.character(ipsc$RNA_snn_res.1.2)

pred.sum <- readRDS(file = paste0("./data/MI_NSC_CorrCyc_v1_label-transfer_", "NSCdor-PAT-Mes", ".pred.sum.rds")) %>%
            group_by(cell) %>%
            mutate(nids = length(unique(subtype2_new))) %>%
            ungroup()

sel_cells <- colnames(ipsc)[!ipsc$RNA_snn_res.1.2 %in% as.character(c(18, 6,7,12,14,15,16,17))]


pred.sum$ident <- pred.sum$subtype2_new
pred.sum$ident[pred.sum$nids > 1] <- "ambiguous"

pmeta <- pred.sum %>%
            filter(cell %in% sel_cells) %>%
            group_by(cell) %>%
            summarize(ident = unique(ident)) %>%
            ungroup() %>%
            mutate(cluster = ipsc$RNA_snn_res.1.2[match(cell, colnames(ipsc))]) %>%
            as.data.frame()


p <- plot_ggsankey(meta = pmeta, 
                    x_col = "cluster", 
                    x_ord = as.character(c(0,1,3,4,11,5,2,9,8,10,13,19,20)), 
                    next_col = "ident", 
                    next_ord = c("PC FGF17", "AntVen NKX2-1", "PC SFRP2", "PC TCF7L2", "PC RSPO3", "PC TTR", "RG early", "RG late", "Mesenchymal", "ambiguous"), 
                    type = "sankey")

pdf(paste0("./report/Anno_Final_NSC-PAT-Mes_macaque_v1_sankey.pdf"), width = 6, height = 10)
print(p)
dev.off()

