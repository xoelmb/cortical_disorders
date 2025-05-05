library(dplyr)
library(Seurat)
source("./sc.fun.R")
source("./preprocess.fun.R")


## Manually curate the data
rawDir <- "/home/sm2726/project2/ASD_iPSC/cellranger_counts/"
sps <- c("D8_Ctrl_290", "D8_Ctrl_311", "D8_Ctrl_317", "D8_ASD_375", "D8_ASD_384", "D8_ASD_494")
ctx_data <- lapply(sps, function(sp) {
	ctx <- Read10X(data.dir = setNames(paste0(rawDir, sp, "/outs/filtered_feature_bc_matrix/"), sp))
	print(paste0("Finish reading sample: ", sp))
	print(dim(ctx))
	return(ctx)
	}) %>%
	do.call(cbind, .)
saveRDS(ctx_data, file = "./data/ASD_iPSC_counts.rds")



## Read the downsampled data
ctx_ds <- Read10X(data.dir = paste0(rawDir, "CtrlASD_D8", "/outs/count/filtered_feature_bc_matrix/"))
ctx_batch_idx <- extract_field(colnames(ctx_ds), 2, "-")
ctx_batch <- c(`1` = "D8_ASD_375", `2` = "D8_ASD_384", `3` = "D8_ASD_494", `4` = "D8_Ctrl_290", `5` = "D8_Ctrl_311", `6` = "D8_Ctrl_317")[ctx_batch_idx]
colnames(ctx_ds) <- gsub("-[0-9]", "-1", colnames(ctx_ds)) %>%
				paste0(ctx_batch, "_", .)
saveRDS(ctx_ds, file = "./data/ASD_iPSC_counts_aggr_default.rds")



## Generate a combined seurat object
ctx_data <- readRDS(file = "./data/ASD_iPSC_counts.rds")
ctx_ds <- readRDS(file = "./data/ASD_iPSC_counts_aggr_default.rds")
seu <- seu_prepare(counts = ctx_data, min.cells = 0, nfeatures = 2500, hvg.method = NULL, assay = "RNA")
seu[["DS"]] <- CreateAssayObject(counts = ctx_ds[, colnames(ctx_data)])


## Add down-sample assay
DefaultAssay(seu) <- "DS"
seu <- NormalizeData(seu, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
DefaultAssay(seu) <- "RNA"


## Set some metadata
seu$cell_origin <- extract_field(colnames(seu), "rm_end", "_")
seu$time <- extract_field(seu$cell_origin, 1, "_")
seu$condition <- extract_field(seu$cell_origin, 2, "_")
seu$sampleid <- extract_field(seu$cell_origin, 3, "_")


## Add doublet scores 
sp_use <- levels(as.factor(seu$cell_origin))
dbmeta <- lapply(sp_use, function(sp) readRDS(file = paste0("./data/scrublet/", sp, "_scrublet_meta.rds"))) %>%
		do.call(rbind, .)
seu@meta.data <- cbind(seu@meta.data, dbmeta[colnames(seu), ])


## Add cell cycle scores (Use DS assay, to avoid batch effects)
DefaultAssay(seu) <- "RNA"
s.genes <- cc.genes$s.genes
g2m.genes <- cc.genes$g2m.genes
cc_scores <- SplitObject(seu, split.by = "cell_origin") %>%
				lapply(., function(xx) {
					xx <- CellCycleScoring(xx, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
					score_data <- xx@meta.data[, c("S.Score", "G2M.Score", "Phase")]
					return(score_data)
					}) %>%
				setNames(., NULL) %>%
				do.call(rbind, .)
seu@meta.data <- cbind(seu@meta.data, cc_scores[colnames(seu), , drop = FALSE])


saveRDS(seu, file = "./data/ASD_iPSC_object.rds")
print(dim(seu))

## Remove doublets
subseu <- seu[, seu$scrublet_assign == 0]
print(dim(subseu))
saveRDS(subseu, file = "./data/ASD_iPSC_object_rmdoublets.rds")



##----------------------------------------------------------------
## Try different set of HVGs

## 1. 
## Set Variable features
hvg_ctrl <- subset(seu, condition %in% "Ctrl") %>%
				SplitObject(., split.by = "cell_origin") %>%
				lapply(., function(x) FindVariableFeatures(x, nfeatures = 2000)) %>%
				SelectIntegrationFeatures(., nfeatures = 2000)
hvg_asd <- subset(seu, condition %in% "ASD") %>%
				SplitObject(., split.by = "cell_origin") %>%
				lapply(., function(x) FindVariableFeatures(x, nfeatures = 2000)) %>%
				SelectIntegrationFeatures(., nfeatures = 2000)
hvg <- union(hvg_ctrl, hvg_asd)
saveRDS(hvg, file = "./data/ASD_iPSC_HVG_set1.rds")


## 2. 
## Set Variable features
hvg <- SplitObject(seu, split.by = "cell_origin") %>%
				lapply(., function(x) FindVariableFeatures(x, nfeatures = 2000)) %>%
				lapply(., function(x) VariableFeatures(x)) %>%
				Reduce("union", .)
saveRDS(hvg, file = "./data/ASD_iPSC_HVG_set2.rds")



## 3. 
## Set Variable features
hvg_ctrl <- subset(seu, condition %in% "Ctrl") %>%
				SplitObject(., split.by = "cell_origin") %>%
				lapply(., function(x) FindVariableFeatures(x, nfeatures = 2500)) %>%
				SelectIntegrationFeatures(., nfeatures = 2000)
hvg_asd <- subset(seu, condition %in% "ASD") %>%
				SplitObject(., split.by = "cell_origin") %>%
				lapply(., function(x) FindVariableFeatures(x, nfeatures = 2500)) %>%
				SelectIntegrationFeatures(., nfeatures = 2000)
hvg_single <- SplitObject(seu, split.by = "cell_origin") %>%
				lapply(., function(x) FindVariableFeatures(x, nfeatures = 1000)) %>%
				lapply(., function(x) VariableFeatures(x)) %>%
				Reduce("union", .)
hvg <- union(hvg_ctrl, hvg_asd) %>% union(., hvg_single)
saveRDS(hvg, file = "./data/ASD_iPSC_HVG_set3.rds")





