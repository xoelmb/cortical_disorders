library(Seurat)
library(dplyr)
source("./mac.fun.R")


##--------------------------------------------------------------------------------
## Organize the macaque data for integration
if (FALSE) {
## load("/gpfs/gibbs/pi/sestan.ycga/sm2726/MonkeyFinal/MF1/overview/load_files/Macaque.developing.seurat.Rdata")
mac <- readRDS(file = "~/myshare2/DataUpload/Micali_Ma_Li_scRNAseq/meta/data/Macaque.dev.seurat.rds")
## Merge subtypes for visualization
clsmerge <- read.table("./data/macaque.cluster.merge.txt", sep = "\t", header = TRUE)
mac$subtype <- clsmerge$subtype[match(mac$cell_subtype, clsmerge$cluster)]
mac$subtype2 <- clsmerge$subtype2[match(mac$cell_subtype, clsmerge$cluster)]

## Define period 
all_ages <- levels(as.factor(mac$age_group))
age2prd <- sapply(all_ages, function(x) {
	yy <- case_when(x %in% c("E37", "E42-43") ~ "E37-43", 
					x %in% c("E54", "E62-64") ~ "E54-64", 
					x %in% c("E77-78") ~ "E77-78", 
					x %in% c("E93") ~ "E93")
	})
mac$period <- age2prd[as.character(mac$age_group)]
saveRDS(mac, file = "./data/Macaque_data_for_integration.rds")
}


mac <- readRDS(file = "./data/Macaque_data_for_integration.rds")
table(mac$samplename)
## [1] "E110_A1C_RMB683"    "E110_A1C_RMB691"    "E110_CGE_RMB683"   
## [4] "E110_CGE_RMB691"    "E110_DFC_RMB683"    "E110_DFC_RMB691"   
## [7] "E110_Insula_RMB683" "E110_Insula_RMB691" "E110_IPC_RMB683"   
##[10] "E110_IPC_RMB691"    "E110_ITC_RMB683"    "E110_ITC_RMB691"   
##[13] "E110_LGE_RMB683"    "E110_LGE_RMB691"    "E110_M1C_RMB683"   
##[16] "E110_M1C_RMB691"    "E110_MFC_RMB683"    "E110_MFC_RMB691"   
##[19] "E110_MGE_RMB683"    "E110_MGE_RMB691"    "E110_OFC_RMB683"   
##[22] "E110_OFC_RMB691"    "E110_PCC_RMB683"    "E110_PCC_RMB691"   
##[25] "E110_S1C_RMB683"    "E110_S1C_RMB691"    "E110_STC_RMB683"   
##[28] "E110_STC_RMB691"    "E110_V1C_RMB683"    "E110_V1C_RMB691"   
##[31] "E110_VFC_RMB683"    "E110_VFC_RMB691"    "E37_FR_0711"       
##[34] "E37_GE_0711"        "E37_OX_0711"        "E42_FR_1018"       
##[37] "E42_GE_1018"        "E42_MS_1018"        "E42_OX_1018"       
##[40] "E43_FR_0926"        "E43_GE_0926"        "E43_OX_0926"       
##[43] "E54_FR_1231"        "E54_GE_1231"        "E54_MS_1231"       
##[46] "E54_OX_1231"        "E54_TP_1231"        "E62_FR_0926"       
##[49] "E62_GE_0926"        "E62_OX_0926"        "E62_TP_0926"       
##[52] "E64_FR_1211"        "E64_GE_1211"        "E64_MS_1211"       
##[55] "E64_OX_1211"        "E64_TP_1211"        "E77_FR_0501"       
##[58] "E77_GE_1222"        "E77_MS_0501"        "E77_MS_1222"       
##[61] "E77_OX_0501"        "E77_OX_1222"        "E77_TP_1222"       
##[64] "E78_FR_0501"        "E78_GE_0501"        "E78_OX_0501"       
##[67] "E93_A1C_1009"       "E93_CGE_1009"       "E93_DFC_1009"      
##[70] "E93_Insula_1009"    "E93_IPC_1009"       "E93_ITC_1009"      
##[73] "E93_LGE_1009"       "E93_M1C_1009"       "E93_MFC_1009"      
##[76] "E93_MGE_1009"       "E93_OFC_1009"       "E93_PCC_1009"      
##[79] "E93_S1C_1009"       "E93_STC_1009"       "E93_V1C_1009"      
##[82] "E93_VFC_1009" 


##--------------------------------------------------------------------------------
## Select representative FR-OC samples for integration
## each age-region pair is treated as a batch
##--------------------------------------------------------------------------------
sel_sps <- c("E37_FR_0711", "E37_OX_0711", "E42_FR_1018", "E42_OX_1018", "E37_GE_0711", "E42_GE_1018", 
        "E54_FR_1231", "E54_OX_1231", "E54_GE_1231", 
        "E64_FR_1211", "E64_OX_1211", "E64_GE_1211",
        "E77_FR_0501", "E77_OX_0501", "E77_GE_1222")
seu <- subset(mac, samplename %in% sel_sps)
seu_new <- OrthoTransfer(object = seu, from = "mmulatta", to = "hsapiens")
seu$period_region <- paste0(seu$period, "_", seu$broad_region) %>% 
                    gsub("E54-64_FC|E54-64_OC", "E54-64_Dorsal", .) %>%
                    gsub("E77-78_FC|E77-78_OC", "E77-78_Dorsal", .)

hvg <- IdentfiyHVG(object = seu, level1.by = "period_region", level2.by = "samplename") %>%
            OrthoHVGs(hvg = .)
VariableFeatures(seu_new) <- hvg
seu_new$inte.batch2 <- paste0(seu_new$age_group, "_", seu_new$broad_region)
saveRDS(seu_new, file = "./data/Macaque_full-FC-GE-OC_orthogene.rds")



sel_sps <- c("E37_FR_0711", "E37_OX_0711", "E42_FR_1018", "E42_OX_1018", 
        "E54_FR_1231", "E54_OX_1231", 
        "E64_FR_1211", "E64_OX_1211", 
        "E77_FR_0501", "E77_OX_0501")
seu <- subset(mac, samplename %in% sel_sps)
seu_new <- OrthoTransfer(object = seu, from = "mmulatta", to = "hsapiens")
seu$period_region <- paste0(seu$period, "_", seu$broad_region) %>% 
                    gsub("E54-64_FC|E54-64_OC", "E54-64", .) %>%
                    gsub("E77-78_FC|E77-78_OC", "E77-78", .) %>%
                    gsub("E93_FC|E93_OC", "E93", .)

hvg <- IdentfiyHVG(object = seu, level1.by = "period_region", level2.by = "samplename") %>%
            OrthoHVGs(hvg = .)
VariableFeatures(seu_new) <- hvg
seu_new$inte.batch2 <- paste0(seu_new$age_group, "_", seu_new$broad_region)
saveRDS(seu_new, file = "./data/Macaque_full-FC-OC_orthogene.rds")


##--------------------------------------------------------------------------------
## Only NSCs (because the cells are much less, intead of treating each age-region pair as
## a batch, i treated each region as a batch)
##--------------------------------------------------------------------------------

## Only NSCs + PAT + Mesenchymal
sel_sps <- c("E37_FR_0711", "E37_OX_0711", "E42_FR_1018", "E42_OX_1018", 
        "E54_FR_1231", "E54_OX_1231", 
        "E64_FR_1211", "E64_OX_1211",
        "E77_FR_0501", "E77_OX_0501")
seu <- subset(mac, broad_region %in% c("FC", "OC") & 
            subtype2 %in% c("RG early", "RG late", "Mesenchymal", "AntVen NKX2-1", "PC FGF17", "PC RSPO3", "PC SFRP2", "PC TCF7L2", "PC TTR") & samplename %in% sel_sps)
seu_new <- OrthoTransfer(object = seu, from = "mmulatta", to = "hsapiens")
hvg <- SplitObject(seu, split.by = "broad_region") %>%
            lapply(., function(mm) FindVariableFeatures(mm, nfeatures = 2000 + 200)) %>%
            SelectIntegrationFeatures(., nfeatures = 2000) %>%
            OrthoHVGs(hvg = .)
VariableFeatures(seu_new) <- hvg
seu_new$inte.batch2 <- seu_new$broad_region
saveRDS(seu_new, file = "./data/Macaque_NSC-PAT-Mes_orthogene.rds")


## only dorsal RGs
sel_sps <- c("E37_FR_0711", "E37_OX_0711", "E42_FR_1018", "E42_OX_1018", 
        "E54_FR_1231", "E54_OX_1231", 
        "E64_FR_1211", "E64_OX_1211",
        "E77_FR_0501", "E77_OX_0501")
seu <- subset(mac, broad_region %in% c("FC", "OC") & 
            subtype2 %in% c("RG early", "RG late") & samplename %in% sel_sps)
seu_new <- OrthoTransfer(object = seu, from = "mmulatta", to = "hsapiens")
hvg <- SplitObject(seu, split.by = "broad_region") %>%
            lapply(., function(mm) FindVariableFeatures(mm, nfeatures = 2000 + 200)) %>%
            SelectIntegrationFeatures(., nfeatures = 2000) %>%
            OrthoHVGs(hvg = .)
VariableFeatures(seu_new) <- hvg
seu_new$inte.batch2 <- seu_new$broad_region
saveRDS(seu_new, file = "./data/Macaque_NSCdor_orthogene.rds")



sel_sps <- c("E37_GE_0711", "E42_GE_1018", 
        "E54_GE_1231", 
        "E64_GE_1211",
        "E77_GE_1222")
seu <- subset(mac, broad_region %in% c("FC", "GE", "OC") & subtype2 %in% c("GE RG HMGA2", "GE RG NKX2-1", "GE RG late") & samplename %in% sel_sps)
seu_new <- OrthoTransfer(object = seu, from = "mmulatta", to = "hsapiens")
hvg <- FindVariableFeatures(seu, nfeatures = 2000) %>%
            VariableFeatures() %>%
            OrthoHVGs(hvg = .)
VariableFeatures(seu_new) <- hvg
seu_new$inte.batch2 <- seu_new$broad_region
saveRDS(seu_new, file = "./data/Macaque_NSC-GE_orthogene.rds")



sel_sps <- c("E37_FR_0711", "E37_OX_0711", "E37_GE_0711")
seu <- subset(mac, broad_region %in% c("FC", "OC") & subtype2 %in% c("AntVen NKX2-1", "PC FGF17", "PC RSPO3", "PC SFRP2", "PC TCF7L2", "PC TTR") & samplename %in% sel_sps)
seu_new <- OrthoTransfer(object = seu, from = "mmulatta", to = "hsapiens")
hvg <- FindVariableFeatures(seu, nfeatures = 2500) %>%
            VariableFeatures() %>%
            OrthoHVGs(hvg = .)
VariableFeatures(seu_new) <- hvg
seu_new$inte.batch2 <- seu_new$age_group
saveRDS(seu_new, file = "./data/Macaque_NSC-PAT_orthogene.rds")




