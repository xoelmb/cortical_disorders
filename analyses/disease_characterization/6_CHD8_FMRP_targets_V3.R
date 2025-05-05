### FMRP CHD8 Targets

library(dplyr)
library(openxlsx)
library(reshape2)
library(ggplot2)
library(biomaRt)

setwd("/Users/gabrielsantpere/Desktop/Neurogenomics/Projectes/CortexMalformation/DATA/DiseaseGenes/FMRP_CHD8/")

chd8_Werling <- read.table("../CHD8_targets.txt", header=F)
chd8_Sugathan <- read.table("Sugathan_et_al_2014.txt", header = TRUE)
chd8_SugathanHC <- read.table("Sugathan_et_al_2014_HC.txt", header = TRUE)
filename <-"41467_2015_BFncomms7404_MOESM1023_ESM.xlsx"
sheets <- openxlsx::getSheetNames(filename)
chd8_Cotney <- lapply(sheets,openxlsx::read.xlsx,xlsxFile=filename)
names(chd8_Cotney) <- sheets

fmrp_Darnell <- read.table("Darnell_2011.txt", header=FALSE)
fmrp_Casingal <- read.table("Casingal_2020.txt", header=FALSE)


convertMouseGeneList <- function(x){
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[,])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

fmrp_Casingal <- convertMouseGeneList(fmrp_Casingal$V1)

fmrp_TOTAL <- unique(c(fmrp_Casingal$HGNC.symbol, fmrp_Darnell$V1))
chd8_TOTAL <- unique( c(chd8_Sugathan$Gene, chd8_SugathanHC$Gene, chd8_Cotney[["hNSC_Chd8_prom"]]$Gene,  
                        chd8_Cotney[["hNSC_specific_Chd8_prom"]]$Gene, chd8_Cotney[["human_brain_Chd8_prom"]]$Gene,
                        chd8_Cotney[["hNSC+human_brain_Chd8_prom"]]$Gene, chd8_Cotney[["hNSC+human+mouse_Chd8_prom"]]$Gene,
                        chd8_Cotney[["hNSC_Chd8_enhancer"]]$Gene, chd8_Cotney[["hNSC+human_brain_Chd8_enhancer"]]$Gene,
                        chd8_Cotney[["human_brain_Chd8_enhancer"]]$Gene, chd8_Cotney[["hNSC+human+mouse_Chd8_enhancer"]]$Gene) )

save(fmrp_TOTAL, chd8_TOTAL, file="../FMRP_CHD8_TOTAL.rda")

load("../FMRP_CHD8_TOTAL.rda")
load("../DiseaseFULL.rda")

#### Aqui la tabla que necesitamos
chd8_table <- data.frame(GENE=disease_full$GENE,
                         Sugathan_2014 = ifelse(disease_full$GENE %in% chd8_Sugathan$Gene, 1,0),
                         Sugathan_2014_HC = ifelse(disease_full$GENE %in% chd8_SugathanHC$Gene, 1,0),
                         hNSC_prom = ifelse(disease_full$GENE %in% chd8_Cotney[["hNSC_Chd8_prom"]]$Gene, 1,0),
                         hNSC_specific_prom = ifelse(disease_full$GENE %in% chd8_Cotney[["hNSC_specific_Chd8_prom"]]$Gene, 1,0),
                         human_brain_prom = ifelse(disease_full$GENE %in% chd8_Cotney[["human_brain_Chd8_prom"]]$Gene, 1,0),
                         hNSC_human_brain_prom = ifelse(disease_full$GENE %in% chd8_Cotney[["hNSC+human_brain_Chd8_prom"]]$Gene, 1,0),
                         hNSC_human_mouse_prom = ifelse(disease_full$GENE %in% chd8_Cotney[["hNSC+human+mouse_Chd8_prom"]]$Gene, 1,0),
                         hNSC_enhancer = ifelse(disease_full$GENE %in% chd8_Cotney[["hNSC_Chd8_enhancer"]]$Gene, 1,0),
                         hNSC_human_brain_enhancer = ifelse(disease_full$GENE %in% chd8_Cotney[["hNSC+human_brain_Chd8_enhancer"]]$Gene, 1,0),
                         human_brain_enhancer = ifelse(disease_full$GENE %in% chd8_Cotney[["human_brain_Chd8_enhancer"]]$Gene, 1,0),
                         hNSC_human_mouse_enhancer = ifelse(disease_full$GENE %in% chd8_Cotney[["hNSC+human+mouse_Chd8_enhancer"]]$Gene, 1,0),
                         FMRP_Casingal = ifelse(disease_full$GENE %in% fmrp_Casingal$HGNC.symbol, 1,0),
                         FMRP_Darnell = ifelse(disease_full$GENE %in% fmrp_Darnell$V1, 1,0)
                         )


disgenet <- read.table("../../DisGENET/curated_gene_disease_associations.tsv", sep="\t", fill=TRUE, quote = "", header = TRUE)
limit <- names(table(disgenet$diseaseId)[which(table(disgenet$diseaseId)>40)])
filterDis <- filter(disgenet, diseaseType == "disease", diseaseId %in% limit)

load("../../MaximumTenporal_231121.rda")
load("../..//CoGAPS_Expression/NicoHsPASS_ExprsLog2.Z_forGabrielHeatmaps.RData")
load("../..//CoGAPS_Expression/NicoHsPASS_ExprsLog2_forGabrielHeatmaps.rdata")
load("../..//CoGAPS_Expression/Maximums.rda")
load("../../../OrderMix_dis_22Nov21.rda")
df2plotDis <- readRDS("Passages2plot.rds")

E <- as.data.frame(NicoHsPASS.ExprsLog2.Z)
E2 <- as.data.frame(NicoHsPASS.ExprsLog2[rownames(E),])

means <- rowMeans(E2)
bottom <- names( head( sort(means), 1000 ) )
top <- names( tail( sort(means), 1000 ) )
oneRPKM <- names( means[means > 1] )
df_bottom <- dfALL[bottom,]
df_bottom$Disease <- "Bottom1000"
df_top <- dfALL[top,]
df_top$Disease <- "Top1000"
df_1RPKM <- dfALL[oneRPKM,]
df_1RPKM$Disease <- ">1 RPKM"

maximums<-as.matrix(t(apply(E, 1, function(x){ (x-min(x))/(max(x)-min(x))   })))
maximums_max2 <- apply(maximums, 1, function(x){ head(which(x==max(x)), 1) })
dfALL <- data.frame(Disease="ALL", Max=maximums_max2)

max_list2 <- c()
for (i in unique(max_list$Disease)){
  g <- disease_full[which(disease_full[,i]==1),]$GENE
  df <- data.frame(Disease=i, Max=maximums_max2[g], Gene=g)
  max_list2 <- rbind(max_list2, df)
}

dfALL$Gene <- rownames(dfALL)
disgenet_EL<-c()
for (i in unique(filterDis$diseaseName)){
  temp <- dfALL[which(rownames(dfALL) %in% filter(filterDis, diseaseName == i)$geneSymbol),]
  temp$Disease <- i
  disgenet_EL <- rbind(disgenet_EL, temp)
}

df_bottom$Gene <- rownames(df_bottom)
df_top$Gene <- rownames(df_top)
df_1RPKM$Gene <- rownames(df_1RPKM)

disgenet_EL <- rbind(disgenet_EL, dfALL, max_list2, df_bottom, df_top, df_1RPKM) 
disgenet_EL$Early <- ifelse( disgenet_EL$Max %in% c(1:4), "P2", "P8" )
disgenet_EL$Early[which(disgenet_EL$Max %in% c(5:8))] <- "P3"
disgenet_EL$Early[which(disgenet_EL$Max %in% c(9:12))] <- "P4"
disgenet_EL$Early[which(disgenet_EL$Max %in% c(13:16))] <- "P6"

disgenet_EL <- left_join(disgenet_EL, chd8_table, by=c("Gene"="GENE"))

chd8_table_disTOTAL <- c()
for (z in colnames(chd8_table)[-1]){
  print(z)
  chd8_table_dis <- melt(table(disgenet_EL[,c("Disease",z)]))  
  chd8_table_dis <- left_join(chd8_table_dis, df2plotDis, by=c("Disease"="disease"))
  chd8_table_dis$Disease <- factor(chd8_table_dis$Disease, levels=order_mix_dis)
  colnames(chd8_table_dis)[2] <- "CHD8"
  mytable <- chd8_table_dis
  fmrp.p <- c()
  fmrp.or <- c()
  abcd <- c()
  for (i in unique(chd8_table_dis$Disease)){
    a <- filter(mytable, CHD8==1, Disease==i)$value
    b <- filter(mytable, CHD8==0, Disease==i)$value
    c <- filter(mytable, CHD8==1, Disease=="ALL")$value - a
    d <- filter(mytable, CHD8==0, Disease=="ALL")$value - b
    f <- fisher.test(matrix(c(a,c,b,d), nrow=2), alternative = "greater")
    fmrp.p <- c(fmrp.p, f$p.value)
    fmrp.or <- c(fmrp.or, f$estimate)
    abcd<-rbind(abcd, c(a,b,c,d))
  }
  names(fmrp.p) <- unique(chd8_table_dis$Disease)
  names(fmrp.or) <- unique(chd8_table_dis$Disease)
  test <- data.frame(Disease = unique(chd8_table_dis$Disease), P=fmrp.p, OR=fmrp.or, PAdj=p.adjust(fmrp.p, method = "BH"))
  test <- cbind(test,abcd)
  chd8_table_dis <- left_join(chd8_table_dis, test)
  chd8_table_dis$DATASET <- z
  chd8_table_disTOTAL <- rbind(chd8_table_disTOTAL,chd8_table_dis)
}

chd8_table_disTOTAL$SET <- ifelse(chd8_table_disTOTAL$DATASET %in% c("FMRP_Casingal","FMRP_Darnell"), "FMRP","CHD8" )
chd8_table_disTOTAL$CHD8[which(chd8_table_disTOTAL$SET=="FMRP" & chd8_table_disTOTAL$CHD8==1)] <- 2

p2 <- ggplot(chd8_table_disTOTAL, aes(y=Disease, x=value, alpha=as.factor(ifelse(PAdj < 0.05,1,0)), fill=factor(CHD8))) + geom_bar(position="fill", stat="identity") + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_alpha_manual(values = c(0.4,1))  + 
  scale_fill_manual(values = c("grey85","tomato", "cadetblue")) +coord_cartesian(xlim=c(0,1))+ facet_grid(Passage2plot~DATASET, space = "free_y", scales = "free_y")


p2
ggsave("AllDiseaseCHD8_FMRP_Complete_V3.pdf", height = 80, width=50,units = "cm")
 