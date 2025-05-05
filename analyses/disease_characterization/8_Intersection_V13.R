
### Intersection of disease genes

require(dplyr)
library("readxl")
library(plyr)
library(biomaRt)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(reshape2)
library(RColorBrewer)


setwd("/Users/gabrielsantpere/Desktop/Neurogenomics/Projectes/CortexMalformation/")

### Set color pallet

c1<-rev(brewer.pal(n = 5, name = "Greens")[-1])
c2<-rev(brewer.pal(n = 5, name = "Blues")[-1])
c3<-rev(brewer.pal(n = 5, name = "Reds")[-1])
c4<-rev(brewer.pal(n = 5, name = "Purples")[-1])
c5<-rev(brewer.pal(n = 5, name = "Greys")[-1])

col <- c(c1,c2,c3,c4,c5)

### Read patterns expression

load("DATA/CoGAPS_Expression/NicoHsPASS_ExprsLog2.Z_forGabrielHeatmaps.RData")
load("DATA/CoGAPS_Expression/NicoDV_ExprsLog2.Z_forGabrielHeatmaps.RData")
colnames(NicoDV.ExprsLog2.Z)=NicoDV.sampleMeta$ids
## Import patterns

out <- read.table("MAGMA/RESULTS/GENE_LEVEL/AllDiseasesMAGMA_GeneLevel_181220.txt",sep = "\t",header=T)
out <- left_join(data.frame(GENE_SYMBOL=NicoHsPASS.geneMeta$GeneSymbol), out)

## Import GWAS genome-wide

gw <- read.table("DATA/DiseaseGenes/GWAS_genes.txt", header = T)
gwm <- data.frame( GENE = out$GENE_SYMBOL  )

for (i in colnames(out)[-c(1:2)]){
  gwm <- cbind( gwm  , ifelse( gwm$GENE %in% dplyr::filter(gw, Dataset == i)$Gene, 1, 0 )  )
}

colnames(gwm)[-1] <- colnames(out)[-c(1:2)]

### Input TFs

tfs <- read.table("DATA/HumanTFs.txt", header = FALSE)$V1

### Input DeCon

load("DATA/Slopes/DeconSlopes_V3_Scaled.rda")
load("DATA/Slopes/telley_slope_V2_Scaled.rda")
decon_slope <- decon_slope[-which(duplicated(decon_slope$GeneHuman)),]
decon_slope[,c(3:5)] <- apply(decon_slope[,c(3:5)],2,as.numeric)
telley_slope <- telley_slope[-which(duplicated(telley_slope$GeneHuman)),]
telley_slope[,c(3:6)] <- apply(telley_slope[,c(3:6)],2,as.numeric)

### Input malformations

mal <- read.table("DATA/DiseaseGenes/CorticalMalformationAndASD.txt", header = T)
fmrp <- read.table("DATA/DiseaseGenes/FMRP_targets.txt")
chd8 <- read.table("DATA/DiseaseGenes/CHD8_targets.txt")

malg <- data.frame( GENE = out$GENE_SYMBOL  )

for (i in unique(mal$DATASET)){
  malg <- cbind( malg, ifelse( malg$GENE %in% filter(mal, DATASET == i)$GENE, 1, 0 )  )
}

colnames(malg)[-1] <- unique(mal$DATASET)
rownames(malg) <- malg$GENE

#### ADD hydrochephaly and SFARI genes

s1 <- read.table("DATA/DiseaseGenes/SFARI_Score1.txt")
s2 <- read.table("DATA/DiseaseGenes/SFARI_Score2.txt")
s3 <- read.table("DATA/DiseaseGenes/SFARI_Score3.txt")
ss <- read.table("DATA/DiseaseGenes/SFARI_Syndromic.txt")
dd <- read.table("DATA/DiseaseGenes/DD_genes.txt")
hc <- read.table("DATA/DiseaseGenes/ASD_HC65.txt")
hy <- read.table("DATA/DiseaseGenes/Hydrocephally.txt")

malg$Hydrocephaly <- ifelse( malg$GENE %in% hy$V1, 1, 0 )
malg$DD <- ifelse( malg$GENE %in% dd$V1, 1, 0 )
malg$SFAR_S1 <- ifelse( malg$GENE %in% s1$V1, 1, 0 )
malg$SFAR_S2 <- ifelse( malg$GENE %in% s2$V1, 1, 0 )
malg$SFAR_S3 <- ifelse( malg$GENE %in% s3$V1, 1, 0 )
malg$SFAR_Synd <- ifelse( malg$GENE %in% ss$V1, 1, 0 )
malg$ASD_HC65 <- ifelse( malg$GENE %in% hc$V1, 1, 0 )

load("DATA/DiseaseGenes/FMRP_CHD8_TOTAL.rda")  

malg$FMRP <- ifelse( malg$GENE %in% fmrp_TOTAL, 1, 0 )
malg$CHD8 <- ifelse( malg$GENE %in% chd8_TOTAL, 1, 0 )

malg <- malg[,-10]
malg <- malg[, c(1,3,2,4:18)]

disease_full <- cbind(malg,gwm[,-1] )
#save(disease_full, file = "DATA/DiseaseGenes/DiseaseFULL.rda")
#write.table(malg, file="DATA/DiseaseGenes/CorticalMalformation_ASD_Marc2021.txt", row.names = FALSE, sep="\t", quote=F)
load("DATA/DiseaseGenes/DiseaseFULL.rda")

### Plot one specific pattern
## 1) Set significance parameters
max_list <- data.frame(Disease=character(), Max=integer())
load("DATA/Slopes/DorsalVentralDeSeq2.rda")
DV_FC <- data.frame(Day8=r8$log2FoldChange, Day17=r17$log2FoldChange, Day30=r30$log2FoldChange, row.names = rownames(r8))
minDV <- min(DV_FC)
maxDV <- max(DV_FC)
DV_FC$GENE <- rownames(DV_FC)

for(disease in colnames(disease_full)[-c(1,17,18)]){
  print(disease)
  if ( sum(disease_full[,disease])==0 ) { next }
  mask <- which(disease_full[,disease] == 1)
  
### 2) Select patterns
  
out_t <- out[mask,]
gwm_t <- gwm[mask,]
malg_t <- malg[mask,]

mat <- out_t[,-c(1,2)]
rownames(mat) <- (out_t$GENE_SYMBOL)
mat <- as.matrix(mat)
gwm_t <- as.matrix(gwm_t[,-1])
rownames(gwm_t) <- (out_t$GENE_SYMBOL)

th <- -log10(0.05/24219)
deconm <- left_join(data.frame(GENE = rownames(gwm_t)), decon_slope, by=c("GENE"="GeneHuman"))
rownames(deconm) <- deconm$GENE 
telleym <- left_join(data.frame(GENE = rownames(gwm_t)), telley_slope, by=c("GENE"="GeneHuman"))
rownames(telleym) <- telleym$GENE


E <- as.data.frame(NicoHsPASS.ExprsLog2.Z)
E$GENE <- rownames(NicoHsPASS.ExprsLog2.Z)
early <- left_join(data.frame(GENE = rownames(gwm_t)), E)[,-1]
early<-t(apply(early, 1, function(x){ (x-min(x))/(max(x)-min(x))   }))

early_max <- apply(early, 1, function(x){ which(x==max(x)) })
df <- data.frame(Disease=disease, Max=early_max)
max_list <- rbind(max_list, df)

dv <- left_join(data.frame(GENE = rownames(gwm_t)), DV_FC)[,-1]
  
### 3) Complex heatmap

col_v <- col
names(col_v) <- c(1:20)
#mapvalues(early_max, from = c(1:20), to = col)
hmax <- HeatmapAnnotation(
  MAX= early_max,
  col = list(MAX = col_v), 
  annotation_name_gp = gpar(col = "black"),
  which =  "row" , gap = unit(0.4, "mm"), border = FALSE,
  gp = gpar(col = "black", lwd = 0.15),
  show_legend = c("MAX" = FALSE)
)

hw <- Heatmap(
  as.matrix(early), row_order = order(early_max),width = unit(9, "cm"),
  cluster_columns  = FALSE,
  col = colorRamp2(c( 0,1 ), c("yellow", "red")),
  rect_gp = gpar(col = "black", lwd = 0.25),
  left_annotation = c(hmax, gap =unit(2, "mm")))

hw2 <- HeatmapAnnotation(
  #H9=as.matrix(early[,-c(5:8, 13:16)]),
  DV=as.matrix(dv),
  col = list(
    DV = colorRamp2(c( -5,-0.5,0,0.5,5) , c("forestgreen", "#90C490" ,"grey97", "#F6D990" ,"#EEB422"))
    #DV = colorRamp2(c( minDV,0,maxDV) , c("darkred","white", "cadetblue"))
    #DV = circlize::colorRamp2(c(minDV,-1.000000001,1,-1,1.000000001,maxDV), 
    #                          c("darkred","pink","grey97", "grey97","cadetblue","midnightblue"))
  ),
  gp = gpar(col = "black", lwd = 0.25),
  which =  "row" , gap = unit(2, "mm"), show_legend = c("DV" = TRUE))



ha <- HeatmapAnnotation(
  TELLEY=as.matrix(telleym[,-c(1,2)]),
  DECON=as.matrix(deconm[,-c(1,2)]),
  col = list(DECON =  colorRamp2(c( min(deconm[,-c(1,2)], na.rm = TRUE),0, max(deconm[,-c(1,2)], na.rm = TRUE) ), c("darkred","grey97", "cadetblue")),
             TELLEY = colorRamp2(c(min(telleym[,-c(1,2)], na.rm = TRUE),0 ,max(telleym[,-c(1,2)], na.rm = TRUE) ), c("darkred","grey97", "cadetblue"))
             ),
  which =  "row" , gap = unit(2, "mm"), show_legend = c("DECON" = TRUE, "TELLEY"=TRUE))

col_bin<-c("0" = "grey97", "1" = "cadetblue")

hb <- HeatmapAnnotation(
  FMRP=  malg_t$FMRP,
  CHD8=  malg_t$CHD8,
  col = list(FMRP = col_bin, CHD8 = col_bin), 
  annotation_name_gp = gpar(col = "black"),
  which =  "row" , gap = unit(0.4, "mm"), border = FALSE,
  gp = gpar(col = "black", lwd = 0.15),
  show_legend = c("FMRP" = FALSE, "CHD8" = FALSE)
)

hc <- HeatmapAnnotation(
  Malformations = as.matrix(malg_t[,-c(1, 17,18)]),
  col = list(Malformations = c("0" = "grey97", "1" = "tomato")),
  annotation_name_gp = gpar(col = "black"),
  which =  "row" , gap = unit(0.4, "mm"), border = FALSE,
  gp = gpar(col = "black", lwd = 0.15),
  show_legend = c("Malformations" = FALSE)
)


ht <- Heatmap(-log10(mat), width = unit(7, "cm"),
        
        col = circlize::colorRamp2(c(0,th-0.000000001,th,30), 
                                   c("grey97","grey97","bisque","midnightblue")),
        row_names_gp = gpar(col = ifelse(rownames(mat) %in% tfs, "red","black")),
        rect_gp = gpar(col = "black", lwd = 0.15), cluster_rows = FALSE, 
        
        layer_fun = function(j, i, x, y, width, height, fill) {
          v = pindex(gwm_t, i, j)
          print(v)
          l = v > 0
          grid.points(x[l], y[l], pch = 16, gp = gpar(col="black"), size = unit(2, "mm"))
        },
        
        right_annotation = c(hb, gap =unit(2, "mm")),
        left_annotation = c(hw2,ha,hc, gap =unit(2, "mm"))
      )


pdf(file = paste0("Figures/Complex/Heatmap_DV_Slope_LIMIT_Rearranged_TARGETS_",disease,".pdf"), width = 20 , height = max(10,nrow(mat)/7))
  draw(hw+ht, heatmap_legend_side = "left", annotation_legend_side = "bottom", show_row_dend = FALSE, 
     heatmap_height=unit(max(10,nrow(mat)/3),"cm"), ht_gap=unit(4, "mm"))     
dev.off()

}


#### Plot distribution of maximums

max_list_early <- max_list
max_list_early$Early <- ifelse( max_list_early$Max %in% c(1:4), "Early", "Late" )
max_prop_early<-as.matrix((table(max_list_early$Disease, max_list_early$Early)))
order2<-rev(names(sort(max_prop_early[,1]/(max_prop_early[,1]+max_prop_early[,2]))))

max_table <- melt(table(max_list))
max_table_prop <- melt(prop.table(table(max_list), "Disease"))
order<-colnames(disease_full)[-c(1,17,18)]
order<-order[c(2,1, 3:15)]
order2<- rev((filter(max_table_prop, Max==1) %>% arrange(value))$Disease)
max_table$Disease <- factor(max_table$Disease, levels=order2)

save(max_list, file="DATA/CoGAPS_Expression/Maximums.rda")

ggplot(max_table, aes(x=Disease, y=value, fill=factor(Max, levels = rev(c(1:20))))) + geom_bar(position="fill", stat="identity") + theme_classic() +
 theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_fill_manual(values = rev(col))

