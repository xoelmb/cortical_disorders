

#### Calculate Dorsal versus ventral gene expression from counts
#### Combine with FOLD_CHANGE

require(dplyr)
library("readxl")
library(plyr)
library(DESeq2)
library(ggplot2)
library(ggrepel)

setwd("/Users/gabrielsantpere/Desktop/Neurogenomics/Projectes/Nicola/")

### 1 ---- Load slopes previously calculated

load("DATA/Slopes/DV_slope.rda")


### 2 ---- Load counts and run Deseq2

load("DATA/CoGAPS_Expression/BRNiPSCnpcNRN_counts.rdata")
load("DATA/CoGAPS_Expression/NicoDV_ExprsLog2_forGabrielHeatmaps.rdata")

ids <- paste0(sampleMETA$SubjectID,".",sampleMETA$LineRep,".day",sampleMETA$Days)
colnames(counts) <- ids

dv_info <- ifelse(grepl("2053.6|2075", ids), "dorsal", "ventral")

coldata <- data.frame( row.names=ids, Day = sampleMETA$Days, DV = dv_info )
coldata$group <- factor(paste0( coldata$DV, "_" , coldata$Day ))
  
dds <- DESeqDataSetFromMatrix(countData=counts, colData=coldata, design=~group)
dds <- DESeq(dds)
resultsNames(dds)

r8 <- results(dds, contrast=c("group", "dorsal_8", "ventral_8"))
r17 <- results(dds, contrast=c("group", "dorsal_17", "ventral_17"))
r30 <- results(dds, contrast=c("group", "dorsal_30", "ventral_30"))

save(r8, r17, r30, file="DATA/Slopes/DorsalVentralDeSeq2.rda")

colnames(NicoDV.ExprsLog2) <- mapvalues(colnames(NicoDV.ExprsLog2), from = sampleMETA$SampleID, to= ids)
col <- rep("black", nrow(r8))
expressed <- which(rowMeans(NicoDV.ExprsLog2[,grep("day8", colnames(NicoDV.ExprsLog2)) ]) > 1 )
col[expressed] <- "red"

par(mfrow=c(1,3))
plot(r8$log2FoldChange, DV_FC$Day8, pch=19, cex=0.8, xlab = "DESEq2FC", ylab="FC from log2RPKM", xlim=c(-6,6), ylim=c(-6,6), col=col, main="Day 8")
plot(r17$log2FoldChange, DV_FC$Day17, pch=19, cex=0.8, xlab = "DESEq2FC", ylab="FC from log2RPKM", xlim=c(-6,6), ylim=c(-6,6), col=col, main="Day 17")
plot(r30$log2FoldChange, DV_FC$Day30, pch=19, cex=0.8, xlab = "DESEq2FC", ylab="FC from log2RPKM", xlim=c(-6,6), ylim=c(-6,6), col=col, main="Day 30")

### 3 --- Volcano plots

load("DATA/DiseaseGenes/DiseaseFULL.rda")

for(disease in colnames(disease_full)[-c(1,17,18)]){
  
  print(disease)
  mask <- rownames(disease_full)[which(disease_full[,disease] == 1)]
  col <- ifelse( rownames(r8) %in% mask, "green", "darkgrey"  )
  m <- which(rownames(r8) %in% mask)
  
  df1 <- data.frame(FC=r8$log2FoldChange, P=-log10(r8$padj), COL=col, DAY="Day 8", GENE=rownames(r8))
  df2 <- data.frame(FC=r17$log2FoldChange, P=-log10(r17$padj), COL=col,DAY="Day 17", GENE=rownames(r8))
  df3 <- data.frame(FC=r30$log2FoldChange, P=-log10(r30$padj), COL=col, DAY="Day 30", GENE=rownames(r8))
  df <- rbind(df1, df2, df3)
  
  df$text<-df$GENE
  df$text[-which(abs(df$FC)>1 & df$COL=="green")] <- ""
  df$DAY <- factor(df$DAY, levels = c("Day 8", "Day 17", "Day 30"))
  
  p <- ggplot( df, aes(FC, P, fill=COL) ) + geom_point(col=df$COL) + 
    scale_fill_manual(values = c("green","grey")) + 
    facet_grid(.~DAY) + theme_bw() + 
    geom_text_repel(aes(x=FC,y=P,label=text), box.padding = 0.5, max.overlaps = Inf) + 
    geom_hline(yintercept=-log10(0.05), linetype="dashed", color = "grey") +
    geom_vline(xintercept = 0) + geom_point(data=df[which(df$COL=="green"),], col="tomato") +
    theme(legend.position = "none") 
  
  ggsave(plot = p, filename =paste0("Figures/DorsoVentral/VOLCANO/", disease, ".pdf"), width = 15, height = 5)

}  
