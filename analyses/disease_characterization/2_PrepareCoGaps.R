
### Filter CoGAPS

require(dplyr)
library("readxl")
library(plyr)

setwd("/Users/gabrielsantpere/Desktop/Neurogenomics/Projectes/Nicola/")

cogaps<-read_excel("DATA/CoGAPS/GSE144156_H9_nicolaFGFnscRNAseqHs_GWCoGAPS_p24.xlsx")
cogaps<-as.data.frame(cogaps)

for(i in paste0("p",c(1:24))){
    temp<-cogaps[,c("gene",i)]
  gene <- temp$gene[ which(temp[,2] > quantile(temp[,2], 0.95)) ]
  write.table(as.data.frame(gene), file=paste0("CoGAPS/",i,"_95.txt"), sep="\t", quote=F, row.names = F, col.names = F)
}

### Select only protein-coding genes for MAGMA analysis

library(biomaRt)

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genesV2 = getBM(filters = "hgnc_symbol", values = cogaps$gene , mart = human, attributes = c("hgnc_symbol","gene_biotype"), uniqueRows=T)

## 3049 genes with no info in BioMart
## 16038 are PCG and 879 ar non-coding

pc<-filter(genesV2, gene_biotype=="protein_coding")$hgnc_symbol
cogaps_clean <- filter(cogaps, gene %in% pc)

geneFilter <- c()
for(i in paste0("p",c(1:24))){
  temp<-cogaps_clean[,c("gene",i)]
  gene <- temp$gene[ which(temp[,2] > quantile(temp[,2], 0.95)) ]
  df<-data.frame(dataset=i,gene=gene)
  geneFilter <- rbind(geneFilter,df)
}


genLoc <- read.table("Prepare_MAGMA/NCBI37.3.gene.loc")
geneFilter <-  filter(geneFilter, gene %in% genLoc$V6)
geneFilter$geneSymbol <- mapvalues( as.character(geneFilter$gene), from = genLoc$V6, to = genLoc$V1  )

write.table(geneFilter[,c(1,3)], col.names = F, row.names = F, quote = F, sep="\t", file="Cogaps_top5.txt")

#for i in $(ls GWAS_MAGMA/RESULTS/*raw); do name=$(basename $i); magma --gene-results ${i} --set-annot GENESETS/Cogaps_top5.txt col=2,1 --out GENESETS/ENRICH_COGAPS/${name}; done


##### CoGAPS form the 6 cell ines for dorsal ventral analysis

cogaps6<-read_excel("DATA/CoGAPS/GSE144157_iPSC6lines_BRNiPSCnpcNRN_GWCoGAPS_p30.xlsx")
cogaps6<-as.data.frame(cogaps6)
colnames(cogaps6)[1] <- "gene"
for(i in paste0("p",c(1:30))){
  temp<-cogaps6[,c("gene",i)]
  gene <- temp$gene[ which(temp[,2] > quantile(temp[,2], 0.95)) ]
  write.table(as.data.frame(gene), file=paste0("CoGAPS/iPSC6lines/",i,"_95_6lines.txt"), sep="\t", quote=F, row.names = F, col.names = F)
}


library(biomaRt)

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genesV2 = getBM(filters = "hgnc_symbol", values = cogaps6$gene , mart = human, attributes = c("hgnc_symbol","gene_biotype"), uniqueRows=T)

## 3049 genes with no info in BioMart
## 16038 are PCG and 879 ar non-coding

pc<-filter(genesV2, gene_biotype=="protein_coding")$hgnc_symbol
cogaps_clean6 <- filter(cogaps6, gene %in% pc)

geneFilter <- c()
for(i in paste0("p",c(1:30))){
  temp<-cogaps_clean6[,c("gene",i)]
  gene <- temp$gene[ which(temp[,2] > quantile(temp[,2], 0.95)) ]
  df<-data.frame(dataset=i,gene=gene)
  geneFilter <- rbind(geneFilter,df)
}


genLoc <- read.table("Prepare_MAGMA/NCBI37.3.gene.loc")
geneFilter <-  filter(geneFilter, gene %in% genLoc$V6)
geneFilter$geneSymbol <- mapvalues( as.character(geneFilter$gene), from = genLoc$V6, to = genLoc$V1  )

write.table(geneFilter[,c(1,3)], col.names = F, row.names = F, quote = F, sep="\t", file="Cogaps6_top5.txt")

##### CoGAPS form the Late Differentiation AZ lines

cogapsL <- read.table("DATA/CoGAPS/AZpilot_ALL.cells_CoGAPS_featureLoadings_139samples_52748genes_1sets_12patterns_50000iterations_12fP_PROP_GeneSymbols.txt", header = T, sep="\t")
cogapsL <- cogapsL[,c(2:14)]

colnames(cogapsL)[1] <- "gene"
for(i in paste0("Pattern_",c(1:12))){
  temp<-cogapsL[,c("gene",i)]
  gene <- temp$gene[ which(temp[,2] > quantile(temp[,2], 0.95)) ]
  write.table(as.data.frame(gene), file=paste0("CoGAPS/AZlines/",i,"_95_AZlines.txt"), sep="\t", quote=F, row.names = F, col.names = F)
}


library(biomaRt)

human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
genesV2 = getBM(filters = "ensembl_gene_id", values = cogapsL$gene , mart = human, attributes = c("ensembl_gene_id","hgnc_symbol","gene_biotype"), uniqueRows=T)

## 3049 genes with no info in BioMart
## 16038 are PCG and 879 ar non-coding

pc<-filter(genesV2, gene_biotype=="protein_coding")$ensembl_gene_id
cogaps_cleanL <- filter(cogapsL, gene %in% pc)

geneFilter <- c()
for(i in paste0("Pattern_",c(1:12))){
  temp<-cogaps_cleanL[,c("gene",i)]
  gene <- temp$gene[ which(temp[,2] > quantile(temp[,2], 0.95)) ]
  df<-data.frame(dataset=i,gene=gene)
  geneFilter <- rbind(geneFilter,df)
}


genLoc <- read.table("Prepare_MAGMA/NCBI37.3.gene.loc")
geneFilter$geneSymbol <- mapvalues( as.character(geneFilter$gene), from = genesV2$ensembl_gene_id, to = genesV2$hgnc_symbol  )
geneFilter <-  filter(geneFilter, geneSymbol %in% genLoc$V6)
geneFilter$geneRef <- mapvalues( as.character(geneFilter$geneSymbol), from = genLoc$V6, to = genLoc$V1  )

write.table(geneFilter[,c(1,4)], col.names = F, row.names = F, quote = F, sep="\t", file="CogapsAZ_top5.txt")
