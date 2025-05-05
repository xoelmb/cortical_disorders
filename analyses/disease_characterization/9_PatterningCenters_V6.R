

#### Dev mouse single-cell data

### 1 --- Run script in cluster to obtain Dotplot data from DevMouse

# module load R/4.0.3-foss-2020b
# module load GLib/2.66.1-GCCcore-10.2.0
# srun --pty --mem 90G -p interactive bash
# sbatch --mem=90G --wrap="module load R/4.0.3-foss-2020b; cd /home/gs574/project/SingleCell/DEVMOUSE; Rscript Prop.R"
## In Ruddle with 90Gb core

#d <- readRDS("/gpfs/ycga/project/sestan/sestanlabShare_from_Jay/publicdata/mouse_dev_atlas/mdev_seu_0503_2021.rds")
d <- readRDS("/gpfs/ycga/project/sestan/sestanlabShare_from_Jay/publicdata/mouse_dev_atlas/mdev_normalized_seurat_09172021.rds")
areas <- c(  "Anteromedial cerebral pole","Cortical hem","Midbrain floor plate","Forebrain glutamatergic","Forebrain GABAergic", "Forebrain astrocyte","Oligodendrocyte", "Dorsal forebrain", "Forebrain", "Neuronal intermediate progenitor")

jj <- subset(x = d, subset = Class %in% c("Neuron","Neuroblast","Radial glia","Oligodendrocyte", "Glioblast") & Subclass %in% areas & Clusters != "214" & Age %in% c("e9.0", "e10.0") )
jj$Myclass <- paste0(jj$Class,"_",jj$Subclass)
Idents(jj) <- "Myclass"

chunks <- split(rownames(jj), ceiling(seq_along(rownames(jj))/150))

pcttable <- c()

for (i in c(1:length(chunks)) ){
  print(i)
  a<-DotPlot(jj, features = chunks[[i]])
  pcttable <- rbind(pcttable, a$data)
}  

save(pcttable, file = "pcttable.Giole.Sept2021.rda")

## Check differential expression

areas <- c(  "Anteromedial cerebral pole","Cortical hem","Midbrain floor plate", "Dorsal forebrain", "Forebrain", "Neuronal intermediate progenitor")
jj <- subset(x = d, subset = Class %in% c("Radial glia") & Subclass %in% areas & Clusters != "214" & (Age %in% c("e9.0", "e10.0") | (Age %in% c("e9.0", "e10.0", "e11.0", "e12.0") & Subclass == "Cortical hem")))
jj$Myclass <- paste0(jj$Class,"_",jj$Subclass)
Idents(jj) <- "Myclass"
metaearly <- jj@meta.data
markers_early <- FindAllMarkers(jj, test.use = "LR", latent.vars="NGenes")
saveRDS(markers_early, file="markers_early.rds")

load("DATA/MOUSE_DEV_ATLAS/GioleMetadata.rda")

########################################################################################3
######## Integrate with Deseq2 dorsal ventral comparisons


require(dplyr)
library(ggplot2)
library(biomaRt)
library(plyr)
library(ggpubr)
library(ComplexHeatmap)
library(RColorBrewer)
library(circlize)

setwd("/Users/gabrielsantpere/Desktop/Neurogenomics/Projectes/CortexMalformation//")

load("DATA/Slopes/DorsalVentralDeSeq2.rda")
#load("DATA/MOUSE_DEV_ATLAS/pcttable.Giole.August2021.rda")
load("DATA/MOUSE_DEV_ATLAS/pcttable.Giole.Sept2021_early.rda")
pearly <- pcttable 
load("DATA/MOUSE_DEV_ATLAS/pcttable.Giole.Sept2021.rda")
load("DATA/DiseaseGenes/DiseaseFULL.rda")

convertMouseGeneList <- function(x){
  human = useEnsembl("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useEnsembl("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[,])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
}

library(Orthology.eg.db)
library(org.Hs.eg.db)
library(org.Mm.eg.db)

fun <- function(gns){
  egs <- mapIds(org.Mm.eg.db, gns, "ENTREZID","SYMBOL")
  mapped <- select(Orthology.eg.db, egs, "Homo.sapiens","Mus.musculus")
  mapped$HS <- mapIds(org.Hs.eg.db, as.character(mapped$Homo.sapiens), "SYMBOL", "ENTREZID")
  mapped$MUS <- names(egs)
  mapped
}

#httr::set_config(httr::config(ssl_verifypeer = FALSE))

convertor <- fun(as.character(pcttable$features.plot))
convertor <- unique(convertor)
convertor$HS <- as.character(as.vector(convertor$HS))
convertor$MUS <- as.character(as.vector(convertor$MUS))
saveRDS(convertor, file = "Convertor_May2024.rds")

pcttable$idHuman <- plyr::mapvalues(pcttable$features.plot, from=convertor$MUS, to=convertor$HS)
pearly$idHuman <- plyr::mapvalues(pearly$features.plot, from=convertor$MUS, to=convertor$HS)

areas <- c(  "Radial glia_Anteromedial cerebral pole","Radial glia_Cortical hem", "Radial glia_Midbrain floor plate","Radial glia_Dorsal forebrain","Radial glia_Forebrain",
             "Neuroblast_Neuronal intermediate progenitor", "Neuroblast_Forebrain" ,"Neuroblast_Forebrain glutamatergic",
             "Neuron_Forebrain glutamatergic", "Neuroblast_Forebrain GABAergic", "Neuron_Forebrain GABAergic","Glioblast_Forebrain","Glioblast_Forebrain astrocyte",
             "Oligodendrocyte_Oligodendrocyte")

areas_early <- c(  "Radial glia_Anteromedial cerebral pole","Radial glia_Cortical hem", "Radial glia_Midbrain floor plate","Radial glia_Dorsal forebrain","Radial glia_Forebrain")

pct <- dplyr::filter(pcttable, id %in% areas)

#### Select dorsal ventral disease genes

r8$GENE <- rownames(r8)
r17$GENE <- rownames(r17)
r30$GENE <- rownames(r30)

FC_DISEASE <- c()
for(disease in colnames(disease_full)[-c(1,17,18,28)]){
  print(disease)
  if ( sum(disease_full[,disease])==0 ) { next }
  mask <- rownames(disease_full)[which(disease_full[,disease] == 1)]
  df8 <- dplyr::filter(as.data.frame(r8), padj < 0.05, GENE %in% mask)
  df17 <- dplyr::filter(as.data.frame(r17), padj < 0.05, GENE %in% mask)
  df30 <- dplyr::filter(as.data.frame(r30), padj < 0.05, GENE %in% mask)
  if(nrow(df8)>0){df8$DAY <- "Day8"}
  if(nrow(df17)>0){df17$DAY <- "Day17"}
  if(nrow(df30)>0){df30$DAY <- "Day30"}
  bigdf <- rbind(df8,df17,df30)
  bigdf$DISEASE=disease
  FC_DISEASE <- rbind(FC_DISEASE,bigdf)
}

FC_DISEASE <- FC_DISEASE[,c(1,2,6,7,8,9)]

### Join disease DV and patterning expression

combined <- left_join(FC_DISEASE, pct, by = c("GENE" = "idHuman" ))
combined$PATTERN <- ifelse(combined$log2FoldChange < 0, "VENTRAL", "DORSAL")
combined <- combined %>% group_by(GENE, DAY, DISEASE) %>% dplyr::mutate(logAvgMinMax = (avg.exp-min(avg.exp))/(max(avg.exp)-min(avg.exp)) )

combined_early <- left_join(FC_DISEASE, pearly, by = c("GENE" = "idHuman" ))
combined_early$PATTERN <- ifelse(combined_early$log2FoldChange < 0, "VENTRAL", "DORSAL")
combined_early <- combined_early %>% group_by(GENE, DAY, DISEASE) %>% dplyr::mutate(logAvgMinMax = (avg.exp-min(avg.exp))/(max(avg.exp)-min(avg.exp)) )


save(combined, file="DATA/MOUSE_DEV_ATLAS/CombinedFCandLaMannoMay2024.rda")
save(combined_early, file="DATA/MOUSE_DEV_ATLAS/CombinedEarlyFCandLaMannoMay2024.rda")

#load("DATA/MOUSE_DEV_ATLAS/CombinedFCandLaMannoSept2021.rda")
#load("DATA/MOUSE_DEV_ATLAS/CombinedEarlyFCandLaMannoSept2021.rda")

### Select only day 8

combined8 <- filter(combined, DAY=="Day8")
combined8<-combined8[!is.na(combined8$id),]
combined8$pct.exp[which(combined8$pct.exp < 5)] <- NA
combined8$id <- factor(combined8$id, levels=areas)

combined8_early <- filter(combined_early, DAY=="Day8")
combined8_early<-combined8_early[!is.na(combined8_early$id),]
combined8_early$pct.exp[which(combined8_early$pct.exp < 5)] <- NA
combined8_early$id <- factor(combined8_early$id, levels=areas_early)

p2 <- ggplot(combined8, aes(x=DAY, y=forcats::fct_reorder(GENE, log2FoldChange), fill=PATTERN) ) + geom_tile() + theme_light() +
  scale_fill_manual(values = c("goldenrod2","forestgreen")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.title.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid = element_blank(), plot.margin = unit(c(0,0.5,0,-0.5),"cm"))

p2 <- ggplot(combined8, aes(x=DAY, y=GENE, fill=PATTERN) ) + geom_tile() + theme_light() +
  scale_fill_manual(values = c("goldenrod2","forestgreen")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  theme(axis.title.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(),
        panel.grid = element_blank(), plot.margin = unit(c(0,0.5,0,-0.5),"cm"))




### Plot results
for(disease in unique(combined$DISEASE)){
  print(disease)
  temp <- filter(combined, DISEASE==disease)
  temp$id <- factor(temp$id, levels = areas)
  temp$DAY <- factor(temp$DAY, levels = c("Day8","Day17","Day30"))
  temp<-temp[!is.na(temp$id),]

  p2 <- ggplot(temp, aes(x=DAY, y=forcats::fct_reorder(GENE, log2FoldChange), fill=PATTERN) ) + geom_tile() + theme_light() +
    scale_fill_manual(values = c("goldenrod2","forestgreen")) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
    theme(axis.title.x = element_blank(), axis.text.y = element_blank(), axis.title.y = element_blank(), axis.ticks.y = element_blank(),
          panel.grid = element_blank(), plot.margin = unit(c(0,0.5,0,-0.5),"cm"))
  
  # p <- ggplot(unique(temp[,-5]), aes(x=id, y=forcats::fct_reorder(GENE, log2FoldChange)) ) + geom_point( aes(size=pct.exp, color=logAvgMinMax )  ) + theme_light() +
  # theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  # scale_color_gradient(low="cadetblue", high = "tomato") +
  #   theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(0,-0.5,0,0.5),"cm"))
  
  temp$pct.exp[which(temp$pct.exp < 5)] <- NA
  p <- ggplot(unique(temp[,-5]), aes(x=id, y=forcats::fct_reorder(GENE, log2FoldChange)) ) + geom_point( aes(size=pct.exp, color=logAvgMinMax )  ) + theme_light() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    scale_color_gradient(low="cadetblue", high = "tomato") + 
    theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(0,-0.5,0,0.5),"cm"))

  if(length(unique(temp$PATTERN)) > 1){
    p <- p + facet_grid(PATTERN~., scales = "free", space = "free_y") + theme(strip.background = element_blank(),strip.text.x = element_blank())
    p2 <- p2 + facet_grid(PATTERN~., scales = "free", space = "free_y") 
  }

  ggarrange(p,p2, align="hv", widths = c(7,4), common.legend = TRUE, nrow = 1)
  ggarrange(p,p2, align="hv", widths = c(5,4), common.legend = TRUE, nrow = 1)
  a <- ggarrange(p,p2, widths = c(15,10), heights = c(20,5), align="hv")
  ggsave(a, filename = paste0("Figures/DorsoVentral/Giole_V2/",disease,".pdf"))
  ggsave(a, filename = paste0("Figures/DorsoVentral/Giole_V2/",disease,".pdf"), width = 7, height = max(5,length(unique(temp$GENE))/4) )
 }

#### Annotations


geneset <- unique(combined8$GENE)
dvsubset <- as.data.frame(unique(combined8[,c(4,12)]))
rownames(dvsubset) <- dvsubset$GENE

c1<-rev(brewer.pal(n = 5, name = "Greens")[-1])
c2<-rev(brewer.pal(n = 5, name = "Blues")[-1])
c3<-rev(brewer.pal(n = 5, name = "Reds")[-1])
c4<-rev(brewer.pal(n = 5, name = "Purples")[-1])
c5<-rev(brewer.pal(n = 5, name = "Greys")[-1])

col <- c(c1,c2,c3,c4,c5)

load("DATA/MaximumTenporal.rda")
early_max <- dfALL[geneset,-1]
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

col_bin<-c("DORSAL" = "forestgreen", "VENTRAL" = "goldenrod2")
hdv <- HeatmapAnnotation(
  DV= dvsubset[,-1],
  col = list(DV = col_bin), 
  annotation_name_gp = gpar(col = "black"),
  which =  "row" , gap = unit(0.4, "mm"), border = FALSE,
  gp = gpar(col = "black", lwd = 0.15),
  show_legend = c("DV" = TRUE)
)

order_disease <- c("Heterotopia", "MDD_2018", "ASD_2019", "ADHD_2019", "DevDyslexia", "FCDandmTOR",
                   "RareMCD", "SFAR_S3", "NEUROT_2018", "SFAR_S2", "SFAR_Synd", "AD_2019", "Lissencephaly",
                   "SFAR_S1", "ASD_HC65", "SCZ_2020", "BD_2019", "Hydrocephaly","AN_2019",
                   "Polymicrogyria","PD_2014","DD","Cobblestone", "Microcephlaly")

malg_t <- disease_full[geneset,order_disease]

hw <- Heatmap(
  as.matrix(malg_t),
  row_order = order(early_max),
  cluster_columns  = FALSE,
  col = colorRamp2(c( 0,1 ), c("grey97", "tomato")),
  rect_gp = gpar(col = "black", lwd = 0.25),
  left_annotation = c(hmax, gap =unit(2, "mm")),
  right_annotation = c(hdv,gap =unit(2, "mm") )
)

draw(hw, heatmap_legend_side = "left", annotation_legend_side = "bottom", show_row_dend = FALSE, ht_gap=unit(4, "mm")) 

pdf(file = paste0("Figures/DorsoVentral/DiseasePart1_2024.pdf"), width = 7 , height = 15)
draw(hw, heatmap_legend_side = "left", annotation_legend_side = "bottom", show_row_dend = FALSE,  ht_gap=unit(4, "mm"))     
dev.off()


order <- geneset[order(early_max)]
combined8$GENE <- factor(combined8$GENE,levels=rev(order))
p1 <- ggplot(combined8, aes(x=id, y=GENE) ) + geom_point( aes(size=pct.exp, color=logAvgMinMax )  ) + theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_color_gradient(low="cadetblue", high = "tomato") + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(0,-0.5,0,0.5),"cm"))

pdf(file = paste0("Figures/DorsoVentral/DiseasePart2_V4_2024.pdf"), width = 5 , height = 15)
p1   
dev.off()

combined8_early$GENE <- factor(combined8_early$GENE,levels=rev(order))
p2 <- ggplot(combined8_early, aes(x=id, y=GENE) ) + geom_point( aes(size=pct.exp, color=log(avg.exp+1) )  ) + theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_color_gradient(low="blue", high = "red") 
 
#+  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(0,-0.5,0,0.5),"cm"))

library(cowplot)
pdf(file = paste0("Figures/DorsoVentral/DiseasePart2_V4_ALL_2024.pdf"), width = 10 , height = 15)
plot_grid(p2,p1)   
dev.off()





p <- ggplot(combined8, aes(x=id, y=GENE) ) + geom_point( aes(size=pct.exp, color=avg.exp.scaled)  ) + theme_light() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_color_gradient(low="cadetblue", high = "tomato") + 
  theme(axis.title.x = element_blank(), axis.title.y = element_blank(), plot.margin = unit(c(0,-0.5,0,0.5),"cm"))
pdf(file = paste0("Figures/DorsoVentral/DiseasePart2_V2log.pdf"), width = 5 , height = 15)
p   
dev.off()

### Include with RG markers of patterning center

markers <- readRDS("DATA/MOUSE_DEV_ATLAS/markers_early.rds")
markers$GENE <- mapvalues(markers$gene, from=convertor$MGI.symbol, to=convertor$HGNC.symbol)
markers_disease <- filter(markers, GENE %in% combined8_early$GENE )
write.table(markers_disease, file="DATA/MOUSE_DEV_ATLAS/Markers_early_disease.tsv", sep="\t",quote=F, row.names = FALSE)
