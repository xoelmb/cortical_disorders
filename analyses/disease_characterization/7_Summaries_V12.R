### Calculate proportions

require(dplyr)
library("readxl")
library(plyr)
library(biomaRt)
library(ComplexHeatmap)
library(circlize)
library(ggplot2)
library(reshape2)
library(RColorBrewer)
library(egg)
library("xlsx")

setwd("/Users/gabrielsantpere/Desktop/Neurogenomics/Projectes/CortexMalformation/")

### Read patterns expression ----

load("DATA/CoGAPS_Expression/NicoHsPASS_ExprsLog2.Z_forGabrielHeatmaps.RData")
load("DATA/CoGAPS_Expression/NicoHsPASS_ExprsLog2_forGabrielHeatmaps.rdata")
load("DATA/DiseaseGenes/DiseaseFULL.rda")
load("DATA/CoGAPS_Expression/Maximums.rda")

E <- as.data.frame(NicoHsPASS.ExprsLog2.Z)
E2 <- as.data.frame(NicoHsPASS.ExprsLog2[rownames(E),])

means <- rowMeans(E2)
bottom <- names( head( sort(means), 1000 ) )
top <- names( tail( sort(means), 1000 ) )
oneRPKM <- names( means[means > 1] )

maximums<-as.matrix(t(apply(E, 1, function(x){ (x-min(x))/(max(x)-min(x))   })))
#maximums_max <- apply(maximums, 1, function(x){ sample(which(x==max(x)), 1) }) ## Cagada guapa
maximums_max2 <- apply(maximums, 1, function(x){ head(which(x==max(x)), 1) })
dfALL <- data.frame(Disease="ALL", Max=maximums_max2)

max_list2 <- c()
for (i in unique(max_list$Disease)){
  g <- disease_full[which(disease_full[,i]==1),]$GENE
  df <- data.frame(Disease=i, Max=maximums_max2[g], Gene=g)
  max_list2 <- rbind(max_list2, df)
}

#save(dfALL, file = "DATA/MaximumTenporal_231121.rda")

load("DATA/MaximumTenporal_231121.rda")
df_bottom <- dfALL[bottom,]
df_bottom$Disease <- "Bottom1000"
df_top <- dfALL[top,]
df_top$Disease <- "Top1000"
df_1RPKM <- dfALL[oneRPKM,]
df_1RPKM$Disease <- ">1 RPKM"

dfALL <- rbind(dfALL, df_bottom, df_top, df_1RPKM, max_list2[,-3])
dfALL$Early <- ifelse( dfALL$Max %in% c(1:4), "P2", "P8" )
dfALL$Early[which(dfALL$Max %in% c(5:8))] <- "P3"
dfALL$Early[which(dfALL$Max %in% c(9:12))] <- "P4"
dfALL$Early[which(dfALL$Max %in% c(13:16))] <- "P6"

max_prop_early<-as.matrix((table(dfALL$Disease, dfALL$Early)))
order_df <- as.data.frame.matrix((prop.table(max_prop_early, margin = 1)))
passage_maximus <- apply(max_prop_early, 1, function(x){ head(which(x==max(x)), 1) })
dfmax <- c()
for(i in names(passage_maximus)){
  dfmax <- c(dfmax, order_df[i,passage_maximus[i]])
}

dfmaxtable <- data.frame(disease=names(passage_maximus), passage=passage_maximus, prop=dfmax)
order2 <- rownames(dfmaxtable[order(dfmaxtable$passage, 1-dfmaxtable$prop), ])

## Use this if you want first P2 and then P3
#order_df <- as.data.frame.matrix((prop.table(max_prop_early, margin = 1)))
#order2 <- rev(rownames(order_df[order(order_df$Early,order_df$Early2),]))
## Use this if you want order by P2+P3
#order2<-rev(names(sort(max_prop_early[,1]/(max_prop_early[,1]+max_prop_early[,2])))) 

max_table <- melt(table(dfALL))
max_table_prop <- melt(prop.table(table(dfALL), "Disease"))

c1<-rev(brewer.pal(n = 5, name = "Greens")[-1])
c2<-rev(brewer.pal(n = 5, name = "Blues")[-1])
c3<-rev(brewer.pal(n = 5, name = "Reds")[-1])
c4<-rev(brewer.pal(n = 5, name = "Purples")[-1])
c5<-rev(brewer.pal(n = 5, name = "Greys")[-1])

col <- c(c1,c2,c3,c4,c5)

df2plot <- left_join(max_table, dfmaxtable,by=c("Disease"="disease"))
df2plot$Passage2plot <- mapvalues(df2plot$passage, from=c(1:5), to=c("P2","P3","P4","P6","P8"))

order_mix <- order2[c(4:1, 13:5, 14:21, 22:29)]
save(order_mix, file="OrderMix_22Nov21.rda")
df2plot$Disease <- factor(df2plot$Disease, levels=(order_mix))

ggplot(df2plot, aes(y=Disease, x=value, fill=factor(Max, levels = rev(c(1:20))))) + geom_bar(position="fill", stat="identity") + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_fill_manual(values = rev(col)) + facet_grid(Passage2plot~., space = "free_y", scales = "free_y")

peaks_disease <- unique(df2plot[,c(1,7)])

  ### Check all disease datasets ----
  disgenet <- read.table("DATA/DisGENET/curated_gene_disease_associations.tsv", sep="\t", fill=TRUE, quote = "", header = TRUE)
  limit <- names(table(disgenet$diseaseId)[which(table(disgenet$diseaseId)>40)])
  filterDis <- filter(disgenet, diseaseType == "disease", diseaseId %in% limit)
  
  load("DATA/MaximumTenporal_231121.rda")
  dfALL$Gene <- rownames(dfALL)
  df_bottom <- dfALL[bottom,]
  df_bottom$Disease <- "Bottom1000"
  df_top <- dfALL[top,]
  df_top$Disease <- "Top1000"
  df_1RPKM <- dfALL[oneRPKM,]
  df_1RPKM$Disease <- ">1 RPKM"
  
  disgenet_EL<-c()
  for (i in unique(filterDis$diseaseName)){
    temp <- dfALL[which(rownames(dfALL) %in% filter(filterDis, diseaseName == i)$geneSymbol),]
    temp$Disease <- i
    disgenet_EL <- rbind(disgenet_EL, temp)
  }

  disgenet_EL <- rbind(disgenet_EL, dfALL, max_list2, df_bottom, df_top, df_1RPKM) 
  disgenet_EL$Early <- ifelse( disgenet_EL$Max %in% c(1:4), "P2", "P8" )
  disgenet_EL$Early[which(disgenet_EL$Max %in% c(5:8))] <- "P3"
  disgenet_EL$Early[which(disgenet_EL$Max %in% c(9:12))] <- "P4"
  disgenet_EL$Early[which(disgenet_EL$Max %in% c(13:16))] <- "P6"
  
  saveRDS(disgenet_EL, file="disgenet_EL.rds")
  
  max_prop_early<-as.matrix((table(disgenet_EL$Disease, disgenet_EL$Early)))
 
  order_df <- as.data.frame.matrix((prop.table(max_prop_early, margin = 1)))
  passage_maximus <- apply(max_prop_early, 1, function(x){ head(which(x==max(x)), 1) })
  dfmax <- c()
  for(i in names(passage_maximus)){
    dfmax <- c(dfmax, order_df[i,passage_maximus[i]])
  }
  
  dfmaxtable <- data.frame(disease=names(passage_maximus), passage=passage_maximus, prop=dfmax)
  dfmaxtable$Passage2plot <- mapvalues(dfmaxtable$passage, from=c(1:5), to=c("P2","P3","P4","P6","P8"))
  saveRDS(unique(dfmaxtable[,c(1,4)]), file="Passages2plot.rds")
  
  order4 <- rownames(dfmaxtable[order(dfmaxtable$passage, 1-dfmaxtable$prop), ])
 
  max_table <- melt(table(disgenet_EL[,-3]))
  
  #max_table_prop <- melt(prop.table(table(disgenet_EL), "Disease"))
  
  df2plotDis <- left_join(max_table, dfmaxtable,by=c("Disease"="disease"))
  df2plotDis$Passage2plot <- mapvalues(df2plotDis$passage, from=c(1:5), to=c("P2","P3","P4","P6","P8"))
  
  #max_table$Disease <- factor(max_table$Disease, levels=order4)
  
  c1<-rev(brewer.pal(n = 5, name = "Greens")[-1])
  c2<-rev(brewer.pal(n = 5, name = "Blues")[-1])
  c3<-rev(brewer.pal(n = 5, name = "Reds")[-1])
  c4<-rev(brewer.pal(n = 5, name = "Purples")[-1])
  c5<-rev(brewer.pal(n = 5, name = "Greys")[-1])
  
  col <- c(c1,c2,c3,c4,c5)

  order_mix_dis <- order4[c(28:1, 59:29, 60:298)]
  save(order_mix_dis, file="OrderMix_dis_22Nov21.rda")
  df2plotDis$Disease <- factor(df2plotDis$Disease, levels=(order_mix_dis))
  
  ggplot(df2plotDis, aes(y=Disease, x=value, fill=factor(Max, levels = rev(c(1:20))))) + geom_bar(position="fill", stat="identity") + theme_classic() +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
    scale_fill_manual(values = rev(col)) + facet_grid(Passage2plot~., space = "free_y", scales = "free_y")
  
  
  ggsave("AllDiseaseV6.pdf", height = 60, width = 30,units = "cm")
  
 ### FMRP CHD8 ----

load("OrderMix_dis_22Nov21.rda")  
load("DATA/DiseaseGenes/FMRP_CHDO_TOTAL.rda")  
disease_full$FMRP <- ifelse(disease_full$GENE %in% fmrp_TOTAL, 1, 0)
disease_full$CHD8 <- ifelse(disease_full$GENE %in% chd8_TOTAL, 1, 0)

frmp<-c()
for(disease in colnames(disease_full)[-c(1,17,18)]){
  if ( sum(disease_full[,disease])==0 ) { next }
  print(disease)
  mask <- which(disease_full[,disease] == 1)
  frmp <- rbind(frmp, data.frame( Disease=disease, FMRP=disease_full[mask,c("FMRP")], CHD8=disease_full[mask,c("CHD8")]))
}  

frmp <- rbind(frmp, data.frame( Disease="ALL", FMRP=disease_full[,c("FMRP")], CHD8=disease_full[,c("CHD8")]))

m <- which(disease_full$GENE %in% bottom)
frmp <- rbind(frmp, data.frame( Disease="Bottom1000", FMRP=disease_full[m,c("FMRP")], CHD8=disease_full[m,c("CHD8")]))
m <- which(disease_full$GENE %in% top)
frmp <- rbind(frmp, data.frame( Disease="Top1000", FMRP=disease_full[m,c("FMRP")], CHD8=disease_full[m,c("CHD8")]))
m <- which(disease_full$GENE %in% oneRPKM)
frmp <- rbind(frmp, data.frame( Disease=">1 RPKM", FMRP=disease_full[m,c("FMRP")], CHD8=disease_full[m,c("CHD8")]))

frmp_table <- melt(table(frmp[,c(1,2)]))
chd8_table <- melt(table(frmp[,c(1,3)]))
frmp_table <- left_join(frmp_table, peaks_disease)
chd8_table <- left_join(chd8_table, peaks_disease)

frmp_table$Disease <- factor(frmp_table$Disease, levels=order_mix)

mytable <- frmp_table
fmrp.p <- c()
fmrp.or <- c()
for (i in unique(frmp_table$Disease)){
    a <- filter(mytable, FMRP==1, Disease==i)$value
    b <- filter(mytable, FMRP==0, Disease==i)$value
    c <- filter(mytable, FMRP==1, Disease=="ALL")$value - a
    d <- filter(mytable, FMRP==0, Disease=="ALL")$value - b
    f <- fisher.test(matrix(c(a,c,b,d), nrow=2), alternative = "greater")
    fmrp.p <- c(fmrp.p, f$p.value)
    fmrp.or <- c(fmrp.or, f$estimate)
}
names(fmrp.p) <- unique(frmp_table$Disease)
names(fmrp.or) <- unique(frmp_table$Disease)
test <- data.frame(Disease = unique(frmp_table$Disease), P=fmrp.p, OR=fmrp.or, PAdj=p.adjust(fmrp.p, method = "BH"))

frmp_table <- left_join(frmp_table, test)

mytable <- chd8_table
fmrp.p <- c()
fmrp.or <- c()
for (i in unique(frmp_table$Disease)){
  a <- filter(mytable, CHD8==1, Disease==i)$value
  b <- filter(mytable, CHD8==0, Disease==i)$value
  c <- filter(mytable, CHD8==1, Disease=="ALL")$value - a
  d <- filter(mytable, CHD8==0, Disease=="ALL")$value - b
  f <- fisher.test(matrix(c(a,c,b,d), nrow=2), alternative = "greater")
  fmrp.p <- c(fmrp.p, f$p.value)
  fmrp.or <- c(fmrp.or, f$estimate)
}
names(fmrp.p) <- unique(frmp_table$Disease)
names(fmrp.or) <- unique(frmp_table$Disease)
test <- data.frame(Disease = unique(frmp_table$Disease), P=fmrp.p, OR=fmrp.or, PAdj=p.adjust(fmrp.p, method = "BH"))

chd8_table <- left_join(chd8_table, test)

p1 <- ggplot(frmp_table, aes(y=Disease, x=value, alpha=as.factor(ifelse(PAdj < 0.05,1,0)), fill=factor(FMRP))) + geom_bar(position="fill", stat="identity") + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_alpha_manual(values = c(0.4,1)) +
  scale_fill_manual(values = c("grey85","cadetblue"))+coord_cartesian(xlim=c(0,1)) + facet_grid(Passage2plot~., space = "free_y", scales = "free_y")

chd8_table$Disease <- factor(chd8_table$Disease, levels=order_mix)
p2 <- ggplot(chd8_table, aes(y=Disease, x=value, alpha=as.factor(ifelse(PAdj < 0.05,1,0)), fill=factor(CHD8))) + geom_bar(position="fill", stat="identity") + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + scale_alpha_manual(values = c(0.4,1))  + 
  scale_fill_manual(values = c("grey85","tomato")) +coord_cartesian(xlim=c(0,1))+ facet_grid(Passage2plot~., space = "free_y", scales = "free_y")

write.table(unique(chd8_table[,-c(2:4)]), file="DATA/Enrichment_CHD8.txt", sep="\t", quote=F, row.names = FALSE)
write.table(unique(frmp_table[,-c(2:4)]), file="DATA/Enrichment_FMRP.txt", sep="\t", quote=F, row.names = FALSE)

write.xlsx(unique(chd8_table[,-c(2:4)]), file="DATA/Enrichment_CHD8_FMRP.xlsx", sheetName = "CHD8", 
           col.names = TRUE, row.names = FALSE, append = TRUE)
write.xlsx(unique(frmp_table[,-c(2:4)]), file="DATA/Enrichment_CHD8_FMRP.xlsx", sheetName = "FMRP", 
           col.names = TRUE, row.names = FALSE, append = TRUE)

ggpubr::ggarrange(p1, p2, labels = c("A", "B"),
          common.legend = TRUE, legend = "bottom")

### Slope Telley ----

marques <- c( "[-1,-0.85]", "(-0.85,-0.7]", "(-0.7,-0.5]", "(-0.5,-0.25]", "(-0.25,0]", 
              "(0,0.25]","(0.25,0.5]", "(0.5,0.7]","(0.7,0.85]","(0.85,1]" )
b <- c( -1.0, -0.85,-0.7,-0.5,-0.25, 0.0,0.25,0.5, 0.7,0.85, 1.0)

load("DATA/Slopes/telley_slope_V2_Scaled.rda")
telley_slope <- telley_slope[-which(duplicated(telley_slope$GeneHuman)),]
telley_slope[,c(3:6)][telley_slope[,c(3:6)]=="NA"] <- 0

breaks <- apply(telley_slope[,-(1:2)], 2, function(x){ cut(as.numeric(x), breaks = b, include.lowest=TRUE)  })
rownames(breaks) <- (telley_slope$GeneHuman)
br <- melt(breaks)
brt <- br %>% group_by(Var2, value) %>% dplyr::summarise(n=n())
brt$Disease="ALL"

brt_bottom <- filter(br, Var1 %in% bottom) %>% group_by(Var2, value) %>% dplyr::summarise(n=n())
brt_bottom$Disease="Bottom1000"
brt_top <- filter(br, Var1 %in% top) %>% group_by(Var2, value) %>% dplyr::summarise(n=n())
brt_top$Disease="Top1000"
brt_one <- filter(br, Var1 %in% oneRPKM) %>% group_by(Var2, value) %>% dplyr::summarise(n=n())
brt_one$Disease=">1 RPKM"

brt <- rbind(brt, brt_bottom, brt_top, brt_one)

for(disease in colnames(disease_full)[-c(1,17,18)]){
  print(disease)
  mask <- which(disease_full[,disease] == 1)
  brt_temp <- filter(br, Var1 %in% disease_full$GENE[mask]) %>% group_by(Var2, value) %>% dplyr::summarise(n=n())
  brt_temp$Disease = disease
  brt <- rbind(brt, brt_temp)
}  

order_earlylate <- filter(brt, Var2=="E12") %>% group_by(Disease) %>% add_tally(n) %>% mutate(freq = n / nn) %>% filter(value %in% c("[-1,-0.85]", "(-0.85,-0.7]", "(-0.7,-0.5]")) %>% dplyr::summarise(freqtotal=sum(freq)) %>% arrange(freqtotal)

brt$Disease <- factor(brt$Disease, levels=order_earlylate$Disease)
brt$value <- factor(brt$value, levels=rev( marques ))

fblue <- colorRampPalette(c("white", "cadetblue"))
c2<-fblue(5)
fred <- colorRampPalette(c("white", "darkred"))
c3<-fred(5)
c2[1:2] <- "grey93"
c3[1:2] <- "grey93"

cols <- c(rev(c2),c3)

#brt <- left_join(brt, peaks_disease)
ggplot(brt, aes(y=Disease, x=n, fill=(value))) + geom_bar(position="fill", stat="identity") + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_fill_manual(values = (cols)) + facet_grid(.~Var2, space = "free_y", scales = "free_y")


### Slope DECON ----

load("DATA/Slopes/DeconSlopes_V3_Scaled.rda")
decon_slope[,-(1:2)][decon_slope[,-(1:2)]=="NA"] <- 0

breaks <- apply(decon_slope[,-(1:2)], 2, function(x){ cut(as.numeric(x), breaks = b, include.lowest=TRUE)  })
rownames(breaks) <- (decon_slope$GeneHuman)
br <- melt(breaks)
brt <- br %>% group_by(Var2, value) %>% dplyr::summarise(n=n())
brt$Disease="ALL"

brt_bottom <- filter(br, Var1 %in% bottom) %>% group_by(Var2, value) %>% dplyr::summarise(n=n())
brt_bottom$Disease="Bottom1000"
brt_top <- filter(br, Var1 %in% top) %>% group_by(Var2, value) %>% dplyr::summarise(n=n())
brt_top$Disease="Top1000"
brt_one <- filter(br, Var1 %in% oneRPKM) %>% group_by(Var2, value) %>% dplyr::summarise(n=n())
brt_one$Disease=">1 RPKM"

brt <- rbind(brt, brt_bottom, brt_top, brt_one)


for(disease in colnames(disease_full)[-c(1,17,18)]){
  print(disease)
  mask <- which(disease_full[,disease] == 1)
  brt_temp <- filter(br, Var1 %in% disease_full$GENE[mask]) %>% group_by(Var2, value) %>% dplyr::summarise(n=n())
  brt_temp$Disease = disease
  brt <- rbind(brt, brt_temp)
}  

brt$Disease <- factor(brt$Disease, levels=order_earlylate$Disease)
brt$value <- factor(brt$value, levels=rev(marques))

fblue <- colorRampPalette(c("white", "cadetblue"))
c2<-fblue(5)
fred <- colorRampPalette(c("white", "darkred"))
c3<-fred(5)

c2[1:2] <- "grey93"
c3[1:2] <- "grey93"

cols <- c(rev(c2), c3)

brt <- left_join(brt, peaks_disease)
ggplot(brt, aes(y=Disease, x=n, fill=(value))) + geom_bar(position="fill", stat="identity") + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_fill_manual(values = (cols), drop=FALSE) + facet_grid(.~Var2, space = "free_y", scales = "free_y")



#### Summary DV ----

load("DATA/Slopes/DorsalVentralDeSeq2.rda")
DV_FC <- data.frame(Day8=r8$log2FoldChange, Day17=r17$log2FoldChange, Day30=r30$log2FoldChange, row.names = rownames(r8))

b <- c(-Inf,-4,-2,-1,-0.5,-0.25,0,0.25,0.5,1,2,4,Inf)
breaks <- apply(DV_FC, 2, function(x){ cut(as.numeric(x), breaks = b)  })
rownames(breaks) <- rownames(DV_FC)
br <- melt(breaks)
brt <- br %>% group_by(Var2, value) %>% dplyr::summarise(n=n())
brt$Disease="ALL"

brt_bottom <- filter(br, Var1 %in% bottom) %>% group_by(Var2, value) %>% dplyr::summarise(n=n())
brt_bottom$Disease="Bottom1000"
brt_top <- filter(br, Var1 %in% top) %>% group_by(Var2, value) %>% dplyr::summarise(n=n())
brt_top$Disease="Top1000"
brt_one <- filter(br, Var1 %in% oneRPKM) %>% group_by(Var2, value) %>% dplyr::summarise(n=n())
brt_one$Disease=">1 RPKM"

brt <- rbind(brt, brt_bottom, brt_top, brt_one)

for(disease in colnames(disease_full)[-c(1,17,18)]){
  print(disease)
  mask <- which(disease_full[,disease] == 1)
  brt_temp <- filter(br, Var1 %in% disease_full$GENE[mask]) %>% group_by(Var2, value) %>% dplyr::summarise(n=n())
  brt_temp$Disease = disease
  brt <- rbind(brt, brt_temp)
}  

order_frontal <- filter(brt, Var2=="Day8", !is.na(value)) %>% group_by(Disease) %>% add_tally(n) %>% mutate(freq = n / nn) %>% filter(value %in% c("(4,7]", "(2,4]", "(1,2]", "(0.5,1]","(0.25,0.5]")) %>% dplyr::summarise(freqtotal=sum(freq)) %>% arrange(freqtotal)

brt$Disease <- factor(brt$Disease, levels=(order_frontal$Disease))
brt$value <- factor(brt$value, levels=rev(c( "(-7,-4]", "(-4,-2]", "(-2,-1]", "(-1,-0.5]", "(-0.5,-0.25]", "(-0.25,0]", "(0,0.25]","(0.25,0.5]", "(0.5,1]", "(1,2]", "(2,4]","(4,7]" ))  )

fblue <- colorRampPalette(c("white", "goldenrod2"))
c2<-fblue(6)
fred <- colorRampPalette(c("white", "forestgreen"))
c3<-fred(6)
c2[1] <- "grey93"
c3[1] <- "grey93"

cols <- c(rev(c3), c2)

ggplot(brt[!is.na(brt$value),], aes(y=Disease, x=n, fill=(value))) + geom_bar(position="fill", stat="identity") + theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + 
  scale_fill_manual(values = (cols), drop=FALSE) + facet_grid(.~Var2)




### Write genes DV biased per disease

# biased <- apply(DV_FC[,-4], 1, function(x) { sum(abs(x) > 0.25)}   )
# 
# for(disease in colnames(disease_full)[-c(1,17,18)]){
#   print(disease)
#   mask <- which(disease_full[,disease] == 1)
#   mask2 <- which(rownames(DV_FC) %in% disease_full[mask,"GENE"] & biased > 0)
#   write.table(round(DV_FC[mask2,], 3), file = paste0( "DATA/Slopes/GenesSlopeDisease/DVBias_",disease,".txt"), quote = FALSE, sep = "\t" )
# }  
#    
  
