

### Prepare SFARI

require(dplyr)

setwd("/Users/gabrielsantpere/Desktop/Neurogenomics/Projectes/Nicola/")

sfari <- read.csv("../DATABASES/SFARI-Gene_genes_01-13-2021release_02-17-2021export.csv", header=TRUE)

table(sfari$gene.score, sfari$syndromic)

write.table(filter(sfari, syndromic==1)$gene.symbol, file="DATA/DiseaseGenes/SFARI_Syndromic.txt", row.names = F, col.names = F, quote = F)
write.table(filter(sfari, gene.score==1)$gene.symbol, file="DATA/DiseaseGenes/SFARI_Score1.txt", row.names = F, col.names = F, quote = F)
write.table(filter(sfari, gene.score==2)$gene.symbol, file="DATA/DiseaseGenes/SFARI_Score2.txt", row.names = F, col.names = F, quote = F)
write.table(filter(sfari, gene.score==3)$gene.symbol, file="DATA/DiseaseGenes/SFARI_Score3.txt", row.names = F, col.names = F, quote = F)

