
### PLOT MAGMA
library(ggplot2)
require(dplyr)
library("readxl")
library(plyr)
library(ComplexHeatmap)
library(reshape2)

setwd("/Users/gabrielsantpere/Desktop/Neurogenomics/Projectes/Nicola/")

whole <- read.table("MAGMA/Prepare_MAGMA/MAGMA_complete_results.txt", header=T, sep="\t")
whole$VARIABLE <- gsub("Pattern_","p", as.character(whole$VARIABLE) )

wholeD <- filter(whole, DATASET=="DISEASE")

diseases <- c( "AD_2019", "SCZ3_2020", "AN_2019", "OCD_2017", "TS", 
              "BD", "NEUROT", "IQ", "MDD", "ADHD","ASD","PD")

wholeD$DISEASE <-  mapvalues( wholeD$DISEASE, from = diseases, 
           to = c("AD_2019","SCZ_2020","AN_2019",
                  "OCD_2018","TS_2019","BD_2019", 
                  "NEUROT_2018", "IQ_2018","MDD_2018",
                  "ADHD_2019", "ASD_2019", "PD_2014"))

order <- c("ASD_2019", "ADHD_2019","SCZ_2020", "BD_2019",
           "OCD_2018", "AN_2019","TS_2019","MDD_2018","NEUROT_2018","IQ_2018", 
           "AD_2019", "PD_2014" )

wholeD <- filter(wholeD, DISEASE %in% order)

wholeD <- wholeD %>% group_by(TEST) %>% dplyr::mutate(P.adj = p.adjust(P, method = "bonferroni"))
wholeD <- wholeD %>% group_by(TEST, DISEASE) %>% dplyr::mutate(P.adjVar = p.adjust(P, method = "bonferroni"))

diseaseDataset <- data.frame(PATTERN = wholeD$VARIABLE, DISEASE = wholeD$DISEASE, P = wholeD$P, P.adj = wholeD$P.adj, P.adj2 = wholeD$P.adjVar, TEST = wholeD$TEST, DATASET = "DISEASE" )

#### Input genes for cortical malformations

raw  <- read.table("DATA/PValuesGenSetTest/CoGAPSpatternEnrichment.NEGlog10pValues_Raw.txt", header=T, sep="\t")
rownames(raw) <- raw$X
rawP <- 10^(-raw[,-1])
rawWhole <- cbind(raw$X, rawP)
colnames(rawWhole)[1] <- "X"  
carlo <- melt(rawWhole[,c(1:9)])
carlo$TEST <- "Early"
carlo$TEST[grep("AZ", carlo$X)] <- "Late"
carlo$TEST[grep("DVlines", carlo$X)] <- "Dorsal_Ventral"

carlo$VARIABLE <- sapply(strsplit(as.character(carlo$X), split="_"),"[" , 2)

carlo <- carlo %>% group_by(TEST) %>% dplyr::mutate (P.adj = p.adjust(value, method = "bonferroni"))
carlo <- carlo %>% group_by(TEST, variable) %>% dplyr::mutate (P.adjVar = p.adjust(value, method = "bonferroni"))

MalformationDataset <- data.frame(PATTERN = carlo$VARIABLE, DISEASE = carlo$variable, P = carlo$value, P.adj = carlo$P.adj, P.adj2 = carlo$P.adjVar, TEST = carlo$TEST, DATASET = "Cortical malformation" )

## Input carlos tests

disease2 <- read.table("DATA/CarloEnrichments/CoGAPSpatternEnrichmentMAS.NEGlog10pValues_Raw.txt", header = T)
disease2 <- disease2[,c(1:11)]
rownames(disease2) <- disease2$X
disease2P <- 10^(-disease2[,-1])
rawWhole2 <- cbind(disease2$X, disease2P)
colnames(rawWhole2)[1] <- "X" 
carlo2 <- melt(rawWhole2)
carlo2$TEST <- "Early"
carlo2$TEST[grep("AZ", carlo2$X)] <- "Late"
carlo2$TEST[grep("DVlines", carlo2$X)] <- "Dorsal_Ventral"

carlo2$VARIABLE <- sapply(strsplit(as.character(carlo2$X), split="_"),"[" , 2)

carlo2 <- carlo2 %>% group_by(TEST) %>% dplyr::mutate (P.adj = p.adjust(value, method = "bonferroni"))
carlo2 <- carlo2 %>% group_by(TEST, variable) %>% dplyr::mutate (P.adjVar = p.adjust(value, method = "bonferroni"))

diseaseDataset2 <- data.frame(PATTERN = carlo2$VARIABLE, DISEASE = carlo2$variable, P = carlo2$value, P.adj = carlo2$P.adj, P.adj2 = carlo2$P.adjVar, TEST = carlo2$TEST, DATASET = "GWAS Hits" )


### Join datasets

df <- rbind(diseaseDataset, MalformationDataset, diseaseDataset2)
df$label <- ifelse(df$P < 0.05, "\U2022","")
df$label[df$P.adj2 < 0.05] <- "\U2020"
df$label[df$P.adj < 0.05] <- "\U2021"

df$PATTERN <- factor(df$PATTERN, levels= rev(paste0("p",c(1:30))))
df$TEST <- factor(df$TEST, levels=c( "Early", "Late", "Dorsal_Ventral" ))

p1 <- ggplot(df, aes(y = PATTERN, x = DISEASE, fill = -log10(P))) + geom_tile( width=0.95, height=0.95, size=0.2, color="darkgrey") + 
  theme(axis.text=element_text(size=6))+ 
  theme(panel.border = element_blank(), panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank(), axis.line = element_line(colour = "white"), 
        panel.background = element_blank()) + 
  scale_fill_gradient(low = "white", high = "tomato") + theme(axis.text.x = element_text(angle = 45, hjust = 1)) + 
  geom_text(aes(label = label), color="black", size=2) + scale_color_manual(values = c("black", "grey")) +
  facet_grid(TEST~DATASET, scales =  "free", space = "free")

ggsave(p1, filename = "CompleteDisease.pdf")
