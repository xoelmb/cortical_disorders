
#### Calculate slopes for Telley et al. and DECON

require(dplyr)
library("readxl")
library(plyr)
library(biomaRt)

setwd("/Users/gabrielsantpere/Desktop/Neurogenomics/Projectes/Nicola/")

## 1- Use Biomart to convert human and mouse genes

convertMouseGeneList <- function(x){
  human = useMart("ensembl", dataset = "hsapiens_gene_ensembl")
  mouse = useMart("ensembl", dataset = "mmusculus_gene_ensembl")
  
  genesV2 = getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = x , mart = mouse, attributesL = c("hgnc_symbol"), martL = human, uniqueRows=T)
  humanx <- unique(genesV2[,])
  
  # Print the first 6 genes found to the screen
  print(head(humanx))
  return(humanx)
  }

## 2- Calculate slopes in DECON

load("DATA/Slopes/DeCon_w_Hs_GeneSymbol.rdata")

decon_s <- p.DeCon.log2 

scp<-c()
cpn<-c()
cth<-c()
for(i in rownames(decon_s)){
  t <- data.frame(time=c(1,1,2,2,3,3,4,4), cpn=unlist(decon_s[i,grep("_cpn", colnames(decon_s))]),
                                   cthl=unlist(decon_s[i,grep("_corticothal", colnames(decon_s))]),
                                   scp=unlist(decon_s[i,grep("_subcereb", colnames(decon_s))]))

  if(length(unique(t$scp)) > 1) { summary(lm(scale(t$scp)~t$time)) ->a } else { 
    a<-c()
    a$coefficients <- data.frame(coefficients=1) 
  }
  scp <- c(scp,ifelse(nrow(a$coefficients)>1,a$coefficients[2,1],"NA"))  
  if(length(unique(t$cthl)) > 1) { summary(lm(scale(t$cthl)~t$time)) ->a } else { 
    a<-c()
    a$coefficients <- data.frame(coefficients=1) 
  }
  cth <- c(cth,ifelse(nrow(a$coefficients)>1,a$coefficients[2,1],"NA"))  
  
  if(length(unique(t$cpn)) > 1) {summary(lm(scale(t$cpn)~t$time)) ->a} else { 
    a<-c()
    a$coefficients <- data.frame(coefficients=1) 
  }
  cpn <- c(cpn,ifelse(nrow(a$coefficients)>1,a$coefficients[2,1],"NA"))   
}  

GeneHuman<-mapvalues( rownames(decon_s), from = map.DeCon$GeneSymbol, to =  map.DeCon$`tbl$HGNC.symbol.human`)
decon_slope <- data.frame(GeneOriginal=rownames(decon_s), GeneHuman=GeneHuman, CPN=cpn,SCP=scp,CTH=cth )

save(decon_slope, file = "DATA/Slopes/DeconSlopes_V3_Scaled.rda")


## 2- Calculate slopes in Telley

load("DATA/Slopes/CtxDevoTempPatApicProgMm_ONLY_CPMlog2.cellAGEpupAGE.means.RData")

#telley<-t(scale(t(CPMlog2.cellAGEpupAGE.means)))
telley <- CPMlog2.cellAGEpupAGE.means

e12<-c()
e13<-c()
e14<-c()
e15<-c()
for(i in rownames(telley)){
  t <- data.frame(time=c(1,2,3), E12=unlist(telley[i,paste0("E12", c(".1H",".24H",".96H"))]),
                                 E13=unlist(telley[i,paste0("E13", c(".1H",".24H",".96H"))]),
                                 E14=unlist(telley[i,paste0("E14", c(".1H",".24H",".96H"))]),
                                 E15=unlist(telley[i,paste0("E15", c(".1H",".24H",".96H"))]))
  
  if(length(unique(t$E12)) > 1) {summary(lm(scale(t$E12)~t$time)) ->a} else { 
    a<-c()
    a$coefficients <- data.frame(coefficients=1) 
  }
  e12 <- c(e12,ifelse(nrow(a$coefficients)>1,a$coefficients[2,1],"NA"))  
  if(length(unique(t$E13)) > 1) {summary(lm(scale(t$E13)~t$time)) ->a} else { 
    a<-c()
    a$coefficients <- data.frame(coefficients=1) 
  }
  e13 <- c(e13,ifelse(nrow(a$coefficients)>1,a$coefficients[2,1],"NA"))  
  if(length(unique(t$E14)) > 1) {summary(lm(scale(t$E14)~t$time)) ->a} else { 
    a<-c()
    a$coefficients <- data.frame(coefficients=1) 
  }
  e14 <- c(e14,ifelse(nrow(a$coefficients)>1,a$coefficients[2,1],"NA"))   
  if(length(unique(t$E15)) > 1) {summary(lm(scale(t$E15)~t$time)) ->a} else { 
    a<-c()
    a$coefficients <- data.frame(coefficients=1) 
  }
  e15 <- c(e15,ifelse(nrow(a$coefficients)>1,a$coefficients[2,1],"NA")) 
}  

telley_slope <- data.frame(GeneOriginal=rownames(telley), E12=e12,E13=e13,E14=e14,E15=e15)
telley_slope<-left_join(telley_slope, GENEinfo[,c(1,5)], by=c("GeneOriginal" = "genes"))
telley_slope<-telley_slope[,c(1,6,2:5)]
colnames(telley_slope)[2]<-c("GeneHuman")

save(telley_slope, file = "DATA/Slopes/telley_slope_V2_Scaled.rda")


### Calculate FC in DV

load("DATA/CoGAPS_Expression/NicoDV_ExprsLog2_forGabrielHeatmaps.rdata")
colnames(NicoDV.ExprsLog2)=NicoDV.sampleMeta$ids

samples <- ifelse(grepl("2053.6|2075", NicoDV.sampleMeta$ids), "dorsal", "ventral")

dv_d <- dv[,grep("2053.6|2075", colnames(dv))]
dv_v <- dv[,grep("2053.2|2063", colnames(dv))]

day8<-NicoDV.ExprsLog2[,grep("day8", colnames(dv))]
day17<-NicoDV.ExprsLog2[,grep("day17", colnames(dv))]
day30<-NicoDV.ExprsLog2[,grep("day30", colnames(dv))]
day32<-NicoDV.ExprsLog2[,grep("day32", colnames(dv))]

d8FC <- apply(day8, 1, function(x) { mean(x[ grep("2053.6|2075", colnames(day8)) ]) -  mean(x[ grep("2053.2|2063", colnames(day8)) ] )  } )
d17FC <- apply(day17, 1, function(x) { mean(x[ grep("2053.6|2075", colnames(day17)) ]) -  mean(x[ grep("2053.2|2063", colnames(day17)) ] )  } )
d30FC <- apply(day30, 1, function(x) { mean(x[ grep("2053.6|2075", colnames(day30)) ]) -  mean(x[ grep("2053.2|2063", colnames(day30)) ] )  } )
d32FC <- apply(day32, 1, function(x) { mean(x[ grep("2053.6|2075", colnames(day32)) ]) -  mean(x[ grep("2053.2|2063", colnames(day32)) ] )  } )

DV_FC <- data.frame(Day8=d8FC, Day17=d17FC, Day30=d30FC, Day32=d32FC, row.names = rownames(NicoDV.ExprsLog2))

save(DV_FC, file = "DATA/Slopes/DV_slope.rda")


d8FC <- apply(day8, 1, function(x) { mean(x[ grep("2075", colnames(day8)) ]) -  mean(x[ grep("2063", colnames(day8)) ] )  } )
d17FC <- apply(day17, 1, function(x) { mean(x[ grep("2075", colnames(day17)) ]) -  mean(x[ grep("2063", colnames(day17)) ] )  } )
d30FC <- apply(day30, 1, function(x) { mean(x[ grep("2075", colnames(day30)) ]) -  mean(x[ grep("2063", colnames(day30)) ] )  } )
d32FC <- apply(day32, 1, function(x) { mean(x[ grep("2075", colnames(day32)) ]) -  mean(x[ grep("2063", colnames(day32)) ] )  } )

DV_FC <- data.frame(Day8=d8FC, Day17=d17FC, Day30=d30FC, Day32=d32FC, row.names = rownames(NicoDV.ExprsLog2))

save(DV_FC, file = "DATA/Slopes/DV_slope_no2053.rda")

