setwd('/users/genomics/xoel/codebases/co_new/results/')
plotdir <- 'Plots/'
dir.create(plotdir, showWarnings = F)

library(ggplot2)
library(ggh4x)
library(ggpubr)
library(ggrepel)
library(tidyr)
library(dplyr)
fig <- function(width, heigth){
 options(repr.plot.width = width, repr.plot.height = heigth)
 }

CellTypeAnnot <- c(
    # 'RGC', 'vRG', 'vtRG', 'oRG', 'RG E', 'RG L',
    'vRG E', 'vRG L', 'tRG', 'oRG E', 'oRG L', 
    
    # 'Neural', 'GluN',
    'nIPC', 'Neu E',
    
    'GluN3', 'GluN1', 'GluN2', 'GluN7',
    'GluN5', 'GluN4',
    'GluN6', 'GluN8',
    
    # 'Glial', 'Diff Glia',
    'mGPC', 'Astro', 'OPC')


CellTypeOrder <- c(
    
    'RGC', 'vRG', 'vtRG', 'oRG', 'RG E', 'RG L',
    'vRG E', 'vRG L', 'tRG', 'oRG E', 'oRG L', 
    
    'Neural', 'GluN',
    'nIPC', 'Neu E',
    
    'GluN3', 'GluN1', 'GluN2', 'GluN7',
    'GluN5', 'GluN4',
    'GluN6', 'GluN8',
    
    'Glial', 'Diff Glia',
    'mGPC', 'Astro', 'OPC')

CellTypeLabels <- CellTypeOrder

CellTypeColor <- setNames(
    c(
        '#bec1d4','#7d87b9','#023fa5','#d6bcc0','#bb7784','#8e063b','#b5bbe3',
        '#8595e1','#4a6fe3','#e6afb9','#e07b91','#d33f6a','#11c638','#c6dec7',
        '#8dd593','#ead3c6','#f0b98d','#ef9708','#0fcfc0','#9cded6','#d5eae7',
        '#f3e1eb','#f6c4e1','#f79cd4','#7f7f7f','#c7c7c7','#1CE6FF','#336600'
    )[1:length(CellTypeOrder)],
    CellTypeOrder)

ctinfo <- data.frame(
    hex=CellTypeColor,
    Name=factor(names(CellTypeColor))
) %>% mutate(`Common label`=Name)

plot.celltype <- ggplot(data = ctinfo, mapping = aes(x = 0, y=Name, fill=`Common label`)) + geom_tile() + scale_fill_manual(values=CellTypeColor)
plot.celltype

write.csv(ctinfo, file = 'CellTypeColorGroups.csv', row.names=FALSE)
write.csv(ctinfo[,c('Name','hex')], file = 'CellTypeColor.csv', row.names=FALSE)

# Orders
RegulonDiseaseOrder <- c(
  'Polymicrogyria',
  'RareMCD',
  'MDD_2018',
  'ASD_2019',
  'ASD_HC65',
  'SFAR_Synd',
  'SFAR_S1',
  'FCDandmTOR',
  'DevDyslexia',
  'SCZ_2020'
)

DiseaseOrder <- c(
  'Microcephaly',
  # 'ASD_2019',
  'Hydrocephaly',
  'Cobblestone',
  'AN_2019',
  'AD_2019',
  'ADHD_2019',
  'BD_2019',
  'NEUROT_2018',
  'IQ_2018',
  'SCZ_2020',
  'MDD_2018',
  'Polymicrogyria',
  'Lissencephaly',
  'Heterotopia',
  # 'ASD_HC65',
  # 'SFAR_Synd',
  'DD',
  # 'SFAR_S1',
  'PD_2014',
  # 'SFAR_S2',
  # 'SFAR_S3',
  'RareMCD',
  'FCDandmTOR',
  'DevDyslexia',
  'ASD'
)



DiseaseColor <- c("Microcephaly"="#1f77b4", "Hydrocephaly"="#aec7e8", "Cobblestone"="#ff7f0e", "AN_2019"="#ffbb78", "AD_2019"="#2ca02c", "ADHD_2019"="#98df8a", "BD_2019"="#d62728", "NEUROT_2018"="#ff9896", "IQ_2018"="#9467bd", "SCZ_2020"="#c5b0d5", "MDD_2018"="#8c564b", "Polymicrogyria"="#c49c94", "Lissencephaly"="#e377c2", "Heterotopia"="#f7b6d2", "DD"="#7f7f7f", "PD_2014"="#c7c7c7", "RareMCD"="#bcbd22", "FCDandmTOR"="#dbdb8d", "DevDyslexia"="#17becf", "ASD"="#9edae5")

write.csv(as.data.frame(DiseaseColor), file = 'DiseaseColor.csv')
DiseaseExclude <- c('CHD8', 'FMRP')

RoleOrder <- c(
  'none',
  "Ultra peripheral",
  "Peripheral",
  "Connector",
  "Kinless" ,
  "Provincial Hub",
  "Connector Hub",
  "Kinless Hub"
)

RoleColors <- c(
  'none'='black',
  "Ultra peripheral" = '#7B7578',
  "Peripheral" = '#F39D85',
  "Connector" = '#B8DC9F',
  "Kinless" = '#8F95C7',
  "Provincial Hub" = '#F7F3A9',
  "Connector Hub" = '#D9C3C2',
  'Kinless Hub' = '#EBEBEB'
)

# INCLUDE INCLUDE INCLUDE
SampleLabels <- c(
    'RGCmaturation' = 'NEP maturation',
    # 'RGCmaturationExtended' = 'NEP maturation',
    'NeuralPCW16' = 'Neurogenesis',
    'NeuralPCW20' = 'Neurogenesis',
    'NeuralPCW21' = 'Neurogenesis',
    'NeuralPCW24' = 'Neurogenesis',
    'Gliogenesis' = 'Gliogenesis'
    # 'GliogenesisExtended' = 'Gliogenesis'
)

SampleOrder <- names(SampleLabels)

SampleLabelOrder <- c('NEP maturation',
                      'Neurogenesis',
                      'Gliogenesis')
# SampleInclude <- c('NeuralPCW20', 'Gliogenesis', 'RGCmaturation')
SampleInclude <- c('NeuralPCW20', 'Gliogenesis')

# LABELS LABELS LABELS
ScoreInclude <- c(
  'degree_all',
  'degree_centrality_all',
  'degree_in',
  'degree_centrality_in',
  'degree_out',
  'degree_centrality_out',
  'betweenness_centrality',
  'eigenvector_centrality',
  'connectivity',
  'participation'
)

ScoreLabels <- c(
  'degree_all' = 'Degree (all)',
  'degree_centrality_all' = 'Degree centrality (all)',
  'degree_in' = 'Degree (in)',
  'degree_centrality_in' = 'Degree centrality (in)',
  'degree_out' = 'Degree (out)',
  'degree_centrality_out' = 'Degree centrality (out)',
  'betweenness_centrality' = 'Betweenness centrality',
  'eigenvector_centrality' = 'Eigenvector centrality',
  'connectivity' = 'Connectivity',
  'participation' = 'Participation'
)

################################################
###################    TFs   ###################
################################################
TranscriptionFactors <- read.csv('~/codebases/cortical_disorders2/raw/HumanTFs/DatabaseExtract_v_1.01.csv', row.names = 1)
tf_genes <- unique(gsub(' ', '', split(TranscriptionFactors, TranscriptionFactors$Is.TF.)[['Yes']]$HGNC.symbol, fixed=T))

################################################
################### REGULONS ###################
################################################

# Read Regulon table
# Regulons <- read.csv('../../cortical_disorders/DATA/Regulons/Diseases/report.geneset=DISEASES.method=iCisTarget.motif_search_range=long.directAnnotation=TRUE.motif_NES=3.disease_list=reduced.motif_db=v10.csv',
DisRegulons.long <- read.csv('~/codebases/cortical_disorders2/results/RCT_diseases/cisTarget_regulons.csv')
DisRegulons.long$Disease <- factor(DisRegulons.long$geneSet, DiseaseOrder)

# Get Disease Regulon genes from Regulons
# Regulons.long <- separate_rows(Regulons, Targets, sep = ', ')
regulon_list_per_disease <- split(DisRegulons.long, DisRegulons.long$Disease)

gene_list_per_disregulon <- lapply(regulon_list_per_disease,
                                function(x){
                                  unique(c(x$Core, x$Target))
                                })
disregulon_genes <- unique(unlist(gene_list_per_disregulon))

# Read Regulon table
# Regulons <- read.csv('../../cortical_disorders/DATA/Regulons/Diseases/report.geneset=DISEASES.method=iCisTarget.motif_search_range=long.directAnnotation=TRUE.motif_NES=3.disease_list=reduced.motif_db=v10.csv',
PeakRegulons.long <- read.csv('~/codebases/cortical_disorders2/results/RCT_peaks/cisTarget_regulons.csv')
PeakRegulons.long$Peak <- factor(PeakRegulons.long$geneSet, unique(PeakRegulons.long$geneSet))

# Get Disease Regulon genes from Regulons
# Regulons.long <- separate_rows(Regulons, Targets, sep = ', ')
regulon_list_per_peak <- split(PeakRegulons.long, PeakRegulons.long$Peak)

gene_list_per_peakregulon <- lapply(regulon_list_per_peak,
                                function(x){
                                  unique(c(x$Core, x$Target))
                                })
peakregulon_genes <- unique(unlist(gene_list_per_peakregulon))

regulon_genes <- unique(c(disregulon_genes, peakregulon_genes))
regulon_cores <- unique(c(DisRegulons.long$Core, PeakRegulons.long$Core))
################################################
################### DISEASES ###################
################################################
# Diseases <- read.csv('../../cortical_disorders/DATA/Regulons/Diseases/Disease_table.csv', row.names=1)
Diseases <- read.csv('~/codebases/cortical_disorders2/data/gene_disorder_associations.original_genes.csv', row.names=1)
Diseases <- Diseases[,!colnames(Diseases) %in% DiseaseExclude]

gene_list_per_disease <- lapply(Diseases, function(x){
  y <- unique(ifelse(x==1, rownames(Diseases), NA))
  y[!is.na(y)]
})
gene_list_per_disease <- gene_list_per_disease[sapply(gene_list_per_disease, length)!=0]
disease_genes <- unique(unlist(gene_list_per_disease))

################################################
###################  SCORES  ###################
################################################

# Read Network Scores
NetworkScores <- readxl::read_excel('NetworkScores.xlsx')#, col_types='text')

NetworkScores$cell.type <- factor(NetworkScores$cell.type,
                                  levels = CellTypeOrder)
NetworkScores$SampleLabel <- factor(SampleLabels[NetworkScores$Sample], levels=SampleLabelOrder)

# Split/duplicate scores by disease
NetworkScoresDiseases <- list()
for (dis in names(gene_list_per_disease)){
  dis_genes <- gene_list_per_disease[[dis]]
  dis_scores <- NetworkScores[NetworkScores$Gene %in% dis_genes,]
  dis_scores$Disease <- dis
  NetworkScoresDiseases[[dis]] <- dis_scores
}
NetworkScoresDiseases <- do.call('rbind', NetworkScoresDiseases)
NetworkScoresDiseases$Disease <- factor(NetworkScoresDiseases$Disease, DiseaseOrder)

################################################
###############  PERTURBATIONS  ################
################################################

options(warn=-1)

# Perturbation Summaries
Perturbations <- readxl::read_excel('PerturbationSummary.xlsx')

Perturbations$log2_ratio <- as.numeric(Perturbations$log2_ratio)
Perturbations$cell.type <- factor(Perturbations$cell.type,
                                  levels = CellTypeOrder)
Perturbations$SampleLabel <- factor(SampleLabels[Perturbations$Sample], levels=SampleLabelOrder)

# pert.dis.annot <- Diseases[Perturbations$Gene,-1]
# Perturbations <- as.data.frame(cbind(Perturbations, pert.dis.annot))


# Split/duplicate scores by disease
PerturbationsDiseases <- list()
for (dis in names(gene_list_per_disease)){
  dis_genes <- gene_list_per_disease[[dis]]
  dis_pert <- Perturbations[Perturbations$Gene %in% dis_genes,]
  if (nrow(dis_pert) < 1){
    next
  }

  dis_pert$Disease <- dis
  PerturbationsDiseases[[dis]] <- dis_pert
}
PerturbationsDiseases <- do.call('rbind', PerturbationsDiseases)
PerturbationsDiseases$Disease <- factor(PerturbationsDiseases$Disease, DiseaseOrder)


################################################
###############  PERTURBATIONS  ################
################################################


# Perturbation Summaries
PerturbationsCT <- readxl::read_excel('PerturbationCT2CT.xlsx')

(unique(PerturbationsCT$CTF) %in% CellTypeOrder)
(unique(PerturbationsCT$CTO) %in% CellTypeOrder)

PerturbationsCT[c('trans.cells', 'init.cells', 'trans.pct', 'Exp')] <- apply(PerturbationsCT[c('trans.cells', 'init.cells', 'trans.pct', 'Exp')], 2, as.numeric)
PerturbationsCT$CTO <- factor(PerturbationsCT$CTO,
                                  levels = CellTypeOrder)
PerturbationsCT$CTF <- factor(PerturbationsCT$CTF,
                                  levels = CellTypeOrder)
PerturbationsCT$SampleLabel <- factor(SampleLabels[PerturbationsCT$Sample], levels=SampleLabelOrder)
########################################


# Perturbation Summaries
PerturbationsCTNJ <- readxl::read_excel('PerturbationCT2CTNJ.xlsx')

(unique(PerturbationsCTNJ$CTF) %in% CellTypeOrder)
(unique(PerturbationsCTNJ$CTO) %in% CellTypeOrder)

PerturbationsCTNJ[c('trans.cells', 'init.cells', 'trans.pct', 'Exp')] <- apply(PerturbationsCTNJ[c('trans.cells', 'init.cells', 'trans.pct', 'Exp')], 2, as.numeric)
PerturbationsCTNJ$CTO <- factor(PerturbationsCTNJ$CTO,
                                  levels = CellTypeOrder)
PerturbationsCTNJ$CTF <- factor(PerturbationsCTNJ$CTF,
                                  levels = CellTypeOrder)
PerturbationsCTNJ$SampleLabel <- factor(SampleLabels[PerturbationsCTNJ$Sample], levels=SampleLabelOrder)

# pert.ct.dis.annot <- Diseases[PerturbationsCT$Gene,-1]
# PerturbationsCT <- as.data.frame(cbind(PerturbationsCT, pert.ct.dis.annot))


# Split/duplicate scores by disease
PerturbationsCTDiseases <- list()
for (dis in names(gene_list_per_disease)){
  dis_genes <- gene_list_per_disease[[dis]]
  dis_pert <- PerturbationsCT[PerturbationsCT$Gene %in% dis_genes,]
  if (nrow(dis_pert) < 1){
    next
  }

  dis_pert$Disease <- dis
  PerturbationsCTDiseases[[dis]] <- dis_pert
}
PerturbationsCTDiseases <- do.call('rbind', PerturbationsCTDiseases)
PerturbationsCTDiseases$Disease <- factor(PerturbationsCTDiseases$Disease, DiseaseOrder)

options(warn=0)

################################################
###############     NETWORKS    ################
################################################


# Network links
Links <- data.table::fread('LinkFiltering.txt')
# (unique(Links$cell.type) %in% CellTypeOrder)

Links$cell.type <- factor(Links$cell.type,
                          levels = CellTypeOrder)
Links$Sample <- factor(Links$Sample, levels=SampleOrder)

Links2 <- Links %>% group_by(Sample, cell.type) %>% mutate(
  p_rank = rank(p, na.last=T),
  p_rank_fr = rank(p, na.last=T)/length(p)
  )

Links2[,c('source.in.reg', 'target.in.reg')] <- apply(
  Links2[,c('source', 'target')], 2, 
  function(x) {x %in% regulon_genes})

Links2$any.in.reg <- apply(Links2[,c('source.in.reg', 'target.in.reg')], 1, any)

Links2[,c('source.is.tf', 'target.is.tf')] <- apply(
  Links2[,c('source', 'target')], 2, 
  function(x) {x %in% tf_genes})
Links2$any.is.tf <- apply(Links2[,c('source.is.tf', 'target.is.tf')],1, any)

Links2[,c('source.is.dis', 'target.is.dis')] <- apply(
  Links2[,c('source', 'target')], 2, 
  function(x) {x %in% disease_genes})
Links2$any.is.dis <- apply(Links2[,c('source.is.dis', 'target.is.dis')],1, any)

### Peaks

Peaks <- read.csv('~/codebases/cortical_disorders2/data/nico_peak_expression.csv', row.names=1)

gene.peak.order <- rownames(Peaks %>% arrange(desc(.)))
peak.per.gene <- apply(Peaks, 1, function(x){colnames(Peaks)[which.max(x)]})


# Tables to subscript
DisTargets <- lapply(split(DisRegulons.long$Target, DisRegulons.long$geneSet), unique)
distargets <- unique(unlist(DisTargets))
DisTargets <- data.frame(lapply(DisTargets, function(x){distargets%in%x}), row.names=distargets)

DisCores <- lapply(split(DisRegulons.long$Core, DisRegulons.long$geneSet), unique)
discores <- unique(unlist(DisCores))
DisCores <- data.frame(lapply(DisCores, function(x){discores%in%x}), row.names=discores)

PeakTargets <- lapply(split(PeakRegulons.long$Target, PeakRegulons.long$geneSet), unique)
peaktargets <- unique(unlist(PeakTargets))
PeakTargets <- data.frame(lapply(PeakTargets, function(x){peaktargets%in%x}), row.names=peaktargets)


PeakCores <- lapply(split(PeakRegulons.long$Core, PeakRegulons.long$geneSet), unique)
peakcores <- unique(unlist(PeakCores))
PeakCores <- data.frame(lapply(PeakCores, function(x){peakcores%in%x}), row.names=peakcores)

## Confusion Disease

dis.conf.cols <- c(
    'None'='white',
    'Risk'='grey',
    'Target'='orange',
    'Core'='darkred',
    'Both'='darkred'
)

get.conf.dis.df <- function(genes){


    conf.risk <- Diseases[genes,]
    rownames(conf.risk) <- genes
    conf.risk[is.na(conf.risk)] <- 0


    conf.distargets <- DisTargets[genes,]
    rownames(conf.distargets) <- genes
    conf.distargets[is.na(conf.distargets)] <- 0
    conf.distargets[,colnames(conf.risk)[!colnames(conf.risk) %in% colnames(conf.distargets)]] <- 0 
    conf.distargets <- conf.distargets[,colnames(conf.risk)]

    conf.discores <- DisCores[genes,]
    rownames(conf.discores) <- genes
    conf.discores[is.na(conf.discores)] <- 0
    conf.discores[,colnames(conf.risk)[!colnames(conf.risk) %in% colnames(conf.discores)]] <- 0 
    conf.discores <- conf.discores[,colnames(conf.risk)]



    conf.risk <- conf.risk %>% mutate(Gene=rownames(.)) %>% reshape2::melt(id.vars='Gene', value.name='Risk', variable.name='Disease')# %>% mutate(dis=as.character(Disease))
    conf.discores <- conf.discores %>% mutate(Gene=rownames(.)) %>% reshape2::melt(id.vars='Gene', value.name='Core', variable.name='Disease') #%>% mutate(Disease=factor(as.character(Disease), levels=levels(conf.risk$Disease)), dis=as.character(Disease))
    conf.distargets <- conf.distargets %>% mutate(Gene=rownames(.)) %>% reshape2::melt(id.vars='Gene', value.name='Target', variable.name='Disease') #%>% mutate(Disease=factor(as.character(Disease), levels=levels(conf.risk$Disease)), dis=as.character(Disease))



    conf.dis.df <- conf.risk %>% merge(y=conf.distargets, by=c('Gene', 'Disease'), all=T) %>%
        merge(y=conf.discores, by=c('Gene', 'Disease'), all=T)
    conf.dis.df[is.na(conf.dis.df)] <- 0

    conf.dis.df$Value <- apply(conf.dis.df, 1, function(x){
        x <- x[-(1:2)]
        x <- setNames(as.logical(as.numeric(x)), names(x))
        if (all(x)) {'Both'} else if (x['Core']){'Core'} else if (x['Target']){'Target'} else if (x['Risk']){'Risk'} else 'None'

    })

    return(conf.dis.df)
}

do.conf.dis.plot <- function(genes){

    conf.dis.df <- get.conf.dis.df(genes)
    conf.dis.df$Gene <- factor(conf.dis.df$Gene, 
                               levels=rev(genes))
    
    
    conf.dis.plot <- ggplot(conf.dis.df, aes(x=Disease, y=Gene, fill=factor(Value, names(dis.conf.cols)))) +

        geom_tile(color='grey', linewidth=0.1) +
        scale_fill_manual('Disease association', values=dis.conf.cols) +

        theme_pubr(legend = 'bottom') +
        labs_pubr(base_family = 'ArialMT') +

        theme(
            panel.spacing = unit(0.2, "lines"),
            strip.background = element_blank(),
            # strip.placement = 'outside',
            panel.background = element_blank(),
            axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5),
            axis.text.y = element_text(),
            # legend.box="vertical",
            # legend.text = element_text(angle=90, hjust=1, vjust=0.5)
        )
    
    return(conf.dis.plot)
    }

## Confusion Peak

peak.conf.cols <- c(
    'None'='white',
    'Peak'='grey',
    'Target'='orange',
    'Core'='darkred',
    'Both'='darkred'
)

get.conf.peak.df <- function(genes){
    conf.peak <- Peaks[genes,]
    rownames(conf.peak) <- genes
    conf.peak[is.na(conf.peak)] <- 0



    conf.peaktargets <- PeakTargets[genes,]
    rownames(conf.peaktargets) <- genes
    conf.peaktargets[is.na(conf.peaktargets)] <- 0
    conf.peaktargets[,colnames(conf.peak)[!colnames(conf.peak) %in% colnames(conf.peaktargets)]] <- 0 
    conf.peaktargets <- conf.peaktargets[,colnames(conf.peak)]

    conf.peakcores <- PeakCores[genes,]
    rownames(conf.peakcores) <- genes
    conf.peakcores[is.na(conf.peakcores)] <- 0
    conf.peakcores[,colnames(conf.peak)[!colnames(conf.peak) %in% colnames(conf.peakcores)]] <- 0 
    conf.peakcores <- conf.peakcores[,colnames(conf.peak)]


    conf.peak <- conf.peak %>% mutate(Gene=rownames(.)) %>% reshape2::melt(id.vars='Gene', value.name='Present', variable.name='Peak')# %>% mutate(peak=as.character(Peak))
    conf.peakcores <- conf.peakcores %>% mutate(Gene=rownames(.)) %>% reshape2::melt(id.vars='Gene', value.name='Core', variable.name='Peak')# %>% mutate(Peak=factor(as.character(Peak), levels=levels(conf.peak$Peak)), peak=as.character(Peak))
    conf.peaktargets <- conf.peaktargets %>% mutate(Gene=rownames(.)) %>% reshape2::melt(id.vars='Gene', value.name='Target', variable.name='Peak')# %>% mutate(Peak=factor(as.character(Peak), levels=levels(conf.peak$Peak)), peak=as.character(Peak))



    conf.peak.df <- conf.peak %>% merge(y=conf.peaktargets, by=c('Gene', 'Peak'), all=T) %>%
        merge(y=conf.peakcores, by=c('Gene', 'Peak'), all=T)
    conf.peak.df[is.na(conf.peak.df)] <- 0

    conf.peak.df$Value <- apply(conf.peak.df, 1, function(x){
        x <- x[-(1:2)]
        x <- setNames(as.logical(as.numeric(x)), names(x))
        if (all(x)) {'Both'} else if (x['Core']){'Core'} else if (x['Target']){'Target'} else if (x['Present']){'Peak'} else 'None'

    })
    
    return(conf.peak.df)

}

do.conf.peak.plot <- function(genes){
    
    conf.peak.df <- get.conf.peak.df(genes)
    conf.peak.df$Gene <- factor(conf.peak.df$Gene, 
                                levels=rev(genes))
    
    conf.peak.plot <- ggplot(conf.peak.df, aes(x=Peak, y=Gene, fill=factor(Value, names(peak.conf.cols)))) +

        geom_tile(color='grey', linewidth=0.1) +
        scale_fill_manual('Disease association', values=peak.conf.cols) +

        theme_pubr(legend = 'bottom') +
        labs_pubr(base_family = 'ArialMT') +

        theme(
            panel.spacing = unit(0.2, "lines"),
            strip.background = element_blank(),
            # strip.placement = 'outside',
            panel.background = element_blank(),
            axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5),
            axis.text.y = element_text(),
            # legend.box="vertical",
            # legend.text = element_text(angle=90, hjust=1, vjust=0.5)
        )
    
    return(conf.peak.plot)
    
}



load('~/codebases/cortical_disorders2/data/graphical.rda', verbose=T)

# save.image('../data/graphical.rda')





### pss
psss <- data.table::fread('PSSs.txt', sep='\t', data.table=F)

psss$groupsample <- paste0(psss$group, ' (', psss$sample, ')')

psss$log1p <- log1p(psss$score)


psss$cell.type <- gsub(pattern = 'not_', replacement = '', fixed=T, psss$group)

groups <- c(sapply(CellTypeOrder, function(x) c(x, paste0('not_', x))), 'Whole data', recursive=T)
groups <- groups[groups %in% unique(psss$group)]

unique(psss$group)[!unique(psss$group) %in% groups]

psss$groupFct <- factor(psss$group, groups)
psss <- psss %>% mutate(sampleFct=factor(sample))

psss <- psss %>% group_by(sample, group) %>% mutate(log1p.scaled=log1p/max(log1p, na.rm=T)) %>% as.data.frame()

psss$Control <- ifelse(grepl('not_', psss$group, fixed=T),'control','test')
                   
rownames(psss) <- 1:nrow(psss)


# psss <- do.call('rbind',
#                 lapply(
#                     c('RGCmaturation', 'NeuralPCW20', 'Gliogenesis'),
#                     function(x){
#                         psss <- data.table::fread(
#                             paste0('../data/',
#                                    x,
#                                    '/perturbations/KO.PerturbationScores.csv'), 
#                             data.table=F) %>% mutate(V1=NULL)
#                         psss$sample <- x
#                         psss$groupsample <- paste0(psss$group, ' (', x, ')')
#                         rownames(psss) <- 1:nrow(psss)
#                         psss
#        }))

# psss$log1p <- log1p(psss$score)

# ## fcts...


# psss$cell.type <- gsub(pattern = 'not_', replacement = '', fixed=T, psss$group)

# groups <- c(sapply(CellTypeOrder, function(x) c(x, paste0('not_', x))), 'Whole data', recursive=T)
# groups <- groups[groups %in% unique(psss$group)]

# unique(psss$group)[!unique(psss$group) %in% groups]

# psss$groupFct <- factor(psss$group, groups)
# psss <- psss %>% mutate(sampleFct=factor(sample))

# psss <- psss %>% group_by(sample, group) %>% mutate(log1p.scaled=log1p/max(log1p, na.rm=T)) 

# psss$Control <- ifelse(grepl('not_', psss$group, fixed=T),'control','test')