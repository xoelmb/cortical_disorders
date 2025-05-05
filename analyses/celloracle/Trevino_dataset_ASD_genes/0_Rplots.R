setwd('/users/genomics/xoel/codebases/co_new/results_foxg1/')
plotdir <- 'Plots/'
dir.create(plotdir, showWarnings = F)

################################################
###################  SCORES  ###################
################################################

# Read Network Scores
NetworkScores <- readxl::read_excel('NetworkScoresFull.xlsx')#, col_types='text')
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
# Perturbations <- readxl::read_excel('PerturbationSummaryFull.xlsx')
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
# PerturbationsCT <- readxl::read_excel('PerturbationCT2CTFull.xlsx')
PerturbationsCT <- readxl::read_excel('PerturbationCT2CT.xlsx')

(unique(PerturbationsCT$CTF) %in% CellTypeOrder)
(unique(PerturbationsCT$CTO) %in% CellTypeOrder)

PerturbationsCT[c('trans.cells', 'init.cells', 'trans.pct', 'Exp')] <- apply(PerturbationsCT[c('trans.cells', 'init.cells', 'trans.pct', 'Exp')], 2, as.numeric)
PerturbationsCT$CTO <- factor(PerturbationsCT$CTO,
                                  levels = CellTypeOrder)
PerturbationsCT$CTF <- factor(PerturbationsCT$CTF,
                                  levels = CellTypeOrder)
PerturbationsCT$SampleLabel <- factor(SampleLabels[PerturbationsCT$Sample], levels=SampleLabelOrder)

# pert.ct.dis.annot <- Diseases[PerturbationsCT$Gene,-1]
# PerturbationsCT <- as.data.frame(cbind(PerturbationsCT, pert.ct.dis.annot))


# # Perturbation Summaries
# PerturbationsCTFull <- readxl::read_excel('PerturbationCT2CTFull.xlsx')

# (unique(PerturbationsCTFull$CTF) %in% CellTypeOrder)
# (unique(PerturbationsCTFull$CTO) %in% CellTypeOrder)

# PerturbationsCTFull[c('trans.cells', 'init.cells', 'trans.pct', 'Exp')] <- apply(PerturbationsCTFull[c('trans.cells', 'init.cells', 'trans.pct', 'Exp')], 2, as.numeric)
# PerturbationsCTFull$CTO <- factor(PerturbationsCTFull$CTO,
#                                   levels = CellTypeOrder)
# PerturbationsCTFull$CTF <- factor(PerturbationsCTFull$CTF,
#                                   levels = CellTypeOrder)
# PerturbationsCTFull$SampleLabel <- factor(SampleLabels[PerturbationsCTFull$Sample], levels=SampleLabelOrder)


# # Perturbation Summaries
# PerturbationsCTNJ <- readxl::read_excel('PerturbationCT2CTNotjoined.xlsx')

# (unique(PerturbationsCTNJ$CTF) %in% CellTypeOrder)
# (unique(PerturbationsCTNJ$CTO) %in% CellTypeOrder)

# PerturbationsCTNJ[c('trans.cells', 'init.cells', 'trans.pct', 'Exp')] <- apply(PerturbationsCTNJ[c('trans.cells', 'init.cells', 'trans.pct', 'Exp')], 2, as.numeric)
# PerturbationsCTNJ$CTO <- factor(PerturbationsCTNJ$CTO,
#                                   levels = CellTypeOrder)
# PerturbationsCTNJ$CTF <- factor(PerturbationsCTNJ$CTF,
#                                   levels = CellTypeOrder)
# PerturbationsCTNJ$SampleLabel <- factor(SampleLabels[PerturbationsCTNJ$Sample], levels=SampleLabelOrder)

# # pert.ct.dis.annot <- Diseases[PerturbationsCT$Gene,-1]
# # PerturbationsCT <- as.data.frame(cbind(PerturbationsCT, pert.ct.dis.annot))



# # Split/duplicate scores by disease
# PerturbationsCTDiseases <- list()
# for (dis in names(gene_list_per_disease)){
#   dis_genes <- gene_list_per_disease[[dis]]
#   dis_pert <- PerturbationsCT[PerturbationsCT$Gene %in% dis_genes,]
#   if (nrow(dis_pert) < 1){
#     next
#   }

#   dis_pert$Disease <- dis
#   PerturbationsCTDiseases[[dis]] <- dis_pert
# }
# PerturbationsCTDiseases <- do.call('rbind', PerturbationsCTDiseases)
# PerturbationsCTDiseases$Disease <- factor(PerturbationsCTDiseases$Disease, DiseaseOrder)

# options(warn=0)

################################################
###############     NETWORKS    ################
################################################


# # Network links
# Links <- data.table::fread('LinkFiltering.txt')
# # (unique(Links$cell.type) %in% CellTypeOrder)

# Links$cell.type <- factor(Links$cell.type,
#                           levels = CellTypeOrder)
# Links$Sample <- factor(Links$Sample, levels=SampleOrder)

# Links2 <- Links %>% group_by(Sample, cell.type) %>% mutate(
#   p_rank = rank(p, na.last=T),
#   p_rank_fr = rank(p, na.last=T)/length(p)
#   )

# Links2[,c('source.in.reg', 'target.in.reg')] <- apply(
#   Links2[,c('source', 'target')], 2, 
#   function(x) {x %in% regulon_genes})

# Links2$any.in.reg <- apply(Links2[,c('source.in.reg', 'target.in.reg')], 1, any)

# Links2[,c('source.is.tf', 'target.is.tf')] <- apply(
#   Links2[,c('source', 'target')], 2, 
#   function(x) {x %in% tf_genes})
# Links2$any.is.tf <- apply(Links2[,c('source.is.tf', 'target.is.tf')],1, any)

# Links2[,c('source.is.dis', 'target.is.dis')] <- apply(
#   Links2[,c('source', 'target')], 2, 
#   function(x) {x %in% disease_genes})
# Links2$any.is.dis <- apply(Links2[,c('source.is.dis', 'target.is.dis')],1, any)


### pss
# psss <- data.table::fread('PSSsFull.txt', sep='\t', data.table=F)
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


                   
# psss.full <- data.table::fread('PSSsFull.txt', sep='\t', data.table=F)
# # psss <- data.table::fread('PSSs.txt', sep='\t', data.table=F)

# psss.full$groupsample <- paste0(psss.full$group, ' (', psss.full$sample, ')')

# psss.full$log1p <- log1p(psss.full$score)


# psss.full$cell.type <- gsub(pattern = 'not_', replacement = '', fixed=T, psss.full$group)

# groups <- c(sapply(CellTypeOrder, function(x) c(x, paste0('not_', x))), 'Whole data', recursive=T)
# groups <- groups[groups %in% unique(psss.full$group)]

# unique(psss.full$group)[!unique(psss.full$group) %in% groups]

# psss.full$groupFct <- factor(psss.full$group, groups)
# psss.full <- psss.full %>% mutate(sampleFct=factor(sample))

# psss.full <- psss.full %>% group_by(sample, group) %>% mutate(log1p.scaled=log1p/max(log1p, na.rm=T)) %>% as.data.frame()

# psss.full$Control <- ifelse(grepl('not_', psss.full$group, fixed=T),'control','test')
                   
# rownames(psss.full) <- 1:nrow(psss.full)


                   
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