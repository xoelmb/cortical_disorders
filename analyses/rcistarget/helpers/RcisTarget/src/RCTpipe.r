Sys.setenv(LOG_LEVEL='TRACE')
Sys.getenv('R_HOME')

PAR_METHOD <- 'iCisTarget'
PAR_SEARCH_RANGE <- 'long'
# PAR_DIRECT_ANNOT <- TRUE
PAR_MOTIF_NES <- 3
PAR_MOTIF_DB_V <- 'V10'
PAR_DB_DIR <- '/users/genomics/xoel/canonades/RcisTarget/cistarget_dbs/'
PAR_RES_DIR <- './cistarget_results/'
PAR_CORES <- 80

PAR_MOTIF_DB_INFO <- list(
    V9 = list(
        MOTIF_TF_URL = 'https://resources.aertslab.org/cistarget/motif2tf/motifs-v9-nr.hgnc-m0.001-o0.0.tbl',
        INDEX_COL = 'features',
        MOTIF_SEARCH_URL = list(
            short = "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__500bp_up_and_100bp_down_tss.mc9nr.feather",
            long = "https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc9nr/gene_based/hg38__refseq-r80__10kb_up_and_down_tss.mc9nr.feather"
    )),
    V10 = list(
        MOTIF_TF_URL = 'https://resources.aertslab.org/cistarget/motif2tf/motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl',
        INDEX_COL = 'motifs',
        MOTIF_SEARCH_URL = list(
            short = 'https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_500bp_up_100bp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather',
            long = 'https://resources.aertslab.org/cistarget/databases/homo_sapiens/hg38/refseq_r80/mc_v10_clust/gene_based/hg38_10kbp_up_10kbp_down_full_tx_v10_clust.genes_vs_motifs.rankings.feather'
    )))


download_dbs <- function(
    motif_db = PAR_MOTIF_DB_V,
    motif_search_range = PAR_SEARCH_RANGE,
    motif_db_dir = PAR_DB_DIR
){
    
    dir.create(motif_db_dir, showWarnings = F)
    
    pars <- PAR_MOTIF_DB_INFO[[motif_db]]    
    
    
    
    # Motif - TF database
    motif_tf_destfile <- file.path(motif_db_dir, basename(pars$MOTIF_TF_URL))
    
    if (!file.exists(motif_tf_destfile)){
        rlog::log_info('Downloading motif-TF DB')
        download.file(
            url = pars$MOTIF_TF_URL,            
            destfile = motif_tf_destfile,
            verbose=T, quiet=F, method='wget'
        )} else {
        rlog::log_info('Motif-TF DB already downloaded')
    }
    
    
    # Motif database
    motif_db_url <- pars$MOTIF_SEARCH_URL[[motif_search_range]]
    motif_db_destfile <- file.path(motif_db_dir, basename(motif_db_url))
    
    if (!file.exists(motif_db_destfile)){
        rlog::log_info('Downloading motif ranking DB')
        download.file(
            url = motif_db_url,            
            destfile = motif_db_destfile,
            verbose=T, quiet=F, method='wget'
        )} else {
        rlog::log_info('Motif ranking DB already downloaded')
    }
    
    
    rlog::log_info('Loading motif annotations...')

    motifRankings <- RcisTarget::importRankings(
      motif_db_destfile, 
      indexCol = pars$INDEX_COL)
    
    motifAnnotations <- RcisTarget::importAnnotations(
        motif_tf_destfile,
        motifsInRanking=RcisTarget::getRanking(motifRankings)[,pars$INDEX_COL,drop=T],
        columnsToKeep = NULL)

    return(list(
        motifRankings = motifRankings,
        motifAnnotations = motifAnnotations
    ))
}

clean_genes_df <- function(df, motif_dbs){
    
    # Obtain list of genes in Ranking DB
    genes.in.ranking <- colnames(motif_dbs$motifRankings@rankings)[-1]
    
    # Clean df
    clean.df <- df[rownames(df) %in% genes.in.ranking,]
    
    # Get number of original genes
    n.genes.original <- colSums(df)
    
    # Get number of genes after cleaning
    n.genes.clean <- colSums(clean.df)
    
    rlog::log_trace('Number of genes NOT present in Motif Ranking DB')
    rlog::log_trace(
        paste(names(n.genes.original), 
              n.genes.original-n.genes.clean,
              '(', round((n.genes.original-n.genes.clean)/n.genes.original*100, 2), '% )'
             ))
    
    
    if (any(n.genes.clean == 0)) {
        
        rlog::log_info('Removing the following gene lists after cleaning')
        rlog::log_info(names(n.genes.clean)[n.genes.clean==0])
        
        clean.df <- clean.df[,n.genes.clean!=0]
    }
    
    
    return(clean.df)
    
}

df_to_genelist <- function(df){
    return(lapply(df, FUN = function(x){rownames(df)[as.logical(x)]}))
}

compute_RCT <- function(
    geneLists,
    motif_dbs,
    ncores=PAR_CORES,
    nes_threshold=PAR_MOTIF_NES,
    method=PAR_METHOD,
    verbose=T
){
  
    
    t0 <- Sys.time()

    rlog::log_info('Running RcisTarget')
    
    ## 1. Compute AUC of motif~geneset
    rlog::log_debug('1/4')
    motifs_AUC <- RcisTarget::calcAUC(
        geneLists, 
        motif_dbs$motifRankings,
        nCores=ncores, 
        verbose=verbose)
    
    ## 2. Select significant motifs and annotate to TFs
    rlog::log_debug('2/4')
    motifEnrichmentTable <- RcisTarget::addMotifAnnotation(
        motifs_AUC,
        nesThreshold=nes_threshold,
        motifAnnot=motif_dbs$motifAnnotations)
    
    # 3. Identify most-enriched genes for each motif in each geneset
    rlog::log_debug('3/4')
    motifEnrichmentTable <- RcisTarget::addSignificantGenes(
        motifEnrichmentTable,
        rankings=motif_dbs$motifRankings,
        geneSets=geneLists,
        method=method,
        nCores=ncores)
      
    # 4. Add logo to results
    rlog::log_debug('4/4')
    motifEnrichmentTable <- RcisTarget::addLogo(motifEnrichmentTable)
    
      
    # 5. Fix "2 genes in one entry" case
    motifEnrichmentTable$TF_highConf <- gsub(').', ');', motifEnrichmentTable$TF_highConf, fixed=T)
    
    rlog::log_debug(capture.output(Sys.time() - t0))
    
    return(motifEnrichmentTable)
    
}

motifEnrTbl_to_long <- function(motifEnrichmentTable, geneLists){
    
    library(dplyr)
    
    # From many-TF-per-row to 1-TF-per-row
    longMotifEnrTbl <- tidyr::separate_rows(motifEnrichmentTable, TF_highConf, sep = ';' )
    
    # Annotate directAnnotation
    longMotifEnrTbl$TF_directAnnotation <- grepl('directAnnotation', longMotifEnrTbl$TF_highConf, fixed=T)
    
    # Delete empty rows
    longMotifEnrTbl <- longMotifEnrTbl[sapply(longMotifEnrTbl$TF_highConf, nchar) != 1,]
    
    # Remote annotation from names
    longMotifEnrTbl$TF_highConf <- gsub("\\s*\\([^\\)]+\\)", '', longMotifEnrTbl$TF_highConf)
    longMotifEnrTbl$TF_highConf <- gsub(" ", '', longMotifEnrTbl$TF_highConf, fixed=T)
    
    # Add presence in original list: regulon or not
    longMotifEnrTbl$TF_inGeneList <- apply(longMotifEnrTbl, 1, function(x){
        x['TF_highConf'] %in% geneLists[[x['geneSet']]]        
    })
    
    # Separate targets
    longMotifEnrTbl <- tidyr::separate_rows(longMotifEnrTbl, enrichedGenes, sep = ';')
    
    # Filter those TFs that are not present in the original list
    longMotifEnrTblSelf <- subset(longMotifEnrTbl, TF_inGeneList)

    
    return(list(
        motEnrTblLong = longMotifEnrTbl,
        motEnrTblLongSelf = longMotifEnrTblSelf   
    ))
}

longSelf_to_summary <- function(motifEnrichmentTableLongSelf){
    
    # Condense (targets, motifs) for summary of cores 
    summaryTbl <- motifEnrichmentTableLongSelf %>% 
        group_by(geneSet, TF_highConf) %>% 
        summarise(
            n.motifs = length(unique(motif)),
            all.motifs = paste(unique(motif), collapse = '; '),
            n.targets = length(unique(enrichedGenes)),
            all.targets = paste(unique(enrichedGenes), collapse = '; ')
        )
    
    return(summaryTbl)
}

longSelf_to_regulon <- function(motifEnrichmentTableLongSelf){
    
    # Condense motifs for regulonTbl
    regulonTbl <- motifEnrichmentTableLongSelf %>% 
        group_by(geneSet, TF_highConf, enrichedGenes) %>% 
        summarise(
            # Target=enrichedGenes,
            n.motifs = length(unique(motif)),
            all.motifs = paste(unique(motif), collapse = '; ')) %>%
        rename(Core = TF_highConf, Target = enrichedGenes)
    
    return(regulonTbl)
}

run_pipe <- function(
    
    genes.df,
    
    motif_search_range = PAR_SEARCH_RANGE,
    # direct_annotation_tf_only = PAR_DIRECT_ANNOT,
    motif_NES_threshold = PAR_MOTIF_NES,
    motif_db = PAR_MOTIF_DB_V,
    motif_db_dir = PAR_DB_DIR,
    motif_method = PAR_METHOD,
    
    ncores = PAR_CORES
    
)
    
{
    
    run.pars <- list(
        motif_search_range = motif_search_range,
        motif_NES_threshold = motif_NES_threshold,
        motif_db = motif_db,
        motif_db_dir = motif_db_dir,
        motif_method = motif_method,
        ncores=ncores
    )
    
    # Download RcisTarget Databases
    motif_dbs <- download_dbs(
        motif_db = motif_db,
        motif_search_range=motif_search_range,
        motif_db_dir=motif_db_dir
    )
    
    # Clean data.frame: remove genes not present in db
    clean.df <- clean_genes_df(genes.df, motif_dbs)
    
        
    # Get gene lists
    geneLists <- df_to_genelist(clean.df)
    
    # Keep original gene.order
    gene.order.original <- rownames(clean.df)
    
    
    # Run RcisTarget on geneLists
    motifEnrichmentTable <- compute_RCT(
        geneLists,
        motif_dbs=motif_dbs,
        ncores=ncores,
        nes_threshold=motif_NES_threshold,
        method=motif_method,
        verbose=T)
    
    
    # Make long data.frame
    longTbls <- motifEnrTbl_to_long(
        motifEnrichmentTable,
        geneLists = geneLists
    )
    motEnrTblLong <- longTbls$motEnrTblLong
    motEnrTblLongSelf <- longTbls$motEnrTblLongSelf
    
    
    # Make regulon table
    summaryTbl <- longSelf_to_summary(motEnrTblLongSelf)
    regulonTbl <- longSelf_to_regulon(motEnrTblLongSelf)
    
    # Get list of targets per (geneSet core)
    targetLists <- split(
        regulonTbl$Target, 
        paste(regulonTbl$geneSet, regulonTbl$Core))
    
    RCT_PKG <- list(
        input.df = genes.df,
        clean.df = clean.df,
        motif_dbs = motif_dbs,
        geneLists = geneLists,
        motifEnrichmentTable = motifEnrichmentTable,
        motifEnrichmentTableLong = motEnrTblLong,
        motifEnrichmentTableLongSelf = motEnrTblLongSelf,
        summaryTbl = summaryTbl,
        regulonTbl = regulonTbl,
        targetLists = targetLists,
        pars = run.pars,
        plots = list()
    )
    

    
    
    return(RCT_PKG)
    
}

save_RCT <- function(RCT_results, save_dir = PAR_RES_DIR){
    
    rlog::log_info(paste('Saving RCT to', save_dir))
    dir.create(save_dir, showWarnings = F, recursive = T)
    
    if ('permutations' %in% names(RCT_results)){

        rlog::log_info('Saving permutated results')

        # Save the regulon table
        write.csv(RCT_results$regulonTblPvals,
                  file = file.path(save_dir, 'cisTarget_regulons.permuted.csv'),
                  row.names=F)
        
        # First save the result list
        saveRDS(RCT_results, 
                file = file.path(save_dir, 'cisTarget_results.permuted.rds'))
    } else {
        
        rlog::log_info('Saving results without permutations')

        # Save the regulon table
        write.csv(RCT_results$regulonTbl,
                  file = file.path(save_dir, 'cisTarget_regulons.csv'),
                  row.names=F)

        # First save the result list
        saveRDS(RCT_results, 
                file = file.path(save_dir, 'cisTarget_results.rds'))
    }
    
    # Save the summary table
    write.csv(RCT_results$summaryTbl,
              file = file.path(save_dir, 'cisTarget_regulons.summary.csv'),
              row.names=F)
    
    # Save bicorrelations
    if ('bicorRes' %in% names(RCT_results)){
        
        rlog::log_info('Saving results of bicorrelations')

        # Save the bicor summary (per core)
        write.csv(RCT_results$bicorRes$bicorSummary,
                  file = file.path(save_dir, 'cisTarget_regulons.bicorsSummary.csv'),
                  row.names=F)
        # Save the bicor full table (core-target)
        write.csv(RCT_results$bicorRes$bicorRes,
                  file = file.path(save_dir, 'cisTarget_regulons.bicorResults.csv'),
                  row.names=F)
    }

    
}

# Function to create e shuffled dataframe with same order of rows
shuffle_df_with_suffix <- function(df, suffix){
    
    set.seed(suffix)
    
    shuf.df <- data.frame(
        apply(df, 2, sample, replace=F),
        row.names=rownames(df))
    
    # shuf.df <- data.frame(as.matrix(df), 
    #                       row.names=sample(rownames(df), replace = F))[rownames(df),]
    
    colnames(shuf.df) <- paste(colnames(df), as.character(suffix), sep='.rep' )
    
    return(shuf.df)
    
}

# Function to permute in parallel a data.frame n times and create a wide one
permute_df <- function(df, ncores=PAR_CORES, n.copies = 100) {
    
    library(doParallel)  

    cl <- makeCluster(min(20,ncores))
    registerDoParallel(cl)  
    on.exit(stopCluster(cl))

    

    perm.df <- do.call('cbind', 
                       foreach(i=1:n.copies,
                               .export = c('shuffle_df_with_suffix')) %dopar% shuffle_df_with_suffix(df, i))

    # stopCluster(cl = cl)
    
    return(perm.df)
}

compute_permutations <- function(
    df, motif_dbs, 
    n.permutations=10,
    ncores=PAR_CORES,
    nes_threshold=PAR_MOTIF_NES,
    method=PAR_METHOD,
    verbose=T
){
    
    # Create permuted lists of genes
    perm.df <- permute_df(df, ncores=ncores, n.copies=n.permutations)
    perm.geneLists <- df_to_genelist(perm.df)
    
    # message('These are the permuted lists:')
    # for (l in names(perm.geneLists)){
    #     print(l)
    #     message(paste(head(perm.geneLists[[l]]), collapse=', '))
    # }
    
    # Run cisTarget
    perm.res <- compute_RCT(
        perm.geneLists,
        motif_dbs=motif_dbs,
        ncores=ncores,
        nes_threshold=nes_threshold,
        method=method,
        verbose=T
    )
    
    # Make long data.frameperm.cores
    perm.longres <- motifEnrTbl_to_long(perm.res, geneLists = perm.geneLists)$motEnrTblLong
    perm.longres[,c('geneSet', 'rep')] <- matrix(unlist(strsplit(perm.longres$geneSet,
                                                                 split = '.rep', fixed = T)),
                                                 ncol = 2, byrow = T)
    
    
    # Get number of observations
    perm.cores <- perm.longres %>% 
        group_by(geneSet, TF_highConf) %>%
        summarise(counts = length(unique(rep)))
    
    # Get nominal p-value
    perm.cores$pval.nom <- perm.cores$counts / n.permutations 
    
    # Get p-value FDR adjusted in all permutations
    perm.cores$pval.adj <- p.adjust(perm.cores$pval.nom, method = 'fdr')
    
    # Adjust p.value per geneSet
    perm.cores <- perm.cores %>% 
        group_by(geneSet) %>% 
        summarize(
            TF_highConf, 
            TFperm_counts = counts, 
            TFperm_pval.nom = pval.nom,
            TFperm_pval.adj = pval.adj,
            TFperm_pval.adj.geneSet = p.adjust(pval.nom, method='fdr')
        )
    
    
    return(list(
        perm.df = perm.df,
        perm.longres = perm.longres,
        perm.cores = perm.cores
    ))
    
}

permute_RCT <- function(RCT_result, 
                        n.permutations = 10,
                        only_regulons = TRUE
                       ){

    if (!only_regulons){
        df.to_permute <- RCT_result$clean.df
    } else {
        df.to_permute <- RCT_result$clean.df[,colnames(RCT_result$clean.df) %in% unique(RCT_result$regulonTbl$geneSet)]
    }
    
    message(paste0('Number of genesets to permute: ', ncol(df.to_permute)))
    
    # Compute permutations
    perm.reslist <- compute_permutations(
        df = df.to_permute,
        motif_dbs = RCT_result$motif_dbs,
        n.permutations = n.permutations,
        ncores = RCT_result$pars$ncores,
        nes_threshold = RCT_result$pars$motif_NES_threshold,
        method = RCT_result$pars$motif_method,
        verbose = T)    
    
    # Annotate p-values
    regulonTblPvals <- merge(
        RCT_result$regulonTbl,
        perm.reslist$perm.cores,
        by.x=c('geneSet', 'Core'),
        by.y=c('geneSet', 'TF_highConf'),
        all.x=T, all.y=F)
    
    # Assign 0 to not observed p-values
    regulonTblPvals[is.na(regulonTblPvals)] <- 0

    # Append results
    RCT_result[['regulonTblPvals']] <- regulonTblPvals
    RCT_result[['permutations']] <- perm.reslist
    
    return(RCT_result)

}

get_probability_matrix <- function(RCT_result, exp.data){

    # Get sample of peak expression per gene
    max_sample <- apply(exp.data, 1, function(x){names(x)[which.max(x)]})
    max_sample.ngenes <- table(max_sample)[colnames(exp.data)]

    # Filter geneLists for those not present in expression data
    geneLists.exp <- lapply(RCT_result$geneLists, function(x){x[x %in% names(max_sample)]})

    # Get number of genes per sample in each disease
    genes.per.sample <- do.call(
        'rbind', 
        lapply(
            geneLists.exp, function(x){
                setNames(table(max_sample[x])[colnames(exp.data)],
                         colnames(exp.data))
    }))
    genes.per.sample[is.na(genes.per.sample)] <- 0

    # Get fractions per sample
    fracs.per.sample <- genes.per.sample / rowSums(genes.per.sample)

    # Get probabilities per sample, divided by n.gene
    probs.per.sample <- data.frame(t(apply(
        fracs.per.sample, 1,
        function(x){ x / max_sample.ngenes }
    )))

    # Get probabilities per gene
    probs.per.gene <- apply(
        probs.per.sample, 1, function(x){
            setNames(x[max_sample], names(max_sample))
        })
    
    return(probs.per.gene)
}

get_bicor <- function(core, targets, exp.data, geneset=NULL){
    
    require(dplyr)
    
    # Check missing genes
    rlog::log_trace(paste('Core in expression data:', core %in%rownames(exp.data)))
    rlog::log_trace(paste('Targets not in expression data:', 
                          sum(!targets %in%rownames(exp.data))))
    
    targets <- targets[targets %in% rownames(exp.data)]
    
    # Retrieve expression data of core and targets
    core.exp <- t(as.matrix(exp.data[core,]))
    targets.exp <- t(exp.data[targets,])

    # Compute bicor
    bicor.list <- WGCNA::bicorAndPvalue(
        x = core.exp,
        y = targets.exp,
        alternative = 'two.sided'
    )

    # Parse and annotate results
    bicor.df <- as.data.frame(t(do.call('rbind', bicor.list)))
    colnames(bicor.df) <- c('bicor', paste('bicor', names(bicor.list)[-1], sep='.'))

    bicor.df$Target <- rownames(bicor.df)
    rownames(bicor.df) <- NULL

    bicor.df$Core <- core
    
    # one-line summary
    summary.df <- bicor.df %>% 
        summarise(
            Core = unique(Core),
            mean.abs = mean(abs(bicor), na.rm=T),
            median.abs = median(abs(bicor), na.rm=T),
            n.sign = sum(bicor.p < 0.05, na.rm=T),
            fr.sign = sum(bicor.p < 0.05, na.rm=T)/length(bicor),
            
            mean.pos = mean(bicor[bicor>0], na.rm=T),
            median.pos = median(bicor[bicor>0], na.rm=T),
            n.sign.pos = sum(bicor.p[bicor>0] < 0.05, na.rm=T),
            fr.sign.pos = sum(bicor.p[bicor>0] < 0.05, na.rm=T)/length(bicor),
            
            mean.neg = mean(bicor[bicor<0], na.rm=T),
            median.neg = median(bicor[bicor<0], na.rm=T),
            n.sign.neg = sum(bicor.p[bicor<0] < 0.05, na.rm=T),
            fr.sign.neg = sum(bicor.p[bicor<0] < 0.05, na.rm=T)/length(bicor)
        )
    
    if (!is.null(geneset)){
        bicor.df$geneSet <- geneset
        summary.df$geneSet <- geneset
    }
    
    return(list(
        bicors = bicor.df,
        summary = summary.df
    ))
}

permute_bicor <- function(exp.data, core, geneset, n.targets,
                          n.permutations = 100, exp.probs=NULL, ncores=PAR_CORES){

    # Init parallel workers
    require(doParallel)  
    
    cl <- makeCluster(PAR_CORES)
    registerDoParallel(cl)  
    on.exit(stopCluster(cl))



    # Run permutations
    rndBicorRes <- foreach(i=1:n.permutations,
                           .export= c('get_bicor')) %dopar% {

        # Get random genes
        random.targets <- sample(
            x = rownames(exp.data), replace = F, 
            size=n.targets, prob = exp.probs)    

        # Compute bicorrelation
        bicorRes <- get_bicor(
            core=core, 
            targets=random.targets,
            exp.data=exp.data,
            geneset=paste(geneset, i, sep='.rep')) 


    }
    # # Stop parallel workers
    # stopCluster(cl = cl)

    # Condense permuted bicorrelations
    longRndBicor <- do.call('rbind', lapply(rndBicorRes, function(x){x[[1]]}))
    # Same for summary
    sumRndBicor <- do.call('rbind', lapply(rndBicorRes, function(x){x[[2]]}))

    # Return both
    return(list(
        fullBicors = longRndBicor,
        summBicors = sumRndBicor
    ))

}


annotate_bicor_summary <- function(bicorRes, bicorResPerm){
    
    annot.summary <- bicorRes$summary %>% 
        summarise(
            Core, 
            geneSet, 
            mean.abs,
            mean.abs.pval = mean(bicorResPerm$summBicors$mean.abs > bicorRes$summary$mean.abs, na.rm=T),
            median.abs,
            median.abs.pval = mean(bicorResPerm$summBicors$median.abs > bicorRes$summary$median.abs, na.rm=T),
            n.sign,
            n.sign.pval = mean(bicorResPerm$summBicors$n.sign > bicorRes$summary$n.sign, na.rm=T),
            fr.sign,
            fr.sign.pval = mean(bicorResPerm$summBicors$fr.sign > bicorRes$summary$fr.sign, na.rm=T),
            
            mean.pos,
            mean.pos.pval = mean(bicorResPerm$summBicors$mean.pos > bicorRes$summary$mean.pos, na.rm=T),
            median.pos,
            median.pos.pval = mean(bicorResPerm$summBicors$median.pos > bicorRes$summary$median.pos, na.rm=T),
            n.sign.pos,
            n.sign.pos.pval = mean(bicorResPerm$summBicors$n.sign.pos > bicorRes$summary$n.sign.pos, na.rm=T),
            fr.sign.pos,
            fr.sign.pos.pval = mean(bicorResPerm$summBicors$fr.sign.pos > bicorRes$summary$fr.sign.pos, na.rm=T),
            
            mean.neg,
            mean.neg.pval = mean(bicorResPerm$summBicors$mean.neg < bicorRes$summary$mean.neg, na.rm=T),
            median.neg,
            median.neg.pval = mean(bicorResPerm$summBicors$median.neg < bicorRes$summary$median.neg, na.rm=T),
            n.sign.neg,
            n.sign.neg.pval = mean(bicorResPerm$summBicors$n.sign.neg > bicorRes$summary$n.sign.neg, na.rm=T),
            fr.sign.neg,
            fr.sign.neg.pval = mean(bicorResPerm$summBicors$fr.sign.neg > bicorRes$summary$fr.sign.neg, na.rm=T)

        )

    return(annot.summary)
}

# Create function to compute and permute bicorrelations from RCT
compute_bicors <- function(
    RCT_result, exp.data, 
    n.permutations=10, ncores=PAR_CORES){


    # Get probability matrix. 
    # we need it to sample genes from expression data following the peak
    # sample distribution of the original geneset (not the regulon one)
    sampling.prob.matrix <- get_probability_matrix(RCT_result = RCT_result, exp.data = exp.data)


    # Init buffers for results
    bicor.summs <- list()
    bicor.dfs <- list()
    perm.summs <- list()
    perm.dfs <- list()

    
    filt.regTbl <- subset(RCT_result$regulonTbl, 
                          Target %in% rownames(sampling.prob.matrix))

    for (geneset in unique(RCT_result$regulonTbl$geneSet)){

        rlog::log_debug(paste('[BICOR] geneSet', geneset))

        # Define vector of probabilities
        geneset.geneprobs <- sampling.prob.matrix[, geneset, drop=T]

        for (core in unique(subset(filt.regTbl, 
                                   geneSet == geneset)$Core)){

            rlog::log_debug(paste('[BICOR]    Core', core))


            # Get list of targets        
            targets <- subset(filt.regTbl,
                              geneSet == geneset & Core == core)$Target
            n.targets <- length(targets)


            # Compute bicorrelations
            bicorRes <- get_bicor(
                core = core,
                targets = targets, 
                exp.data = exp.data,
                geneset = geneset) 

            
            # Compute perm bicors
            bicorResPerm <- permute_bicor(
                exp.data = exp.data,
                core=core,
                geneset=geneset,
                n.targets=n.targets,
                n.permutations=n.permutations, 
                exp.probs = geneset.geneprobs, 
                ncores=ncores)


            # Compute pvalues
            annot.bicor.summ <- annotate_bicor_summary(bicorRes=bicorRes, bicorResPerm=bicorResPerm)


            # Append
            bicor.summs[[paste(geneset, core)]] <- annot.bicor.summ
            bicor.dfs[[paste(geneset, core)]] <- bicorRes$bicors

            # Append perms
            perm.summs[[paste(geneset, core)]] <- bicorResPerm$summBicors
            perm.dfs[[paste(geneset, core)]] <- bicorResPerm$fullBicors


        }

    }

    # Condense
    bicor.summs.dense <- do.call('rbind', bicor.summs)
    bicor.dfs.dense <- do.call('rbind', bicor.dfs)

    perm.summs.dense <- do.call('rbind', perm.summs)
    perm.dfs.dense <- do.call('rbind', perm.dfs)

    # Return
    return(list(
        bicorSummary = bicor.summs.dense,
        bicorRes = bicor.dfs.dense,
        permSummary = perm.summs.dense,
        permRes = perm.dfs.dense,
        n.permutations = n.permutations
    ))

}

fisher.reg.test <- function(geneset.genes, reg.genes, test.genes, haldane_correction=T){
    
    require(dplyr)
    
    # Define levels
    lvls <- c('REG TEST', 'REG NOT TEST', 'NOT REG TEST', 'NOT REG NOT TEST')
    
    # Get genes in regulon and test set
    in.reg <- ifelse(geneset.genes %in% reg.genes, 'REG', 'NOT REG')
    in.test <- ifelse(geneset.genes %in% test.genes, 'TEST', 'NOT TEST')
    
    # Create table
    ot <- table(paste(in.reg, in.test, sep = ' '))
    # Fill 0
    # if zero, correct haldane
    correct_haldane <- F
    for (cc in lvls[!lvls %in% names(ot)]){
        ot[cc] <- 0
        correct_haldane <- T
    }
    # Reorder
    ot <- ot[lvls]
    
    # Test
    fish <- unlist(fisher.test(matrix(ot, ncol=2, nrow=2, byrow = T), alternative = 'greater'))
    # Rename
    names(fish) <- names(fish) %>% gsub(pattern = 'conf.int2', replacement = 'conf.int.max') %>% gsub(pattern = 'conf.int1', replacement = 'conf.int.min')
    fish <- fish[1:(length(fish)-2)]   
    
    # Join
    fish <- c(ot, fish)
    
    # Correct Haldane
    if (haldane_correction & correct_haldane){
        
        message('Using Haldane correction')
        ot <- ot + 0.5
        
        fish['estimate.odds ratio'] <- unname((ot[1] * ot[4]) / (ot[2] * ot[3]))
        fish['alternative'] <- 'haldane'
        
    } 
    
    # data.frame
    fish <- as.data.frame(t(as.data.frame(fish)))

    
    return(fish)
}

odds.ratio.RCT <- function(RCT_result, df_to_test, haldane_correction=T){
    
    regulonTbl <- RCT_result$regulonTbl

    # Remove uncomputed genes
    untested.ngenes <- sum(!rownames(df_to_test) %in% rownames(RCT_result$clean.df)) 
    # Info
    rlog::log_info(paste0('Number of genes not tested in RcisTarget: ', untested.ngenes))
    df_to_test <- df_to_test[rownames(df_to_test) %in% rownames(RCT_result$clean.df),]


    # Result buffer
    fishes <- list()

    for (geneset in unique(regulonTbl$geneSet)){

        genesetTbl <- subset(regulonTbl, geneSet == geneset)
        genesetGenes <- RCT_result$geneLists[[geneset]]

        rlog::log_info(paste0(geneset, ': ', length(genesetGenes)))

        for (core in unique(genesetTbl$Core)){

            coreTbl <- subset(genesetTbl, Core == core)
            regGenes <- unique(c(coreTbl$Core, coreTbl$Target))

            rlog::log_info(paste0('Regulon core ', core, ': ', length(regGenes)))

            for (test.gs in colnames(df_to_test)){

                testGenes <- rownames(df_to_test)[as.logical(df_to_test[,test.gs])]

                rlog::log_info(paste0(
                    '\tTesting ', test.gs, ': ', length(testGenes)))

                fish <- fisher.reg.test(geneset.genes = genesetGenes, reg.genes = regGenes, test.genes = testGenes, haldane_correction=haldane_correction)

                fish$geneSet <- geneset
                fish$Core <- core
                fish$testSet <- test.gs

                fishes[[length(fishes)+1]] <- fish
            }
        }
    }
    

    fishes <- do.call('rbind', fishes)
    rownames(fishes) <- 1:nrow(fishes)
    fishes[,1:9] <- apply(fishes[,1:9], 2, as.numeric)
    
    # Adjust p-values
    # Per Regulon
    fishes <- fishes %>% group_by(geneSet, Core) %>% mutate(
        p.value.adj.core = p.adjust(p.value, method='fdr')
    )
    # Per geneSet
    fishes <- fishes %>% group_by(geneSet) %>% mutate(
        p.value.adj.geneset = p.adjust(p.value, method='fdr')
    )
    # Per testSet
    fishes <- fishes %>% group_by(testSet) %>% mutate(
        p.value.adj.testset = p.adjust(p.value, method='fdr')
    )
    # All
    fishes$p.value.adj <- p.adjust(fishes$p.value, method='fdr')
    # Log2 odd-ratio
    fishes$OR.log2 <- log2(fishes$`estimate.odds ratio`)    

    return(fishes)
    
}

PAR_EXPORT <- FALSE

if (PAR_EXPORT) {
    message('Exporting code...')
    system(
        command = 'jupyter nbconvert --to script --output RCTpipe --output-dir ./src/ ./RCTpipe.ipynb',
        intern=TRUE)
}


