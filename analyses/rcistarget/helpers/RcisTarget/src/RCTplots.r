figsize <- function(x,y){
    options(repr.plot.width=x, repr.plot.height=y)
}

upset_rct <- function(
    RCT_result, 
    original_list=TRUE, 
    target_list=TRUE, 
    geneSet=NULL){
    
    
    if (original_list){
        
        p <- UpSetR::upset(UpSetR::fromList(RCT_result$geneLists),
                           order.by = "freq", nintersects = 25, nsets = 30)
        RCT_result[['plots']][['Input lists']] <- p
        plot(p)
    }
    
    if (target_list){
        
        p <- UpSetR::upset(UpSetR::fromList(RCT_result$targetLists),
                           order.by = "freq", nintersects = 25, nsets = 30)
        RCT_result[['plots']][['Target lists']] <- p
        plot(p)
    }
    
    if (!is.null(geneSet)){
        
        gsLists <- RCT_result$targetLists[grep(geneSet, names(RCT_result$targetLists), fixed=T, value=T)]
        
        p <- UpSetR::upset(UpSetR::fromList(gsLists),
                           order.by = "freq", nintersects = 25, nsets = 30)
        RCT_result[['plots']][[geneSet]] <- p
        plot(p)
    }
    
    return(RCT_result)
}











PAR_EXPORT <- FALSE

if (PAR_EXPORT) {
    message('Exporting code...')
    system(
        command = 'jupyter nbconvert --verbose --to script --output RCTplots --output-dir ./src/ ./RCTplots.ipynb',
        intern=TRUE)
}
