source('~/canonades/bioinforgalician/src/R/utils.R', echo=T)

library(cicero)
library(monocle3)

# Load data and metadata
load_data <- function(atac_path, peaks_bed_path, cell_metadata_path, atac_sep='\t', verbose=F){
    
    if (verbose){message('Loading metadata')}
    ## Peak metadata
    peakinfo <- data.table::fread(peaks_bed_path, select = c(1:3), data.table=F, verbose=F)
    # Format peak info
    names(peakinfo) <- c("chr", "bp1", "bp2")
    peakinfo$site_name <- paste(peakinfo$chr, peakinfo$bp1, peakinfo$bp2, sep="_")
    row.names(peakinfo) <- peakinfo$site_name

    ## Cell metadata
    cellinfo <- data.table::fread(cell_metadata_path, header=T, data.table=F, verbose=F)
    # Format cell info
    rownames(cellinfo) <- cellinfo[,1]
    
    
    if (verbose){message('Loading  atac data')}
    indata <- load_sparse_csv(atac_path, 
                              sep=atac_sep, 
                              header=T,
                              # nrows=300,
                              # select=c(1:3),
                              verbose=verbose)

    peakinfo <- peakinfo[rownames(indata),]
    cellinfo <- cellinfo[colnames(indata),]

    return(list(
        atac_data=indata,
        peakinfo=peakinfo,
        cellinfo=cellinfo
    ))
}
# Preprocess scATAC-Seq Data

## Create CDS object 

create_cds <- function(indata, cellinfo, peakinfo, binarize=TRUE){
    
    cds <- indata

    if (binarize){
        # Binarize the matrix
        cds@x[cds@x > 0] <- 1
    }
    
    cds <- suppressWarnings(
        monocle3::new_cell_data_set(
            indata,
            cell_metadata = cellinfo,
            gene_metadata = peakinfo))
    
    return(cds)
}

## Quality check and Filtering
qc_cds <- function(input_cds, min_count=2000, max_count=Inf){
    
    # Ensure there are no peaks included with zero reads
    cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,] 

    
    # Visualize peak_count_per_cell
    hist(Matrix::colSums(exprs(cds)), breaks = 30)

    # Filter cells by peak_count
    # Please set an appropriate threshold values according to your data 
    cds <- cds[,Matrix::colSums(exprs(cds)) >= min_count] 
    cds <- cds[,Matrix::colSums(exprs(cds)) <= max_count] 

    # Visualize peak_count_per_cell
    hist(Matrix::colSums(exprs(cds)), breaks = 30)
    
    # Ensure there are no peaks included with zero reads after removing cells
    cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,] 

    return(cds)
}

## Prepare CDS
prepare_cds <- function(input_cds, verbose=TRUE){
    
    if (verbose) {message('[1] Detecting genes')}
    # Data preprocessing
    cds <- monocle3::detect_genes(input_cds)
    
    if (verbose) {message('[2] Estimating size factors')}
    cds <- monocle3::estimate_size_factors(cds)
    
    if (verbose) {message('[3] Preprocessing')}
    cds <- monocle3::preprocess_cds(cds, method = "LSI")
    
    if (verbose) {message('[4] UMAP')}
    cds <- monocle3::reduce_dimension(cds,
                                      reduction_method = 'UMAP', 
                                      preprocess_method = "LSI")
    return(cds)
}




# Cicero

analysis_cicero <- function(input_cds, chromosome_length,
                            return_cicero_cds = FALSE, 
                            verbose=TRUE){

    if (verbose) {message('[1] Making cicero CDS')}
    # Make cicero')
    cicero_cds <- cicero::make_cicero_cds(input_cds,
                                  reduced_coordinates = reducedDims(input_cds)$UMAP)

    if (verbose) {message('[2] Running cicero analysis')}
    # Run the main function 
    # Takes long to run
    conns <- cicero::run_cicero(
        cicero_cds, 
        chromosome_length)
    
    
    if (return_cicero_cds){
        return(list('cicero_cds'=cicero_cds,
                    'conns'=conns))
    } else{
        return(conns)
    }

}



# PIPELINE

pipe_cicero_conns <- function(
    
    atac_path, peaks_bed_path, cell_metadata_path, 
    plot_by=NULL, save_dir=NULL, prefix=NULL, verbose=TRUE){
    
    v <- verbose
    
    
    if (v) {message('[1] Loading input data')}
    input_data <- load_data(
        atac_path = atac_path,
        peaks_bed_path=peaks_bed_path,
        cell_metadata_path=cell_metadata_path,
        verbose=verbose
            )
    indata <- input_data$atac_data
    cellinfo <- input_data$cellinfo
    peakinfo <- input_data$peakinfo
    
    if (v) {message('[2] Sanitizing input')}
    # Exclude 0-sum columns and rows
    indata <- indata[Matrix::rowSums(indata) != 0,Matrix::colSums(indata) != 0] 
    peakinfo <- peakinfo[rownames(indata),]
    cellinfo <- cellinfo[colnames(indata),]
    
    if (v) {message('[3] Getting chromosome lengths')}
    # Load chromosome length
    chromosome_length <- load_chr_length()
    
    if (v) {message('[4] Creating CDS object')}
    
    # Create CDS for analysis
    cds <- create_cds(indata, binarize=TRUE, 
                      cellinfo=cellinfo, peakinfo=peakinfo)
    
    if (v) {message('[5] Detecting genes with Monocle3')}
    # Detect genes with monocle3
    cds <- monocle3::detect_genes(cds)
    
    if (v) {message('[6] Running basic QC')}
    # Remove low counts
    cds <- qc_cds(cds)
    
    if (v) {message('[7] Normalizing data')}
    # Normalize
    cds <- prepare_cds(cds, verbose=v)
    
    if (v) {message('[8] Plotting cells')}
    # Plot cells
    p <- monocle3::plot_cells(
        cds, reduction_method = 'UMAP',
        show_trajectory_graph = F,
        color_cells_by = plot_by,
        label_cell_groups=F)
    plot(p)
    
    if (v) {message('[9] Running Cicero analysis')}
    # Run cicero analysis
    conns <- analysis_cicero(cds, chromosome_length = chromosome_length, verbose=v)

    
    if (!is.null(save_dir)){
        
        dir.create(save_dir, showWarnings = F)
        
        if (v) {message('[SAVE] Saving connections & peaks')}
        fname <- 'cicero'
        
        if (!is.null(prefix)){
            fname <- paste(prefix, fname, sep='.')
        }
        
        all_peaks <- row.names(exprs(cds))
        write.csv(x = all_peaks, file = file.path(save_dir, paste0(fname, '_peaks.csv')))
        
        write.csv(x = conns, file = file.path(save_dir, paste0(fname, '_connections.csv')))
        saveRDS(conns,  file.path(save_dir, paste0(fname, '_connections.rds')))
        
        pdf(file.path(save_dir, paste0(fname, '_umap.pdf')))
        plot(p)
        dev.off()
        
    }
                
    return(conns)
    
}