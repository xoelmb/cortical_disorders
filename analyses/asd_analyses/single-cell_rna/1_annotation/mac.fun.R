OrthoTransfer <- function(object, from = "mmulatta", to = "hsapiens") {
    ## Get the updated count and normalized data
    ctx_new <- orthogene::convert_orthologs(gene_df = object$RNA@counts,
                                            gene_input = "rownames", 
                                            gene_output = "rownames", 
                                            input_species = from,
                                            output_species = to,
                                            non121_strategy = "drop_both_species",
                                            method = "gprofiler")
    data_new <- orthogene::convert_orthologs(gene_df = object$RNA@data,
                                            gene_input = "rownames", 
                                            gene_output = "rownames", 
                                            input_species = from,
                                            output_species = to,
                                            non121_strategy = "drop_both_species",
                                            method = "gprofiler")


    ## Create a new seurat object
    newobj <- CreateSeuratObject(ctx_new, meta.data = object@meta.data, assay = "RNA", min.cells = 0, min.features = 0, names.field = 1, names.delim = "_")
    newobj[["RNA"]] <- CreateAssayObject(data = data_new)
    return(newobj)
}

IdentfiyHVG <- function(object, level1.by = "period", level2.by = "inte.batch", nf = 1500) {
    all_prds <- levels(as.factor(object@meta.data[, level1.by]))
    hvgs <- lapply(all_prds, function(xx) {
        ss <- object[, object@meta.data[, level1.by] %in% xx]
        if (length(unique(ss@meta.data[, level2.by])) == 1){
            hvg <- FindVariableFeatures(ss, nfeatures = nf) %>%
                    VariableFeatures()
        } else {
            hvg <- SplitObject(ss, split.by = level2.by) %>%
                    lapply(., function(mm) FindVariableFeatures(mm, nfeatures = nf + 200)) %>%
                    SelectIntegrationFeatures(., nfeatures = nf)
        }
        return(hvg)
        }) %>%
        unlist() %>%
        unique()
    return(hvgs)
}

OrthoHVGs <- function(hvg){
    df <- data.frame(row.names = hvg, value = rep(1, length(hvg)), stringsAsFactor = FALSE)
    df_new <- orthogene::convert_orthologs(gene_df = df,
                                            gene_input = "rownames", 
                                            gene_output = "rownames", 
                                            input_species = "mmulatta",
                                            output_species = "hsapiens",
                                            non121_strategy = "drop_both_species",
                                            method = "gprofiler")
    hvgs_new <- rownames(df_new)
    return(hvgs_new)
}


