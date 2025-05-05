RunScrublet <- function(counts, file_name, input_dir = "./data/", output_dir = "./report/", script_path = "/home/sm2726/project2/ASD_iPSC/preprocess/scrublet.doublets.py") {
	##Execute the python. scrublet scripts to detect the doublets
	writeMM(counts, file = paste0(input_dir, file_name, "_matrix.mtx"))
	data.frame(ids = rownames(counts), name = rownames(counts), stringsAsFactors = FALSE) %>% write.table(., file = paste0(input_dir, file_name, "_genes.tsv"), sep = "\t", col.names = FALSE, row.names = FALSE, quote= FALSE)
	paste0("python ", script_path, " ", input_dir, file_name, "_matrix.mtx ", input_dir, file_name, "_genes.tsv ", output_dir, " ", file_name, " > ", input_dir, file_name, "_scrublet_out.txt") %>% system(.)
	file.remove(paste0(input_dir, file_name, "_matrix.mtx"));
	file.remove(paste0(input_dir, file_name, "_genes.tsv"));


	meta <- read.table(paste0(output_dir, file_name, "_scrublet_res.tsv"), sep = "\t", header = FALSE) %>% 
				t() %>%
				as.data.frame(., check.names = FALSE)
	rownames(meta) <- colnames(counts)
	colnames(meta) <- c("scrublet_score", "scrublet_assign")
    saveRDS(meta, file = paste0(input_dir, file_name, "_scrublet_meta.rds"))
    file.remove(paste0(output_dir, file_name, "_scrublet_res.tsv"))
}



##------------------------------------------------------------------------------------------------
## Prepare seurat object
seu_prepare <- function(counts, min.cells = 5, nfeatures = 2500, hvg.method = "vst", assay = "RNA") {
    inseu <- CreateSeuratObject(counts, meta.data = NULL, assay = "RNA", min.cells = min.cells, min.features = 0, names.field = 1, names.delim = "_");
    quality_genes <- get_genes(input_genes = rownames(inseu$RNA@data), gene_type = c("mito","ribo"), return_list = TRUE, revised = FALSE)
    inseu[["percent.mt"]] <- PercentageFeatureSet(inseu, pattern = NULL, features = quality_genes[["mito"]], col.name = NULL)
    inseu[["percent.ribo"]] <- PercentageFeatureSet(inseu, pattern = NULL, features = quality_genes[["ribo"]], col.name = NULL)


    #normalize the dataset if needed
    inseu <- NormalizeData(inseu, normalization.method = "LogNormalize", scale.factor = 10000, verbose = FALSE)
    if (!is.null(hvg.method)){
        inseu <- FindVariableFeatures(inseu, selection.method = hvg.method, nfeatures = nfeatures, verbose = FALSE)
    }
    return(inseu)
}

QuickCluster <- function(object, features, dims = 1:30, resolution = 1.2, do_clustering = TRUE) {
	object <- ScaleData(object, features = features) %>%
					RunPCA(., features = features, npcs = 50)

	if (do_clustering){
		object <- FindNeighbors(object, dims = dims) %>%
			FindClusters(., resolution = resolution, algorithm = 3)
	}

	set.seed(42)
	object <- RunUMAP(object, dims = dims, umap.method = "umap-learn", metric = "correlation")
	return(object)
}



#---------------------------------------------------------------------------------------------
## Function to perform the Seurat integration 
## Seurat integration (with out cell cycle correction)
Integratelist.seurat <- function(obj.list, hvg, file_name, input_dir = inputdir, inte.dims = 1:30, cluster.dims = 1:30, reference = NULL, do.cluster = FALSE) {
    if (length(obj.list) == 2){
        newseu <- merge(x = obj.list[[1]], y = obj.list[[2]])
    } else {
        newseu <- merge(x = obj.list[[1]], y = obj.list[2:length(obj.list)])
    }


    inte.slim.file <- paste0(input_dir, file_name, ".seurat.rds")
    if (!file.exists(inte.slim.file)){
        dg.anchors <- FindIntegrationAnchors(object.list = obj.list, dims = inte.dims, assay = NULL, anchor.features = hvg, reference = reference)
        seuinte <- IntegrateData(anchorset = dg.anchors, dims = inte.dims)
        DefaultAssay(seuinte) <- "integrated"
        seuinte <- ScaleData(seuinte, verbose = FALSE) %>%
                                RunPCA(., npcs = 50, verbose = FALSE)
        
        newseu[["pca"]] <- CreateDimReducObject(embeddings = seuinte$pca@cell.embeddings[colnames(newseu), ], loadings = seuinte$pca@feature.loadings, stdev = seuinte$pca@stdev, key = "PC_", assay = "RNA")
        rm(seuinte)
        newseu <- RunUMAP(newseu, dims = cluster.dims, umap.method = "umap-learn", metric = "correlation")
        if (do.cluster){
            newseu <- FindNeighbors(newseu, dims = cluster.dims,  k.param = 25) %>%
                    FindClusters(., resolution = 1.2, n.iter = 20, algorithm = 3)
        } else {
            newseu@meta.data$seurat_clusters <- "empty"
        }
        saveRDS(newseu, file = inte.slim.file)
    }else {
        newseu <- readRDS(file = inte.slim.file)
    }
    return(newseu)
}


## Seurat integration (with cell cycle correction)
Integratelist.seurat.cellcycle <- function(obj.list, hvg, file_name, input_dir = inputdir, inte.dims = 1:30, cluster.dims = 1:30, reference = NULL, do.cluster = FALSE) {
    if (length(obj.list) == 2){
        newseu <- merge(x = obj.list[[1]], y = obj.list[[2]])
    } else {
        newseu <- merge(x = obj.list[[1]], y = obj.list[2:length(obj.list)])
    }


    inte.slim.file <- paste0(input_dir, file_name, ".seurat.rds")
    if (!file.exists(inte.slim.file)){
        dg.anchors <- FindIntegrationAnchors(object.list = obj.list, dims = inte.dims, assay = NULL, anchor.features = hvg, reference = reference)
        seuinte <- IntegrateData(anchorset = dg.anchors, dims = inte.dims)
        DefaultAssay(seuinte) <- "integrated"
        seuinte <- ScaleData(seuinte, verbose = FALSE, vars.to.regress = c("S.Score", "G2M.Score")) %>%
                                RunPCA(., npcs = 50, verbose = FALSE)
        
        newseu[["pca"]] <- CreateDimReducObject(embeddings = seuinte$pca@cell.embeddings[colnames(newseu), ], loadings = seuinte$pca@feature.loadings, stdev = seuinte$pca@stdev, key = "PC_", assay = "RNA")
        rm(seuinte)
        newseu <- RunUMAP(newseu, dims = cluster.dims, umap.method = "umap-learn", metric = "correlation")
        if (do.cluster){
            newseu <- FindNeighbors(newseu, dims = cluster.dims,  k.param = 25) %>%
                    FindClusters(., resolution = 1.2, n.iter = 20, algorithm = 3)
        } else {
            newseu@meta.data$seurat_clusters <- "empty"
        }
        saveRDS(newseu, file = inte.slim.file)
    }else {
        newseu <- readRDS(file = inte.slim.file)
    }
    return(newseu)
}


LabelTransfer.seurat <- function(ref.list, query.list, hvg, inte.dims = 1:30, ref.anno) {
    ## For both the reference and the query dataset, perform data integration if needed
    FastIntegration <- function(object.list, inte.dims, hvg){
        if (length(object.list) == 1){
            seuinte <- object.list[[1]]
            DefaultAssay(seuinte) <- "RNA"
        } else {
            dg.anchors <- FindIntegrationAnchors(object.list = object.list, dims = inte.dims, assay = NULL, anchor.features = hvg)
            seuinte <- IntegrateData(anchorset = dg.anchors, dims = inte.dims)
            DefaultAssay(seuinte) <- "integrated"
        }
        
        return(seuinte)
    }


    ## For reference data
    ref <- FastIntegration(object.list = ref.list, inte.dims = inte.dims, hvg = hvg)
    query <- FastIntegration(object.list = query.list, inte.dims = inte.dims, hvg = hvg)


    lt.anchors <- FindTransferAnchors(reference = ref, query = query, dims = inte.dims, reduction = "cca",
                    reference.assay = "integrated", query.assay = "RNA")
    print(names(lt.anchors))
    predictions <- TransferData(anchorset = lt.anchors, refdata = ref@meta.data[, ref.anno], dims = inte.dims, weight.reduction = "cca")
    ##query <- AddMetaData(query, metadata = predictions)
    return(predictions)
}


## Transfer identities from reference data to query data
LabelTransfer.custom <- function(object, reduction = "pca", dims.use = 30, transfer_cols = "label", k = 20){
	meta_use <- object@meta.data
	data.use <- object[[reduction]]@cell.embeddings[, 1:dims.use]

	## The function to get the neighoring information
	get_most <- function(x){
		return(names(sort(table(x), decreasing = TRUE))[1])
	}


	##Assume the NA values are the inquiry cell
	anno_list <- list()
	for (anno in transfer_cols){
		## Here, we assume the NA annotaions represent query cells while non-NA annotations are reference cells
		ref_cells <- rownames(meta_use)[!is.na(meta_use[, anno])]
		inquiry_cells <- setdiff(rownames(meta_use), ref_cells)

		## Get the KNN cells for all the inquiry cells
		knn_cells <- FNN::get.knnx(data = data.use[ref_cells, ,drop = FALSE], query = data.use[inquiry_cells, ,drop = FALSE], k = k) %>%
						.$nn.index

		## Get the annotation for reference cells
		ref_anno <- meta_use[ref_cells, anno] %>% setNames(., ref_cells)

		## Do the annotation transfer
		if (class(ref_anno) %in% "character"){
			new_label <- apply(knn_cells, 1, function(x) get_most(ref_anno[x])) %>%
							setNames(., inquiry_cells) %>%
							c(., ref_anno)
		} else if (class(ref_anno) %in% "numeric") {
			new_label <- apply(knn_cells, 1, function(x) mean(ref_anno[x])) %>%
							setNames(., inquiry_cells) %>%
							c(., ref_anno)
		} else {
			stop(paste0("column ", column, " has unsupported object type"))
		}
		anno_list[[anno]] <- new_label[rownames(meta_use)]
	}
	## Add the label-transferred the metadata to the seurat object
	new_meta <- anno_list %>%
					as.data.frame(., stringsAsFactors = FALSE) %>%
					setNames(., paste0(names(anno_list), "_new"))
    return(new_meta)
}


LabelTransfer.parallel <- function(object, subset.list, reduction = "pca", dims.use = 30, transfer_cols = "label", k = 20, nreps = 100, nCores = 4){
    print(Sys.time())
    cl = makeCluster(nCores, outfile="")
    doParallel::registerDoParallel(cl);

    print("Start parallel label transfer")
    tem_idx <- NULL    

    label.all <- foreach(tem_idx = 1:length(subset.list), 
		                .combine = rbind,
                        .export = "LabelTransfer.custom",
		                .packages = c("FNN", "dplyr")) %dopar% {
				        
				        ## Some computing
                        obj <- object[, subset.list[[tem_idx]]]
                        lb_df <- LabelTransfer.custom(object = obj, reduction = reduction, dims.use = dims.use, transfer_cols = transfer_cols, k = k) %>%
                                    tibble::rownames_to_column("cell") %>%
                                    mutate(replicate = tem_idx)
                        print(paste0("Finish ", tem_idx, "/", length(subset.list), " iterations;"))

				return(lb_df)
		    	}
	stopCluster(cl)
	print("Finish parallel label transfer")
	print(Sys.time())
    return(label.all)
}



library(ggplot2)
library(ggsankey)
## meta: meta data containing x_col, next_col
## x_col: the first column name
## x_ord: the order of labels in the first column
## next_col: the second column name
## next_ord: the order of labels in the second column
## type: "sankey"
## space: space between nodes
## width: width of nodes

plot_ggsankey <- function(meta, x_col, next_col, x_ord = NULL, next_ord = NULL, type = "sankey", ...) {
    let_ord <- paste0(rep(letters, each = 9), rep(1:9, times = length(letters)))
    if (is.null(x_ord)){
        x_ord <- levels(as.factor(as.character(meta[, x_col])))
    }
    x_trans <- paste0(let_ord[1:length(x_ord)], x_ord) %>%
                    setNames(., x_ord)
    meta$newx <- x_trans[as.character(meta[, x_col])]



    if (is.null(next_ord)){
        next_ord <- levels(as.factor(as.character(meta[, next_col])))
    }
    next_trans <- paste0(let_ord[1:length(next_ord)], next_ord) %>%
                    setNames(., next_ord)
    meta$newnext <- next_trans[as.character(meta[, next_col])]


    df <- meta %>%
                make_long(newx, newnext) %>%
                mutate(nlabel = substring(node, 3))

    cols <- c(gg_color_hue(length(x_ord)), rep("#bdbdbd", length(next_trans))) %>%
                setNames(., c(setNames(x_trans, NULL), setNames(next_trans, NULL)))
    p <- ggplot(df, aes(x = x, next_x = next_x, node = node, next_node = next_node, fill = factor(node), label = nlabel)) +
            geom_sankey(flow.alpha = .6, node.color = "gray30", type = type, ...) +
            geom_sankey_text(size = 3, color = "black") +
            scale_fill_manual(values = cols) +
            theme_sankey(base_size = 18) +
            labs(x = NULL) +
            theme(legend.position = "none",  plot.title = element_text(hjust = .5), axis.text = element_blank())
    return(p)
}

