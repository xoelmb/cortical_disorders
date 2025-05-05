## Some commonly used functions
gg_color_hue <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
}



map_gene <- function(gene_names, input_genes,ignore_case=TRUE){
    input_genes <- unique(input_genes)
  
    if (sum(grepl("\\|",gene_names))==length(gene_names)){
        if (sum(grepl("\\|",input_genes))==length(input_genes)){
              gene_id <- input_genes[input_genes %in% gene_names]
          }else{
            input_genes <- extract_field(input_genes=input_genes, 2, "|")
              gene_id <- unlist(parallel::mclapply(input_genes, function(i) ifelse(length(grep(paste0("\\|",i,"$"),gene_names, ignore.case= ignore_case))==1,gene_names[grep(paste0("\\|",i,"$"),gene_names, ignore.case= ignore_case)],"empty")))
              gene_id <- gene_id[gene_id != "empty"]
          }
    } else if(sum(grepl("\\|",gene_names))==0){
        input_genes <- extract_field(input_genes=input_genes, 2, "|")
          gene_id <- unlist(parallel::mclapply(input_genes, function(i) ifelse(length(grep(paste0("^",i,"$"),gene_names, ignore.case= ignore_case))==1,gene_names[grep(paste0("^",i,"$"),gene_names, ignore.case= ignore_case)],"empty")))
          gene_id <- gene_id[gene_id != "empty"]
    } else {
        stop("Inconsistent gene name format")
    }
    return(gene_id)
}



#Get the mito genes based on the inquiry genes
get_genes <- function(input_genes, gene_type = c("mito","ribo", "cc")[1], return_list = FALSE, revised = FALSE){
    gene_use <- list()
    if ("mito" %in% gene_type){
        mito.known <- map_gene(gene_names=input_genes, input_genes=c("ND1","ND2","COX1","COX2","ATP8","ATP6","COX3","ND3","ND4L","ND4","ND5","ND6","CYTB")) #Refseq annotation 103 (Macaque)
        mito.combine <- grep(pattern = "\\|MT-", x = input_genes, value = TRUE, ignore.case=TRUE) #All human mitochondria genes
        mito.single <- grep(pattern = "^MT-", x = input_genes, value = TRUE, ignore.case=TRUE) #All human mitochondria genes
        gene_use[["mito"]] <- unique(c(mito.known, mito.combine, mito.single))
    }

    if ("ribo" %in% gene_type){
        ribo.combine <- c(grep("\\|RPL", input_genes, value =TRUE, ignore.case=TRUE), grep("\\|RPS", input_genes, value =TRUE, ignore.case=TRUE), grep("\\|MRPL", input_genes, value =TRUE, ignore.case=TRUE), grep("\\|MRPS", input_genes, value =TRUE, ignore.case=TRUE))
        ribo.single <- c(grep("^RPL", input_genes, value =TRUE, ignore.case=TRUE), grep("^RPS", input_genes, value =TRUE, ignore.case=TRUE), grep("^MRPL", input_genes, value =TRUE, ignore.case=TRUE), grep("^MRPS", input_genes, value =TRUE, ignore.case=TRUE))
        gene_use[["ribo"]] <- c(ribo.combine, ribo.single)
    }

    if ("cc" %in% gene_type){
        if (!revised){
            regev.s.genes <- read.table(file=paste0("/gpfs/ycga/project/ysm/sestan/sm2726/public_data/CellCycle_correction/",'s_genes.txt'), header=FALSE, stringsAsFactors=FALSE)$V1
            regev.g2m.genes <- read.table(file=paste0("/gpfs/ycga/project/ysm/sestan/sm2726/public_data/CellCycle_correction/",'g2m_genes.txt'), header=FALSE, stringsAsFactors=FALSE)$V1
        } else {
            regev.s.genes <- read.table(file=paste0("/gpfs/ycga/project/ysm/sestan/sm2726/public_data/CellCycle_correction/",'s_genes_revised_04202019.txt'), header=FALSE, stringsAsFactors=FALSE)$V1
            regev.g2m.genes <- read.table(file=paste0("/gpfs/ycga/project/ysm/sestan/sm2726/public_data/CellCycle_correction/",'g2m_genes_revised_04202019.txt'), header=FALSE, stringsAsFactors=FALSE)$V1
        }
        
        gene_use[["s"]] <- map_gene(gene_names=input_genes, input_genes=regev.s.genes,ignore_case=TRUE)
        gene_use[["g2m"]] <- map_gene(gene_names=input_genes, input_genes=regev.g2m.genes,ignore_case=TRUE)
    }

    if (return_list){
        return(gene_use)
    } else {
        all_genes <- setNames(unlist(gene_use), NULL)
        return(all_genes)
    }
}


extract_field <- function(input_genes, field = 1, split_by = "_") {
    split_strings <- strsplit(input_genes, split_by, fixed=TRUE)
    if (is.numeric(field)){
        if (all(field > 0)){
            if (length(field) == 1){
                idx <- ifelse(sapply(split_strings, length) == 1, 1, field)
                output_strings <- sapply(1:length(split_strings), function(x) split_strings[[x]][idx[x]])
            } else {
                idx <- sapply(split_strings, length) == 1
                output_strings <- sapply(1:length(split_strings), function(x) ifelse(idx[x], input_genes[x], paste(split_strings[[x]][field], collapse = split_by)))
            }
        } else {
            if (length(field) == 1) {
                field = as.integer(abs(field))
                idx <- ifelse(sapply(split_strings, length) == 1, 1, field)
                output_strings <- sapply(1:length(split_strings), function(x) rev(split_strings[[x]])[idx[x]])
            } else {
                stop("currently doesnot support field with length >= 2")
            }
        }
    } else if (is.character(field)) {
        if (field == "rm_start"){
            idx <- ifelse(sapply(split_strings, length) == 1, 1, field)
            output_strings <- sapply(1:length(split_strings), function(x) paste(split_strings[[x]][-1],collapse=split_by))
        } else if (field == "rm_end") {
            output_strings <- sapply(1:length(split_strings), function(x) paste(rev(rev(split_strings[[x]])[-1]),collapse=split_by))
        } else {
            stop("Currently only support rm_start, rm_end for the character field")
        }
    }
    output_strings <- output_strings %>% setNames(., input_genes)
    return(output_strings)
}


#Function in Seurat package, useful
MinMax <- function (data, min, max) {
    data2 <- data
    data2[data2 > max] <- max
    data2[data2 < min] <- min
    return(data2)
}


## AUCell enrichment

library(AUCell)
GetModuleScore <- function (assay.data, features, nbin = 24, ctrl = 100, k = FALSE, seed = 42, method = c("seurat","aucell")[2], input_dir = new_inputdir, file_name, output_dir = outputdir, rethreshold_list = NULL, cellbin.size = 8000) {
     if (is.null(x = features)) {
        stop("Missing input feature list")
    }
    features <- lapply(X = features, FUN = function(x) {
        return(intersect(x = x, y = rownames(x = assay.data)))
    })
    cluster.length <- length(x = features) #number of feature list

    if (method == "seurat"){
        set.seed(seed = seed)
        pool <- rownames(x = assay.data)
        data.avg <- Matrix::rowMeans(x = assay.data[pool, ])
        data.avg <- data.avg[order(data.avg)]
        data.cut <- cut_number(x = data.avg + rnorm(n = length(data.avg))/1e+30, 
            n = nbin, labels = FALSE, right = FALSE)
        names(x = data.cut) <- names(x = data.avg)
        ctrl.use <- vector(mode = "list", length = cluster.length)
        for (i in 1:cluster.length) {
            features.use <- features[[i]]
            for (j in 1:length(x = features.use)) {
                ctrl.use[[i]] <- c(ctrl.use[[i]], names(x = sample(x = data.cut[which(x = data.cut == 
                    data.cut[features.use[j]])], size = ctrl, replace = FALSE)))
            }
        }
        ctrl.use <- lapply(X = ctrl.use, FUN = unique)
        ctrl.scores <- matrix(data = numeric(length = 1L), nrow = length(x = ctrl.use), 
            ncol = ncol(x = assay.data))
        for (i in 1:length(ctrl.use)) {
            features.use <- ctrl.use[[i]]
            ctrl.scores[i, ] <- Matrix::colMeans(x = assay.data[features.use, 
               ])
        }
        features.scores <- matrix(data = numeric(length = 1L), nrow = cluster.length, 
            ncol = ncol(x = assay.data))
        for (i in 1:cluster.length) {
            features.use <- features[[i]]
            data.use <- assay.data[features.use, , drop = FALSE]
            features.scores[i, ] <- Matrix::colMeans(x = data.use)
        }
        features.scores.use <- features.scores - ctrl.scores
        features.scores.use <- as.data.frame(x = t(x = features.scores.use))
        rownames(x = features.scores.use) <- colnames(x = assay.data)
        colnames(features.scores.use) <- names(features)

        return(features.scores.use)
    } else if (method == "aucell"){
        library(AUCell)
        if (!file.exists(paste0(input_dir, file_name, "_modulescore_auc_res.rds"))){
            #Split the cells to bins[Sometimes, the matrix is too large for rankings]
            if (ncol(assay.data) < cellbin.size){
                cellbin.size <- ceiling(ncol(assay.data)/2)
            }
            bin.ind <- ceiling(c(1:ncol(assay.data))/cellbin.size)
            max.bin <- max(bin.ind)

            auc_list <- lapply(1:max.bin, function(tembin) {
                tem_matrix <- assay.data[, bin.ind == tembin]
                tem_rankings <- AUCell_buildRankings(tem_matrix, nCores=1, plotStats=FALSE) 
                tem_auc <- AUCell_calcAUC(features, tem_rankings)#, aucMaxRank = 500)
                tem_aucmatrix <- t(as.matrix(getAUC(tem_auc)))
                rm(tem_matrix, tem_rankings)
                return(tem_aucmatrix)
                })

            hauc_matrix <- do.call(rbind, auc_list)

            if (length(auc_list) == 1){ #When the input is not a named list, the colnames will be "Geneset"
                colnames(hauc_matrix) <- names(features)
            }
            saveRDS(hauc_matrix, file = paste0(input_dir, file_name, "_modulescore_auc_res.rds"))
        } else {
            hauc_matrix <- readRDS(paste0(input_dir, file_name, "_modulescore_auc_res.rds"))
        }
        
        set.seed(seed)
        pdf(paste0(output_dir, file_name, "_modulescore_auc_auto_assignment.pdf"), paper="letter")
        par(mfrow=c(2,2)); cells_assignment <- AUCell_exploreThresholds(t(hauc_matrix), plotHist=TRUE, assign=TRUE) 
        dev.off()

        #Generate assignment matrix (rownames as cells)
        default_assign <- hauc_matrix * 0 #build an empty one
        for (gset in colnames(hauc_matrix)){
            default_assign[cells_assignment[[gset]]$assignment, gset] <- 1
        }

        outdata <- list(auc = hauc_matrix, auto = default_assign, custom = default_assign)

        if (!is.null(rethreshold_list)){
            pdf(paste0(output_dir, file_name, "_modulescore_auc_custom_assignment.pdf"), paper="letter")
            for (j in names(rethreshold_list)){
                AUCell_plotHist(t(hauc_matrix)[j,,drop = FALSE], aucThr=rethreshold_list[[j]])
                abline(v=rethreshold_list[[j]])
                cells_assignment[[j]]$assignment <- rownames(hauc_matrix)[hauc_matrix[, j] >= rethreshold_list[[j]]]
            }
            dev.off()

            ##default_assign <- hauc_matrix * 0 #build an empty one
            for (gset in colnames(hauc_matrix)){
                default_assign[, gset] <- 0
                default_assign[cells_assignment[[gset]]$assignment, gset] <- 1
            }
            outdata$custom <- default_assign
        }
        return(outdata)
    }
}



##------------------------------------------------------------------------------------------
## Visualization
##------------------------------------------------------------------------------------------
DimFig <- function(object, group.by, pt.size = 0.01, raster = TRUE, ncol = 3, output_dir = "./report/", file_name, plot.scale = 1, ...) {
	p_dim <- DimPlot(object = object, group.by = group.by, pt.size = pt.size, raster = raster, ncol = ncol, ...) &
			theme(plot.title = element_text(size = 8, face = "bold"), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())

	nrow <- ceiling(length(group.by)/ncol)
	if (raster){
		pdf(paste0(output_dir, file_name, ".pdf"), width = 6.5 * ncol * plot.scale, height = 5 * nrow * plot.scale)
	} else {
		jpeg(paste0(output_dir, file_name, ".jpeg"), width = 6.5 * ncol * plot.scale, height = 5 * nrow * plot.scale, unit = "in", res = 300)
	}
	print(p_dim)
	dev.off()
}



FeatureFig <- function(object, features, cols = c("lightgrey", "red"), output_dir = "./report/", file_name, plot.scale = 1, ncol = 4, raster = FALSE, ...) {
	plist <- lapply(features, function(gg) {
		p <- FeaturePlot(object = object, features = gg, cols = cols, ...) +
				theme(plot.title = element_text(size = 8, face = "bold"), axis.line = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), axis.title = element_blank())
		return(p)
		})

	nrow <- ceiling(length(features)/ncol)
	jpeg(paste0(output_dir, file_name, ".jpeg"), width = 5 * ncol * plot.scale, height = 5 * nrow * plot.scale, unit = "in", res = 300)
	pp <- patchwork::wrap_plots(plist, nrow = nrow, ncol = ncol) & theme(legend.position = "none")
	print(pp)
	dev.off()
}



##---------------------------------------------------------------------------------------------
## Function to perform the harmony integration
library(harmony)
InteAllSp.harmony <- function(object, split.by, hvg, file_name, input_dir = "./data/", dims = 1:30, theta = 2, lambda = 1, sigma = 0.1, resolution = 1.2, do.cluster = TRUE) {
    inte.slim.file <- paste0(input_dir, file_name, ".harmony.rds")
    if (!file.exists(inte.slim.file)){
        ## Do the Harmony integration
        object <- ScaleData(object, split.by = split.by, do.center = FALSE, features = hvg)%>%
                    RunPCA(., features = hvg, verbose = FALSE) %>%
                    RunHarmony(., group.by.vars = split.by, lambda = lambda, theta = theta, dims.use = dims, sigma = sigma)
        set.seed(0)
        object <- RunUMAP(object, dims = dims, reduction = "harmony", umap.method = "umap-learn", metric = "correlation")
        if (do.cluster){
        	object <- FindNeighbors(object, dims = dims, reduction = "harmony") %>%
                    FindClusters(., resolution = resolution, algorithm = 3)
        } else {
        	object@meta.data$seurat_clusters <- "empty"
        }
        saveRDS(object, file = inte.slim.file)
    }else {
        object <- readRDS(file = inte.slim.file)
    }
    return(object)
}


##------------------------------------------------------------------------------------------
## scHex visualization
##------------------------------------------------------------------------------------------
ifelse_check <- function(test, yes, no){if (test){return(yes)} else{ return(no) }}

make_hexbin_byregion <- function(sce, nbins = 80, dimension_reduction = "UMAP", split.by = "species"){
    all_regions <- levels(as.factor(sce@meta.data[, split.by]))
    sce_list <- list()
    for (region in all_regions){
        sub_sce <- sce[, rownames(sce@meta.data)[sce@meta.data[, split.by] == region]]
        sce_list[[region]] <- make_hexbin(sub_sce, nbins = nbins, dimension_reduction = dimension_reduction)
    }
    return(sce_list)
} 


plot_hex_byregion <- function(sce, feature, action = "mean", cols = colorRampPalette(c("lightgrey", "red"))(3), split.order = c("Human", "Chimpanzee", "Rhesus", "Marmoset"), assay = "RNA"){
    if (!feature %in% rownames(sce[[1]][[assay]])) {
        stop("Specify a valid assay.")
    }

    ##First prepare the data list
    all_regions <- ifelse_check(is.null(split.order), names(sce), split.order)
    exp_list <- lapply(all_regions, function(reg) as.numeric(GetAssayData(sce[[reg]], slot = "data")[feature, ])) %>% setNames(., all_regions)
    plot.data <- lapply(all_regions, function(reg) {
        x <- exp_list[[reg]]
        out <- sce[[reg]]@misc$hexbin[[2]]
        cID <- sce[[reg]]@misc$hexbin[[1]]
        df <- data.frame(out, 
                exp = schex:::.make_hexbin_function(x, action, cID), 
                region = reg, 
                stringsAsFactors = FALSE, check.names = FALSE)
        return(df)
        }) %>% do.call(rbind, .)
    

    maxexp <- max(plot.data$exp) %>% max(., 1)

    regp_list <- list()
    for (reg in all_regions){
        sub.data <- subset(plot.data, region == reg)
        ##cur.cols <- colorRampPalette(colors = cols)(nbks)[1:max(sub.data$scale.exp)]
        regp_list[[reg]] <- ggplot(data = sub.data) +
                    geom_hex(mapping = aes_string(x = "x", y = "y", fill = "exp"), stat = "identity") + 
                    theme_classic() + 
                    labs(title = reg) + 
                    coord_fixed() + 
                    scale_fill_gradientn(colors = cols, limits = c(0, maxexp)) + 
                    theme(legend.position = "bottom", axis.line = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(size = 12, hjust = 0.5, face = "bold"), legend.title = element_blank())
    }
    ##p <- plot_grid(plotlist = regp_list, nrow = 1, ncol = length(all_regions))
    
    return(regp_list)
}




plot_hex_cbn <- function(sce, feature, action = "mean", cols = viridis(3)){
    if (!feature %in% rownames(sce$RNA)) {
        stop("Specify a valid assay.")
    }

    ##First prepare the data list
    xx <- as.numeric(GetAssayData(sce, slot = "data")[feature, ])
    out <- sce@misc$hexbin[[2]]
    cID <- sce@misc$hexbin[[1]]
    plot.data <- data.frame(out, 
                exp = schex:::.make_hexbin_function(xx, action, cID), 
                stringsAsFactors = FALSE, check.names = FALSE)

    scale_exp <- function(x) {
        oo <- (x - mean(x))/sd(x)
        return(oo)
    }
    nbks <- 30
    plot.data <- plot.data %>%
                    mutate(scale.exp = scale_exp(exp)) %>%
                    mutate(scale.exp = MinMax(exp, -2.5, 2.5)) %>%
                    mutate(scale.exp = as.numeric(cut(scale.exp, nbks)))
    plot.data$scale.exp[is.na(plot.data$scale.exp)] <- 1


    p <- ggplot(data = plot.data) +
                    geom_hex(mapping = aes_string(x = "x", y = "y", fill = "scale.exp"), stat = "identity") + 
                    theme_classic() + 
                    labs(title = feature) + 
                    coord_fixed() + 
                    scale_fill_gradientn(colors = colorRampPalette(colors = cols)(nbks)) + 
                    theme(legend.position = "bottom", axis.line = element_blank(), axis.title = element_blank(), axis.ticks = element_blank(), axis.text = element_blank(), plot.title = element_text(size = 12, hjust = 0.5, face = "bold"), legend.title = element_blank())
    return(p)
}












