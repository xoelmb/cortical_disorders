library(dplyr)
library(tibble)
library(lisi)
library(doParallel)
library(parallel)
library(foreach)

## Balance the number of the cells between datasets and calculate the lisi scores per cell & clusters
## Identify the non-aligned cell types via LISI scores

## Functions to calculate the LISI scores
## meta: cell meta data, a data.frame, rows as cells and columns as meta
## emb: low dimensional embeddings, rows as cells and columns as dimensions (assuming all dimensions will be used). Cells should have the same order as the meta input
## split.by: batch column name in the meta object
## perplexity: k for KNN
## equal.size: whether to balance the cell numbers acrosss different batches when computing the the LISI scores
## seed.use: this is only used when equal.size = TRUE
CalcLISI.default <- function(meta, emb, split.by, perplexity = 20, equal.size = TRUE, seed.use = 0) {
    ## check if the meta data contains split.by column
    if (!split.by %in% colnames(meta)){
        stop(paste0("The meta data doesn't contain the column: ", split.by))
    }

    ## To have balanced cell number between all batches
    if (equal.size) {
        meta$index_temp <- 1:nrow(meta)
        batch.ind <- split(meta$index_temp, as.character(meta[, split.by]))
        batch.size <- sapply(batch.ind, length)

        set.seed(seed.use)
        new_idx <- c()
        for (bat in names(batch.ind)){
            if (batch.size[bat] > min(batch.size)){
                subcells <- sample(batch.ind[[bat]], min(batch.size), replace = FALSE)
            } else {
                subcells <- batch.ind[[bat]]
            }
            new_idx <- c(new_idx, subcells)
        }

        ## update the meta and emb
        emb <- emb[new_idx, , drop = FALSE]
        meta <- meta[new_idx, , drop = FALSE]
    }

    ## Compute LISI scores using the updated index
	lisi_score <- compute_lisi(X = emb, meta_data = meta, label_colnames = split.by, perplexity = perplexity, nn_eps = 0) %>%
					setNames(., "LISI") %>%
                    tibble::rownames_to_column("cell")
	return(lisi_score)
}


## Functions to calculate the LISI scores in parallel
## meta: cell meta data, a data.frame, rows as cells and columns as meta
## emb: low dimensional embeddings, rows as cells and columns as dimensions (assuming all dimensions will be used). Cells should have the same order as the meta input
## split.by: batch column name in the meta object
## perplexity: k for KNN
## n.reps: number of parallel replicates, default 50
## nCores: number of cores for parallel computing, default 4
CalcLISI.parallel <- function(meta, emb, split.by, perplexity = 15, n.reps = 50, nCores = 4, query.name = "query") {
	## Generate a list of random seeds
	set.seed(42)
	seeds_use <- sample(1:10000, n.reps, replace = FALSE) 

	print(Sys.time())
    cl = makeCluster(nCores, outfile="")
	doParallel::registerDoParallel(cl);

    print("Start LISI score calculation")
	tem_idx <- NULL
	scores <- foreach(tem_idx = 1:length(seeds_use), 
		                .combine = rbind,
                        .export = "CalcLISI.default",
		                .packages = c("lisi", "dplyr")) %dopar% {
				        
				        ## Some computing
						seed.use <- seeds_use[tem_idx]
                        lscores <- CalcLISI.default(meta = meta, emb = emb, split.by = split.by, perplexity = perplexity, equal.size = TRUE, seed.use = seed.use) %>%
									mutate(replicate = paste0("rep_", tem_idx))
                        print(paste0("Finish ", tem_idx, "/", length(seeds_use), " iterations;"))

				return(lscores)
		    	}
	stopCluster(cl)
	print("Finish LISI score calculation")
	print(Sys.time())

    ## Add the split.by information back
    scores[[split.by]] <- meta[[split.by]][match(scores[["cell"]], rownames(meta))]

    ## Filter to query cells
    scores <- filter(scores, !!sym(split.by) %in% query.name)
    ## Average lisi scores for each cell
    scores.avg <- scores %>%
            group_by(cell) %>%
            summarize(avgLISI = mean(LISI)) %>%
            ungroup()
	return(list(raw = scores, avg = scores.avg))
}




LISIFlag <- function(data, LISI.col, group.by, cell.thre = 0.075, cluster.thre = 0.125){
    ## Check fisher's exact p values
    all.cls <- unique(data[[group.by]])
    min.cell.lisi <- quantile(data[[LISI.col]], probs = cell.thre) %>%
                        max(., 1.02) %>% ## in case the quantile give rise to very small value such as 1.00001
                        min(., 1.05) ## incase the min.cell.lisi is too big
                        
    test.data <- filter(data, !!sym(LISI.col) < min.cell.lisi)
    bg.data <- filter(data, !!sym(LISI.col) >= min.cell.lisi)

    pvals <- sapply(all.cls, function(cls){
        a1 <- sum(test.data[[group.by]] %in% cls) ## Small LISI (Yes) + Celltype cls (yes)
        b1 <- nrow(test.data) - a1 ## Small LISI (Yes) + Celltype cls (No)
        c1 <- sum(bg.data[[group.by]] %in% cls) ## Small LISI (NO) + Celltype cls (yes)
        d1 <- nrow(bg.data) - c1 ## Small LISI (NO) + Celltype cls (NO)
        
        mat <- matrix(c(a1, b1, c1, d1), nrow = 2, ncol = 2, byrow = FALSE)
        pp <- fisher.test(x = mat, alternative = "greater")$p.val
        return(pp)
    }) %>%
        sort() %>%
        p.adjust(., method = "bonferroni")
    
    ## average LISI per cluster
    out.cluster <- data %>%
                group_by(!!sym(group.by)) %>%
                summarize(avgLISI = mean(!!sym(LISI.col))) %>%
                ungroup() %>%
                mutate(pval = pvals[match(!!sym(group.by), names(pvals))]) %>%
                arrange(avgLISI, pval)
    
    ## cluster to flag
    min.cluster.lisi <- quantile(data[[LISI.col]], probs = cluster.thre) %>%
                        max(., 1.05) %>%
                        min(., 1.1) ## in case the threshold is too bit or small
    flag.cls <- out.cluster[[group.by]][out.cluster[["avgLISI"]] < min.cluster.lisi & 
                                        out.cluster[["pval"]] < 0.01]
    ## cellular level output
    data[["cell.lisi.flag"]] <- ifelse(data[[LISI.col]] < min.cell.lisi, 1, 0)
    data[["cluster.lisi.flag"]] <- ifelse(data[[group.by]] %in% flag.cls, 1, 0)
    return(list(cell = data, cluster = out.cluster))
}


LISI.wrapper <- function(meta = meta, emb = pca, split.by = "dataset", query.name, group.by,
            query.split.by = NULL,
            perplexity = 20, n.reps = 100, nCores = 4, 
            cell.thre = 0.075, cluster.thre = 0.125){
    ## First obtain the lisi results
    ## The CalcLISI function assumes more reference cells than query cells
    ## So here, we split the query meta to multiple subsets if there are more query than refs.
    query.cells <- which(meta[[split.by]] %in% query.name)
    ref.cells <- setdiff(1:nrow(meta), query.cells)

    if (length(query.cells) <= (length(ref.cells) * 0.95)){
        lisi <- CalcLISI.parallel(meta = meta, emb = emb, split.by = split.by, perplexity = perplexity, n.reps = n.reps, nCores = nCores) 
        lisi.avg <- lisi[["avg"]]
    } else {
        ## Randomly split the query cells to multiple groups
        if (is.null(query.split.by)) {
            set.seed(0)
            nbins <- ceiling(length(query.cells)/(length(ref.cells) * 0.95))
            bin.size <- ceiling(length(query.cells)/nbins)
            query.bins <- ceiling(sample(1:length(query.cells), replace = FALSE)/bin.size)
        } else {## Split the query cells based on the query split.by
            query.meta <- meta[query.cells, , drop = FALSE]
            query.bins <- as.numeric(as.factor(query.meta[["query.split.by"]]))
        }


        lisi.avg.list <- list()
        for (ii in seq_len(nbins)){
            new.cells <- union(ref.cells, query.cells[query.bins == ii])
            new.meta <- meta[new.cells, , drop = FALSE]
            new.emb <- emb[new.cells, , drop = FALSE]
            new.lisi <- CalcLISI.parallel(meta = new.meta, emb = new.emb, split.by = split.by, perplexity = perplexity, n.reps = n.reps, nCores = nCores, query.name = query.name) 
            lisi.avg.list[[ii]] <- new.lisi[["avg"]]
        }
        lisi.avg <- do.call(rbind, lisi.avg.list)
    }


    ## Add cluster identity 
    lisi.avg <- lisi.avg %>%
            mutate(!!sym(group.by) := meta[[group.by]][match(cell, rownames(meta))])
    return(lisi.avg)
    #lisi.summary <- LISIFlag(data = lisi.avg, 
    #        LISI.col = "avgLISI", group.by = group.by, 
    #        cell.thre = 0.075, cluster.thre = 0.125)
    #return(list(avg = lisi.avg, summary = lisi.summary))
}

