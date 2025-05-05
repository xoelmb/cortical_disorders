args <- commandArgs(trailingOnly = TRUE)
library(dplyr)
library(Seurat)
library(Matrix)
source("./sc.fun.R")
source("./preprocess.fun.R")



## Manually curate the data
rawDir <- "/home/sm2726/project2/ASD_iPSC/cellranger_counts/"
sps <- c("D8_Ctrl_290", "D8_Ctrl_311", "D8_Ctrl_317", "D8_ASD_375", "D8_ASD_384", "D8_ASD_494")


sp <- sps[as.numeric(args[1])]
ctx <- Read10X(data.dir = setNames(paste0(rawDir, sp, "/outs/filtered_feature_bc_matrix/"), sp))
RunScrublet(counts = ctx, file_name = sp)


##sps <- c("D8_Ctrl_290", "D8_Ctrl_311", "D8_Ctrl_317", "D8_ASD_375", "D8_ASD_384", "D8_ASD_494")
##codes <- paste0("source ~/.bashrc; conda activate scmultiome; Rscript db.scrublet.R ", seq_along(sps), " > doublet_", sps, ".out; conda deactivate;")
##writeLines(codes, con = "db.scrublet.task.txt")

## Generate dSQ scripts
##module load dSQ
##dsq --job-file db.scrublet.task.txt --batch-file db.scrublet.sh -c 2 -p day --mem-per-cpu=15G -J db --max-jobs 6 -N 1





