# Install and load required package
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("CoGAPS")
# Load the library
library("CoGAPS")
Npats=12# Npats=12 for Burke 2020;  Npats=11 for Ziller 2015
cells <- colnames(DATi)
genes <- rownames(DATi)
params <- new("CogapsParams")
params <- setParam(params, "sparseOptimization", FALSE)
params <- setParam(params, "nIterations", Niter)
params <- setParam(params, "nPatterns", Npats)
params <- setDistributedParams(params, nSets=1)
xxCoGAPS<-CoGAPS(data=DATi,params=params,geneNames=genes,sampleNames=cells, messages=TRUE,transposeData=FALSE,nPatterns=Npats,nIterations=50000,nThreads=1)
