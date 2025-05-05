# Select default mirror
options(repos=structure(c(CRAN="http://cran.r-project.org")))

# Install BiocManager for RcisTarget
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager", Ncpus=8)

# Install main package: RcisTarget
BiocManager::install("RcisTarget")

# Install extra packages
install.packages('doParallel', Ncpus=8)  # For anything and everything
install.packages('tidyverse', Ncpus=8)  # For anything and everything
install.packages('devtools', Ncpus=8)  # To install logger
devtools::install_github("sellorm/rlog")  # Logger
BiocManager::install("WGCNA") # For bicorrelation
BiocManager::install("DESeq2") # For DV data

# For plots
install.packages("UpSetR", Ncpus=8)
install.packages("ggpubr", Ncpus=8)
install.packages("ggh4x", Ncpus=8)
install.packages("igraph", Ncpus=8)
install.packages("ggvenn", Ncpus=8)
install.packages("ggsci", Ncpus=8)
devtools::install_github("jokergoo/ComplexHeatmap") 

# install.packages("openxlsx", Ncpus=8)


install.packages('IRkernel', Ncpus=8)  # To use in JupyterLab
IRkernel::installspec(name='cano_rcistarget', displayname='RcisTarget (R 4.2.2)', user=T, verbose=T)
