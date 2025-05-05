#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -J label
#SBATCH -p day
#SBATCH -t 1-00:00:00
#SBATCH --mem-per-cpu=50G
#SBATCH -o label_out.txt
#SBATCH -e label_err.txt


source ~/.bashrc;
conda activate multiome;
###Rscript nsc.seurat.transfer.R
Rscript nsc.label.transfer.R
conda deactivate

