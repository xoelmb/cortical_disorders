#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 4
#SBATCH -J ctm
#SBATCH -p day
#SBATCH -t 1-00:00:00
#SBATCH --mem-per-cpu=30G
#SBATCH -o ctm_out.txt
#SBATCH -e ctm_err.txt


source ~/.bashrc;
conda activate multiome;
Rscript counts.hvg.extract.R
conda deactivate




