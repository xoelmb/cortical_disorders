#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 16
#SBATCH -J lisi
#SBATCH -p day
#SBATCH -t 1-00:00:00
#SBATCH --mem-per-cpu=25G
#SBATCH -o lisi_out.txt
#SBATCH -e lisi_err.txt


source ~/.bashrc;
conda activate multiome;
Rscript nsc.inte.lisi.R
conda deactivate

