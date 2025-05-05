#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 1
#SBATCH -J VIS
#SBATCH -p scavenge
#SBATCH -t 1-00:00:00
#SBATCH --mem-per-cpu=50G
#SBATCH -o vis_out.txt
#SBATCH -e vis_err.txt


source ~/.bashrc;
conda activate multiome;
Rscript vis.integration.R
conda deactivate

