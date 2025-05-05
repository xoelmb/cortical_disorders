#!/bin/bash
#SBATCH -N 1
#SBATCH -n 1
#SBATCH -c 20
#SBATCH -J label
#SBATCH -p pi_sestan
#SBATCH -t 1-00:00:00
#SBATCH --mem-per-cpu=50G
#SBATCH -o label_out.txt
#SBATCH -e label_err.txt


source ~/.bashrc;
conda activate bigdata;
Rscript nsc-pat-mes.label.transfer.R
conda deactivate

