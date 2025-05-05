#!/bin/bash
#SBATCH --output dsq-nsc.inte.sub.task-%A_%1a-%N.out
#SBATCH --array 0-4%5
#SBATCH --job-name umap
#SBATCH -c 4 -p day --mem-per-cpu=20G -N 1 -t 1-00:00:00

# DO NOT EDIT LINE BELOW
/vast/palmer/apps/avx2/software/dSQ/1.05/dSQBatch.py --job-file /gpfs/gibbs/pi/sestan/sm2726/ASD_iPSC/qualitycontrol/nsc.inte.sub.task.txt --status-dir /gpfs/gibbs/pi/sestan/sm2726/ASD_iPSC/qualitycontrol

