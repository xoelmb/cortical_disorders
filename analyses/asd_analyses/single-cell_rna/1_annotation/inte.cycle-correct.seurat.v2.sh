#!/bin/bash
#SBATCH --output dsq-inte.cycle-correct.seurat.task-%A_%1a-%N.out
#SBATCH --array 0-2%3
#SBATCH --job-name seurat2
#SBATCH -c 6 -p day --mem-per-cpu=80G -N 1 -t 1-00:00:00

# DO NOT EDIT LINE BELOW
/vast/palmer/apps/avx2/software/dSQ/1.05/dSQBatch.py --job-file /gpfs/gibbs/pi/sestan/sm2726/ASD_iPSC/qualitycontrol/inte.cycle-correct.seurat.task.v2.txt --status-dir /gpfs/gibbs/pi/sestan/sm2726/ASD_iPSC/qualitycontrol

