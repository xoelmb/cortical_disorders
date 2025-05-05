#!/bin/bash
#SBATCH --output dsq-mac.inte.full.task-%A_%1a-%N.out
#SBATCH --array 0-1%2
#SBATCH --job-name full
#SBATCH -c 16 -p bigmem --mem-per-cpu=90G -N 1 -t 1-00:00:00

# DO NOT EDIT LINE BELOW
/vast/palmer/apps/avx2/software/dSQ/1.05/dSQBatch.py --job-file /gpfs/gibbs/pi/sestan/sm2726/ASD_iPSC/qualitycontrol/mac.inte.full.task.txt --status-dir /gpfs/gibbs/pi/sestan/sm2726/ASD_iPSC/qualitycontrol

