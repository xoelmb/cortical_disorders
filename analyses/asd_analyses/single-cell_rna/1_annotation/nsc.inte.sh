#!/bin/bash
#SBATCH --output dsq-nsc.inte.task-%A_%1a-%N.out
#SBATCH --array 0-3%4
#SBATCH --job-name nsc
#SBATCH -c 8 -p day --mem-per-cpu=40G -N 1 -t 1-00:00:00

# DO NOT EDIT LINE BELOW
/vast/palmer/apps/avx2/software/dSQ/1.05/dSQBatch.py --job-file /gpfs/gibbs/pi/sestan/sm2726/ASD_iPSC/QC/nsc.inte.task.txt --status-dir /gpfs/gibbs/pi/sestan/sm2726/ASD_iPSC/QC

