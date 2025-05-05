#!/bin/bash
#SBATCH --output dsq-db.scrublet.task-%A_%1a-%N.out
#SBATCH --array 0-5%6
#SBATCH --job-name db
#SBATCH -c 2 -p day --mem-per-cpu=15G -N 1

# DO NOT EDIT LINE BELOW
/vast/palmer/apps/avx2/software/dSQ/1.05/dSQBatch.py --job-file /gpfs/gibbs/pi/sestan/sm2726/ASD_iPSC/qualitycontrol/db.scrublet.task.txt --status-dir /gpfs/gibbs/pi/sestan/sm2726/ASD_iPSC/qualitycontrol

