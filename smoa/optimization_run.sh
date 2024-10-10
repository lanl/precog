#!/bin/bash -l

#SBATCH -J murph_opt
#SBATCH -p general
#SBATCH -c 100
#SBATCH -t 02-00:00:00
#SBATCH --mem-per-cpu=64G
#SBATCH --qos=long
#SBATCH --output=SLURMOUT/murph_%A-%a.out

# Load software and activate environment (I loaded the modules into this env already):
module load R
R CMD BATCH "--no-save" R/murph_optimize_smoa.R ./logfiles/murph_optimize_smoa_$SLURM_ARRAY_TASK_ID.Rout

