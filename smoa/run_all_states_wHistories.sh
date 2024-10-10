#!/bin/bash -l

#SBATCH -J runSnH
#SBATCH -p general
#SBATCH -c 1
#SBATCH --array=1-1
#SBATCH -t 02-00:00:00
#SBATCH --mem-per-cpu=64G
#SBATCH --qos=long
#SBATCH --output=SLURMOUT/history_%A-%a.out

# Load software and activate environment (I loaded the modules into this env already):
module load R
R CMD BATCH "--no-save" R/run_smoa_w_online_history.R ./logfiles/history_and_synthetic_smoa_by_state_$SLURM_ARRAY_TASK_ID.Rout

