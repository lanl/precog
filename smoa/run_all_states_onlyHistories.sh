#!/bin/bash -l

#SBATCH -J runOnlyH
#SBATCH -p general
#SBATCH -c 1
#SBATCH --array=1-1
#SBATCH -t 04:00:00
#SBATCH --mem=120G
#SBATCH --output=SLURMOUT/onlyhistory_%A-%a.out

# Load software and activate environment (I loaded the modules into this env already):
module load R
R CMD BATCH "--no-save" R/run_smoa_w_only_history.R ./logfiles/onlyHistory_smoa_by_state_$SLURM_ARRAY_TASK_ID.Rout

