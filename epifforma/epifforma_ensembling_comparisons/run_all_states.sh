#!/bin/bash -l

#SBATCH -J runStates
#SBATCH -p ccs6
#SBATCH -c 53
#SBATCH -t 05:00:00
#SBATCH --array=1-1
#SBATCH --mem=150G
#SBATCH --output=SLURMOUT/run_states_%A-%a.out

# Load software and activate environment (I loaded the modules into this env already):
module load R
R CMD BATCH "--no-save" R/run_epifforma.R ./logfiles/epifforma_by_state_$SLURM_ARRAY_TASK_ID.Rout

