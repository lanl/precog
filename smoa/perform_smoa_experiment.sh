#!/bin/bash -l

#SBATCH -J runStates
#SBATCH -p general
#SBATCH -c 53
#SBATCH -t 05:00:00
#SBATCH --array=1-1
#SBATCH --mem=150G
#SBATCH --output=SLURMOUT/run_states_%A-%a.out

# Change to the directory where sbatch was submitted from
cd "$SLURM_SUBMIT_DIR"

# Create output directories if they don't exist
mkdir -p SLURMOUT
mkdir -p logfiles

# Cleanup old SLURMOUT files from previous runs (keep only current job ID)
if [ -n "$SLURM_ARRAY_JOB_ID" ]; then
  echo "Cleaning up old SLURMOUT files (keeping job ID: $SLURM_ARRAY_JOB_ID)"
  find SLURMOUT/ -type f -name "run_states_*.out" ! -name "run_states_${SLURM_ARRAY_JOB_ID}-*.out" -delete
elif [ -n "$SLURM_JOB_ID" ]; then
  echo "Cleaning up old SLURMOUT files (keeping job ID: $SLURM_JOB_ID)"
  find SLURMOUT/ -type f -name "run_states_*.out" ! -name "run_states_${SLURM_JOB_ID}-*.out" -delete
fi

# Load software and activate environment (I loaded the modules into this env already):
module load R

R CMD BATCH "--no-save" R/run_smoa.R ./logfiles/smoa_by_state_$SLURM_ARRAY_TASK_ID.Rout

