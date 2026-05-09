#!/bin/bash -l

#SBATCH -J runMutAnti
#SBATCH -p ccs6
#SBATCH --qos=work
#SBATCH -c 100
#SBATCH -t 07-00:00:00
#SBATCH --array=1-1
#SBATCH --mem=150G
#SBATCH --output=SLURMOUT/run_states_%A-%a.out

# Load software and activate environment (I loaded the modules into this env already):
module load R
R CMD BATCH "--no-save" R/run_mutantigen_parallel_wOrigYML.R ./logfiles/mutantigen_$SLURM_ARRAY_TASK_ID.Rout


