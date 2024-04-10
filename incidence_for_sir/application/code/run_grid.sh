#!/bin/bash -l

#SBATCH -J runGrid
#SBATCH -p general
#SBATCH -c 100
#SBATCH -t 10:00:00
#SBATCH --array=1-100
#SBATCH --output=SLURMOUT/featQuantiles_%A-%a.out

# Load software and activate environment (I loaded the modules into this env already):
module add R/4.1.2
R CMD BATCH "--no-save" make_data_grid.R ./logfiles/grid_logs_$SLURM_ARRAY_TASK_ID.Rout
