#!/bin/bash -l
#SBATCH --error=./features2weights/arraylog/PP_%a_%N_%j.err.txt
#SBATCH --output=./features2weights/arraylog/PP_%a_%N_%j.log
#SBATCH --job-name=PP
#SBATCH --time=540
#SBATCH --mem=5G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-1
module load R
bash -c "$(head -n $SLURM_ARRAY_TASK_ID ./features2weights/PP_20251119_175433.cmd.txt | tail -n 1)"
