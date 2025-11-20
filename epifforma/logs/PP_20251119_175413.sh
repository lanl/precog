#!/bin/bash -l
#SBATCH --error=./logs/arraylog/PP_%a_%N_%j.err.txt
#SBATCH --output=./logs/arraylog/PP_%a_%N_%j.log
#SBATCH --job-name=PP
#SBATCH --time=540
#SBATCH --mem=5G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-1
module load R
bash -c "$(head -n $SLURM_ARRAY_TASK_ID ./logs/PP_20251119_175413.cmd.txt | tail -n 1)"
