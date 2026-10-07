#!/bin/bash
#SBATCH --job-name=OSFD
##SBATCH -p general
#SBATCH -p ccs6
#SBATCH --qos=work
#SBATCH -c 100
#SBATCH -t 03-00:00:00
#SBATCH --array=1-1
#SBATCH --mem=150G
#SBATCH --output=SLURMOUT/run_states_%A-%a.out

# --- Try custom R first, fall back to module if unavailable ---
# Try to use custom R installation
if [ -d "$HOME/R-4.3.2/lib64/R" ] && [ -x "$HOME/R-4.3.2/lib64/R/bin/R" ]; then
    echo "Using custom R installation"
    module purge  # Avoid conflicts with module-provided R

    export R_HOME="$HOME/R-4.3.2/lib64/R"
    export PATH="$HOME/R-4.3.2/lib64/R/bin:$PATH"
    export LD_LIBRARY_PATH="$R_HOME/lib:$LD_LIBRARY_PATH"
    export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.3"
else
    # Fall back to module-provided R
    echo "Custom R not found, using module load R"
    module load R
fi

# --- Pre-run cleanup: delete all files under data/embeddings_gam_real/ ---
EMB_DIR="data/embeddings_gam_real"

if [[ -d "$EMB_DIR" ]]; then
  echo "[$(date)] Cleaning embedding dir: $EMB_DIR"
  # Delete files only (leave subdirs intact). Remove -print if too chatty.
  find "$EMB_DIR" -type f -print -delete
else
  echo "[$(date)] WARNING: Directory not found: $EMB_DIR"
fi

# Run the script with the custom R
R CMD BATCH --no-save R/run_experiment_smoa.R "./logfiles/bayesian_fit_${SLURM_ARRAY_TASK_ID}.Rout"
