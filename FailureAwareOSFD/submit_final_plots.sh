#!/bin/bash
#SBATCH --job-name=plots
##SBATCH -p general
#SBATCH -p ccs6
#SBATCH --qos=work
#SBATCH -c 100
#SBATCH -t 01:00:00
#SBATCH --mem=16G
#SBATCH --output=SLURMOUT/final_plots_%j.out

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

R CMD BATCH --no-save R/plot_regime_filled_sir_experiment.R "./logfiles/plot_regime_filled_sir_experiment.Rout"
R CMD BATCH --no-save R/plot_smoa_disease_comparisons.R "./logfiles/plot_smoa_disease_comparisons.Rout"
R CMD BATCH --no-save R/make_lauren_plots.R "./logfiles/lauren_plots.Rout"