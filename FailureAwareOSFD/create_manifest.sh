#!/bin/bash
#SBATCH --job-name=create_manifest
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -p ccs6
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --mem=8G
#SBATCH --output=SLURMOUT/create_manifest_%j.out
#SBATCH --error=SLURMOUT/create_manifest_%j.err

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

Rscript R/create_work_manifest.R

# Verify the manifest was created
if [ ! -f data/work_manifest.csv ]; then
    echo "ERROR: work_manifest.csv was not created!"
    exit 1
fi

echo "Successfully created data/work_manifest.csv"
