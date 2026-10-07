#!/bin/bash
#SBATCH --job-name=prep_embeddings
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH -p ccs6
#SBATCH --cpus-per-task=100
#SBATCH --time=04:00:00
#SBATCH --output=SLURMOUT/prep_embeddings_%j.out
#SBATCH --error=SLURMOUT/prep_embeddings_%j.err

module load R/4.4.0-fosscuda-2020b

Rscript R/prepare_synthetic_embeddings.R

# Verify the output file was created
if [ ! -f data/synthetic_embeddings.RData ]; then
    echo "ERROR: synthetic_embeddings.RData was not created!"
    exit 1
fi

echo "Successfully created data/synthetic_embeddings.RData"

# Clean up the stochastic_sir directory
echo "Removing data/stochastic_sir directory..."
rm -rf data/stochastic_sir
echo "Cleanup complete."
