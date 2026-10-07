#!/bin/bash
#SBATCH --job-name=neural_sbi
#SBATCH --output=logs/slurm/neural_sbi_%j.log
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --time=144:00:00

# Neural SBI Pipeline - SLURM Submission Script for LinusBlack
# This script manages the entire snakemake workflow with cluster execution

echo "=========================================="
echo "Neural SBI Pipeline - Cluster Execution"
echo "Job ID: $SLURM_JOB_ID"
echo "Node: $SLURM_NODELIST"
echo "Started: $(date)"
echo "=========================================="

# Create log directory
mkdir -p logs/slurm

# Load conda/mamba
# Adjust this path if your conda installation is elsewhere
source ~/miniforge3/etc/profile.d/conda.sh
conda activate neural_sbi

# Verify environment
echo "Python: $(which python)"
echo "Snakemake: $(which snakemake)"
echo "R: $(which R)"
echo ""
echo "Configuration files:"
echo "  - Loaded internally by Snakefile (4 domain-specific config files)"
echo "  - config/cluster_config.yaml (resource specifications for SLURM)"
echo ""

# Run snakemake with cluster execution
# Each rule will be submitted as a separate SLURM job
# IMPORTANT: --cluster-config specifies per-rule resources (time, threads)
# NOTE: Memory allocation removed - not supported on this cluster
# NOTE: Config files loaded by Snakefile's configfile: directives (simulation, training, evaluation, system)
snakemake \
    --cluster-config config/cluster_config.yaml \
    --cluster "sbatch \
        --parsable \
        --cpus-per-task={cluster.threads} \
        --time={cluster.time} \
        --output=logs/slurm/{rule}_{wildcards}_%j.log" \
    --jobs 100 \
    --latency-wait 60 \
    --restart-times 3 \
    --keep-going \
    --rerun-incomplete \
    --printshellcmds \
    --reason \
    --use-conda \
    --default-resources \
        time=\"04:00:00\" \
        threads=1 \
    2>&1 | tee logs/snakemake_$(date +%Y%m%d_%H%M%S).log

EXIT_CODE=$?

echo "=========================================="
echo "Pipeline completed with exit code: $EXIT_CODE"
echo "Finished: $(date)"
echo "=========================================="

# Summary of results
if [ $EXIT_CODE -eq 0 ]; then
    echo "✓ SUCCESS: All jobs completed successfully"
    echo ""
    echo "Output files:"
    echo "  - results/summaries/train_aggregated.npz"
    echo "  - results/summaries/test_aggregated.npz"
    echo "  - results/models/model_y_only/trained_model.pkl"
    echo "  - results/models/model_y_z/trained_model.pkl"
    echo "  - results/evaluation/*.csv"
else
    echo "✗ FAILED: Pipeline encountered errors"
    echo "Check logs in logs/slurm/ for details"
fi

exit $EXIT_CODE
