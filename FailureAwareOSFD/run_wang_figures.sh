#!/bin/bash
#SBATCH --job-name=WangFigs
#SBATCH -p ccs6
#SBATCH --qos=work
#SBATCH -c 60
#SBATCH -t 48:00:00
#SBATCH --mem=200G
#SBATCH --output=SLURMOUT/wang_figures_%j.out

# ==============================================================================
# run_wang_figures.sh
# ==============================================================================
#
# Purpose:
#   Submit all three Wang extension experiment scripts to SLURM to reproduce
#   Figures 1, 2, and 8 from the paper.
#
# Scripts executed in order:
#   01_reproduce_wang_inverse_radius.R      - Reproduces Wang et al. (2024)
#   03_failure_stochastic_inverse_radius.R  - Extends with failures + stochasticity
#   04_plot_wang_extension_experiment.R     - Creates publication figures
#
# Output figures:
#   viz/failure_stochastic_diagnostic_scatter.png  (Figure 1)
#   viz/failure_stochastic_extension_metrics.png   (Figure 2)
#   viz/wang_reproduction_metrics.png              (Figure 8, appendix)
#
# Usage:
#   sbatch run_wang_figures.sh
#
# Author: AC Murph
# Date: Oct 2026
# ==============================================================================

# --- Setup R environment ---
# Try to use custom R installation first, fall back to module if unavailable
if [ -d "$HOME/R-4.3.2/lib64/R" ] && [ -x "$HOME/R-4.3.2/lib64/R/bin/R" ]; then
    echo "Using custom R installation"
    export R_HOME=$HOME/R-4.3.2/lib64/R
    export PATH=$R_HOME/bin:$PATH
    export LD_LIBRARY_PATH=$R_HOME/lib:$LD_LIBRARY_PATH
else
    echo "Custom R not found, loading R module"
    module load R/4.3.2
fi

# Verify R is available
which R
R --version

# --- Create output directories ---
mkdir -p data/wang_extension
mkdir -p viz
mkdir -p SLURMOUT

# --- Navigate to repository root ---
cd "$(dirname "$0")"
echo "Working directory: $(pwd)"

# ==============================================================================
# Step 1: Reproduce Wang et al. (2024) inverse-radius example
# ==============================================================================
echo ""
echo "========================================"
echo "Step 1/3: Running 01_reproduce_wang_inverse_radius.R"
echo "========================================"
echo "Start time: $(date)"
echo ""

Rscript R/01_reproduce_wang_inverse_radius.R

if [ $? -ne 0 ]; then
    echo "ERROR: Script 01 failed. Aborting pipeline."
    exit 1
fi

echo ""
echo "Step 1 completed successfully at $(date)"
echo ""

# ==============================================================================
# Step 2: Run failure + stochastic extension
# ==============================================================================
echo ""
echo "========================================"
echo "Step 2/3: Running 03_failure_stochastic_inverse_radius.R"
echo "========================================"
echo "Start time: $(date)"
echo ""

Rscript R/03_failure_stochastic_inverse_radius.R

if [ $? -ne 0 ]; then
    echo "ERROR: Script 03 failed. Aborting pipeline."
    exit 1
fi

echo ""
echo "Step 2 completed successfully at $(date)"
echo ""

# ==============================================================================
# Step 3: Generate all figures
# ==============================================================================
echo ""
echo "========================================"
echo "Step 3/3: Running 04_plot_wang_extension_experiment.R"
echo "========================================"
echo "Start time: $(date)"
echo ""

Rscript R/04_plot_wang_extension_experiment.R

if [ $? -ne 0 ]; then
    echo "ERROR: Script 04 failed. Figures may not be complete."
    exit 1
fi

echo ""
echo "Step 3 completed successfully at $(date)"
echo ""

# ==============================================================================
# Pipeline complete
# ==============================================================================
echo ""
echo "========================================"
echo "Pipeline completed successfully!"
echo "========================================"
echo "End time: $(date)"
echo ""
echo "Output figures (paper numbering):"
echo "  viz/failure_stochastic_diagnostic_scatter.png  (Figure 1)"
echo "  viz/failure_stochastic_extension_metrics.png   (Figure 2)"
echo "  viz/wang_reproduction_metrics.png              (Figure 8, appendix)"
echo ""
