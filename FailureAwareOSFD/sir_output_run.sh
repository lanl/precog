#!/bin/bash
#SBATCH --job-name=SIRout
#SBATCH -p general
##SBATCH -p ccs6
##SBATCH --qos=work
#SBATCH -c 100
#SBATCH -t 08:00:00
#SBATCH --mem=150G
#SBATCH --output=SLURMOUT/run_states_%j.out

# Set nsuccesses value (can be overridden via: sbatch --export=NSUCCESSES=10000 sir_output_run.sh)
NSUCCESSES=${NSUCCESSES:-1000}

# Set number of replicates to run
NREP=20

# ==============================================================================
# Clean up old logfiles and SLURMOUT files
# ==============================================================================
echo "==========================================="
echo "Cleaning up old files..."
echo "==========================================="

# Remove all logfiles from previous runs
if [ -d "logfiles" ]; then
    echo "Removing old logfiles..."
    rm -f logfiles/sir_part*.Rout
    rm -f logfiles/coverage_*.Rout
    rm -f logfiles/sir_plot_*.Rout
    echo "  Logfiles cleaned"
fi

# Remove SLURMOUT files from other jobs (keep current job only)
if [ -d "SLURMOUT" ] && [ ! -z "${SLURM_JOB_ID}" ]; then
    echo "Removing old SLURMOUT files (keeping current job ${SLURM_JOB_ID})..."
    find SLURMOUT/ -type f -name "run_states_*.out" ! -name "run_states_${SLURM_JOB_ID}.out" -delete 2>/dev/null || true
    echo "  SLURMOUT files cleaned"
fi

echo "Cleanup complete"
echo ""

# ==============================================================================
# --- Try custom R first, fall back to module if unavailable ---
# ==============================================================================
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

# Run all four parts sequentially for each replicate
echo "Using NSUCCESSES=${NSUCCESSES}"
echo "Running ${NREP} replicates (each replicate runs Parts 1-4)"
echo ""

# Initialize timing CSV file for all parts
TIMING_FILE="./data/runtimes_all_parts_n${NSUCCESSES}.csv"
echo "replicate,method,runtime_seconds,runtime_minutes,nsuccesses" > ${TIMING_FILE}

# Loop through replicates
for ((rep=1; rep<=NREP; rep++)); do
    echo ""
    echo "==========================================="
    echo "Starting replicate ${rep} of ${NREP}"
    echo "==========================================="
    echo ""

    # --- Part 1: Inverse SIR Mapping ---
    echo "[Part 1/4] Running Inverse SIR Mapping (rep ${rep})..."
    START_TIME=$(date +%s)
    Rscript R/run_sir_part1_inverse_mapping.R ${NSUCCESSES} ${rep} > "./logfiles/sir_part1_rep${rep}_${SLURM_JOB_ID}.Rout" 2>&1
    END_TIME=$(date +%s)
    RUNTIME_SEC=$((END_TIME - START_TIME))
    RUNTIME_MIN=$(awk "BEGIN {printf \"%.2f\", ${RUNTIME_SEC}/60}")
    echo "  Part 1 completed in ${RUNTIME_SEC} seconds (${RUNTIME_MIN} minutes)"
    echo "${rep},SIR_Maps,${RUNTIME_SEC},${RUNTIME_MIN},${NSUCCESSES}" >> ${TIMING_FILE}

    # --- Part 2: Basic LHS ---
    echo "[Part 2/4] Running Basic LHS (rep ${rep})..."
    START_TIME=$(date +%s)
    Rscript R/run_sir_part2_basic_lhs.R ${NSUCCESSES} ${rep} > "./logfiles/sir_part2_rep${rep}_${SLURM_JOB_ID}.Rout" 2>&1
    END_TIME=$(date +%s)
    RUNTIME_SEC=$((END_TIME - START_TIME))
    RUNTIME_MIN=$(awk "BEGIN {printf \"%.2f\", ${RUNTIME_SEC}/60}")
    echo "  Part 2 completed in ${RUNTIME_SEC} seconds (${RUNTIME_MIN} minutes)"
    echo "${rep},Basic_LHS,${RUNTIME_SEC},${RUNTIME_MIN},${NSUCCESSES}" >> ${TIMING_FILE}

    # --- Part 3: OSFD Basic ---
    echo "[Part 3/4] Running OSFD Basic with max_budget=${NSUCCESSES} (rep ${rep})..."
    START_TIME=$(date +%s)
    Rscript R/run_sir_part3_osfd_basic.R ${NSUCCESSES} ${NSUCCESSES} ${rep} > "./logfiles/sir_part3_rep${rep}_${SLURM_JOB_ID}.Rout" 2>&1
    END_TIME=$(date +%s)
    RUNTIME_SEC=$((END_TIME - START_TIME))
    RUNTIME_MIN=$(awk "BEGIN {printf \"%.2f\", ${RUNTIME_SEC}/60}")
    echo "  Part 3 completed in ${RUNTIME_SEC} seconds (${RUNTIME_MIN} minutes)"
    echo "${rep},OSFD_Basic,${RUNTIME_SEC},${RUNTIME_MIN},${NSUCCESSES}" >> ${TIMING_FILE}

    # Note: This last comparison model became computationally infeasible. 
    # # --- Part 4: Wang Original OSFD ---
    # echo "[Part 4/4] Running Wang Original OSFD with max_budget=${NSUCCESSES} (rep ${rep})..."
    # START_TIME=$(date +%s)
    # Rscript R/run_sir_WangOrig_osfd_basic.R ${NSUCCESSES} ${NSUCCESSES} ${rep} > "./logfiles/sir_part4_rep${rep}_${SLURM_JOB_ID}.Rout" 2>&1
    # END_TIME=$(date +%s)
    # RUNTIME_SEC=$((END_TIME - START_TIME))
    # RUNTIME_MIN=$(awk "BEGIN {printf \"%.2f\", ${RUNTIME_SEC}/60}")
    # echo "  Part 4 completed in ${RUNTIME_SEC} seconds (${RUNTIME_MIN} minutes)"
    # echo "${rep},OSFD_WangOrig,${RUNTIME_SEC},${RUNTIME_MIN},${NSUCCESSES}" >> ${TIMING_FILE}

    # --- Calculate grid coverage after all parts complete ---
    echo ""
    echo "All 4 parts complete for replicate ${rep}. Calculating grid coverage..."
    Rscript R/calculate_grid_coverage.R ${NSUCCESSES} ${rep} > "./logfiles/coverage_rep${rep}_${SLURM_JOB_ID}.Rout" 2>&1
    echo "Grid coverage calculated."

    echo ""
    echo "Replicate ${rep} complete (Parts 1-4 + coverage)"
    echo "-------------------------------------------"
done

echo ""
echo "==========================================="
echo "All ${NREP} replicates complete!"
echo "Timing results saved to: ${TIMING_FILE}"
echo "==========================================="

# Run plotting script for the final results
echo ""
echo "Running plotting script with NSUCCESSES=${NSUCCESSES}..."
Rscript R/plot_basic_sir_experiment.R ${NSUCCESSES} > "./logfiles/sir_plot_${SLURM_JOB_ID}.Rout" 2>&1
echo "Plotting complete."

echo ""
echo "All tasks complete!"
