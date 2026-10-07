#!/bin/bash
#SBATCH --job-name=stoSIR_all
#SBATCH -p ccs6
#SBATCH --qos=work
#SBATCH -c 100
#SBATCH -t 20:00:00
#SBATCH --array=1-4
#SBATCH --mem=150G
#SBATCH --output=SLURMOUT/stochastic_sir_all_%A-%a.out

# Top-level script to run all four stochastic SIR experiments
# Array jobs: 1=LHS no rep, 2=LHS w/rep, 3=OSFD w/rep, 4=OSFD no rep

# Default nsuccesses (can be overridden by setting NSUCCESSES environment variable)
NSUCCESSES=${NSUCCESSES:-5000}

echo "========================================="
echo "Running stochastic SIR experiment array"
echo "Job ID: ${SLURM_ARRAY_JOB_ID}"
echo "Task ID: ${SLURM_ARRAY_TASK_ID}"
echo "nsuccesses: ${NSUCCESSES}"
echo "========================================="

# --- Adaptive cleanup on first task only ---
if [ "${SLURM_ARRAY_TASK_ID}" -eq 1 ]; then
    echo "Task 1: Cleaning up old job outputs (keeping current job: ${SLURM_ARRAY_JOB_ID})"

    # Remove old SLURMOUT files (not from current job)
    if [ -d "SLURMOUT" ] && [ -n "$SLURM_ARRAY_JOB_ID" ]; then
        find SLURMOUT -type f -name "stochastic_sir_all_*.out" ! -name "stochastic_sir_all_${SLURM_ARRAY_JOB_ID}-*.out" -delete
        echo "Removed old SLURMOUT files"
    fi

    # Remove old logfiles
    if [ -d "logfiles" ]; then
        find logfiles -type f -name "stochastic_sir_*.Rout" -delete
        echo "Removed old logfiles"
    fi

    # Remove intermediate data files
    echo "Removing old data files..."
    rm -f \
        data/inputs_and_outputs_sirStochasticLHS_noReplicates.RData \
        data/inputs_and_outputs_sirStochasticLHS_wReplicates.RData \
        data/inputs_and_outputs_StochasticSirOutputFilling_wReplicates.RData \
        data/inputs_and_outputs_StochasticSirOutputFilling_noReplicates.RData
    echo "Cleanup complete"
fi

# --- Set up R environment ---
if [ -d "$HOME/R-4.3.2/lib64/R" ] && [ -x "$HOME/R-4.3.2/lib64/R/bin/R" ]; then
    echo "Using custom R installation"
    module purge
    export R_HOME="$HOME/R-4.3.2/lib64/R"
    export PATH="$HOME/R-4.3.2/lib64/R/bin:$PATH"
    export LD_LIBRARY_PATH="$R_HOME/lib:$LD_LIBRARY_PATH"
    export R_LIBS_USER="$HOME/R/x86_64-pc-linux-gnu-library/4.3"
else
    echo "Custom R not found, using module load R"
    module load R
fi

# --- Run the appropriate script based on array task ID ---
case ${SLURM_ARRAY_TASK_ID} in
    1)
        echo "Running: LHS without replicates"
        SCRIPT="R/run_stochasticSIR_LHS_noReplicates.R"
        LOGFILE="logfiles/stochastic_sir_LHS_noRep.Rout"
        ;;
    2)
        echo "Running: LHS with replicates"
        SCRIPT="R/run_stochasticSIR_LHS_wReplicates.R"
        LOGFILE="logfiles/stochastic_sir_LHS_wRep.Rout"
        ;;
    3)
        echo "Running: OSFD with replicates"
        SCRIPT="R/run_stochasticSIR_OSFD_wReplicates.R"
        LOGFILE="logfiles/stochastic_sir_OSFD_wRep.Rout"
        ;;
    4)
        echo "Running: OSFD without replicates"
        SCRIPT="R/run_stochasticSIR_OSFD_noReplicates.R"
        LOGFILE="logfiles/stochastic_sir_OSFD_noRep.Rout"
        ;;
    *)
        echo "ERROR: Invalid SLURM_ARRAY_TASK_ID: ${SLURM_ARRAY_TASK_ID}"
        exit 1
        ;;
esac

echo "Script: ${SCRIPT}"
echo "Log file: ${LOGFILE}"
echo "Starting at: $(date)"
echo "========================================="

# Run the R script with nsuccesses as argument
Rscript ${SCRIPT} ${NSUCCESSES} > ${LOGFILE} 2>&1

EXIT_CODE=$?
echo "========================================="
echo "Finished at: $(date)"
echo "Exit code: ${EXIT_CODE}"
echo "========================================="

exit ${EXIT_CODE}
