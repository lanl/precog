#!/bin/bash
#SBATCH --job-name=mutantigen_space_filling
#SBATCH -p ccs6
#SBATCH --qos=work
#SBATCH -c 100
#SBATCH -t 03-00:00:00
#SBATCH --array=1-1
#SBATCH --output=SLURMOUT/mutantigen_space_filling_%j.out
#SBATCH --error=SLURMOUT/mutantigen_space_filling_%j.err

# Print job information
echo "Job started at: $(date)"
echo "Job ID: $SLURM_JOB_ID"
echo "Running on node: $(hostname)"
echo "CPUs allocated: $SLURM_CPUS_PER_TASK"
echo ""
# Set working directory to the repository root
# Change to script directory
cd "$(dirname "$(readlink -f "${BASH_SOURCE[0]}")")"

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

# Clean up previous run files
echo "Cleaning up previous run files..."
rm -f mutantigen_parallel/input_files/parameters_load_temp_*.yml
rm -f mutantigen_parallel/outfiles/out_*.*
echo "Cleanup complete"
echo ""

# Clean up old SLURM output files (keep only current job)
echo "Cleaning up old SLURMOUT files (keeping current job $SLURM_JOB_ID)..."
find SLURMOUT -type f \( -name "*.out" -o -name "*.err" \) ! -name "*_${SLURM_JOB_ID}.*" -delete 2>/dev/null || true
echo "SLURMOUT cleanup complete"
echo ""

# Clean up old logfiles (keep only files from current run, which will be created during execution)
echo "Cleaning up old logfiles..."
rm -f logfiles/java_*.log 2>/dev/null || true
rm -f logfiles/osfd_algorithm_*.log 2>/dev/null || true
rm -f mutantigen_parallel/logfiles/*.log 2>/dev/null || true
echo "Logfiles cleanup complete"
echo ""

# Run the experiment script
Rscript R/run_mutantigen_space_filling.R

# Print completion time
echo ""
echo "Job completed at: $(date)"
