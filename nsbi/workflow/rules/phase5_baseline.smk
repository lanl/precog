"""
Phase 5: Baseline Methods for Comparison

Fit standard SIR model via maximum likelihood as a baseline to compare against NSBI.
Uses the same test simulations and daily binned case counts from Phase 2.

Uses r0rt parameterization (fits R0 and recovery_time directly) with prior mean starting values.

Also includes Bayesian baseline using Stan/CmdStanR with MCMC inference.
"""

import os
import time

# ============================================================================
# SIR MLE Baseline - r0rt parameterization with prior mean starting values
# ============================================================================

rule sir_mle_batch:
    """Fit SIR model via MLE to a batch of test simulations."""
    input:
        params="results/phase1_simulation/parameters/test_design.csv",
        simulations=expand(
            "results/batch_markers/beast_batch_{batch_id}.done",
            batch_id=range(N_BEAST_BATCHES)
        )
    output:
        temp("results/phase5_baseline/sir_mle/batch_{batch_id}.csv")
    params:
        sim_ids=lambda wildcards: BASELINE_BATCHES[int(wildcards.batch_id)][1],
        starting_value_strategy="r0rt_prior_mean",  # Hardcoded strategy
        R0_min=config["R0_min"],
        R0_max=config["R0_max"],
        recovery_time_min=config["recovery_time_min"],
        recovery_time_max=config["recovery_time_max"],
        S_fixed=config["S_fixed"]
    log:
        "logs/sir_mle_batch/batch_{batch_id}.log"
    threads: 1
    conda:
        "../envs/r-baseline.yaml"
    script:
        "../scripts/phase5_baseline/01_sir_mle_batch_modular.R"

rule aggregate_sir_mle:
    """Aggregate all SIR MLE batch results."""
    input:
        batch_files=expand(
            "results/phase5_baseline/sir_mle/batch_{batch_id}.csv",
            batch_id=range(N_BASELINE_BATCHES)
        ),
        params="results/phase1_simulation/parameters/test_design.csv"
    output:
        results="results/phase5_baseline/sir_mle/results.csv",
        summary="results/phase5_baseline/sir_mle/summary.yaml"
    params:
        expected_sims=N_TEST_SIMS
    log:
        "logs/aggregate_sir_mle.log"
    conda:
        "../envs/r-baseline.yaml"
    script:
        "../scripts/phase5_baseline/02_aggregate_sir_mle.R"

# ============================================================================
# SIR Bayesian Baseline - Stan/CmdStanR with MCMC
# ============================================================================

rule sir_bayes_batch:
    """Fit SIR model via Bayesian MCMC to a batch of test simulations."""
    input:
        params="results/phase1_simulation/parameters/test_design.csv",
        stan_model="workflow/scripts/phase5_baseline/sir_bayes.stan",
        simulations=expand(
            "results/batch_markers/beast_batch_{batch_id}.done",
            batch_id=range(N_BEAST_BATCHES)
        )
    output:
        temp("results/phase5_baseline/sir_bayes/batch_{batch_id}.csv")
    params:
        sim_ids=lambda wildcards: BAYES_BATCHES[int(wildcards.batch_id)][1],
        sim_ids_r=lambda wildcards: ', '.join(f'"{sid}"' for sid in BAYES_BATCHES[int(wildcards.batch_id)][1]),
        R0_min=config["R0_min"],
        R0_max=config["R0_max"],
        recovery_time_min=config["recovery_time_min"],
        recovery_time_max=config["recovery_time_max"],
        S_fixed=config["S_fixed"],
        prior_buffer=config["bayesian_baseline"]["prior_buffer_fraction"],
        n_chains=config["bayesian_baseline"]["mcmc"]["n_chains"],
        n_warmup=config["bayesian_baseline"]["mcmc"]["n_warmup"],
        n_iter=config["bayesian_baseline"]["mcmc"]["n_iter"],
        adapt_delta=config["bayesian_baseline"]["mcmc"]["adapt_delta"],
        max_treedepth=config["bayesian_baseline"]["mcmc"]["max_treedepth"]
    log:
        "logs/sir_bayes_batch/batch_{batch_id}.log"
    threads: 3  # Number of parallel chains
    conda:
        "../envs/r-baseline.yaml"  # Use baseline R environment (now includes cmdstanr)
    shell:
        """
        # Force conda R to be used (macOS PATH workaround)
        export R_HOME="$CONDA_PREFIX/lib/R"
        export R_LIBS="$CONDA_PREFIX/lib/R/library"
        export R_LIBS_USER="$CONDA_PREFIX/lib/R/library"
        export PATH="$CONDA_PREFIX/bin:$PATH"
        
        # Stan/CmdStanR cluster compatibility fixes
        export STAN_NUM_THREADS={threads}  # Match allocated threads
        export TBB_CXX_TYPE=gcc  # Use GCC-compatible TBB (cluster compatibility)
        export CMDSTAN_OUTPUT_DIR=/tmp/cmdstan_output_$$  # Use local temp (avoid NFS issues)
        mkdir -p $CMDSTAN_OUTPUT_DIR
        
        # Run R script
        "$CONDA_PREFIX/bin/Rscript" --vanilla workflow/scripts/phase5_baseline/03_sir_bayes_batch_shell.R \
            {wildcards.batch_id} \
            {input.params} \
            {input.stan_model} \
            {output[0]} \
            {log[0]} \
            'c({params.sim_ids_r})' \
            {params.R0_min} {params.R0_max} \
            {params.recovery_time_min} {params.recovery_time_max} \
            {params.S_fixed} {params.prior_buffer} \
            {params.n_chains} {params.n_warmup} {params.n_iter} \
            {params.adapt_delta} {params.max_treedepth}
        
        # Cleanup temp directory
        rm -rf $CMDSTAN_OUTPUT_DIR
        """

rule aggregate_sir_bayes:
    """Aggregate all SIR Bayesian batch results."""
    input:
        batch_files=expand(
            "results/phase5_baseline/sir_bayes/batch_{batch_id}.csv",
            batch_id=range(N_BAYES_BATCHES)
        ),
        params="results/phase1_simulation/parameters/test_design.csv"
    output:
        results="results/phase5_baseline/sir_bayes/results.csv",
        summary="results/phase5_baseline/sir_bayes/summary.yaml"
    params:
        expected_sims=len(BAYES_TEST_IDS)
    log:
        "logs/aggregate_sir_bayes.log"
    conda:
        "../envs/r-baseline.yaml"
    script:
        "../scripts/phase5_baseline/04_aggregate_sir_bayes.R"

# ============================================================================
# MLE vs Bayes Comparison Visualization
# ============================================================================

rule identify_comparison_candidates:
    """Identify simulations that converged in both MLE and Bayes, select subset for viz."""
    input:
        mle="results/phase5_baseline/sir_mle/results.csv",
        bayes="results/phase5_baseline/sir_bayes/results.csv"
    output:
        converged="results/phase5_baseline/comparison/converged_sims.txt",
        selected="results/phase5_baseline/comparison/selected_sims.txt"
    params:
        n_examples=lambda wildcards: config.get("bayesian_baseline", {}).get("comparison", {}).get("n_examples", 3),
        seed=lambda wildcards: config.get("bayesian_baseline", {}).get("comparison", {}).get("seed") if config.get("bayesian_baseline", {}).get("comparison", {}).get("seed") is not None else int(time.time() * 1000) % 2**31
    log:
        "logs/identify_comparison_candidates.log"
    conda:
        "../envs/python-core.yaml"
    script:
        "../scripts/phase5_baseline/identify_converged_for_snakemake.py"

rule refit_comparison_simulations:
    """Re-fit selected simulations with full covariance matrices and posterior draws."""
    input:
        selected="results/phase5_baseline/comparison/selected_sims.txt",
        params="results/phase1_simulation/parameters/test_design.csv",
        stan_model="workflow/scripts/phase5_baseline/sir_bayes.stan"
    output:
        done="results/phase5_baseline/comparison/refit.done"
    params:
        sim_ids=lambda wildcards, input: open(input.selected).read().strip().split() if os.path.getsize(input.selected) > 0 else []
    log:
        "logs/refit_comparison_simulations.log"
    threads: 4
    conda:
        "../envs/r-baseline.yaml"
    shell:
        """
        # Check if we have simulations to refit
        if [ -s {input.selected} ]; then
            # Force conda R to be used
            export R_HOME="$CONDA_PREFIX/lib/R"
            export R_LIBS="$CONDA_PREFIX/lib/R/library"
            export R_LIBS_USER="$CONDA_PREFIX/lib/R/library"
            export PATH="$CONDA_PREFIX/bin:$PATH"
            
            # Stan/CmdStanR settings
            export STAN_NUM_THREADS={threads}
            export TBB_CXX_TYPE=gcc
            export CMDSTAN_OUTPUT_DIR=/tmp/cmdstan_comparison_$$
            mkdir -p $CMDSTAN_OUTPUT_DIR
            
            # Run refitting
            Rscript workflow/scripts/phase5_baseline/refit_for_comparison.R {params.sim_ids} > {log} 2>&1
            
            # Cleanup
            rm -rf $CMDSTAN_OUTPUT_DIR
        else
            echo "No converged simulations found, skipping refit" > {log}
        fi
        
        touch {output.done}
        """

rule plot_comparison:
    """Generate MLE vs Bayes comparison visualization."""
    input:
        refit_done="results/phase5_baseline/comparison/refit.done",
        selected="results/phase5_baseline/comparison/selected_sims.txt"
    output:
        plot="reports/mle_vs_bayes_comparison.pdf"
    log:
        "logs/plot_comparison.log"
    conda:
        "../envs/python-core.yaml"
    script:
        "../scripts/phase5_baseline/plot_mle_bayes_comparison.py"



