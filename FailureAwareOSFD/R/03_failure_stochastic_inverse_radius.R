# ==============================================================================
# 03_failure_stochastic_inverse_radius.R
# ==============================================================================
#
# Purpose:
#   Extend the inverse-radius example to a simulator with both failures and
#   stochastic outputs (output noise).
#
#   This combines the challenges from script 02 (stochastic failures) with
#   stochastic output variation. This is NOT a replication experiment - each
#   input is evaluated once with noisy outputs.
#
#   Comparison:
#     - Blind Wang-like OSFD uses EI only and discards failures (no feasibility learning)
#     - Failure-aware OSFD uses the same single-point EI structure but adds
#       a learned success-probability adjustment (feasibility-weighted acquisition)
#
# Important:
#   This script does NOT assess batching or replication capabilities.
#   Our implementation is always called with:
#       batch_size = 1        (sequential, one point at a time)
#       n_replicates = 1      (no replication, single noisy evaluation per input)
#
# Outputs:
#   data/wang_extension/03_failure_stochastic_inverse_radius.rds
#
# Author: AC Murph
# Date: Sep 2026
# ==============================================================================

# ------------------------------------------------------------------------------
# Load required packages
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(lhs)       # Latin hypercube sampling
  library(dplyr)     # Data manipulation
  library(tidyr)     # Data tidying
  library(purrr)     # Functional programming
  library(tibble)    # Modern data frames
  library(FNN)       # Fast nearest neighbor search
  library(here)      # Path management
  library(foreach)   # Parallel foreach loops
  library(doParallel) # Parallel backend for foreach
  library(magrittr)  # Pipe operator
})

# ------------------------------------------------------------------------------
# Source our OSFD implementation and helper functions
# ------------------------------------------------------------------------------

# Load our constrained OSFD implementation and experiment setup
source(here::here("R", "sir_experiment_setup.R"))
source(here::here("R", "source_helpers.R"))
source_helpers()  # Automatically loads all helpers including wang_helper_functions.R

# ------------------------------------------------------------------------------
# Experiment settings
# ------------------------------------------------------------------------------

# Define output directory for results (using data/ instead of results/)
OUT_DIR <- here::here("data", "wang_extension")
# Create directory if it doesn't exist
dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)

# Dimensions
P <- 2  # Input dimension (2D input space)
Q <- 2  # Output dimension (2D output space)

n_cores <- 60
# Inverse-radius function parameter
EPSILON <- 0.1  # Smoothing parameter to avoid division by zero at origin

# Design parameters
N_INI <- NA                          # going with half the number of attempts rounded down
N_GRID <- c(25, 50, 75, 100, 150, 200, 300, 500, 1000)   # Attempt budgets to test (total evaluations allowed)
N_REPS <- 20                         # Number of replicate experiments per size

# Candidate set and reference grid sizes
CAND_SIZE <- 50000  # Size of candidate input set (large to account for failures during initialization)
REF_N <- 75000      # Number of stochastic reference points to generate
GRID_K <- 100       # Grid resolution for coverage metric (100x100 = 10,000 cells)

# Random seed for reproducibility
BASE_SEED <- 139834

# ------------------------------------------------------------------------------
# Failure surface parameters
# ------------------------------------------------------------------------------

# Success probability is lowest near the origin, which is also the region
# that generates large inverse-radius outputs. This stresses the original OSFD
# objective because the desired output region is expensive to sample.
#
# This failure surface uses a smooth logistic transition with no hard cutoff:
#   - Near origin: low success probability (but never exactly 0)
#   - Far from origin: high success probability
#   - Smooth transition controlled by A_FAIL steepness
R_HARD_FAIL <- 0.10    # No hard-failure core (set to 0 for smooth transition)
R_SOFT_FAIL <- 0.2   # Midpoint of soft transition where p_success ≈ 50%
A_FAIL <- 20          # Steepness of logistic transition (lower = smoother)
P_MIN_TRUE <- 0.00    # Minimum success probability (asymptote near origin)

# ------------------------------------------------------------------------------
# Stochastic output noise parameters
# ------------------------------------------------------------------------------

# Independent Gaussian noise added to each output dimension after deterministic evaluation
# Noise is added on the raw output scale (not standardized)
SIGMA_Y1 <- 0.15  # Standard deviation for first output dimension (inverse radius)
SIGMA_Y2 <- 0.03  # Standard deviation for second output dimension (angle)

# ------------------------------------------------------------------------------
# Failure-aware OSFD settings
# ------------------------------------------------------------------------------

# These parameters control the feasibility-weighted acquisition function
GAMMA_FEAS <- 2     # Feasibility weighting exponent (beta in code): higher values more strongly favor high-probability regions
P_FLOOR <- 0.1     # Minimum feasibility probability floor (prevents complete rejection of low-prob regions)

# ------------------------------------------------------------------------------
# Helper functions loaded from wang_helper_functions.R via source_helpers():
#   - inverse_radius(): Wang test function (inverse radius + angle)
#   - p_success(): Spatially-varying success probability function
#   - failure_inverse_radius(): Inverse radius with stochastic failures
#   - failure_stochastic_inverse_radius(): Inverse radius with failures AND output noise
#   - make_reference_inputs(): Generate dense regular grid
#   - eval_matrix(): Evaluate function over matrix rows
#   - eval_reference_stochastic(): Generate stochastic reference set
#   - scale_to_reference(): Scale outputs to [0,1] based on reference
#   - nearest_distances(): Compute nearest neighbor distances
#   - grid_cells(): Map outputs to discrete grid cells
#   - compute_metrics(): Comprehensive space-filling metrics
#   - safe_eval_once(): Safely evaluate function with error handling
#   - run_lhs_success_budget(): LHS until reaching success budget
# ------------------------------------------------------------------------------

# Note: The failure and stochastic parameters (P_MIN_TRUE, R0_FAIL, A_FAIL,
# SIGMA_Y1, SIGMA_Y2) defined above will be passed to the functions when needed.

# ------------------------------------------------------------------------------
# Method wrappers
# ------------------------------------------------------------------------------

#' Run LHS baseline with attempt budget
#'
#' Generates LHS inputs and evaluates until reaching max_attempts total attempts.
#' This is the naive baseline that doesn't optimize output space filling or learn
#' from failures. Outputs are stochastic (noisy) and subject to failures.
#'
#' @param max_attempts Maximum number of attempts (total evaluations)
#' @param seed Random seed for reproducibility
#' @return List with D (successful inputs), Y (successful outputs), attempts,
#'         successes, success_rate, X_all (all inputs), success (logical vector)
run_lhs_attempt_budget_wrapper <- function(max_attempts, seed) {
  # Create a wrapper function that uses our failure and stochastic parameters
  # IMPORTANT: Capture parameter values in local variables so they are
  # available when the function is evaluated in parallel workers
  epsilon_val <- EPSILON
  p_min_val <- P_MIN_TRUE
  r_hard_val <- R_HARD_FAIL
  r_soft_val <- R_SOFT_FAIL
  a_val <- A_FAIL
  sigma_y1_val <- SIGMA_Y1
  sigma_y2_val <- SIGMA_Y2

  f_with_params <- function(x) {
    failure_stochastic_inverse_radius(x, epsilon = epsilon_val,
                                     p_min_true = p_min_val,
                                     r_hard_fail = r_hard_val,
                                     r_soft_fail = r_soft_val,
                                     a_fail = a_val,
                                     sigma_y1 = sigma_y1_val,
                                     sigma_y2 = sigma_y2_val)
  }

  # Force evaluation to capture values in closure
  force(f_with_params)

  # Call the helper function from wang_helper_functions.R with attempt budget
  # Set n_success very high since we're using max_attempts to control budget
  run_lhs_success_budget(
    n_success = 1e9,  # Effectively unlimited success target
    seed = seed,
    f = f_with_params,
    p = P,
    q = Q,
    max_attempts = max_attempts
  )
}

#' Run blind Wang-like OSFD (no feasibility learning) with attempt budget
#'
#' Uses the basic Wang-style EI criterion with safe retry but does NOT learn
#' from failures. When a failure occurs, it simply retries without adjusting
#' the acquisition function. This is the baseline OSFD approach for failure-prone
#' stochastic functions.
#'
#' @param max_attempts Maximum number of attempts (total evaluations)
#' @param seed Random seed
#' @param CAND Candidate input matrix
#' @return List with D, Y, attempts, successes, success_rate, X_all, success, full
run_blind_wang_like <- function(max_attempts, seed, CAND) {
  set.seed(seed)  # Set random seed

  # Create wrapper function with failure and stochastic parameters
  # IMPORTANT: Use force() to capture parameter values in function closure
  # so they are available in parallel workers
  epsilon_val <- EPSILON
  p_min_val <- P_MIN_TRUE
  r_hard_val <- R_HARD_FAIL
  r_soft_val <- R_SOFT_FAIL
  a_val <- A_FAIL
  sigma_y1_val <- SIGMA_Y1
  sigma_y2_val <- SIGMA_Y2

  f_with_params <- function(x) {
    failure_stochastic_inverse_radius(x, epsilon = epsilon_val,
                                     p_min_true = p_min_val,
                                     r_hard_fail = r_hard_val,
                                     r_soft_fail = r_soft_val,
                                     a_fail = a_val,
                                     sigma_y1 = sigma_y1_val,
                                     sigma_y2 = sigma_y2_val)
  }

  # Force evaluation to capture values in closure
  force(f_with_params)

  # Call our implementation with Wang-like settings (NO feasibility weighting)
  # This is intentionally blind to failures - uses EI only
  fit <- constrained_osfd_ei_with_retry(
    f = f_with_params,         # Test function with failures and noise
    p = P,                     # Input dimension
    q = Q,                     # Output dimension
    CAND = CAND,               # Candidate inputs
    n_success = 1e9,           # Effectively unlimited (using max_attempts instead)
    n_ini_success = floor(max_attempts / 2),     # Initial LHS sample size

    # Sequential acquisition (one point at a time)
    batch_size = 1,
    n_replicates = 1,          # No replication (single noisy evaluation per input)
    mc.cores = 1,              # Single core (sequential)

    # Score the full remaining candidate set at each iteration
    cand_batch = nrow(CAND),

    # DISABLE feasibility weighting (blind to failures, only uses EI)
    beta = 0,                  # No feasibility multiplier
    p_floor = 1,               # No probability floor (disabled)
    tau_hard = NULL,           # No hard constraint threshold
    update_feas_every = 0,     # Never update feasibility model

    # Use attempt budget instead of success target
    max_attempts = max_attempts,

    # Diversity parameter (irrelevant when batch_size = 1)
    repel_radius = 0.05,

    verbose = FALSE            # Suppress iteration messages
  )

  # Return results with extended diagnostics
  list(
    D = fit$D_success,            # Successful inputs
    Y = fit$Y_success,            # Successful outputs
    attempts = nrow(fit$X_all),   # Total attempts (including failures)
    successes = sum(fit$success), # Number of successes
    success_rate = mean(fit$success), # Empirical success rate
    X_all = fit$X_all,            # All attempted inputs
    success = fit$success,        # Success indicator for each attempt
    full = fit                    # Full fit object for diagnostics
  )
}

#' Run failure-aware OSFD (with feasibility learning) with attempt budget
#'
#' Uses single-point EI with feasibility-weighted acquisition. The only
#' methodological difference from blind Wang-like is that this learns a GP model
#' of success probability and weights the acquisition function accordingly.
#' This allows it to adaptively avoid low-probability regions even with noisy outputs.
#'
#' @param max_attempts Maximum number of attempts (total evaluations)
#' @param seed Random seed
#' @param CAND Candidate input matrix
#' @return List with D, Y, attempts, successes, success_rate, X_all, success, full
run_failure_aware <- function(max_attempts, seed, CAND) {
  set.seed(seed)  # Set random seed

  # Create wrapper function with failure and stochastic parameters
  # IMPORTANT: Use local variables to capture parameter values in function closure
  # so they are available in parallel workers
  epsilon_val <- EPSILON
  p_min_val <- P_MIN_TRUE
  r_hard_val <- R_HARD_FAIL
  r_soft_val <- R_SOFT_FAIL
  a_val <- A_FAIL
  sigma_y1_val <- SIGMA_Y1
  sigma_y2_val <- SIGMA_Y2

  f_with_params <- function(x) {
    failure_stochastic_inverse_radius(x, epsilon = epsilon_val,
                                     p_min_true = p_min_val,
                                     r_hard_fail = r_hard_val,
                                     r_soft_fail = r_soft_val,
                                     a_fail = a_val,
                                     sigma_y1 = sigma_y1_val,
                                     sigma_y2 = sigma_y2_val)
  }

  # Force evaluation to capture values in closure
  force(f_with_params)

  # Call our implementation with feasibility-weighted acquisition
  # This learns from failures and adjusts the acquisition function
  fit <- constrained_osfd_ei_with_retry(
    f = f_with_params,         # Test function with failures and noise
    p = P,                     # Input dimension
    q = Q,                     # Output dimension
    CAND = CAND,               # Candidate inputs
    n_success = 1e9,           # Effectively unlimited (using max_attempts instead)
    n_ini_success = floor(max_attempts / 2),     # Initial LHS sample size

    # Sequential acquisition (one point at a time)
    batch_size = 1,
    n_replicates = 1,          # No replication (single noisy evaluation per input)
    mc.cores = 1,              # Single core (sequential)

    # Score the full remaining candidate set at each iteration
    cand_batch = nrow(CAND),

    # ENABLE feasibility weighting (learns from failures)
    beta = GAMMA_FEAS,         # Feasibility weighting exponent (higher = more conservative)
    p_floor = P_FLOOR,         # Minimum feasibility probability floor
    tau_hard = NULL,           # No hard constraint threshold (soft weighting only)
    update_feas_every = 1,     # Update feasibility model after each evaluation

    # Use attempt budget instead of success target
    max_attempts = max_attempts,

    # Diversity parameter (irrelevant when batch_size = 1)
    repel_radius = 0.05,

    verbose = FALSE            # Suppress iteration messages
  )

  # Return results with extended diagnostics
  list(
    D = fit$D_success,            # Successful inputs
    Y = fit$Y_success,            # Successful outputs
    attempts = nrow(fit$X_all),   # Total attempts (including failures)
    successes = sum(fit$success), # Number of successes
    success_rate = mean(fit$success), # Empirical success rate
    X_all = fit$X_all,            # All attempted inputs
    success = fit$success,        # Success indicator for each attempt
    full = fit                    # Full fit object for diagnostics
  )
}

# ------------------------------------------------------------------------------
# Run experiment
# ------------------------------------------------------------------------------

# Generate stochastic reference set for computing metrics
# Unlike scripts 01 and 02, we can't use a deterministic reference grid here
# because the outputs are stochastic (noisy). Instead, we sample a large number
# of successful stochastic evaluations to estimate the output distribution.
message("Generating stochastic reference set...")

# Create a wrapper function with the correct parameter values for reference generation
f_ref <- function(x) {
  failure_stochastic_inverse_radius(
    x,
    epsilon = EPSILON,
    p_min_true = P_MIN_TRUE,
    r_hard_fail = R_HARD_FAIL,
    r_soft_fail = R_SOFT_FAIL,
    a_fail = A_FAIL,
    sigma_y1 = SIGMA_Y1,
    sigma_y2 = SIGMA_Y2
  )
}

ref <- eval_reference_stochastic(
  n_ref = REF_N,                       # Attempt 75,000 evaluations
  seed = BASE_SEED + 999,              # Fixed seed for reference set
  f = f_ref,                           # Stochastic function with our parameters
  p = P,
  q = Q
)
X_ref <- ref$X  # Successful reference inputs
Y_ref <- ref$Y  # Successful stochastic reference outputs
message(sprintf("Generated %d successful reference outputs", nrow(Y_ref)))

# ------------------------------------------------------------------------------
# Set up parallel cluster for outer loop
# ------------------------------------------------------------------------------

# Determine number of cores to use (leave 2 cores free for system)
message(sprintf("Setting up parallel cluster with %d cores", n_cores))

# Create parallel cluster
cl <- parallel::makeCluster(n_cores, type = "PSOCK")
doParallel::registerDoParallel(cl)

# Load required packages and compile C++ on each worker FIRST
# This must happen BEFORE exporting to avoid serializing NULL C++ pointers
parallel::clusterEvalQ(cl, {
  suppressPackageStartupMessages({
    library(lhs)
    library(dplyr)
    library(tibble)
    library(FNN)
    library(here)
    library(magrittr)  # Explicitly load pipe operator
    library(Rcpp)
  })

  # Load helper functions on each worker (includes C++ compilation)
  # This creates valid C++ function pointers in each worker
  source(here::here("R", "sir_experiment_setup.R"))
  source(here::here("R", "source_helpers.R"))
  source_helpers()

  NULL
})

# Create grid of all (rep_id, n) combinations before exporting
param_grid <- expand.grid(
  rep_id = seq_len(N_REPS),
  n = N_GRID
)

# NOW export parameters, data, param_grid, AND the wrapper functions from this script
# Wrapper functions are defined in this script, not in helper files
message("Exporting parameters, data, and wrapper functions to workers...")
parallel::clusterExport(cl, varlist = c(
  "BASE_SEED", "N_INI", "N_GRID", "N_REPS", "CAND_SIZE", "P", "Q", "GRID_K",
  "EPSILON", "P_MIN_TRUE", "R_HARD_FAIL", "R_SOFT_FAIL", "A_FAIL",
  "SIGMA_Y1", "SIGMA_Y2", "GAMMA_FEAS", "P_FLOOR",
  "X_ref", "Y_ref",
  "param_grid",
  "run_lhs_attempt_budget_wrapper", "run_blind_wang_like", "run_failure_aware"
), envir = environment())

# Check if C++ functions were actually compiled in workers
message("Checking C++ function availability in workers...")
cpp_check <- parallel::clusterEvalQ(cl, {
  list(
    approx_gen_exists = exists("approx_gen"),
    approx_gen_is_function = exists("approx_gen") && is.function(approx_gen)
  )
})
message(sprintf("Worker 1 - approx_gen exists: %s, is function: %s",
                cpp_check[[1]]$approx_gen_exists,
                cpp_check[[1]]$approx_gen_is_function))

# ------------------------------------------------------------------------------
# Run parallelized experiment
# ------------------------------------------------------------------------------

message(sprintf("Running %d parameter combinations in parallel...", nrow(param_grid)))

# Prevent foreach from automatically re-exporting objects from parent environment
# Workers already have everything they need (loaded fresh with valid C++ pointers)
objects_already_available <- ls(envir = .GlobalEnv)

# Run all combinations in parallel using foreach
# Returns a list where each element contains: list(metrics = list(...), diagnostics = ...)
all_results_raw <- foreach(
  i = seq_len(nrow(param_grid)),
  .combine = "c",                    # Combine results into a single list
  .multicombine = TRUE,
  .inorder = FALSE,                  # Don't need to maintain order
  .errorhandling = "stop",           # Stop on first error
  .verbose = FALSE,
  .noexport = objects_already_available  # Don't re-export; workers already have everything
) %dopar% {
  # Extract parameters for this iteration
  rep_id <- param_grid$rep_id[i]
  n_attempts <- param_grid$n[i]  # Now this is the attempt budget, not success target

  # Generate unique seed for this (replicate, n_attempts) combination
  seed <- BASE_SEED + 30000 * rep_id + n_attempts

  # Generate candidate set for OSFD methods (same for both OSFD variants)
  # Use seed+1 to make it independent of method seeds
  set.seed(seed + 1)
  CAND <- lhs::randomLHS(CAND_SIZE, P)

  # Print progress (will be somewhat interleaved due to parallel execution)
  message(sprintf("Failure + stochastic experiment: rep %d / %d, max_attempts = %d", rep_id, N_REPS, n_attempts))

  # Run all three methods with different sub-seeds and capture timing
  # All methods get the same attempt budget for fair comparison
  # Note: Each method creates its own internal 1-core cluster (mc.cores = 1)
  # This is safe because we're already in a worker process
  fits <- list()
  timings <- list()

  # Method 1: LHS
  timing_lhs <- system.time({
    fits$LHS <- run_lhs_attempt_budget_wrapper(max_attempts = n_attempts, seed = seed + 10)
  })
  timings$LHS <- as.numeric(timing_lhs["elapsed"])

  # Method 2: Blind Wang-like OSFD
  timing_blind <- system.time({
    fits$Blind_Wang_like_OSFD <- run_blind_wang_like(max_attempts = n_attempts, seed = seed + 20, CAND = CAND)
  })
  timings$Blind_Wang_like_OSFD <- as.numeric(timing_blind["elapsed"])

  # Method 3: Failure-aware OSFD
  timing_aware <- system.time({
    fits$Failure_aware_OSFD <- run_failure_aware(max_attempts = n_attempts, seed = seed + 30, CAND = CAND)
  })
  timings$Failure_aware_OSFD <- as.numeric(timing_aware["elapsed"])

  # Create union of all successful outputs from all methods for this replicate/budget
  # This represents the total achievable output space with this attempt budget
  Y_union <- do.call(rbind, lapply(fits, function(fit) fit$Y))

  # Compute metrics for each method and collect diagnostics if needed
  method_results <- list()
  method_diagnostics <- list()

  for (method_name in names(fits)) {
    fit <- fits[[method_name]]

    # Compute output space-filling metrics against the UNION of all methods' outputs
    # This measures how well each method covers the total achievable output space
    metrics <- compute_metrics(fit$Y, Y_union, k = GRID_K)

    # Add metadata including timing
    method_results[[method_name]] <- metrics %>%
      mutate(
        experiment = "failure_stochastic_extension",     # Experiment type
        method = method_name,                 # Method name
        rep_id = rep_id,                      # Replicate ID
        max_attempts = n_attempts,            # Attempt budget used
        attempts = fit$attempts,              # Total attempts actually made
        successes = fit$successes,            # Successful evaluations
        success_rate = fit$success_rate,      # Empirical success rate
        attempts_per_success = fit$attempts / fit$successes,  # Efficiency metric
        elapsed_time_sec = timings[[method_name]]   # Wall-clock time in seconds
      )

    # Save one representative diagnostic object for plotting attempted points
    # (only for first replicate at n_attempts = 500 for visualization)
    if (rep_id == 1 && n_attempts == 500) {
      method_diagnostics[[method_name]] <- list(
        method = method_name,   # Method name for plotting
        D = fit$D,              # Successful inputs
        Y = fit$Y,              # Successful outputs
        X_all = fit$X_all,      # All attempted inputs (including failures)
        success = fit$success,  # Success indicator for each attempt
        max_attempts = n_attempts  # Store the attempt count for filtering in plotting
      )
    }
  }

  # Return both metrics and diagnostics
  # The outer list preserves one complete result per iteration
  list(list(
    metrics = method_results,
    diagnostics = method_diagnostics
  ))
}

# Stop the parallel cluster and clean up connections
parallel::stopCluster(cl)

# Force garbage collection to clean up any lingering connections
gc()

# Separate metrics from diagnostics
all_results <- list()
diagnostics <- list()

for (item in all_results_raw) {
  # Extract metrics (list of tibbles, one per method)
  all_results <- c(all_results, item$metrics)

  # Extract diagnostics (only non-empty for rep_id=1, n=max(N_GRID))
  if (length(item$diagnostics) > 0) {
    diagnostics <- c(diagnostics, item$diagnostics)
  }
}

# Combine all results into a single tibble
metrics_df <- bind_rows(all_results)

# ------------------------------------------------------------------------------
# Save results
# ------------------------------------------------------------------------------

# Save as RDS with all settings, reference data, metrics, and diagnostics
saveRDS(
  list(
    settings = list(
      epsilon = EPSILON,
      n_ini = N_INI,
      n_grid = N_GRID,
      n_reps = N_REPS,
      cand_size = CAND_SIZE,
      ref_n = REF_N,
      grid_k = GRID_K,
      base_seed = BASE_SEED,
      p_min_true = P_MIN_TRUE,
      r_hard_fail = R_HARD_FAIL,
      r_soft_fail = R_SOFT_FAIL,
      a_fail = A_FAIL,
      sigma_y1 = SIGMA_Y1,
      sigma_y2 = SIGMA_Y2,
      gamma_feas = GAMMA_FEAS,
      p_floor = P_FLOOR
    ),
    X_ref = X_ref,              # Stochastic reference input grid
    Y_ref = Y_ref,              # Stochastic reference output grid
    metrics = metrics_df,       # All computed metrics
    diagnostics = diagnostics   # Diagnostic info for visualization
  ),
  file = file.path(OUT_DIR, "03_failure_stochastic_inverse_radius.rds")
)

# Print confirmation message
message("Saved: data/wang_extension/03_failure_stochastic_inverse_radius.rds")