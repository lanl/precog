# ==============================================================================
# 01_reproduce_wang_inverse_radius.R
# ==============================================================================
#
# Purpose:
#   Reproduce the canonical Wang et al. (2024) inverse-radius OSFD example
#   and compare three methods:
#     1. Input-space LHS (baseline random sampling)
#     2. CRAN OSFD::OSFD() (Wang's published implementation)
#     3. Our constrained_osfd_ei_with_retry() run in Wang-like mode
#
# Important:
#   This script does NOT assess batching or replication capabilities.
#   Our implementation is called with:
#       batch_size = 1        (sequential, one point at a time)
#       n_replicates = 1      (deterministic, no replication)
#
#   In the Wang-like mode, feasibility modeling is disabled:
#       beta = 0              (no feasibility weighting)
#       p_floor = 1           (no probability floor)
#       update_feas_every = 0 (never update feasibility model)
#
# Outputs:
#   data/wang_extension/01_reproduce_wang_inverse_radius.rds
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
  library(OSFD)      # Wang's CRAN package
  library(here)      # Path management
  library(foreach)   # Parallel foreach loops
  library(doParallel) # Parallel backend for foreach
  library(magrittr)  # Pipe operator
})

# ------------------------------------------------------------------------------
# Source our OSFD implementation and helper functions
# ------------------------------------------------------------------------------

# Load our constrained OSFD implementation
source(here::here("R", "sir_experiment_setup.R"))
source(here::here("R", "source_helpers.R"))
source_helpers()

# Load shared helper functions for Wang experiments
# source(here::here("R", "helpers", "wang_helper_functions.R"))

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
N_INI <- 10                          # Number of initial LHS points for OSFD
N_GRID <- c(25, 50, 75, 100, 150, 200, 300, 500, 1000)    # Target successful output sizes to test
N_REPS <- 20                         # Number of replicate experiments per size

# Candidate set and reference grid sizes
CAND_SIZE <- 5000  # Size of candidate input set for OSFD algorithms
REF_SIDE <- 300    # Resolution of reference grid (300x300 = 90,000 points)
GRID_K <- 100      # Grid resolution for coverage metric (100x100 = 10,000 cells)

# Random seed for reproducibility
BASE_SEED <- 13

# ------------------------------------------------------------------------------
# Helper functions loaded from wang_helper_functions.R:
#   - inverse_radius(): Wang test function (inverse radius + angle)
#   - make_reference_inputs(): Generate dense regular grid
#   - eval_matrix(): Evaluate function over matrix rows
#   - scale_to_reference(): Scale outputs to [0,1] based on reference
#   - nearest_distances(): Compute nearest neighbor distances
#   - grid_cells(): Map outputs to discrete grid cells
#   - compute_metrics(): Comprehensive space-filling metrics
# ------------------------------------------------------------------------------

# ------------------------------------------------------------------------------
# Method wrappers
# ------------------------------------------------------------------------------

#' Run input-space LHS baseline
#'
#' Generates n random LHS inputs and evaluates the deterministic function.
#' This is the naive baseline that doesn't optimize output space filling.
#'
#' @param n Number of design points to generate
#' @param seed Random seed for reproducibility
#' @return List with D (inputs), Y (outputs), attempts, successes
run_lhs <- function(n, seed) {
  set.seed(seed)  # Set random seed

  # Generate n-point LHS in input space
  D <- lhs::randomLHS(n, P)

  # Evaluate inverse_radius for each input (deterministic, no failures)
  Y <- eval_matrix(D, inverse_radius)

  # Return design and metadata (for deterministic case, attempts = successes = n)
  list(
    D = D,           # Input design matrix (n x P)
    Y = Y,           # Output matrix (n x Q)
    attempts = n,    # Number of function evaluations
    successes = n    # Number of successful evaluations
  )
}

#' Run CRAN OSFD::OSFD() from Wang et al. (2024)
#'
#' Uses the published Wang OSFD implementation with Expected Improvement.
#' This is the reference method we are comparing against.
#'
#' @param n Target number of outputs
#' @param seed Random seed
#' @param CAND Candidate input matrix for acquisition function
#' @return List with D (inputs), Y (outputs), attempts, successes
run_cran_osfd <- function(n, seed, CAND) {
  set.seed(seed)  # Set random seed

  # Call CRAN OSFD::OSFD with Expected Improvement acquisition
  # The function directly evaluates the deterministic test function
  fit <- OSFD::OSFD(
    f = inverse_radius,   # Test function
    p = P,                # Input dimension
    q = Q,                # Output dimension
    n_ini = N_INI,        # Initial LHS sample size
    n = n,                # Target output size
    scale = TRUE,         # Standardize outputs for GP modeling
    method = "EI",        # Expected Improvement acquisition function
    CAND = CAND,          # Candidate inputs for optimization
    rand_out = FALSE,  # Use quasi-random/twinning points, not random points, for output-space approximation
    rand_in  = FALSE   # Use quasi-random/twinning points, not random points, for input-space candidate generation
  )

  # Return results in standard format
  list(
    D = fit$D,                # Final input design (n x P)
    Y = fit$Y,                # Final outputs (n x Q)
    attempts = nrow(fit$D),   # Number of function evaluations
    successes = nrow(fit$Y)   # Number of successful evaluations
  )
}

#' Run our implementation in Wang-like mode
#'
#' Calls constrained_osfd_ei_with_retry with settings that match Wang's
#' approach: sequential (batch_size=1), deterministic (n_replicates=1),
#' and no feasibility modeling (beta=0).
#'
#' @param n Target number of successful outputs
#' @param seed Random seed
#' @param CAND Candidate input matrix
#' @return List with D, Y, attempts, successes, and full fit object
run_ours_wang_like <- function(n, seed, CAND) {
  set.seed(seed)  # Set random seed

  # Call our implementation with Wang-like settings
  # This is intentionally single-point, non-replicated, and feasibility-free.
  fit <- constrained_osfd_ei_with_retry(
    f = inverse_radius,        # Test function
    p = P,                     # Input dimension
    q = Q,                     # Output dimension
    CAND = CAND,               # Candidate inputs
    n_success = n,             # Target number of successful outputs
    n_ini_success = N_INI,     # Initial LHS sample size

    # Sequential acquisition (one point at a time, like Wang)
    batch_size = 1,
    n_replicates = 1,          # No replication (deterministic)
    mc.cores = 1,              # Single core (sequential)

    # Score the full remaining candidate set at each iteration
    cand_batch = nrow(CAND),

    # Disable feasibility weighting (this is a deterministic function, no failures)
    beta = 0,                  # No feasibility multiplier
    p_floor = 1,               # No probability floor (disabled)
    tau_hard = NULL,           # No hard constraint threshold
    update_feas_every = 0,     # Never update feasibility model

    # Diversity parameter (irrelevant when batch_size = 1, but kept for completeness)
    repel_radius = 0.05,

    verbose = FALSE            # Suppress iteration messages
  )

  # Return results in standard format
  list(
    D = fit$D_success,          # Successful inputs (n x P)
    Y = fit$Y_success,          # Successful outputs (n x Q)
    attempts = nrow(fit$X_all), # Total function evaluations attempted
    successes = sum(fit$success), # Number of successful evaluations
    full = fit                  # Save full fit object for diagnostics
  )
}

# ------------------------------------------------------------------------------
# Run experiment
# ------------------------------------------------------------------------------

# Generate dense reference grid for computing metrics
# This gives us 90,000 output points covering the output space
X_ref <- make_reference_inputs(REF_SIDE)  # 300x300 grid in [0,1]^2
Y_ref <- eval_matrix(X_ref, inverse_radius)  # Evaluate inverse_radius on grid

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
    library(OSFD)
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

# NOW export parameters, data, param_grid, and wrapper functions
# Do NOT export C++ functions - they're already compiled on workers
message("Exporting parameters, data, and wrapper functions to workers...")
parallel::clusterExport(cl, varlist = c(
  "BASE_SEED", "N_INI", "N_GRID", "N_REPS", "CAND_SIZE", "P", "Q", "GRID_K",
  "X_ref", "Y_ref", "param_grid",
  "run_lhs", "run_cran_osfd", "run_ours_wang_like"
), envir = environment())

# ------------------------------------------------------------------------------
# Run parallelized experiment
# ------------------------------------------------------------------------------

message(sprintf("Running %d parameter combinations in parallel...", nrow(param_grid)))

# Prevent foreach from automatically re-exporting objects from parent environment
# Workers already have everything they need (loaded fresh with valid C++ pointers)
objects_already_available <- ls(envir = .GlobalEnv)

# Run all combinations in parallel using foreach
all_results <- foreach(
  i = seq_len(nrow(param_grid)),
  .combine = "c",                    # Combine results into a single list
  .multicombine = TRUE,
  .inorder = FALSE,                  # Don't need to maintain order
  .errorhandling = "stop",           # Stop on first error
  .verbose = FALSE,
  .noexport = objects_already_available  # Don't re-export; workers already have everything
) %dopar% {
  # Create grid of all (rep_id, n) combinations before exporting
  param_grid <- expand.grid(
    rep_id = seq_len(N_REPS),
    n = N_GRID
  )
  # Extract parameters for this iteration
  rep_id <- param_grid$rep_id[i]
  n <- param_grid$n[i]

  # Generate unique seed for this (replicate, n) combination
  seed <- BASE_SEED + 10000 * rep_id + n

  # Generate candidate set for OSFD methods (same for CRAN and ours)
  # Use seed+1 to make it independent of method seeds
  set.seed(seed + 1)
  CAND <- lhs::randomLHS(CAND_SIZE, P)

  # Print progress (will be somewhat interleaved due to parallel execution)
  message(sprintf("Reproduction: rep %d / %d, n = %d", rep_id, N_REPS, n))

  # Run all three methods with different sub-seeds and capture timing
  # Note: Each method creates its own internal 1-core cluster (mc.cores = 1)
  # This is safe because we're already in a worker process
  fits <- list()
  timings <- list()

  # Method 1: LHS
  timing_lhs <- system.time({
    fits$LHS <- run_lhs(n = n, seed = seed + 10)
  })
  timings$LHS <- as.numeric(timing_lhs["elapsed"])

  # Method 2: CRAN OSFD
  timing_cran <- system.time({
    fits$CRAN_OSFD <- run_cran_osfd(n = n, seed = seed + 20, CAND = CAND)
  })
  timings$CRAN_OSFD <- as.numeric(timing_cran["elapsed"])

  # Method 3: Ours Wang-like
  timing_ours <- system.time({
    fits$Ours_Wang_like <- run_ours_wang_like(n = n, seed = seed + 30, CAND = CAND)
  })
  timings$Ours_Wang_like <- as.numeric(timing_ours["elapsed"])

  # Compute metrics for each method and return as list
  method_results <- list()
  for (method_name in names(fits)) {
    fit <- fits[[method_name]]

    # Compute output space-filling metrics against reference grid
    metrics <- compute_metrics(fit$Y, Y_ref, k = GRID_K)

    # Add metadata including timing
    method_results[[method_name]] <- metrics %>%
      mutate(
        experiment = "deterministic_reproduction",  # Experiment type
        method = method_name,                       # Method name
        rep_id = rep_id,                            # Replicate ID
        n_target = n,                               # Target output size
        attempts = fit$attempts,                    # Function evaluations
        successes = fit$successes,                  # Successful evaluations
        elapsed_time_sec = timings[[method_name]]   # Wall-clock time in seconds
      )
  }

  # Return list of results for this (rep_id, n) combination
  method_results
}

# Stop the parallel cluster
parallel::stopCluster(cl)

# Combine all results into a single tibble
# all_results is now a list of lists (one per method per combination)
metrics_df <- bind_rows(all_results)

# ------------------------------------------------------------------------------
# Save results
# ------------------------------------------------------------------------------

# Save as RDS with all settings, reference data, and metrics
saveRDS(
  list(
    settings = list(
      epsilon = EPSILON,
      n_ini = N_INI,
      n_grid = N_GRID,
      n_reps = N_REPS,
      cand_size = CAND_SIZE,
      ref_side = REF_SIDE,
      grid_k = GRID_K,
      base_seed = BASE_SEED
    ),
    X_ref = X_ref,        # Reference input grid
    Y_ref = Y_ref,        # Reference output grid
    metrics = metrics_df  # All computed metrics
  ),
  file = file.path(OUT_DIR, "01_reproduce_wang_inverse_radius.rds")
)

# Print confirmation message
message("Saved: data/wang_extension/01_reproduce_wang_inverse_radius.rds")