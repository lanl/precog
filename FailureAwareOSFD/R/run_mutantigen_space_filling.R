# Output-space filling experiment for MutAntiGen simulation
# Compares LHS on input space vs OSFD (output-space filling design)
## Author: AC Murph
## Date: May 2026

library(parallel)
library(doParallel)
library(here)
# Load helper functions
source(here::here("R", "source_helpers.R"))
source_helpers()

setwd(here::here())

# Load the wrapper function
source(here::here("R", "run_mutantigen_and_extract.R"))

# -------------------------------
# Experiment parameters
# -------------------------------
nsuccesses = 500      # Total number of successful simulations to collect
p = 10                 # Input dimension (10 MutAntiGen parameters)
q = 10                 # Output dimension (10 summary statistics)
mc.cores = 99          # Number of parallel cores (adjust based on SLURM allocation)
cand_batch = 10000     # Candidates to score per iteration in OSFD
CAND_size = 100000     # Size of candidate pool
pkg_location = this.path::here()
sim_length_days = 1095 # 3 years

# Check if running on SLURM and adjust cores if needed
if (!is.na(Sys.getenv("SLURM_CPUS_PER_TASK", unset = NA))) {
  slurm_cores <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK"))
  mc.cores <- min(slurm_cores - 1, mc.cores)  # Use n-1 cores, save 1 for overhead
  cat("Detected SLURM allocation:", slurm_cores, "cores\n")
  cat("Using:", mc.cores, "cores for parallelization\n\n")
}

# Parameter bounds - NARROWED around defaults from parameters_load.yml
# Defaults: initialNs=40M, beta=0.5627, nu=0.25, lambda=0.10, mutCost=0.008,
#          lambdaAntigenic=0.00075, meanAntigenicSize=0.012, epsilon=0.16
# NOTE: We sample R0 and mutCost directly, then calculate derived parameters:
#   nu = beta / R0 (ensures R0 >= 1.01 for transmission)
# FIXED: Previously incorrectly calculated mutCost = mutationLoad / lambda which inverted the relationship!
param_ranges <- list(
  initialNs = c(20000000, 60000000),        # ±50% around 40M default
  demeAmplitudes = c(0, 0.1),               # Narrowed from 0.2
  lambdaAntigenic = c(0.0003, 0.0015),      # ±2x around 0.00075 default
  meanAntigenicSize = c(0.006, 0.024),      # ±2x around 0.012 default
  lambda = c(0.05, 0.20),                   # ±2x around 0.10 default
  mutCost = c(0.004, 0.016),                # ±2x around 0.008 default - RENAMED from mutationLoad!
  beta = c(0.28, 1.13),                     # ±2x around 0.5627 default
  R0 = c(1.5, 3.5),                         # Narrowed from 1.01-5.0, centered around 2.25
  epsilon_mut = c(0.08, 0.32),              # ±2x around 0.16 default
  initialI_prop = c(0.001, 0.003)           # Narrowed slightly
)

# -------------------------------
# Create wrapper function that expects unit hypercube input
# -------------------------------
make_mutantigen_f <- function(param_ranges, sim_length_days, method_label) {
  force(param_ranges); force(sim_length_days); force(method_label)

  function(x) {
    # x is on unit hypercube [0,1]^10
    # Transform to actual parameter ranges

    # initialNs: log-uniform scaling (as in lhs.R line 44-47)
    initialNs_range <- param_ranges$initialNs
    log_initialNs <- log10(initialNs_range[1]) + x[1] * (log10(initialNs_range[2]) - log10(initialNs_range[1]))
    initialNs <- 10^log_initialNs

    # Other parameters: linear scaling
    demeAmplitudes <- param_ranges$demeAmplitudes[1] + x[2] * diff(param_ranges$demeAmplitudes)
    lambdaAntigenic <- param_ranges$lambdaAntigenic[1] + x[3] * diff(param_ranges$lambdaAntigenic)
    meanAntigenicSize <- param_ranges$meanAntigenicSize[1] + x[4] * diff(param_ranges$meanAntigenicSize)
    lambda <- param_ranges$lambda[1] + x[5] * diff(param_ranges$lambda)
    mutCost <- param_ranges$mutCost[1] + x[6] * diff(param_ranges$mutCost)
    beta <- param_ranges$beta[1] + x[7] * diff(param_ranges$beta)
    R0 <- param_ranges$R0[1] + x[8] * diff(param_ranges$R0)
    epsilon_mut <- param_ranges$epsilon_mut[1] + x[9] * diff(param_ranges$epsilon_mut)
    initialI_prop <- param_ranges$initialI_prop[1] + x[10] * diff(param_ranges$initialI_prop)

    # Calculate derived parameters from biological constraints
    nu <- beta / R0  # Ensures R0 >= 1.01 for sustained transmission

    # Create parameter vector
    params <- c(
      initialNs = initialNs,
      demeAmplitudes = demeAmplitudes,
      lambdaAntigenic = lambdaAntigenic,
      meanAntigenicSize = meanAntigenicSize,
      lambda = lambda,
      mutCost = mutCost,
      beta = beta,
      nu = nu,
      epsilon_mut = epsilon_mut,
      initialI_prop = initialI_prop
    )

    # Generate unique NUM for this simulation
    NUM <- as.integer((as.numeric(Sys.time()) * 1000) %% 1000000) + sample(1:10000, 1)

    # Run simulation and extract outputs
    run_mutantigen_and_extract(
      x = params,
      NUM = NUM,
      sim_length_days = sim_length_days,
      max_runtime = 3600*4,  # 4 hour timeout
      method = method_label,
      verbose = FALSE
    )
  }
}

# -------------------------------
# -------------------------------
## Method 1: Regular LHS on input space
# -------------------------------
# -------------------------------
cat("========================================\n")
cat("METHOD 1: LHS on Input Space\n")
cat("========================================\n")
cat("Generating", nsuccesses, "samples via LHS...\n\n")

X_lhs <- lhs::randomLHS(nsuccesses, p)

# Create LHS-specific wrapper
f_lhs <- make_mutantigen_f(param_ranges, sim_length_days, method_label = "lhs")

cl <- parallel::makeCluster(mc.cores)
doParallel::registerDoParallel(cl)

# Send variables to workers
parallel::clusterExport(cl, varlist = c("f_lhs", "q", "param_ranges", "sim_length_days"),
                       envir = environment())
parallel::clusterEvalQ(cl, {
  source(here::here("R", "source_helpers.R")); source_helpers()
  source(here::here("R", "run_mutantigen_and_extract.R"))
  NULL
})

cat("Running", nsuccesses, "MutAntiGen simulations in parallel...\n")
cat("This will take a while (~", ceiling(nsuccesses * 150 / mc.cores / 60), "minutes)...\n\n")

results_lhs <- foreach::foreach(
  i = seq_len(nrow(X_lhs)),
  .combine = 'rbind',
  .multicombine = TRUE,
  .inorder = TRUE,
  .verbose = TRUE
) %dopar% {
  x <- X_lhs[i, ]  # Extract as vector, not matrix

  # Transform inputs for saving
  initialNs_range <- param_ranges$initialNs
  log_initialNs <- log10(initialNs_range[1]) + x[1] * (log10(initialNs_range[2]) - log10(initialNs_range[1]))
  initialNs <- 10^log_initialNs

  demeAmplitudes <- param_ranges$demeAmplitudes[1] + x[2] * diff(param_ranges$demeAmplitudes)
  lambdaAntigenic <- param_ranges$lambdaAntigenic[1] + x[3] * diff(param_ranges$lambdaAntigenic)
  meanAntigenicSize <- param_ranges$meanAntigenicSize[1] + x[4] * diff(param_ranges$meanAntigenicSize)
  lambda <- param_ranges$lambda[1] + x[5] * diff(param_ranges$lambda)
  mutCost <- param_ranges$mutCost[1] + x[6] * diff(param_ranges$mutCost)
  beta <- param_ranges$beta[1] + x[7] * diff(param_ranges$beta)
  R0 <- param_ranges$R0[1] + x[8] * diff(param_ranges$R0)
  epsilon_mut <- param_ranges$epsilon_mut[1] + x[9] * diff(param_ranges$epsilon_mut)
  initialI_prop <- param_ranges$initialI_prop[1] + x[10] * diff(param_ranges$initialI_prop)

  # Calculate derived parameters
  nu <- beta / R0

  # Run simulation
  templist <- safe_eval_f(f_lhs, x, q)

  if (templist$success) {
    # Return inputs and outputs
    data.frame(
      # Outputs
      out1 = templist$y[1],
      out2 = templist$y[2],
      out3 = templist$y[3],
      out4 = templist$y[4],
      out5 = templist$y[5],
      out6 = templist$y[6],
      out7 = templist$y[7],
      out8 = templist$y[8],
      out9 = templist$y[9],
      out10 = templist$y[10],
      # Inputs (sampled parameters)
      initialNs = initialNs,
      demeAmplitudes = demeAmplitudes,
      lambdaAntigenic = lambdaAntigenic,
      meanAntigenicSize = meanAntigenicSize,
      lambda = lambda,
      mutCost = mutCost,            # Sampled parameter (FIXED: was mutationLoad)
      beta = beta,
      R0 = R0,                      # Sampled parameter
      epsilon_mut = epsilon_mut,
      initialI_prop = initialI_prop
    )
  } else {
    NULL  # Skip failed simulations
  }
}

parallel::stopCluster(cl)

cat("\n\nLHS completed!\n")
cat("Successful runs:", nrow(results_lhs), "/", nsuccesses, "\n")

# Save results
save(results_lhs, X_lhs,
     file = here::here("data", "mutantigen_lhs_results.RData"))
cat("Results saved to data/mutantigen_lhs_results.RData\n\n")

# -------------------------------
# -------------------------------
## Method 2: Output-Space Filling Design (OSFD)
# -------------------------------
# -------------------------------
cat("========================================\n")
cat("METHOD 2: Output-Space Filling Design\n")
cat("========================================\n")
cat("Generating", nsuccesses, "samples via OSFD...\n\n")

# Create OSFD-specific wrapper
f_osfd <- make_mutantigen_f(param_ranges, sim_length_days, method_label = "osfd")

# Generate candidate pool
CAND <- lhs::randomLHS(CAND_size, p)

cat("Running OSFD with:\n")
cat("  Target successes:", nsuccesses, "\n")
cat("  Initial successes:", floor(2 * nsuccesses / 4), "\n")
cat("  Candidate pool:", CAND_size, "\n")
cat("  Batch size:", mc.cores, "\n\n")

# Start logging OSFD output to file in logfiles/
logfiles_dir <- here::here("logfiles")
if (!dir.exists(logfiles_dir)) {
  dir.create(logfiles_dir, recursive = TRUE)
}
osfd_log <- file.path(logfiles_dir, paste0("osfd_algorithm_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".log"))
cat("Logging OSFD algorithm output to:", osfd_log, "\n\n")
sink(osfd_log, split = TRUE)  # split=TRUE means output goes to both file and console

res <- constrained_osfd_ei_with_retry(
  f_osfd, p, q,
  CAND,
  nsuccesses,
  n_ini_success = floor(2 * nsuccesses / 4),
  batch_size = mc.cores,
  n_replicates = 1,           # No replicates for MutAntiGen (too expensive)
  cand_batch = cand_batch,
  beta = 2.0,
  p_floor = 0.05,
  tau_hard = NULL,
  repel_radius = 0.05,        # Input-space diversity scale
  update_feas_every = 1,
  mc.cores = mc.cores,
  verbose = TRUE,
  pkg_location = pkg_location
)

cat("\n\nOSFD completed!\n")

# Extract results
X_osfd <- res$D_success
Y_osfd <- res$Y_success

cat("Successful runs:", nrow(X_osfd), "\n")

# Save results
save(X_osfd, Y_osfd,
     file = here::here("data", "mutantigen_osfd_results.RData"))
cat("Results saved to data/mutantigen_osfd_results.RData\n\n")

# -------------------------------
# Summary
# -------------------------------
cat("========================================\n")
cat("EXPERIMENT COMPLETE\n")
cat("========================================\n")
cat("LHS results:", nrow(results_lhs), "simulations\n")
cat("OSFD results:", nrow(X_osfd), "simulations\n")
cat("\nData files saved to data/\n")
cat("Timing logs saved to mutantigen_parallel/timing_logs/\n")
