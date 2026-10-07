# Part 4: Wang Original OSFD Settings
# Runs constrained_osfd_ei_with_retry with original Wang 2024 paper settings
# Author: AC Murph
# Date: Sep 2026

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript run_sir_WangOrig_osfd_basic.R <nsuccesses> [budget] [rep]")
}
nsuccesses <- as.numeric(args[1])

# Optional: computational budget (max_attempts)
budget <- NULL
if (length(args) >= 2) {
  budget <- as.numeric(args[2])
  cat(sprintf("Running Part 4 (Wang Original Settings) with computational budget = %d\n", budget))
} else {
  cat(sprintf("Running Part 4 (Wang Original Settings) with nsuccesses = %d\n", nsuccesses))
}

# Optional: replicate number
rep <- 1
if (length(args) >= 3) {
  rep <- as.numeric(args[3])
  cat(sprintf("Replicate number: %d\n", rep))
}

# Load shared setup
source(here::here("R", "sir_experiment_setup.R"))

# -------------------------------
# Wang Original OSFD Settings
# -------------------------------

cat("\n=== PART 4: Wang Original OSFD Settings ===\n")
cat("Using Wang 2024 paper original settings:\n")
cat("  - batch_size = 1 (sequential acquisition)\n")
cat("  - n_replicates = 1 (deterministic)\n")
cat("  - cand_batch = full candidate set\n")
cat("  - beta = 0 (no feasibility modification)\n")
cat("  - p_floor = 1 (no feasibility floor)\n")
cat("  - tau_hard = NULL (no hard threshold)\n")
cat("  - update_feas_every = 0 (never update feasibility)\n\n")

# Set parameters based on whether budget is provided
if (!is.null(budget)) {
  # Budget mode: run until max_attempts reached
  n_target <- 1e10  # Very large number to effectively disable success-based stopping
  n_ini_success <- floor(budget / 2)
  max_attempts <- budget
  cat(sprintf("Budget mode: max_attempts = %d, n_ini_success = %d\n",
              max_attempts, n_ini_success))
} else {
  # Success mode: run until nsuccesses reached
  n_target <- nsuccesses
  n_ini_success <- floor(2 * nsuccesses / 4)
  max_attempts <- NULL
  cat(sprintf("Success mode: nsuccesses = %d, n_ini_success = %d\n",
              n_target, n_ini_success))
}

# LHS on output space for candidate set
CAND <- lhs::randomLHS(CAND_size, p)

# Define the forward mapping function
make_f <- function(alpha_bounds, reproduction_number_bounds, s0_bounds) {
  force(alpha_bounds); force(reproduction_number_bounds); force(s0_bounds)
  function(x) {
    calculate_sir_for_outputSpaceFilling(
      x[1], x[2], x[3],
      alpha_bounds = alpha_bounds,
      reproduction_number_bounds = reproduction_number_bounds,
      s0_bounds = s0_bounds
    )
  }
}
f <- make_f(alpha_bounds, reproduction_number_bounds, s0_bounds)

print("About to enter OSFD function with Wang original settings.")
res <- constrained_osfd_ei_with_retry(
  f = f,
  p = p,
  q = q,
  CAND = CAND,
  n_success = n_target,
  n_ini_success = n_ini_success,

  # Wang-like settings
  batch_size = min(mc.cores, n_ini_success),
  n_replicates = 1,
  cand_batch = nrow(CAND),

  # Turn off feasibility modification
  beta = 0,
  p_floor = 1,
  tau_hard = NULL,
  update_feas_every = 0,

  # Irrelevant when batch_size = 1
  repel_radius = 0.05,

  max_attempts = max_attempts,  # NULL or budget value
  mc.cores = mc.cores,
  verbose = TRUE,
  pkg_location = pkg_location
)
print("We have now exited the OSFD function.")

D <- res$D_success
Y <- res$Y_success

# Determine actual number of successes obtained
actual_n <- nrow(D)

# Save results with appropriate suffix (unique naming for Wang Original)
if (!is.null(budget)) {
  # Budget mode: use budget in filename with replicate number
  save(D, file = here::here("data", sprintf("inputs_sirOutputFilling_WangOrig_budget%d_n%d_rep%d.RData", budget, actual_n, rep)))
  save(Y, file = here::here("data", sprintf("outputs_sirOutputFilling_WangOrig_budget%d_n%d_rep%d.RData", budget, actual_n, rep)))
  save(res, file = here::here("data", sprintf("fullOutput_OSFD_WangOrig_budget%d_n%d_rep%d.RData", budget, actual_n, rep)))
  cat(sprintf("Results saved with WangOrig_budget=%d_rep%d suffix (obtained %d successes)\n", budget, rep, actual_n))
} else {
  # Success mode: use nsuccesses in filename with replicate number
  save(D, file = here::here("data", sprintf("inputs_sirOutputFilling_WangOrig_n%d_rep%d.RData", nsuccesses, rep)))
  save(Y, file = here::here("data", sprintf("outputs_sirOutputFilling_WangOrig_n%d_rep%d.RData", nsuccesses, rep)))
  save(res, file = here::here("data", sprintf("fullOutput_OSFD_WangOrig_n%d_rep%d.RData", nsuccesses, rep)))
  cat(sprintf("Results saved with WangOrig_n=%d_rep%d suffix\n", nsuccesses, rep))
}

cat("\n=== PART 4 COMPLETE ===\n")
