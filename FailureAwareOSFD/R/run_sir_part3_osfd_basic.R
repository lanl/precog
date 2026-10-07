# Part 3: Direct OSFD (Output Space Filling Design)
# Runs constrained_osfd_ei_with_retry from scratch
# Author: AC Murph
# Date: Feb 2026

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript run_sir_part3_osfd_basic.R <nsuccesses> [budget] [rep]")
}
nsuccesses <- as.numeric(args[1])

# Optional: computational budget (max_attempts)
budget <- NULL
if (length(args) >= 2) {
  budget <- as.numeric(args[2])
  cat(sprintf("Running Part 3 with computational budget = %d\n", budget))
} else {
  cat(sprintf("Running Part 3 with nsuccesses = %d\n", nsuccesses))
}

# Optional: replicate number
rep_num <- if (length(args) >= 3) as.numeric(args[3]) else NULL
rep_suffix <- if (!is.null(rep_num)) sprintf("_rep%d", rep_num) else ""

if (!is.null(rep_num)) {
  cat(sprintf("  Replicate = %d\n", rep_num))
}

# Load shared setup
source(here::here("R", "sir_experiment_setup.R"))

# -------------------------------
# Direct OSFD
# -------------------------------

cat("\n=== PART 3: Direct OSFD ===\n")

# Set parameters based on whether budget is provided
if (!is.null(budget)) {
  # Budget mode: run until max_attempts reached
  n_target <- 1e10  # Very large number to effectively disable success-based stopping
  n_ini_success <- floor(budget / 2)
  max_attempts <- budget
  cand_batch <- CAND_size #floor(budget / 2)
  cat(sprintf("Budget mode: max_attempts = %d, n_ini_success = %d, cand_batch = %d\n",
              max_attempts, n_ini_success, cand_batch))
} else {
  # Success mode: run until nsuccesses reached
  n_target <- nsuccesses
  n_ini_success <- floor(2 * nsuccesses / 4)
  max_attempts <- NULL
  cand_batch <- CAND_size #floor(nsuccesses / 2)
  cat(sprintf("Success mode: nsuccesses = %d, n_ini_success = %d, cand_batch = %d\n",
              n_target, n_ini_success, cand_batch))
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

print("About to enter OSFD function.")
res <- constrained_osfd_ei_with_retry(
  f, p, q,
  CAND,
  n_success = n_target,
  n_ini_success = n_ini_success,
  batch_size = min(mc.cores, n_ini_success),
  n_replicates = 1,      # deterministic this time, i think.
  cand_batch = cand_batch,     # score this many candidates per iteration
  beta = 1,
  p_floor = 0,
  tau_hard = NULL,
  repel_radius = 0.05,    # input-space diversity scale
  update_feas_every = 1,
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

# Save results with appropriate suffix
if (!is.null(budget)) {
  # Budget mode: use budget in filename
  save(D, file = here::here("data", sprintf("inputs_sirOutputFilling_budget%d_n%d%s.RData", budget, actual_n, rep_suffix)))
  save(Y, file = here::here("data", sprintf("outputs_sirOutputFilling_budget%d_n%d%s.RData", budget, actual_n, rep_suffix)))
  save(res, file = here::here("data", sprintf("fullOutput_OSFD_basic_budget%d_n%d%s.RData", budget, actual_n, rep_suffix)))
  cat(sprintf("Results saved with budget=%d suffix (obtained %d successes)%s\n", budget, actual_n, rep_suffix))
} else {
  # Success mode: use nsuccesses in filename (original behavior)
  save(D, file = here::here("data", sprintf("inputs_sirOutputFilling_n%d%s.RData", nsuccesses, rep_suffix)))
  save(Y, file = here::here("data", sprintf("outputs_sirOutputFilling_n%d%s.RData", nsuccesses, rep_suffix)))
  save(res, file = here::here("data", sprintf("fullOutput_OSFD_basic_n%d%s.RData", nsuccesses, rep_suffix)))
  cat(sprintf("Results saved with n=%d suffix%s\n", nsuccesses, rep_suffix))
}

cat("\n=== PART 3 COMPLETE ===\n")
