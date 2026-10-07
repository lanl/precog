# Output-space filling design (no replicates) for stochastic SIR
# Author: AC Murph
# Date: Mar 2026

library(parallel)
library(doParallel)
library(doSNOW)
library(ggplot2)
library(patchwork)

# Load helper functions
source(here::here("R", "sir_experiment_setup.R"))
source(here::here("R", "source_helpers.R"))
source_helpers()

setwd(here::here())
pkg_location <- this.path::here()

# Get nsuccesses from command line argument (default to 100000 if not provided)
args <- commandArgs(trailingOnly = TRUE)
nsuccesses <- if (length(args) > 0) as.numeric(args[1]) else 100000

cat(sprintf("\n=== Running OSFD (no replicates) with nsuccesses = %d ===\n", nsuccesses))

p <- 3
q <- 2
mc.cores <- 99
cand_batch <- floor(nsuccesses / 2)

# Create objective function
make_f <- function(alpha_bounds, reproduction_number_bounds, s0_bounds) {
  force(alpha_bounds); force(reproduction_number_bounds); force(s0_bounds)
  function(x) {
    calculate_stochastic_sir_for_outputSpaceFilling(
      x[1], x[2], x[3],
      alpha_bounds = alpha_bounds,
      reproduction_number_bounds = reproduction_number_bounds,
      s0_bounds = s0_bounds
    )
  }
}
f <- make_f(alpha_bounds, reproduction_number_bounds, s0_bounds)

# Generate candidate set
CAND <- lhs::randomLHS(CAND_size, p)

# Run OSFD
print("About to enter OSFD function.")
res <- constrained_osfd_ei_with_retry(
  f, p, q,
  CAND,
  nsuccesses,
  n_ini_success = floor(2 * nsuccesses / 4),
  batch_size = mc.cores,
  n_replicates = 1,
  cand_batch = cand_batch,
  beta = 2.0,
  p_floor = 0.05,
  tau_hard = NULL,
  repel_radius = 0.05,
  update_feas_every = 1,
  mc.cores = mc.cores,
  verbose = TRUE,
  pkg_location = pkg_location
)
print("We have now exited the OSFD function.")

# Extract and save results (both in one save to avoid race condition)
D <- res$D_success
Y <- res$Y_success
save(D, Y, file = here::here("data", "inputs_and_outputs_StochasticSirOutputFilling_noReplicates.RData"))

cat(sprintf("\nSuccessfully saved D and Y to inputs_and_outputs_StochasticSirOutputFilling_noReplicates.RData\n"))
cat(sprintf("D dimensions: %d x %d\n", nrow(D), ncol(D)))
cat(sprintf("Y dimensions: %d x %d\n\n", nrow(Y), ncol(Y)))
