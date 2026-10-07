# Regular LHS on input space (with replicates) for stochastic SIR
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

cat(sprintf("\n=== Running LHS (with replicates) with nsuccesses = %d ===\n", nsuccesses))

p <- 3
q <- 2
n_replicates <- 5
mc.cores <- 99

# Generate LHS design with replicates
X_isfd <- lhs::randomLHS(ceiling(nsuccesses / n_replicates), p)
X_isfd_full <- NULL
for (tmp_idx in 1:n_replicates) {
  X_isfd_full <- rbind(X_isfd_full, X_isfd)
}
X_isfd <- X_isfd_full

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

# Set up parallel cluster
cl <- parallel::makeCluster(mc.cores)
doSNOW::registerDoSNOW(cl)

# Send variables to workers
parallel::clusterExport(cl, varlist = c("f", "q"), envir = environment())
parallel::clusterEvalQ(cl, {
  suppressPackageStartupMessages(library(deSolve))
  source(here::here("R", "sir_experiment_setup.R"))
  source(here::here("R", "source_helpers.R"))
  source_helpers()
  NULL
})

# Set up progress bar
cat("Running LHS on input space (with replicates)...\n")
pb <- txtProgressBar(max = nrow(X_isfd), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Run parallel evaluation
results <- foreach::foreach(
  i = seq_len(nrow(X_isfd)),
  .combine = 'rbind',
  .multicombine = TRUE,
  .inorder = TRUE,
  .options.snow = opts
) %dopar% {
  x <- X_isfd[i, , drop = FALSE]

  alpha <- alpha_bounds[1] + x[1] * diff(alpha_bounds)
  rho <- reproduction_number_bounds[1] + x[2] * diff(reproduction_number_bounds)
  beta <- alpha / rho

  templist <- safe_eval_f(f, x, q)
  if (templist$success) {
    data.frame(piv = templist$y[1], pit = templist$y[2], alpha = alpha, beta = beta)
  }
}
close(pb)
parallel::stopCluster(cl)

# Save results
save(results, file = here::here("data", "inputs_and_outputs_sirStochasticLHS_wReplicates.RData"))
cat(sprintf("\nSuccessfully saved results to inputs_and_outputs_sirStochasticLHS_wReplicates.RData\n"))
cat(sprintf("Total successful evaluations: %d\n\n", nrow(results)))
