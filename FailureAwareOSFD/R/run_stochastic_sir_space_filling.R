# For this script -- I am going to output-space fill on PIT and PIV to SIR curves. Stochastic version.
# This will use three methods.  1) basic LHS, 2) my SIR paper stuff, 3) this output-space filling design stuff
## Author: AC Murph
## Date Mar 2026
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
pkg_location = this.path::here()

p = 3
q = 2

nsuccesses = 50000 
n_replicates = 5
mc.cores = 99
cand_batch <- floor(nsuccesses / 2)

# alpha_bounds = c(0.001, 100)
# reproduction_number_bounds = c(1.001, 100)
# pit_bounds = c(0, 50)
# piv_bounds = c(0.005, 0.75)
# s0_bounds = c(0.95, 0.999)

# -------------------------------
# -------------------------------
## Regular LHS on input space (alpha/beta)
X_isfd <- lhs::randomLHS(nsuccesses, p) 
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

cl <- parallel::makeCluster(mc.cores)
doSNOW::registerDoSNOW(cl)

# send the variable to workers
parallel::clusterExport(cl, varlist = c("f", "q"), envir = environment())
parallel::clusterEvalQ(cl, {
  suppressPackageStartupMessages(library(deSolve))
  source(here::here("R", "sir_experiment_setup.R")); source(here::here("R", "source_helpers.R")); source_helpers()
  NULL
})

# Set up progress bar for first loop (no replicates)
cat("Running LHS on input space (no replicates)...\n")
pb <- txtProgressBar(max = nrow(X_isfd), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

results <- foreach::foreach(
  i = seq_len(nrow(X_isfd)),
  .combine = 'rbind',
  .multicombine = TRUE,
  .inorder = TRUE,
  .options.snow = opts
) %dopar% {
  x <- X_isfd[i, , drop = FALSE]

  alpha = alpha_bounds[1] + x[1] * diff(alpha_bounds)
  rho = reproduction_number_bounds[1] + x[2] * diff(reproduction_number_bounds)
  beta = alpha / rho

  templist = safe_eval_f(f, x, q)
  if(templist$success){
    data.frame(piv = templist$y[1], pit = templist$y[2], alpha = alpha, beta = beta)
  }
}
close(pb)
parallel::stopCluster(cl)
save(results, file=here::here("data", "inputs_and_outputs_sirStochasticLHS_noReplicates.RData"))

# -------------------------------
# -------------------------------
# Now run with replicates
X_isfd <- lhs::randomLHS(ceiling(nsuccesses/n_replicates), p) 
X_isfd_full = NULL
for(tmp_idx in 1:(n_replicates)){
  X_isfd_full = rbind(X_isfd_full, X_isfd)
}
X_isfd = X_isfd_full

cl <- parallel::makeCluster(mc.cores)
doSNOW::registerDoSNOW(cl)

# send the variable to workers
parallel::clusterExport(cl, varlist = c("f", "q"), envir = environment())
parallel::clusterEvalQ(cl, {
  suppressPackageStartupMessages(library(deSolve))
  source(here::here("R", "sir_experiment_setup.R")); source(here::here("R", "source_helpers.R")); source_helpers()
  NULL
})

# Set up progress bar for second loop (with replicates)
cat("Running LHS on input space (with replicates)...\n")
pb <- txtProgressBar(max = nrow(X_isfd), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

results <- foreach::foreach(
  i = seq_len(nrow(X_isfd)),
  .combine = 'rbind',
  .multicombine = TRUE,
  .inorder = TRUE,
  .options.snow = opts
) %dopar% {
  x <- X_isfd[i, , drop = FALSE]

  alpha = alpha_bounds[1] + x[1] * diff(alpha_bounds)
  rho = reproduction_number_bounds[1] + x[2] * diff(reproduction_number_bounds)
  beta = alpha / rho

  templist = safe_eval_f(f, x, q)
  if(templist$success){
    data.frame(piv = templist$y[1], pit = templist$y[2], alpha = alpha, beta = beta)
  }
}
close(pb)
parallel::stopCluster(cl)
save(results, file=here::here("data", "inputs_and_outputs_sirStochasticLHS_wReplicates.RData"))



# -------------------------------
# -------------------------------
## LHS on output space (with replicates):
CAND <- lhs::randomLHS(CAND_size, p)
print("About to enter OSFD function.")
res <- constrained_osfd_ei_with_retry(
  f, p, q,
  CAND,
  nsuccesses,
  n_ini_success = floor(2 * nsuccesses / 4),
  batch_size = mc.cores,
  n_replicates = n_replicates,      # deterministic this time, i think. 
  cand_batch = cand_batch,     # score this many candidates per iteration
  beta = 1,
  p_floor = 0,
  tau_hard = NULL,
  repel_radius = 0.05,    # input-space diversity scale
  update_feas_every = 1,
  mc.cores = mc.cores,
  verbose = TRUE, 
  pkg_location = pkg_location
)
print("We have now exited the OSFD function.")

D <- res$D_success
Y <- res$Y_success
save(D, file=here::here("data", "inputs_StochasticSirOutputFilling_wReplicates.RData"))
save(Y, file=here::here("data", "outputs_StochasticSirOutputFilling_wReplicates.RData"))

# -------------------------------
# -------------------------------
## LHS on output space (no replicates):
CAND <- lhs::randomLHS(CAND_size, p)
print("About to enter OSFD function.")
res <- constrained_osfd_ei_with_retry(
  f, p, q,
  CAND,
  nsuccesses,
  n_ini_success = floor(2 * nsuccesses / 4),
  batch_size = mc.cores,
  n_replicates = 1,      # deterministic this time, i think. 
  cand_batch = cand_batch,     # score this many candidates per iteration
  beta = 2.0,
  p_floor = 0.05,
  tau_hard = NULL,
  repel_radius = 0.05,    # input-space diversity scale
  update_feas_every = 1,
  mc.cores = mc.cores,
  verbose = TRUE, 
  pkg_location = pkg_location
)
print("We have now exited the OSFD function.")

D <- res$D_success
Y <- res$Y_success
save(D, file=here::here("data", "inputs_StochasticSirOutputFilling_noReplicates.RData"))
save(Y, file=here::here("data", "outputs_StochasticSirOutputFilling_noReplicates.RData"))