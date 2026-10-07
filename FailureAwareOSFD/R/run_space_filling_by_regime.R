# Run space-filling samples on SIR space by regime
## Author: AC Murph
## Date: February 2026
library(parallel)
library(doParallel)
library(ggplot2)
library(patchwork)
library(zoo)
library(roll)
# Load helper functions
source(here::here("R", "sir_experiment_setup.R"))
source(here::here("R", "source_helpers.R"))
source_helpers()

setwd(here::here())

nsuccesses = 300
cand_batch = 300
n_replicates = 1
p = 3
mc.cores = 99
pkg_location = this.path::here()
# alpha_bounds = c(0.001, 10)
# reproduction_number_bounds = c(1.001, 30)
# s0_bounds = c(0.95, 0.999)
k = 5
h = 4
snippet_length = k + h
q = snippet_length

regimes_to_fill = c("Inc", "Dec", "Surge", "Near Peak", "Crash")

print(paste("Running space-filling for", length(regimes_to_fill), "regimes:", paste(regimes_to_fill, collapse=", ")))


# Loop through regimes to fill
for (regime_name in regimes_to_fill) {
  print(paste("========================================"))
  print(paste("Processing regime:", regime_name))
  print(paste("========================================"))

  # -------------------------------
  # -------------------------------
  ## LHS on output space:
  make_f <- function(alpha_bounds, reproduction_number_bounds, s0_bounds, regime_name, snippet_length, k) {
    force(alpha_bounds); force(reproduction_number_bounds); force(s0_bounds); force(regime_name); force(snippet_length); force(k)
    function(x) {
      calculate_sir_for_outputSpaceFilling_by_regime(
        x[1], x[2], x[3], x[4],
        regime_name = regime_name,
        snippet_length = snippet_length,
        alpha_bounds = alpha_bounds,
        reproduction_number_bounds = reproduction_number_bounds,
        s0_bounds = s0_bounds,
        # epsilon_bounds = epsilon_bounds,
        match_length = k,
        method = "OSFD"
      )
    }
  }
  f <- make_f(alpha_bounds, reproduction_number_bounds, s0_bounds, regime_name, snippet_length, k)

  CAND <- lhs::randomLHS(10000, p)
  print("About to enter OSFD function.")
  res <- constrained_osfd_ei_with_retry(
    f, p, q,
    CAND,
    nsuccesses,
    n_ini_success = floor(2 * nsuccesses / 4),
    batch_size = mc.cores,
    n_replicates = n_replicates,      # deterministic this time, i think.
    cand_batch = cand_batch,     # score this many candidates per iteration
    beta = 1.0,
    p_floor = 0,
    tau_hard = NULL,
    repel_radius = 0.05,    # input-space diversity scale
    update_feas_every = 1,
    mc.cores = mc.cores,
    verbose = TRUE,
    pkg_location = pkg_location
  )
  print("We have now exited the OSFD function.")

  # Based on the number of attempts used in the OSFD function, run ISFD the same number of times
  # Use X_all to match total attempts (not just successes), since OSFD includes retry logic
  nsuccesses_isfd = nrow(res$X_all)
  cat(sprintf("OSFD made %d total attempts (X_all) to get %d successes (D_success)\n",
              nrow(res$X_all), nrow(res$D_success)))
  cat(sprintf("ISFD will make %d attempts\n", nsuccesses_isfd))

  # -------------------------------
  # -------------------------------
  ## LHS directly on PIV/PIT space:
  X_isfd <- lhs::randomLHS(nsuccesses_isfd, p)
  make_f <- function(alpha_bounds, reproduction_number_bounds, s0_bounds, regime_name, snippet_length, k) {
    force(alpha_bounds); force(reproduction_number_bounds); force(s0_bounds); force(regime_name); force(snippet_length); force(k)
    function(x) {
      calculate_sir_for_outputSpaceFilling_by_regime(
        x[1], x[2], x[3], x[4],
        regime_name = regime_name,
        snippet_length = snippet_length,
        alpha_bounds = alpha_bounds,
        reproduction_number_bounds = reproduction_number_bounds,
        s0_bounds = s0_bounds,
        # epsilon_bounds = epsilon_bounds,
        match_length = k,
        method = 'ISFD'
      )
    }
  }
  f <- make_f(alpha_bounds, reproduction_number_bounds, s0_bounds, regime_name, snippet_length, k)

  cl <- parallel::makeCluster(mc.cores)
  doParallel::registerDoParallel(cl)

  # send the variable to workers
  parallel::clusterExport(cl, varlist = c("f", "q", "alpha_bounds", "reproduction_number_bounds", "s0_bounds"), envir = environment())
  parallel::clusterEvalQ(cl, {
    suppressPackageStartupMessages(library(deSolve))
    source(here::here("R", "source_helpers.R")); source_helpers()
    NULL
  })
  results <- foreach::foreach(
    i = seq_len(nrow(X_isfd)),
    .multicombine = TRUE,
    .combine = 'rbind',
    .inorder = TRUE
  ) %dopar% {
    x <- X_isfd[i, , drop = FALSE]
    templist = safe_eval_f(f, x, q)
    if(templist$success){
      xx = data.frame(matrix(templist$y[1:snippet_length], nrow = 1))
      names(xx) = as.character(1:snippet_length)
      return(xx)
    }
  }
  parallel::stopCluster(cl)
  # save(results, file=here::here("data", paste0("inputs_and_outputs_sirMappings_byRegime_", gsub("\\s+", "", regime_name), ".RData")))

  print(paste("Completed regime:", regime_name))
}

print("All regimes completed successfully!")





