# Part 2: Regular LHS on input space (alpha/beta)
# Standard LHS directly on input space with forward simulation to get outputs
# Author: AC Murph
# Date: Feb 2026

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript run_sir_part2_basic_lhs.R <nsuccesses> [rep]")
}
nsuccesses <- as.numeric(args[1])
rep_num <- if (length(args) >= 2) as.numeric(args[2]) else NULL
rep_suffix <- if (!is.null(rep_num)) sprintf("_rep%d", rep_num) else ""

if (!is.null(rep_num)) {
  cat(sprintf("Running Part 2 with nsuccesses = %d, replicate = %d\n", nsuccesses, rep_num))
} else {
  cat(sprintf("Running Part 2 with nsuccesses = %d\n", nsuccesses))
}

# Load shared setup
source(here::here("R", "sir_experiment_setup.R"))

# -------------------------------
# Regular LHS on input space (alpha/beta)
# -------------------------------

cat("\n=== PART 2: Regular LHS on input space (alpha/beta) ===\n")

X_isfd <- lhs::randomLHS(nsuccesses, p)

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

cl <- parallel::makeCluster(mc.cores)
doSNOW::registerDoSNOW(cl)

# send the variable to workers
parallel::clusterExport(cl, varlist = c("f", "q"), envir = environment())
parallel::clusterEvalQ(cl, {
  suppressPackageStartupMessages(library(deSolve))
  source(here::here("R", "source_helpers.R")); source_helpers()
  NULL
})

# Set up progress bar
n_tasks <- nrow(X_isfd)
cat(sprintf("Processing %d tasks across %d cores...\n", n_tasks, mc.cores))
pb <- txtProgressBar(max = n_tasks, style = 3)
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
  s0 = s0_bounds[1] + x[3] * diff(s0_bounds)

  templist = safe_eval_f(f, x, q)
  if(templist$success){
    data.frame(piv = templist$y[1], pit = templist$y[2], alpha = alpha, beta = beta, s0 = s0)
  }
}

close(pb)
parallel::stopCluster(cl)

cat(sprintf("Completed %d tasks\n\n", n_tasks))

# Save results with nsuccesses in filename
output_file <- here::here("data", sprintf("inputs_and_outputs_sirBasicLHS_n%d%s.RData", nsuccesses, rep_suffix))
save(results, file = output_file)
cat(sprintf("Results saved to: %s\n", output_file))

cat("\n=== PART 2 COMPLETE ===\n")
