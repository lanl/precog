# Part 4: OSFD with SIR-based initialization
# Pre-calculates initial samples using inverse SIR mappings and runs OSFD with those as initialization
# Author: AC Murph
# Date: Feb 2026

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript run_sir_part4_osfd_with_sir.R <nsuccesses>")
}
nsuccesses <- as.numeric(args[1])
cat(sprintf("Running Part 4 with nsuccesses = %d\n", nsuccesses))

# Load shared setup
source(here::here("R", "sir_experiment_setup.R"))

# -------------------------------
# OSFD using initial sample from INVERSE SIR maps
# -------------------------------

cat("\n=== PART 4: OSFD using initial sample from INVERSE SIR maps ===\n")

# Set cand_batch based on nsuccesses
cand_batch <- floor(nsuccesses / 2)
cat(sprintf("Using cand_batch = %d (floor(nsuccesses/2))\n", cand_batch))

# Pre-calculate LHS initial
n_ini_success = floor(2 * nsuccesses / 4)

# Create cache directory if it doesn't exist
cache_dir <- here::here("data", "feasible_LHSs")
if (!dir.exists(cache_dir)) {
  dir.create(cache_dir, recursive = TRUE)
  cat("Created cache directory:", cache_dir, "\n")
}

# Check if cached feasible LHS exists for the initial sample size
cache_file <- file.path(cache_dir, sprintf("feasible_LHS_%d.RData", n_ini_success))

if (file.exists(cache_file)) {
  cat(sprintf("Loading cached feasible LHS from: %s\n", cache_file))
  load(cache_file)  # This loads X_isfd
  cat(sprintf("Loaded %d cached feasible points\n", nrow(X_isfd)))
} else {
  cat(sprintf("No cached file found. Generating new feasible LHS with %d points...\n", n_ini_success))

  # Generate 10x oversample
  n_oversample <- 10 * n_ini_success
  X_candidate <- lhs::randomLHS(n_oversample, p)  # p=3: PIV, PIT, s0

  # Scale to bounds
  PIV_candidate <- piv_bounds[1] + X_candidate[, 1] * diff(piv_bounds)
  PIT_candidate <- pit_bounds[1] + X_candidate[, 2] * diff(pit_bounds)
  s0_candidate <- s0_bounds[1] + X_candidate[, 3] * diff(s0_bounds)

  # Test feasibility in parallel
  cl <- parallel::makeCluster(mc.cores)
  doSNOW::registerDoSNOW(cl)

  parallel::clusterExport(cl,
                         varlist = c("sir_piv_pit_feasible", "PIV_candidate",
                                     "PIT_candidate", "s0_candidate"),
                         envir = environment())

  cat(sprintf("Testing feasibility for %d candidate points...\n", n_oversample))
  pb <- txtProgressBar(max = n_oversample, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  feasible_flags <- foreach::foreach(
    i = 1:n_oversample,
    .combine = 'c',
    .inorder = TRUE,
    .options.snow = opts
  ) %dopar% {
    i0 <- 0.0001
    sir_piv_pit_feasible(PIV_candidate[i], PIT_candidate[i], s0_candidate[i], i0)
  }

  close(pb)
  parallel::stopCluster(cl)

  # Subset to feasible points
  X_feasible <- X_candidate[feasible_flags, ]
  n_feasible <- sum(feasible_flags)

  cat(sprintf("\nFeasible points: %d / %d (%.1f%%)\n",
              n_feasible, n_oversample, 100 * n_feasible / n_oversample))

  if (n_feasible < n_ini_success) {
    warning(sprintf("Only found %d feasible points but need %d. Consider increasing oversample factor.",
                    n_feasible, n_ini_success))
    X_isfd <- X_feasible
  } else {
    # Subsample to exactly n_ini_success
    keep_idx <- sample(n_feasible, n_ini_success)
    X_isfd <- X_feasible[keep_idx, ]
  }

  # Save X_isfd to cache for future runs
  save(X_isfd, file = cache_file)
  cat(sprintf("Saved feasible LHS (%d points) to cache: %s\n", nrow(X_isfd), cache_file))
}

# Always compute SIR mapping from X_isfd (whether loaded from cache or just generated)
cat("\n=== Computing SIR mapping from feasible LHS ===\n")

# Pre-calculate the LHS using the SIR mapping stuff:
make_f <- function(pit_bounds, piv_bounds, s0_bounds, reproduction_number_bounds) {
  force(pit_bounds); force(piv_bounds); force(s0_bounds); force(reproduction_number_bounds)
  function(x) {
    calculate_sir_from_PIVPIT(
      x[1], x[2], x[3],
      pit_bounds = pit_bounds,
      piv_bounds = piv_bounds,
      s0_bounds = s0_bounds,
      reproduction_number_bounds = reproduction_number_bounds
    )
  }
}
f2 <- make_f(pit_bounds, piv_bounds, s0_bounds, reproduction_number_bounds)

# Use doSNOW for progress bar support
cl <- parallel::makeCluster(mc.cores)
doSNOW::registerDoSNOW(cl)

# send the variable to workers
parallel::clusterExport(cl, varlist = c("f2", "q", "alpha_bounds", "reproduction_number_bounds", "s0_bounds"), envir = environment())
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

# Use a list to collect all results including failures
reses_list <- foreach::foreach(
  i = seq_len(nrow(X_isfd)),
  .inorder = TRUE,
  .options.snow = opts
) %dopar% {
  x <- X_isfd[i, , drop = FALSE]
  templist = safe_eval_f(f2, x, q+5)

  result <- list(
    index = i,
    input = x,
    success = templist$success,
    error_msg = templist$msg,
    data = NULL
  )

  if(templist$success){
    tmp_rho = templist$y[3] /  templist$y[4]
    tmp_alpha = templist$y[3]

    # Check bounds
    if(tmp_rho < min(reproduction_number_bounds)) {
      result$success <- FALSE
      result$error_msg <- sprintf("rho out of bounds: %.4f < %.4f", tmp_rho, min(reproduction_number_bounds))
    } else if(tmp_rho > max(reproduction_number_bounds)) {
      result$success <- FALSE
      result$error_msg <- sprintf("rho out of bounds: %.4f > %.4f", tmp_rho, max(reproduction_number_bounds))
    # } else if(tmp_alpha < min(alpha_bounds)) {
    #   result$success <- FALSE
    #   result$error_msg <- sprintf("alpha out of bounds: %.4f < %.4f", tmp_alpha, min(alpha_bounds))
    # } else if(tmp_alpha > max(alpha_bounds)) {
    #   result$success <- FALSE
    #   result$error_msg <- sprintf("alpha out of bounds: %.4f > %.4f", tmp_alpha, max(alpha_bounds))
    } else {
      result$data <- data.frame(
        piv = templist$y[1],
        pit = templist$y[2],
        alpha = templist$y[3],
        beta = templist$y[4],
        s0 = templist$y[7]
      )
    }
  }

  result
}

# Analyze failures
n_success <- sum(sapply(reses_list, function(x) x$success && !is.null(x$data)))
n_fail <- length(reses_list) - n_success

cat(sprintf("\n\nResults: %d successes, %d failures (%.1f%% failure rate)\n",
            n_success, n_fail, 100 * n_fail / length(reses_list)))

# Categorize error messages
if(n_fail > 0) {
  error_msgs <- sapply(reses_list, function(x) {
    if(!x$success || is.null(x$data)) x$error_msg else NA_character_
  })
  error_msgs <- error_msgs[!is.na(error_msgs)]

  cat("\nError summary:\n")
  error_table <- sort(table(error_msgs), decreasing = TRUE)
  for(i in seq_along(error_table)) {
    cat(sprintf("  [%d occurrences] %s\n", error_table[i], names(error_table)[i]))
  }
}

# Extract successful results for the rest of the script
reses <- do.call(rbind, lapply(reses_list, function(x) x$data))

close(pb)
parallel::stopCluster(cl)
cat(sprintf("\nCompleted %d tasks\n", n_tasks))

# Build initial data structures for OSFD
X_all = NULL
succ_flag = c()
fail_msg = c()
X_succ = NULL
Y_succ = NULL
remaining <- seq_len(CAND_size)  # Use CAND_size, not nrow(CAND) since we haven't created CAND yet

for(ii in 1:nrow(reses)){
  x_tmp = c(
    (reses[ii,]$alpha - alpha_bounds[1]) / diff(alpha_bounds),
    (reses[ii,]$alpha/reses[ii,]$beta - reproduction_number_bounds[1]) / diff(reproduction_number_bounds),
    (reses[ii,]$s0 - s0_bounds[1]) / diff(s0_bounds)
  )

  X_all <- rbind(X_all, matrix(x_tmp, nrow = 1) )
  succ_flag <- c(succ_flag, TRUE)
  fail_msg <- c(fail_msg, "")

  X_succ <- rbind(X_succ, matrix(x_tmp, nrow = 1))
  Y_succ <- rbind(Y_succ, matrix(c(reses[ii,]$piv, reses[ii,]$pit) , nrow = 1))
}

# Save pre-computed initialization data with nsuccesses in filename
save(X_all, file = paste0(here::here(), "/data/X_all.RData") )
save(succ_flag, file = paste0(here::here(), "/data/succ_flag.RData") )
save(fail_msg, file = paste0(here::here(), "/data/fail_msg.RData") )
save(X_succ, file = paste0(here::here(), "/data/X_succ.RData") )
save(Y_succ, file = paste0(here::here(), "/data/Y_succ.RData") )
save(remaining, file = paste0(here::here(), "/data/remaining.RData") )

# Now generate CAND for the OSFD run
CAND <- lhs::randomLHS(CAND_size, p)

# Define forward mapping function
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
  floor(nsuccesses - nrow(X_all)),
  n_ini_success = nrow(X_all),
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
  pkg_location = here::here(),
  pre_calculate_LHS = TRUE
)
print("We have now exited the OSFD function.")

D <- res$D_success
Y <- res$Y_success

# Save results with nsuccesses in filename
save(D, file = here::here("data", sprintf("inputs_sirOutputFilling_wSIRStuff_n%d.RData", nsuccesses)))
save(Y, file = here::here("data", sprintf("outputs_sirOutputFilling_wSIRStuff_n%d.RData", nsuccesses)))
save(res, file = here::here("data", sprintf("fullOutput_OSFD_wSIRStuff_n%d.RData", nsuccesses)))
cat(sprintf("Results saved with n=%d suffix\n", nsuccesses))

cat("\n=== PART 4 COMPLETE ===\n")
