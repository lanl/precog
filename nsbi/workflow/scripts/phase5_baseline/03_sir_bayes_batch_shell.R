#!/usr/bin/env Rscript
# Bayesian Batch SIR Fitting with CmdStanR - Shell wrapper

suppressPackageStartupMessages({
  library(cmdstanr)
  library(posterior)
  library(dplyr)
})

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 17) {
  stop("Expected 17 arguments")
}

batch_id <- args[1]
params_file <- args[2]
stan_model_file <- args[3]
output_file <- args[4]
log_file <- args[5]
sim_ids_str <- args[6]
R0_min <- as.numeric(args[7])
R0_max <- as.numeric(args[8])
rt_min <- as.numeric(args[9])
rt_max <- as.numeric(args[10])
S_fixed <- as.numeric(args[11])
prior_buffer <- as.numeric(args[12])
n_chains <- as.integer(args[13])
n_warmup <- as.integer(args[14])
n_iter <- as.integer(args[15])
adapt_delta <- as.numeric(args[16])
max_treedepth <- as.integer(args[17])

# Parse sim_ids (passed as R vector string like 'c("test_p00000_r0", ...)')
sim_ids <- eval(parse(text = sim_ids_str))

# Configure CmdStan
conda_prefix <- Sys.getenv("CONDA_PREFIX")
if (conda_prefix != "") {
  cmdstan_path <- file.path(conda_prefix, "bin", "cmdstan")
  if (dir.exists(cmdstan_path)) {
    set_cmdstan_path(cmdstan_path)
  }
}

cat("\n========================================\n")
cat("SIR Bayesian Batch Fitting\n")
cat("========================================\n\n")

sink(log_file, split = TRUE)

# Load parameters
params_df <- read.csv(params_file)

# Extract param_id and replicate from sim_id
# sim_id format: "test_p00123_r0"
parse_sim_id <- function(sim_id) {
  parts <- strsplit(sim_id, "_")[[1]]
  param_id <- as.integer(sub("p", "", parts[2]))
  replicate <- as.integer(sub("r", "", parts[3]))
  list(param_id = param_id, replicate = replicate)
}

# Compile model
cat("Compiling Stan model...\n")
cat("Stan model file:", stan_model_file, "\n")
cat("File exists:", file.exists(stan_model_file), "\n")
cat("CmdStan path:", cmdstan_path(), "\n\n")

stan_model <- tryCatch({
  cmdstan_model(stan_model_file)
}, error = function(e) {
  cat("\n!!! STAN MODEL COMPILATION FAILED !!!\n")
  cat("Error:", conditionMessage(e), "\n")
  cat("\nChecking environment:\n")
  cat("CONDA_PREFIX:", Sys.getenv("CONDA_PREFIX"), "\n")
  cat("CMDSTAN_OUTPUT_DIR:", Sys.getenv("CMDSTAN_OUTPUT_DIR"), "\n")
  cat("STAN_NUM_THREADS:", Sys.getenv("STAN_NUM_THREADS"), "\n")
  cat("Working directory:", getwd(), "\n")
  stop("Stan model compilation failed")
})

cat("Model compiled successfully\n")
cat("Model executable:", stan_model$exe_file(), "\n\n")

# Fit function
fit_one_simulation <- function(sim_id) {
  cat(sprintf("\n--- %s ---\n", sim_id))
  
  # Parse sim_id to get param_id and replicate
  parsed <- parse_sim_id(sim_id)
  param_id <- parsed$param_id
  replicate <- parsed$replicate
  
  # Initialize result with metadata
  result <- data.frame(
    sim_id = sim_id,
    param_id = param_id,
    replicate = replicate,
    converged = FALSE,
    error_msg = NA_character_,
    stringsAsFactors = FALSE
  )
  
  # Try-catch wrapper for error handling
  tryCatch({
    # Get true parameters
    true_params <- params_df[params_df$param_id == param_id, ]
    if (nrow(true_params) == 0) {
      stop(sprintf("Parameters not found for param_id %d", param_id))
    }
    
    true_R0 <- true_params$R0
    true_recovery_time <- true_params$recovery_time
    
    # Store true parameters in result immediately (for error handling)
    result$R0_true <- true_R0
    result$D_true <- true_recovery_time
    
    cat(sprintf("True: R0=%.3f, recovery_time=%.3f\n", true_R0, true_recovery_time))
    
    # Load case data using binned case counts from phase2
    binned_file <- file.path("results", "phase2_processing", "case_counts", 
                             paste0(sim_id, "_binned_case_counts.csv"))
    if (!file.exists(binned_file)) {
      stop(sprintf("Case file not found: %s", binned_file))
    }
    cases <- read.csv(binned_file)
    
    # Stan data
    # FIXED: Pass actual day values [0, 1, 2, ...]
    # Stan model extends this to [0, 1, 2, ..., n+1] for diff(R) alignment
    stan_data <- list(
      n_obs = nrow(cases),
      y = as.integer(cases$case_count),  # Stan expects integer array
      ts = cases$day,                     # Day values: 0, 1, 2, ...
      t0 = 0.0,                           # Initial time
      S0_fixed = S_fixed,                 # Initial susceptible
      I0_fixed = 1.0,                     # Initial infected (always 1)
      R0_min = R0_min,
      R0_max = R0_max,
      recovery_time_min = rt_min,
      recovery_time_max = rt_max,
      prior_buffer = prior_buffer
    )
    
    # MCMC
    cat(sprintf("Running MCMC (%d chains, %d warmup, %d samples)...\n", 
                n_chains, n_warmup, n_iter))
    
    # Run MCMC sampling
    fit <- stan_model$sample(
      data = stan_data,
      chains = n_chains,
      parallel_chains = n_chains,
      iter_warmup = n_warmup,
      iter_sampling = n_iter - n_warmup,
      adapt_delta = adapt_delta,
      max_treedepth = max_treedepth,
      refresh = 100,
      show_messages = TRUE,
      show_exceptions = TRUE
    )
    
    # Check if chains finished successfully
    if (!is.null(fit$output_files())) {
      cat("Stan output files created:\n")
      print(fit$output_files())
      cat("\n")
    }
    
    # Try to extract draws and catch failures
    draws <- tryCatch({
      fit$draws(variables = c("R0", "recovery_time", "beta", "gamma"), format = "df")
    }, error = function(e) {
      cat("\n!!! MCMC CHAINS FAILED !!!\n")
      cat("Error when extracting draws:", conditionMessage(e), "\n\n")
      
      # Read Stan output files to get actual error
      cat("Stan model executable:", stan_model$exe_file(), "\n")
      cat("Executable exists:", file.exists(stan_model$exe_file()), "\n\n")
      
      # Get and read output CSV files
      csv_files <- fit$output_files()
      if (!is.null(csv_files) && length(csv_files) > 0) {
        cat("Reading Stan chain outputs for error details...\n\n")
        for (i in seq_along(csv_files)) {
          csv_file <- csv_files[i]
          if (file.exists(csv_file)) {
            cat(sprintf("=== Chain %d: %s ===\n", i, basename(csv_file)))
            # Read first 100 lines to capture error messages
            lines <- tryCatch({
              readLines(csv_file, n = 100, warn = FALSE)
            }, error = function(e2) {
              c(sprintf("ERROR reading file: %s", conditionMessage(e2)))
            })
            cat(paste(lines, collapse = "\n"), "\n\n")
          } else {
            cat(sprintf("=== Chain %d: FILE NOT FOUND: %s ===\n\n", i, csv_file))
          }
        }
      } else {
        cat("No Stan output files found!\n")
      }
      
      return(NULL)
    })
    
    # Check if draws extraction failed
    if (is.null(draws)) {
      cat("\nSkipping", sim_id, "due to MCMC failure\n\n")
      result$converged <- FALSE
      result$error_msg <- "MCMC chains failed - no draws extracted"
      # True params already populated above; fill estimates with NA
      result$R0_point_estimate <- NA_real_
      result$D_point_estimate <- NA_real_
      result$R0_lower_bound <- NA_real_
      result$R0_upper_bound <- NA_real_
      result$D_lower_bound <- NA_real_
      result$D_upper_bound <- NA_real_
      result$nll <- NA_real_
      result$rhat_R0 <- NA_real_
      result$rhat_recovery_time <- NA_real_
      result$rhat_beta <- NA_real_
      result$rhat_gamma <- NA_real_
      result$ess_bulk_R0 <- NA_real_
      result$ess_bulk_recovery_time <- NA_real_
      result$ess_bulk_beta <- NA_real_
      result$ess_bulk_gamma <- NA_real_
      result$ess_tail_R0 <- NA_real_
      result$ess_tail_recovery_time <- NA_real_
      result$n_divergences <- NA_real_
      result$n_max_treedepth <- NA_real_
      return(result)
    }
    
    # Extract summary statistics and diagnostics
    summary_stats <- fit$summary(variables = c("R0", "recovery_time", "beta", "gamma"))
    diagnostics <- fit$diagnostic_summary()
    
    # Get statistics for each parameter
    R0_stats <- summary_stats[summary_stats$variable == "R0", ]
    rt_stats <- summary_stats[summary_stats$variable == "recovery_time", ]
    beta_stats <- summary_stats[summary_stats$variable == "beta", ]
    gamma_stats <- summary_stats[summary_stats$variable == "gamma", ]
    
    # Point estimates (use posterior mean to match MLE expectation)
    est_R0 <- mean(draws$R0)
    est_recovery_time <- mean(draws$recovery_time)
    
    # 95% credible intervals (HPD-style, using quantiles)
    ci_R0 <- quantile(draws$R0, c(0.025, 0.975))
    ci_rt <- quantile(draws$recovery_time, c(0.025, 0.975))
    
    # Convergence: Rhat < 1.1 AND ESS > 400 for both main parameters
    rhat_ok <- (R0_stats$rhat < 1.1) && (rt_stats$rhat < 1.1)
    ess_ok <- (R0_stats$ess_bulk > 400) && (rt_stats$ess_bulk > 400)
    converged <- rhat_ok && ess_ok
    
    cat(sprintf("Posterior R0: %.3f [%.3f, %.3f]\n", est_R0, ci_R0[1], ci_R0[2]))
    cat(sprintf("Posterior RT: %.3f [%.3f, %.3f]\n", est_recovery_time, ci_rt[1], ci_rt[2]))
    cat(sprintf("Converged: %s (Rhat OK: %s, ESS OK: %s)\n", converged, rhat_ok, ess_ok))
    
    # Populate result with standard schema matching MLE output
    result$converged <- converged
    result$error_msg <- NA_character_
    
    # True values (already populated above, but keep for clarity)
    result$R0_true <- true_R0
    result$D_true <- true_recovery_time
    
    # Point estimates (standardized names)
    result$R0_point_estimate <- est_R0
    result$D_point_estimate <- est_recovery_time
    
    # Credible intervals (standardized names)
    result$R0_lower_bound <- as.numeric(ci_R0[1])
    result$R0_upper_bound <- as.numeric(ci_R0[2])
    result$D_lower_bound <- as.numeric(ci_rt[1])
    result$D_upper_bound <- as.numeric(ci_rt[2])
    
    # NLL approximation (use -mean(log_lik) if available, else NA)
    # For now, set to NA as Stan model doesn't compute log_lik by default
    result$nll <- NA_real_
    
    # Bayesian-specific diagnostics (keep for aggregation script)
    result$rhat_R0 <- R0_stats$rhat
    result$rhat_recovery_time <- rt_stats$rhat
    result$rhat_beta <- beta_stats$rhat
    result$rhat_gamma <- gamma_stats$rhat
    result$ess_bulk_R0 <- R0_stats$ess_bulk
    result$ess_bulk_recovery_time <- rt_stats$ess_bulk
    result$ess_bulk_beta <- beta_stats$ess_bulk
    result$ess_bulk_gamma <- gamma_stats$ess_bulk
    result$ess_tail_R0 <- R0_stats$ess_tail
    result$ess_tail_recovery_time <- rt_stats$ess_tail
    result$n_divergences <- sum(diagnostics$num_divergent, na.rm = TRUE)
    result$n_max_treedepth <- sum(diagnostics$num_max_treedepth, na.rm = TRUE)
    
  }, error = function(e) {
    cat(sprintf("ERROR: %s\n", e$message))
    result$error_msg <<- e$message
    result$converged <<- FALSE
    # Only set to NA if not already populated (early errors vs late errors)
    # True params populated after line 113, so early errors get NA, late errors keep values
    if (!"R0_true" %in% names(result)) result$R0_true <<- NA_real_
    if (!"D_true" %in% names(result)) result$D_true <<- NA_real_
    # Estimates are always NA on error
    if (!"R0_point_estimate" %in% names(result)) result$R0_point_estimate <<- NA_real_
    if (!"D_point_estimate" %in% names(result)) result$D_point_estimate <<- NA_real_
    if (!"R0_lower_bound" %in% names(result)) result$R0_lower_bound <<- NA_real_
    if (!"R0_upper_bound" %in% names(result)) result$R0_upper_bound <<- NA_real_
    if (!"D_lower_bound" %in% names(result)) result$D_lower_bound <<- NA_real_
    if (!"D_upper_bound" %in% names(result)) result$D_upper_bound <<- NA_real_
    if (!"nll" %in% names(result)) result$nll <<- NA_real_
    if (!"rhat_R0" %in% names(result)) result$rhat_R0 <<- NA_real_
    if (!"rhat_recovery_time" %in% names(result)) result$rhat_recovery_time <<- NA_real_
    if (!"rhat_beta" %in% names(result)) result$rhat_beta <<- NA_real_
    if (!"rhat_gamma" %in% names(result)) result$rhat_gamma <<- NA_real_
    if (!"ess_bulk_R0" %in% names(result)) result$ess_bulk_R0 <<- NA_real_
    if (!"ess_bulk_recovery_time" %in% names(result)) result$ess_bulk_recovery_time <<- NA_real_
    if (!"ess_bulk_beta" %in% names(result)) result$ess_bulk_beta <<- NA_real_
    if (!"ess_bulk_gamma" %in% names(result)) result$ess_bulk_gamma <<- NA_real_
    if (!"ess_tail_R0" %in% names(result)) result$ess_tail_R0 <<- NA_real_
    if (!"ess_tail_recovery_time" %in% names(result)) result$ess_tail_recovery_time <<- NA_real_
    if (!"n_divergences" %in% names(result)) result$n_divergences <<- NA_real_
    if (!"n_max_treedepth" %in% names(result)) result$n_max_treedepth <<- NA_real_
  })
  
  return(result)
}

results_list <- lapply(sim_ids, fit_one_simulation)

# Combine results with better error handling
# Use dplyr::bind_rows which handles mismatched columns gracefully
results_df <- tryCatch({
  bind_rows(results_list)
}, error = function(e) {
  cat("\nWARNING: bind_rows failed, trying rbind:\n")
  cat(e$message, "\n")
  do.call(rbind, results_list)
})

# Check if any results succeeded
n_success <- sum(results_df$converged, na.rm = TRUE)
n_failed <- sum(!results_df$converged, na.rm = TRUE)

write.csv(results_df, output_file, row.names = FALSE)

cat(sprintf("\nBatch complete: %d simulations (%d succeeded, %d failed)\n", 
            nrow(results_df), n_success, n_failed))
sink()
