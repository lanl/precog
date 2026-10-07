#!/usr/bin/env Rscript

#' Modular Batch SIR MLE Fitting
#' 
#' Fit SIR model to binned daily case counts via maximum likelihood estimation.
#' 
#' Supports composite strategies: {parameterization}_{starting_values}
#' - Parameterizations: "betagamma" (original), "r0rt" (new, constrained R0)
#' - Starting values: "prior_mean", "true_values"
#' 
#' Examples: "betagamma_prior_mean", "r0rt_true_values"

suppressPackageStartupMessages({
  library(deSolve)
  library(bbmle)
  library(mvtnorm)
})

# Load parameterization module
source(file.path(snakemake@scriptdir, "parameterizations.R"))

# ============================================================================
# Setup
# ============================================================================

cat("\n========================================\n")
cat("SIR MLE Batch Fitting (Modular)\n")
cat("========================================\n\n")

# Get snakemake parameters
sim_ids <- snakemake@params$sim_ids
batch_id <- snakemake@wildcards$batch_id
composite_strategy <- snakemake@params$starting_value_strategy
params_file <- snakemake@input$params
output_file <- snakemake@output[[1]]
log_file <- snakemake@log[[1]]

# Prior bounds for parameter constraints and starting values
R0_min <- snakemake@params$R0_min
R0_max <- snakemake@params$R0_max
rt_min <- snakemake@params$recovery_time_min
rt_max <- snakemake@params$recovery_time_max
S_fixed <- snakemake@params$S_fixed

# Parse composite strategy: {parameterization}_{starting_values}
# Examples: "betagamma_prior_mean", "r0rt_true_values"
strategy_parts <- strsplit(composite_strategy, "_", fixed = TRUE)[[1]]
if (length(strategy_parts) >= 2) {
  parameterization_name <- strategy_parts[1]
  starting_value_strategy <- paste(strategy_parts[-1], collapse = "_")
} else {
  # Backward compatibility: if no underscore, assume betagamma
  parameterization_name <- "betagamma"
  starting_value_strategy <- composite_strategy
}

# Load parameterization
parameterization <- get_parameterization(parameterization_name)

cat(sprintf("Batch ID: %s\n", batch_id))
cat(sprintf("Composite strategy: %s\n", composite_strategy))
cat(sprintf("  Parameterization: %s\n", parameterization_name))
cat(sprintf("  Starting values: %s\n", starting_value_strategy))
cat(sprintf("Processing %d simulations\n", length(sim_ids)))
cat(sprintf("S fixed at: %d\n", S_fixed))
cat(sprintf("Prior ranges:\n"))
cat(sprintf("  R0: [%.2f, %.2f]\n", R0_min, R0_max))
cat(sprintf("  Recovery time: [%.2f, %.2f]\n", rt_min, rt_max))
cat("\n")

# Redirect output to log file
log_conn <- file(log_file, open = "wt")
sink(log_conn, type = "output")
sink(log_conn, type = "message")

# ============================================================================
# SIR Model Definition
# ============================================================================

sir_ode <- function(t, y, params) {
  S <- y[1]
  I <- y[2]
  R <- y[3]
  beta <- params["beta"]
  gamma <- params["gamma"]
  N <- S + I + R  # Total population (constant)
  
  # Frequency-dependent transmission (matches BEAST2)
  dS <- -(beta/N) * S * I
  dI <- (beta/N) * S * I - gamma * I
  dR <- gamma * I
  
  list(c(dS, dI, dR))
}

# ============================================================================
# Helper Functions
# ============================================================================

load_case_counts <- function(sim_id, S0) {
  #' Load pre-binned case counts from Phase 2 output
  
  binned_file <- file.path("results", "phase2_processing", "case_counts", paste0(sim_id, "_binned_case_counts.csv"))
  
  if (!file.exists(binned_file)) {
    stop(sprintf("Binned case counts file not found: %s", binned_file))
  }
  
  df <- read.csv(binned_file)
  
  # Standard initial conditions: I=1, R=0 at t=0
  I0 <- 1
  R0 <- 0
  
  list(
    daily_cases = df$case_count,
    times = df$day,
    S0 = S0,
    I0 = I0,
    R0_init = R0
  )
}


fit_sir_mle <- function(sim_id, true_beta, true_gamma, true_R0, true_recovery_time,
                       starting_value_strategy, S0, parameterization) {
  #' Fit SIR model via maximum likelihood using specified parameterization
  
  result <- list(
    sim_id = sim_id,
    converged = FALSE,
    error_msg = NA,
    true_beta = true_beta,  # Removed in aggregation
    true_gamma = true_gamma,  # Removed in aggregation
    R0_true = true_R0,  # Standardized name
    D_true = true_recovery_time,  # Standardized name
    starting_value_strategy = composite_strategy,  # Removed in aggregation
    est_beta = NA,  # Removed in aggregation
    est_gamma = NA,  # Removed in aggregation
    R0_point_estimate = NA,  # Standardized name
    D_point_estimate = NA,  # Standardized name
    ci_lower_beta = NA,  # Removed in aggregation
    ci_upper_beta = NA,  # Removed in aggregation
    ci_lower_gamma = NA,  # Removed in aggregation
    ci_upper_gamma = NA,  # Removed in aggregation
    R0_lower_bound = NA,  # Standardized name
    R0_upper_bound = NA,  # Standardized name
    D_lower_bound = NA,  # Standardized name
    D_upper_bound = NA,  # Standardized name
    nll = NA,  # Kept
    aic = NA,  # Removed in aggregation
    fit_time_sec = NA,  # Removed in aggregation
    n_samples = NA,  # Removed in aggregation
    total_cases = NA  # Removed in aggregation
  )
  
  tryCatch({
    # Load data
    data <- load_case_counts(sim_id, S0)
    observed_cases <- data$daily_cases
    times <- data$times
    init_state <- c(S = data$S0, I = data$I0, R = data$R0_init)
    
    # Remove any days with no data (all zeros at the end)
    last_nonzero <- max(which(observed_cases > 0))
    if (last_nonzero < length(observed_cases)) {
      observed_cases <- observed_cases[1:last_nonzero]
      times <- times[1:last_nonzero]
    }
    
    # Create NLL function using parameterization
    nll <- parameterization$create_nll(sir_ode, init_state, times, observed_cases)
    
    # Compute starting values using parameterization
    start_vals <- parameterization$compute_starting_values(
      strategy = starting_value_strategy,
      true_R0 = true_R0,
      true_recovery_time = true_recovery_time,
      R0_min = R0_min,
      R0_max = R0_max,
      rt_min = rt_min,
      rt_max = rt_max
    )
    
    # Compute bounds using parameterization
    bounds <- parameterization$compute_bounds(
      R0_min = R0_min,
      R0_max = R0_max,
      rt_min = rt_min,
      rt_max = rt_max,
      buffer = 0.1
    )
    
    cat(sprintf("\n--- %s ---\n", sim_id))
    cat(sprintf("Parameterization: %s\n", parameterization$name))
    cat(sprintf("Starting values: %s\n", starting_value_strategy))
    cat(sprintf("True: R0=%.4f, rt=%.4f (beta=%.6f, gamma=%.6f)\n",
                true_R0, true_recovery_time, true_beta, true_gamma))
    cat(sprintf("Start: %s\n", paste(names(start_vals), "=", 
                sprintf("%.6f", unlist(start_vals)), collapse=", ")))
    cat(sprintf("Bounds: [%s] to [%s]\n",
                paste(sprintf("%.4f", bounds$lower), collapse=", "),
                paste(sprintf("%.4f", bounds$upper), collapse=", ")))
    
    # Fit model with bounds
    fit_start_time <- Sys.time()
    fit <- mle2(nll,
                start = start_vals,
                method = "L-BFGS-B",
                lower = bounds$lower,
                upper = bounds$upper,
                control = list(maxit = 1000, trace = 0))
    fit_elapsed <- as.numeric(difftime(Sys.time(), fit_start_time, units = "secs"))
    
    # Transform estimates to standard format (beta, gamma, R0, recovery_time)
    estimates <- parameterization$transform_estimates(coef(fit))
    
    cat(sprintf("Estimated: R0=%.4f, rt=%.4f (beta=%.6f, gamma=%.6f)\n",
                estimates$R0, estimates$recovery_time, estimates$beta, estimates$gamma))
    cat(sprintf("Fit time: %.2f sec\n", fit_elapsed))

    
    # Compute 95% confidence intervals via Monte Carlo sampling from Hessian
    vcov_mat <- vcov(fit)
    n_samples_ci <- 10000
    
    # Sample from multivariate normal in optimization parameter space
    param_samples <- rmvnorm(n_samples_ci, mean = coef(fit), sigma = vcov_mat)
    
    # Transform all samples to (beta, gamma, R0, recovery_time)
    all_estimates <- t(apply(param_samples, 1, function(row) {
      named_row <- setNames(row, parameterization$param_names)
      est <- parameterization$transform_estimates(named_row)
      c(est$beta, est$gamma, est$R0, est$recovery_time)
    }))
    
    beta_samples <- all_estimates[, 1]
    gamma_samples <- all_estimates[, 2]
    R0_samples <- all_estimates[, 3]
    rt_samples <- all_estimates[, 4]
    
    # 95% CI from quantiles
    ci_beta <- quantile(beta_samples, c(0.025, 0.975))
    ci_gamma <- quantile(gamma_samples, c(0.025, 0.975))
    ci_R0 <- quantile(R0_samples, c(0.025, 0.975))
    ci_rt <- quantile(rt_samples, c(0.025, 0.975))
    
    # Store results
    result$converged <- TRUE
    result$est_beta <- estimates$beta
    result$est_gamma <- estimates$gamma
    result$R0_point_estimate <- estimates$R0
    result$D_point_estimate <- estimates$recovery_time
    result$ci_lower_beta <- as.numeric(ci_beta[1])
    result$ci_upper_beta <- as.numeric(ci_beta[2])
    result$ci_lower_gamma <- as.numeric(ci_gamma[1])
    result$ci_upper_gamma <- as.numeric(ci_gamma[2])
    result$R0_lower_bound <- as.numeric(ci_R0[1])
    result$R0_upper_bound <- as.numeric(ci_R0[2])
    result$D_lower_bound <- as.numeric(ci_rt[1])
    result$D_upper_bound <- as.numeric(ci_rt[2])
    result$nll <- as.numeric(-logLik(fit))
    result$aic <- AIC(fit)
    result$fit_time_sec <- fit_elapsed
    result$n_samples <- length(observed_cases)
    result$total_cases <- sum(observed_cases)
    
  }, error = function(e) {
    cat(sprintf("ERROR: %s\n", e$message))
    result$error_msg <<- e$message
  })
  
  return(result)
}


# ============================================================================
# Main Processing
# ============================================================================

# Load true parameters
params_df <- read.csv(params_file)

# Initialize results list
results_list <- list()

# Process each simulation
for (i in seq_along(sim_ids)) {
  sim_id <- sim_ids[i]
  
  # Parse sim_id to get param_id
  match <- regexec("(train|test)_p(\\d{5})_r(\\d+)", sim_id)
  matches <- regmatches(sim_id, match)[[1]]
  param_id <- as.integer(matches[3])
  replicate <- as.integer(matches[4])
  
  # Get true parameters
  param_row <- params_df[params_df$param_id == param_id, ]
  if (nrow(param_row) == 0) {
    cat(sprintf("WARNING: No parameters found for param_id=%d\n", param_id))
    next
  }
  
  true_beta <- param_row$beta[1]
  true_gamma <- param_row$gamma[1]
  true_R0 <- param_row$R0[1]
  true_recovery_time <- param_row$recovery_time[1]
  S0 <- param_row$S[1]
  
  # Fit model with specified parameterization and starting value strategy
  result <- fit_sir_mle(sim_id, true_beta, true_gamma, true_R0, true_recovery_time,
                       starting_value_strategy, S0, parameterization)
  result$param_id <- param_id
  result$replicate <- replicate
  
  results_list[[i]] <- result
}

# Convert to data frame
results_df <- do.call(rbind, lapply(results_list, as.data.frame))

# Reorder columns
col_order <- c(
  "sim_id", "param_id", "replicate", "converged", "error_msg",
  "starting_value_strategy",
  "true_beta", "true_gamma", "R0_true", "D_true",
  "est_beta", "est_gamma", "R0_point_estimate", "D_point_estimate",
  "ci_lower_beta", "ci_upper_beta", 
  "ci_lower_gamma", "ci_upper_gamma",
  "R0_lower_bound", "R0_upper_bound",
  "D_lower_bound", "D_upper_bound",
  "nll", "aic", "fit_time_sec", "n_samples", "total_cases"
)
results_df <- results_df[, col_order[col_order %in% names(results_df)]]

# Save results
write.csv(results_df, output_file, row.names = FALSE)

# Print summary
cat("\n========================================\n")
cat("Batch Summary\n")
cat("========================================\n")
cat(sprintf("Composite strategy: %s\n", composite_strategy))
cat(sprintf("  Parameterization: %s\n", parameterization_name))
cat(sprintf("  Starting values: %s\n", starting_value_strategy))
cat(sprintf("Total: %d\n", nrow(results_df)))
cat(sprintf("Converged: %d (%.1f%%)\n", 
            sum(results_df$converged, na.rm = TRUE),
            100 * mean(results_df$converged, na.rm = TRUE)))
cat(sprintf("\nOutput: %s\n", output_file))

# Close log
sink(type = "message")
sink(type = "output")
close(log_conn)

