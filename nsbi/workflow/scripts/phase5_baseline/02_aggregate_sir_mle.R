#!/usr/bin/env Rscript
# Aggregate SIR MLE Batch Results

suppressPackageStartupMessages({ library(yaml) })

cat("\n========================================\n")
cat("Aggregate SIR MLE Results\n")
cat("========================================\n\n")

# Get snakemake inputs
batch_files <- snakemake@input$batch_files
params_file <- snakemake@input$params
output_results <- snakemake@output$results
output_summary <- snakemake@output$summary
expected_sims <- snakemake@params$expected_sims

cat(sprintf("Loading %d batch files...\n", length(batch_files)))
batch_list <- lapply(batch_files, read.csv)
all_results <- do.call(rbind, batch_list)

cat(sprintf("Total: %d (expected: %d)\n", nrow(all_results), expected_sims))

# Add metadata columns
all_results$method <- "mle"
alpha <- 0.05  # 95% CI
all_results$alpha <- alpha
all_results$nominal_coverage <- 1 - alpha

# Compute diagnostic columns for ALL rows
# For non-converged rows, these will be NA (semantically correct - diagnostics are undefined)
all_results$R0_error <- ifelse(all_results$converged, 
                                all_results$R0_point_estimate - all_results$R0_true, 
                                NA)
all_results$D_error <- ifelse(all_results$converged,
                                all_results$D_point_estimate - all_results$D_true,
                                NA)
all_results$R0_abs_error <- ifelse(all_results$converged,
                                    abs(all_results$R0_error),
                                    NA)
all_results$D_abs_error <- ifelse(all_results$converged,
                                    abs(all_results$D_error),
                                    NA)

# Coverage (logical columns)
all_results$R0_coverage <- ifelse(all_results$converged,
                                   (all_results$R0_true >= all_results$R0_lower_bound) & 
                                   (all_results$R0_true <= all_results$R0_upper_bound),
                                   NA)
all_results$D_coverage <- ifelse(all_results$converged,
                                   (all_results$D_true >= all_results$D_lower_bound) & 
                                   (all_results$D_true <= all_results$D_upper_bound),
                                   NA)

# Interval widths
all_results$R0_interval_width <- ifelse(all_results$converged,
                                         all_results$R0_upper_bound - all_results$R0_lower_bound,
                                         NA)
all_results$D_interval_width <- ifelse(all_results$converged,
                                         all_results$D_upper_bound - all_results$D_lower_bound,
                                         NA)

# Interval Scores (proper scoring rule for confidence intervals)
# IS = width + (2/alpha) * penalties_if_not_covered
# Lower is better; balances sharpness vs coverage
all_results$R0_interval_score <- ifelse(all_results$converged,
  all_results$R0_interval_width + 
    (2.0 / alpha) * pmax(0, all_results$R0_lower_bound - all_results$R0_true) +
    (2.0 / alpha) * pmax(0, all_results$R0_true - all_results$R0_upper_bound),
  NA)

all_results$D_interval_score <- ifelse(all_results$converged,
  all_results$D_interval_width +
    (2.0 / alpha) * pmax(0, all_results$D_lower_bound - all_results$D_true) +
    (2.0 / alpha) * pmax(0, all_results$D_true - all_results$D_upper_bound),
  NA)

# Remove columns not needed for comparison
all_results$starting_value_strategy <- NULL
all_results$true_beta <- NULL
all_results$true_gamma <- NULL
all_results$est_beta <- NULL
all_results$est_gamma <- NULL
all_results$ci_lower_beta <- NULL
all_results$ci_upper_beta <- NULL
all_results$ci_lower_gamma <- NULL
all_results$ci_upper_gamma <- NULL
all_results$aic <- NULL
all_results$fit_time_sec <- NULL
all_results$n_samples <- NULL
all_results$total_cases <- NULL

# Reorder columns (27 total)
col_order <- c(
  "sim_id", "param_id", "replicate", "converged", "error_msg",
  "method", "alpha", "nominal_coverage",
  "R0_true", "D_true",
  "R0_point_estimate", "D_point_estimate",
  "R0_lower_bound", "R0_upper_bound", "D_lower_bound", "D_upper_bound",
  "R0_interval_width", "D_interval_width",
  "R0_coverage", "D_coverage",
  "R0_error", "D_error",
  "R0_abs_error", "D_abs_error",
  "R0_interval_score", "D_interval_score",
  "nll"
)
all_results <- all_results[, col_order[col_order %in% names(all_results)]]

# Save consolidated results (27 columns total)
write.csv(all_results, output_results, row.names = FALSE)

# Calculate summary statistics (only for converged rows)
n_converged <- sum(all_results$converged, na.rm = TRUE)
convergence_rate <- n_converged / nrow(all_results)

cat(sprintf("Converged: %d (%.1f%%)\n", n_converged, 100*convergence_rate))

# Summary statistics
summary_stats <- list(
  n_simulations = nrow(all_results),
  n_converged = n_converged,
  convergence_rate = convergence_rate,
  method = "mle",
  alpha = alpha,
  nominal_coverage = 1 - alpha
)

if (n_converged > 0) {
  # Use only converged rows for summary metrics (na.rm handles filtering)
  summary_stats$R0 <- list(
    mae = mean(all_results$R0_abs_error, na.rm = TRUE),
    rmse = sqrt(mean(all_results$R0_error^2, na.rm = TRUE)),
    bias = mean(all_results$R0_error, na.rm = TRUE),
    coverage_95 = mean(all_results$R0_coverage, na.rm = TRUE),
    mean_interval_width = mean(all_results$R0_interval_width, na.rm = TRUE),
    mean_interval_score = mean(all_results$R0_interval_score, na.rm = TRUE)
  )
  
  summary_stats$D <- list(
    mae = mean(all_results$D_abs_error, na.rm = TRUE),
    rmse = sqrt(mean(all_results$D_error^2, na.rm = TRUE)),
    bias = mean(all_results$D_error, na.rm = TRUE),
    coverage_95 = mean(all_results$D_coverage, na.rm = TRUE),
    mean_interval_width = mean(all_results$D_interval_width, na.rm = TRUE),
    mean_interval_score = mean(all_results$D_interval_score, na.rm = TRUE)
  )
}

write_yaml(summary_stats, output_summary)
cat("\n✓ Aggregation complete!\n")
cat(sprintf("Output: %s\n", output_results))

