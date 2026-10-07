#!/usr/bin/env Rscript

#' Extract tsfeatures in batch with parallel processing (Snakemake version)
#'
#' Processes multiple BINNED case count files in parallel using mclapply.
#' Each series takes ~0.04 seconds for 1000 time points.
#' With 4 cores, 1000 series take ~10 seconds.
#'
#' This script extracts 55 time series features from BINNED daily case count trajectories
#' using fast, fixed-length tsfeatures functions.
#'
#' Note: This reads BINNED case counts (daily aggregated) from *_binned_case_counts.csv,
#' NOT the raw event-level data from *_case_counts.csv.
#'
#' Note: hw_parameters and station_features are excluded as they return NaN
#' on non-seasonal epidemic data.

# Load required packages
library(tsfeatures)
library(parallel)

# Source the optimized calc_fast_tsfeatures function
source("workflow/scripts/utils/calc_fast_tsfeatures.R")

#' Extract tsfeatures for a single simulation
#' Reads case count CSV and calculates time series features
extract_tsfeatures_one_file <- function(sim_id) {
  tryCatch({
    # Input: BINNED case counts CSV from phase2 step 4b
    case_counts_file <- paste0("results/phase2_processing/case_counts/", sim_id, "_binned_case_counts.csv")
    output_stats <- paste0("results/phase2_processing/features/temp/tsfeatures/", sim_id, "_tsfeatures.csv")
    
    # Create output directory
    dir.create(dirname(output_stats), recursive = TRUE, showWarnings = FALSE)
    
    # Read BINNED case counts CSV
    df <- read.csv(case_counts_file)
    
    # Extract case_count column (the binned daily case count time series)
    case_counts <- df$case_count
    
    # Calculate tsfeatures using optimized function
    # Use frequency=1 for daily data (non-seasonal)
    features <- calc_fast_tsfeatures(case_counts, frequency = 1)
    
    # Convert to data frame (single row with multiple columns)
    features_df <- as.data.frame(t(features))
    
    # Save to CSV
    write.csv(features_df, output_stats, row.names = FALSE)
    
    return(list(sim_id = sim_id, success = TRUE, error = NULL))
    
  }, error = function(e) {
    error_msg <- paste0(class(e)[1], ": ", e$message)
    return(list(sim_id = sim_id, success = FALSE, error = error_msg))
  })
}

# Main execution
if (exists("snakemake")) {
  # Extract parameters from Snakemake
  sim_ids <- snakemake@params$sim_ids
  batch_id <- snakemake@wildcards$batch_id
  log_file <- snakemake@log[[1]]
  output_marker <- snakemake@output[[1]]
  n_cores <- snakemake@threads
  
  cat("tsfeatures Extraction - Batch", batch_id, "\n")
  cat("Processing", length(sim_ids), "files with", n_cores, "cores\n")
  cat("Output: 49 time series features per simulation\n")
  cat("(FFORMA arch/garch features already in heterogeneity())\n")
  cat("(hw_parameters and station_features excluded - return NaN on epidemic data)\n")
  cat("(stl_features_nperiods and stl_features_seasonal_period excluded - constant for non-seasonal data)\n\n")
  
  # Process files in parallel
  results <- mclapply(
    sim_ids, 
    extract_tsfeatures_one_file,
    mc.cores = n_cores,
    mc.preschedule = TRUE  # Fast uniform processing
  )
  
  # Analyze results
  successful <- sapply(results, function(x) x$success)
  failed_results <- results[!successful]
  
  n_success <- sum(successful)
  n_failed <- sum(!successful)
  
  # Write detailed log
  log_dir <- dirname(log_file)
  dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
  
  sink(log_file)
  cat("tsfeatures Extraction - Batch", batch_id, "\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat("Configuration:\n")
  cat("  Frequency: 1 (daily data)\n")
  cat("  Expected features: 49\n")
  cat("  Note: FFORMA arch/garch features already computed by heterogeneity()\n")
  cat("  (hw_parameters and station_features excluded - return NaN on epidemic data)\n")
  cat("  (stl_features_nperiods and stl_features_seasonal_period excluded - constant for non-seasonal data)\n\n")
  cat("Total files:", length(sim_ids), "\n")
  cat("Successful:", n_success, "\n")
  cat("Failed:", n_failed, "\n")
  cat("Cores used:", n_cores, "\n\n")
  
  if (n_failed > 0) {
    cat("Failed extractions:\n")
    for (result in failed_results) {
      cat("  -", result$sim_id, ":", result$error, "\n")
    }
    
    # Write failures to separate file
    failure_file <- file.path(log_dir, paste0("batch_", batch_id, "_tsfeatures_failures.txt"))
    writeLines(sapply(failed_results, function(x) x$sim_id), failure_file)
    cat("\nFailure list saved to:", failure_file, "\n")
  }
  
  sink()
  
  # Print summary to stdout
  cat("Batch", batch_id, "complete:\n")
  cat("  Successful:", n_success, "/", length(sim_ids), "\n")
  cat("  Failed:", n_failed, "\n")
  
  # Create output marker file (touched by Snakemake)
  # The marker file is created by Snakemake's touch() directive
  
} else {
  stop("This script must be run via Snakemake")
}
