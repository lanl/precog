#!/usr/bin/env Rscript

#' Extract tree features in batch with parallel processing (Snakemake version)
#' 
#' Processes multiple tree files in parallel using mclapply.
#' Each tree takes ~5-30 seconds depending on size.
#' With 4 cores, 200 trees take ~10-15 minutes.
#' 
#' This script extracts ALL 51 tree features (size_agnostic=FALSE).
#' It excludes LTT and NTT files - only uses the aggregated tree features CSV.
#' 
#' Note: treestats and dependencies are pre-installed by setup_r_treestats rule
#' to avoid race conditions during parallel batch processing.

# Load required packages (already installed via conda + CRAN setup)
library(ape)
library(parallel)

# Verify treestats is available
if (!requireNamespace("treestats", quietly = TRUE)) {
  stop("treestats package not found. The setup_r_treestats rule should have installed it.")
}
library(treestats)

# Configuration for tree features extraction
# Use ALL 51 features (size_agnostic=FALSE)
size_agnostic <- FALSE
normalization <- "both"  # Ignored when size_agnostic=FALSE

# Source the optimized calc_fast_stats function
source("workflow/scripts/utils/calc_fast_stats.R")

#' Extract tree features for a single simulation
#' Each simulation produces exactly one tree in the .trees file
extract_treefeatures_one_file <- function(sim_id, size_agnostic, normalization) {
  tryCatch({
    trees_file <- paste0("results/phase1_simulation/simulations/", sim_id, ".trees")
    output_stats <- paste0("results/phase2_processing/features/temp/treefeatures/", sim_id, "_treefeatures.csv")
    
    # Create output directory
    dir.create(dirname(output_stats), recursive = TRUE, showWarnings = FALSE)
    
    # Read tree file (BEAST format - Nexus)
    # Each simulation produces exactly one tree
    tree <- read.nexus(trees_file)
    
    # Calculate tree features using optimized function
    stats <- calc_fast_stats(
      tree, 
      size_agnostic = size_agnostic,
      normalization = normalization
    )
    
    # Convert to data frame (single row with multiple columns)
    # stats is a named vector, need to convert to single-row dataframe
    stats_df <- as.data.frame(t(stats))
    
    # Save to CSV
    write.csv(stats_df, output_stats, row.names = FALSE)
    
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
  
  cat("Tree Features Extraction - Batch", batch_id, "\n")
  cat("Processing", length(sim_ids), "files with", n_cores, "cores\n")
  cat("Settings: size_agnostic =", size_agnostic, ", normalization =", normalization, "\n")
  cat("Output: ALL 51 tree features per simulation\n\n")
  
  # Process files in parallel
  results <- mclapply(
    sim_ids, 
    extract_treefeatures_one_file,
    size_agnostic = size_agnostic,
    normalization = normalization,
    mc.cores = n_cores,
    mc.preschedule = FALSE  # Better load balancing for variable tree sizes
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
  cat("Tree Features Extraction - Batch", batch_id, "\n")
  cat(paste(rep("=", 80), collapse = ""), "\n")
  cat("Configuration:\n")
  cat("  size_agnostic =", size_agnostic, "\n")
  cat("  normalization =", normalization, "\n")
  cat("  Expected features: 51\n\n")
  cat("Total files:", length(sim_ids), "\n")
  cat("Successful:", n_success, "\n")
  cat("Failed:", n_failed, "\n")
  cat("Cores used:", n_cores, "\n\n")
  
  # No additional summary needed - one tree per simulation
  
  if (n_failed > 0) {
    cat("Failed extractions:\n")
    for (result in failed_results) {
      cat("  -", result$sim_id, ":", result$error, "\n")
    }
    
    # Write failures to separate file
    failure_file <- file.path(log_dir, paste0("batch_", batch_id, "_treefeatures_failures.txt"))
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

