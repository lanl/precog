#!/usr/bin/env Rscript

#' Benchmark tsfeatures functions with intelligent timeouts and variable-length detection
#' Identifies fast, fixed-length features for case count time series
#' 
#' Features:
#' - Timeout mechanism to skip extremely slow features
#' - Detection of variable-length outputs (excluded from fast features)
#' - Progressive benchmarking from short to long time series
#' - Skip slow features on longer series to save time
#' - **Per-feature time threshold**: Fairly evaluates multi-feature functions

library(tsfeatures)

cat("=== tsfeatures Benchmarking Suite ===\n\n")

# Configuration
SERIES_LENGTHS <- c(100, 500, 1000)  # Progressive lengths
TIMEOUT_SECONDS <- 30  # Max time per feature
SLOW_THRESHOLD <- 1.0   # Features slower than this won't be tested on longer series
TIME_THRESHOLD_FAST <- 0.1  # Classify as fast (PER FEATURE)
TIME_THRESHOLD_SLOW <- 1.0  # Classify as slow (PER FEATURE)
N_REPS <- 3
FREQUENCY <- 1  # Daily data (non-seasonal)

# Track slow features to skip on longer series
slow_features <- character(0)

# Generate synthetic epidemic-like time series
generate_epidemic_ts <- function(n, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  
  # Generate epidemic curve (gamma-like distribution)
  # Simulate realistic case count trajectory
  t <- seq(0, 1, length.out = n)
  
  # Peak somewhere in middle
  peak_time <- runif(1, 0.3, 0.7)
  
  # Generate smooth epidemic curve
  curve <- dgamma((t - peak_time + 0.5) * 20, shape = 3, rate = 1) * runif(1, 500, 2000)
  
  # Add some noise
  noise <- rnorm(n, 0, max(curve) * 0.05)
  
  # Case counts should be non-negative integers
  cases <- pmax(0, round(curve + noise))
  
  # Convert to ts object
  ts(cases, frequency = FREQUENCY)
}

# Get all available tsfeatures functions
# Based on package documentation and exports
get_all_tsfeatures <- function() {
  list(
    # Original 18 functions
    acf_features = tsfeatures::acf_features,
    arch_stat = tsfeatures::arch_stat,
    crossing_points = tsfeatures::crossing_points,
    entropy = tsfeatures::entropy,
    flat_spots = tsfeatures::flat_spots,
    heterogeneity = tsfeatures::heterogeneity,
    holt_parameters = tsfeatures::holt_parameters,
    hw_parameters = tsfeatures::hw_parameters,
    lumpiness = tsfeatures::lumpiness,
    max_kl_shift = tsfeatures::max_kl_shift,
    max_level_shift = tsfeatures::max_level_shift,
    max_var_shift = tsfeatures::max_var_shift,
    nonlinearity = tsfeatures::nonlinearity,
    pacf_features = tsfeatures::pacf_features,
    stl_features = tsfeatures::stl_features,
    stability = tsfeatures::stability,
    unitroot_kpss = tsfeatures::unitroot_kpss,
    unitroot_pp = tsfeatures::unitroot_pp,
    
    # Previously missed functions (5 basic ones)
    autocorr_features = tsfeatures::autocorr_features,
    binarize_mean = tsfeatures::binarize_mean,
    firstzero_ac = tsfeatures::firstzero_ac,
    hurst = tsfeatures::hurst,
    zero_proportion = tsfeatures::zero_proportion,
    
    # Additional feature sets
    compengine = tsfeatures::compengine,
    dist_features = tsfeatures::dist_features,
    pred_features = tsfeatures::pred_features,
    scal_features = tsfeatures::scal_features,
    station_features = tsfeatures::station_features
  )
}

# Safely execute with timeout
safe_execute_with_timeout <- function(expr, timeout_sec) {
  tryCatch({
    result <- NULL
    timed_out <- FALSE
    
    # Use setTimeLimit for timeout
    setTimeLimit(cpu = timeout_sec, elapsed = timeout_sec, transient = TRUE)
    
    result <- expr
    
    # Reset time limit
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
    
    list(result = result, timed_out = FALSE, error = NULL)
  }, error = function(e) {
    # Reset time limit on error
    setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
    
    if (grepl("reached elapsed time limit|reached CPU time limit", e$message)) {
      list(result = NULL, timed_out = TRUE, error = "TIMEOUT")
    } else {
      list(result = NULL, timed_out = FALSE, error = as.character(e$message))
    }
  })
}

# Check if feature returns fixed-length output
check_fixed_length <- function(feature_func, feature_name) {
  tryCatch({
    # Test on three different lengths
    test_lengths <- c(50, 100, 200)
    output_lengths <- sapply(test_lengths, function(n) {
      ts_data <- generate_epidemic_ts(n, seed = 42)
      result <- feature_func(ts_data)
      length(unlist(result))
    })
    
    # Check if all outputs have same length
    is_fixed <- length(unique(output_lengths)) == 1
    
    list(
      is_fixed_length = is_fixed,
      output_length = output_lengths[1],
      test_lengths = test_lengths,
      output_lengths = output_lengths
    )
  }, error = function(e) {
    list(
      is_fixed_length = FALSE,
      output_length = NA,
      test_lengths = test_lengths,
      output_lengths = rep(NA, length(test_lengths)),
      error = as.character(e$message)
    )
  })
}


# Benchmark a single feature
benchmark_feature <- function(feature_name, feature_func, series_length, n_reps) {
  
  # Skip if marked as slow from previous iteration
  if (feature_name %in% slow_features) {
    return(data.frame(
      feature = feature_name,
      series_length = series_length,
      mean_time = NA,
      sd_time = NA,
      n_reps = 0,
      timed_out = FALSE,
      skipped = TRUE,
      error = NA,
      is_fixed_length = NA,
      output_length = NA
    ))
  }
  
  cat(sprintf("  Testing %s (n=%d)... ", feature_name, series_length))
  
  # First check if it's fixed-length (only on first size)
  is_fixed_length <- TRUE
  output_length <- NA
  if (series_length == SERIES_LENGTHS[1]) {
    length_check <- check_fixed_length(feature_func, feature_name)
    is_fixed_length <- length_check$is_fixed_length
    output_length <- length_check$output_length
    
    if (!is_fixed_length) {
      cat(sprintf("VARIABLE-LENGTH (lengths: %s)\n", 
                  paste(length_check$output_lengths, collapse = ", ")))
      return(data.frame(
        feature = feature_name,
        series_length = series_length,
        mean_time = NA,
        sd_time = NA,
        n_reps = 0,
        timed_out = FALSE,
        skipped = FALSE,
        error = "VARIABLE_LENGTH",
        is_fixed_length = FALSE,
        output_length = NA
      ))
    }
  }
  
  # Generate test time series
  ts_data <- generate_epidemic_ts(series_length, seed = 12345)
  
  # Try with timeout first
  test_result <- safe_execute_with_timeout({
    feature_func(ts_data)
  }, TIMEOUT_SECONDS)
  
  if (test_result$timed_out) {
    cat("TIMEOUT\n")
    return(data.frame(
      feature = feature_name,
      series_length = series_length,
      mean_time = NA,
      sd_time = NA,
      n_reps = 0,
      timed_out = TRUE,
      skipped = FALSE,
      error = "TIMEOUT",
      is_fixed_length = is_fixed_length,
      output_length = output_length
    ))
  }
  
  if (!is.null(test_result$error)) {
    cat(sprintf("ERROR: %s\n", test_result$error))
    return(data.frame(
      feature = feature_name,
      series_length = series_length,
      mean_time = NA,
      sd_time = NA,
      n_reps = 0,
      timed_out = FALSE,
      skipped = FALSE,
      error = test_result$error,
      is_fixed_length = is_fixed_length,
      output_length = output_length
    ))
  }
  
  # Feature works - now benchmark it
  times <- numeric(n_reps)
  for (i in 1:n_reps) {
    # Generate new series each time
    ts_data <- generate_epidemic_ts(series_length, seed = 12345 + i)
    
    start_time <- Sys.time()
    result <- feature_func(ts_data)
    end_time <- Sys.time()
    
    times[i] <- as.numeric(difftime(end_time, start_time, units = "secs"))
  }
  
  mean_time <- mean(times)
  sd_time <- sd(times)
  
  cat(sprintf("%.4f sec (±%.4f)\n", mean_time, sd_time))
  
  # Mark as slow if over threshold
  if (mean_time > SLOW_THRESHOLD) {
    slow_features <<- c(slow_features, feature_name)
  }
  
  return(data.frame(
    feature = feature_name,
    series_length = series_length,
    mean_time = mean_time,
    sd_time = sd_time,
    n_reps = n_reps,
    timed_out = FALSE,
    skipped = FALSE,
    error = NA,
    is_fixed_length = is_fixed_length,
    output_length = output_length
  ))
}

# Main benchmarking loop
cat("Configuration:\n")
cat(sprintf("  Series lengths: %s\n", paste(SERIES_LENGTHS, collapse = ", ")))
cat(sprintf("  Timeout: %d seconds\n", TIMEOUT_SECONDS))
cat(sprintf("  Fast threshold: %.2f seconds\n", TIME_THRESHOLD_FAST))
cat(sprintf("  Slow threshold: %.2f seconds\n", TIME_THRESHOLD_SLOW))
cat(sprintf("  Replicates: %d\n", N_REPS))
cat(sprintf("  Frequency: %d (daily)\n\n", FREQUENCY))

# Get all features
all_features <- get_all_tsfeatures()
cat(sprintf("Total features to test: %d\n\n", length(all_features)))

# Initialize results
results <- list()

# Progressive benchmarking
for (series_length in SERIES_LENGTHS) {
  cat(paste(rep("=", 70), collapse = ""), "\n")
  cat(sprintf("Testing with series length: %d time points\n", series_length))
  cat(paste(rep("=", 70), collapse = ""), "\n\n")
  
  for (feature_name in names(all_features)) {
    result <- benchmark_feature(
      feature_name = feature_name,
      feature_func = all_features[[feature_name]],
      series_length = series_length,
      n_reps = N_REPS
    )
    results[[length(results) + 1]] <- result
  }
  
  cat("\n")
}

# Combine results
results_df <- do.call(rbind, results)

# Save results
output_dir <- "results"
if (!dir.exists(output_dir)) {
  dir.create(output_dir, recursive = TRUE)
}

write.csv(results_df, file.path(output_dir, "benchmark_results.csv"), row.names = FALSE)


# Analyze results for longest series
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("ANALYSIS: Fast, Fixed-Length Features\n")
cat(paste(rep("=", 70), collapse = ""), "\n\n")

largest_size <- max(SERIES_LENGTHS)
large_series_results <- results_df[results_df$series_length == largest_size & 
                                   !is.na(results_df$mean_time), ]

# Filter for fixed-length features
fixed_length_results <- large_series_results[
  large_series_results$is_fixed_length == TRUE | 
  is.na(large_series_results$is_fixed_length), 
]

# Get output_length from 100-point results (since 1000 has NA)
results_100 <- results_df[results_df$series_length == 100, ]
output_lengths <- sapply(fixed_length_results$feature, function(f) {
  row <- results_100[results_100$feature == f, ]
  if (nrow(row) > 0 && !is.na(row$output_length[1])) {
    return(row$output_length[1])
  }
  return(NA)
})

# Add output_length and calculate per-feature time
fixed_length_results$output_length_actual <- output_lengths
fixed_length_results$time_per_feature <- ifelse(
  !is.na(fixed_length_results$output_length_actual) & fixed_length_results$output_length_actual > 0,
  fixed_length_results$mean_time / fixed_length_results$output_length_actual,
  NA
)

# Classify by per-feature speed (primary metric)
fast_features <- fixed_length_results[!is.na(fixed_length_results$time_per_feature) & 
                                       fixed_length_results$time_per_feature < TIME_THRESHOLD_FAST, ]
moderate_features <- fixed_length_results[!is.na(fixed_length_results$time_per_feature) &
                                           fixed_length_results$time_per_feature >= TIME_THRESHOLD_FAST & 
                                           fixed_length_results$time_per_feature < TIME_THRESHOLD_SLOW, ]
slow_features_analyzed <- fixed_length_results[!is.na(fixed_length_results$time_per_feature) &
                                                fixed_length_results$time_per_feature >= TIME_THRESHOLD_SLOW, ]

# Count problem features
timeout_count <- sum(results_df$series_length == largest_size & results_df$timed_out, na.rm = TRUE)
skipped_count <- sum(results_df$series_length == largest_size & results_df$skipped, na.rm = TRUE)
variable_length_count <- sum(results_df$series_length == largest_size & 
                              results_df$error == "VARIABLE_LENGTH", na.rm = TRUE)

cat(sprintf("Results for series with %d time points:\n", largest_size))
cat(sprintf("  Fast features (< %.1fs PER FEATURE): %d\n", TIME_THRESHOLD_FAST, nrow(fast_features)))
cat(sprintf("  Moderate features (%.1f-%.1fs PER FEATURE): %d\n", 
            TIME_THRESHOLD_FAST, TIME_THRESHOLD_SLOW, nrow(moderate_features)))
cat(sprintf("  Slow features (> %.1fs PER FEATURE): %d\n", TIME_THRESHOLD_SLOW, nrow(slow_features_analyzed)))
cat(sprintf("  Variable-length features (excluded): %d\n", variable_length_count))
cat(sprintf("  Timed out features (>%ds): %d\n", TIMEOUT_SECONDS, timeout_count))
cat(sprintf("  Skipped features (slow on shorter series): %d\n", skipped_count))

# Calculate time savings
if (nrow(fast_features) > 0) {
  total_all_stats <- sum(large_series_results$mean_time, na.rm = TRUE)
  total_fast_only <- sum(fast_features$mean_time, na.rm = TRUE)
  speedup <- if (total_fast_only > 0) total_all_stats / total_fast_only else 1
  
  cat(sprintf("\n--- Time Estimates for %d time points ---\n", largest_size))
  cat(sprintf("Time to compute ALL successful features: %.2f seconds\n", total_all_stats))
  cat(sprintf("Time to compute FAST features only: %.2f seconds\n", total_fast_only))
  cat(sprintf("Speedup: %.1fx faster\n", speedup))
  cat(sprintf("Time saved per series: %.2f seconds\n", total_all_stats - total_fast_only))
  
  # Project for large datasets
  cat(sprintf("\n--- Projections for Large Datasets ---\n"))
  for (n_series in c(10000, 20000, 50000)) {
    all_time_hrs <- (total_all_stats * n_series) / 3600
    fast_time_hrs <- (total_fast_only * n_series) / 3600
    cat(sprintf("For %d series:\n", n_series))
    cat(sprintf("  All features: %.1f hours\n", all_time_hrs))
    cat(sprintf("  Fast only: %.1f hours (%.1fx faster)\n", fast_time_hrs, speedup))
  }
}

# Save fast features list
if (nrow(fast_features) > 0) {
  fast_feature_names <- fast_features$feature
  save(fast_feature_names, file = file.path(output_dir, "fast_features.RData"))
  writeLines(fast_feature_names, file.path(output_dir, "fast_features.txt"))
  cat(sprintf("\n✓ Fast features list saved to: %s\n", 
              file.path(output_dir, "fast_features.txt")))
  
  # Show detailed lists
  cat(sprintf("\n--- Fast Features (fixed-length, < %.1fs PER FEATURE) ---\n", TIME_THRESHOLD_FAST))
  fast_features_sorted <- fast_features[order(fast_features$time_per_feature), ]
  for (i in 1:nrow(fast_features_sorted)) {
    cat(sprintf("%2d. %s: %.4f sec total / %s features = %.4f sec/feature\n", 
                i, 
                fast_features_sorted$feature[i], 
                fast_features_sorted$mean_time[i],
                if(is.na(fast_features_sorted$output_length_actual[i])) "?" else as.character(fast_features_sorted$output_length_actual[i]),
                fast_features_sorted$time_per_feature[i]))
  }
}

# Save variable-length features list
variable_length_features <- results_df[results_df$error == "VARIABLE_LENGTH" & 
                                       !is.na(results_df$error), ]
if (nrow(variable_length_features) > 0) {
  variable_length_names <- unique(variable_length_features$feature)
  writeLines(variable_length_names, file.path(output_dir, "variable_length_features.txt"))
  cat(sprintf("\n--- Variable-Length Features (excluded) ---\n"))
  for (fname in variable_length_names) {
    cat(sprintf("  %s\n", fname))
  }
}

# Show problematic features
if (nrow(slow_features_analyzed) > 0 || timeout_count > 0) {
  cat(sprintf("\n--- Slow/Problematic Features (AVOID) ---\n"))
  
  # Show timed out ones
  timeout_features <- results_df[results_df$series_length == largest_size & 
                                 results_df$timed_out, ]
  if (nrow(timeout_features) > 0) {
    for (i in 1:nrow(timeout_features)) {
      cat(sprintf("  %s: TIMEOUT (>%ds)\n", timeout_features$feature[i], TIMEOUT_SECONDS))
    }
  }
  
  # Show slow but completed
  if (nrow(slow_features_analyzed) > 0) {
    slow_sorted <- slow_features_analyzed[order(-slow_features_analyzed$mean_time), ]
    for (i in 1:nrow(slow_sorted)) {
      cat(sprintf("  %s: %.3f sec\n", slow_sorted$feature[i], slow_sorted$mean_time[i]))
    }
  }
}

cat("\n=== Benchmark complete! ===\n")

cat(sprintf("\n✓ Results saved to: %s\n\n", file.path(output_dir, "benchmark_results.csv")))

    output_length = output_length
  ))
}

}
