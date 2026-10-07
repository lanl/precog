#!/usr/bin/env Rscript

#' Calculate Fast Time Series Features
#' 
#' Optimized function that computes only fast, fixed-length tsfeatures
#' (< 0.1 seconds for 1000 time points).
#' 
#' Based on comprehensive benchmarking:
#' - 25 fast tsfeatures functions (< 0.1s on 1000 points)
#' - All return fixed-length output (independent of time series length)
#' - Total computation time: ~0.07s per 1000-point series
#' 
#' Features included (25 functions → ~49 scalar features):
#' - acf_features: 6 ACF-based features
#' - arch_stat: 1 ARCH LM statistic  
#' - crossing_points: 1 crossing point count
#' - entropy: 1 spectral entropy
#' - flat_spots: 1 flat spot indicator
#' - heterogeneity: 4 heterogeneity measures (includes arch_acf, garch_acf, arch_r2, garch_r2)
#' - holt_parameters: 2 Holt's linear trend parameters
#' - lumpiness: 1 lumpiness measure
#' - max_kl_shift: 2 maximum KL shift values
#' - max_level_shift: 2 maximum level shift values
#' - max_var_shift: 2 maximum variance shift values
#' - nonlinearity: 1 nonlinearity measure
#' - pacf_features: 3 PACF-based features
#' - stl_features: 6 STL decomposition features (nperiods and seasonal_period excluded - constant for non-seasonal data)
#' - stability: 1 stability measure
#' - unitroot_kpss: 1 KPSS unit root test statistic
#' - unitroot_pp: 1 Phillips-Perron test statistic
#' - autocorr_features: 7 autocorrelation features
#' - firstzero_ac: 1 first zero crossing of ACF
#' - hurst: 1 Hurst exponent
#' - zero_proportion: 1 proportion of zeros
#' - dist_features: 2 distribution features
#' - total_case_count: 1 total outbreak size (sum of all cases)
#' 
#' Note: total_case_count is equivalent to tree statistic "number_of_lineages"
#' (total tips in phylogeny), representing cumulative outbreak size.
#' 
#' EXCLUDED (return NaN on epidemic data):
#' - hw_parameters: Requires seasonal data (epidemic curves are non-seasonal)
#' - station_features: spreadrandomlocal_meantaul_50 fails on non-stationary data
#' 
#' Excluded (slow or variable-length):
#' - binarize_mean: variable-length output
#' - compengine: 1.9s (too slow, 16 features)
#' - pred_features: 1.3s (too slow, 3 features)
#' - scal_features: 0.5s (moderate, 1 feature)
#' 
#' Total: ~49 scalar features (after removing 2 constant STL features)
#' 
#' Performance:
#' - ~0.07s per 1000-point series
#' - For 20,000 series: ~0.39 hours total
#' 
#' @param case_counts Numeric vector of case counts (can include zeros)
#' @param frequency Time series frequency (default: 1 for daily data)
#' @return Named numeric vector of time series features
#' @export

library(tsfeatures)

calc_fast_tsfeatures <- function(case_counts, frequency = 1) {
  
  # Helper function to safely compute a feature
  try_feature <- function(feature_name, func_call) {
    tryCatch({
      result <- func_call
      # If result is a list/data.frame, flatten to named vector
      if (is.list(result) || is.data.frame(result)) {
        result <- unlist(result)
      }
      # Add feature name prefix to avoid collisions
      if (!is.null(names(result))) {
        names(result) <- paste0(feature_name, "_", names(result))
      } else {
        names(result) <- feature_name
      }
      result
    }, error = function(e) {
      # Return NA with appropriate name
      na_val <- NA_real_
      names(na_val) <- feature_name
      na_val
    })
  }
  
  # Initialize feature list
  features <- list()
  
  # ===== BASIC SUMMARY STATISTIC =====
  # Compute total case count (sum of all sampling events)
  # This is equivalent to tree statistic "number_of_lineages" (total tips)
  # Represents the total number of sampled individuals (I→R transitions)
  # Computed directly from raw case_counts (doesn't require ts conversion)
  # 
  # Note: case_counts should already include all sampled individuals (I→R transitions)
  # including the initial seed infection, so no +1 adjustment needed
  total_case_count <- sum(case_counts, na.rm = TRUE)
  features$total_case_count <- total_case_count
  
  # Convert to ts object for other features
  ts_data <- ts(case_counts, frequency = frequency)
  
  # Compute all fast tsfeatures
  # Listed in order of speed (fastest first)
  
  # Fast features that work on epidemic data (52 scalar features)
  features$zero_proportion <- try_feature("zero_proportion", zero_proportion(ts_data))
  # EXCLUDED: hw_parameters - requires seasonal data, returns NaN on epidemic curves
  # features$hw_parameters <- try_feature("hw_parameters", hw_parameters(ts_data))
  features$flat_spots <- try_feature("flat_spots", flat_spots(ts_data))
  features$crossing_points <- try_feature("crossing_points", crossing_points(ts_data))
  features$max_level_shift <- try_feature("max_level_shift", max_level_shift(ts_data))
  features$max_var_shift <- try_feature("max_var_shift", max_var_shift(ts_data))
  features$unitroot_kpss <- try_feature("unitroot_kpss", unitroot_kpss(ts_data))
  features$firstzero_ac <- try_feature("firstzero_ac", firstzero_ac(ts_data))
  features$stability <- try_feature("stability", stability(ts_data))
  features$lumpiness <- try_feature("lumpiness", lumpiness(ts_data))
  features$arch_stat <- try_feature("arch_stat", arch_stat(ts_data))
  features$nonlinearity <- try_feature("nonlinearity", nonlinearity(ts_data))
  features$acf_features <- try_feature("acf_features", acf_features(ts_data))
  features$entropy <- try_feature("entropy", entropy(ts_data))
  features$pacf_features <- try_feature("pacf_features", pacf_features(ts_data))
  features$hurst <- try_feature("hurst", hurst(ts_data))
  features$autocorr_features <- try_feature("autocorr_features", autocorr_features(ts_data))
  features$unitroot_pp <- try_feature("unitroot_pp", unitroot_pp(ts_data))
  features$stl_features <- try_feature("stl_features", stl_features(ts_data))
  
  # Remove constant STL features for non-seasonal data
  # stl_features_nperiods is always 0 (no seasonal periods detected)
  # stl_features_seasonal_period is always 1 (frequency parameter)
  # These provide no information for epidemic time series with frequency=1
  if (!is.null(features$stl_features)) {
    stl_names <- names(features$stl_features)
    keep_stl <- !grepl("stl_features_nperiods$|stl_features_seasonal_period$", stl_names)
    features$stl_features <- features$stl_features[keep_stl]
  }

  features$heterogeneity <- try_feature("heterogeneity", heterogeneity(ts_data))
  features$dist_features <- try_feature("dist_features", dist_features(ts_data))
  features$max_kl_shift <- try_feature("max_kl_shift", max_kl_shift(ts_data))
  features$holt_parameters <- try_feature("holt_parameters", holt_parameters(ts_data))
  # EXCLUDED: station_features - spreadrandomlocal_meantaul_50 returns NaN on epidemic data
  # features$station_features <- try_feature("station_features", station_features(ts_data))
  
  # Flatten to named vector
  features <- unlist(features)
  
  # Sort by name for consistency
  features <- features[order(names(features))]
  
  return(features)
}


# If running as script
if (!interactive()) {
  cat("========================================\n")
  cat("calc_fast_tsfeatures() function loaded\n")
  cat("========================================\n\n")
  cat("Usage:\n")
  cat("  source('calc_fast_tsfeatures.R')\n")
  cat("  case_counts <- c(0, 5, 12, 45, 89, 120, ...)  # Your case count data\n\n")
  cat("  # Calculate fast tsfeatures:\n")
  cat("  features <- calc_fast_tsfeatures(case_counts)\n")
  cat("  # Returns named vector of ~49 features\n\n")
  cat("  # For daily data (default):\n")
  cat("  features <- calc_fast_tsfeatures(case_counts, frequency=1)\n\n")
  cat("Performance:\n")
  cat("  ~0.07s per 1000-point series\n")
  cat("  For 20,000 series: ~0.39 hours total\n\n")
  cat("Features computed (25 functions):\n")
  cat("  total_case_count, acf_features, arch_stat, autocorr_features,\n")
  cat("  crossing_points, dist_features, entropy, firstzero_ac, flat_spots,\n")
  cat("  heterogeneity, holt_parameters, hurst, lumpiness, max_kl_shift,\n")
  cat("  max_level_shift, max_var_shift, nonlinearity, pacf_features,\n")
  cat("  stl_features, stability, unitroot_kpss, unitroot_pp,\n")
  cat("  zero_proportion\n\n")
  cat("Note: heterogeneity() includes arch_acf, garch_acf, arch_r2, garch_r2\n\n")
  cat("Note: total_case_count equals tree statistic 'number_of_lineages'\n")
  cat("      (total outbreak size, cumulative infections)\n\n")
  cat("Note: STL features nperiods and seasonal_period excluded (constant for non-seasonal data)\n\n")
  cat("Excluded (slow/variable): binarize_mean, compengine, pred_features, scal_features\n")
  cat("Excluded (NaN on epidemic data): hw_parameters, station_features\n\n")
  cat("Total: ~49 scalar features (after removing 2 constant STL features)\n")
  cat("========================================\n")
}
