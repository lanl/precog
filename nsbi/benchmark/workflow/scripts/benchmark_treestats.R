#!/usr/bin/env Rscript

#' Benchmark all treestats summary statistics with intelligent timeouts
#' Identifies fast vs slow statistics for large phylogenetic trees
#' 
#' Features:
#' - Timeout mechanism to skip extremely slow statistics
#' - Skip slow statistics on larger trees to save time
#' - Progressive benchmarking from small to large trees

library(ape)
library(treestats)

cat("=== TreeStats Benchmarking Suite ===\n\n")

# Configuration
REAL_TREE_PATH <- "../remaster/sim_sir_trees.full.trees"
TREE_SIZES <- c(1000, 5000, 10000)  # Progressive sizes
TIMEOUT_SECONDS <- 30  # Max time per statistic
SLOW_THRESHOLD <- 1.0   # Statistics slower than this won't be tested on larger trees
TIME_THRESHOLD_FAST <- 0.1  # Classify as fast
TIME_THRESHOLD_SLOW <- 1.0  # Classify as slow
N_REPS <- 3

# Track slow statistics to skip on larger trees
slow_stats <- character(0)

# Function name mappings: output_name -> actual_function_name
# Some statistics have different function names than their output names
get_function_name_mapping <- function() {
  list(
    beta = "beta_statistic",
    gamma = "gamma_statistic",
    il_number = "ILnumber",
    symmetry_nodes = "sym_nodes",
    mpd = "mean_pair_dist",
    vpd = "var_pair_dist",
    phylogenetic_div = "phylogenetic_diversity",
    var_depth = "var_leaf_depth",
    j_stat = "entropy_j",
    tot_path = "tot_path_length",
    i_stat = "mean_i",
    nltt_base = "nLTT_base"
  )
}

# Special parameters needed for certain statistics
get_stat_parameters <- function(stat_name) {
  params <- list(
    il_number = list(normalization = "none"),
    symmetry_nodes = list(normalization = "none"),
    max_depth = list(normalization = "none"),
    max_width = list(normalization = "none"),
    average_leaf_depth = list(normalization = "none"),
    max_del_width = list(normalization = "none"),
    tot_coph = list(normalization = "none"),
    area_per_pair = list(normalization = "none"),
    wiener = list(normalization = FALSE),
    mpd = list(normalization = "none"),
    var_depth = list(normalization = "none"),
    max_betweenness = list(normalization = "none"),
    max_closeness = list(weight = FALSE, normalization = "none"),
    psv = list(normalization = "none"),
    imbalance_steps = list(normalization = FALSE)
  )
  
  if (stat_name %in% names(params)) {
    return(params[[stat_name]])
  } else {
    return(list())
  }
}

# Statistics that require ultrametric trees
ultrametric_only_stats <- function() {
  c("gamma", "imbalance_steps", "nltt_base")
}

# Statistics that need complex extraction (not simple function calls)
get_complex_stat_handler <- function(stat_name) {
  handlers <- list(
    laplace_spectrum_a = function(tree) laplacian_spectrum(tree)$asymmetry,
    laplace_spectrum_e = function(tree) laplacian_spectrum(tree)$eigengap,
    laplace_spectrum_g = function(tree) laplacian_spectrum(tree)$principal_eigenvalue,
    laplace_spectrum_p = function(tree) laplacian_spectrum(tree)$peakedness,
    max_laplace = function(tree) max(laplacian_spectrum(tree)$eigenvalues),
    min_laplace = function(tree) min(laplacian_spectrum(tree)$eigenvalues),
    max_adj = function(tree) minmax_adj(tree)$max,
    min_adj = function(tree) minmax_adj(tree)$min,
    eigen_centralityW = function(tree) max(eigen_centrality(tree, weight = TRUE)[[1]]),
    max_closenessW = function(tree) max_closeness(tree, weight = TRUE, normalization = "none")
  )
  
  if (stat_name %in% names(handlers)) {
    return(handlers[[stat_name]])
  } else {
    return(NULL)
  }
}

# Helper function to time a single statistic with timeout
time_statistic_with_timeout <- function(tree, tree_ultra, stat_name, timeout_sec = 30) {
  
  # Check if this is a complex stat that needs special handling
  complex_handler <- get_complex_stat_handler(stat_name)
  
  if (!is.null(complex_handler)) {
    # Use the complex handler directly
    result <- tryCatch({
      setTimeLimit(cpu = timeout_sec, elapsed = timeout_sec, transient = TRUE)
      start_time <- Sys.time()
      complex_handler(tree)
      end_time <- Sys.time()
      setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
      as.numeric(difftime(end_time, start_time, units = "secs"))
    }, error = function(e) {
      setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
      if (grepl("time limit", e$message, ignore.case = TRUE)) {
        return("TIMEOUT")
      }
      return(NA_real_)
    })
  } else {
    # Standard function call
    
    # Check if ultrametric-only and use appropriate tree
    use_tree <- if (stat_name %in% ultrametric_only_stats()) tree_ultra else tree
    
    # Map output name to actual function name if needed
    func_mappings <- get_function_name_mapping()
    func_name <- if (stat_name %in% names(func_mappings)) {
      func_mappings[[stat_name]]
    } else {
      stat_name
    }
    
    # Get the function for this statistic
    stat_func <- tryCatch({
      get(func_name, envir = asNamespace("treestats"))
    }, error = function(e) {
      return(NULL)
    })
    
    if (is.null(stat_func)) {
      return(list(mean_time = NA, min_time = NA, max_time = NA, 
                  failed = TRUE, timed_out = FALSE))
    }
    
    # Get any special parameters needed for this statistic
    stat_params <- get_stat_parameters(stat_name)
    
    # Try to run with timeout using setTimeLimit
    result <- tryCatch({
      setTimeLimit(cpu = timeout_sec, elapsed = timeout_sec, transient = TRUE)
      start_time <- Sys.time()
      do.call(stat_func, c(list(use_tree), stat_params))
      end_time <- Sys.time()
      setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
      as.numeric(difftime(end_time, start_time, units = "secs"))
    }, error = function(e) {
      setTimeLimit(cpu = Inf, elapsed = Inf, transient = FALSE)
      if (grepl("time limit", e$message, ignore.case = TRUE)) {
        return("TIMEOUT")
      }
      return(NA_real_)
    })
  }
  
  if (identical(result, "TIMEOUT")) {
    return(list(mean_time = timeout_sec, min_time = timeout_sec, 
                max_time = timeout_sec, failed = FALSE, timed_out = TRUE))
  }
  
  if (is.na(result)) {
    return(list(mean_time = NA, min_time = NA, max_time = NA, 
                failed = TRUE, timed_out = FALSE))
  }
  
  # If first run succeeded and was fast enough, do more reps
  if (result < timeout_sec / 3) {
    times <- numeric(N_REPS)
    times[1] <- result
    
    for (i in 2:N_REPS) {
      if (!is.null(complex_handler)) {
        # Use complex handler
        times[i] <- tryCatch({
          start_time <- Sys.time()
          complex_handler(tree)
          end_time <- Sys.time()
          as.numeric(difftime(end_time, start_time, units = "secs"))
        }, error = function(e) {
          NA_real_
        })
      } else {
        # Use standard function
        use_tree <- if (stat_name %in% ultrametric_only_stats()) tree_ultra else tree
        times[i] <- tryCatch({
          start_time <- Sys.time()
          do.call(stat_func, c(list(use_tree), stat_params))
          end_time <- Sys.time()
          as.numeric(difftime(end_time, start_time, units = "secs"))
        }, error = function(e) {
          NA_real_
        })
      }
    }
    
    return(list(
      mean_time = mean(times, na.rm = TRUE),
      min_time = min(times, na.rm = TRUE),
      max_time = max(times, na.rm = TRUE),
      failed = all(is.na(times)),
      timed_out = FALSE
    ))
  } else {
    # First run was slow, just use that measurement
    return(list(
      mean_time = result,
      min_time = result,
      max_time = result,
      failed = FALSE,
      timed_out = FALSE
    ))
  }
}

# Generate synthetic trees of specified size
generate_trees <- function(n_tips) {
  cat(sprintf("\nGenerating synthetic trees with %d tips...\n", n_tips))
  
  # Non-ultrametric tree (standard)
  tree <- rtree(n_tips)
  tree$edge.length <- abs(tree$edge.length)
  
  # Ultrametric tree (coalescent)
  tree_ultra <- rcoal(n_tips)
  
  return(list(tree = tree, tree_ultra = tree_ultra))
}

# Main benchmarking function
benchmark_all_stats <- function() {
  
  # Get list of all available statistics
  cat("Getting list of all available treestats statistics...\n")
  all_stats <- treestats::list_statistics()
  cat(sprintf("Found %d statistics to benchmark\n\n", length(all_stats)))
  
  # Initialize results data frame
  results <- data.frame(
    statistic = character(),
    tree_size = integer(),
    mean_time = numeric(),
    min_time = numeric(),
    max_time = numeric(),
    fast = logical(),
    failed = logical(),
    timed_out = logical(),
    skipped = logical(),
    stringsAsFactors = FALSE
  )
  
  # Benchmark on progressively larger trees
  for (size in TREE_SIZES) {
    cat(sprintf("\n========================================\n"))
    cat(sprintf("=== Benchmarking on %d tip tree ===\n", size))
    cat(sprintf("========================================\n"))
    
    trees <- generate_trees(size)
    tree <- trees$tree
    tree_ultra <- trees$tree_ultra
    n_to_test <- length(all_stats) - length(slow_stats)
    tested <- 0
    
    for (i in seq_along(all_stats)) {
      stat_name <- all_stats[i]
      
      # Skip if identified as slow on smaller tree
      if (stat_name %in% slow_stats) {
        cat(sprintf("[SKIP] %s (too slow on smaller trees)\n", stat_name))
        results <- rbind(results, data.frame(
          statistic = stat_name,
          tree_size = size,
          mean_time = NA,
          min_time = NA,
          max_time = NA,
          fast = FALSE,
          failed = FALSE,
          timed_out = FALSE,
          skipped = TRUE,
          stringsAsFactors = FALSE
        ))
        next
      }
      
      tested <- tested + 1
      cat(sprintf("[%d/%d] Testing %s... ", tested, n_to_test, stat_name))
      flush.console()
      
      timing <- time_statistic_with_timeout(tree, tree_ultra, stat_name, 
                                           timeout_sec = TIMEOUT_SECONDS)
      
      if (timing$timed_out) {
        cat(sprintf("TIMEOUT (>%ds)\n", TIMEOUT_SECONDS))
        slow_stats <<- c(slow_stats, stat_name)
      } else if (timing$failed) {
        cat("FAILED\n")
      } else {
        cat(sprintf("%.4f sec", timing$mean_time))
        if (timing$mean_time > SLOW_THRESHOLD) {
          cat(" [SLOW - will skip on larger trees]")
          slow_stats <<- c(slow_stats, stat_name)
        }
        cat("\n")
      }
      
      results <- rbind(results, data.frame(
        statistic = stat_name,
        tree_size = size,
        mean_time = timing$mean_time,
        min_time = timing$min_time,
        max_time = timing$max_time,
        fast = !timing$failed && !timing$timed_out && timing$mean_time < TIME_THRESHOLD_FAST,
        failed = timing$failed,
        timed_out = timing$timed_out,
        skipped = FALSE,
        stringsAsFactors = FALSE
      ))
    }
    
    cat(sprintf("\nSummary for %d tips: %d slow stats identified, %d will be skipped on larger trees\n", 
                size, length(slow_stats), length(slow_stats)))
  }
  
  return(results)
}

# Run benchmarking
cat("Starting intelligent benchmark with timeout protection...\n")
cat(sprintf("- Maximum time per statistic: %d seconds\n", TIMEOUT_SECONDS))
cat(sprintf("- Statistics slower than %.1fs will be skipped on larger trees\n", SLOW_THRESHOLD))
cat(sprintf("- Testing tree sizes: %s tips\n\n", paste(TREE_SIZES, collapse = ", ")))

start_total <- Sys.time()
results <- benchmark_all_stats()
end_total <- Sys.time()

total_time <- as.numeric(difftime(end_total, start_total, units = "secs"))
cat(sprintf("\n========================================\n"))
cat(sprintf("Total benchmarking time: %.2f seconds (%.2f minutes)\n", 
            total_time, total_time/60))
cat(sprintf("========================================\n\n"))

# Save raw results
output_file <- "benchmark_results.csv"
write.csv(results, output_file, row.names = FALSE)
cat(sprintf("Raw results saved to: %s\n", output_file))

# Analyze results for largest tree size
cat("\n=== ANALYSIS ===\n\n")
largest_size <- max(TREE_SIZES)
large_tree_results <- results[results$tree_size == largest_size & 
                               !results$failed & 
                               !results$skipped & 
                               !results$timed_out, ]

if (nrow(large_tree_results) > 0) {
  large_tree_results <- large_tree_results[order(large_tree_results$mean_time), ]
  
  fast_stats <- large_tree_results[large_tree_results$mean_time < TIME_THRESHOLD_FAST, ]
  moderate_stats <- large_tree_results[large_tree_results$mean_time >= TIME_THRESHOLD_FAST & 
                                       large_tree_results$mean_time < TIME_THRESHOLD_SLOW, ]
  slow_stats_analyzed <- large_tree_results[large_tree_results$mean_time >= TIME_THRESHOLD_SLOW, ]
  
  # Also count timed out and skipped stats
  timeout_count <- sum(results$tree_size == largest_size & results$timed_out)
  skipped_count <- sum(results$tree_size == largest_size & results$skipped)
  
  cat(sprintf("Results for trees with %d tips:\n", largest_size))
  cat(sprintf("  Fast statistics (< %.1fs): %d\n", TIME_THRESHOLD_FAST, nrow(fast_stats)))
  cat(sprintf("  Moderate statistics (%.1f-%.1fs): %d\n", 
              TIME_THRESHOLD_FAST, TIME_THRESHOLD_SLOW, nrow(moderate_stats)))
  cat(sprintf("  Slow statistics (> %.1fs): %d\n", TIME_THRESHOLD_SLOW, nrow(slow_stats_analyzed)))
  cat(sprintf("  Timed out statistics (>%ds): %d\n", TIMEOUT_SECONDS, timeout_count))
  cat(sprintf("  Skipped statistics (slow on smaller trees): %d\n", skipped_count))
  
  # Calculate time savings
  total_all_stats <- sum(large_tree_results$mean_time)
  total_fast_only <- sum(fast_stats$mean_time)
  speedup <- if (total_fast_only > 0) total_all_stats / total_fast_only else 1
  
  cat(sprintf("\n--- Time Estimates for %d tips ---\n", largest_size))
  cat(sprintf("Time to compute ALL successful statistics: %.2f seconds\n", total_all_stats))
  cat(sprintf("Time to compute FAST statistics only: %.2f seconds\n", total_fast_only))
  cat(sprintf("Speedup: %.1fx faster\n", speedup))
  cat(sprintf("Time saved per tree: %.2f seconds\n", total_all_stats - total_fast_only))
  
  # Project for tens of thousands of trees
  cat(sprintf("\n--- Projections for Large Datasets ---\n"))
  for (n_trees in c(10000, 50000, 100000)) {
    all_time_hrs <- (total_all_stats * n_trees) / 3600
    fast_time_hrs <- (total_fast_only * n_trees) / 3600
    cat(sprintf("For %d trees:\n", n_trees))
    cat(sprintf("  All stats: %.1f hours\n", all_time_hrs))
    cat(sprintf("  Fast only: %.1f hours (%.1fx faster)\n", fast_time_hrs, speedup))
  }
  
  # Save fast statistics list
  fast_stat_names <- fast_stats$statistic
  save(fast_stat_names, file = "fast_statistics.RData")
  writeLines(fast_stat_names, "fast_statistics.txt")
  cat(sprintf("\nFast statistics list saved to: fast_statistics.txt\n"))
  
  # Show detailed lists
  cat(sprintf("\n--- Top 10 Fastest Statistics ---\n"))
  n_show <- min(10, nrow(fast_stats))
  for (i in 1:n_show) {
    cat(sprintf("%2d. %s: %.4f sec\n", i, fast_stats$statistic[i], fast_stats$mean_time[i]))
  }
  
  if (nrow(slow_stats_analyzed) > 0 || timeout_count > 0) {
    cat(sprintf("\n--- Slowest/Problematic Statistics (AVOID) ---\n"))
    # Show timed out ones
    timeout_stats <- results[results$tree_size == largest_size & results$timed_out, ]
    if (nrow(timeout_stats) > 0) {
      for (i in 1:nrow(timeout_stats)) {
        cat(sprintf("  %s: TIMEOUT (>%ds)\n", timeout_stats$statistic[i], TIMEOUT_SECONDS))
      }
    }
    # Show slow but completed
    n_show <- min(10, nrow(slow_stats_analyzed))
    if (n_show > 0) {
      for (i in n_show:1) {
        idx <- nrow(slow_stats_analyzed) - i + 1
        cat(sprintf("  %s: %.3f sec\n", 
                    slow_stats_analyzed$statistic[idx], 
                    slow_stats_analyzed$mean_time[idx]))
      }
    }
  }
  
} else {
  cat("Warning: No successful benchmarks on largest tree size!\n")
}

cat("\n=== Benchmark complete! ===\n")
