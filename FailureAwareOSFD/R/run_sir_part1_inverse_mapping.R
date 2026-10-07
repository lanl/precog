# Part 1: LHS on PIV/PIT output space using inverse SIR mappings
# This script filters feasible PIV/PIT/s0 combinations and maps from output space back to input space
# Author: AC Murph
# Date: Feb 2026

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 1) {
  stop("Usage: Rscript run_sir_part1_inverse_mapping.R <nsuccesses> [rep]")
}
nsuccesses <- as.numeric(args[1])
rep_num <- if (length(args) >= 2) as.numeric(args[2]) else NULL
rep_suffix <- if (!is.null(rep_num)) sprintf("_rep%d", rep_num) else ""

if (!is.null(rep_num)) {
  cat(sprintf("Running Part 1 with nsuccesses = %d, replicate = %d\n", nsuccesses, rep_num))
} else {
  cat(sprintf("Running Part 1 with nsuccesses = %d\n", nsuccesses))
}

# Load shared setup
source(here::here("R", "sir_experiment_setup.R"))

# -------------------------------
# Filter LHS samples to only feasible PIV/PIT/s0 combinations
# (with caching to avoid recomputation)
# -------------------------------

cat("\n=== PART 1: Filtering LHS samples for feasibility ===\n")

# Create cache directory if it doesn't exist
cache_dir <- here::here("data", "feasible_LHSs")
if (!dir.exists(cache_dir)) {
  dir.create(cache_dir, recursive = TRUE)
  cat("Created cache directory:", cache_dir, "\n")
}

# Check if cached feasible LHS exists
cache_file <- file.path(cache_dir, sprintf("feasible_LHS_%d.RData", nsuccesses))

if (file.exists(cache_file)) {
  cat(sprintf("Loading cached feasible LHS from: %s\n", cache_file))
  load(cache_file)  # This loads X_isfd
  cat(sprintf("Loaded %d cached feasible points\n", nrow(X_isfd)))
} else {
  cat(sprintf("No cached file found. Generating new feasible LHS...\n"))

  # Generate 10x oversample
  n_oversample <- 10 * nsuccesses
  X_candidate <- lhs::randomLHS(n_oversample, 3)  # 3D: PIV, PIT, s0

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

  if (n_feasible < nsuccesses) {
    warning(sprintf("Only found %d feasible points but need %d. Consider increasing oversample factor.",
                    n_feasible, nsuccesses))
    X_isfd <- X_feasible
  } else {
    # Subsample to exactly nsuccesses
    X_isfd <- X_feasible[1:nsuccesses, ]
    cat(sprintf("Using first %d feasible points\n", nsuccesses))
  }

  # Save to cache for future runs
  save(X_isfd, file = cache_file)
  cat(sprintf("Saved feasible LHS to cache: %s\n", cache_file))
}

# -------------------------------
# LHS directly on PIV/PIT space
# -------------------------------

## LHS directly on PIV/PIT space:
X_isfd <- lhs::randomLHS(nsuccesses, p)

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
f <- make_f(pit_bounds, piv_bounds, s0_bounds, reproduction_number_bounds)

cat("\n=== LHS directly on PIV/PIT space ===\n")
cl <- parallel::makeCluster(mc.cores)
doSNOW::registerDoSNOW(cl)

# send the variable to workers
parallel::clusterExport(cl, varlist = c("f", "q", "alpha_bounds", "reproduction_number_bounds", "s0_bounds"), envir = environment())
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
results_list <- foreach::foreach(
  i = seq_len(nrow(X_isfd)),
  .inorder = TRUE,
  .options.snow = opts
) %dopar% {
  x <- X_isfd[i, , drop = FALSE]

  # Call the function directly to get the raw output
  # DISABLED: This evaluates f() twice, doubling computation time
  # raw_output <- tryCatch({
  #   f(x)
  # }, error = function(e) {
  #   list(error = conditionMessage(e))
  # })
  raw_output <- NULL  # Set to NULL to avoid breaking downstream code

  templist = safe_eval_f(f, x, q+5)

  result <- list(
    index = i,
    input = x,
    success = templist$success,
    error_msg = templist$msg,
    raw_output = raw_output,  # Store raw output for diagnosis (currently NULL)
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

close(pb)
parallel::stopCluster(cl)

# Analyze failures
n_success <- sum(sapply(results_list, function(x) x$success && !is.null(x$data)))
n_fail <- length(results_list) - n_success

cat(sprintf("\n\nResults: %d successes, %d failures (%.1f%% failure rate)\n",
            n_success, n_fail, 100 * n_fail / length(results_list)))

# Categorize error messages
if(n_fail > 0) {
  error_msgs <- sapply(results_list, function(x) {
    if(!x$success || is.null(x$data)) x$error_msg else NA_character_
  })
  error_msgs <- error_msgs[!is.na(error_msgs)]

  cat("\nError summary:\n")
  error_table <- sort(table(error_msgs), decreasing = TRUE)
  for(i in seq_along(error_table)) {
    cat(sprintf("  [%d occurrences] %s\n", error_table[i], names(error_table)[i]))
  }

  # Detailed diagnosis for "wrong length or non-finite" errors
  # DISABLED: raw_output evaluation is disabled to avoid double computation
  # wrong_length_cases <- which(error_msgs == "Returned wrong length or non-finite values.")
  # if(length(wrong_length_cases) > 0) {
  #   cat("\n\n=== Detailed diagnosis of 'wrong length or non-finite' errors ===\n")
  #
  #   # Sample first 10 cases
  #   sample_indices <- head(wrong_length_cases, 10)
  #
  #   for(idx in sample_indices) {
  #     result_item <- results_list[[idx]]
  #     raw <- result_item$raw_output
  #
  #     cat(sprintf("\nCase %d (input: PIT_norm=%.3f, PIV_norm=%.3f, s0_norm=%.3f):\n",
  #                 result_item$index,
  #                 result_item$input[1],
  #                 result_item$input[2],
  #                 result_item$input[3]))
  #
  #     if(is.null(raw)) {
  #       cat(sprintf("  Raw output not available (double evaluation disabled)\n"))
  #     } else if(is.list(raw) && !is.null(raw$error)) {
  #       cat(sprintf("  Error during function call: %s\n", raw$error))
  #     } else if(is.list(raw)) {
  #       cat(sprintf("  Returned list with names: %s\n", paste(names(raw), collapse=", ")))
  #       cat(sprintf("  List elements:\n"))
  #       for(nm in names(raw)) {
  #         val <- raw[[nm]]
  #         if(is.numeric(val)) {
  #           cat(sprintf("    %s: %.6f (finite=%s)\n", nm, val, is.finite(val)))
  #         } else {
  #           cat(sprintf("    %s: %s (class=%s)\n", nm, as.character(val), class(val)[1]))
  #         }
  #       }
  #     } else {
  #       cat(sprintf("  Returned type: %s, length: %d, expected: %d\n",
  #                   class(raw)[1], length(raw), q+5))
  #       if(is.numeric(raw) && length(raw) <= 20) {
  #         cat(sprintf("  Values: %s\n", paste(sprintf("%.4f", raw), collapse=", ")))
  #         cat(sprintf("  Finite? %s\n", paste(is.finite(raw), collapse=", ")))
  #       }
  #     }
  #   }
  # }

  # Save detailed error log for later analysis
  error_log <- data.frame(
    index = sapply(results_list, function(x) x$index),
    success = sapply(results_list, function(x) x$success),
    error_msg = sapply(results_list, function(x) if(!x$success || is.null(x$data)) x$error_msg else NA_character_),
    stringsAsFactors = FALSE
  )

  # Save the full results_list for deep dive analysis
  save(results_list, file = here::here("data", sprintf("full_results_list_sirMappings_n%d%s.RData", nsuccesses, rep_suffix)))
  save(error_log, file = here::here("data", sprintf("error_log_sirMappings_n%d%s.RData", nsuccesses, rep_suffix)))
  cat(sprintf("\n\nDetailed error log saved to: data/error_log_sirMappings_n%d%s.RData", nsuccesses, rep_suffix))
  cat(sprintf("\nFull results list saved to: data/full_results_list_sirMappings_n%d%s.RData\n\n", nsuccesses, rep_suffix))
}

# Extract successful results
results <- do.call(rbind, lapply(results_list, function(x) x$data))
cat(sprintf("Completed with %d valid results\n\n", nrow(results)))

# Save results with nsuccesses in filename
output_file <- here::here("data", sprintf("inputs_and_outputs_sirMappings_n%d%s.RData", nsuccesses, rep_suffix))
save(results, file = output_file)
cat(sprintf("Results saved to: %s\n", output_file))

cat("\n=== PART 1 COMPLETE ===\n")
