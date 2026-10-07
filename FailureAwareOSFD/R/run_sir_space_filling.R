# For this script -- I am going to output-space fill on PIT and PIV to SIR curves.
# This will use three methods.  1) basic LHS, 2) my SIR paper stuff, 3) this output-space filling design stuff
## Author: AC Murph
## Date Feb 2026
library(parallel)
library(doParallel)
library(doSNOW)
library(ggplot2)
library(patchwork)
# Load helper functions
source(here::here("R", "source_helpers.R"))
source_helpers()

setwd(here::here())

nsuccesses = 50000 
p = 3
q = 2
mc.cores = 99
cand_batch = 30000
CAND_size = 500000
pkg_location = here::here()
alpha_bounds = c(0.001, 100)
reproduction_number_bounds = c(1.001, 100)
pit_bounds = c(3, 50)
piv_bounds = c(0.005, 1)
s0_bounds = c(0.95, 0.999)

# -------------------------------
# -------------------------------
## Filter LHS samples to only feasible PIV/PIT/s0 combinations, to then use for the SIR inverse maps LHS on the output space
source(here::here("R", "sir_piv_pit_feasible.R"))

cat("\n=== Filtering LHS samples for feasibility ===\n")
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

## LHS directly on PIV/PIT space:
X_isfd <- lhs::randomLHS(nsuccesses, p) 
# f <- function(x) {
#   calculate_sir_from_PIVPIT(x[1], x[2], x[3],pit_bounds = pit_bounds,
#                                       piv_bounds = piv_bounds,
#                                       s0_bounds = s0_bounds, max_reproduction_number = max(reproduction_number_bounds))
# }
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
  raw_output <- tryCatch({
    f(x)
  }, error = function(e) {
    list(error = conditionMessage(e))
  })

  templist = safe_eval_f(f, x, q+5)

  result <- list(
    index = i,
    input = x,
    success = templist$success,
    error_msg = templist$msg,
    raw_output = raw_output,  # Store raw output for diagnosis
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
        beta = templist$y[4]
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
  wrong_length_cases <- which(error_msgs == "Returned wrong length or non-finite values.")
  if(length(wrong_length_cases) > 0) {
    cat("\n\n=== Detailed diagnosis of 'wrong length or non-finite' errors ===\n")

    # Sample first 10 cases
    sample_indices <- head(wrong_length_cases, 10)

    for(idx in sample_indices) {
      result_item <- results_list[[idx]]
      raw <- result_item$raw_output

      cat(sprintf("\nCase %d (input: PIT_norm=%.3f, PIV_norm=%.3f, s0_norm=%.3f):\n",
                  result_item$index,
                  result_item$input[1],
                  result_item$input[2],
                  result_item$input[3]))

      if(is.list(raw) && !is.null(raw$error)) {
        cat(sprintf("  Error during function call: %s\n", raw$error))
      } else if(is.list(raw)) {
        cat(sprintf("  Returned list with names: %s\n", paste(names(raw), collapse=", ")))
        cat(sprintf("  List elements:\n"))
        for(nm in names(raw)) {
          val <- raw[[nm]]
          if(is.numeric(val)) {
            cat(sprintf("    %s: %.6f (finite=%s)\n", nm, val, is.finite(val)))
          } else {
            cat(sprintf("    %s: %s (class=%s)\n", nm, as.character(val), class(val)[1]))
          }
        }
      } else {
        cat(sprintf("  Returned type: %s, length: %d, expected: %d\n",
                    class(raw)[1], length(raw), q+5))
        if(is.numeric(raw) && length(raw) <= 20) {
          cat(sprintf("  Values: %s\n", paste(sprintf("%.4f", raw), collapse=", ")))
          cat(sprintf("  Finite? %s\n", paste(is.finite(raw), collapse=", ")))
        }
      }
    }
  }

  # Save detailed error log for later analysis
  error_log <- data.frame(
    index = sapply(results_list, function(x) x$index),
    success = sapply(results_list, function(x) x$success),
    error_msg = sapply(results_list, function(x) if(!x$success || is.null(x$data)) x$error_msg else NA_character_),
    stringsAsFactors = FALSE
  )

  # Save the full results_list for deep dive analysis
  save(results_list, file = here::here("data", "full_results_list_sirMappings.RData"))
  save(error_log, file = here::here("data", "error_log_sirMappings.RData"))
  cat(sprintf("\n\nDetailed error log saved to: data/error_log_sirMappings.RData"))
  cat(sprintf("\nFull results list saved to: data/full_results_list_sirMappings.RData\n\n"))
}

# Extract successful results
results <- do.call(rbind, lapply(results_list, function(x) x$data))
cat(sprintf("Completed with %d valid results\n\n", nrow(results)))
save(results, file=here::here("data", "inputs_and_outputs_sirMappings.RData"))

# -------------------------------
# -------------------------------
# Regular LHS on input space (alpha/beta)
X_isfd <- lhs::randomLHS(nsuccesses, p) 
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

cat("\n=== Regular LHS on input space (alpha/beta) ===\n")
cl <- parallel::makeCluster(mc.cores)
doSNOW::registerDoSNOW(cl)

# send the variable to workers
parallel::clusterExport(cl, varlist = c("f", "q"), envir = environment())
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

results <- foreach::foreach(
  i = seq_len(nrow(X_isfd)),
  .combine = 'rbind',
  .multicombine = TRUE,
  .inorder = TRUE,
  .options.snow = opts
) %dopar% {
  x <- X_isfd[i, , drop = FALSE]

  alpha = alpha_bounds[1] + x[1] * diff(alpha_bounds)
  rho = reproduction_number_bounds[1] + x[2] * diff(reproduction_number_bounds)
  beta = alpha / rho

  templist = safe_eval_f(f, x, q)
  if(templist$success){
    data.frame(piv = templist$y[1], pit = templist$y[2], alpha = alpha, beta = beta)
  }
}

close(pb)
parallel::stopCluster(cl)
cat(sprintf("Completed %d tasks\n\n", n_tasks)) 
save(results, file=here::here("data", "inputs_and_outputs_sirBasicLHS.RData"))

# -------------------------------
# -------------------------------
# Direct OSFD
## LHS on output space:
CAND <- lhs::randomLHS(CAND_size, p)
print("About to enter OSFD function.")
res <- constrained_osfd_ei_with_retry(
  f, p, q,
  CAND,
  nsuccesses,
  n_ini_success = floor(2 * nsuccesses / 4),
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
  pkg_location = pkg_location
)
print("We have now exited the OSFD function.")

D <- res$D_success
Y <- res$Y_success
save(D, file=here::here("data", "inputs_sirOutputFilling.RData"))
save(Y, file=here::here("data", "outputs_sirOutputFilling.RData"))
save(res, file=here::here("data", "fullOutput_OSFD_basic.RData"))


# -------------------------------
# -------------------------------
# Doing OSFD using an initial sample from the INVERSE SIR maps.

## Pre-calculate LHS initial
n_ini_success = floor(2 * nsuccesses / 4)
X_isfd <- lhs::randomLHS(nsuccesses, p) 
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

cat("\n=== Pre-calculate LHS initial using SIR mapping ===\n")
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
    } else if(tmp_alpha < min(alpha_bounds)) {
      result$success <- FALSE
      result$error_msg <- sprintf("alpha out of bounds: %.4f < %.4f", tmp_alpha, min(alpha_bounds))
    } else if(tmp_alpha > max(alpha_bounds)) {
      result$success <- FALSE
      result$error_msg <- sprintf("alpha out of bounds: %.4f > %.4f", tmp_alpha, max(alpha_bounds))
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

X_all = NULL
succ_flag = c()
fail_msg = c()
X_succ = NULL
Y_succ = NULL
remaining <- seq_len(nrow(CAND))

for(ii in 1:nrow(reses)){
  # if(is.null(reses[[ii]])) next
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

save(X_all, file = paste0(pkg_location, "/data/X_all.RData") )
save(succ_flag, file = paste0(pkg_location, "/data/succ_flag.RData") )
save(fail_msg, file = paste0(pkg_location, "/data/fail_msg.RData") )
save(X_succ, file = paste0(pkg_location, "/data/X_succ.RData") )
save(Y_succ, file = paste0(pkg_location, "/data/Y_succ.RData") )
save(remaining, file = paste0(pkg_location, "/data/remaining.RData") )


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
  pkg_location = pkg_location,
  pre_calculate_LHS = TRUE
)
print("We have now exited the OSFD function.")

D <- res$D_success
Y <- res$Y_success
save(D, file=here::here("data", "inputs_sirOutputFilling_wSIRStuff.RData"))
save(Y, file=here::here("data", "outputs_sirOutputFilling_wSIRStuff.RData"))
save(res, file=here::here("data", "fullOutput_OSFD_wSIRStuff.RData"))