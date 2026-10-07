# depends: safe_eval_f, fit_feas_model_glm, predict_feas_glm, osfd_scale_Y, estimate_hi_from_outputs, osfd_ei_acquisition
#' Output-Space Filling Design with Expected Improvement and Feasibility Constraints
#'
#' This function implements a sequential design algorithm that fills the output space
#' of a computationally expensive black-box function f: R^p -> R^q. It balances three
#' competing objectives:
#'   1. Filling holes in the observed output space (via Expected Improvement)
#'   2. Avoiding infeasible regions of the input space (via feasibility modeling)
#'   3. Maintaining diversity in the selected input points (via repulsion)
#'
#' The algorithm is particularly useful when:
#'   - The function f is expensive to evaluate
#'   - The function may fail on some inputs (constrained/feasible region)
#'   - You want comprehensive coverage of the output space for emulation or uncertainty quantification
#'
#' @param f Function to evaluate: f(x) returns a length-q numeric vector (the outputs)
#' @param p Integer: dimension of input space
#' @param q Integer: dimension of output space
#' @param CAND Matrix (n_candidates x p): pool of candidate input points to choose from
#' @param n_success Integer: target number of successful function evaluations (ignored if max_attempts is set)
#' @param n_ini_success Integer: number of initial successful points to collect randomly
#' @param batch_size Integer: number of points to evaluate in parallel per iteration
#' @param n_replicates Integer: number of stochastic replicates to run at each input point
#' @param cand_batch Integer: size of candidate pool to score per iteration (memory/speed tradeoff)
#' @param beta Numeric: exponent for feasibility weighting (higher = more conservative)
#' @param p_floor Numeric: minimum feasibility probability (prevents over-confidence)
#' @param tau_hard Numeric or NULL: hard threshold for feasibility (reject if p_feas < tau_hard)
#' @param repel_radius Numeric: radius for input-space repulsion (promotes diversity)
#' @param update_feas_every Integer: refit feasibility model every N attempts (0 = never)
#' @param max_attempts Integer or NULL: maximum number of function evaluations (including failures).
#'   If set, algorithm stops after this many attempts regardless of n_success. If NULL (default),
#'   algorithm runs until n_success successful evaluations are obtained.
#' @param mc.cores Integer: number of parallel cores to use
#' @param verbose Logical: print progress messages
#' @param pkg_location Character: path to package root (for sourcing helper functions)
#' @param pre_calculate_LHS Logical: whether to load pre-calculated LHS initialization
#'
#' @return List containing:
#'   - D_success: matrix of successful input points (n_success x p)
#'   - Y_success: matrix of successful output points (n_success x q)
#'   - X_all: matrix of all attempted input points
#'   - success: logical vector indicating which attempts succeeded
#'   - fail_msg: character vector of failure messages
#'   - remaining_idx: indices of unused candidates from CAND
#'   - feas_model: final feasibility model (GLM)
#'
#' @export
constrained_osfd_ei_with_retry <- function(
  f, p, q,
  CAND,
  n_success,
  n_ini_success = 20,
  batch_size = 8,
  n_replicates = 100,
  cand_batch = 20000,
  beta = 2.0,
  p_floor = 0.05,
  tau_hard = NULL,
  repel_radius = 0.05,
  update_feas_every = 1,
  max_attempts = NULL,
  mc.cores = batch_size,
  verbose = TRUE,
  pkg_location = this.path::here(),
  pre_calculate_LHS = FALSE
) {

  suppressPackageStartupMessages({
    library(parallel) 
    library(doParallel) 
  })

  # ==============================================================================
  # Input validation and setup
  # ==============================================================================

  # Ensure CAND is a matrix with the correct number of columns
  stopifnot(is.matrix(CAND), ncol(CAND) == p)

  # Set up file locations for helper functions (not used currently but kept for compatibility)
  save_eval_location <- here::here("R", "safe_eval_f.R")
  save_generate_location <- here::here("R", "generate_smoa_synthetic.R")
  pkg_location <- here::here()

  # Use PSOCK socket type for parallel cluster (works across platforms)
  sockettype <- "PSOCK"

  # Track which candidates from CAND have not yet been evaluated
  remaining <- seq_len(nrow(CAND))

  # ==============================================================================
  # Initialize storage for results
  # ==============================================================================

  # Successful evaluations:
  X_succ <- matrix(numeric(0), ncol = p)  # inputs that succeeded
  Y_succ <- matrix(numeric(0), ncol = q)  # corresponding outputs

  # All evaluation attempts (successes and failures):
  X_all <- matrix(numeric(0), ncol = p)   # all attempted inputs
  succ_flag <- logical(0)                  # TRUE if evaluation succeeded
  fail_msg <- character(0)                 # error message if failed, "" if succeeded

  # Feasibility model (GLM predicting P(success | input))
  feas_model <- NULL

  # Counter for total evaluation attempts
  attempt <- 0

  # ==============================================================================
  # Set up parallel cluster (only if mc.cores > 1)
  # ==============================================================================

  # Only create a cluster if we actually need parallelism
  # When mc.cores = 1, we'll use serial evaluation (no cluster creation)
  # This avoids nested cluster issues when called from already-parallel code
  use_cluster <- (mc.cores > 1)
  cl <- NULL

  if (use_cluster) {
    # Create parallel cluster with specified number of cores
    cl <- parallel::makeCluster(spec = mc.cores, type = sockettype)
    parallel::setDefaultCluster(cl)
    doParallel::registerDoParallel(cl)

    # Export necessary objects to worker processes
    print(paste("pkg location is:", pkg_location))
    parallel::clusterExport(cl, varlist = c("pkg_location", "f", "q"), envir = environment())

    # Load helper functions on each worker
    parallel::clusterEvalQ(cl, {
      # Load C++ functions first (avoid RcppExports.R namespace issues)
      cpp_file <- here::here("src", "OSFD.cpp")
      if (!exists("approx_gen") || !is.function(approx_gen)) {
        # Add inst/include to the include path for nanoflann.hpp
        Sys.setenv(PKG_CPPFLAGS = paste0("-I", here::here("inst/include")))
        Rcpp::sourceCpp(cpp_file)
      }
      # Load R helper functions (will skip RcppExports.R)
      source(here::here("R", "source_helpers.R"))
      source_helpers()
      NULL
    })
  } else {
    # Serial mode: no cluster, just print message
    print(paste("pkg location is:", pkg_location, "(serial mode, no cluster)"))
  }

  # ==============================================================================
  # Phase 1: Initial space-filling phase
  # Collect n_ini_success successful evaluations using random sampling from CAND
  # This provides a baseline for the feasibility model and output-space coverage
  # ==============================================================================

  if (!pre_calculate_LHS) {
    # Standard initialization: randomly sample until we get n_ini_success successes

    while (nrow(X_succ) < n_ini_success) {
      # Check if we've exhausted the candidate pool
      if (length(remaining) == 0) {
        stop("Candidate pool exhausted during initialization.")
      }

      # Determine how many more successes we need
      num_to_try <- n_ini_success - nrow(X_succ)
      attempt <- attempt + num_to_try

      # Randomly sample candidate indices to try
      js <- sample(remaining, num_to_try)

      # Current number of successes (used for seeding)
      xx <- nrow(X_succ)

      # Export current state to workers (only if using cluster)
      if (use_cluster) {
        parallel::clusterExport(cl, varlist = c("pkg_location", "f", "q", "xx"), envir = environment())
      }

      # Evaluate candidates (parallel if cluster, serial if not)
      if (use_cluster) {
        reses <- foreach::foreach(
          i = js,
          .verbose = FALSE
        ) %dopar% {
          # Set seed for reproducibility (different for each point)
          set.seed(xx * length(js) + i)

          # Get the candidate input point
          x <- CAND[i, ]

          # Safely evaluate f(x), catching any errors
          res <- safe_eval_f(f, x, q)

          # Return both input and result
          list(x = x, res = res)
        }
      } else {
        # Serial evaluation using regular for loop (avoid foreach nesting issues)
        reses <- vector("list", length(js))
        for (idx in seq_along(js)) {
          i <- js[idx]

          # Set seed for reproducibility (different for each point)
          set.seed(xx * length(js) + i)

          # Get the candidate input point
          x <- CAND[i, ]

          # Safely evaluate f(x), catching any errors
          res <- safe_eval_f(f, x, q)

          # Return both input and result
          reses[[idx]] <- list(x = x, res = res)
        }
      }

      # Remove attempted candidates from the remaining pool
      remaining <- setdiff(remaining, js)

      # Process results from this batch
      for (ii in 1:length(reses)) {
        # Store the input point in X_all
        X_all <- rbind(X_all, matrix(reses[[ii]]$x, nrow = 1))

        # Store success status
        succ_flag <- c(succ_flag, reses[[ii]]$res$success)

        # Store failure message (empty string if successful)
        fail_msg <- c(fail_msg, ifelse(reses[[ii]]$res$success, "", reses[[ii]]$res$msg))

        # If successful, store in the success matrices
        if (reses[[ii]]$res$success) {
          X_succ <- rbind(X_succ, matrix(reses[[ii]]$x, nrow = 1))
          Y_succ <- rbind(Y_succ, matrix(reses[[ii]]$res$y, nrow = 1))
        }
      }

      # Print progress
      if (verbose) {
        cat("  attempts:", attempt, " successes:", nrow(X_succ), "\n")
      }
    }

    if (verbose) {
      cat("Sequential (batched) phase: targeting", n_success, "successful points...\n")
    }

  } else {
    # Alternative: load pre-calculated initialization from saved files
    # This is useful for resuming interrupted runs or testing with fixed initial conditions

    load(file = paste0(pkg_location, "/data/X_all.RData"))
    load(file = paste0(pkg_location, "/data/succ_flag.RData"))
    load(file = paste0(pkg_location, "/data/fail_msg.RData"))
    load(file = paste0(pkg_location, "/data/X_succ.RData"))
    load(file = paste0(pkg_location, "/data/Y_succ.RData"))
    load(file = paste0(pkg_location, "/data/remaining.RData"))
  }

  # ==============================================================================
  # Phase 2: Sequential adaptive phase
  # Iteratively select batches of points that fill holes in the output space
  # while avoiding infeasible regions and maintaining input diversity
  # ==============================================================================

  # Determine stopping criterion based on max_attempts or n_success
  use_attempt_budget <- !is.null(max_attempts)

  while (TRUE) {
    # Check stopping criteria
    if (use_attempt_budget) {
      # Stop if we've reached the attempt budget
      if (attempt >= max_attempts) {
        if (verbose) cat("Reached maximum attempts budget:", max_attempts, "\n")
        break
      }
    } else {
      # Stop if we've reached the success target
      if (nrow(X_succ) >= n_success) {
        break
      }
    }

    # Check if we've exhausted the candidate pool
    if (length(remaining) == 0) {
      if (use_attempt_budget) {
        warning("Candidate pool exhausted before reaching max_attempts.")
      } else {
        warning("Candidate pool exhausted before reaching n_success.")
      }
      break
    }

    # Increment attempt counter
    attempt <- attempt + batch_size

    # --------------------------------------------------------------------------
    # Step 1: Update feasibility model
    # Fit a GLM to predict P(success | input) based on all attempts so far
    # --------------------------------------------------------------------------

    if (update_feas_every > 0 && (attempt %% update_feas_every == 0)) {
      feas_model <- fit_feas_model_glm(X_all, succ_flag)
    }

    # --------------------------------------------------------------------------
    # Step 2: Compute hole-filling metric in output space
    # --------------------------------------------------------------------------

    # Scale outputs to [0,1]^q for distance calculations
    t0 <- proc.time()[["elapsed"]]
    Y_scaled <- osfd_scale_Y(Y_succ)
    t1 <- proc.time()[["elapsed"]]
    cat(sprintf(
      "osfd_scale_Y(Y_succ) took %.3f sec (%.2f min; %.2f hr)\n",
      (t1 - t0), (t1 - t0) / 60, (t1 - t0) / 3600
    ))

    # Estimate the "hole index" hi for each observed output point
    # hi[i] measures how much of a "hole" exists at Y_scaled[i,]
    # Higher hi = more of a hole = prioritize filling nearby regions
    t2 <- proc.time()[["elapsed"]]
    hi <- estimate_hi_from_outputs(Y_scaled)
    hi[is.na(hi)] <- 0  # Replace any NAs with 0 (no hole)
    t3 <- proc.time()[["elapsed"]]
    cat(sprintf(
      "estimate_hi_from_outputs(Y_scaled) took %.3f sec (%.2f min; %.2f hr)\n",
      (t3 - t2), (t3 - t2) / 60, (t3 - t2) / 3600
    ))

    # --------------------------------------------------------------------------
    # Step 3: Sample a pool of candidates and compute acquisition scores
    # --------------------------------------------------------------------------

    # Randomly sample up to cand_batch candidates from the remaining pool
    # (we don't score all remaining candidates for computational efficiency)
    k <- min(cand_batch, length(remaining))
    pool_idx <- sample(remaining, k, replace = FALSE)
    X_pool <- CAND[pool_idx, , drop = FALSE]

    # Predict feasibility probability for each candidate
    p_feas <- predict_feas_glm(feas_model, X_pool, fallback_p = mean(succ_flag))

    # Enforce floor and ceiling on feasibility probabilities
    # Floor prevents over-confidence in infeasibility
    p_feas <- pmin(1, pmax(p_floor, p_feas))

    # Apply hard threshold if specified (reject candidates below threshold)
    if (!is.null(tau_hard)) {
      ok <- p_feas >= tau_hard
      # If no candidates pass, ignore the threshold for this iteration
      if (!any(ok)) ok <- rep(TRUE, length(p_feas))
    } else {
      ok <- rep(TRUE, length(p_feas))
    }

    # Filter to candidates that pass the threshold
    X_eval <- X_pool[ok, , drop = FALSE]
    idx_eval <- pool_idx[ok]
    p_eval <- p_feas[ok]

    # Compute Expected Improvement for filling output space holes
    t_ei0 <- proc.time()[["elapsed"]]
    EI <- osfd_ei_acquisition(X_obs = X_succ, hi = hi, X_cand = X_eval)
    t_ei1 <- proc.time()[["elapsed"]]
    cat(sprintf(
      "osfd_ei_acquisition(X_succ, hi, X_eval) took %.3f sec (%.2f min; %.2f hr)\n",
      (t_ei1 - t_ei0), (t_ei1 - t_ei0) / 60, (t_ei1 - t_ei0) / 3600
    ))

    # Compute final acquisition score combining EI and feasibility
    # Higher beta = more conservative (stronger penalization of low-feasibility regions)
    score <- EI * (p_eval ^ beta)

    # --------------------------------------------------------------------------
    # Step 4: Greedily select a diverse batch using repulsion
    # --------------------------------------------------------------------------

    chosen_idx <- integer(0)          # indices in CAND of chosen points
    chosen_X <- matrix(numeric(0), ncol = p)  # actual chosen input points
    score_work <- score                # working copy of scores (modified during selection)

    # Greedily select batch_size points
    for (b in seq_len(batch_size)) {
      # Stop if no valid candidates remain
      if (length(score_work) == 0 || all(!is.finite(score_work)) || max(score_work, na.rm = TRUE) <= 0) {
        break
      }

      # Select point with highest score
      jloc <- which.max(score_work)
      chosen_idx <- c(chosen_idx, idx_eval[jloc])
      chosen_X <- rbind(chosen_X, X_eval[jloc, , drop = FALSE])

      # Apply repulsion: downweight scores of nearby points in input space
      # This promotes diversity in the selected batch
      dx <- sqrt(rowSums((X_eval - matrix(X_eval[jloc, ], nrow(X_eval), p, byrow = TRUE))^2))
      score_work <- score_work * pmin(1, dx / repel_radius)

      # Blacklist the just-selected point (set score to -Inf)
      score_work[jloc] <- -Inf
    }

    # Check if we successfully selected any points
    if (length(chosen_idx) == 0) {
      warning("Could not propose any batch points with positive score; relaxing feasibility gate might help.")
      next
    }

    # --------------------------------------------------------------------------
    # Step 5: Evaluate the selected batch in parallel with replicates
    # --------------------------------------------------------------------------

    if (verbose) cat("  proposing batch of", length(chosen_idx), "...\n")

    # Create replicate indices: each point in chosen_X is evaluated n_replicates times
    temp_list <- rep(seq_len(nrow(chosen_X)), each = n_replicates)

    # Export batch to workers (only if using cluster)
    if (use_cluster) {
      parallel::clusterExport(
        cl,
        varlist = c("pkg_location", "f", "q", "chosen_X", "n_replicates", "temp_list"),
        envir = environment()
      )
    }

    # Time the evaluation
    t_par0 <- proc.time()[["elapsed"]]

    # Evaluate each replicate (parallel if cluster, serial if not)
    if (use_cluster) {
      results <- foreach::foreach(
        ii = 1:length(temp_list),
        .multicombine = TRUE,
        .inorder = TRUE
      ) %dopar% {
        # Get the input point index
        i <- temp_list[ii]

        # Set unique seed for this replicate
        set.seed(nrow(chosen_X) * length(n_replicates) + ii)

        # Evaluate f(x) with error handling
        suppressWarnings(
          suppressMessages({
            x <- chosen_X[i, , drop = FALSE]
            safe_eval_f(f, x, q)
          })
        )
      }
    } else {
      # Serial evaluation using regular for loop (avoid foreach nesting issues)
      results <- vector("list", length(temp_list))
      for (ii in seq_along(temp_list)) {
        # Get the input point index
        i <- temp_list[ii]

        # Set unique seed for this replicate
        set.seed(nrow(chosen_X) * length(n_replicates) + ii)

        # Evaluate f(x) with error handling
        results[[ii]] <- suppressWarnings(
          suppressMessages({
            x <- chosen_X[i, , drop = FALSE]
            safe_eval_f(f, x, q)
          })
        )
      }
    }

    t_par1 <- proc.time()[["elapsed"]]
    cat(sprintf(
      "foreach %%dopar%% over %d points took %.3f sec (%.2f min; %.2f hr)\n",
      nrow(chosen_X),
      (t_par1 - t_par0), (t_par1 - t_par0) / 60, (t_par1 - t_par0) / 3600
    ))

    # --------------------------------------------------------------------------
    # Step 6: Process results and update state
    # --------------------------------------------------------------------------

    # Map each result back to its corresponding row in chosen_X
    x_rows <- rep(seq_len(nrow(chosen_X)), each = n_replicates)

    for (i in 1:length(x_rows)) {
      x <- chosen_X[x_rows[i], ]
      res <- results[[i]]

      # Store the attempt in X_all
      X_all <- rbind(X_all, matrix(x, nrow = 1))
      succ_flag <- c(succ_flag, res$success)
      fail_msg <- c(fail_msg, ifelse(res$success, "", res$msg))

      # Blacklist: remove this candidate from remaining pool
      # This prevents re-evaluating the same input in future iterations
      remaining <- setdiff(remaining, chosen_idx[x_rows[i]])

      # If successful, add to success matrices
      if (res$success) {
        X_succ <- rbind(X_succ, matrix(x, nrow = 1))
        Y_succ <- rbind(Y_succ, matrix(res$y, nrow = 1))
      }
    }

    # Print progress
    if (verbose) {
      cat("   successes so far:", nrow(X_succ), "/", n_success,
          "   (remaining candidates:", length(remaining), ")\n")
    }
  }

  # ==============================================================================
  # Clean up and return results
  # ==============================================================================

  # Shut down parallel cluster (only if we created one)
  if (use_cluster) {
    parallel::stopCluster(cl)
  }

  # Return all results
  list(
    D_success = X_succ,      # Successful input points (n_success x p)
    Y_success = Y_succ,      # Successful output points (n_success x q)
    X_all = X_all,           # All attempted inputs
    success = succ_flag,     # Success indicators for all attempts
    fail_msg = fail_msg,     # Failure messages for all attempts
    remaining_idx = remaining,  # Unused candidate indices
    feas_model = feas_model  # Final feasibility model
  )
}
