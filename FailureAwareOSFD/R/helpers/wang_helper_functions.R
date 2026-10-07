# Wang Helper Functions
#
# Purpose:
#   Shared helper functions for Wang OSFD extension experiments (scripts 01-04)
#   These functions support:
#     - Reference grid generation
#     - Matrix evaluation
#     - Output space filling metrics (scaling, distances, grid coverage)
#     - Safe function evaluation
#     - LHS-based success budget sampling
#
# Author: AC Murph
# Date: Sep 2026

# =============================================================================
# Test Functions
# =============================================================================

#' Inverse radius test function
#'
#' Computes the inverse of the Euclidean distance from the origin, plus
#' the angle in polar coordinates.
#'
#' @param x Numeric vector of length 2 (input point in [0,1]^2)
#' @param epsilon Small constant added to denominator to avoid division by zero
#' @return Numeric vector of length 2: (inverse_radius, angle)
inverse_radius <- function(x, epsilon = 0.1) {
  # Ensure x is a numeric vector
  x <- as.numeric(x)

  # y1: inverse of distance from origin (with epsilon to prevent div/0)
  y1 <- 1 / sqrt(x[1]^2 + x[2]^2 + epsilon^2)

  # y2: angle in polar coordinates (atan2 handles x1=0 cleanly)
  y2 <- atan2(x[2], x[1])

  # Return 2D output
  c(y1, y2)
}

#' Success probability for failure experiments
#'
#' Defines a spatially-varying success probability with a hard-failure core
#' near the origin. Inputs within r <= r_hard_fail ALWAYS fail (simulator crashes).
#' Outside the core, success probability increases with distance via logistic transition.
#'
#' @param x Numeric vector of length 2 (input point)
#' @param p_min_true Minimum success probability outside hard-failure core (default 0.05)
#' @param r_hard_fail Radius of hard-failure core where p_success = 0 (default 0.08)
#' @param r_soft_fail Midpoint of soft transition where p_success ≈ 50% (default 0.22)
#' @param a_fail Steepness of logistic transition (default 35)
#' @return Scalar success probability in [0, 1]
p_success <- function(x, p_min_true = 0.05, r_hard_fail = 0.08, r_soft_fail = 0.22, a_fail = 35) {
  # Ensure x is numeric
  x <- as.numeric(x)

  # Compute distance from origin
  r <- sqrt(sum(x^2))

  # Hard simulator-failure region: inputs with r <= r_hard_fail always fail
  # These inputs never return usable outputs (simulator crashes)
  if (r <= r_hard_fail) {
    return(0)
  }

  # Outside the hard-failure core, success probability increases with distance
  # from the origin. The region just outside the hard core is still difficult,
  # but not impossible.
  p_min_true + (1 - p_min_true) * plogis(a_fail * (r - r_soft_fail))
}

#' Inverse radius with simulated failures
#'
#' Wraps inverse_radius with a stochastic failure mechanism.
#' Throws an error if the evaluation fails.
#'
#' @param x Input point
#' @param epsilon Epsilon parameter for inverse_radius
#' @param p_min_true Minimum success probability outside hard-failure core
#' @param r_hard_fail Radius of hard-failure core (always fails)
#' @param r_soft_fail Midpoint of soft transition
#' @param a_fail Logistic steepness parameter
#' @return Output vector if successful, throws error if failed
failure_inverse_radius <- function(x, epsilon = 0.1,
                                   p_min_true = 0.05, r_hard_fail = 0.08,
                                   r_soft_fail = 0.22, a_fail = 35) {
  # Compute success probability for this input
  ps <- p_success(x, p_min_true, r_hard_fail, r_soft_fail, a_fail)

  # Random Bernoulli trial: fail with probability (1 - ps)
  if (runif(1) > ps) {
    stop("simulated failure")
  }

  # If successful, return the output
  inverse_radius(x, epsilon)
}

#' Inverse radius with failures and stochastic output noise
#'
#' Adds independent Gaussian noise to each output dimension.
#'
#' @param x Input point
#' @param epsilon Epsilon parameter for inverse_radius
#' @param p_min_true Minimum success probability outside hard-failure core
#' @param r_hard_fail Radius of hard-failure core (always fails)
#' @param r_soft_fail Midpoint of soft transition
#' @param a_fail Logistic steepness parameter
#' @param sigma_y1 Standard deviation of noise for first output (default 0.15)
#' @param sigma_y2 Standard deviation of noise for second output (default 0.03)
#' @return Noisy output vector if successful, throws error if failed
failure_stochastic_inverse_radius <- function(x, epsilon = 0.1,
                                              p_min_true = 0.05, r_hard_fail = 0.08,
                                              r_soft_fail = 0.22, a_fail = 35,
                                              sigma_y1 = 0.15, sigma_y2 = 0.03) {
  # Check for failure first
  ps <- p_success(x, p_min_true, r_hard_fail, r_soft_fail, a_fail)
  if (runif(1) > ps) {
    stop("simulated failure")
  }

  # Compute deterministic output
  y <- inverse_radius(x, epsilon)

  # Add independent Gaussian noise to each output dimension
  y + c(
    rnorm(1, mean = 0, sd = sigma_y1),
    rnorm(1, mean = 0, sd = sigma_y2)
  )
}

# =============================================================================
# Reference Grid and Evaluation
# =============================================================================

#' Generate a dense regular grid of inputs for reference evaluation
#'
#' Creates a uniform grid over [0,1]^2 for computing reference outputs
#' and evaluating space-filling metrics.
#'
#' @param side Number of points along each dimension (default 300)
#' @param p Input dimension (default 2)
#' @return Matrix of dimension (side^p) x p
make_reference_inputs <- function(side = 300, p = 2) {
  # Create uniform grid over [0,1]^p
  if (p == 2) {
    grid <- expand.grid(
      x1 = seq(0, 1, length.out = side),
      x2 = seq(0, 1, length.out = side)
    )
  } else {
    stop("make_reference_inputs currently only supports p=2")
  }

  # Return as matrix
  as.matrix(grid)
}

#' Evaluate a function f over all rows of a matrix X
#'
#' Applies function f to each row of X and returns a matrix of outputs.
#'
#' @param X Input matrix (n x p)
#' @param f Function that takes a vector of length p and returns vector of length q
#' @return Output matrix (n x q) with column names y1, y2, ..., yq
eval_matrix <- function(X, f) {
  # Apply f to each row, transpose to get n x q matrix
  Y <- t(apply(X, 1, f))

  # Add column names
  colnames(Y) <- paste0("y", seq_len(ncol(Y)))

  Y
}

#' Generate stochastic reference set by evaluating failure-prone function
#'
#' For stochastic experiments, we need a reference set drawn from the
#' actual stochastic output distribution (not the deterministic version).
#'
#' @param n_ref Number of reference points to attempt (default 75000)
#' @param seed Random seed for reproducibility
#' @param f Stochastic function to evaluate (should throw error on failure)
#' @param p Input dimension (default 2)
#' @param q Output dimension (default 2)
#' @return List with X (successful inputs) and Y (successful outputs)
eval_reference_stochastic <- function(n_ref = 75000, seed = 20260902 + 999,
                                     f = failure_stochastic_inverse_radius,
                                     p = 2, q = 2) {
  set.seed(seed)

  # Generate LHS sample of candidate inputs
  X <- lhs::randomLHS(n_ref, p)

  # Initialize storage for successful evaluations
  Y <- matrix(numeric(0), ncol = q)
  X_success <- matrix(numeric(0), ncol = p)

  # Evaluate each input, catching failures
  for (i in seq_len(nrow(X))) {
    res <- tryCatch(
      {
        y <- f(X[i, ])
        list(success = TRUE, y = y)
      },
      error = function(e) {
        list(success = FALSE, y = rep(NA_real_, q))
      }
    )

    # Only keep finite successful outputs
    if (res$success && all(is.finite(res$y))) {
      X_success <- rbind(X_success, matrix(X[i, ], nrow = 1))
      Y <- rbind(Y, matrix(res$y, nrow = 1))
    }
  }

  # Return successful inputs and outputs
  list(X = X_success, Y = Y)
}

# =============================================================================
# Output Space Scaling and Distance Metrics
# =============================================================================

#' Scale outputs Y relative to reference outputs Y_ref
#'
#' Transforms Y to [0,1]^q based on the min/max of Y_ref in each dimension.
#' This normalization is essential for comparing distances across dimensions.
#'
#' @param Y Output matrix to scale (n x q)
#' @param Y_ref Reference outputs defining scale (m x q)
#' @return Scaled matrix (n x q) where each column is in [0,1] based on Y_ref range
scale_to_reference <- function(Y, Y_ref) {
  # Compute min and max for each output dimension from reference
  mins <- apply(Y_ref, 2, min, na.rm = TRUE)
  maxs <- apply(Y_ref, 2, max, na.rm = TRUE)

  # Compute range, with machine epsilon floor to avoid division by zero
  den <- pmax(maxs - mins, .Machine$double.eps)

  # Center by min, then scale by range: (Y - min) / (max - min)
  sweep(sweep(Y, 2, mins, "-"), 2, den, "/")
}

#' Compute nearest neighbor distances from reference to design points
#'
#' For each reference point, finds the distance to the nearest design point.
#' Used to quantify how well the design fills the output space.
#'
#' @param Y_ref_scaled Scaled reference outputs (m x q)
#' @param Y_design_scaled Scaled design outputs (n x q)
#' @return Vector of length m with nearest neighbor distances
nearest_distances <- function(Y_ref_scaled, Y_design_scaled) {
  # Use FNN package for efficient k-nearest neighbors (k=1)
  nn <- FNN::get.knnx(
    data = Y_design_scaled,    # Search among design points
    query = Y_ref_scaled,       # Query from reference points
    k = 1                       # Find 1 nearest neighbor
  )

  # Extract distances (column 1 since k=1)
  as.numeric(nn$nn.dist[, 1])
}

#' Map scaled outputs to discrete grid cells
#'
#' Discretizes the output space into a k x k grid and assigns each point
#' to a cell index. Used for grid coverage metric.
#'
#' @param Y_scaled Scaled outputs in [0,1]^2 (n x 2)
#' @param k Grid resolution (default 100, giving 10,000 cells)
#' @return Vector of cell indices (1 to k^2)
grid_cells <- function(Y_scaled, k = 100) {
  # Convert to matrix if needed and ensure proper dimensions
  Y_scaled <- as.matrix(Y_scaled)

  # Handle case where it's been converted to a vector (single point case)
  # Check if we lost the 2-column structure
  if (ncol(Y_scaled) != 2) {
    # If it's a vector of length 2, make it a 1x2 matrix
    if (length(Y_scaled) == 2) {
      Y_scaled <- matrix(Y_scaled, nrow = 1, ncol = 2)
    } else {
      stop("Y_scaled must have 2 columns (Q=2)")
    }
  }

  # Extract columns explicitly to prevent dimension dropping
  # Using drop=TRUE converts to vector, which is safe for pmin/pmax
  col1 <- Y_scaled[, 1, drop = TRUE]
  col2 <- Y_scaled[, 2, drop = TRUE]

  # Clamp each column to [0,1] to handle numerical issues
  col1 <- pmin(1, pmax(0, col1))
  col2 <- pmin(1, pmax(0, col2))

  # Map each dimension to grid index [1, k]
  # floor(col * k) + 1 maps [0,1) to [1,k], with col=1 clamped to k
  ix <- pmin(k, pmax(1, floor(col1 * k) + 1))
  iy <- pmin(k, pmax(1, floor(col2 * k) + 1))

  # Convert 2D grid index to 1D cell index (column-major ordering)
  ix + k * (iy - 1)
}

# =============================================================================
# Composite Metrics
# =============================================================================

#' Compute comprehensive output space-filling metrics
#'
#' Evaluates how well a design Y_design fills the output space defined by Y_ref.
#' Metrics include fill distance, mean/median/quantile distances, and grid coverage.
#'
#' @param Y_design Design outputs (n x q)
#' @param Y_ref Reference outputs defining target space (m x q)
#' @param k Grid resolution for coverage metric (default 100)
#' @return Tibble with metrics: fill_distance, mean/median/q90/q95_target_distance,
#'         grid_coverage, n_design_outputs
compute_metrics <- function(Y_design, Y_ref, k = 100) {
  # Scale both design and reference to [0,1]^q based on reference range
  Y_ref_scaled <- scale_to_reference(Y_ref, Y_ref)
  Y_design_scaled <- scale_to_reference(Y_design, Y_ref)

  # Compute nearest neighbor distances: for each ref point, distance to nearest design point
  d <- nearest_distances(Y_ref_scaled, Y_design_scaled)

  # Compute grid cells covered by reference and design
  ref_cells <- unique(grid_cells(Y_ref_scaled, k = k))
  design_cells <- unique(grid_cells(Y_design_scaled, k = k))

  # Return metrics as tibble
  tibble::tibble(
    fill_distance = max(d),                                    # Worst-case coverage (lower is better)
    mean_target_distance = mean(d),                            # Average coverage (lower is better)
    median_target_distance = median(d),                        # Median coverage (lower is better)
    q90_target_distance = as.numeric(quantile(d, 0.90)),      # 90th percentile (lower is better)
    q95_target_distance = as.numeric(quantile(d, 0.95)),      # 95th percentile (lower is better)
    grid_coverage = length(intersect(design_cells, ref_cells)) / length(ref_cells),  # Proportion of ref grid covered
    n_design_outputs = nrow(Y_design)                          # Number of design points
  )
}

# =============================================================================
# Safe Function Evaluation
# =============================================================================

#' Safely evaluate a potentially-failing function once
#'
#' Wraps function evaluation in tryCatch to handle errors gracefully.
#' Returns success status, output (or NA if failed), and error message.
#'
#' @param f Function to evaluate
#' @param x Input point
#' @param q Expected output dimension (default 2)
#' @return List with success (logical), y (numeric vector), msg (character)
safe_eval_once <- function(f, x, q = 2) {
  tryCatch(
    {
      # Attempt evaluation
      y <- f(x)

      # Validate output: must be numeric, correct length, and finite
      if (!is.numeric(y) || length(y) != q || any(!is.finite(y))) {
        list(success = FALSE, y = rep(NA_real_, q), msg = "invalid output")
      } else {
        list(success = TRUE, y = as.numeric(y), msg = "")
      }
    },
    error = function(e) {
      # Catch any error and return failure with error message
      list(success = FALSE, y = rep(NA_real_, q), msg = conditionMessage(e))
    }
  )
}

# =============================================================================
# LHS with Success Budget
# =============================================================================

#' Run LHS until reaching a target number of successful evaluations
#'
#' Generates a large LHS stream and evaluates points sequentially until
#' n_success successful outputs are obtained. Useful for failure-prone functions.
#'
#' @param n_success Target number of successful evaluations
#' @param seed Random seed
#' @param f Function to evaluate (may throw errors)
#' @param p Input dimension (default 2)
#' @param q Output dimension (default 2)
#' @param max_attempts Maximum total evaluations to prevent infinite loops (default 100000)
#' @return List with D (successful inputs), Y (successful outputs),
#'         attempts, successes, success_rate, X_all, success vector
run_lhs_success_budget <- function(n_success, seed, f, p = 2, q = 2, max_attempts = 100000) {
  set.seed(seed)

  # Generate a large LHS stream upfront
  X_stream <- lhs::randomLHS(max_attempts, p)

  # Initialize storage
  X_all <- matrix(numeric(0), ncol = p)      # All attempted inputs
  success <- logical(0)                       # Success status for each attempt
  Y_success <- matrix(numeric(0), ncol = q)   # Successful outputs
  X_success <- matrix(numeric(0), ncol = p)   # Successful inputs

  # Iterate through stream until success budget is met
  for (i in seq_len(nrow(X_stream))) {
    x <- X_stream[i, ]

    # Safely evaluate
    res <- safe_eval_once(f, x, q)

    # Record attempt
    X_all <- rbind(X_all, matrix(x, nrow = 1))
    success <- c(success, res$success)

    # Record success if applicable
    if (res$success) {
      X_success <- rbind(X_success, matrix(x, nrow = 1))
      Y_success <- rbind(Y_success, matrix(res$y, nrow = 1))
    }

    # Check if target reached
    if (nrow(Y_success) >= n_success) {
      break
    }
  }

  # Warn if we exhausted the stream without reaching target
  if (nrow(Y_success) < n_success) {
    warning("LHS stream ended before reaching n_success.")
  }

  # Return complete evaluation record
  list(
    D = X_success,
    Y = Y_success,
    attempts = nrow(X_all),
    successes = nrow(Y_success),
    success_rate = nrow(Y_success) / nrow(X_all),
    X_all = X_all,
    success = success
  )
}

# =============================================================================
# Plotting Helper Functions
# =============================================================================

#' Clean method names for publication-quality labels
#'
#' Maps short method identifiers to human-readable names for plot legends.
#' Names are chosen to be consistent with other plots in the paper.
#'
#' @param x Character vector of method names
#' @return Character vector of cleaned names
clean_method_names <- function(x) {
  method_labels <- c(
    LHS = "Basic LHS",                             # Input-space Latin Hypercube Sampling
    CRAN_OSFD = "CRAN OSFD",                       # Wang's published implementation
    Ours_Wang_like = "Our Wang-like OSFD",         # Our implementation matching Wang (deterministic)
    Blind_Wang_like_OSFD = "Our Wang-like OSFD",   # Also maps to Our Wang-like (stochastic experiments)
    Failure_aware_OSFD = "Failure-aware OSFD"      # Our main method (with feasibility learning)
  )
  dplyr::recode(x, !!!method_labels)
}

#' Summarize metric across replicates for plotting
#'
#' Computes mean and 5th/95th percentile bands (90% interval) for a given
#' metric across replicates at each target size n.
#'
#' @param df Metrics data frame
#' @param metric Name of metric column to summarize
#' @return Tibble with method, n_target, mean, lo, hi, method_label
summarise_metric <- function(df, metric) {
  df %>%
    dplyr::group_by(method, n_target) %>%
    dplyr::summarise(
      mean = mean(.data[[metric]], na.rm = TRUE),     # Mean across replicates
      lo = as.numeric(quantile(.data[[metric]], 0.05, na.rm = TRUE)),  # 5th percentile
      hi = as.numeric(quantile(.data[[metric]], 0.95, na.rm = TRUE)),  # 95th percentile
      .groups = "drop"
    ) %>%
    dplyr::mutate(method_label = clean_method_names(method))  # Add human-readable labels
}

#' Plot metric vs. target output size n
#'
#' Creates a line plot with ribbon showing metric performance as a function
#' of target output size. Used for deterministic experiments where n fully
#' determines computational cost.
#'
#' @param df Metrics data frame
#' @param metric Name of metric column to plot
#' @param ylab Y-axis label
#' @param title Plot title
#' @param lower_is_better Logical: is lower metric value better? (not used currently)
#' @return ggplot object
plot_metric_by_n <- function(df, metric, ylab, title, lower_is_better = TRUE) {
  # Summarize metric across replicates
  sm <- summarise_metric(df, metric)

  # Define consistent color mapping across all plots
  # Using viridis palette with 5 distinct colors for all methods
  viridis_colors <- scales::viridis_pal(option = "D")(5)
  method_colors <- c(
    "Basic LHS" = viridis_colors[1],                 # Purple/dark blue (color 1)
    "Failure-aware OSFD" = viridis_colors[2],        # Blue/cyan (color 2)
    "Inverse SIR Maps" = viridis_colors[3],          # Teal/green (color 3)
    "Our Wang-like OSFD" = viridis_colors[4],        # Yellow/green (color 4)
    "CRAN OSFD" = viridis_colors[5]                  # Yellow (color 5)
  )

  # Create line plot with ribbon for uncertainty bands
  # Use manual color scheme for consistency across all plots in the paper
  ggplot2::ggplot(sm, ggplot2::aes(x = n_target, y = mean, color = method_label, group = method_label)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi, fill = method_label), alpha = 0.25, color = NA) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_color_manual(values = method_colors, name = "Method") +
    ggplot2::scale_fill_manual(values = method_colors, name = "Method") +
    ggplot2::labs(
      title = title,
      x = "Simulator attempt budget",   # X-axis: target number of successes
      y = ylab
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )
}

#' Plot metric vs. mean total attempts
#'
#' Creates a line plot showing metric performance as a function of
#' mean simulator attempts. Used for failure-prone experiments where
#' different methods require different numbers of attempts to reach
#' the same number of successes.
#'
#' @param df Metrics data frame
#' @param metric Name of metric column to plot
#' @param ylab Y-axis label
#' @param title Plot title
#' @return ggplot object
plot_metric_by_attempts <- function(df, metric, ylab, title) {
  # Summarize metric and attempts across replicates
  sm <- df %>%
    dplyr::group_by(method, n_target) %>%
    dplyr::summarise(
      attempts = mean(attempts, na.rm = TRUE),  # Mean attempts to reach n_target successes
      mean = mean(.data[[metric]], na.rm = TRUE),
      lo = as.numeric(quantile(.data[[metric]], 0.05, na.rm = TRUE)),  # 5th percentile
      hi = as.numeric(quantile(.data[[metric]], 0.95, na.rm = TRUE)),  # 95th percentile
      .groups = "drop"
    ) %>%
    dplyr::mutate(method_label = clean_method_names(method))

  # Define consistent color mapping across all plots
  # Using viridis palette with 5 distinct colors for all methods
  viridis_colors <- scales::viridis_pal(option = "D")(5)
  method_colors <- c(
    "Basic LHS" = viridis_colors[1],                 # Purple/dark blue (color 1)
    "Failure-aware OSFD" = viridis_colors[2],        # Blue/cyan (color 2)
    "Inverse SIR Maps" = viridis_colors[3],          # Teal/green (color 3)
    "Our Wang-like OSFD" = viridis_colors[4],        # Yellow/green (color 4)
    "CRAN OSFD" = viridis_colors[5]                  # Yellow (color 5)
  )

  # Create line plot with X-axis = attempts (not n_target)
  # Use manual color scheme for consistency across all plots in the paper
  ggplot2::ggplot(sm, ggplot2::aes(x = attempts, y = mean, color = method_label, group = method_label)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lo, ymax = hi, fill = method_label), alpha = 0.25, color = NA) +
    ggplot2::geom_line(linewidth = 1.2) +
    ggplot2::geom_point(size = 3) +
    ggplot2::scale_color_manual(values = method_colors, name = "Method") +
    ggplot2::scale_fill_manual(values = method_colors, name = "Method") +
    ggplot2::labs(
      title = title,
      x = "Mean total simulator attempts",  # X-axis: computational cost
      y = ylab
    ) +
    ggplot2::theme_bw(base_size = 14) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      legend.position = "bottom"
    )
}
