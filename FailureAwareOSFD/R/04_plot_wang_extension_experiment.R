# ==============================================================================
# 04_plot_wang_extension_experiment.R
# ==============================================================================
#
# Purpose:
#   Visualize the two Wang-extension experiments used in the paper:
#
#     1. Deterministic inverse-radius reproduction
#        - establishes that our Wang-like implementation reproduces the
#          behavior of the original CRAN OSFD implementation.
#
#     2. Failure + stochastic inverse-radius extension
#        - compares input-space LHS, blind Wang-like OSFD, and failure-aware OSFD
#          when simulator calls can fail and successful outputs are noisy.
#
# Figures produced:
#
#   viz/wang_reproduction_metrics.png
#     Three-panel validation figure:
#       - output-space fill distance;
#       - mean nearest-target distance;
#       - output grid coverage.
#
#   viz/failure_stochastic_extension_metrics.png
#     Three-panel extension figure:
#       - output-space fill distance;
#       - mean nearest-target distance;
#       - attempts per success.
#
#   viz/failure_stochastic_diagnostic_scatter.png
#     Two-row diagnostic figure:
#       - attempted input locations under the true success-probability surface;
#       - successful realized output locations.
#
# Notes:
#   These figures deliberately focus on single-point sequential design.
#   They do not compare batching and do not compare replication.
#
# Author: AC Murph
# Date: Sep 2026
# ==============================================================================

# ------------------------------------------------------------------------------
# Load required packages
# ------------------------------------------------------------------------------

suppressPackageStartupMessages({
  library(dplyr)      # Data manipulation
  library(tidyr)      # Data tidying
  library(purrr)      # Functional programming
  library(tibble)     # Modern data frames
  library(ggplot2)    # Plotting
  library(patchwork)  # Plot composition
  library(here)       # Path management
})

# ------------------------------------------------------------------------------
# Load helper functions
# ------------------------------------------------------------------------------

# Load plotting helper functions from wang_helper_functions.R.
# Expected helper functions include:
#   - plot_metric_by_n()
#   - plot_metric_by_attempts()
#   - clean_method_names()
source(here::here("R", "source_helpers.R"))
source_helpers()

# ------------------------------------------------------------------------------
# Define directories
# ------------------------------------------------------------------------------

OUT_DIR <- here::here("data", "wang_extension")
FIG_DIR <- here::here("viz")
dir.create(FIG_DIR, recursive = TRUE, showWarnings = FALSE)

# ------------------------------------------------------------------------------
# Load experiment results
# ------------------------------------------------------------------------------

# Script 01: deterministic reproduction of Wang inverse-radius example.
repro <- readRDS(file.path(OUT_DIR, "01_reproduce_wang_inverse_radius.rds"))

# Script 03: failure + stochastic inverse-radius extension.
stoch <- readRDS(file.path(OUT_DIR, "03_failure_stochastic_inverse_radius.rds"))

repro_metrics <- repro$metrics
stoch_metrics <- stoch$metrics

# ------------------------------------------------------------------------------
# Optional method ordering for diagnostic panels
# ------------------------------------------------------------------------------

method_order <- c(
  "Basic LHS",
  "Our Wang-like OSFD",
  "Failure-aware OSFD"
)

# ==============================================================================
# Print timing summary for all methods and replicates
# ==============================================================================

cat("\n=== Timing Summary: Deterministic Wang Reproduction ===\n\n")

# Calculate timing statistics by method and n_target
timing_summary <- repro_metrics %>%
  group_by(method, n_target) %>%
  summarise(
    mean_time_sec = mean(elapsed_time_sec, na.rm = TRUE),
    sd_time_sec = sd(elapsed_time_sec, na.rm = TRUE),
    min_time_sec = min(elapsed_time_sec, na.rm = TRUE),
    max_time_sec = max(elapsed_time_sec, na.rm = TRUE),
    n_reps = n(),
    .groups = "drop"
  ) %>%
  mutate(
    mean_time_min = mean_time_sec / 60,
    sd_time_min = sd_time_sec / 60
  )

# Print detailed timing by method and n
cat("Timing by method and n_target (mean ± sd across replicates):\n")
timing_summary %>%
  arrange(method, n_target) %>%
  mutate(
    timing_str = sprintf("%.2f ± %.2f sec (%.3f ± %.3f min)",
                         mean_time_sec, sd_time_sec,
                         mean_time_min, sd_time_min)
  ) %>%
  select(method, n_target, timing_str, n_reps) %>%
  print(n = Inf)

# Print timing for each individual replicate
cat("\n\nDetailed timing for all methods and replicates:\n")
repro_metrics %>%
  select(method, rep_id, n_target, elapsed_time_sec) %>%
  arrange(method, n_target, rep_id) %>%
  mutate(elapsed_time_min = elapsed_time_sec / 60) %>%
  print(n = Inf)

cat("\n==========================================================\n\n")

# ==============================================================================
# Figure 1: Deterministic Wang reproduction
# ==============================================================================

# This figure establishes that our Wang-like implementation agrees with the
# original CRAN OSFD implementation when failures, stochasticity, batching,
# replication, and feasibility weighting are not part of the problem.

# Panel 1: Fill distance.
#
# Fill distance is the worst-case nearest-neighbor distance from the reference
# output set to the successful design outputs. Lower is better.
p_repro_fill <- plot_metric_by_n(
  repro_metrics,
  metric = "fill_distance",
  ylab = "Estimated output-space fill distance",
  title = "Deterministic inverse-radius reproduction"
)

# Panel 2: Mean nearest-target distance.
#
# Mean distance measures typical closeness to the reference output set.
# Lower is better.
p_repro_mean <- plot_metric_by_n(
  repro_metrics,
  metric = "mean_target_distance",
  ylab = "Mean nearest-target distance",
  title = "Mean target distance"
)

# Panel 3: Output grid coverage.
#
# Grid coverage is an output-space occupancy metric. Higher is better, but this
# metric is less directly aligned with the OSFD objective than fill distance.
p_repro_grid <- plot_metric_by_n(
  repro_metrics,
  metric = "grid_coverage",
  ylab = "Grid coverage",
  title = "Output grid coverage",
  lower_is_better = FALSE
)

# Panel 4: Timing (elapsed time in seconds).
#
# Shows wall-clock runtime for each method across sample sizes.
# Lower is better (faster).
p_repro_timing <- plot_metric_by_n(
  repro_metrics,
  metric = "elapsed_time_sec",
  ylab = "Elapsed time (seconds)",
  title = "Computational time",
  lower_is_better = TRUE
)

fig_repro <- (p_repro_fill | p_repro_mean | p_repro_grid | p_repro_timing) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  filename = file.path(FIG_DIR, "wang_reproduction_metrics.png"),
  plot = fig_repro,
  width = 20,
  height = 5,
  dpi = 300
)

# ==============================================================================
# Print timing summary for failure + stochastic experiment
# ==============================================================================

cat("\n=== Timing Summary: Failure + Stochastic Extension ===\n\n")

# Calculate timing statistics by method and max_attempts
timing_summary_stoch <- stoch_metrics %>%
  group_by(method, max_attempts) %>%
  summarise(
    mean_time_sec = mean(elapsed_time_sec, na.rm = TRUE),
    sd_time_sec = sd(elapsed_time_sec, na.rm = TRUE),
    min_time_sec = min(elapsed_time_sec, na.rm = TRUE),
    max_time_sec = max(elapsed_time_sec, na.rm = TRUE),
    n_reps = n(),
    .groups = "drop"
  ) %>%
  mutate(
    mean_time_min = mean_time_sec / 60,
    sd_time_min = sd_time_sec / 60
  )

# Print detailed timing by method and max_attempts
cat("Timing by method and max_attempts (mean ± sd across replicates):\n")
timing_summary_stoch %>%
  arrange(method, max_attempts) %>%
  mutate(
    timing_str = sprintf("%.2f ± %.2f sec (%.3f ± %.3f min)",
                         mean_time_sec, sd_time_sec,
                         mean_time_min, sd_time_min)
  ) %>%
  select(method, max_attempts, timing_str, n_reps) %>%
  print(n = Inf)

# Print timing for each individual replicate
cat("\n\nDetailed timing for all methods and replicates:\n")
stoch_metrics %>%
  select(method, rep_id, max_attempts, elapsed_time_sec) %>%
  arrange(method, max_attempts, rep_id) %>%
  mutate(elapsed_time_min = elapsed_time_sec / 60) %>%
  print(n = Inf)

cat("\n==========================================================\n\n")

# ==============================================================================
# Figure 2: Failure + stochastic extension metrics
# ==============================================================================

# This figure is the main extension comparison. It asks how the methods behave
# when simulator calls can fail and successful outputs are stochastic.

# Panel 1: Fill distance vs. max_attempts.
#
# NOTE: Script 03 uses ATTEMPT BUDGETS. All methods get the same number of
# attempts (max_attempts), and we compare their output space filling performance.
# X-axis is the attempt budget, not number of successes.
p_stoch_fill <- plot_metric_by_n(
  stoch_metrics %>% rename(n_target = max_attempts),  # Rename for compatibility with plot function
  metric = "fill_distance",
  ylab = "Estimated output-space fill distance",
  title = "Failure + stochastic extension:\nfill distance"
)

# Panel 2: Mean nearest-target distance vs. max_attempts.
#
# This metric emphasizes typical output locations and can favor methods that
# sample common regions of the stochastic output distribution well.
p_stoch_mean <- plot_metric_by_n(
  stoch_metrics %>% rename(n_target = max_attempts),  # Rename for compatibility with plot function
  metric = "mean_target_distance",
  ylab = "Mean nearest-target distance",
  title = "Failure + stochastic extension:\nmean distance"
)

# Panel 3: Number of successes vs. max_attempts (higher is better, more efficient)
# This shows how many successes each method achieves with the same attempt budget
#
# NOTE: Script 03 uses ATTEMPT BUDGETS. The x-axis is max_attempts (the budget),
# and we plot how many successes each method achieves.
p_stoch_eff <- plot_metric_by_n(
  stoch_metrics %>% rename(n_target = max_attempts),  # Rename for compatibility with plot function
  metric = "successes",
  ylab = "Number of successes achieved",
  title = "Failure + stochastic extension:\nsuccess efficiency",
  lower_is_better = FALSE  # Higher successes = better
)

fig_stoch <- (p_stoch_fill | p_stoch_mean | p_stoch_eff) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

ggsave(
  filename = file.path(FIG_DIR, "failure_stochastic_extension_metrics.png"),
  plot = fig_stoch,
  width = 16,
  height = 5,
  dpi = 300
)

# ==============================================================================
# Figure 3: Diagnostic scatter plots for failure + stochastic experiment
# ==============================================================================

# This figure uses the representative replicate saved from script 03.
#
# Top row:
#   Attempted input locations in x-space. Points are colored by the true
#   success probability at that input and shaped by whether the simulator call
#   succeeded.
#
# Bottom row:
#   Successful realized output locations in y-space. Failed calls are not shown
#   in the output row because they do not return usable simulator outputs.
#
# Because this is the stochastic-output experiment, the bottom row shows realized
# noisy outputs, not latent deterministic inverse-radius outputs.

failure_settings <- stoch$settings
# Filter diagnostics to only use max_attempts == 500
# The experiment script (03) saves diagnostics for the first replicate at n_attempts = 500
diagnostics <- stoch$diagnostics %>%
  purrr::keep(~{
    # Check if max_attempts is stored in each diagnostic entry
    if (!is.null(.x$max_attempts)) {
      .x$max_attempts == 500
    } else {
      # Fallback: keep all if max_attempts field doesn't exist (backward compatibility)
      TRUE
    }
  })

# If no entries match, use all diagnostics (backward compatibility with older saved results)
if (length(diagnostics) == 0) {
  diagnostics <- stoch$diagnostics
  warning("No diagnostics found with max_attempts == 500. Using all available diagnostics.")
}

# ------------------------------------------------------------------------------
# Robust settings helper
# ------------------------------------------------------------------------------

# Some saved experiment versions use a soft-only failure surface with r0_fail.
# Newer versions use a hard-failure core plus a soft transition, with
# r_hard_fail and r_soft_fail. This helper supports both.
get_setting <- function(settings, names, default = NULL) {
  for (nm in names) {
    if (!is.null(settings[[nm]])) {
      return(settings[[nm]])
    }
  }
  default
}

p_min_true <- get_setting(
  failure_settings,
  c("p_min_true", "P_MIN_TRUE"),
  default = 0.05
)

a_fail <- get_setting(
  failure_settings,
  c("a_fail", "A_FAIL"),
  default = 35
)

r_hard_fail <- get_setting(
  failure_settings,
  c("r_hard_fail", "R_HARD_FAIL"),
  default = NULL
)

r_soft_fail <- get_setting(
  failure_settings,
  c("r_soft_fail", "R_SOFT_FAIL"),
  default = NULL
)

r0_fail <- get_setting(
  failure_settings,
  c("r0_fail", "R0_FAIL"),
  default = NULL
)

# ------------------------------------------------------------------------------
# True success-probability surface
# ------------------------------------------------------------------------------

#' Compute true success probability at input location(s)
#'
#' This function is vectorized so that it can be used inside dplyr::mutate().
#'
#' @param x1 First input dimension.
#' @param x2 Second input dimension.
#' @return Success probability in [0, 1].
p_success_surface <- function(x1, x2) {
  r <- sqrt(x1^2 + x2^2)

  # Preferred current form: hard-failure core plus soft transition.
  if (!is.null(r_hard_fail) && !is.null(r_soft_fail)) {
    soft_prob <-
      p_min_true +
      (1 - p_min_true) *
      plogis(a_fail * (r - r_soft_fail))

    return(ifelse(r <= r_hard_fail, 0, soft_prob))
  }

  # Backward-compatible form: soft logistic failure surface only.
  if (!is.null(r0_fail)) {
    return(
      p_min_true +
        (1 - p_min_true) *
        plogis(a_fail * (r - r0_fail))
    )
  }

  stop("Could not identify failure-surface parameters in stoch$settings.")
}

# Generate grid of success probabilities for the input-space background.
surface_df <- expand.grid(
  x1 = seq(0, 1, length.out = 200),
  x2 = seq(0, 1, length.out = 200)
) %>%
  mutate(p_success = p_success_surface(x1, x2))

# ------------------------------------------------------------------------------
# Diagnostic data
# ------------------------------------------------------------------------------

# All attempted input locations.
diag_input_df <- purrr::map_dfr(diagnostics, function(z) {
  tibble(
    method = clean_method_names(z$method),
    x1 = z$X_all[, 1],
    x2 = z$X_all[, 2],
    success = as.logical(z$success)
  ) %>%
    mutate(
      p_success = p_success_surface(x1, x2),
      method = factor(method, levels = method_order),
      success = factor(success, levels = c(FALSE, TRUE))
    )
})

# Successful realized output locations.
diag_output_df <- purrr::map_dfr(diagnostics, function(z) {
  tibble(
    method = clean_method_names(z$method),
    y1 = z$Y[, 1],
    y2 = z$Y[, 2]
  ) %>%
    mutate(method = factor(method, levels = method_order))
})

# ------------------------------------------------------------------------------
# Top row: attempted input locations
# ------------------------------------------------------------------------------

p_attempts <- ggplot() +
  geom_raster(
    data = surface_df,
    aes(x = x1, y = x2, fill = p_success),
    alpha = 0.25
  ) +
  geom_point(
    data = diag_input_df,
    aes(x = x1, y = x2, color = p_success, shape = success),
    size = 1.5,
    alpha = 0.85
  ) +
  facet_wrap(~ method, nrow = 1) +
  scale_fill_gradient(
    name = "True success\nprobability",
    low = "#2166AC",
    high = "#F4A582"
  ) +
  scale_color_gradient(
    name = "True success\nprobability",
    low = "#2166AC",
    high = "#F4A582"
  ) +
  scale_shape_manual(
    name = "Succeeded",
    values = c("FALSE" = 1, "TRUE" = 16),  # 1 = open circle, 16 = closed circle
    labels = c("FALSE" = "Failed", "TRUE" = "Completed")
  ) +
  labs(
    title = "Attempted input locations under simulated failures and stochastic outputs",
    x = expression(x[1]),
    y = expression(x[2])
  ) +
  guides(
    fill = "none",
    color = guide_colorbar(order = 1),
    shape = guide_legend(order = 2)
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold"),
    legend.position = "right"
  )

# ------------------------------------------------------------------------------
# Bottom row: successful realized output locations
# ------------------------------------------------------------------------------

p_outputs <- ggplot(diag_output_df, aes(x = y1, y = y2)) +
  geom_point(size = 1, alpha = 0.7) +
  facet_wrap(~ method, nrow = 1) +
  labs(
    title = "Successful realized output locations under simulated failures and stochastic outputs",
    x = expression(y[1]),
    y = expression(y[2])
  ) +
  theme_bw(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold")
  )

fig_diag <- p_attempts / p_outputs

ggsave(
  filename = file.path(FIG_DIR, "failure_stochastic_diagnostic_scatter.png"),
  plot = fig_diag,
  width = 16,
  height = 9,
  dpi = 300
)

# ------------------------------------------------------------------------------
# Print confirmation messages
# ------------------------------------------------------------------------------

message("Saved figures:")
message("  viz/wang_reproduction_metrics.png")
message("  viz/failure_stochastic_extension_metrics.png")
message("  viz/failure_stochastic_diagnostic_scatter.png")