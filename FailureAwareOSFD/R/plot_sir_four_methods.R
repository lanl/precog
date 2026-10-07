# Script to compare all four methods (SIR_Maps, Basic_LHS, OSFD_Basic, OSFD_WangOrig)
# Creates separate plots for n=1000-10000
# Author: AC Murph
library(ggplot2)
library(patchwork)
library(tidyr)
library(dplyr)
library(latex2exp)
setwd(here::here())
source(here::here("R", "sir_experiment_setup.R"))

alpha_lo <- alpha_bounds[1];  alpha_hi <- alpha_bounds[2]
rho_lo   <- reproduction_number_bounds[1];  rho_hi   <- reproduction_number_bounds[2]
s0_lo    <- s0_bounds[1];  s0_hi    <- s0_bounds[2]

# ============================================================================
# Helper Functions
# ============================================================================

#' Convert normalized design matrix D to physical parameters
rescale_design <- function(D, Y = NULL) {
  alpha <- alpha_bounds[1] + D[, 1] * diff(alpha_bounds)
  rho   <- reproduction_number_bounds[1] + D[, 2] * diff(reproduction_number_bounds)
  s0    <- s0_bounds[1] + D[, 3] * diff(s0_bounds)
  beta  <- alpha / rho

  result <- data.frame(alpha = alpha, beta = beta, s0 = s0, rho = rho)

  if (!is.null(Y)) {
    result$piv <- Y[, 1]
    result$pit <- Y[, 2]
  }

  return(result)
}

#' Calculate grid coverage with exclusions
#' @param u Normalized x coordinates (0-1)
#' @param v Normalized y coordinates (0-1)
#' @param k Grid resolution
#' @param total_cells Total possible cells (for calculating percentage)
grid_coverage <- function(u, v, k, total_cells = NULL) {
  # Clamp to [0,1]
  u <- pmin(1, pmax(0, u))
  v <- pmin(1, pmax(0, v))

  # Convert to grid indices
  ix <- pmin(k, pmax(1, floor(u * k) + 1))
  iy <- pmin(k, pmax(1, floor(v * k) + 1))

  # Calculate unique cells
  identifiers <- ix + k * (iy - 1)
  n_unique <- length(unique(identifiers))

  if (is.null(total_cells)) {
    return(n_unique / k^2)
  } else {
    return(n_unique / total_cells)
  }
}

#' Create input space plot
make_input_plot <- function(data, title) {
  n <- nrow(data)

  ggplot(data, aes(x = alpha, y = log(rho))) +
    geom_point(size = 1, alpha = 0.9) +
    labs(
      title = title,
      subtitle = sprintf("n = %d", n),
      x = TeX(r'($\alpha$)'),
      y = TeX(r'(log($\rho = \alpha/\beta$))')
    ) +
    theme_bw() +
    theme(
      plot.subtitle = element_text(size = 20)
    )
}

#' Create output space plot with coverage annotation
make_output_plot <- function(data, title, coverage_pct) {
  label_txt <- sprintf("Grid (of possible) \n coverage = %.1f%%", coverage_pct)

  ggplot(data, aes(piv, pit)) +
    geom_point(size = 1, alpha = 0.9) +
    labs(title = title, x = "PIV", y = "PIT") +
    theme_bw() +
    xlim(0, piv_bounds[2]) +
    ylim(0, pit_bounds[2]) +
    annotate(
      "text",
      x = piv_bounds[2] * 0.98,
      y = pit_bounds[2] * 0.98,
      hjust = 1, vjust = 1,
      label = label_txt,
      size = 12,
      fontface = "bold"
    )
}

# ============================================================================
# Loop through sample sizes and create plots
# ============================================================================

NSAMPLES <- seq(1000, 10000, 1000)
k <- 50  # Grid resolution for coverage calculation

for (n in NSAMPLES) {
  cat(sprintf("Creating plots for n=%d...\n", n))

  # Method 1: SIR Mappings
  sirmap_file <- here::here('data', sprintf('inputs_and_outputs_sirMappings_n%d.RData', n))
  if (file.exists(sirmap_file)) {
    load(sirmap_file)
    data1 <- results %>%
      as_tibble() %>%
      filter(alpha > 0, beta > 0) %>%
      mutate(rho = alpha / beta)
  } else {
    cat(sprintf("  Warning: SIR_Maps file not found for n=%d\n", n))
    data1 <- data.frame(alpha = numeric(0), beta = numeric(0), rho = numeric(0),
                        piv = numeric(0), pit = numeric(0))
  }

  # Method 2: Basic LHS
  isfd_file <- here::here('data', sprintf('inputs_and_outputs_sirBasicLHS_n%d.RData', n))
  if (file.exists(isfd_file)) {
    load(isfd_file)
    data2 <- results %>%
      as_tibble() %>%
      mutate(rho = alpha / beta)
  } else {
    cat(sprintf("  Warning: Basic_LHS file not found for n=%d\n", n))
    data2 <- data.frame(alpha = numeric(0), beta = numeric(0), rho = numeric(0),
                        piv = numeric(0), pit = numeric(0))
  }

  # Method 3: OSFD Basic
  osfd_file <- Sys.glob(here::here('data', sprintf('fullOutput_OSFD_basic_budget%d_n*_rep1.RData', n)))
  if (length(osfd_file) > 0) {
    load(osfd_file[1])
    D3 <- res$D_success
    Y3 <- res$Y_success
    data3 <- rescale_design(D3, Y3)
  } else {
    cat(sprintf("  Warning: OSFD_Basic file not found for n=%d\n", n))
    data3 <- data.frame(alpha = numeric(0), beta = numeric(0), rho = numeric(0),
                        piv = numeric(0), pit = numeric(0))
  }

  # Method 4: OSFD WangOrig
  wangOrig_file <- Sys.glob(here::here('data', sprintf('fullOutput_OSFD_WangOrig_budget%d_n*.RData', n)))
  if (length(wangOrig_file) > 0) {
    load(wangOrig_file[1])
    D4 <- res$D_success
    Y4 <- res$Y_success
    data4 <- rescale_design(D4, Y4)
  } else {
    cat(sprintf("  Warning: OSFD_WangOrig file not found for n=%d\n", n))
    data4 <- data.frame(alpha = numeric(0), beta = numeric(0), rho = numeric(0),
                        piv = numeric(0), pit = numeric(0))
  }

  # Calculate total possible grid cells from all methods combined
  all_data <- bind_rows(
    mutate(data1, source = "SIR_Maps"),
    mutate(data2, source = "Basic_LHS"),
    mutate(data3, source = "OSFD_Basic"),
    mutate(data4, source = "OSFD_WangOrig")
  )

  # Normalize output space coordinates (PIV, PIT) to [0,1]
  u_all <- (all_data$piv - piv_bounds[1]) / diff(piv_bounds)
  v_all <- (all_data$pit - pit_bounds[1]) / diff(pit_bounds)

  # Count unique grid cells occupied by ALL methods combined
  total_possible_cells <- grid_coverage(u_all, v_all, k) * k^2

  # Helper to calculate coverage for a dataset
  calc_coverage <- function(data) {
    if (nrow(data) == 0) return(0)
    u <- (data$piv - piv_bounds[1]) / diff(piv_bounds)
    v <- (data$pit - pit_bounds[1]) / diff(pit_bounds)
    grid_coverage(u, v, k, total_possible_cells) * 100
  }

  # Calculate coverage for each method
  cov1 <- calc_coverage(data1)
  cov2 <- calc_coverage(data2)
  cov3 <- calc_coverage(data3)
  cov4 <- calc_coverage(data4)

  # Create plots for all four methods
  p1 <- make_input_plot(data1, "Input Space - SIR Maps")
  p2 <- make_output_plot(data1, "Output Space - SIR Maps", cov1)

  p3 <- make_input_plot(data2, "Input Space - Basic LHS")
  p4 <- make_output_plot(data2, "Output Space - Basic LHS", cov2)

  p5 <- make_input_plot(data3, "Input Space - OSFD Basic")
  p6 <- make_output_plot(data3, "Output Space - OSFD Basic", cov3)

  p7 <- make_input_plot(data4, "Input Space - OSFD WangOrig")
  p8 <- make_output_plot(data4, "Output Space - OSFD WangOrig", cov4)

  # Combine into 2x4 grid
  grid <- (p1 | p3 | p5 | p7) / (p2 | p4 | p6 | p8)

  grid_big <- grid & theme(
    axis.title = element_text(size = 24),
    axis.text  = element_text(size = 20),
    plot.title = element_text(size = 28),
    plot.subtitle = element_text(size = 24),
    legend.title = element_text(size = 22),
    legend.text  = element_text(size = 20)
  )

  # Save plot
  scale <- 0.7
  output_file <- here::here("viz", sprintf("sir_four_methods_n%d.png", n))
  ggsave(
    filename = output_file,
    plot = grid_big,
    width = 50 * scale,
    height = 20 * scale,
    dpi = 300
  )

  cat(sprintf("  Saved: %s\n", output_file))
  cat(sprintf("  Coverage: SIR Maps=%.1f%%, Basic LHS=%.1f%%, OSFD Basic=%.1f%%, OSFD WangOrig=%.1f%%\n",
              cov1, cov2, cov3, cov4))
}

cat("\nAll plots complete!\n")
