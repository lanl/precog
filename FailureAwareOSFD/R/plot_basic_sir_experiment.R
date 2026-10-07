# Script to plot the SIR curve space-filling experiment.
# Author: AC Murph
#
# IMPORTANT NAMING CONVENTION:
# =============================
# This codebase uses a naming convention that differs from some epidemiology papers.
#
# In this code:
#   - "alpha" = transmission rate (often called β in epidemiology literature)
#   - "beta" = recovery rate (often called γ in epidemiology literature)
#   - "rho" = basic reproduction number = alpha/beta (often called R₀)
#
# In some epidemiology papers (and in the plots we create):
#   - β (beta) = transmission rate = our "alpha"
#   - γ (gamma) = recovery rate = our "beta"
#
# Therefore, when PLOTTING, we show:
#   - x-axis: our "alpha" parameter, but LABELED as "β" to match paper convention
#   - y-axis: log(rho) = log(α/β) in our notation, which is log(β/γ) = log(R₀) in papers
#
# This allows the code to remain consistent internally while presenting results
# in the conventional notation expected by readers familiar with the epidemiology literature.
# =============================

library(ggplot2)
library(patchwork)
library(tidyr)
library(dplyr)
library(latex2exp)
setwd(here::here())
source(here::here("R", "sir_experiment_setup.R"))

# Define parameter bounds
# These bounds are in our internal naming convention where:
#   - alpha = transmission rate (β in papers)
#   - beta = recovery rate (γ in papers) [computed as alpha/rho]
#   - rho = reproduction number R₀ (alpha/beta in our code)
k = 50
# alpha_bounds = c(0.001, 100)  # Transmission rate bounds (β in epidemiology notation)
# reproduction_number_bounds = c(1.001, 100)  # R₀ bounds
# pit_bounds = c(3, 50)  # Peak incidence time bounds
# piv_bounds = c(0.005, 1)  # Peak incidence value bounds
# s0_bounds = c(0.95, 0.999)  # Initial susceptible fraction bounds

# ============================================================================
# Helper Functions
# ============================================================================

#' Convert normalized design matrix D to physical parameters
#' @param D Design matrix with columns [alpha_norm, rho_norm, s0_norm]
#'          where values are in [0,1] (normalized)
#' @param Y Output matrix with columns [piv, pit] (optional)
#' @return Data frame with columns: alpha, beta, s0, rho (and piv, pit if Y provided)
#'
#' NAMING: In this code, alpha = transmission rate, beta = recovery rate
#'         This differs from epidemiology papers where β = transmission, γ = recovery
rescale_design <- function(D, Y = NULL) {
  # Convert normalized [0,1] inputs to physical parameter ranges
  alpha <- alpha_bounds[1] + D[, 1] * diff(alpha_bounds)  # Transmission rate (our alpha)
  rho   <- reproduction_number_bounds[1] + D[, 2] * diff(reproduction_number_bounds)  # R₀
  s0    <- s0_bounds[1] + D[, 3] * diff(s0_bounds)  # Initial susceptible fraction

  # Calculate recovery rate (our beta) from alpha and rho
  # rho = alpha/beta, so beta = alpha/rho
  beta  <- alpha / rho  # Recovery rate (our beta)

  result <- data.frame(alpha = alpha, beta = beta, s0 = s0, rho = rho)

  # Add output variables if provided
  if (!is.null(Y)) {
    result$piv <- Y[, 1]  # Peak incidence value
    result$pit <- Y[, 2]  # Peak incidence time
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
#'
#' IMPORTANT: This function plots the "alpha" parameter on the x-axis, but LABELS it as β
#'            to match epidemiology paper convention where β = transmission rate.
#'            In our code, alpha = transmission rate, but papers call this β.
#'
#' @param data Data frame with columns: alpha, beta, rho
#' @param title Plot title
make_input_plot <- function(data, title) {
  n <- nrow(data)

  ggplot(data, aes(x = alpha, y = log(rho))) +  # Plotting alpha (transmission rate)
    geom_point(size = 1, alpha = 0.9) +
    labs(
      title = title,
      subtitle = sprintf("n = %s", format(n, big.mark = ",", scientific = FALSE)),
      # KEY: x-axis plots "alpha" but is LABELED as β to match paper convention
      x = TeX(r'($\beta$)'),  # Label as β (paper convention) even though plotting alpha
      y = TeX(r'(log($R_0 = \beta/\gamma$))')  # R₀ in paper notation
    ) +
    theme_bw()
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
# Load All Data
# ============================================================================
# Each dataset contains parameters in our internal convention:
#   - alpha = transmission rate (labeled as β in plots)
#   - beta = recovery rate (labeled as γ in plots)
#   - rho = R₀ = alpha/beta

# ============================================================================
# SELECT BUDGET HERE
# ============================================================================
BUDGET <- 5000  # Change this value to select different budget

cat(sprintf("\n=== Loading data for budget n = %d ===\n", BUDGET))

# Dataset 1: SIR Mappings (inverse problem: PIV/PIT → alpha/beta)
# Load first replicate for this budget
sirmap_file <- list.files(
  here::here('data'),
  pattern = sprintf("^inputs_and_outputs_sirMappings_n%d_rep1\\.RData$", BUDGET),
  full.names = TRUE
)

if (length(sirmap_file) == 0) {
  stop(sprintf("No SIR Mappings file found for budget %d (rep 1)", BUDGET))
}

load(sirmap_file[1])
data1 <- results %>%
  as_tibble() %>%
  filter(alpha > 0, beta > 0) %>%  # Remove invalid parameter combinations
  mutate(rho = alpha / beta)  # Calculate R₀

cat(sprintf("  Loaded SIR Mappings: %d points\n", nrow(data1)))

# Dataset 2: Basic LHS (Latin Hypercube Sampling in input space)
# Directly samples alpha, rho, s0 and computes beta = alpha/rho
lhs_file <- list.files(
  here::here('data'),
  pattern = sprintf("^inputs_and_outputs_sirBasicLHS_n%d_rep1\\.RData$", BUDGET),
  full.names = TRUE
)

if (length(lhs_file) == 0) {
  stop(sprintf("No Basic LHS file found for budget %d (rep 1)", BUDGET))
}

load(lhs_file[1])
data2 <- results %>%
  as_tibble() %>%
  mutate(rho = alpha / beta)  # Calculate R₀

cat(sprintf("  Loaded Basic LHS: %d points\n", nrow(data2)))

# Dataset 3: Failure-aware OSFD (Output-Space Filling Design)
# Adaptively samples to achieve uniform coverage in OUTPUT space (PIV, PIT)
osfd_file <- list.files(
  here::here('data'),
  pattern = sprintf("^fullOutput_OSFD_basic_budget%d_n\\d+_rep1\\.RData$", BUDGET),
  full.names = TRUE
)

if (length(osfd_file) == 0) {
  stop(sprintf("No OSFD file found for budget %d (rep 1)", BUDGET))
}

load(osfd_file[1])
# OSFD files contain res$Y_success (outputs) and res$D_success (normalized inputs)
Y <- res$Y_success
D <- res$D_success
data3 <- rescale_design(D, Y)  # Convert normalized inputs to physical parameters

cat(sprintf("  Loaded Failure-aware OSFD: %d points\n", nrow(data3)))


# ============================================================================
# Calculate Total Possible Grid Cells (for coverage calculation)
# ============================================================================
# We estimate the "theoretically achievable" output space by pooling all methods.
# This gives a more realistic coverage metric than using the full k×k grid,
# since some (PIV, PIT) combinations may be impossible given the parameter constraints.

# Combine all datasets to find which cells are theoretically achievable
all_data <- bind_rows(
  mutate(data1, source = "SIR_Maps"),
  mutate(data2, source = "Basic_LHS"),
  mutate(data3, source = "Failure-aware OSFD")
)

# Normalize output space coordinates (PIV, PIT) to [0,1]
u_all <- (all_data$piv - piv_bounds[1]) / diff(piv_bounds)
v_all <- (all_data$pit - pit_bounds[1]) / diff(pit_bounds)

# Count unique grid cells occupied by ALL methods combined
total_possible_cells <- grid_coverage(u_all, v_all, k) * k^2

# ============================================================================
# Create Plots for Each Method
# ============================================================================
# For each sampling method, we create:
#   1. INPUT SPACE plot: Shows distribution in (alpha, log(rho)) space
#      Note: x-axis plots alpha but is LABELED as β (paper convention)
#   2. OUTPUT SPACE plot: Shows distribution in (PIV, PIT) space
#      with grid coverage percentage relative to achievable cells

# Helper to calculate coverage for a dataset
# Coverage = (unique cells occupied by this method) / (total achievable cells) × 100%
calc_coverage <- function(data) {
  u <- (data$piv - piv_bounds[1]) / diff(piv_bounds)  # Normalize PIV to [0,1]
  v <- (data$pit - pit_bounds[1]) / diff(pit_bounds)  # Normalize PIT to [0,1]
  grid_coverage(u, v, k, total_possible_cells) * 100  # Return percentage
}

# --- Method 1: SIR Mappings (inverse problem approach) ---
cov1 <- calc_coverage(data1)
p1 <- make_input_plot(data1, "Input Space - SIR Mappings")  # x-axis = alpha, labeled as β
p2 <- make_output_plot(data1, "Output Space - SIR Mappings", cov1)

# --- Method 2: Basic LHS (standard input-space sampling) ---
cov2 <- calc_coverage(data2)
p3 <- make_input_plot(data2, "Input Space - Basic LHS")  # x-axis = alpha, labeled as β
p4 <- make_output_plot(data2, "Output Space - Basic LHS", cov2)

# --- Method 3: Failure-aware OSFD (output-space filling design) ---
cov3 <- calc_coverage(data3)
p5 <- make_input_plot(data3, "Input Space - Failure-aware OSFD")  # x-axis = alpha, labeled as β
p6 <- make_output_plot(data3, "Output Space - Failure-aware OSFD", cov3)


# ============================================================================
# Combine and Save
# ============================================================================
# Create a 2×3 grid layout:
#   Top row: Input space plots for all three methods
#   Bottom row: Output space plots for all three methods
#
# REMINDER: Input space x-axis shows alpha (transmission rate) but is
#           LABELED as β to match epidemiology paper convention

grid <- (p1 | p3 | p5) / (p2 | p4 | p6)

# Apply consistent theme with large text for publication
grid_big <- grid & theme(
  axis.title = element_text(size = 24),
  axis.text  = element_text(size = 20),
  plot.title = element_text(size = 28),
  plot.subtitle = element_text(size = 24),
  legend.title = element_text(size = 22),
  legend.text  = element_text(size = 20)
)

# Save with high resolution
scale <- 0.7
output_filename <- sprintf("sir_experiment_n%d.png", BUDGET)
ggsave(
  filename = here::here("viz", output_filename),
  plot = grid_big,
  width = 40 * scale,
  height = 20 * scale,
  dpi = 300
)

cat(sprintf("\nPlot saved to viz/%s\n", output_filename))
cat(sprintf("Coverage: SIR Maps=%.1f%%, Basic LHS=%.1f%%, OSFD=%.1f%%\n",
            cov1, cov2, cov3))
cat("\nNOTE: Input space x-axis plots 'alpha' (transmission rate) but is labeled as 'β' (paper convention)\n")
