# Script to plot the STOCHASTIC SIR curve space-filling experiment.
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
library(latex2exp)
setwd(here::here())
source(here::here('R', 'sir_experiment_setup.R'))

# Grid size parameter: divides output space into k×k grid cells for coverage calculation
k = 50

# ==============================================================================
# Load all experimental results
# ==============================================================================
# All datasets contain parameters in our internal convention:
#   - alpha = transmission rate (will be labeled as β in plots)
#   - beta = recovery rate (will be labeled as γ in plots)
#   - rho = R₀ = alpha/beta

# Population size N for converting PIV counts to proportions
N <- 10000

# Load basic LHS with replicates
# LHS = Latin Hypercube Sampling in input space (alpha, rho, s0)
load(file = here::here('data', 'inputs_and_outputs_sirStochasticLHS_wReplicates.RData'))
lhs_w_replicates <- results
# Convert PIV from absolute count to proportion (stochastic model uses counts, we plot proportions)
lhs_w_replicates$piv <- lhs_w_replicates$piv / N

# Load basic LHS without replicates
load(file = here::here('data', 'inputs_and_outputs_sirStochasticLHS_noReplicates.RData'))
lhs_no_replicates <- results
# Convert PIV from absolute count to proportion
lhs_no_replicates$piv <- lhs_no_replicates$piv / N

# Load Failure-aware OSFD with replicates (D and Y saved together)
# Failure-aware OSFD = Output-Space Filling Design, adaptively samples to achieve uniform coverage in output (PIV, PIT) space
load(file = here::here("data", "inputs_and_outputs_StochasticSirOutputFilling_wReplicates.RData"))
D_w_replicates <- D  # Normalized input design matrix [0,1]
Y_w_replicates <- Y  # Output observations (PIV in counts, PIT in time units)

# Transform Failure-aware OSFD with replicates from normalized [0,1] to physical parameters
alpha_w_rep <- alpha_bounds[1] + D_w_replicates[, 1] * diff(alpha_bounds)  # Transmission rate (our alpha)
s0_w_rep    <- s0_bounds[1]    + D_w_replicates[, 3] * diff(s0_bounds)     # Initial susceptible fraction
rho_w_rep   <- reproduction_number_bounds[1] + D_w_replicates[, 2] * diff(reproduction_number_bounds)  # R₀
beta_w_rep  <- alpha_w_rep / rho_w_rep  # Recovery rate (our beta) = alpha/rho

osfd_w_replicates <- data.frame(
  alpha = alpha_w_rep,  # Transmission rate (will be labeled as β in plots)
  beta  = beta_w_rep,   # Recovery rate (will be labeled as γ in plots)
  piv   = Y_w_replicates[, 1] / N,  # Convert PIV from count to proportion
  pit   = Y_w_replicates[, 2]       # Peak incidence time (already in time units)
)

# Load Failure-aware OSFD without replicates (D and Y saved together)
load(file = here::here("data", "inputs_and_outputs_StochasticSirOutputFilling_noReplicates.RData"))
D_no_replicates <- D  # Normalized input design matrix [0,1]
Y_no_replicates <- Y  # Output observations (PIV in counts, PIT in time units)

# Transform Failure-aware OSFD without replicates from normalized [0,1] to physical parameters
alpha_no_rep <- alpha_bounds[1] + D_no_replicates[, 1] * diff(alpha_bounds)  # Transmission rate (our alpha)
s0_no_rep    <- s0_bounds[1]    + D_no_replicates[, 3] * diff(s0_bounds)     # Initial susceptible fraction
rho_no_rep   <- reproduction_number_bounds[1] + D_no_replicates[, 2] * diff(reproduction_number_bounds)  # R₀
beta_no_rep  <- alpha_no_rep / rho_no_rep  # Recovery rate (our beta) = alpha/rho

osfd_no_replicates <- data.frame(
  alpha = alpha_no_rep,  # Transmission rate (will be labeled as β in plots)
  beta  = beta_no_rep,   # Recovery rate (will be labeled as γ in plots)
  piv   = Y_no_replicates[, 1] / N,  # Convert PIV from count to proportion
  pit   = Y_no_replicates[, 2]       # Peak incidence time (already in time units)
)

# Combine all results for overall grid coverage calculation
all_results <- rbind(lhs_w_replicates, lhs_no_replicates, osfd_w_replicates, osfd_no_replicates)

# For stochastic SIR, PIT (time to peak) has a different range than deterministic SIR
# Stochastic outbreaks peak faster (median ~3 days) but some take longer (up to 400 days)
# Override pit_bounds to match the stochastic data distribution
stochastic_pit_bounds <- c(0, 100)  # Captures 99.5% of data
cat("Using stochastic PIT bounds:", stochastic_pit_bounds, "(captures 99.5% of data)\n")
cat("PIV bounds (proportions):", piv_bounds, "\n")

# ==============================================================================
# Calculate baseline grid coverage (from all experiments combined)
# This identifies which grid cells actually contain observations, used to
# calculate the percentage of *feasible* grid cells covered by each method
# ==============================================================================

# Normalize output space to [0,1]
u <- (all_results$piv - piv_bounds[1]) / diff(piv_bounds)
v <- (all_results$pit - stochastic_pit_bounds[1]) / diff(stochastic_pit_bounds)
u <- pmin(1, pmax(0, u))
v <- pmin(1, pmax(0, v))

# Map to grid indices (1 to k)
ix <- pmin(k, pmax(1, floor(u * k) + 1))
iy <- pmin(k, pmax(1, floor(v * k) + 1))

# Create unique identifier for each grid cell
exclusion_identifiers <- (ix + k * (iy - 1))
num_grid_w_obs <- length(unique(exclusion_identifiers))

# ==============================================================================
# Helper function to calculate grid coverage
# ==============================================================================
# This function calculates what percentage of "achievable" grid cells are
# covered by a given dataset.
#
# The "grid of possible" is defined as all grid cells occupied by ANY method
# (num_grid_w_obs). We measure what fraction of these achievable cells each
# individual method covers.
#
# Why not use k^2? Because some grid cells may be infeasible (no parameter
# combinations lead to outputs in those cells). Using the union of all methods
# as the baseline gives a more fair comparison.

grid_coverage_with_exclusions <- function(u, v, k) {
  # Clamp to [0, 1] to handle boundary cases
  u <- pmin(1, pmax(0, u))
  v <- pmin(1, pmax(0, v))

  # Convert to grid indices (1 to k)
  ix <- pmin(k, pmax(1, floor(u * k) + 1))
  iy <- pmin(k, pmax(1, floor(v * k) + 1))

  # Create unique identifier for each grid cell
  # Formula: converts 2D grid position (ix, iy) to 1D index
  identifiers <- ix + k * (iy - 1)

  # Return: (# unique cells this method covers) / (# achievable cells total)
  length(unique(identifiers)) / num_grid_w_obs
}

# ==============================================================================
# Plot 1 & 2: Basic LHS with replicates
# ==============================================================================
# IMPORTANT: Plot 1 shows INPUT SPACE
#   - x-axis plots "alpha" (transmission rate) but is LABELED as β (paper convention)
#   - y-axis plots log(rho) = log(alpha/beta) which is log(R₀) = log(β/γ) in paper notation

p1 <- ggplot(lhs_w_replicates, aes(x = alpha, y = log(alpha / beta))) +  # Plotting alpha on x-axis
  geom_point(size = 1, alpha = 0.9) +
  labs(
    title = "Input Space - Basic LHS \nwith Replicates",
    subtitle = sprintf("N = %s", format(nrow(lhs_w_replicates), big.mark = ",")),
    # KEY: x-axis plots "alpha" but is LABELED as β to match paper convention
    x = TeX(r'($\beta$)'),  # Label as β (paper convention) even though plotting alpha
    y = TeX(r'(log($R_0 = \beta/\gamma$))')  # R₀ in paper notation
  ) +
  theme_bw()

# Calculate grid coverage for this method
# Normalize outputs to [0,1] then compute what fraction of achievable cells are covered
u <- (lhs_w_replicates$piv - piv_bounds[1]) / diff(piv_bounds)
v <- (lhs_w_replicates$pit - stochastic_pit_bounds[1]) / diff(stochastic_pit_bounds)
cov <- grid_coverage_with_exclusions(u, v, k = k)
label_txt <- sprintf("Grid (of possible) \n coverage = %.1f%%", 100 * cov)

# Plot 2 shows OUTPUT SPACE (PIV vs PIT)
p2 <- ggplot(lhs_w_replicates, aes(piv, pit)) +
  geom_point(size = 1, alpha = 0.9) +
  labs(title = "Output Space - Basic LHS \nwith Replicates", x = "PIV", y = "PIT") +
  theme_bw() + xlim(piv_bounds[1], piv_bounds[2]) + ylim(stochastic_pit_bounds[1], stochastic_pit_bounds[2]) +
  annotate(
    "text",
    x = piv_bounds[1] + diff(piv_bounds) * 0.02,
    y = stochastic_pit_bounds[2] * 0.98,
    hjust = 0, vjust = 1,
    label = label_txt,
    size = 12,
    fontface = "bold"
  )

# ==============================================================================
# Plot 3 & 4: Basic LHS without replicates
# ==============================================================================
# IMPORTANT: Plot 3 shows INPUT SPACE
#   - x-axis plots "alpha" (transmission rate) but is LABELED as β (paper convention)
#   - y-axis plots log(rho) = log(alpha/beta) which is log(R₀) = log(β/γ) in paper notation

p3 <- ggplot(lhs_no_replicates, aes(x = alpha, y = log(alpha / beta))) +  # Plotting alpha on x-axis
  geom_point(size = 1, alpha = 0.9) +
  labs(
    title = "Input Space - Basic LHS, \nno Replicates",
    subtitle = sprintf("N = %s", format(nrow(lhs_no_replicates), big.mark = ",")),
    # KEY: x-axis plots "alpha" but is LABELED as β to match paper convention
    x = TeX(r'($\beta$)'),  # Label as β (paper convention) even though plotting alpha
    y = TeX(r'(log($R_0 = \beta/\gamma$))')  # R₀ in paper notation
  ) +
  theme_bw()

# Calculate grid coverage for this method
u <- (lhs_no_replicates$piv - piv_bounds[1]) / diff(piv_bounds)
v <- (lhs_no_replicates$pit - stochastic_pit_bounds[1]) / diff(stochastic_pit_bounds)
cov <- grid_coverage_with_exclusions(u, v, k = k)
label_txt <- sprintf("Grid (of possible) \n coverage = %.1f%%", 100 * cov)

# Plot 4 shows OUTPUT SPACE (PIV vs PIT)
p4 <- ggplot(lhs_no_replicates, aes(piv, pit)) +
  geom_point(size = 1, alpha = 0.9) +
  labs(title = "Output Space - Basic LHS, \nno Replicates", x = "PIV", y = "PIT") +
  theme_bw() + xlim(piv_bounds[1], piv_bounds[2]) + ylim(stochastic_pit_bounds[1], stochastic_pit_bounds[2]) +
  annotate(
    "text",
    x = piv_bounds[1] + diff(piv_bounds) * 0.02,
    y = stochastic_pit_bounds[2] * 0.98,
    hjust = 0, vjust = 1,
    label = label_txt,
    size = 12,
    fontface = "bold"
  )

# ==============================================================================
# Plot 5 & 6: OSFD with replicates
# ==============================================================================
# IMPORTANT: Plot 5 shows INPUT SPACE
#   - x-axis plots "alpha" (transmission rate) but is LABELED as β (paper convention)
#   - y-axis plots log(rho) = log(alpha/beta) which is log(R₀) = log(β/γ) in paper notation

p5 <- ggplot(osfd_w_replicates, aes(x = alpha, y = log(alpha / beta))) +  # Plotting alpha on x-axis
  geom_point(size = 1, alpha = 0.9) +
  labs(
    title = "Input Space - Failure-aware OSFD \nwith Replicates",
    subtitle = sprintf("N = %s", format(nrow(osfd_w_replicates), big.mark = ",")),
    # KEY: x-axis plots "alpha" but is LABELED as β to match paper convention
    x = TeX(r'($\beta$)'),  # Label as β (paper convention) even though plotting alpha
    y = TeX(r'(log($R_0 = \beta/\gamma$))')  # R₀ in paper notation
  ) +
  theme_bw()

# Calculate grid coverage for this method
u <- (osfd_w_replicates$piv - piv_bounds[1]) / diff(piv_bounds)
v <- (osfd_w_replicates$pit - stochastic_pit_bounds[1]) / diff(stochastic_pit_bounds)
cov <- grid_coverage_with_exclusions(u, v, k = k)
label_txt <- sprintf("Grid (of possible) \n coverage = %.1f%%", 100 * cov)

# Plot 6 shows OUTPUT SPACE (PIV vs PIT)
p6 <- ggplot(osfd_w_replicates, aes(piv, pit)) +
  geom_point(size = 1, alpha = 0.9) +
  labs(title = "Output Space - Failure-aware OSFD \nwith Replicates", x = "PIV", y = "PIT") +
  theme_bw() + xlim(piv_bounds[1], piv_bounds[2]) + ylim(stochastic_pit_bounds[1], stochastic_pit_bounds[2]) +
  annotate(
    "text",
    x = piv_bounds[1] + diff(piv_bounds) * 0.02,
    y = stochastic_pit_bounds[2] * 0.98,
    hjust = 0, vjust = 1,
    label = label_txt,
    size = 12,
    fontface = "bold"
  )

# ==============================================================================
# Plot 7 & 8: Failure-aware OSFD without replicates
# ==============================================================================
# IMPORTANT: Plot 7 shows INPUT SPACE
#   - x-axis plots "alpha" (transmission rate) but is LABELED as β (paper convention)
#   - y-axis plots log(rho) = log(alpha/beta) which is log(R₀) = log(β/γ) in paper notation

p7 <- ggplot(osfd_no_replicates, aes(x = alpha, y = log(alpha / beta))) +  # Plotting alpha on x-axis
  geom_point(size = 1, alpha = 0.9) +
  labs(
    title = "Input Space - Failure-aware OSFD, \nno Replicates",
    subtitle = sprintf("N = %s", format(nrow(osfd_no_replicates), big.mark = ",")),
    # KEY: x-axis plots "alpha" but is LABELED as β to match paper convention
    x = TeX(r'($\beta$)'),  # Label as β (paper convention) even though plotting alpha
    y = TeX(r'(log($R_0 = \beta/\gamma$))')  # R₀ in paper notation
  ) +
  theme_bw()

# Calculate grid coverage for this method
u <- (osfd_no_replicates$piv - piv_bounds[1]) / diff(piv_bounds)
v <- (osfd_no_replicates$pit - stochastic_pit_bounds[1]) / diff(stochastic_pit_bounds)
cov <- grid_coverage_with_exclusions(u, v, k = k)
label_txt <- sprintf("Grid (of possible) \n coverage = %.1f%%", 100 * cov)

# Plot 8 shows OUTPUT SPACE (PIV vs PIT)
p8 <- ggplot(osfd_no_replicates, aes(piv, pit)) +
  geom_point(size = 1, alpha = 0.9) +
  labs(title = "Output Space - Failure-aware OSFD, \nno Replicates", x = "PIV", y = "PIT") +
  theme_bw() + xlim(piv_bounds[1], piv_bounds[2]) + ylim(stochastic_pit_bounds[1], stochastic_pit_bounds[2]) +
  annotate(
    "text",
    x = piv_bounds[1] + diff(piv_bounds) * 0.02,
    y = stochastic_pit_bounds[2] * 0.98,
    hjust = 0, vjust = 1,
    label = label_txt,
    size = 12,
    fontface = "bold"
  )

# ==============================================================================
# Combine plots into final grid and save
# ==============================================================================
# Create a 2×4 grid layout:
#   Top row: Input space plots for all four methods
#   Bottom row: Output space plots for all four methods
#
# REMINDER: Input space x-axis shows alpha (transmission rate) but is
#           LABELED as β to match epidemiology paper convention

# Arrange plots: top row = input spaces, bottom row = output spaces
# Columns: Basic LHS w/ rep, Failure-aware OSFD w/ rep, Basic LHS no rep, Failure-aware OSFD no rep
grid <- (p1 | p5 | p3 | p7) / (p2 | p6 | p4 | p8)

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
ggsave(
  filename = here::here("viz", "stochastic_sir_experiment.png"),
  plot = grid_big,
  width = 40 * scale,
  height = 20 * scale,
  dpi = 300
)

cat("Plot saved to viz/stochastic_sir_experiment.png\n")
cat("\nNOTE: Input space x-axis plots 'alpha' (transmission rate) but is labeled as 'β' (paper convention)\n")
