# Calculate grid coverage and distance metrics for all three methods
# This script loads the output files from each method and calculates:
#   - Grid coverage percentage
#   - Fill distance (max nearest neighbor distance)
#   - Mean distance (average nearest neighbor distance)
# Usage: Rscript calculate_grid_coverage.R <nsuccesses> <replicate_number>
# Author: AC Murph
# Date: Sep 2026

library(dplyr)
library(here)
setwd(here::here())
source(here::here("R", "sir_experiment_setup.R"))

# Load helper functions for metrics calculation
source(here::here("R", "helpers", "wang_helper_functions.R"))

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  stop("Usage: Rscript calculate_grid_coverage.R <nsuccesses> <replicate_number>")
}

nsuccesses <- as.numeric(args[1])
rep_num <- as.numeric(args[2])
rep_suffix <- sprintf("_rep%d", rep_num)

cat(sprintf("Calculating metrics for nsuccesses=%d, replicate=%d\n", nsuccesses, rep_num))

# ============================================================================
# Grid Coverage Function (from plot_basic_sir_experiment.R)
# Note: Other metric functions are loaded from wang_helper_functions.R
# ============================================================================

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

# ============================================================================
# Load Data for All Three Methods
# ============================================================================

k <- 100  # Grid resolution
n <- nsuccesses

# Method 1: SIR Mappings (inverse problem)
cat("Loading Method 1: SIR Mappings...\n")
load(here::here('data', paste0('inputs_and_outputs_sirMappings_n', format(n, scientific = FALSE), rep_suffix, '.RData')))
data1 <- results %>%
  as_tibble() %>%
  filter(alpha > 0, beta > 0) %>%
  mutate(rho = alpha / beta)

# Method 2: Basic LHS
cat("Loading Method 2: Basic LHS...\n")
load(here::here('data', paste0('inputs_and_outputs_sirBasicLHS_n', format(n, scientific = FALSE), rep_suffix, '.RData')))
data2 <- results %>%
  as_tibble() %>%
  mutate(rho = alpha / beta)

# Method 3: OSFD Basic (budget mode)
cat("Loading Method 3: OSFD Basic...\n")
# Try budget filename first, fall back to old naming if not found
osfd_file <- here::here("data", sprintf("outputs_sirOutputFilling_budget%d_n*%s.RData", n, rep_suffix))
osfd_files <- Sys.glob(osfd_file)
if (length(osfd_files) > 0) {
  # Budget mode files exist
  load(osfd_files[1])  # Y output
  Y3 <- Y
  # Load corresponding input file
  input_file <- gsub("outputs_sirOutputFilling", "inputs_sirOutputFilling", osfd_files[1])
  load(input_file)
  D3 <- D
} else {
  # Fall back to old naming
  load(here::here("data", paste0("inputs_sirOutputFilling_n", format(n, scientific = FALSE), rep_suffix, ".RData")))
  D3 <- D
  load(here::here("data", paste0("outputs_sirOutputFilling_n", format(n, scientific = FALSE), rep_suffix, ".RData")))
  Y3 <- Y
}
data3 <- rescale_design(D3, Y3)

# ============================================================================
# Calculate Total Possible Grid Cells
# ============================================================================
cat("Calculating total possible cells...\n")

all_data <- bind_rows(
  mutate(data1, source = "SIR_Maps"),
  mutate(data2, source = "Basic_LHS"),
  mutate(data3, source = "OSFD_Basic")
)

# Normalize output space coordinates to [0,1]
u_all <- (all_data$piv - piv_bounds[1]) / diff(piv_bounds)
v_all <- (all_data$pit - pit_bounds[1]) / diff(pit_bounds)

# Count unique grid cells occupied by ALL methods combined
total_possible_cells <- grid_coverage(u_all, v_all, k) * k^2

cat(sprintf("Total possible cells: %.0f out of %d\n", total_possible_cells, k^2))

# ============================================================================
# Calculate Metrics for Each Method
# ============================================================================
cat("Calculating metrics for each method...\n")

# Create Y_union from all methods (this is the reference space)
Y_union <- rbind(
  as.matrix(data1[, c("piv", "pit")]),
  as.matrix(data2[, c("piv", "pit")]),
  as.matrix(data3[, c("piv", "pit")])
)

# Function to calculate all metrics for a method
calc_metrics <- function(data, method_name) {
  # Extract output matrix
  Y_method <- as.matrix(data[, c("piv", "pit")])

  # Scale both method and union to [0,1] based on union range
  Y_union_scaled <- scale_to_reference(Y_union, Y_union)
  Y_method_scaled <- scale_to_reference(Y_method, Y_union)

  # Compute nearest neighbor distances: for each union point, distance to nearest method point
  d <- nearest_distances(Y_union_scaled, Y_method_scaled)

  # Compute grid coverage
  u <- (data$piv - piv_bounds[1]) / diff(piv_bounds)
  v <- (data$pit - pit_bounds[1]) / diff(pit_bounds)
  coverage_pct <- grid_coverage(u, v, k, total_possible_cells) * 100

  # Return all metrics
  list(
    method = method_name,
    coverage_pct = coverage_pct,
    fill_distance = max(d),
    mean_distance = mean(d),
    median_distance = median(d),
    q90_distance = as.numeric(quantile(d, 0.90)),
    q95_distance = as.numeric(quantile(d, 0.95)),
    n_points = nrow(data)
  )
}

# Calculate metrics for each method
metrics_sir_maps <- calc_metrics(data1, "SIR_Maps")
metrics_basic_lhs <- calc_metrics(data2, "Basic_LHS")
metrics_osfd_basic <- calc_metrics(data3, "OSFD_Basic")

cat(sprintf("\nMethod: SIR Maps\n"))
cat(sprintf("  Coverage: %.2f%%, Fill Distance: %.4f, Mean Distance: %.4f\n",
            metrics_sir_maps$coverage_pct, metrics_sir_maps$fill_distance, metrics_sir_maps$mean_distance))

cat(sprintf("\nMethod: Basic LHS\n"))
cat(sprintf("  Coverage: %.2f%%, Fill Distance: %.4f, Mean Distance: %.4f\n",
            metrics_basic_lhs$coverage_pct, metrics_basic_lhs$fill_distance, metrics_basic_lhs$mean_distance))

cat(sprintf("\nMethod: OSFD Basic\n"))
cat(sprintf("  Coverage: %.2f%%, Fill Distance: %.4f, Mean Distance: %.4f\n",
            metrics_osfd_basic$coverage_pct, metrics_osfd_basic$fill_distance, metrics_osfd_basic$mean_distance))

# ============================================================================
# Save Results
# ============================================================================

coverage_results <- data.frame(
  nsuccesses = n,
  replicate = rep_num,
  method = c("SIR_Maps", "Basic_LHS", "OSFD_Basic"),
  coverage_pct = c(metrics_sir_maps$coverage_pct, metrics_basic_lhs$coverage_pct,
                   metrics_osfd_basic$coverage_pct),
  fill_distance = c(metrics_sir_maps$fill_distance, metrics_basic_lhs$fill_distance,
                    metrics_osfd_basic$fill_distance),
  mean_distance = c(metrics_sir_maps$mean_distance, metrics_basic_lhs$mean_distance,
                    metrics_osfd_basic$mean_distance),
  median_distance = c(metrics_sir_maps$median_distance, metrics_basic_lhs$median_distance,
                      metrics_osfd_basic$median_distance),
  q90_distance = c(metrics_sir_maps$q90_distance, metrics_basic_lhs$q90_distance,
                   metrics_osfd_basic$q90_distance),
  q95_distance = c(metrics_sir_maps$q95_distance, metrics_basic_lhs$q95_distance,
                   metrics_osfd_basic$q95_distance),
  n_points = c(metrics_sir_maps$n_points, metrics_basic_lhs$n_points,
               metrics_osfd_basic$n_points),
  total_possible_cells = total_possible_cells,
  grid_resolution = k
)

# Save to data directory
output_file <- here::here("data", sprintf("grid_coverage_n%d_rep%d.RData", n, rep_num))
save(coverage_results, file = output_file)

# Also save as CSV for easy reading
csv_file <- here::here("data", sprintf("grid_coverage_n%d_rep%d.csv", n, rep_num))
write.csv(coverage_results, file = csv_file, row.names = FALSE)

cat(sprintf("\nResults saved to:\n"))
cat(sprintf("  - %s\n", output_file))
cat(sprintf("  - %s\n", csv_file))
cat("\nMetrics calculation complete.\n")
cat("Metrics include: coverage_pct, fill_distance, mean_distance, median_distance, q90_distance, q95_distance\n")
