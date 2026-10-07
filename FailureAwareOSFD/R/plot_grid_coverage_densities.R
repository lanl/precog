# Script to plot grid coverage metrics as densities by method
# Author: AC Murph
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(here)
library(patchwork)
library(fields)
library(parallel)
library(doParallel)
library(doSNOW)
library(foreach)

setwd(here::here())
source(here::here('R', 'sir_experiment_setup.R'))

# ============================================================================
# Load all grid coverage CSV files
# ============================================================================

csv_files <- list.files(
  path = here::here("data"),
  pattern = "^grid_coverage_n5000_rep\\d+\\.csv$",
  full.names = TRUE
)

cat(sprintf("Found %d CSV files\n", length(csv_files)))

# Read and combine all files
all_data <- csv_files %>%
  map_dfr(read.csv) %>%
  as_tibble()

cat(sprintf("Total rows: %d\n", nrow(all_data)))
cat(sprintf("Unique methods: %s\n", paste(unique(all_data$method), collapse = ", ")))
cat(sprintf("Replicates per method: %d\n", n_distinct(all_data$replicate)))

# ============================================================================
# Reshape data for faceting
# ============================================================================

# Select metrics and rename for display
plot_data <- all_data %>%
  select(method, coverage_pct, fill_distance, mean_distance) %>%
  pivot_longer(
    cols = c(coverage_pct, fill_distance, mean_distance),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    # Rename methods for display
    method = case_when(
      method == "Basic_LHS" ~ "LHS",
      method == "OSFD_Basic" ~ "OSFD",
      method == "SIR_Maps" ~ "Inverse SIR Maps",
      TRUE ~ method
    ),
    metric_label = case_when(
      metric == "coverage_pct" ~ "Grid Coverage",
      metric == "fill_distance" ~ "Fill Distance",
      metric == "mean_distance" ~ "Mean Distance"
    ),
    metric_label = factor(
      metric_label,
      levels = c("Grid Coverage", "Fill Distance", "Mean Distance")
    )
  )

# ============================================================================
# Create density plot
# ============================================================================

# Use viridis colorblind-friendly palette
method_colors <- scales::viridis_pal(option = "D")(n_distinct(plot_data$method))

p <- ggplot(plot_data, aes(x = value, fill = method)) +
  geom_density(alpha = 0.6) +
  facet_wrap(~ metric_label, scales = "free", ncol = 1, strip.position = "top") +
  scale_fill_viridis_d(option = "D", name = "Method") +
  labs(x = "", y = "Density") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 14),
    strip.text = element_text(size = 16),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    legend.position = "bottom"
  )

# ============================================================================
# Save plot
# ============================================================================

ggsave(
  filename = here::here("viz", "grid_coverage_densities.png"),
  plot = p,
  width = 10,
  height = 12,
  dpi = 300
)

cat("Plot saved to viz/grid_coverage_densities.png\n")

# Print summary statistics
summary_stats <- plot_data %>%
  group_by(method, metric_label) %>%
  summarize(
    mean = mean(value),
    sd = sd(value),
    median = median(value),
    .groups = "drop"
  )

print(summary_stats)

# ============================================================================
# Part 2: Grid coverage metrics vs sample size
# ============================================================================

cat("\n============================================\n")
cat("Computing metrics vs sample size...\n")
cat("============================================\n")

# Define parameter bounds
alpha_lo <- alpha_bounds[1];  alpha_hi <- alpha_bounds[2]
rho_lo   <- reproduction_number_bounds[1];  rho_hi   <- reproduction_number_bounds[2]
s0_lo    <- s0_bounds[1];  s0_hi    <- s0_bounds[2]
piv_lo   <- piv_bounds[1];  piv_hi   <- piv_bounds[2]
pit_lo   <- pit_bounds[1];  pit_hi   <- pit_bounds[2]

# Find available budget sizes based on actual files (only rep files)
osfd_files <- list.files(here::here('data'), pattern = "^fullOutput_OSFD_basic_budget\\d+_n\\d+_rep\\d+\\.RData$", full.names = TRUE)
sirmap_files <- list.files(here::here('data'), pattern = "^inputs_and_outputs_sirMappings_n\\d+_rep\\d+\\.RData$", full.names = TRUE)
isfd_files <- list.files(here::here('data'), pattern = "^inputs_and_outputs_sirBasicLHS_n\\d+_rep\\d+\\.RData$", full.names = TRUE)

# Extract budgets from files (accounting for _rep## suffix)
osfd_budgets <- as.numeric(gsub(".*budget(\\d+)_n\\d+_rep\\d+\\.RData", "\\1", basename(osfd_files)))
sirmap_n <- as.numeric(gsub(".*_n(\\d+)_rep\\d+\\.RData", "\\1", basename(sirmap_files)))
isfd_n <- as.numeric(gsub(".*_n(\\d+)_rep\\d+\\.RData", "\\1", basename(isfd_files)))

# Use budgets that exist across methods, excluding n=100000 (too large for memory)
BUDGETS <- sort(unique(c(osfd_budgets, sirmap_n, isfd_n)))
BUDGETS <- BUDGETS[BUDGETS < 6000]
cat(sprintf("Found budgets: %s\n", paste(BUDGETS, collapse = ", ")))

k <- 100 # Grid resolution

# Helper function to find all replicate files with given budget
find_files_by_budget <- function(pattern, budget) {
  files <- list.files(here::here('data'), pattern = pattern, full.names = TRUE)

  # Try budget pattern first (for OSFD)
  matching <- grep(sprintf("budget%d_n\\d+_rep\\d+\\.RData$", budget), files, value = TRUE)
  if (length(matching) > 0) return(matching)

  # For ISFD/SIRMap files, try direct n match with rep
  matching <- grep(sprintf("_n%d_rep\\d+\\.RData$", budget), files, value = TRUE)
  if (length(matching) > 0) return(matching)

  return(NULL)
}

# Helper function to calculate grid coverage
calc_grid_coverage <- function(piv, pit, piv_bounds, pit_bounds, k, total_cells) {
  u <- (piv - piv_bounds[1]) / diff(piv_bounds)
  v <- (pit - pit_bounds[1]) / diff(pit_bounds)

  u <- pmin(1, pmax(0, u))
  v <- pmin(1, pmax(0, v))

  ix <- pmin(k, pmax(1, floor(u * k) + 1))
  iy <- pmin(k, pmax(1, floor(v * k) + 1))

  identifiers <- ix + k * (iy - 1)
  n_unique <- length(unique(identifiers))

  return(n_unique / total_cells * 100)
}

# Helper function to calculate fill distance
calc_fill_distance <- function(Y_method, Y_union) {
  dists <- fields::rdist(Y_union, Y_method)
  min_dists <- apply(dists, 1, min)
  return(max(min_dists))
}

# Helper function to calculate mean distance
calc_mean_distance <- function(Y_method, Y_union) {
  dists <- fields::rdist(Y_union, Y_method)
  min_dists <- apply(dists, 1, min)
  return(mean(min_dists))
}

# Create grid of all budget x replicate combinations
budget_rep_grid <- expand.grid(
  budget = BUDGETS,
  rep = 1:20,
  stringsAsFactors = FALSE
)

cat(sprintf("Processing %d budget x replicate combinations in parallel...\n", nrow(budget_rep_grid)))

# Set up parallel cluster with doSNOW for progress bar
# Limit to reasonable number of cores to avoid memory issues
# Use fewer cores to prevent serialization/memory errors
n_cores <- min(10, parallel::detectCores() - 2)
cl <- parallel::makeCluster(n_cores, type = "PSOCK", outfile = "")
doSNOW::registerDoSNOW(cl)

cat(sprintf("Using %d parallel workers\n", n_cores))

# Set up progress bar
pb <- txtProgressBar(max = nrow(budget_rep_grid), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

# Run parallel loop - one iteration per budget-replicate combination
metrics_list <- foreach::foreach(
  i = 1:nrow(budget_rep_grid),
  .packages = c("here", "fields"),
  .combine = rbind,
  .inorder = FALSE,  # Don't require order for better performance
  .errorhandling = "pass",  # Continue even if some iterations fail
  .options.snow = opts
) %dopar% {

  # Wrap in tryCatch to handle errors gracefully
  tryCatch({
    budget <- budget_rep_grid$budget[i]
    rep_num <- budget_rep_grid$rep[i]

    # Initialize result rows
    result_rows <- data.frame()

    # Clear any leftover objects from previous iteration
    gc(verbose = FALSE)

  # Load OSFD Basic for this budget and replicate (rep can be 1-digit or 2-digit)
  osfd_file <- list.files(
    here::here('data'),
    pattern = sprintf("^fullOutput_OSFD_basic_budget%d_n\\d+_rep%d\\.RData$", budget, rep_num),
    full.names = TRUE
  )

  if (length(osfd_file) > 0 && file.exists(osfd_file[1])) {
    load(osfd_file[1])
    y_OSFD <- res$Y_success
    actual_n_osfd <- nrow(y_OSFD)
  } else {
    y_OSFD <- NULL
    actual_n_osfd <- NA
  }

  # Load ISFD (Basic LHS) for this budget and replicate
  isfd_file <- list.files(
    here::here('data'),
    pattern = sprintf("^inputs_and_outputs_sirBasicLHS_n%d_rep%d\\.RData$", budget, rep_num),
    full.names = TRUE
  )

  if (length(isfd_file) > 0 && file.exists(isfd_file[1])) {
    load(isfd_file[1])
    Y_ISFD <- as.matrix(results[, c('piv', 'pit')])
    actual_n_isfd <- nrow(Y_ISFD)
  } else {
    Y_ISFD <- NULL
    actual_n_isfd <- NA
  }

  # Load SIR Maps for this budget and replicate
  sirmap_file <- list.files(
    here::here('data'),
    pattern = sprintf("^inputs_and_outputs_sirMappings_n%d_rep%d\\.RData$", budget, rep_num),
    full.names = TRUE
  )

  if (length(sirmap_file) > 0 && file.exists(sirmap_file[1])) {
    load(sirmap_file[1])
    Y_SIRMap <- as.matrix(results[, c('piv', 'pit')])
    actual_n_sirmap <- nrow(Y_SIRMap)
  } else {
    Y_SIRMap <- NULL
    actual_n_sirmap <- NA
  }

  # Create union of all outputs for this replicate
  # Only include non-NULL methods
  Y_union <- rbind(
    if (!is.null(y_OSFD)) y_OSFD else matrix(ncol = 2, nrow = 0),
    if (!is.null(Y_ISFD)) Y_ISFD else matrix(ncol = 2, nrow = 0),
    if (!is.null(Y_SIRMap)) Y_SIRMap else matrix(ncol = 2, nrow = 0)
  )

  # Skip if union is empty (no methods found for this budget-rep combo)
  if (nrow(Y_union) == 0) return(NULL)

  # Calculate total possible cells from union
  u_all <- (Y_union[, 1] - piv_bounds[1]) / diff(piv_bounds)
  v_all <- (Y_union[, 2] - pit_bounds[1]) / diff(pit_bounds)
  u_all <- pmin(1, pmax(0, u_all))
  v_all <- pmin(1, pmax(0, v_all))
  ix_all <- pmin(k, pmax(1, floor(u_all * k) + 1))
  iy_all <- pmin(k, pmax(1, floor(v_all * k) + 1))
  identifiers_all <- ix_all + k * (iy_all - 1)
  total_possible_cells <- length(unique(identifiers_all))

  # Calculate metrics for OSFD for this replicate (only if file was found)
  if (!is.null(y_OSFD)) {
    result_rows <- rbind(result_rows, data.frame(
      budget = budget,
      n = actual_n_osfd,
      replicate = rep_num,
      method = "OSFD_Basic",
      coverage_pct = calc_grid_coverage(y_OSFD[, 1], y_OSFD[, 2], piv_bounds, pit_bounds, k, total_possible_cells),
      fill_distance = calc_fill_distance(y_OSFD, Y_union),
      mean_distance = calc_mean_distance(y_OSFD, Y_union),
      error_message = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  # Calculate metrics for ISFD for this replicate (only if file was found)
  if (!is.null(Y_ISFD)) {
    result_rows <- rbind(result_rows, data.frame(
      budget = budget,
      n = actual_n_isfd,
      replicate = rep_num,
      method = "Basic_LHS",
      coverage_pct = calc_grid_coverage(Y_ISFD[, 1], Y_ISFD[, 2], piv_bounds, pit_bounds, k, total_possible_cells),
      fill_distance = calc_fill_distance(Y_ISFD, Y_union),
      mean_distance = calc_mean_distance(Y_ISFD, Y_union),
      error_message = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  # Calculate metrics for SIR Maps for this replicate (only if file was found)
  if (!is.null(Y_SIRMap)) {
    result_rows <- rbind(result_rows, data.frame(
      budget = budget,
      n = actual_n_sirmap,
      replicate = rep_num,
      method = "SIR_Maps",
      coverage_pct = calc_grid_coverage(Y_SIRMap[, 1], Y_SIRMap[, 2], piv_bounds, pit_bounds, k, total_possible_cells),
      fill_distance = calc_fill_distance(Y_SIRMap, Y_union),
      mean_distance = calc_mean_distance(Y_SIRMap, Y_union),
      error_message = NA_character_,
      stringsAsFactors = FALSE
    ))
  }

  # Clean up large intermediate objects before returning
  # (result_rows will be preserved as it's the return value)
  rm(y_OSFD, Y_ISFD, Y_SIRMap, Y_union)
  gc(verbose = FALSE)

  result_rows

  }, error = function(e) {
    # On error, return error info
    return(data.frame(
      budget = budget,
      n = NA,
      replicate = rep_num,
      method = "ERROR",
      coverage_pct = NA,
      fill_distance = NA,
      mean_distance = NA,
      error_message = as.character(e$message),
      stringsAsFactors = FALSE
    ))
  })
}

# Close progress bar and cluster
close(pb)
parallel::stopCluster(cl)

# Clean up results
# Use dplyr::bind_rows to properly combine list of data frames
# This handles column types correctly and preserves data frame structure
metrics_by_n <- dplyr::bind_rows(metrics_list)

# Ensure all columns are proper types
metrics_by_n$budget <- as.numeric(metrics_by_n$budget)
metrics_by_n$n <- as.numeric(metrics_by_n$n)
metrics_by_n$replicate <- as.numeric(metrics_by_n$replicate)
metrics_by_n$method <- as.character(metrics_by_n$method)
metrics_by_n$coverage_pct <- as.numeric(metrics_by_n$coverage_pct)
metrics_by_n$fill_distance <- as.numeric(metrics_by_n$fill_distance)
metrics_by_n$mean_distance <- as.numeric(metrics_by_n$mean_distance)
metrics_by_n$error_message <- as.character(metrics_by_n$error_message)

# Separate errors from valid data
errors <- metrics_by_n[metrics_by_n$method == "ERROR" & !is.na(metrics_by_n$method), , drop = FALSE]
metrics_by_n <- metrics_by_n[metrics_by_n$method != "ERROR" & !is.na(metrics_by_n$method), , drop = FALSE]

# Remove any duplicate rows (shouldn't happen, but just in case)
metrics_by_n <- unique(metrics_by_n)

cat("Metrics by n calculation complete\n")
cat(sprintf("Total rows: %d replicates across %d budgets and %d methods\n",
            nrow(metrics_by_n),
            n_distinct(metrics_by_n$budget),
            n_distinct(metrics_by_n$method)))

# Report any errors
if (nrow(errors) > 0) {
  cat(sprintf("\n*** WARNING: %d budget-replicate combinations encountered errors ***\n", nrow(errors)))
  cat("\nError details:\n")
  print(errors[, c("budget", "replicate", "error_message")])
  cat("\n")
} else {
  cat("No errors encountered during processing\n")
}

# Reshape for plotting
metrics_by_n_long <- metrics_by_n %>%
  pivot_longer(
    cols = c(coverage_pct, fill_distance, mean_distance),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    # Rename methods for display (matching wang_helper_functions.R naming)
    method = case_when(
      method == "Basic_LHS" ~ "LHS",
      method == "OSFD_Basic" ~ "Failure-aware OSFD",
      method == "SIR_Maps" ~ "Inverse SIR Maps",
      TRUE ~ method
    ),
    # Factor method to ensure consistent ordering (and thus consistent colors)
    method = factor(method, levels = c("LHS", "Failure-aware OSFD", "Inverse SIR Maps")),
    metric_label = case_when(
      metric == "coverage_pct" ~ "Grid Coverage",
      metric == "fill_distance" ~ "Fill Distance",
      metric == "mean_distance" ~ "Mean Distance"
    ),
    metric_label = factor(
      metric_label,
      levels = c("Grid Coverage", "Fill Distance", "Mean Distance")
    )
  )

# Calculate central 90% region (5th to 95th percentile) and median for each method-budget-metric combo
# Group by budget (not n) since n might vary slightly across replicates
metrics_summary <- metrics_by_n_long %>%
  group_by(budget, method, metric_label) %>%
  summarize(
    n = median(n),  # Use median n for x-axis positioning
    median_value = median(value),
    lower_90 = quantile(value, 0.05),
    upper_90 = quantile(value, 0.95),
    .groups = "drop"
  )

cat(sprintf("Computed central 90%% regions for %d method-budget-metric combinations\n", nrow(metrics_summary)))

# Define consistent color mapping across all plots
# Using viridis palette with 5 distinct colors, but only using first 3 positions
# Match the color positions used in other plot scripts
viridis_colors <- scales::viridis_pal(option = "D")(5)
method_colors <- c(
  "LHS" = viridis_colors[1],                    # Purple/dark blue (color 1)
  "Failure-aware OSFD" = viridis_colors[2],     # Blue/cyan (color 2)
  "Inverse SIR Maps" = viridis_colors[3]        # Teal/green (color 3)
)

# Create line plots for each metric with shaded central 90% region
p_n <- ggplot() +
  # Add shaded ribbon for central 90% region
  geom_ribbon(
    data = metrics_summary,
    aes(x = n, ymin = lower_90, ymax = upper_90, fill = method, group = method),
    alpha = 0.25
  ) +
  # Add median line
  geom_line(
    data = metrics_summary,
    aes(x = n, y = median_value, color = method, group = method),
    linewidth = 1.2
  ) +
  # Add points at median
  geom_point(
    data = metrics_summary,
    aes(x = n, y = median_value, color = method, group = method),
    size = 3
  ) +
  facet_wrap(~ metric_label, scales = "free_y", ncol = 1, strip.position = "top") +
  scale_x_log10() +
  scale_color_manual(values = method_colors, name = "Method") +
  scale_fill_manual(values = method_colors, name = "Method") +
  labs(x = "Sample Size (n)", y = "Metric Value") +
  theme_bw() +
  theme(
    axis.title = element_text(size = 22),      # Axis titles (x and y) - a bit larger
    axis.text = element_text(size = 20),       # Axis tick labels - MUCH bigger
    strip.text = element_text(size = 18),      # Facet labels - slightly bigger
    legend.title = element_text(size = 20),    # Legend title - MUCH bigger
    legend.text = element_text(size = 18),     # Legend text - MUCH bigger
    legend.position = "bottom"
  )

# Save sample size plot separately
ggsave(
  filename = here::here("viz", "grid_coverage_vs_n.png"),
  plot = p_n,
  width = 10,
  height = 12,
  dpi = 300
)

cat("Sample size plot saved to viz/grid_coverage_vs_n.png\n")


# ============================================================================
# Combine density and sample size plots
# ============================================================================

combined_plot <- p | p_n

combined_plot <- combined_plot & theme(
  axis.title = element_text(size = 18),
  axis.text = element_text(size = 14),
  strip.text = element_text(size = 16),
  legend.title = element_text(size = 16),
  legend.text = element_text(size = 14),
  legend.position = "bottom"
)

ggsave(
  filename = here::here("viz", "grid_coverage_combined.png"),
  plot = combined_plot,
  width = 20,
  height = 12,
  dpi = 300
)

cat("\nCombined plot saved to viz/grid_coverage_combined.png\n")
cat("\n============================================\n")
cat("Summary Statistics (Central 90% Region)\n")
cat("============================================\n")
print(metrics_summary)

cat("\n============================================\n")
cat("All Replicate Data\n")
cat("============================================\n")
print(metrics_by_n)
