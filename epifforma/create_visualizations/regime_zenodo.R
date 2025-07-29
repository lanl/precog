################################################################################
# Shapelet-Based Regime Classification for epiFFORMA Forecasting Evaluation
#
# This script implements shapelet-based methods for classifying epidemic curves
# into different regimes (increasing, decreasing, surge, etc.) and evaluates
# model performance across these regimes.
#
# Adapted from: https://github.com/Satwant-Singh-ADS/Shapelet_Methods/
# Original Paper: https://ieeexplore.ieee.org/document/10020895
################################################################################

# Required packages
library(tidyverse)  # For data manipulation and visualization
library(zoo)        # For rolling window operations
library(roll)       # For additional rolling window functions
library(this.path)
setwd(paste0(this.path::here(), '/../'))

################################################################################
# FUNCTIONS
################################################################################

#' Assign epidemic regime to a time series
#'
#' This function takes a time series and determines the best matching
#' epidemic regime based on shapelet similarity.
#'
#' @param trial A dataframe containing a column 'truth' with time series values
#' @return Character string indicating the assigned regime
assign_regime2 <- function(trial) {
  vec <- trial$truth
  regime <- return_best_shapelet_pearson(vec)
  return(regime)
}

#' Calculate similarity between two vectors using Pearson correlation
#'
#' @param vector1 First vector for comparison
#' @param vector2 Second vector for comparison
#' @return Numeric value representing similarity (-1 to 1)
similarity_matrix <- function(vector1, vector2) {
  # Handle near-constant vectors to avoid correlation calculation issues
  if (sd(vector1) < 1e-100 || sd(vector2) < 1e-100) {
    similarity_value <- 0
  } else {
    similarity_value <- cor(vector1, vector2, method = 'pearson')
  }
  similarity_value
}

#' Find the best matching shapelet for a time series
#'
#' Determines which predefined shapelet pattern best matches the input vector
#' based on Pearson correlation, accounting for slope characteristics.
#'
#' @param vector Time series vector to classify
#' @param slope_thres Threshold for determining flatness, default 0.0005
#' @return Character string indicating the best matching shapelet name
return_best_shapelet_pearson <- function(vector, slope_thres = 0.0005) {
  corrs <- return_all_shapelet_pearson(vector, slope_thres)
  scenario <- which.max(corrs)
  shapelet_standard_names[scenario]
}

#' Compute similarity scores for all standard shapelet patterns
#'
#' Calculates correlation between input vector and each predefined shapelet,
#' adjusting for flatness of the input vector.
#'
#' @param vector Time series vector to compare against shapelets
#' @param slope_thres Threshold for determining flatness, default 0.0005
#' @return Numeric vector of similarity scores for each shapelet
return_all_shapelet_pearson <- function(vector, slope_thres = 0.0005) {
  corrs <- c()
  beta <- -log(0.1) / slope_thres
  m0 <- 0
  slope <- mean(abs(diff(vector)))
  
  # Calculate flatness measure based on slope
  if (slope < m0) {
    flatness <- 1
  } else {
    flatness <- exp(-beta * (slope - m0))
  }
  
  # Calculate correlation with each shapelet, adjusting for flatness
  for (i in 1:length(shapelet_standard_array)) {
    if (all(shapelet_standard_array[[i]] == 0)) {
      score <- 2 * flatness - 1
    } else {
      score <- (1 - flatness) * similarity_matrix(shapelet_standard_array[[i]], vector)
    }
    corrs[i] <- score
  }
  corrs
}

#' Apply LOESS smoothing to time series data
#'
#' @param truth Vector of observed values to smooth
#' @param x Vector of time points corresponding to truth values
#' @param span Controls the degree of smoothing (0-1)
#' @return Vector of smoothed values
do_smooth <- function(truth, x, span = 0.25) { 
  # Use larger span for shorter time series to avoid overfitting
  if(length(truth) < 40) {
    span <- 0.55
  }
  lo_tmp <- loess(truth ~ x, span = span)
  smooth_cases <- predict(lo_tmp, x)
  # Ensure non-negative values (cases can't be negative)
  smooth_cases[smooth_cases < 0] <- 0 
  return(smooth_cases)
}

################################################################################
# SHAPELET CONFIGURATION
################################################################################

# Hyperparameters for shapelet definition
vector_length <- c(0, 4)   # Look 4 weeks ahead in future while defining shapelet
number_of_shapelets <- 6
shapelet_length <- vector_length[1] + vector_length[2]

# Initialize shapelet arrays
shapelet_standard_array <- vector("list", number_of_shapelets)
shapelet_standard_array <- lapply(1:number_of_shapelets, function(x) rep(0, shapelet_length))

# Define names for each shapelet pattern
shapelet_standard_names <- c("Inc", "Dec", "Surge", "Near Peak", "Past Peak", "Flat")

# Uncomment to verify configuration consistency
#stopifnot(length(shapelet_standard_names) == Number_of_shapelets, 
#          "Size of array mismatch for shapelet_standard_names and Number_of_shapelets")

# Define shapelet patterns
shapelet_standard_array[[1]] <- c(1, 2, 3, 4)            # Increasing
shapelet_standard_array[[2]] <- c(4, 3, 2, 1)            # Decreasing
shapelet_standard_array[[3]] <- c(1, 2, 4, 8)            # Surge
shapelet_standard_array[[4]] <- c(-1, -0.5, -0.25, -0.125)  # Near peak
shapelet_standard_array[[5]] <- c(-1, -2, -4, -8)        # Past peak
# Note: Flat pattern (index 6) remains as zeros

################################################################################
# DATA PROCESSING AND ANALYSIS
################################################################################

# Find all evaluation files, excluding synthetic data
list_rds <- list.files("evaluate_model/evaluation/", 
                       full.names = TRUE, pattern = ".RDS")
diseases <- list_rds[str_detect(list_rds, negate = TRUE, pattern = "synthetic")]

# Initialize list to store error metrics for each disease
error_list <- vector(mode = "list", length = length(diseases))

# Process each disease dataset
for(i in 1:length(diseases)) { 
  
  # Read disease data
  disease <- diseases[i]
  disease_list <- readRDS(disease)
  
  # Extract plot data 
  df <- disease_list$plot_df
  
  # Extract and smooth observed data
  cases_smooth <- df %>% 
    filter(type == "obs") %>% 
    group_by(geography, last_obs_time) %>% 
    filter(x > (last_obs_time - 40)) %>% 
    mutate(smoothed_cases = do_smooth(truth = truth, x = x)) 
  
  # Extract forecast/future data
  cases_smooth_sub <- cases_smooth %>% 
    group_by(geography, last_obs_time) %>%
    filter(x > last_obs_time)
  
  # Assign regime to each forecast period
  trial_regime <- plyr::ddply(cases_smooth_sub, 
                              .variables = c("geography", "last_obs_time"), 
                              .fun = assign_regime2)
  
  # Calculate error metrics for each forecast
  error <- df %>% 
    left_join(trial_regime) %>% 
    filter(type == "fcst") %>% 
    rename("regime" = V1) %>% 
    group_by(geography, model, h, regime) %>% 
    mutate(mae = abs(fcst - truth),
           se = abs(fcst - truth)^2) 
  
  # Summarize error metrics
  error <- error %>% 
    summarize(mean_mae = mean(mae, na.rm = TRUE),
              rmse = sqrt(mean(se))) 
  
  # Calculate relative error metrics
  rel_error <- error %>% 
    group_by(geography, regime, h) %>% 
    mutate(rel_mae = (mean_mae - min(mean_mae)) / (max(mean_mae) - min(mean_mae)),
           rel_rmse = (rmse - min(rmse)) / (max(rmse) - min(rmse))) 
  
  error_list[[i]] <- rel_error
}

# Name the list elements with disease names
disease_names <- str_remove(basename(diseases), pattern = "_order_multierror.RDS")
names(error_list) <- disease_names

# Combine all disease results into a single dataframe
error_df <- bind_rows(error_list, .id = "disease")

################################################################################
# VISUALIZATION
################################################################################

# Heatmap of relative MAE by disease, model, horizon, and regime
error_df %>% 
  ggplot(aes(x = disease, y = model, fill = rel_mae)) + 
  geom_tile() + 
  facet_grid(h ~ regime) + 
  scale_fill_gradient2(midpoint = 0.5, high = "red", low = "blue") + 
  theme(axis.text.x = element_text(angle = 60, hjust = 1, vjust = 1),
        text = element_text(size = 16)) +
  labs(title = "Relative MAE by Disease, Model, Horizon and Regime",
       x = "Disease", 
       y = "Model",
       fill = "Relative MAE")

# Summarize errors across horizons
error_df_summary <- error_df %>% 
  group_by(disease, regime, model) %>% 
  dplyr::summarize(mean_rel_mae = mean(rel_mae),
                   mean_rel_rmse = mean(rel_rmse))

# Create boxplot visualization of error metrics by regime
metric_labels <- c(mean_rel_rmse = "Normalized RMSE", mean_rel_mae = "Normalized MAE")
error_df_summary$regime <- factor(error_df_summary$regime, 
                                  levels = c("Surge", "Inc", "Near Peak", "Past Peak", "Dec", "Flat"))

plot_error_boxplot <- error_df_summary %>% 
  mutate(color = ifelse(model == "epifforma", "blue", "grey")) %>% 
  pivot_longer(names_to = "metric", values_to = "value", cols = c("mean_rel_rmse", "mean_rel_mae")) %>% 
  ggplot(aes(x = value, y = model, fill = color)) + 
  geom_boxplot(alpha = 0.5) + 
  geom_point(alpha = 0.2) +
  facet_grid(metric ~ regime, labeller = labeller(metric = metric_labels)) + 
  theme_bw(base_size = 16) +
  scale_fill_manual(values = c("blue", "grey")) +
  guides(fill = "none") + 
  labs(x = "Normalized Error Value", 
       y = "Model",
       title = "Model Performance by Epidemic Regime")
cowplot::save_plot(plot_error_boxplot, filename = "create_visualizations/error_summary.pdf", base_height = 8)
