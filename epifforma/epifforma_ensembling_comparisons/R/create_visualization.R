
# Fitting Synthetic Method of Analogues and Testing it against CovidHub Ensemble Models
## Authors: GC Gibson and AC Murph
library(lhs)
library(mgcv)
library(ggplot2)
library(plyr)
library(data.table)
library(LearnBayes)
library(LaplacesDemon)
library(dplyr)
library(spatstat)
library(BASS)
library(GPfit)
library(nnet)
library(dplyr)
library(KernelKnn)
library(ggridges)
library(utils)
library(parallel)
library(doParallel)
library(grid)
library(forecast)
library(ggplot2)
library(dplyr)
library(tidyr)
library(purrr)
library(patchwork)

source(here::here("epifforma","epifforma_ensembling_comparisons", "R", "smoa_helpers.R"))
code_path = here::here("epifforma","epifforma_ensembling_comparisons", "R")
data_path = here::here("epifforma","epifforma_ensembling_comparisons", "data")
source(here::here("epifforma","process_data", "epi_functions.R"))
source(paste0(code_path,"/alternative_ensemblers.R"))

ncores <- 51
sim_idx <- 1
h       <- 4

# These are all from the bayesian optimization.
num_curves <- 18387
k <- 5
closest <- 4422
lower_CI_scale <- 1
upper_CI_scale <- 1
dispersion_forecast <- 1
mle_lower_bound <- 1

###########
# Make a directory to hold the results for each state, using this run's particular set of hyperparameters.
name_of_change <- paste("k", k, "num_curves", num_curves, "closest", closest, "dispersion", round(dispersion_forecast*10000), 'mlebound', round(mle_lower_bound*10000), sep = "_")
state_log_directory <- paste(name_of_change, "_state_records", sep = "")


############################################
### Generate Synthetic List for Training ###
############################################
sockettype <- "PSOCK"


quantiles <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
### read in stored truth data with as_of

# Coverage Data:
coverage_data <- NULL

# We pre-built this data file to cut on api calls to github.
truth_as_of_tot                 <- read.csv(paste0(data_path,"/tdat_list_tot_weekly.csv"))

### Iterate through the states and calculate the MAE and WIS for the sMOA forecast.
mse_df_list                     <- list()
count_list                      <- 1
curr_state <- "California"
mle_start_value <- 20

location = curr_state
#print(paste("working on location", location))  

##### true weekly data as of a final data of reporting
state_truth_data_file <- paste(data_path, '/state_truths/', gsub(" ", "", location), '.RData', sep = "")
if(!file.exists(state_truth_data_file)){
  truth_weekly                  <- load_truth(truth_source = "JHU",hub = "US", target_variable = "inc case", as_of= "2023-03-04",locations =location,data_location="covidData")
  save(truth_weekly, file = state_truth_data_file)
}else{
  load(state_truth_data_file)
}

#### get the list of forecast dates
fcast_dates_to_match          <- unique(truth_weekly$target_end_date)

# We need to collect the one-step-ahead forecasts for all available data (potentially not just
# from august 15th onwards).
one_step_ahead_forecasts <- c()
two_step_ahead_forecasts <- c()
three_step_ahead_forecasts <- c()
four_step_ahead_forecasts <- c()

#### subset to "2020-08-15" when the revised API started up
fcast_dates_to_match          <- fcast_dates_to_match[which(fcast_dates_to_match > as.Date("2020-08-15"))]

#### subset big as of data frame to this specific location
truth_as_of_tot_loc           <- truth_as_of_tot[truth_as_of_tot$location_name == location ,]

  
state_coverage_location = paste(data_path, "/", state_log_directory, sep = "")
curr_state = "California"
load(paste0(state_coverage_location, "/", curr_state, "weights_list.RData"))
load(paste0(state_coverage_location, "/", curr_state, "preds_list.RData"))
load(paste0(state_coverage_location, "/", curr_state, "true_values.RData"))

###############################################
# Helper: compute MAE over time
###############################################

mae_over_time <- function(pred, truth) {
  abs(pred - truth)
}

###############################################
# Calculate MAE time series for all ensembles
###############################################
preds_list$ewa <- NULL
# preds_list$equal_wt <- NULL

plots <- vector("list", h)  # to store one ggplot per horizon

for (horizon_num in 1:h) {
  idx <- seq(from = horizon_num,
             to   = length(true_values),
             by   = h)
  
  truth <- true_values[idx]
  horizon_dates <- fcast_dates_to_match  # <- dates for this horizon
  
  # 1. Compute cumulative / rolling-average MAE for each ensemble
  cum_mae_list <- list()
  
  for (nm in names(preds_list)) {
    pred <- preds_list[[nm]][, horizon_num]
    abs_err <- abs(pred - truth)
    
    # cumulative mean MAE up to each time t
    cum_mae <- cumsum(abs_err) / seq_along(abs_err)
    
    cum_mae_list[[nm]] <- cum_mae
  }
  
  # 2. Turn cum_mae_list into a tidy data frame for ggplot
  cum_mae_df <- imap_dfr(
    cum_mae_list,
    ~ data.frame(
      date   = horizon_dates,      # use dates instead of time index
      cum_mae = .x,
      method  = .y,
      horizon = horizon_num
    )
  )
  
  # 3. Build the ggplot for this horizon
  p <- ggplot(cum_mae_df, aes(x = date, y = cum_mae, color = method)) +
    geom_line(linewidth = 1) +
    labs(
      x = "Date",
      y = "Cumulative Average MAE",
      title = paste0(curr_state, ": Cumulative (Rolling Average) MAE"),
      subtitle = paste("Horizon", horizon_num)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "top"
    ) +
    scale_x_date(date_breaks = "3 months", date_labels = "%Y-%m")  # tweak as you like
  
  plots[[horizon_num]] <- p
}

combined_plot <- wrap_plots(plots, ncol = 2)
combined_plot
