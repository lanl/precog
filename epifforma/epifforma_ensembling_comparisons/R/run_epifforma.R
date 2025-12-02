# Fitting Synthetic Method of Analogues and Testing it against CovidHub Ensemble Models
## Authors: GC Gibson and AC Murph
library(lhs)
library(mgcv)
library(ggplot2)
library(plyr)
library(data.table)
library(LearnBayes)
library(LaplacesDemon)
#library(covidHubUtils)
library(dplyr)
library(spatstat)
library(BASS)
#library(covidcast)
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

### Read in embedding matrices
embed_mat_X <- data.table::fread(file=paste0(data_path,"/embed_mat/embed_mat_X.csv"))
embed_mat_y <- data.table::fread(file=paste0(data_path,"/embed_mat/embed_mat_y.csv"))

embed_mat_X_deriv <- data.table::fread(file=paste0(data_path,"/embed_mat/embed_mat_X_deriv.csv"))
embed_mat_y_deriv <- data.table::fread(file=paste0(data_path,"/embed_mat/embed_mat_y_deriv.csv"))

quantiles <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
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

parfctn = function(x){
  library(covidHubUtils)
  library(mgcv)
  library(collapse)
  library(dplyr)
  quantiles <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
  ### read in stored truth data with as_of
  
  # Coverage Data:
  coverage_data <- NULL
  
  # We pre-built this data file to cut on api calls to github.
  truth_as_of_tot                 <- read.csv(paste0(data_path,"/tdat_list_tot_weekly.csv"))
  
  ### Iterate through the states and calculate the MAE and WIS for the sMOA forecast.
  mse_df_list                     <- list()
  count_list                      <- 1
  curr_state <- state.name[x]
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
  
  ### hold mse 
  mse_df                        <- NULL 
  
  
  # Set up for all the other ensemble weighting techinques
  num_components = 9
  T_max = length(fcast_dates_to_match)
  
  ewa = list()
  rls = list()
  ridge = list()
  kalman = list()
  roll = list()
  
  for(horizon_num in 1:h){
    ewa[[horizon_num]]    <- make_ewa_ensemble(num_components, eta = 0.05)
    rls[[horizon_num]]    <- make_rls_ensemble(num_components, lambda = 0.99, delta = 1e3)
    ridge[[horizon_num]]  <- make_ridge_rls_ensemble(num_components, tau = 1.0, delta = 1.0)
    kalman[[horizon_num]] <- make_kalman_ensemble(num_components, Q_scale = 1e-4, R = 1.0)
    roll[[horizon_num]]   <- make_rolling_reg_ensemble(num_components, window = 40, lambda = 0.1)
  }
  
  # ewa    <- make_ewa_ensemble(num_components, eta = 0.05)
  # rls    <- make_rls_ensemble(num_components, lambda = 0.99, delta = 1e3)
  # ridge  <- make_ridge_rls_ensemble(num_components, tau = 1.0, delta = 1.0)
  # kalman <- make_kalman_ensemble(num_components, Q_scale = 1e-4, R = 1.0)
  # roll   <- make_rolling_reg_ensemble(num_components, window = 40, lambda = 0.1)
  
  # Storage for weights and predictions
  ewa_w    <- array(NA_real_, dim = c( T_max, num_components, h) )
  rls_w    <- array(NA_real_, dim = c( T_max, num_components, h))
  ridge_w  <- array(NA_real_, dim = c( T_max, num_components, h))
  kalman_w <- array(NA_real_, dim = c( T_max, num_components, h))
  roll_w   <- array(NA_real_, dim = c( T_max, num_components, h))
  epifforma_w   <- array(NA_real_, dim = c( T_max, num_components, h))
  eq_w   <- array(1/9, dim = c( T_max, num_components, h))
  
  ewa_pred    <- matrix(NA_real_, T_max, h)
  rls_pred    <- matrix(NA_real_, T_max, h)
  ridge_pred  <- matrix(NA_real_, T_max, h)
  kalman_pred <- matrix(NA_real_, T_max, h)
  roll_pred   <- matrix(NA_real_, T_max, h)
  epifforma_pred   <- matrix(NA_real_, T_max, h)
  eqw_pred   <- matrix(NA_real_, T_max, h)
  
  true_values = c()
  
  
  ### iterate through forecast dates
  for (fcast_date_idx in 1:(length(fcast_dates_to_match))){
    print(paste(fcast_date_idx, 'of', length(fcast_dates_to_match)))
    #truth_weekly                  <- load_truth(truth_source = "JHU",hub = "US", target_variable = "inc case", 
    #					  as_of= as.character(as.Date(fcast_dates_to_match[fcast_date_idx])), locations =location,data_location="covidData")

    #### grab the date formatted
    fcast_date                  <- fcast_dates_to_match[fcast_date_idx]
    
    #### subset to data that was available at forecast date
    truth_as_of                 <- truth_as_of_tot_loc[truth_as_of_tot_loc$as_of == fcast_date ,]
    #### double check its ordered properly
    truth_as_of                 <- truth_as_of[order(truth_as_of$target_end_date),]
    
    #### make sure that the target end date is less than or equal to fcast_date
    data_till_now               <- truth_as_of[truth_as_of$target_end_date <= fcast_date,]
    data_till_now$value         <- pmax(1,data_till_now$value)
    data_till_now$t             <- 1:nrow(data_till_now)
    
    #### light smoothing and differencing and get last k 
    data_till_now_smoothed      <- gam(value~ s(t,k=round(nrow(data_till_now)/2)),data=data_till_now)$fitted.values
    
    #### could use gam smoother 
    to_match_in_moa             <- tail((data_till_now_smoothed),k+1)
   
    ###############################################
    #### THIS IS WHERE I GET THE EPIFFORMA FORECASTS 

    data_future                 <- tail(truth_weekly[truth_weekly$target_end_date <= (fcast_date + h*7),]$value,h) # Note that that '7' is for the length of a week (7 days).

    ######################################
    ### Step 1: Simulate a time series ###
    ######################################
    ts_scaled = data_till_now_smoothed
    info_packet = list(ts = ts_scaled,
                       ts_time_cadence = 'weekly',
                       ts_scale = 'counts')
    # plot(ts_scaled, xlab = 'Weeks', ylab = 'Cases')
    
    
    ################################################
    ### Step 2: Estimate features and components ###
    ################################################
    
    ## handle outliers
    info_packet$ts = as.numeric(handle_outliers(info_packet))
    
    ## distribute zeros
    info_packet$ts = as.numeric(distribute_zeros(info_packet))
    
    ## compute the features
    features <- make_features(info_packet, h=4)
    
    ## compute the components
    components <- make_components(info_packet, h=4)
    
    ## compute the components
    component_intervals <- make_component_intervals(info_packet, h=4)
    
    
    ############################################################
    ### Step 3: Generate epiFFORMA predictions and intervals ###
    ############################################################
    
    
    ### Get Point Estimates
    feature2wt <- readRDS(paste0(data_path,"/fitted_lgb_models_order_multierror.RDS"))
    pred_wts_list <- list()
    for(j in 1:length(feature2wt)){
      pred_wts_list[[j]] <- predict(feature2wt[[j]], newdata = as.matrix(features))
    }
    pred_wts <- apply(simplify2array(pred_wts_list), 1:2, mean)
    pred_wts_sd <- apply(simplify2array(pred_wts_list), 1:2, sd)
    pred_wts[pred_wts<0.01] = 0 
    pred_wts = t(apply(pred_wts, 1, FUN = function(x){x/sum(x)}))
    components_reordered = NULL
    for(j in 1:nrow(components)){
      components_reordered = rbind(components_reordered, sort(unlist(components[j,])))
    }
    fcst_df <- data.frame(h = 1:4,epifforma = rowSums(pred_wts*components_reordered))
    
    for(horizon_num in 1:h){
      x_t <- components_reordered[horizon_num, ]
      y_t <- data_future[horizon_num]
      true_values = c(true_values, y_t)
      t = fcast_date_idx
      
      # 0. Grab epifforma stuff
      epifforma_pred[t, horizon_num] <- fcst_df[horizon_num,2]
      epifforma_w[t, , horizon_num] <- pred_wts[horizon_num,]
      
      # 1. EWA
      ewa_pred[t,horizon_num] <- sum(ewa[[horizon_num]]$weights * x_t)
      ewa[[horizon_num]] <- update_ewa_ensemble(ewa[[horizon_num]], x_t, y_t)
      ewa_w[t, ,horizon_num] <- ewa[[horizon_num]]$weights
      
      # 2. RLS
      rls_pred[t, horizon_num] <- sum(rls[[horizon_num]]$w * x_t)
      rls[[horizon_num]] <- update_rls_ensemble(rls[[horizon_num]], x_t, y_t)
      rls_w[t, , horizon_num] <- rls[[horizon_num]]$w
      
      # 3. Ridge RLS
      ridge_pred[t, horizon_num] <- sum(ridge[[horizon_num]]$w * x_t)
      ridge[[horizon_num]] <- update_ridge_rls_ensemble(ridge[[horizon_num]], x_t, y_t)
      ridge_w[t, , horizon_num] <- ridge[[horizon_num]]$w
      
      # 5. Kalman
      kalman_pred[t, horizon_num] <- sum(kalman[[horizon_num]]$w * x_t)
      kalman[[horizon_num]] <- update_kalman_ensemble(kalman[[horizon_num]], x_t, y_t)
      kalman_w[t, , horizon_num] <- kalman[[horizon_num]]$w
      
      # 6. Rolling window regression
      roll_pred[t, horizon_num] <- sum(roll[[horizon_num]]$w * x_t)
      roll[[horizon_num]] <- update_rolling_reg_ensemble(roll[[horizon_num]], x_t, y_t)
      roll_w[t, ,horizon_num] <- roll[[horizon_num]]$w
      
      eqw_pred[t, horizon_num] = mean(x_t)
    }
  }
  
  #### create the data frame and return
  weights_list = list()
  weights_list[["ewa"]] = ewa_w
  weights_list[["rls"]] = rls_w
  weights_list[["ridge"]] = ridge_w
  weights_list[["kalman"]] = kalman_w
  weights_list[["roll"]] = roll_w
  weights_list[["epifforma"]] = epifforma_w
  
  preds_list = list()
  preds_list[["ewa"]] = ewa_pred
  preds_list[["rls"]] = rls_pred
  preds_list[["ridge"]] = ridge_pred
  preds_list[["kalman"]] = kalman_pred
  preds_list[["roll"]] = roll_pred
  preds_list[["epifforma"]] = epifforma_pred
  preds_list[["equal_wt"]] = eqw_pred
  
  state_coverage_location = paste("data/", state_log_directory, sep = "")
  save(weights_list, file = paste0(state_coverage_location, "/", curr_state, "weights_list.RData") )
  save(preds_list, file = paste0(state_coverage_location, "/", curr_state, "preds_list.RData") )
  save(true_values, file = paste0(state_coverage_location, "/", curr_state, "true_values.RData") )
}


# If the directory exists, change that here.
directory_name = paste(data_path, "/", state_log_directory, sep = "")
if(file.exists(directory_name)) unlink(directory_name, recursive = TRUE)

dir.create(directory_name)
sockettype <- "PSOCK"

## Uncomment this to work with a simple example (one run).
# parfctn(3)

cl <- parallel::makeCluster(spec = ncores,type = sockettype) #, outfile=""
setDefaultCluster(cl)
registerDoParallel(cl)
sim_ts <- foreach(i=1:50,
                  .verbose = T)%dopar%{
		    print(i)
                    parfctn(i)
                  }
stopCluster(cl)


state_coverage_location = paste("data/", state_log_directory, sep = "")
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
  truth <- true_values[ seq(from = horizon_num, 
                            to   = length(true_values), 
                            by   = h) ]
  
  # 1. Compute cumulative / rolling-average MAE for each ensemble
  cum_mae_list <- list()
             # keep your original exclusion
  
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
      time = seq_along(.x),
      cum_mae = .x,
      method = .y,
      horizon = horizon_num
    )
  )
  
  # 3. Build the ggplot for this horizon
  p <- ggplot(cum_mae_df, aes(x = time, y = cum_mae, color = method)) +
    geom_line(linewidth = 1) +
    labs(
      x = "Time",
      y = "Cumulative Average MAE",
      title = paste0(curr_state, ": Cumulative (Rolling Average) MAE"),
      subtitle = paste("Horizon", horizon_num)
    ) +
    theme_minimal() +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "top"
    )
  
  plots[[horizon_num]] <- p
}

# 4. Combine with patchwork in a 2x2 grid
combined_plot <- wrap_plots(plots, ncol = 2)
combined_plot