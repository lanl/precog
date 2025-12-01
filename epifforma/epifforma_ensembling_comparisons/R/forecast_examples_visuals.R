# Visualizations for the sMOA Paper; example on early-epidemic California.
## Author: Alexander C. Murph
## Date: September 2024
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
library(tidybayes)
library(ggridges)
library(utils)
library(parallel)
library(doParallel)
library(covidHubUtils)
library(mgcv)
library(collapse)
library(dplyr)
setwd("~/GitLab/smoa")
source("R/smoa_helpers.R")
source('R/vecchia_scaled.R')

ncores <- 51
lopez_models <- c("BPagano-RtDriven", 'CEID-Walk', 'CovidAnalytics-DELPHI', 'COVIDhub-baseline', 'COVIDhub-4_week_ensemble',
                  'COVIDhub-trained_ensemble', 'Covid19Sim-Simulator', 'CU-select', 'FAIR-NRAR', 'FRBSF_Wilson-Econometric',
                  'IEM_MED-CovidProject', 'IowaStateLW-STEM', 'JHUAPL-Bucky', 'JHU_CSSE-DECOM', 'JHU_IDD-CovidSP',
                  'Karlen-pypm', 'LNQ-ens1', 'LANL-GrowthRate', 'Microsoft-DeepSTIA', 'MOBS-GLEAM_COVID',
                  'RobertWalraven-ESG', 'UCLA-SuEIR', 'USC-SI_kJalpha', 'UMass-MechBayes', 'UVA-Ensemble')

###########
# The following is only for when we run the grid search to test several different
# hyperparameters:
sim_idx                         <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(paste("Running slurm job:", sim_idx))
if(is.na(sim_idx)){
  sim_idx                       <- 1
}
curr_state                      <- state.name[sim_idx]


num_of_each_curve <- c(20000, 30000, 40000, 50000)
ks <- c(4, 6,7,8)
epsilons <- c(5e-4, 5e-3, 5e-5, 5e-2)
closest_ids_vec <- c(1000, 2000, 3000, 5000)

dispersion_forecast_scaling <- seq(from=0.1, to=5, length.out = 64)
dispersion_forecast <- dispersion_forecast_scaling[sim_idx]
h                               <- 4

# These are all from the bayesian optimization.
upper_mle_bound <- 2000
num_curves <- 18387
k <- 7
closest <- 4422
lower_CI_scale <- 1
upper_CI_scale <- 1
# dispersion_forecast <- 91.32496 / 3
# mle_lower_bound <- 8.156253

dispersion_forecast <- 1
mle_lower_bound <- 1
lower_mle_bound <- mle_lower_bound

###########
# Make a directory to hold the results for each state, using this run's particular set of hyperparameters.
name_of_change <- paste("k", k, "num_curves", num_curves, "closest", closest, "dispersion", round(dispersion_forecast*10000), sep = "_")
state_log_directory <- paste(name_of_change, "_state_records", sep = "")

############################################
### Generate Synthetic List for Training ###
############################################
types_of_curves <- c("sir_rollercoaster", "sir_rollercoaster_wiggle", 'seasonal')
N <- num_curves*length(types_of_curves)

sockettype <- "PSOCK"

load(paste("data/synthetic_logs/synthetic_simidx_", sim_idx, "_num_curves_", num_curves, ".RData", sep = ""))
print("Creating the embedding matrix.")

### read in embedding mat design matrix and outcome matrix
embed_mat                       <- create_embed_matrix(sim_ts,h=(h+1),k=(k+1))
X                               <- embed_mat[[1]]
y                               <- embed_mat[[2]]

### convert to differences
X_diff                          <- t(apply(X,1,diff))
y_diff                          <- t(apply(y,1,diff))

quantiles <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
## read in stored truth data with as_of

# We pre-built this data file to cut on api calls to github.
truth_as_of_tot                 <- read.csv("data/tdat_list_tot_weekly.csv")


### Iterate through the states and calculate the MAE and WIS for the sMOA forecast.
mse_df_list                     <- list()
count_list                      <- 1
curr_state                      <- state.name[c(5,9,32,43)]# 
count <- 1
graph_data_log_95s = NULL
graph_data_log_50s = NULL
true_data = NULL
mle_start_value <- 30

for (location in c(curr_state)){
  print(paste("working on location", location))  
  
  ##### true weekly data as of a final data of reporting
  state_truth_data_file <- paste('data/state_truths/', gsub(" ", "", location), '.RData', sep = "")
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
  fcast_dates_to_match          <- fcast_dates_to_match[which(fcast_dates_to_match > as.Date("2020-02-28"))]
  
  #### subset big as of data frame to this specific location
  truth_as_of_tot_loc           <- truth_as_of_tot[truth_as_of_tot$location_name == location ,]
  
  ### hold mse 
  mse_df                        <- NULL
  ### iterate through forecast dates
  future_dates_count = 1
  for(fcast_date_idx in 1:(length(fcast_dates_to_match))){
    #truth_weekly                  <- load_truth(truth_source = "JHU",hub = "US", target_variable = "inc case", 
    #					  as_of= as.character(as.Date(fcast_dates_to_match[fcast_date_idx])), locations =location,data_location="covidData")
    
    
    #### grab the date formatted
    fcast_date                  <- fcast_dates_to_match[fcast_date_idx]
    
    actual_dates <- truth_weekly[truth_weekly$target_end_date <= (fcast_date + h*7),]$target_end_date
    
    #### subset to data that was available at forecast date
    truth_as_of                 <- truth_as_of_tot_loc[truth_as_of_tot_loc$as_of == (as.Date("2020-08-15")+(future_dates_count)*7) ,]
    #### double check its ordered properly
    truth_as_of                 <- truth_as_of[order(truth_as_of$target_end_date),]
    
    #### make sure that the target end date is less than or equal to fcast_date
    data_till_now               <- truth_as_of[truth_as_of$target_end_date <= fcast_date,]
    data_till_now$value         <- pmax(1,data_till_now$value)
    data_till_now$t             <- 1:nrow(data_till_now)
    
    if(nrow(data_till_now)<=k) next
    
    #### light smoothing and differencing and get last k 
    data_till_now_smoothed      <- gam(value~ s(t,k=round(nrow(data_till_now)/2)),data=data_till_now)$fitted.values
    actual_data_temp      <- data_till_now$value
    
    if(actual_dates[(length(actual_data_temp))]>"2020-08-15"){
      # I need at least four more weeks of data from here in the true data, then leave.
      true_data <- rbind(true_data, data.frame(date = actual_dates[(length(actual_data_temp))],
                                               value = actual_data_temp[(length(actual_data_temp))], 
                                               location = location,
                                               location_name = location)
      )
      future_dates_count = future_dates_count + 1
      if(future_dates_count == 5) break
      next
    } 
    
    #### could use gam smoother 
    to_match_in_moa             <- tail((data_till_now_smoothed),k+1)
    #### This is where the sMOA calculation actually happens! 
    #### take our matrix to match (X_diff) and efficiently compute distances from to_match_in_moa to X_diff
    gc()
    dist_to_test <- rowSums(abs(X_diff %r-% tail(diff(to_match_in_moa),k)))
    gc()
    
    #### get the closest ids in the X_diff mat
    closest_ids                 <- sort(dist_to_test,index.return = TRUE)$ix[1:closest]
    dists_of_closest            <- sort(dist_to_test)[1:closest]
    fcast                       <- head(apply(y_diff[closest_ids,],2,median),h)
    
    #### convert from difference back to raw cases
    fcast                       <- tail(data_till_now$value,1) + cumsum(fcast)
    fcast                       <- pmax(1,fcast)
    data_future                 <- tail(truth_weekly[truth_weekly$target_end_date <= (fcast_date + h*7),]$value,h) # Note that that '7' is for the length of a week (7 days).
    
    point <- pmax(1,tail(data_till_now$value,1) +cumsum(head(apply(y_diff[closest_ids,],2,median),h)))
    vars <- pmax(1,tail(data_till_now$value,1) +cumsum(head(apply(y_diff[closest_ids,],2,sd),h)))
    
    # If there are only 5 or less forecasts so far, get the error bounds from NBR models with pre-learned dispersion parameters
    # according to a bayesian optimization scheme on a hold-out set of synthetic data.
    if (length(one_step_ahead_forecasts) < 7){
      sim_nb1 <- matrix(rnbinom(5000*length(point),mu = (point+1), size = 1*dispersion_forecast),ncol=length(point),byrow = T)
      sim_nb2 <- matrix(rnbinom(5000*length(point),mu = (point+1), size = 0.5*dispersion_forecast),ncol=length(point),byrow = T)
      sim_nb3 <- matrix(rnbinom(5000*length(point),mu = (point+1), size = 0.25*dispersion_forecast),ncol=length(point),byrow = T)
      sim_nb4 <- matrix(rnbinom(5000*length(point),mu = (point+1), size = 0.125*dispersion_forecast),ncol=length(point),byrow = T)
      
      print("no optimal ks yet")
      
      # Otherwise, we learn the dispersion parameters for each forecast horizon using an MLE on the observations so far.
      # This is an online update for the MLE of the dispersion parameters.
    } else{
      mle_memory <- Inf
      
      fit_nb_function <- function(k){
        mu <- one_step_ahead_forecasts
        shifted_data <- tail(data_till_now_smoothed[2:length(data_till_now_smoothed)],length(one_step_ahead_forecasts))
        one_step_ahead_forecasts[which(one_step_ahead_forecasts<100)] = shifted_data[which(one_step_ahead_forecasts<100)]
        if(length(shifted_data)>mle_memory){
          one_step_ahead_forecasts <- one_step_ahead_forecasts[(length(shifted_data)-mle_memory+1):length(shifted_data)]
          shifted_data <- shifted_data[(length(shifted_data)-mle_memory+1):length(shifted_data)]
        }
        
        -sum(dnbinom(round(shifted_data),mu = one_step_ahead_forecasts,size=k,log = T))
      }
      optimal_k_1 <- optim(mle_start_value,fit_nb_function,method="Brent",lower = lower_mle_bound,upper = upper_mle_bound)
      
      fit_nb_function <- function(k){
        mu <- two_step_ahead_forecasts
        two_step_subset <- two_step_ahead_forecasts[1:(length(two_step_ahead_forecasts)-1)]
        shifted_data <- tail(data_till_now_smoothed[2:length(data_till_now_smoothed)],length(two_step_subset))
        if(length(shifted_data)>mle_memory){
          two_step_subset <- two_step_subset[(length(shifted_data)-mle_memory+1):length(shifted_data)]
          shifted_data <- shifted_data[(length(shifted_data)-mle_memory+1):length(shifted_data)]
        }
        two_step_subset[which(two_step_subset<100)] = shifted_data[which(two_step_subset<100)]
        -sum(dnbinom(round(shifted_data),mu = two_step_subset,size=k,log = T))
      }
      optimal_k_2 <- optim(mle_start_value,fit_nb_function,method="Brent",lower = lower_mle_bound,upper = upper_mle_bound)
      
      fit_nb_function <- function(k){
        mu <- three_step_ahead_forecasts
        three_step_subset <- two_step_ahead_forecasts[1:(length(three_step_ahead_forecasts)-2)]
        shifted_data <- tail(data_till_now_smoothed[2:length(data_till_now_smoothed)],length(three_step_subset))
        if(length(shifted_data)>mle_memory){
          three_step_subset <- three_step_subset[(length(shifted_data)-mle_memory+1):length(shifted_data)]
          shifted_data <- shifted_data[(length(shifted_data)-mle_memory+1):length(shifted_data)]
        }
        three_step_subset[which(three_step_subset<100)] = shifted_data[which(three_step_subset<100)]
        -sum(dnbinom(round(shifted_data),mu = three_step_subset,size=k,log = T))
      }
      optimal_k_3 <- optim(mle_start_value,fit_nb_function,method="Brent",lower = lower_mle_bound,upper = upper_mle_bound)
      
      fit_nb_function <- function(k){
        mu <- four_step_ahead_forecasts
        four_step_subset <- two_step_ahead_forecasts[1:(length(four_step_ahead_forecasts)-3)]
        shifted_data <- tail(data_till_now_smoothed[2:length(data_till_now_smoothed)],length(four_step_subset))
        if(length(shifted_data)>mle_memory){
          four_step_subset <- four_step_subset[(length(shifted_data)-mle_memory+1):length(shifted_data)]
          shifted_data <- shifted_data[(length(shifted_data)-mle_memory+1):length(shifted_data)]
        }
        four_step_subset[which(four_step_subset<100)] = shifted_data[which(four_step_subset<100)]
        -sum(dnbinom(round(shifted_data),mu = four_step_subset,size=k,log = T))
      }
      optimal_k_4 <- optim(mle_start_value,fit_nb_function,method="Brent",lower = lower_mle_bound,upper = upper_mle_bound)
      
      dispersion_estimates = c(optimal_k_1$par,
                               optimal_k_2$par,
                               optimal_k_3$par,
                               optimal_k_4$par)
      for(horizon_num in 2:h){
        if(dispersion_estimates[horizon_num] >= upper_mle_bound-5){
          dispersion_estimates[horizon_num] <- dispersion_estimates[horizon_num-1]*0.5
        }
      }
      
      temp_point <- point
      # for(horizon_idx in 2:h){
      #   if(temp_point[horizon_idx]<200) temp_point[horizon_idx] <- temp_point[horizon_idx-1]
      # }
      sim_nb1 <- matrix(rnbinom(5000*length(point),mu = (temp_point), size = 1*dispersion_estimates[1]*dispersion_forecast),ncol=length(point),byrow = T)
      sim_nb2 <- matrix(rnbinom(5000*length(point),mu = (temp_point), size = 0.5*dispersion_estimates[2]*dispersion_forecast),ncol=length(point),byrow = T)
      sim_nb3 <- matrix(rnbinom(5000*length(point),mu = (temp_point), size = 0.25*dispersion_estimates[3]*dispersion_forecast),ncol=length(point),byrow = T)
      sim_nb4 <- matrix(rnbinom(5000*length(point),mu = (temp_point), size = 0.125*dispersion_estimates[4]*dispersion_forecast),ncol=length(point),byrow = T)
      
      print("optimal ks are")
      # print(optimal_k_1$par)
      # print(optimal_k_2$par)
      # print(optimal_k_3$par)
      # print(optimal_k_4$par)
      print(dispersion_estimates)
      print("and the forecasts are")
      print(point)
      print(fcast)
      print("-----------")
      print(actual_dates[(length(data_till_now_smoothed))])
      print(one_step_ahead_forecasts)
      print(two_step_ahead_forecasts)
      print(three_step_ahead_forecasts)
      print(four_step_ahead_forecasts)
      print("-----------")
    }
    
    # These must be collected for MLE calculations on future dates. 
    one_step_ahead_forecasts <- c(one_step_ahead_forecasts,fcast[1])
    two_step_ahead_forecasts <- c(two_step_ahead_forecasts,fcast[2])
    three_step_ahead_forecasts <- c(three_step_ahead_forecasts,fcast[3])
    four_step_ahead_forecasts <- c(four_step_ahead_forecasts,fcast[4])
    
    # In the following, we calculate the WIS for each forecast horizon. 
    lower_95s <- c()
    upper_95s <- c()
    lower_50s <- c()
    upper_50s <- c()
    medians   <- c()
    
    lower_CI_scale <- 1
    upper_CI_scale <- 1
    
    est_intervals <- quantile(sim_nb1[,1],probs = quantiles)
    est_intervals[quantiles<0.5] <- est_intervals[quantiles<0.5] * lower_CI_scale
    est_intervals[quantiles>0.5] <- est_intervals[quantiles>0.5] * upper_CI_scale
    lower_95s <- c(lower_95s, est_intervals[2])
    upper_95s <- c(upper_95s, est_intervals[22])
    lower_50s <- c(lower_50s, est_intervals[7])
    upper_50s <- c(upper_50s, est_intervals[17])
    medians   <- c(medians, est_intervals[12])
    
    wis_tmp_1 <- weighted_interval_score(quantiles,value = est_intervals, actual_value = data_future[1])
    
    temp_est_intervals <- est_intervals
    est_intervals <- quantile(sim_nb2[,2],probs = quantiles)
    # if(any(est_intervals<=10)){
    #   est_intervals[est_intervals<=10] <- temp_est_intervals[est_intervals<=10]
    # }
    est_intervals[quantiles<0.5] <- est_intervals[quantiles<0.5] * lower_CI_scale
    est_intervals[quantiles>0.5] <- est_intervals[quantiles>0.5] * upper_CI_scale
    wis_tmp_2 <- weighted_interval_score(quantiles,value = est_intervals, actual_value = data_future[2])
    lower_95s <- c(lower_95s, est_intervals[2])
    upper_95s <- c(upper_95s, est_intervals[22])
    lower_50s <- c(lower_50s, est_intervals[7])
    upper_50s <- c(upper_50s, est_intervals[17])
    medians   <- c(medians, est_intervals[12])
    
    temp_est_intervals <- est_intervals
    est_intervals <- quantile(sim_nb3[,3],probs = quantiles)
    # if(any(est_intervals<=10)){
    #   est_intervals[est_intervals<=10] <- temp_est_intervals[est_intervals<=10]
    # }
    est_intervals[quantiles<0.5] <- est_intervals[quantiles<0.5] * lower_CI_scale
    est_intervals[quantiles>0.5] <- est_intervals[quantiles>0.5] * upper_CI_scale
    wis_tmp_3 <- weighted_interval_score(quantiles,value = est_intervals, actual_value = data_future[3])
    lower_95s <- c(lower_95s, est_intervals[2])
    upper_95s <- c(upper_95s, est_intervals[22])
    lower_50s <- c(lower_50s, est_intervals[7])
    upper_50s <- c(upper_50s, est_intervals[17])
    medians   <- c(medians, est_intervals[12])
    
    temp_est_intervals <- est_intervals
    est_intervals <- quantile(sim_nb4[,4],probs = quantiles)
    # if(any(est_intervals<=10)){
    #   est_intervals[est_intervals<=10] <- temp_est_intervals[est_intervals<=10]
    # }
    est_intervals[quantiles<0.5] <- est_intervals[quantiles<0.5] * lower_CI_scale
    est_intervals[quantiles>0.5] <- est_intervals[quantiles>0.5] * upper_CI_scale
    wis_tmp_4 <- weighted_interval_score(quantiles,value = est_intervals, actual_value = data_future[4])
    lower_95s <- c(lower_95s, est_intervals[2])
    upper_95s <- c(upper_95s, est_intervals[22])
    lower_50s <- c(lower_50s, est_intervals[7])
    upper_50s <- c(upper_50s, est_intervals[17])
    medians   <- c(medians, est_intervals[12])
    
    print("-----------")
    print(actual_dates[(length(data_till_now_smoothed))])
    print(upper_50s)
    print(lower_50s)
    print("-----------")
    
    # # Now plot all of this for Casey:
    all_data = c(data_till_now_smoothed, data_future)
    actual_dates = truth_weekly[truth_weekly$target_end_date <= (fcast_date + h*7),]$target_end_date
    graph_data = data.frame(x = actual_dates, 
                            value = all_data,
                            type = rep('true_data', times = length(all_data)),
                            ymin = all_data,
                            ymax = all_data)
    graph_data = rbind(graph_data, data.frame(x = actual_dates[(length(data_till_now_smoothed)):length(all_data)], 
                                              value = c(data_till_now_smoothed[length(data_till_now_smoothed)], point),
                                              type = rep('sMOA Prediction', times = (h+1)),
                                              ymax = c(data_till_now_smoothed[length(data_till_now_smoothed)],upper_95s),
                                              ymin = c(data_till_now_smoothed[length(data_till_now_smoothed)], lower_95s)))
    
    forecast_horizon = length(data_till_now_smoothed)
    
    
    temp_data = data.frame(x = (length(data_till_now_smoothed)):length(all_data),
                           upper95 = c(data_till_now_smoothed[length(data_till_now_smoothed)],upper_95s),
                           lower95 = c(data_till_now_smoothed[length(data_till_now_smoothed)], lower_95s),
                           value = c(data_till_now_smoothed[length(data_till_now_smoothed)], medians),
                           type = rep('error_region', times = (h+1)))
    dates_to_forecast <- as.Date(c("2020-03-21", "2020-04-18",
                           "2020-05-23", "2020-06-20", 
                           "2020-07-18","2020-08-15"))
    if((actual_dates[(length(data_till_now_smoothed))]+7) %in% dates_to_forecast){ 
      # if(actual_dates[(length(data_till_now_smoothed))] == "2020-07-18") browser()
      # Add in this data for the final graph.
      graph_data_log_95s <- rbind(graph_data_log_95s, data.frame(reference_date = rep(actual_dates[(length(data_till_now_smoothed))+1]-5, times = h),
                                                                 x = actual_dates[(length(data_till_now_smoothed)+1):length(all_data)], 
                                                                 y = c(medians),
                                                                 model = rep(paste('sMOA'), times = (h)),
                                                                 location = rep(location, times = h),
                                                                 .upper = c(upper_95s),
                                                                 .lower = c(lower_95s),
                                                                 .width = rep(0.95, times = h),
                                                                 .interval = rep('qi', times = h),
                                                                 .point = rep('median', times = h))
                                  )
      graph_data_log_50s <- rbind(graph_data_log_50s, data.frame(reference_date = rep(actual_dates[(length(data_till_now_smoothed))+1]-5, times = h),
                                                                 x = actual_dates[(length(data_till_now_smoothed)+1):length(all_data)], 
                                                                 y = c(medians),
                                                                 model = rep(paste('sMOA'), times = (h)),
                                                                 location = rep(location, times = h),
                                                                 .upper = c(upper_50s),
                                                                 .lower = c(lower_50s),
                                                                 .width = rep(0.5, times = h),
                                                                 .interval = rep('qi', times = h),
                                                                 .point = rep('median', times = h))
                                  )
      true_data <- rbind(true_data, data.frame(date = actual_dates[(length(data_till_now_smoothed))],
                                               value = actual_data_temp[(length(data_till_now_smoothed))], 
                                               location = location,
                                               location_name = location)
      )
      count = count + 1
    }
    
    wis_tot <- c(wis_tmp_1,wis_tmp_2,wis_tmp_3,wis_tmp_4)
    
    ### compute mse and append
    ### could undo horizon mean and include horizon 
    tmp_df                      <- data.frame(location=location,target_end_date=fcast_date,horizon=1:h,abs_error_moa=abs(fcast-data_future), wis_error_moa=wis_tot)
    mse_df                      <- rbind(mse_df,tmp_df)
    
    
    one_step_ahead_forecasts <- one_step_ahead_forecasts[one_step_ahead_forecasts>100]
    two_step_ahead_forecasts <- two_step_ahead_forecasts[two_step_ahead_forecasts>100]
    three_step_ahead_forecasts <- three_step_ahead_forecasts[three_step_ahead_forecasts>100]
    four_step_ahead_forecasts <- four_step_ahead_forecasts[four_step_ahead_forecasts>100]
    
  }
  
  #### create the data frame and return
  mse_df_list[[count_list]]     <- mse_df
  count_list                    <- count_list + 1
}

graph_data_log = rbind(graph_data_log_50s, graph_data_log_95s)
graph_data_log = as_tibble(graph_data_log)

true_data = as_tibble(true_data)  

# Now get the ForecastHub forecasts and put these on the tibble in the same way.
forecasts_case = readRDS(file = 'data/forecasts.rds')
# forecasts_case = subset(forecasts_case, subset = (temporal_resolution == 'wk')&(target_variable == 'inc case'))
forecast_dates <- dates_to_forecast - 5

attach(forecasts_case)
hub_forecasts <- NULL
for(temp_date in forecast_dates){
  count = 1
  for(state in curr_state){
    for(horizon in 1:4){
      
      for(model_name in unique(forecasts_case$model)){
        temp_df = forecasts_case[(target_end_date <= (as.Date(temp_date)+horizon*7))&
                                   (target_end_date >= (as.Date(temp_date)+(horizon-1)*7))&
                                   (forecast_date == (as.Date(temp_date)))&
                                   (location_name == state)&
                                   (model == model_name),]
        if(nrow(temp_df) == 0) next
        temp_df$quantile[is.na(temp_df$quantile)]<- 0
        temp_new_row = data.frame(model = model_name,
                                  location = state,
                                  reference_date = as.Date(temp_date),
                                  y = temp_df[temp_df$type == 'point',]$value,
                                  x = temp_df$target_end_date[1],
                                  .lower = temp_df[temp_df$quantile == 0.025,]$value,
                                  .upper = temp_df[temp_df$quantile == 0.975,]$value,
                                  .width = 0.95,
                                  .interval = 'qi',
                                  .point = 'median') 
        hub_forecasts = rbind(hub_forecasts, temp_new_row)
        temp_new_row = data.frame(model = model_name,
                                  location = state,
                                  reference_date = as.Date(temp_date),
                                  y = temp_df[temp_df$type == 'point',]$value,
                                  x = temp_df$target_end_date[1],
                                  .lower = temp_df[temp_df$quantile == 0.25,]$value,
                                  .upper = temp_df[temp_df$quantile == 0.75,]$value,
                                  .width = 0.5,
                                  .interval = 'qi',
                                  .point = 'median') 
        hub_forecasts = rbind(hub_forecasts, temp_new_row)
      }
    }
  }
}
hub_forecasts = as_tibble(hub_forecasts)
detach(forecasts_case)

full_data = rbind(graph_data_log, hub_forecasts)
full_data$location_name = full_data$location

ggplot() +
  geom_lineribbon(
    mapping = aes(x = x, y = y, ymin = .lower, ymax = .upper, group = reference_date),
    data = full_data 
  ) +
  geom_line(
    mapping = aes(x = date, y = value),
    linewidth = 0.75,
    color = "orange",
    data = true_data)  +
  scale_linetype("Reported Data") +
  scale_fill_brewer("Interval Level", labels = c("95%", "50%")) +
  scale_x_date("Forecast target date", date_breaks = "2 month", date_labels = "%b %Y",
               # limits = c(as.Date("2020-07-25"), as.Date("2021-11-15")),
               expand = expansion()) +
  # scale_y_continuous("Cases", labels = comma) +
  facet_grid(rows = vars(location_name), cols = vars(model), scales = "free_y") +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 20, hjust=1, vjust=1)
  )+ 
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 15),
        title = element_text(size=15),strip.text.x = element_text(size = 15),strip.text.y = element_text(size = 15)) +
  theme(panel.spacing = unit(0, "cm"),
        plot.margin = margin(0, 0, 0, 0, "cm"), 
        plot.caption = element_blank()) + theme(legend.text=element_text(size=15))




