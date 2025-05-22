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
library(this.path)
setwd(paste0(this.path::here(),"/../"))
source("R/smoa_helpers.R")
ncores <- 51
lopez_models <- c("BPagano-RtDriven", 'CEID-Walk', 'CovidAnalytics-DELPHI', 'COVIDhub-baseline', 'COVIDhub-4_week_ensemble',
                  'COVIDhub-trained_ensemble', 'Covid19Sim-Simulator', 'CU-select', 'FAIR-NRAR', 'FRBSF_Wilson-Econometric',
                  'IEM_MED-CovidProject', 'IowaStateLW-STEM', 'JHUAPL-Bucky', 'JHU_CSSE-DECOM', 'JHU_IDD-CovidSP',
                  'Karlen-pypm', 'LNQ-ens1', 'LANL-GrowthRate', 'Microsoft-DeepSTIA', 'MOBS-GLEAM_COVID',
                  'RobertWalraven-ESG', 'UCLA-SuEIR', 'USC-SI_kJalpha', 'UMass-MechBayes', 'UVA-Ensemble')

sim_idx                       <- 1
curr_state                      <- state.name[sim_idx]

dispersion_forecast_scaling <- seq(from=0.1, to=5, length.out = 64)
dispersion_forecast <- dispersion_forecast_scaling[sim_idx]
h                               <- 4

# These are all from the bayesian optimization.
num_curves <- 18387
k <- 5
closest <- 4422
lower_CI_scale <- 1
upper_CI_scale <- 1
dispersion_forecast <- 1
mle_lower_bound <- 1

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

quantiles <- c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)
## read in stored truth data with as_of

# dates_to_forecast <- as.Date(c("2021-05-22", "2021-07-03", "2021-08-14", "2021-09-25", "2021-11-13", "2021-12-25"))
dates_to_forecast <- as.Date(c("2021-11-20", "2021-12-25", "2022-01-22", "2022-03-12", "2022-04-23", "2022-06-04"))-14

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
mle_start_value <- 0.5

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
  targetend_dates_to_match          <- unique(truth_weekly$target_end_date)
  
  # We need to collect the one-step-ahead forecasts for all available data (potentially not just
  # from august 15th onwards).
  # if(length(which(targetend_dates_to_match <= as.Date("2020-02-28")))>=6){
  #   forecasts_data = get_early_pandemic_errors(targetend_dates_to_match, location, k, 
  #                                              truth_as_of_tot, X_diff, y_diff,
  #                                              truth_weekly, current_end_date = "2020-02-28")
  # } else {
  #   forecasts_data = NULL
  # }
  if(length(which(targetend_dates_to_match <= as.Date("2020-08-15")))>=5){
    forecasts_data = get_early_pandemic_errors(targetend_dates_to_match, location, k, 
                                               truth_as_of_tot, X_diff, y_diff,
                                               truth_weekly)
  } else {
    forecasts_data = NULL
  }
  
  #### subset to "2020-08-15" when the revised API started up
  targetend_dates_to_match          <- targetend_dates_to_match[which((targetend_dates_to_match > as.Date("2020-08-15")))]
  
  #### subset big as of data frame to this specific location
  truth_as_of_tot_loc           <- truth_as_of_tot[truth_as_of_tot$location_name == location ,]
  
  ### hold mse 
  mse_df                        <- NULL
  ### iterate through forecast dates
  future_dates_count = 0
  for(last_targetend_date_idx in 1:(length(targetend_dates_to_match))){
    paste0("calculating forecast instance ", last_targetend_date_idx, " of ", length(targetend_dates_to_match))
    #truth_weekly                  <- load_truth(truth_source = "JHU",hub = "US", target_variable = "inc case", 
    #					  as_of= as.character(as.Date(targetend_dates_to_match[last_targetend_date_idx])), locations =location,data_location="covidData")
    
    #### grab the date formatted
    curr_targetend_date                  <- targetend_dates_to_match[last_targetend_date_idx]
    
    actual_dates <- truth_weekly[truth_weekly$target_end_date <= (curr_targetend_date + h*7),]$target_end_date
    
    #### subset to data that was available at forecast date
    truth_as_of                 <- truth_as_of_tot_loc[truth_as_of_tot_loc$as_of == curr_targetend_date ,]
    #### double check its ordered properly
    truth_as_of                 <- truth_as_of[order(truth_as_of$target_end_date),]
    #### make sure that the target end date is less than or equal to curr_targetend_date
    data_till_now               <- truth_as_of[truth_as_of$target_end_date <= curr_targetend_date,]
    data_till_now$value         <- pmax(1,data_till_now$value)
    data_till_now$t             <- 1:nrow(data_till_now)
    
    if(nrow(data_till_now)<=k) next
    
    #### light smoothing and differencing and get last k 
    data_till_now_smoothed      <- pmax(1,gam(value~ s(t,k=round(nrow(data_till_now)/2)),data=data_till_now)$fitted.values)
    actual_data_temp            <- data_till_now$value
    
    if(curr_targetend_date>="2022-06-04"){
      # I need at least four more weeks of data from here in the true data, then leave.
      true_data <- rbind(true_data, data.frame(date = actual_dates[(length(actual_data_temp))],
                                               value = actual_data_temp[(length(actual_data_temp))], 
                                               location = location,
                                               location_name = location)
      )
      future_dates_count = future_dates_count + 1
      if(future_dates_count == 4) break
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
    data_future                 <- tail(truth_weekly[truth_weekly$target_end_date <= (curr_targetend_date + h*7),]$value,h) # Note that that '7' is for the length of a week (7 days).
    
    ##
    # data_till_now is the vector we used to get the bit that we 'match' in sMOA.  This is where I grab the most recent value.
    # This is used if we decide to impose 'guard rails' on the sMOA forecast so it cannot explode upwards to high.
    most_recent_value = data_till_now$value[length(data_till_now$value)]
    ##
    
    mle_memory <- Inf
    
    point <- pmax(1,tail(data_till_now$value,1) +cumsum(head(apply(y_diff[closest_ids,],2,median),h)))
    vars <- pmax(1,tail(data_till_now$value,1) +cumsum(head(apply(y_diff[closest_ids,],2,sd),h)))
    
    # # If there are only 5 or less forecasts so far, get the error bounds from NBR models with pre-learned dispersion parameters
    # # according to a bayesian optimization scheme on a hold-out set of synthetic data.
    # if (length(one_step_ahead_forecasts) < 7){
    #   sim_nb1 <- matrix(rnbinom(5000*length(point),mu = (point+1), size = 1*dispersion_forecast),ncol=length(point),byrow = T)
    #   sim_nb2 <- matrix(rnbinom(5000*length(point),mu = (point+1), size = 0.5*dispersion_forecast),ncol=length(point),byrow = T)
    #   sim_nb3 <- matrix(rnbinom(5000*length(point),mu = (point+1), size = 0.25*dispersion_forecast),ncol=length(point),byrow = T)
    #   sim_nb4 <- matrix(rnbinom(5000*length(point),mu = (point+1), size = 0.125*dispersion_forecast),ncol=length(point),byrow = T)
    #   
    #   print("no optimal ks yet")
    #   
    #   # Otherwise, we learn the dispersion parameters for each forecast horizon using an MLE on the observations so far.
    #   # This is an online update for the MLE of the dispersion parameters.
    # } 
    if (is.null(nrow(forecasts_data))){
      sim_nb1 <- matrix(rnbinom(5000*length(point),mu = (point+1), size = 20*dispersion_forecast),ncol=length(point),byrow = T)
      sim_nb2 <- matrix(rnbinom(5000*length(point),mu = (point+1), size = 20*0.5*dispersion_forecast),ncol=length(point),byrow = T)
      sim_nb3 <- matrix(rnbinom(5000*length(point),mu = (point+1), size = 20*0.3*dispersion_forecast),ncol=length(point),byrow = T)
      sim_nb4 <- matrix(rnbinom(5000*length(point),mu = (point+1), size = 20*0.3*dispersion_forecast),ncol=length(point),byrow = T)
      
      print(curr_targetend_date)
      # Otherwise, we learn the dispersion parameters for each forecast horizon using an MLE on the observations so far.
      # This is an online update for the MLE of the dispersion parameters.
    } else if ((nrow(forecasts_data) < 6*4)){
      sim_nb1 <- matrix(rnbinom(5000*length(point),mu = (point+1), size = 20*dispersion_forecast),ncol=length(point),byrow = T)
      sim_nb2 <- matrix(rnbinom(5000*length(point),mu = (point+1), size = 20*0.5*dispersion_forecast),ncol=length(point),byrow = T)
      sim_nb3 <- matrix(rnbinom(5000*length(point),mu = (point+1), size = 20*0.3*dispersion_forecast),ncol=length(point),byrow = T)
      sim_nb4 <- matrix(rnbinom(5000*length(point),mu = (point+1), size = 20*0.3*dispersion_forecast),ncol=length(point),byrow = T)
      
      print(curr_targetend_date)
      # Otherwise, we learn the dispersion parameters for each forecast horizon using an MLE on the observations so far.
      # This is an online update for the MLE of the dispersion parameters.
    } else {
      print("WE ARE ACTUALLY DOING THE ONLINE STUFF")
      mle_memory <- Inf
      
      fit_nb_function <- function(k_){
        subset_forecasts = forecasts_data %>% subset(horizon == 1)
        if(nrow(subset_forecasts)>mle_memory){
          subset_forecasts = subset_forecasts[(nrow(subset_forecasts)-mle_memory):nrow(subset_forecasts),]
        } 
        mu <- subset_forecasts$true_values
        horizon_subset <- subset_forecasts$forecasts
        horizon_subset[which(horizon_subset<100)] = mu[which(horizon_subset<100)]
        xx=-sum(dnbinom(round(mu),mu = horizon_subset,size=1/k_,log = T))
        return(xx)
      }
      optimal_k_1 <- optim(mle_start_value,fit_nb_function,method="Brent",lower = 0, upper = 1)
      
      fit_nb_function <- function(k_){
        subset_forecasts = forecasts_data %>% subset(horizon == 2)
        if(nrow(subset_forecasts)>mle_memory){
          subset_forecasts = subset_forecasts[(nrow(subset_forecasts)-mle_memory):nrow(subset_forecasts),]
        } 
        mu <- subset_forecasts$true_values
        horizon_subset <- subset_forecasts$forecasts
        horizon_subset[which(horizon_subset<100)] = mu[which(horizon_subset<100)]
        xx=-sum(dnbinom(round(mu),mu = horizon_subset,size=1/k_,log = T))
        return(xx)
      }
      optimal_k_2 <- optim(mle_start_value,fit_nb_function,method="Brent",lower = 0, upper = 1)
      
      fit_nb_function <- function(k_){
        subset_forecasts = forecasts_data %>% subset(horizon == 3)
        if(nrow(subset_forecasts)>mle_memory){
          subset_forecasts = subset_forecasts[(nrow(subset_forecasts)-mle_memory):nrow(subset_forecasts),]
        } 
        mu <- subset_forecasts$true_values
        horizon_subset <- subset_forecasts$forecasts
        horizon_subset[which(horizon_subset<100)] = mu[which(horizon_subset<100)]
        xx=-sum(dnbinom(round(mu),mu = horizon_subset,size=1/k_,log = T))
        return(xx)
      }
      optimal_k_3 <- optim(mle_start_value,fit_nb_function,method="Brent",lower = 0, upper = 1)
      
      fit_nb_function <- function(k_){
        subset_forecasts = forecasts_data %>% subset(horizon == 4)
        if(nrow(subset_forecasts)>mle_memory){
          subset_forecasts = subset_forecasts[(nrow(subset_forecasts)-mle_memory):nrow(subset_forecasts),]
        } 
        mu <- subset_forecasts$true_values
        horizon_subset <- subset_forecasts$forecasts
        horizon_subset[which(horizon_subset<100)] = mu[which(horizon_subset<100)]
        xx=-sum(dnbinom(round(mu),mu = horizon_subset,size=1/k_,log = T))
        return(xx)
      }
      optimal_k_4 <- optim(mle_start_value,fit_nb_function,method="Brent",lower = 0, upper = 1)
      
      
      dispersion_estimates = c(1/optimal_k_1$par,
                               1/optimal_k_2$par,
                               1/optimal_k_3$par,
                               1/optimal_k_4$par)
      
      temp_point <- point
      sim_nb1 <- matrix(rnbinom(5000*length(point),mu = (temp_point), size = dispersion_estimates[1]*dispersion_forecast),ncol=length(point),byrow = T)
      sim_nb2 <- matrix(rnbinom(5000*length(point),mu = (temp_point), size = dispersion_estimates[2]*dispersion_forecast),ncol=length(point),byrow = T)
      sim_nb3 <- matrix(rnbinom(5000*length(point),mu = (temp_point), size = dispersion_estimates[3]*dispersion_forecast),ncol=length(point),byrow = T)
      sim_nb4 <- matrix(rnbinom(5000*length(point),mu = (temp_point), size = dispersion_estimates[4]*dispersion_forecast),ncol=length(point),byrow = T)
    }
    
    # In the following, we calculate the WIS for each forecast horizon. 
    lower_95s <- c()
    upper_95s <- c()
    lower_50s <- c()
    upper_50s <- c()
    medians   <- c()
    
    lower_CI_scale <- 1
    upper_CI_scale <- 1
    
    ## Horizon == 1/
    curr_horizon <- 1
    est_intervals <- quantile(sim_nb1[,1],probs = quantiles)
    est_intervals <- point[curr_horizon] + (est_intervals-est_intervals[quantiles==0.5]) 
    lower_95s <- c(lower_95s, est_intervals[1])
    upper_95s <- c(upper_95s, est_intervals[7])
    lower_50s <- c(lower_50s, est_intervals[3])
    upper_50s <- c(upper_50s, est_intervals[5])
    medians   <- c(medians, est_intervals[4])
    wis_tmp_1 <- weighted_interval_score(quantiles,value = est_intervals, actual_value = data_future[1])
    temp_est_intervals <- est_intervals
    
    ## Horizon == 2
    curr_horizon <- 2
    est_intervals <- quantile(sim_nb2[,2],probs = quantiles)
    est_intervals <- point[curr_horizon] + (est_intervals-est_intervals[quantiles==0.5]) 
    lower_95s <- c(lower_95s, est_intervals[1])
    upper_95s <- c(upper_95s, est_intervals[7])
    lower_50s <- c(lower_50s, est_intervals[3])
    upper_50s <- c(upper_50s, est_intervals[5])
    medians   <- c(medians, est_intervals[4])
    wis_tmp_2 <- weighted_interval_score(quantiles,value = est_intervals, actual_value = data_future[2])
    temp_est_intervals <- est_intervals
    
    ## Horizon == 3
    curr_horizon <- 3
    est_intervals <- quantile(sim_nb3[,3],probs = quantiles)
    est_intervals <- point[curr_horizon] + (est_intervals-est_intervals[quantiles==0.5]) 
    lower_95s <- c(lower_95s, est_intervals[1])
    upper_95s <- c(upper_95s, est_intervals[7])
    lower_50s <- c(lower_50s, est_intervals[3])
    upper_50s <- c(upper_50s, est_intervals[5])
    medians   <- c(medians, est_intervals[4])
    wis_tmp_3 <- weighted_interval_score(quantiles,value = est_intervals, actual_value = data_future[3])
    temp_est_intervals <- est_intervals
    
    ## Horizon == 4
    curr_horizon <- 4
    est_intervals <- quantile(sim_nb4[,4],probs = quantiles)
    est_intervals <- point[curr_horizon] + (est_intervals-est_intervals[quantiles==0.5]) 
    lower_95s <- c(lower_95s, est_intervals[1])
    upper_95s <- c(upper_95s, est_intervals[7])
    lower_50s <- c(lower_50s, est_intervals[3])
    upper_50s <- c(upper_50s, est_intervals[5])
    medians   <- c(medians, est_intervals[4])
    wis_tmp_4 <- weighted_interval_score(quantiles,value = est_intervals, actual_value = data_future[4])
    temp_est_intervals <- est_intervals
    
    print("-----------")
    print(actual_dates[(length(data_till_now_smoothed))])
    print(upper_50s)
    print(lower_50s)
    print("-----------")
    
    # # Now plot all of this for Casey:
    all_data = c(data_till_now_smoothed, data_future)
    actual_dates = truth_weekly[truth_weekly$target_end_date <= (curr_targetend_date + h*7),]$target_end_date
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
    # dates_to_forecast <- as.Date(c("2020-03-21", "2020-04-18",
    #                        "2020-05-23", "2020-06-20", 
    #                        "2020-07-18","2020-08-15"))
    # dates_to_forecast <- as.Date(c("2020-09-04", "2020-10-17",
    #                                "2020-11-21", "2021-01-02", 
    #                                "2021-02-13","2021-03-20"))
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
    tmp_df                      <- data.frame(location=location,fcast_date=curr_targetend_date,forecasts=fcast, true_values = data_future,
                                              horizon=1:h,abs_error_moa=abs(fcast-data_future), 
                                              wis_error_moa=wis_tot,
                                              target_end_date = c(curr_targetend_date + 7,
                                                                  curr_targetend_date + 14,
                                                                  curr_targetend_date + 21,
                                                                  curr_targetend_date + 28),
                                              lower_95s = lower_95s,
                                              upper_95s = upper_95s
    )
    mse_df                      <- rbind(mse_df,tmp_df)
    
    forecasts_data = rbind(forecasts_data, data.frame(forecasts = fcast,
                                                      horizon = c(1:4),
                                                      target_end_date = c(curr_targetend_date + 7,
                                                                          curr_targetend_date + 14,
                                                                          curr_targetend_date + 21,
                                                                          curr_targetend_date + 28), 
                                                      true_values = data_future
                          )
    )
    
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
# dates_to_forecast <- as.Date(c("2020-03-21", "2020-04-18",
#                                "2020-05-23", "2020-06-20", 
#                                "2020-07-18","2020-08-15"))
# dates_to_forecast <- as.Date(c("2020-09-04", "2020-10-17",
#                                "2020-11-21", "2021-01-02", 
#                                "2021-02-13","2021-03-20"))
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

full_data$model = as.factor(full_data$model)
levels(full_data$model) = c('4-week ensemble', 'baseline', 'trained ensemble', 'sMOA')

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
  facet_grid(cols = vars(location_name), rows = vars(model)) + #, scales = "free_y"
  theme_bw() +
  ylab('Weekly Incident COVID Cases') +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(angle = 20, hjust=1, vjust=1)
  )+ 
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 15),
        title = element_text(size=15),strip.text.x = element_text(size = 15),strip.text.y = element_text(size = 15)) +
  theme(panel.spacing = unit(0, "cm"),
        plot.margin = margin(0, 0, 0, 0, "cm"), 
        plot.caption = element_blank()) + theme(legend.text=element_text(size=15))




