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
library(this.path)
setwd(paste0(this.path::here(), '/../'))
source("R/smoa_helpers.R")

ncores <- 51

###########
# The following is only for when we run the grid search to test several different
# hyperparameters:
sim_idx                         <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(paste("Running slurm job:", sim_idx))
if(is.na(sim_idx)){
  sim_idx                       <- 3
}
curr_state                      <- state.name[sim_idx]


num_of_each_curve <- c(20000, 30000, 40000, 50000)
ks <- c(4, 6,7,8)
epsilons <- c(5e-4, 5e-3, 5e-5, 5e-2)
closest_ids_vec <- c(1000, 2000, 3000, 5000)

dispersion_forecast_scaling <- seq(from=1, to=40, length.out = 8)
dispersion_forecast <- dispersion_forecast_scaling[sim_idx%%8+1]
mle_lower_bound_scaling <- seq(from=1, to=0.001, length.out = 8)
mle_lower_bound <- mle_lower_bound_scaling[floor((sim_idx-1)/8)+1]

h                               <- 4
optim_method <- 'Brent'

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
name_of_change <- paste("k", k, "num_curves", num_curves, "closest", closest, "dispersion", 
                        round(dispersion_forecast*10000), 'mlebound', round(mle_lower_bound*10000), sep = "_")
state_log_directory <- paste(name_of_change, "_state_records", sep = "")

############################################
### Generate Synthetic List for Training ###
############################################
types_of_curves <- c("sir_rollercoaster", "sir_rollercoaster_wiggle", 'seasonal')
N <- num_curves*length(types_of_curves)

sockettype <- "PSOCK"

if(!file.exists(paste("data/synthetic_logs/synthetic_simidx_", sim_idx, "_num_curves_", num_curves, ".RData", sep = ""))){
  cl <- parallel::makeCluster(spec = ncores,type = sockettype)
  setDefaultCluster(cl)
  registerDoParallel(cl)
  sim_ts <- foreach(i=1:N, #added extra 300 to compensate for extra seasonality sims
                    .errorhandling = "pass",
                    .verbose = F)%dopar%{
                      #curve_type <- types_of_curves[rep(c(1:length(types_of_curves)), each = num_curves)][i]
                      library(data.table)
                      library(LearnBayes)
                      library(LaplacesDemon)
                      source("R/smoa_helpers.R")
                      
                      curve_type <- types_of_curves[rep(c(1:length(types_of_curves)), each = num_curves)][i]
                      templist <- gen_curve(curve_type)
                      templist
                    }
  stopCluster(cl)
  A = lapply(sim_ts,function(x){ length(x)})
  unique(A)
  sim_ts = sim_ts[A > 2]
  save(sim_ts, file = paste("data/synthetic_logs/synthetic_simidx_", sim_idx, "_num_curves_", num_curves, ".RData", sep = ""))
} else {
  load(paste("data/synthetic_logs/synthetic_simidx_", sim_idx, "_num_curves_", num_curves, ".RData", sep = ""))
}
print("Creating the embedding matrix.")

#### read in embedding mat design matrix and outcome matrix
embed_mat                       <- create_embed_matrix(sim_ts,h=(h+1),k=(k+1))
X                               <- embed_mat[[1]]
y                               <- embed_mat[[2]]

#save(X, file = 'synthetic_X.RData')
#save(y, file = 'synthetic_y.RData')
# load(file = 'synthetic_X.RData')
# load(file = 'synthetic_y.RData')

### convert to differences
X_diff                          <- t(apply(X,1,diff))
y_diff                          <- t(apply(y,1,diff))

#save(X_diff, file = 'synthetic_X_diff.RData')
#save(y_diff, file = 'synthetic_y_diff.RData')
# load(file = 'synthetic_X_diff.RData')
# load(file = 'synthetic_y_diff.RData')

parfctn = function(x){
  set.seed(13)
  library(covidHubUtils)
  library(mgcv)
  library(collapse)
  library(dplyr)
  #quantiles <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
  quantiles <- c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975)
  ### read in stored truth data with as_of
  
  # Coverage Data:
  coverage_data <- NULL
  dispersion_logs <- NULL
  
  # We pre-built this data file to cut on api calls to github.
  truth_as_of_tot                 <- read.csv("data/tdat_list_tot_weekly.csv")
  
  ### Iterate through the states and calculate the MAE and WIS for the sMOA forecast.
  mse_df_list                     <- list()
  count_list                      <- 1
  curr_state <- state.name[x]
  mle_start_value <- 0.5
  
  for (location in c(curr_state)){
    #print(paste("working on location", location))  
    
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
    if(length(which(targetend_dates_to_match <= as.Date("2020-08-15")))>=5){
      forecasts_data = get_early_pandemic_errors(targetend_dates_to_match, location, k, 
                                                 truth_as_of_tot, X_diff, y_diff,
                                                 truth_weekly)
    } else {
      forecasts_data = NULL
    }
    
    #### subset to "2020-08-15" when the revised API started up
    targetend_dates_to_match          <- targetend_dates_to_match[which(targetend_dates_to_match >=  as.Date("2020-08-15"))]
    
    #### subset big as of data frame to this specific location
    truth_as_of_tot_loc           <- truth_as_of_tot[truth_as_of_tot$location_name == location ,]
    
    ### hold mse 
    mse_df                        <- NULL 
    
    state_X = NULL
    state_y = NULL

    write.csv(forecasts_data, file = paste("data/early_pandemic_forecasts/",location, ".csv",sep = ""))
    
    ### iterate through forecast dates
    for ( last_targetend_date_idx in 1:(length(targetend_dates_to_match)) ){ 
      paste0("calculating forecast instance ", last_targetend_date_idx, " of ", length(targetend_dates_to_match))
      #truth_weekly                  <- load_truth(truth_source = "JHU",hub = "US", target_variable = "inc case", 
      #					  as_of= as.character(as.Date(targetend_dates_to_match[last_targetend_date_idx])), locations =location,data_location="covidData")
      
      #### grab the date formatted
      curr_targetend_date                  <- targetend_dates_to_match[last_targetend_date_idx]
      
      #### subset to data that was available at forecast date
      truth_as_of                 <- truth_as_of_tot_loc[truth_as_of_tot_loc$as_of == curr_targetend_date ,]
      #### double check its ordered properly
      truth_as_of                 <- truth_as_of[order(truth_as_of$target_end_date),]
      #### make sure that the target end date is less than or equal to curr_targetend_date
      data_till_now               <- truth_as_of[truth_as_of$target_end_date <= curr_targetend_date,]
      data_till_now$value         <- pmax(1,data_till_now$value)
      data_till_now$t             <- 1:nrow(data_till_now)
      
      #### light smoothing and differencing and get last k 
      data_till_now_smoothed      <- pmax(1,gam(value~ s(t,k=round(nrow(data_till_now)/2)),data=data_till_now)$fitted.values)
      #data_till_now_smoothed      <- data_till_now$value 
      
      #### could use gam smoother 
      to_match_in_moa             <- tail((data_till_now_smoothed),k+1)
      
      ##
      # The only thing that's happening here is that I'm collecting all the segments on which we forecast for this
      # state.  This is for an analysis afterwards and has nothing to do with the model fit here. 
      state_X = rbind(state_X, to_match_in_moa)
      state_y = rbind(state_y, tail(truth_weekly[truth_weekly$target_end_date <= (curr_targetend_date + h*7),]$value,h))
      ##
      
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
      
      point <- pmax(1,tail(data_till_now$value,1) +cumsum(head(apply(y_diff[closest_ids,],2,median),h)))
      vars <- pmax(1,tail(data_till_now$value,1) +cumsum(head(apply(y_diff[closest_ids,],2,sd),h)))
      
      mle_memory <- Inf
      
      # If there are only 5 or less forecasts so far, get the error bounds from NBR models with pre-learned dispersion parameters
      # according to a bayesian optimization scheme on a hold-out set of synthetic data.
      if ((nrow(forecasts_data) < 6*4)){
        sim_nb1 <- matrix(rnbinom(5000*length(point),mu = (point+1), size = 20*dispersion_forecast),ncol=length(point),byrow = T)
        sim_nb2 <- matrix(rnbinom(5000*length(point),mu = (point+1), size = 20*0.5*dispersion_forecast),ncol=length(point),byrow = T)
        sim_nb3 <- matrix(rnbinom(5000*length(point),mu = (point+1), size = 20*0.3*dispersion_forecast),ncol=length(point),byrow = T)
        sim_nb4 <- matrix(rnbinom(5000*length(point),mu = (point+1), size = 0.3*dispersion_forecast),ncol=length(point),byrow = T)
        
        print(curr_targetend_date)
        # Otherwise, we learn the dispersion parameters for each forecast horizon using an MLE on the observations so far.
        # This is an online update for the MLE of the dispersion parameters.
      } else{
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
      
      lower_CI_scale <- 1
      upper_CI_scale <- 1
      
      ## Horizon == 1/
      curr_horizon <- 1
      est_intervals <- quantile(sim_nb1[,1],probs = quantiles)
      est_intervals <- point[curr_horizon] + (est_intervals-est_intervals[quantiles==0.5]) 
      lower_95s <- c(lower_95s, est_intervals[1])
      upper_95s <- c(upper_95s, est_intervals[7])
      wis_tmp_1 <- weighted_interval_score(quantiles,value = est_intervals, actual_value = data_future[1])
      overprediction_tmp_1 <- overprediction(quantiles,value = est_intervals, actual_value = data_future[1])
      underprediction_tmp_1 <- underprediction(quantiles,value = est_intervals, actual_value = data_future[1])
      sharpness_tmp_1 <- sharpness(quantiles,value = est_intervals, actual_value = data_future[1])
      
      # Gather coverage data for this horizon:
      for(quantile_idx in 0:2){
        lower_quan <- est_intervals[1+quantile_idx]
        upper_quan <- est_intervals[length(quantiles)-quantile_idx]
        coverage <- as.numeric((data_future[curr_horizon]>lower_quan)&
                                 (data_future[curr_horizon]<=upper_quan))
        
        temp_row <- data.frame(State = location, Date = curr_targetend_date, 
                               horizon = curr_horizon, 
                               Lower_Quantile = lower_quan, 
                               Upper_Quantile = upper_quan, 
                               Coverage = coverage)
        coverage_data <- rbind(coverage_data, temp_row)
      }
      
      
      ## Horizon == 2
      curr_horizon <- 2
      temp_est_intervals <- est_intervals
      est_intervals <- quantile(sim_nb2[,2],probs = quantiles)
      est_intervals <- point[curr_horizon] + (est_intervals-est_intervals[quantiles==0.5])
      
      wis_tmp_2 <- weighted_interval_score(quantiles,value = est_intervals, actual_value = data_future[2])
      overprediction_tmp_2 <- overprediction(quantiles,value = est_intervals, actual_value = data_future[2])
      underprediction_tmp_2 <- underprediction(quantiles,value = est_intervals, actual_value = data_future[2])
      sharpness_tmp_2 <- sharpness(quantiles,value = est_intervals, actual_value = data_future[2])
      lower_95s <- c(lower_95s, est_intervals[1])
      upper_95s <- c(upper_95s, est_intervals[7])
      
      # Gather coverage data for this horizon:
      for(quantile_idx in 0:2){
        lower_quan <- est_intervals[1+quantile_idx]
        upper_quan <- est_intervals[length(quantiles)-quantile_idx]
        coverage <- as.numeric((data_future[curr_horizon]>lower_quan)&
                                 (data_future[curr_horizon]<=upper_quan))
        
        temp_row <- data.frame(State = location, Date = curr_targetend_date, 
                               horizon = curr_horizon, 
                               Lower_Quantile = lower_quan, 
                               Upper_Quantile = upper_quan, 
                               Coverage = coverage)
        coverage_data <- rbind(coverage_data, temp_row)
      }
      
      
      ## Horizon == 3
      curr_horizon <- 3
      temp_est_intervals <- est_intervals
      est_intervals <- quantile(sim_nb3[,3],probs = quantiles)
      est_intervals <- point[curr_horizon] + (est_intervals-est_intervals[quantiles==0.5])
      wis_tmp_3 <- weighted_interval_score(quantiles,value = est_intervals, actual_value = data_future[3])
      overprediction_tmp_3 <- overprediction(quantiles,value = est_intervals, actual_value = data_future[3])
      underprediction_tmp_3 <- underprediction(quantiles,value = est_intervals, actual_value = data_future[3])
      sharpness_tmp_3 <- sharpness(quantiles,value = est_intervals, actual_value = data_future[3])
      
      lower_95s <- c(lower_95s, est_intervals[1])
      upper_95s <- c(upper_95s, est_intervals[7])
      
      # Gather coverage data for this horizon:
      for(quantile_idx in 0:2){
        lower_quan <- est_intervals[1+quantile_idx]
        upper_quan <- est_intervals[length(quantiles)-quantile_idx]
        coverage <- as.numeric((data_future[curr_horizon]>lower_quan)&
                                 (data_future[curr_horizon]<=upper_quan))
        
        temp_row <- data.frame(State = location, Date = curr_targetend_date, 
                               horizon = curr_horizon, 
                               Lower_Quantile = lower_quan, 
                               Upper_Quantile = upper_quan, 
                               Coverage = coverage)
        coverage_data <- rbind(coverage_data, temp_row)
      }
      
      
      ## Horizon == 4
      curr_horizon <- 4
      temp_est_intervals <- est_intervals
      est_intervals <- quantile(sim_nb4[,4],probs = quantiles)
      est_intervals <- point[curr_horizon] + (est_intervals-est_intervals[quantiles==0.5])
      wis_tmp_4 <- weighted_interval_score(quantiles,value = est_intervals, actual_value = data_future[4])
      overprediction_tmp_4 <- overprediction(quantiles,value = est_intervals, actual_value = data_future[4])
      underprediction_tmp_4 <- underprediction(quantiles,value = est_intervals, actual_value = data_future[4])
      sharpness_tmp_4 <- sharpness(quantiles,value = est_intervals, actual_value = data_future[4])
      
      lower_95s <- c(lower_95s, est_intervals[1])
      upper_95s <- c(upper_95s, est_intervals[7])
      
      # Gather coverage data for this horizon:
      for(quantile_idx in 0:2){
        lower_quan <- est_intervals[1+quantile_idx]
        upper_quan <- est_intervals[length(quantiles)-quantile_idx]
        coverage <- as.numeric((data_future[curr_horizon]>lower_quan)&
                                 (data_future[curr_horizon]<=upper_quan))
        
        temp_row <- data.frame(State = location, Date = curr_targetend_date, 
                               horizon = curr_horizon, 
                               Lower_Quantile = lower_quan, 
                               Upper_Quantile = upper_quan, 
                               Coverage = coverage)
        coverage_data <- rbind(coverage_data, temp_row)
      }
      
      
      # # Now plot all of this for Casey:
      all_data <- c(data_till_now_smoothed, data_future)
      
      wis_tot <- c(wis_tmp_1,wis_tmp_2,wis_tmp_3,wis_tmp_4)
      sharpness_tot <- c(sharpness_tmp_1,sharpness_tmp_2,sharpness_tmp_3,sharpness_tmp_4)
      overprediction_tot <- c(overprediction_tmp_1,overprediction_tmp_2,overprediction_tmp_3,overprediction_tmp_4)
      underprediction_tot <- c(underprediction_tmp_1,underprediction_tmp_2,underprediction_tmp_3,underprediction_tmp_4)
      
      ### compute mse and append
      ### could undo horizon mean and include horizon 
      tmp_df                      <- data.frame(location=location,fcast_date=curr_targetend_date,forecasts=fcast, true_values = data_future,
                                                horizon=1:h,abs_error_moa=abs(fcast-data_future), 
                                                wis_error_moa=wis_tot,
                                                sharpness_moa = sharpness_tot,
                                                overprediction_moa = overprediction_tot,
                                                underprediction_moa = underprediction_tot,
                                                target_end_date = c(curr_targetend_date + 7,
                                                                    curr_targetend_date + 14,
                                                                    curr_targetend_date + 21,
                                                                    curr_targetend_date + 28),
                                                lower_95s = lower_95s,
                                                upper_95s = upper_95s
      )
      mse_df                      <- rbind(mse_df,tmp_df)
      
      dispersion_logs = rbind(dispersion_logs, data.frame(dispersion = dispersion_estimates,
                                                          horizon = c(1:4),
                                                          target_end_date = c(curr_targetend_date + 7,
                                                                              curr_targetend_date + 14,
                                                                              curr_targetend_date + 21,
                                                                              curr_targetend_date + 28)
      )
      )
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
  
  # The following will write the scores_tot.csv file if it has not already been written.
  # source("R/get_model_scores.R")
  scores                          <- read.csv("data/scores_tot_w_wis_components.csv")
  
  write.csv(dispersion_logs, file = paste('data/dispersion_logs/', location, '.csv', sep = ""))
  
  save(state_X, file = paste0('data/embeddings_by_location/', location, '_X.RData'))
  save(state_y, file = paste0('data/embeddings_by_location/', location, '_y.RData'))
  
  #### formatting output of moa
  mse_df_tot                      <- do.call(rbind,mse_df_list)
  
  #### format MOA results to match covidHub forecasts
  mse_df_tot$location_name        <- mse_df_tot$location
  mse_df_tot$location             <- unlist(lapply(mse_df_tot$location,function(x){substr(name_to_fips(x),1,2)}))
  mse_df_tot$location             <- as.integer(mse_df_tot$location)
  
  # The following makes every forecast_date the MONDAY AFTER the target end data that we have
  # data up through.
  mse_df_tot$forecast_date      <- as.Date(mse_df_tot$fcast_date) + 2
  
  ### begin scoring stuff
  scores$target_end_date          <- as.Date(scores$target_end_date)
  scores$forecast_date          <- as.Date(scores$forecast_date)
  
  
  #### keep track of model wins
  results_list                    <- list()
  wins                            <- c()
  wins_wis                        <- c()
  
  # Iterate through every model to compare its performance to sMOA.
  days_data = c()
  for (model_ in unique(scores$model)){
    #### take big precomputed scores data frame exclude the us and subset to current model 
    scores_subset                 <- scores[scores$location !="US" & scores$model == model_,] %>% dplyr::select(horizon,location,abs_error,
                                                                                                                forecast_date,model,wis, 
                                                                                                                dispersion, overprediction,
                                                                                                                underprediction, target_end_date)
    days_data = c(days_data, weekdays(scores_subset$forecast_date))
    #### convert horizon to int
    scores_subset$horizon         <- as.integer(scores_subset$horizon)
    scores_subset$forecast_date         <- as.Date(scores_subset$forecast_date)
    scores_subset$location       <- as.integer(scores_subset$location)
    mse_df_tot$horizon         <- as.integer(mse_df_tot$horizon)
    mse_df_tot$location       <- as.integer(mse_df_tot$location)
    #### if a model submits < 10 forecasts ignore
    if (nrow(scores_subset) > 10){
      ### merge scores by model with moa
      mse_df_tot_gg               <- mse_df_tot %>% left_join(scores_subset,by=c("location","target_end_date","horizon"))
      results_list[[model_]]      <- mse_df_tot_gg
    }
  }
  print(paste("finished", location))
  directory_name = paste("data/", state_log_directory, sep = "")
  save(results_list, file = paste(directory_name, "/", gsub(" ", "", curr_state), ".RData", sep = ""))
  state_coverage_location = paste("data/coverage_data/", location, '.csv', sep = '')
  write.csv(coverage_data, file = state_coverage_location)
  
}


# If the directory exists, change that here.
directory_name = paste("data/", state_log_directory, sep = "")
if(file.exists(directory_name)) unlink(directory_name, recursive = TRUE)


dir.create(directory_name)
sockettype <- "PSOCK"

## Uncomment this to work with a simple example (one run).
#parfctn(3)

cl <- parallel::makeCluster(spec = ncores,type = sockettype) #, outfile=""
setDefaultCluster(cl)
registerDoParallel(cl)
sim_ts <- foreach(i=1:50,
                  .verbose = T)%dopar%{
                    print(i)
                    parfctn(i)
                  }
stopCluster(cl)    





