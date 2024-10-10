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
  sim_idx                       <- 3
}
curr_state                      <- state.name[sim_idx]


num_of_each_curve <- c(20000, 30000, 40000, 50000)
ks <- c(4, 6,7,8)
epsilons <- c(5e-4, 5e-3, 5e-5, 5e-2)
closest_ids_vec <- c(1000, 2000, 3000, 5000)


h                               <- 4

# These are all from the bayesian optimization.
num_curves <- 18387
k <- 7
closest <- 4422
lower_CI_scale <- 1.038942204
upper_CI_scale <- 0.81145339130
dispersion_forecast <- 1
mle_lower_bound <- 1

###########
# Make a directory to hold the results for each state, using this run's particular set of hyperparameters.
name_of_change <- paste("wHistories_k", k, "num_curves", num_curves, "closest", closest, "dispersion", round(dispersion_forecast*10000), 'mlebound', round(mle_lower_bound*10000), sep = "_")
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
                      #templist = gen_curve(curve_type)
                      library(plyr)
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

### convert to differences
X_diff                          <- as.data.frame(t(apply(X,1,diff)))
y_diff                          <- as.data.frame(t(apply(y,1,diff)))

original_y_diff <- y_diff
original_X_diff <- X_diff

parfctn = function(x){
  library(covidHubUtils)
  library(mgcv)
  library(collapse)
  library(dplyr)
  quantiles <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
  ### read in stored truth data with as_of
  
  # We pre-built this data file to cut on api calls to github.
  truth_as_of_tot                 <- read.csv("data/tdat_list_tot_weekly.csv")
  
  ### Iterate through the states and calculate the MAE and WIS for the sMOA forecast.
  mse_df_list                     <- list()
  count_list                      <- 1
  curr_state <- state.name[x]
  mle_start_value <- 20

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
    fcast_dates_to_match          <- fcast_dates_to_match[which(fcast_dates_to_match > as.Date("2020-08-15"))]
    
    #### subset big as of data frame to this specific location
    truth_as_of_tot_loc           <- truth_as_of_tot[truth_as_of_tot$location_name == location ,]
    
    ### hold mse 
    mse_df                        <- NULL 
    
    ### Reset the X_diff y_diff (we added in histories in this version)
    X_diff <- original_X_diff
    y_diff <- original_y_diff

    ### iterate through forecast dates
    added_first_history_flag <- FALSE
    for (fcast_date_idx in 1:(length(fcast_dates_to_match))){
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

      # At this point, we are assuming that we have 'new data.' Tack it onto the data library:
      if((!added_first_history_flag)&(length(data_till_now$value)>(h+k+1))){
      	# In this case, I need to segment all the data that we observe at the beginning.
      	temp_histories <- list(list(ts = data_till_now$value))
      	embed_mat_histories                       <- create_embed_matrix(temp_histories,h=(h+1),k=(k+1))
      	X_temp                               <- embed_mat_histories[[1]]
      	y_temp                               <- embed_mat_histories[[2]]
      
      	### convert to differences
      	X_diff_temp                          <- data.frame(t(apply(X_temp,1,diff)))
      	y_diff_temp                          <- data.frame(t(apply(y_temp,1,diff)))
      	
      	gc()
      	X_diff <- collapse::rowbind(X_diff, X_diff_temp, use.names=FALSE)
      	y_diff <- collapse::rowbind(y_diff, y_diff_temp, use.names=FALSE) 
      	gc()
      
      	added_first_history_flag <- TRUE
      } else if((length(data_till_now$value)>(h+k+1))) {
      	# In this case, I only need to take on a single new 'segment' onto the data matrices.
      	new_X_vals <- data_till_now$value[(length(data_till_now$value)-h-k-1):(length(data_till_now$value)-h-1)]
              new_y_vals <- data_till_now$value[(length(data_till_now$value)-h):(length(data_till_now$value))]
      	
      	### convert to differences
              X_diff_temp                          <- data.frame(t(as.matrix(diff(new_X_vals))))
              y_diff_temp                          <- data.frame(t(as.matrix(diff(new_y_vals))))
      
      	gc()
              X_diff                               <- collapse::rowbind(X_diff, X_diff_temp, use.names=FALSE)
              y_diff  			     <- collapse::rowbind(y_diff, y_diff_temp, use.names=FALSE)
      	gc()
      }

      #### light smoothing and differencing and get last k 
      data_till_now_smoothed      <- gam(value~ s(t,k=round(nrow(data_till_now)/2)),data=data_till_now)$fitted.values
      
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
      if ((length(one_step_ahead_forecasts) < 7)){
	      sim_nb1 <- matrix(rnbinom(5000*length(point),mu = (point+1), size = 1*dispersion_forecast),ncol=length(point),byrow = T)
        sim_nb2 <- matrix(rnbinom(5000*length(point),mu = (point+1), size = 0.5*dispersion_forecast),ncol=length(point),byrow = T)
        sim_nb3 <- matrix(rnbinom(5000*length(point),mu = (point+1), size = 0.25*dispersion_forecast),ncol=length(point),byrow = T)
        sim_nb4 <- matrix(rnbinom(5000*length(point),mu = (point+1), size = 0.175*dispersion_forecast),ncol=length(point),byrow = T)

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
        optimal_k_1 <- optim(mle_start_value,fit_nb_function,method="Brent",lower = mle_lower_bound, upper = 2000)
        
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
        optimal_k_2 <- optim(mle_start_value,fit_nb_function,method="Brent",lower = mle_lower_bound, upper = 2000)
        
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
        optimal_k_3 <- optim(mle_start_value,fit_nb_function,method="Brent",lower = mle_lower_bound, upper = 2000)
        
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
        optimal_k_4 <- optim(mle_start_value,fit_nb_function,method="Brent",lower = mle_lower_bound, upper = 2000)
      
	dispersion_estimates = c(optimal_k_1$par,
                               optimal_k_2$par,
                               optimal_k_3$par,
                               optimal_k_4$par)
        for(horizon_num in 2:h){
          if(dispersion_estimates[horizon_num] >= 1990){
            dispersion_estimates[horizon_num] <- dispersion_estimates[horizon_num-1]*0.5
          }
        }

	temp_point <- point
        for(horizon_idx in 2:h){
          if(temp_point[horizon_idx]<200) temp_point[horizon_idx] <- temp_point[horizon_idx-1]
        }
        sim_nb1 <- matrix(rnbinom(5000*length(point),mu = (temp_point), size = dispersion_estimates[1]*dispersion_forecast),ncol=length(point),byrow = T)
        sim_nb2 <- matrix(rnbinom(5000*length(point),mu = (temp_point), size = dispersion_estimates[2]*dispersion_forecast),ncol=length(point),byrow = T)
        sim_nb3 <- matrix(rnbinom(5000*length(point),mu = (temp_point), size = dispersion_estimates[3]*dispersion_forecast),ncol=length(point),byrow = T)
        sim_nb4 <- matrix(rnbinom(5000*length(point),mu = (temp_point), size = dispersion_estimates[4]*dispersion_forecast),ncol=length(point),byrow = T)


        print("optimal ks are")
        print(optimal_k_1$par)
        print(optimal_k_2$par)
        print(optimal_k_3$par)
        print(optimal_k_4$par)
        print("and the forecasts are")
        print(point)
        print(fcast)
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
      
      lower_CI_scale <- 1
      upper_CI_scale <- 1
      
      est_intervals <- quantile(sim_nb1[,1],probs = quantiles)
      est_intervals[quantiles<0.5] <- est_intervals[quantiles<0.5] * lower_CI_scale
      est_intervals[quantiles>0.5] <- est_intervals[quantiles>0.5] * upper_CI_scale
      lower_95s <- c(lower_95s, est_intervals[2])
      upper_95s <- c(upper_95s, est_intervals[22])
      wis_tmp_1 <- weighted_interval_score(quantiles,value = est_intervals, actual_value = data_future[1])
      
      temp_est_intervals <- est_intervals
      est_intervals <- quantile(sim_nb2[,2],probs = quantiles)
      if(any(est_intervals<=10)){
        est_intervals[est_intervals<=10] <- temp_est_intervals[est_intervals<=10]
      }
      est_intervals[quantiles<0.5] <- est_intervals[quantiles<0.5] * lower_CI_scale
      est_intervals[quantiles>0.5] <- est_intervals[quantiles>0.5] * upper_CI_scale
      wis_tmp_2 <- weighted_interval_score(quantiles,value = est_intervals, actual_value = data_future[2])
      lower_95s <- c(lower_95s, est_intervals[2])
      upper_95s <- c(upper_95s, est_intervals[22])
      
      temp_est_intervals <- est_intervals
      est_intervals <- quantile(sim_nb3[,3],probs = quantiles)
      if(any(est_intervals<=10)){
        est_intervals[est_intervals<=10] <- temp_est_intervals[est_intervals<=10]
      }
      est_intervals[quantiles<0.5] <- est_intervals[quantiles<0.5] * lower_CI_scale
      est_intervals[quantiles>0.5] <- est_intervals[quantiles>0.5] * upper_CI_scale
      wis_tmp_3 <- weighted_interval_score(quantiles,value = est_intervals, actual_value = data_future[3])
      lower_95s <- c(lower_95s, est_intervals[2])
      upper_95s <- c(upper_95s, est_intervals[22])
      
      temp_est_intervals <- est_intervals
      est_intervals <- quantile(sim_nb4[,4],probs = quantiles)
      if(any(est_intervals<=10)){
        est_intervals[est_intervals<=10] <- temp_est_intervals[est_intervals<=10]
      }
      est_intervals[quantiles<0.5] <- est_intervals[quantiles<0.5] * lower_CI_scale
      est_intervals[quantiles>0.5] <- est_intervals[quantiles>0.5] * upper_CI_scale
      wis_tmp_4 <- weighted_interval_score(quantiles,value = est_intervals, actual_value = data_future[4])
      lower_95s <- c(lower_95s, est_intervals[2])
      upper_95s <- c(upper_95s, est_intervals[22])
      
      # # Now plot all of this for Casey:
      all_data = c(data_till_now_smoothed, data_future)
      graph_data = data.frame(x = 1:length(all_data), 
                              value = all_data,
                              type = rep('true_data', times = length(all_data)))
      graph_data = rbind(graph_data, data.frame(x = (length(data_till_now_smoothed)):length(all_data), 
                                                value = c(data_till_now_smoothed[length(data_till_now_smoothed)], point),
                                                type = rep('mean_value', times = (h+1))))
      graph_data = rbind(graph_data, data.frame(x = (length(data_till_now_smoothed)):length(all_data), 
                                                value = c(data_till_now_smoothed[length(data_till_now_smoothed)],upper_95s),
                                                type = rep('upper 95 CI', times = (h+1))))
      graph_data = rbind(graph_data, data.frame(x = (length(data_till_now_smoothed)):length(all_data), 
                                                value = c(data_till_now_smoothed[length(data_till_now_smoothed)], lower_95s),
                                                type = rep('lower 95 CI', times = (h+1))))
      
      forecast_horizon = length(data_till_now_smoothed)
      
      
      temp_data = data.frame(x = (length(data_till_now_smoothed)):length(all_data),
                             upper95 = c(data_till_now_smoothed[length(data_till_now_smoothed)],upper_95s),
                             lower95 = c(data_till_now_smoothed[length(data_till_now_smoothed)], lower_95s),
                             value = c(data_till_now_smoothed[length(data_till_now_smoothed)], point),
                             type = rep('error_region', times = (h+1)))
      
      p1 <- ggplot(graph_data, aes(x = x, y = value, color = type)) + geom_line() + geom_vline(xintercept = forecast_horizon) +
        geom_ribbon(data = temp_data, aes(x = x, y = value, ymax = upper95, ymin = lower95), alpha = 0.2)
      print(p1)
      Sys.sleep(0.1)
      
      wis_tot <- c(wis_tmp_1,wis_tmp_2,wis_tmp_3,wis_tmp_4)
      
      ### compute mse and append
      ### could undo horizon mean and include horizon 
      tmp_df                      <- data.frame(location=location,target_end_date=fcast_date,horizon=1:h,abs_error_moa=abs(fcast-data_future), wis_error_moa=wis_tot)
      mse_df                      <- rbind(mse_df,tmp_df)

      # Remove abnormally low errors:
      one_step_ahead_forecasts <- one_step_ahead_forecasts[one_step_ahead_forecasts>100]
      two_step_ahead_forecasts <- two_step_ahead_forecasts[two_step_ahead_forecasts>100]
      three_step_ahead_forecasts <- three_step_ahead_forecasts[three_step_ahead_forecasts>100]
      four_step_ahead_forecasts <- four_step_ahead_forecasts[four_step_ahead_forecasts>100]
    }
    
    #### create the data frame and return
    mse_df_list[[count_list]]     <- mse_df
    count_list                    <- count_list + 1
  }
  
  # The following will write the scores_tot.csv file if it has not already been written.
  # source("R/get_model_scores.R")
  scores                          <- read.csv("data/scores_tot.csv")
  
  
  #### formatting output of moa
  mse_df_tot                      <- do.call(rbind,mse_df_list)
  
  
  #### format MOA results to match covidHub forecasts
  mse_df_tot$location_name        <- mse_df_tot$location
  mse_df_tot$location             <- unlist(lapply(mse_df_tot$location,function(x){substr(name_to_fips(x),1,2)}))
  mse_df_tot$model                <- "moa"
  mse_df_tot$location             <- as.integer(mse_df_tot$location)
  mse_df_tot$target_end_date      <- as.Date(mse_df_tot$target_end_date)
  
  ### begin scoring stuff
  scores$target_end_date          <- as.Date(scores$target_end_date)
  
  
  #### keep track of model wins
  results_list                    <- list()
  wins                            <- c()
  wins_wis                        <- c()
  
  # Iterate through every model to compare its performance to sMOA.
  for (model_ in unique(scores$model)){
    
    #### take big precomputed scores data frame exclude the us and subset to current model 
    scores_subset                 <- scores[scores$location !="US" & scores$model == model_,] %>% select(horizon,location,abs_error,target_end_date,model,wis)
    
    #### convert horizon to int
    scores_subset$horizon         <- as.integer(scores_subset$horizon)
    #### if a model submits < 10 forecasts ignore
    if (nrow(scores_subset) > 10){
      ### merge scores by model with moa
      mse_df_tot_gg               <- mse_df_tot %>% left_join(scores_subset,by=c("location","target_end_date","horizon"))
      #mse_df_tot_gg               <- mse_df_tot_gg[complete.cases(mse_df_tot_gg),]
      results_list[[model_]]      <- mse_df_tot_gg
      
      # We should really only be making comparisons for places where both models have predictions.
      temp_df <- mse_df_tot_gg[which(!is.na(mse_df_tot_gg$wis)),]

      if(nrow(temp_df) > 0){
      	moa_median_wis <- mean(temp_df$wis_error_moa)
      	covid_median_wis <- mean(temp_df$wis)
      	print('testing:')
      	print(moa_median_wis)
      	print(covid_median_wis)
      	if((length(covid_median_wis) >0)&(length(moa_median_wis) >0)){
      	  if ((moa_median_wis) < (covid_median_wis)){
      	    wins_wis                    <- c(wins_wis,1)
      	  } else{
      	    wins_wis                    <- c(wins_wis,0)
      	  }
      	}
      } else {
	moa_median_wis <- NA
        covid_median_wis <- NA
      }
      
      # We should really only be making comparisons for places where both models have predictions.
      temp_df <- mse_df_tot_gg[which(!is.na(mse_df_tot_gg$abs_error)),]
      if(nrow(temp_df) > 0){
      	moa_median_abs <- mean(temp_df$abs_error_moa)
      	covid_median_abs <- mean(temp_df$abs_error)
      	
      	if((length(covid_median_abs) >0)&(length(moa_median_abs)>0)){
      	  if (moa_median_abs < covid_median_abs){
      	    wins                    <- c(wins,1)
      	  } else{
      	    wins                    <- c(wins,0)
      	  }
      	}
      } else {
	moa_median_abs <- NA
        covid_median_abs <- NA
      }
      
      results_list[[paste("justresults_", model_, sep = "")]] <- data.frame(moa_abs = moa_median_abs, covidhub_abs = covid_median_abs,
                                                                            moa_wis = moa_median_wis, covidhub_wis = covid_median_wis)
    }
  }
  
  print(mean(wins))
  print(mean(wins_wis))
  results_list[['proportion_wins_abs']] <- mean(wins)
  results_list[['wins_vector_abs']] <- wins
  results_list[['proportion_wins_wis']] <- mean(wins_wis)
  results_list[['wins_vector_wis']] <- wins_wis
  results_list[['name_of_change']] <- name_of_change
  
  print(paste("finished", location))
  directory_name = paste("data/", state_log_directory, sep = "")
  save(results_list, file = paste(directory_name, "/", gsub(" ", "", curr_state), ".RData", sep = ""))    
}


# If the directory exists, change that here.
directory_name = paste("data/", state_log_directory, sep = "")
if(file.exists(directory_name)) unlink(directory_name, recursive = TRUE)

dir.create(directory_name)
sockettype <- "PSOCK"

parfctn(3)
stop("just testing.")

## Uncomment this to work with a simple example (one run).
cl <- parallel::makeCluster(spec = ncores,type = sockettype, outfile="")
setDefaultCluster(cl)
registerDoParallel(cl)
sim_ts <- foreach(i=1:50,
                  .verbose = T)%do%{
                    
                    parfctn(i)
                  }
stopCluster(cl)    


print(paste("The results for the change described as",name_of_change , "are as follows:"))
source("R/compile_results.R")

xx = compile_results(directory_name)
xx$k = k
#xx$epsilon = epsilon
xx$how_many_closest = closest
xx$num_curves = num_curves
xx$length_synthetic = nrow(X_diff)
xx$upper_CI_scale = upper_CI_scale
xx$lower_CI_scale = lower_CI_scale
xx$dispersion_scale = dispersion_forecast
xx$mle_lower_bound = mle_lower_bound

print(directory_name)

write.csv(xx, file=paste("data/results/", name_of_change, ".csv", sep = ""))

