# Optimize input parameters for sMOA using a Bayesian optimization scheme.
## Author: Alexander C. Murph
## Date: August 2024
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
library(collapse)
library(ParBayesianOptimization)
setwd("~/GitLab/smoa")
source("R/smoa_helpers.R")
source('R/vecchia_scaled.R')
source("R/parallel_bayesian_optimization.R")
set.seed(13)

num_cores                            <- 99
curve_range                          <- c(10000L, 30000L)
k_range 			     <- c(2L, 7L)
closest_ids_range                    <- c(1000L, 6000L)
dispersion_scale_range               <- c(10, 40)
mle_lower_bound_range	             <- c(0.05, 5)

param_bounds = list(num_curves = curve_range, 
                    k = k_range, 
                    closest = closest_ids_range, 
                    dispersion_scale = dispersion_scale_range,
		    mle_lower_bound_range = mle_lower_bound_range)

is_discrete = c(1,1,1,0,0)

fit_model_wrapper = function(input_vector){
  num_curves                         <- input_vector[1]
  k                                  <- input_vector[2]
  closest                            <- input_vector[3]
  dispersion_forecast                <- input_vector[4]
  lower_mle_value 		     <- input_vector[5]
  lower_CI_scale                     <- 1
  upper_CI_scale                     <- 1
  library(collapse)
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
  library(deSolve)
  source("R/smoa_helpers.R")
  source('R/vecchia_scaled.R')
  h                                  <- 4
  types_of_curves                    <- c("sir_rollercoaster", "sir_rollercoaster_wiggle", 'seasonal')
  N                                  <- num_curves*length(types_of_curves)
  quantiles                          <- c(0.01, 0.025, seq(0.05, 0.95, by = 0.05), 0.975, 0.99)
  
  ############################################
  ### Generate Synthetic List for Training ###
  ############################################
  sim_ts                             <- foreach(i=1:N,
                    .errorhandling = "pass",
                     .verbose = F)%do%{
                       library(plyr)
                       library(data.table)
                       library(LearnBayes)
                       library(LaplacesDemon)
                       source("R/smoa_helpers.R")
                       
                       curve_type     <- types_of_curves[rep(c(1:length(types_of_curves)), each = num_curves)][i]
                       templist       <- gen_curve(curve_type)
                       templist
                     }
   A = lapply(sim_ts,function(x){ length(x)})
   unique(A)
   sim_ts = sim_ts[A > 2]
   
   print("Creating the embedding matrix.")
   embed_mat                          <- create_embed_matrix(sim_ts,h=(h+1),k=(k+1))#I want k to represent the number we compare on the diff scale
   
   #### read in embedding mat design matrix and outcome matrix
   X                                  <- embed_mat[[1]]
   y                                  <- embed_mat[[2]]
   
   ### convert to differences
   X_diff                             <- t(apply(X,1,diff))
   y_diff                             <- t(apply(y,1,diff))
   
   # Create the test synthetic data.
   xx                                 <- length(types_of_curves[rep(c(1:length(types_of_curves)), each = floor(num_curves/100))])
   sim_ts_test                        <- foreach(i=1:xx,
                     .errorhandling = "pass",
                     .verbose = F)%do%{
                       library(plyr)
                       library(data.table)
                       library(LearnBayes)
                       library(LaplacesDemon)
                       source("R/smoa_helpers.R")
                       
                       curve_type     <- types_of_curves[rep(c(1:length(types_of_curves)), each = floor(num_curves/100) )][i]
                       templist       <- gen_curve(curve_type)
                       templist
                     }
   A = lapply(sim_ts_test,function(x){ length(x)})
   unique(A)
   sim_ts_test <- sim_ts_test[A > 2]
   sim_ts_test_new <- list()
   list_count <- 1
   shuffled_indices <- sample(1:length(sim_ts_test), size = length(sim_ts_test))
   for(list_idx in shuffled_indices){
     sim_ts_test_new[[list_count]] <- sim_ts_test[[list_idx]]
     list_count <- list_count + 1
   }
   sim_ts_test <- sim_ts_test_new
   sim_ts_test_new <- NULL
 
  #save(X, file = "X.RData")
  #save(y, file = "y.RData")
  #save(X_test, file = "X_test.RData")
  #save(y_test, file = "y_test.RData")
  #save(X_diff, file = "X_diff.RData")
  #save(y_diff, file = "y_diff.RData")
  #save(sim_ts_test, file = 'sim_ts_test.RData')
  #save(X_diff_test, file = "X_diff_test.RData")
  #save(y_diff_test, file = "y_diff_test.RData")

  #load(file = "X.RData")
  #load(file = "y.RData")
  #load(file = "X_test.RData")
  #load(file = "y_test.RData")
  #load(file = "X_diff.RData")
  #load(file = "y_diff.RData")
  #load(file = 'sim_ts_test.RData')
  #load(file = "X_diff_test.RData")
  #load(file = "y_diff_test.RData")

  wis_dump                           <- c()
  mae_dump                           <- c()
  full_testing_count                 <- 1

  one_step_ahead_forecasts                                  <- c()
  two_step_ahead_forecasts                                  <- c()
  three_step_ahead_forecasts                                <- c()
  four_step_ahead_forecasts                                 <- c()

  for (test_idx in 1:length(sim_ts_test)){
    print(paste("starting calculation", test_idx))
    if(full_testing_count > 300) break
    
    full_testing_count <- full_testing_count + 1

    if((length(sim_ts_test[[test_idx]])-h-1)>(k+1)){
      sim_ts_test[[test_idx]]$ts <- sim_ts_test[[test_idx]]$ts + 100
      
      for(end_date_index in (k+1):(length(sim_ts_test[[test_idx]]$ts)-h-1)){
        #### make sure that the target end date is less than or equal to fcast_date
        data_till_now                    <- data.frame(value = sim_ts_test[[test_idx]]$ts[1:end_date_index])
        data_till_now$value              <- pmax(0,data_till_now$value)
        data_till_now$t                  <- 1:nrow(data_till_now)
    	
	if(all(data_till_now$value<=1e-8)) next

        #### light smoothing and differencing and get last k 
        data_till_now_smoothed           <- gam(value~ s(t,k=round(nrow(data_till_now)/2)),data=data_till_now)$fitted.values
    
        #### could use gam smoother 
        to_match_in_moa                  <- tail((data_till_now_smoothed),k+1)
        
        ##### take our matrix to match (X_diff) and efficiently compute distances from to_match_in_moa to X_diff
        ##### have to get the scale right 
        gc()
        dist_to_test                     <- rowSums(abs(X_diff %r-% diff(to_match_in_moa)))
        gc()
        
        #### get the closest ids in the X_diff mat
        closest_ids                      <- sort(dist_to_test,index.return = TRUE)$ix[1:closest]
        dists_of_closest                 <- sort(dist_to_test)[1:closest]
        
        #### forecast is then just y_diff mat subsetted to the closest ids and then take median over first h horizons
        # decay_weights <-  1 / dists_of_closest
        # decay_weights <- decay_weights / sum(decay_weights)
        # closest_ids[decay_weights < epsilon] <- 0
        fcast                            <- head(apply(y_diff[closest_ids,],2,median),h)
        
        #### convert from difference back to raw cases
        fcast                            <- tail(data_till_now$value,1) + cumsum(fcast)
        fcast                            <- pmax(0,fcast)
        data_future                      <- sim_ts_test[[test_idx]]$ts[(end_date_index+1):(end_date_index+h)] # Note that that '7' is for the length of a week (7 days).
     
    mae_dump = c(mae_dump, abs(fcast-data_future))

    #### UQ 
    point                            <- pmax(0,tail(data_till_now$value,1) +cumsum(head(apply(y_diff[closest_ids,],2,median),h)))
    vars                             <- pmax(0,tail(data_till_now$value,1) +cumsum(head(apply(y_diff[closest_ids,],2,sd),h)))
  
   
    if(length(one_step_ahead_forecasts) < 8){
    	sim_nb1 <- matrix(rnbinom(5000*length(point),mu = (point+1), size = 4*dispersion_forecast),ncol=length(point),byrow = T)
    	sim_nb2 <- matrix(rnbinom(5000*length(point),mu = (point+1), size = 3*dispersion_forecast),ncol=length(point),byrow = T)
    	sim_nb3 <- matrix(rnbinom(5000*length(point),mu = (point+1), size = 2*dispersion_forecast),ncol=length(point),byrow = T)
    	sim_nb4 <- matrix(rnbinom(5000*length(point),mu = (point+1), size = 1*dispersion_forecast),ncol=length(point),byrow = T)
    } else {
	
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
        optimal_k_1 <- optim(30,fit_nb_function,method="Brent",lower = lower_mle_value, upper = 2000)
        
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
        optimal_k_2 <- optim(30,fit_nb_function,method="Brent",lower = lower_mle_value, upper = 2000)
        
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
        optimal_k_3 <- optim(30,fit_nb_function,method="Brent",lower = lower_mle_value, upper = 2000)
        
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
        optimal_k_4 <- optim(30,fit_nb_function,method="Brent",lower = lower_mle_value, upper = 2000)
       
	sim_nb1 <- matrix(rnbinom(5000*length(point),mu = (point+1), size = optimal_k_1$par*dispersion_forecast),ncol=length(point),byrow = T)
        sim_nb2 <- matrix(rnbinom(5000*length(point),mu = (point+1), size = optimal_k_2$par*dispersion_forecast),ncol=length(point),byrow = T)
        sim_nb3 <- matrix(rnbinom(5000*length(point),mu = (point+1), size = optimal_k_3$par*dispersion_forecast),ncol=length(point),byrow = T)
        sim_nb4 <- matrix(rnbinom(5000*length(point),mu = (point+1), size = optimal_k_4$par*dispersion_forecast),ncol=length(point),byrow = T)

    }
    #### UQ 
    #point <- pmax(0,tail(data_till_now$value,1) +cumsum(head(apply(y_diff[closest_ids,],2,median),h)))
    #vars <- pmax(0,tail(data_till_now$value,1) +cumsum(head(apply(y_diff[closest_ids,],2,sd),h)))
    
    est_intervals <- quantile(sim_nb1[,1],probs = quantiles)
    est_intervals[quantiles<0.5] <- est_intervals[quantiles<0.5] * lower_CI_scale
    est_intervals[quantiles>0.5] <- est_intervals[quantiles>0.5] * upper_CI_scale
    wis_tmp_1 <- weighted_interval_score(quantiles,value = est_intervals, actual_value = data_future[1])
    
    est_intervals <- quantile(sim_nb2[,2],probs = quantiles)
    est_intervals[quantiles<0.5] <- est_intervals[quantiles<0.5] * lower_CI_scale
    est_intervals[quantiles>0.5] <- est_intervals[quantiles>0.5] * upper_CI_scale
    wis_tmp_2 <- weighted_interval_score(quantiles,value = est_intervals, actual_value = data_future[2])
    
    est_intervals <- quantile(sim_nb3[,3],probs = quantiles)
    est_intervals[quantiles<0.5] <- est_intervals[quantiles<0.5] * lower_CI_scale
    est_intervals[quantiles>0.5] <- est_intervals[quantiles>0.5] * upper_CI_scale
    wis_tmp_3 <- weighted_interval_score(quantiles,value = est_intervals, actual_value = data_future[3])
    
    est_intervals <- quantile(sim_nb4[,4],probs = quantiles)
    est_intervals[quantiles<0.5] <- est_intervals[quantiles<0.5] * lower_CI_scale
    est_intervals[quantiles>0.5] <- est_intervals[quantiles>0.5] * upper_CI_scale
    wis_tmp_4 <- weighted_interval_score(quantiles,value = est_intervals, actual_value = data_future[4])
    
    wis_dump <- c(wis_dump, wis_tmp_1,wis_tmp_2,wis_tmp_3,wis_tmp_4)
    one_step_ahead_forecasts                                  <- c(one_step_ahead_forecasts, fcast[1])
    two_step_ahead_forecasts                                  <- c(two_step_ahead_forecasts, fcast[2])
    three_step_ahead_forecasts                                <- c(three_step_ahead_forecasts, fcast[3])
    four_step_ahead_forecasts                                 <- c(four_step_ahead_forecasts, fcast[4])
  }
  }
  }
  return(list(Score = mean(wis_dump, na.rm = T)))
}

num_cores = 99

print("starting example version")
example_version = c(param_bounds[[1]][1],param_bounds[[2]][1], param_bounds[[3]][1], param_bounds[[4]][1], param_bounds[[5]][1])
print(fit_model_wrapper(example_version))
print('finished example version')

stop("just testing.")


opt_model = parallel_bayesian_optimization(fit_model_wrapper,
                     parameter_bounds = param_bounds,
                     initial_lhs = num_cores,
                     number_of_iters = 5, 
                     is_discrete = is_discrete,
		     num_cores = num_cores)

print(opt_model)

write.csv(data.frame(best_params=opt_model) , file = paste('murph_best_params', '.csv', sep = ""))
