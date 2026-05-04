
################################
### Get Global SLURM Options ###
################################

options(warn = -1)

suppressMessages(library("optparse"))

option_list <- list(
  make_option("--total_run", type="character", default=""),
  make_option("--cur_run", type="character", default=""),
  make_option("--path_to_dir", type="character", default="")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options

total_run = as.numeric(opt$total_run) 
cur_run = as.numeric(opt$cur_run) 
path_to_dir = as.character(opt$path_to_dir) 








######################
### Body of Script ###
######################


## load libraries
library(ggplot2)
library(data.table)
library(plyr)
library(gridExtra)
library(lubridate)
library(parallel)
library(doParallel)
library(grid)
library(plotly)
library(GGally)
theme_set(theme_bw())


output_path = gsub('evaluations','uq',path_to_dir)
FILES_ALL = list.files(path_to_dir)
FILES_CUR = FILES_ALL[seq(cur_run,length(FILES_ALL), by = total_run)]

library(mgcv)
library(collapse)
library(qgam)


source('/vast/home/lvandervort/GitLab/frodo/uq/quantile_baseline_func.R')


source('/vast/home/lvandervort/GitLab/frodo/uq/covidhub_baseline_custom.R')

get_quantiles = function(x,ts){
  ts_short = ts[1:x]
  
  simple_ts_fit <- fit_quantile_baseline(incidence = ts_short)
  
  forecast_inc_trajectories <- predict.quantile_baseline(
    simple_ts_fit,
    inc_data = ts_short,
    quantiles = c(0.025,0.1,0.25, 0.75,0.9,0.975),
    horizon = 4,
    num_samples = 100000
  )
  forecast_inc_trajectories = forecast_inc_trajectories[forecast_inc_trajectories$type == 'inc',]
  forecast_inc_trajectories[forecast_inc_trajectories < 0] <- 0.0
  forecast_inc_trajectories$last_obs_time = x
  return(forecast_inc_trajectories)
}



get_quantiles_dave = function(x,ts){
  ts_short = ts[1:x]
  
  forecast_inc_trajectories <- covidhub_baseline(y = ts_short,
                                     quantiles = c(0.025,0.1,0.25, 0.75,0.9,0.975),
                                     h = 4,
                                     N = 100000
                                     )

  forecast_inc_trajectories[forecast_inc_trajectories < 0] <- 0.0
  forecast_inc_trajectories$last_obs_time = x
  return(forecast_inc_trajectories)
}

library(dplyr)


## define number of cores
ncores <- pmin(floor(.5*detectCores()),10)

## define socket type
sockettype <- "PSOCK"

## set up cluster
cl <- parallel::makeCluster(spec = ncores,
                            type = sockettype)
setDefaultCluster(cl)
registerDoParallel(cl)

print(Sys.time())
train_data <- foreach(i=1:length(FILES_CUR), 
                      .errorhandling = "pass",
                      .verbose = T,
                      .packages = c('dplyr', 'tsfeatures', 'timeDate', 'lubridate','mgcv','zoo','forecast','collapse','qgam'))%dopar%{  
  eval_key = gsub('.csv','',FILES_CUR[i])
  eval_key = gsub('real_eval_mat_','',eval_key)
  
  
  dat <- read.csv(paste0(path_to_dir,FILES_CUR[i]))
  if('fcst_sm' %in% colnames(dat)){
    dat$fcst = dat$fcst_sm
  }
  if('row_num' %in% colnames(dat)){
    dat$last_obs_time = dat$row_num
  }
  if(!('last_obs_time' %in% colnames(dat))){
    DATES = sort(unique(dat$date))
    dat = merge(dat, data.frame(date = DATES, last_obs_time = c(1:length(DATES))), by = 'date', all.x = T, all.y = F)
  }
  
  dat$fcst = pmax(0,dat$fcst)
  g_inds = unique(dat$last_obs_time)
  dat$x = dat$last_obs_time + dat$h
  dat$scaling_factor = ifelse(dat$fcst==0,rep(0.0001,length(dat$fcst)),dat$fcst)
  g_inds = sort(g_inds)
  g_inds = g_inds[-c(1:10)] #have at least 10 obs
  
  if(!('row_num' %in% colnames(dat))){
    dat$row_num = dat$last_obs_time
  }
  
  if(!file.exists(paste0(output_path,"real_uq_mat_",eval_key, ".csv"))){
    RESULTS = NULL
    for(j in 1:length(g_inds)){
      
      dat_train = dat[dat$x <= g_inds[j], ]
      dat_train$error = abs(dat_train$fcst - dat_train$truth)
      dat_test = dat[dat$last_obs_time == g_inds[j],]
      
      ############################
      ### Negative Binomial UQ ###
      ############################
      h=4
      dispersion_forecast <- 1
      mle_lower_bound <- 1
      mle_start_value <- 0.5
      SUBSET2 = reshape(data.frame(dat_train[,c('error', 'last_obs_time','h')]), direction = 'wide', v.names = 'error', timevar = 'h', idvar = 'last_obs_time')
      SUBSET2 = SUBSET2[,-1]
      if(nrow(SUBSET2) < 6*4){
        sim_nb1 <- rnbinom(5000,mu = (dat_test$fcst[1]+1), size = 20*dispersion_forecast)
        sim_nb2 <- rnbinom(5000,mu = (dat_test$fcst[2]+1), size = 20*0.5*dispersion_forecast)
        sim_nb3 <- rnbinom(5000,mu = (dat_test$fcst[3]+1), size = 20*0.3*dispersion_forecast)
        sim_nb4 <- rnbinom(5000,mu = (dat_test$fcst[4]+1), size = 0.3*dispersion_forecast)
      }else{
        mle_memory <- Inf
        fit_nb_function <- function(k_, h){
          subset_forecasts = dat_train %>% subset(h == h)
          if(nrow(subset_forecasts)>mle_memory){
            subset_forecasts = subset_forecasts[(nrow(subset_forecasts)-mle_memory):nrow(subset_forecasts),]
          }
          mu <- subset_forecasts$truth
          horizon_subset <- subset_forecasts$fcst
          horizon_subset[which(horizon_subset<100)] = mu[which(horizon_subset<100)]
          xx=-sum(dnbinom(round(mu),mu = horizon_subset,size=1/k_,log = T))
          return(xx)
        }
        optimal_k_1 <- optim(mle_start_value,fit_nb_function,method="Brent",lower = 0, upper = 1, h=1)
        optimal_k_2 <- optim(mle_start_value,fit_nb_function,method="Brent",lower = 0, upper = 1, h=2)
        optimal_k_3 <- optim(mle_start_value,fit_nb_function,method="Brent",lower = 0, upper = 1, h=3)
        optimal_k_4 <- optim(mle_start_value,fit_nb_function,method="Brent",lower = 0, upper = 1, h=4)
        dispersion_estimates = c(1/optimal_k_1$par,
                                 1/optimal_k_2$par,
                                 1/optimal_k_3$par,
                                 1/optimal_k_4$par)
        sim_nb1 <- rnbinom(5000,mu = dat_test$fcst[1], size = dispersion_estimates[1]*dispersion_forecast)
        sim_nb2 <- rnbinom(5000,mu = dat_test$fcst[2], size = dispersion_estimates[2]*dispersion_forecast)
        sim_nb3 <- rnbinom(5000,mu = dat_test$fcst[3], size = dispersion_estimates[3]*dispersion_forecast)
        sim_nb4 <- rnbinom(5000,mu = dat_test$fcst[4], size = dispersion_estimates[4]*dispersion_forecast)
      }
      se_preds_nb = c(sd(sim_nb1), sd(sim_nb2), sd(sim_nb3), sd(sim_nb4))

      DAT = expand.grid(last_obs_time = j, row_num = j, h = c(1:4), quant = c(0.025, 0.1, 0.25, 0.5, 0.75, 0.9, 0.975))
      DAT = merge(DAT, data.frame(h=c(1:4), se = se_preds_nb), by = 'h', all.x = T, all.y = T)
      DAT = merge(DAT, dat_test[,c('fcst', 'obs', 'truth','date','h')], by = 'h', all.x = T, all.y = T)
      DAT$value = pmax(0,DAT$fcst + qnorm(p=DAT$quant, lower.tail = T)*DAT$se)
      
      
      
      ######################
      ### Random Walk UQ ###
      ######################
      
      dat_test$fcst_rw = dat_test$obs
      dat_train = dat_train[!duplicated(dat_train$last_obs_time),]
      dat_train = dat_train[order(dat_train$last_obs_time),]
      QUANT = get_quantiles(nrow(dat_train), dat_train$obs)
      #QUANT_DAVE = get_quantiles_dave(nrow(dat_train), dat_train$obs)
      
      DAT_RW = data.frame(last_obs_time = j, row_num = j,  h = QUANT$horizon,
                       quant = QUANT$quantile,
                       value_rw = QUANT$value)
      DAT_RW = rbind(DAT_RW, data.frame(last_obs_time = j, row_num = j,  h = c(1:4),
                                        quant = rep(0.5,4),
                                        value_rw =dat_test$obs))
      DAT = merge(DAT, DAT_RW[,c('h','quant','value_rw')], by = c('h','quant'), all.x = T, all.y = T)
      
      RESULTS <- rbind(RESULTS, DAT)
      print(j)
    }
    write.csv(RESULTS, file = paste0(output_path,"real_uq_mat_",eval_key, ".csv"), quote = F, row.names = F)
  }
  rm('RESULTS')
  xx <- 1
  xx
}
stopCluster(cl)



