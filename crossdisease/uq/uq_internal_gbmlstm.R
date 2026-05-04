
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
FILES_ALL = FILES_ALL[grepl('_0.5', FILES_ALL)]
FILES_CUR = FILES_ALL[seq(cur_run,length(FILES_ALL), by = total_run)]

library(mgcv)
library(collapse)
library(qgam)
library(dplyr)

source('/vast/home/lvandervort/GitLab/frodo/uq/quantile_baseline_func.R')


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
  
  
  dat_0.5 <- read.csv(paste0(path_to_dir,FILES_CUR[i]))

  dat_0.025 <- read.csv(paste0(path_to_dir,gsub('_0.5','_0.025',FILES_CUR[i])))
  dat_0.1 <- read.csv(paste0(path_to_dir,gsub('_0.5','_0.1',FILES_CUR[i])))
  dat_0.25 <- read.csv(paste0(path_to_dir,gsub('_0.5','_0.25',FILES_CUR[i])))
  dat_0.75 <- read.csv(paste0(path_to_dir,gsub('_0.5','_0.75',FILES_CUR[i])))
  dat_0.9 <- read.csv(paste0(path_to_dir,gsub('_0.5','_0.9',FILES_CUR[i])))
  dat_0.975 <- read.csv(paste0(path_to_dir,gsub('_0.5','_0.975',FILES_CUR[i])))
  
  if(!('last_obs_time' %in% colnames(dat_0.5))){
    DATES = sort(unique(dat_0.5$date))
    dat_0.5 = merge(dat_0.5, data.frame(date = DATES, last_obs_time = c(1:length(DATES))), by = 'date', all.x = T, all.y = F)
    dat_0.025 = merge(dat_0.025, data.frame(date = DATES, last_obs_time = c(1:length(DATES))), by = 'date', all.x = T, all.y = F)
    dat_0.1 = merge(dat_0.1, data.frame(date = DATES, last_obs_time = c(1:length(DATES))), by = 'date', all.x = T, all.y = F)
    dat_0.25 = merge(dat_0.25, data.frame(date = DATES, last_obs_time = c(1:length(DATES))), by = 'date', all.x = T, all.y = F)
    dat_0.75 = merge(dat_0.75, data.frame(date = DATES, last_obs_time = c(1:length(DATES))), by = 'date', all.x = T, all.y = F)
    dat_0.9 = merge(dat_0.9, data.frame(date = DATES, last_obs_time = c(1:length(DATES))), by = 'date', all.x = T, all.y = F)
    dat_0.975 = merge(dat_0.975, data.frame(date = DATES, last_obs_time = c(1:length(DATES))), by = 'date', all.x = T, all.y = F)
  }
  
  
  PREDS_SORTED = cbind(dat_0.025$fcst, dat_0.1$fcst, dat_0.25$fcst, dat_0.5$fcst, dat_0.75$fcst, dat_0.9$fcst, dat_0.975$fcst )
  PREDS_SORTED2 = t(apply(PREDS_SORTED,1,sort))
  dat_0.025$fcst = pmax(0,PREDS_SORTED2[,1])
  dat_0.1$fcst = pmax(0,PREDS_SORTED2[,2])
  dat_0.25$fcst = pmax(0,PREDS_SORTED2[,3])
  dat_0.5$fcst = pmax(0,PREDS_SORTED2[,4])
  dat_0.75$fcst = pmax(0,PREDS_SORTED2[,5])
  dat_0.9$fcst = pmax(0,PREDS_SORTED2[,6])
  dat_0.975$fcst = pmax(0,PREDS_SORTED2[,7])
  
  g_inds = unique(dat_0.5$last_obs_time)
  dat_0.5$x = dat_0.5$last_obs_time + dat_0.5$h
  g_inds = sort(g_inds)
  #g_inds = g_inds[-c(1:10)] #have at least 10 obs
  
  
  
  if(!file.exists(paste0(output_path,"real_uq_mat_",eval_key, ".csv"))){
    RESULTS = NULL
    for(j in 1:length(g_inds)){
      
      dat_train = dat_0.5[dat_0.5$x <= g_inds[j], ]
      dat_train$error = abs(dat_train$fcst - dat_train$truth)
      dat_test_0.5 = dat_0.5[dat_0.5$last_obs_time == g_inds[j],]

      DAT = data.frame(quant = 0.5, value = dat_0.5[dat_0.5$last_obs_time == g_inds[j],'fcst'], dat_test_0.5[,c('fcst', 'obs', 'truth','date','h')])
      DAT = rbind(DAT,data.frame(quant = 0.025, value = dat_0.025[dat_0.025$last_obs_time == g_inds[j],'fcst'], dat_test_0.5[,c('fcst', 'obs', 'truth','date','h')]))
      DAT = rbind(DAT,data.frame(quant = 0.1, value = dat_0.1[dat_0.1$last_obs_time == g_inds[j],'fcst'], dat_test_0.5[,c('fcst', 'obs', 'truth','date','h')]))
      DAT = rbind(DAT,data.frame(quant = 0.25, value = dat_0.25[dat_0.25$last_obs_time == g_inds[j],'fcst'], dat_test_0.5[,c('fcst', 'obs', 'truth','date','h')]))
      DAT = rbind(DAT,data.frame(quant = 0.75, value = dat_0.75[dat_0.75$last_obs_time == g_inds[j],'fcst'], dat_test_0.5[,c('fcst', 'obs', 'truth','date','h')]))
      DAT = rbind(DAT,data.frame(quant = 0.9, value = dat_0.9[dat_0.9$last_obs_time == g_inds[j],'fcst'], dat_test_0.5[,c('fcst', 'obs', 'truth','date','h')]))
      DAT = rbind(DAT,data.frame(quant = 0.975, value = dat_0.975[dat_0.975$last_obs_time == g_inds[j],'fcst'], dat_test_0.5[,c('fcst', 'obs', 'truth','date','h')]))
      DAT$last_obs_time = j
      DAT$row_num = j

      
      
      ######################
      ### Random Walk UQ ###
      ######################
      
      # dat_train = dat_train[!duplicated(dat_train$last_obs_time),]
      # dat_train = dat_train[order(dat_train$last_obs_time),]
      # QUANT = get_quantiles(nrow(dat_train), dat_train$obs)
      # 
      # DAT_RW = data.frame(last_obs_time = j, row_num = j,  h = QUANT$horizon,
      #                  quant = QUANT$quantile,
      #                  value_rw = QUANT$value)
      # 
      # DAT = merge(DAT, DAT_RW[,c('h','quant','value_rw')], by = c('h','quant'), all.x = T, all.y = T)
      # DAT$value_rw[is.na(DAT$value_rw)] = dat_test_0.5$obs
      
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



