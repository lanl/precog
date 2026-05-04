
################################
### Get Global SLURM Options ###
################################

options(warn = -1)

suppressMessages(library("optparse"))

option_list <- list(
  make_option("--alpha", type="character", default=""),
  make_option("--include_disease_source", type="character", default=""),
  make_option("--warm_start", type="character", default=""),
  make_option("--scale_ts", type="character", default="")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options

alpha = as.numeric(opt$alpha) 
include_disease_source = as.character(opt$include_disease_source) 
warm_start = as.character(opt$warm_start) 
scale_ts = as.character(opt$scale_ts)




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

## define paths
my_path = "~/GitLab/frodo/models"
setwd(my_path)
savetrainpath <- my_path

## define number of cores
ncores <- floor(.9*detectCores())

## define socket type
sockettype <- "PSOCK"


lut_path = '~/GitLab/frodo/lookup_tables/'
DISEASE_SOURCES = readxl::read_xlsx(paste0(lut_path, 'Evaluations_ML.xlsx'))
DISEASE_SOURCES = data.frame(DISEASE_SOURCES)
DISEASE_SOURCES$FILES = gsub('.RDS','',DISEASE_SOURCES$FILES)
DISEASE_SOURCES = DISEASE_SOURCES[DISEASE_SOURCES$FILES == include_disease_source, ]
YEAR_RUNS = c(DISEASE_SOURCES$ML_start:DISEASE_SOURCES$ML_end)


if(alpha == -1){
  ALPHA_GRID = c(0.025, 0.1, 0.25, 0.75, 0.9, 0.975)
}else{
  ALPHA_GRID = c(alpha)
}

for(year_cutoff in YEAR_RUNS){
  
  for(alpha in ALPHA_GRID){
    
    if(scale_ts == 'no'){
      if(warm_start == 'yes'){
        FILENAME = paste0(my_path,'/results_sourcespecific/fit_',year_cutoff,'_', alpha,'_',include_disease_source,'_warm_noscale.RDS')
      }else{
        FILENAME = paste0(my_path,'/results_sourcespecific/fit_',year_cutoff,'_', alpha,'_',include_disease_source,'_cold_noscale.RDS')
      }
    }else{
      if(warm_start == 'yes'){
        FILENAME = paste0(my_path,'/results_sourcespecific/fit_',year_cutoff,'_', alpha,'_',include_disease_source,'_warm.RDS')
      }else{
        FILENAME = paste0(my_path,'/results_sourcespecific/fit_',year_cutoff,'_', alpha,'_',include_disease_source,'_cold.RDS')
      }
    }
    
    if(!file.exists(FILENAME)){
      ################################
      ### Read in Training Dataset ###
      ################################
      
      library(data.table)
      
      PMIN = 0
      
      if(scale_ts == 'yes'){
        real_data = data.table::fread(file =paste0(savetrainpath,"/../features/results_sourcespecific/",include_disease_source,"_train_mat.csv"))
      }else{
        real_data = data.table::fread(file =paste0(savetrainpath,"/../features/results_sourcespecific_unscaled/",include_disease_source,"_train_mat.csv"))
      }
      
      
      if(scale_ts == 'yes'){
        real_data$scaling_factor = as.numeric(real_data$scaling_factor)
        PMAX = 3
      }else{
        real_data$scaling_factor = 1
        PMAX = Inf
      }    
      real_data$truth = pmax(as.numeric(real_data$truth),0)
      
      if('truth_sm' %in% names(real_data)){
        real_data = subset(real_data, select = -c(truth_sm))
      }
      if('truth_before_boot' %in% names(real_data)){
        real_data = subset(real_data, select = -c(truth_before_boot))
      }
      
      
      real_data$obs_smoothed = pmax(as.numeric(real_data$obs_smoothed),0)
      real_data$obs = pmax(real_data$obs,0)
      real_data = real_data[real_data$obs_smoothed > 0,]
      
      real_data$y = pmax(pmin(real_data$truth / real_data$scaling_factor,PMAX),PMIN)
      real_data$obs = pmax(pmin(real_data$obs / real_data$scaling_factor,PMAX),PMIN)
      real_data$obs_smoothed = pmax(pmin(real_data$obs_smoothed / real_data$scaling_factor,PMAX),PMIN)
      real_data$obs_smoothed_minus1 = pmax(pmin(real_data$obs_smoothed_minus1 / real_data$scaling_factor,PMAX),PMIN)
      real_data$sloperoll02 = as.numeric(real_data$sloperoll02)
      real_data$sloperoll12 = as.numeric(real_data$sloperoll12)
      
      YEAR = substr(real_data$date,1,4)
      
      real_data = real_data[YEAR < year_cutoff,]
      real_data = subset(real_data, select = -c(truth,date,obs))
      real_data = real_data[!is.na(real_data$y),]
      
      FILENAME_TEMP = paste0(my_path,'/results/cat_to_numeric.RDS')
      LEVELS_STORAGE <- readRDS(FILENAME_TEMP)
      
      
      real_data$transmission = as.numeric(factor(as.character(real_data$transmission), levels = c(LEVELS_STORAGE$transmission)))
      real_data$viral_bacterial_fungal = as.numeric(factor(as.character(real_data$viral_bacterial_fungal), levels = c(LEVELS_STORAGE$viral_bacterial_fungal)))
      real_data$ts_time_cadence = as.numeric(factor(as.character(real_data$ts_time_cadence), levels = c(LEVELS_STORAGE$ts_time_cadence)))
      real_data$ts_scale = as.numeric(factor(as.character(real_data$ts_scale), levels = c(LEVELS_STORAGE$ts_scale)))
      real_data$ts_measurement_type = as.numeric(factor(as.character(real_data$ts_measurement_type), levels = c(LEVELS_STORAGE$ts_measurement_type)))
      real_data$disease = as.numeric(factor(as.character(real_data$disease), levels = c(LEVELS_STORAGE$disease)))
      real_data$source = as.numeric(factor(as.character(real_data$source), levels = c(LEVELS_STORAGE$source)))
      
      
      COLS_REORDERED = colnames(real_data)
      COLS_REORDERED = sort(COLS_REORDERED)
      real_data = data.frame(real_data)[,COLS_REORDERED]
      
      #########################
      ### Fitting Functions ###
      #########################
      
      lgbmfit <- function(alpha){
        library(lightgbm)
        ## fit nmodels
        model_list <- list()
        for(h in 1:4){
          real_data_h = real_data[real_data$h == h,]
          ## divide df into train and validate
          df_train_ids <- sample(1:nrow(real_data_h), .8*nrow(real_data_h), replace=F)
          df_valid_id <- setdiff(1:nrow(real_data_h),df_train_ids)
          
          CAT_VARS = c('transmission','viral_bacterial_fungal','ts_time_cadence','ts_scale','ts_measurement_type','disease','source')
          
          ## convert the data to LightGBM dataset format: training data
          lgbm_train <- lightgbm::lgb.Dataset(data = as.matrix(subset(real_data_h[df_train_ids,], select=setdiff(names(real_data_h),c("y")))),
                                              label = real_data_h$y[df_train_ids],
                                              categorical_feature = CAT_VARS)
          ## convert the data to LightGBM dataset format: validation data
          lgbm_valid <- lightgbm::lgb.Dataset(data = as.matrix(subset(real_data_h[df_valid_id,], select=setdiff(names(real_data_h),c("y")))),
                                              label = real_data_h$y[df_valid_id],
                                              reference = lgbm_train,
                                              categorical_feature = CAT_VARS)
          
          
          #https://lightgbm.readthedocs.io/en/latest/Parameters.html
          ## define parameters for LightGBM model
          params = list(objective = "quantile",
                        metric = c("quantile"),
                        alpha = alpha,
                        num_leaves = pmin(63,round(nrow(real_data_h)^(1/2),0)), #added pmin, 2^max_depth - 1 common choice
                        max_depth = 6, #from default
                        learning_rate = .1,
                        feature_fraction = 0.9,
                        verbose = -1)
          ## train the LightGBM model
          
          if(warm_start == 'yes'){
            if(scale_ts == 'no'){
              FILENAME_WARM = paste0(my_path,'/results/fit_',year_cutoff,'_', alpha,'_noscale.RDS')
            }else{
              FILENAME_WARM = paste0(my_path,'/results/fit_',year_cutoff,'_', alpha,'.RDS')
            }
            model_warm_start <- readRDS(FILENAME_WARM)
            model_list[[h]] <- lightgbm::lgb.train(params = params,
                                                   data = lgbm_train,
                                                   valids = list(train = lgbm_train,
                                                                 valid = lgbm_valid),
                                                   early_stopping_rounds = 10,
                                                   categorical_feature = CAT_VARS,
                                                   nrounds = 5000,
                                                   init_model = model_warm_start)
          }else{
            model_list[[h]] <- lightgbm::lgb.train(params = params,
                                                   data = lgbm_train,
                                                   valids = list(train = lgbm_train,
                                                                 valid = lgbm_valid),
                                                   early_stopping_rounds = 10,
                                                   categorical_feature = CAT_VARS,
                                                   nrounds = 5000)
          }
          
        }
        ## get outta here
        return(model_list)
      }
      
      
      #######################
      ### Perform Fitting ###
      #######################
      
      mod = lgbmfit(alpha = alpha) 
      
      
      saveRDS(mod, file = FILENAME)
      gc()
    }
  }
  
}
