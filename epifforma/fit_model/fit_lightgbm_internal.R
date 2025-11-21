##########################
##########################
### epiFFORMA Training ###
##########################
##########################


################################
### Get Global SLURM Options ###
################################

options(warn = -1)

suppressMessages(library("optparse"))

option_list <- list(
  make_option("--fit_type", type="character", default=""),
  make_option("--fit_num", type="character", default=""),
  make_option("--param_type", type="character", default="")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options

fit_type = as.character(opt$fit_type) 
fit_num = as.character(opt$fit_num) 
param_type = as.character(opt$param_type) 

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
my_path = here::here("epifforma", "fit_model") #this.path::here()
# setwd(my_path)
savetrainpath <- paste0(my_path, '/../process_data/')
syntheticpath <- paste0(my_path, '/../raw_data/synthetic/output/') #EDITED 3/25/24
figsavepath <- paste0(savetrainpath,"figs/")
savemodelpath <- paste0(my_path, '/features2weights/') 

## source in function
source(paste0(my_path, "/../process_data/epi_functions.R"))

#### set up computer for parallelization
## define number of cores
ncores <- floor(.9*detectCores())

## define socket type
sockettype <- "PSOCK"


## forecast horizon
h = 4

## read in data
  #############
if(fit_type == 'order'){ #class = model ORDER (lowest = 1, highest = 7)
  #############
  lf <- list.files(paste0(savetrainpath,"training_data_packets/"))[grep("features_se_",list.files(paste0(savetrainpath,"training_data_packets/")))]
  FILENAME = paste0(savemodelpath,"fitted_lgb_models_order_",fit_num, ".RDS")
  train_data <- NULL
  for(i in 1:length(lf)){
    train_data <- rbind(train_data,fread(paste0(paste0(savetrainpath,"training_data_packets/"),lf[i])))
  }
  train_data = train_data[!duplicated(paste0(train_data$ts_length, '_', train_data$ts_id,'_', train_data$class, '_', train_data$h)),] #added 7/7/24
  
  comp_ex <- list.files(paste0(savetrainpath,"training_data_packets/"))[grepl("component",list.files(paste0(savetrainpath,"training_data_packets/"))) & !grepl("components_minus",list.files(paste0(savetrainpath,"training_data_packets/")))]
  component_data <- NULL
  for(i in 1:length(comp_ex)){
    component_data <- rbind(component_data,fread(paste0(paste0(savetrainpath,"training_data_packets/"),comp_ex[i])))
  } 
  COMPONENTS = setdiff(names(component_data),c('truth','h','ts_length','ts_id'))
  ORDERS = apply(as.matrix(data.frame(component_data)[,COMPONENTS]),1,FUN = order)
  ORDERS = data.frame(t(ORDERS))
  ORDER_BY_COMPONENTS = matrix(NA, ncol = length(COMPONENTS), nrow = length(ORDERS[,1]))
  for(i in 1:length(COMPONENTS)){
    ORDER_BY_COMPONENTS[,i] = apply(ORDERS,1,FUN = function(x,i){which(x==i)}, i = i)
  }
  ORDER_BY_COMPONENTS = data.frame(ORDER_BY_COMPONENTS)
  colnames(ORDER_BY_COMPONENTS) = COMPONENTS
  components_long = NULL
  for(i in 1:length(COMPONENTS)){
    dat_temp = data.frame(data.frame(component_data)[,c('truth','h','ts_length','ts_id')], 
                          ordering = data.frame(ORDER_BY_COMPONENTS)[,paste0(COMPONENTS[i])],
                          model = COMPONENTS[i])
    components_long = rbind(components_long, dat_temp)
  }
  train_data$model = COMPONENTS[train_data$class]
  train_data = merge(train_data, components_long[,c('h','ts_length','ts_id','model','ordering')],
                     by = c('h','ts_length','ts_id','model'), all.x = T, all.y = F)
  train_data$class = train_data$ordering
  train_data = subset(train_data, select=-c(ordering, model))
  #############  
}else{
  #############
  stop('Invalid fit_type')
  #############
}




input = subset(train_data, select=setdiff(names(train_data),c("ts_length","ts_id")))
df <- na.omit(input)
df$class <- as.factor(df$class)
  

if(param_type == 'multilogloss'){
  FILENAME = gsub('.RDS','_multilogloss.RDS',FILENAME)
  calibrated_params = read.csv(paste(my_path, "/best_params/best_param_", fit_num, ".csv", sep = ""))
  ## define parameters for LightGBM model
  params = list(objective = "multiclass",
                    metric = c("multi_logloss"),
                    num_class = length(unique(df$class)),
                    num_leaves = calibrated_params[which(calibrated_params$X=='num_leaves'),2],
                    learning_rate = calibrated_params[which(calibrated_params$X=='learning_rate'),2],
                    feature_fraction = calibrated_params[which(calibrated_params$X=='feature_fraction'),2],
                    max_depth = calibrated_params[which(calibrated_params$X=='max_depth'),2],
                    verbose = -1)
}else if(param_type == 'multierror'){
  FILENAME = gsub('.RDS','_multierror.RDS',FILENAME)
  calibrated_params = read.csv(paste(my_path,"/best_params_for_multierror/best_param_", fit_num, ".csv", sep = ""))
  ## define parameters for LightGBM model
  params = list(objective = "multiclass",
                metric = c("multi_error"),
                num_class = length(unique(df$class)),
                num_leaves = calibrated_params[which(calibrated_params$X=='num_leaves'),2],
                learning_rate = calibrated_params[which(calibrated_params$X=='learning_rate'),2],
                feature_fraction = calibrated_params[which(calibrated_params$X=='feature_fraction'),2],
                max_depth = calibrated_params[which(calibrated_params$X=='max_depth'),2],
                verbose = -1)
}else{
  stop('Invalid param_type')
}

## divide df into train and validate
df_train_ids <- sample(1:nrow(df), calibrated_params[which(calibrated_params$X=='prop_holdout'),2]*nrow(df), replace=F)
df_train <- df[df_train_ids,]
df_valid_id <- setdiff(1:nrow(df),df_train_ids)
df_valid <- df[df_valid_id,]

## convert the data to LightGBM dataset format: training data
lgbm_train <- lightgbm::lgb.Dataset(data = as.matrix(subset(df_train, select=setdiff(names(df_train),c("class","class_wt")))),
                                        label = as.integer(df_train$class)-1,
                                        weight = df_train$class_wt)

## convert the data to LightGBM dataset format: validation data
lgbm_valid <- lightgbm::lgb.Dataset(data = as.matrix(subset(df_valid, select=setdiff(names(df_valid),c("class","class_wt")))),
                                        label = as.integer(df_valid$class)-1,
                                        weight = df_valid$class_wt,
                                        reference = lgbm_train)


## train the LightGBM model
lgb_model <- lightgbm::lgb.train(params = params,
                                     data = lgbm_train,
                                     valids = list(train = lgbm_train,
                                                   valid = lgbm_valid),
                                     early_stopping_rounds = calibrated_params[which(calibrated_params$X=='early_stopping_rounds'),2],
                                     nrounds = 10000)

## save the fitted models
saveRDS(lgb_model, file = FILENAME)




