###########################################
###########################################
### Calibrate epiFFORMA Hyperparameters ###
###########################################
###########################################

################################
### Get Global SLURM Options ###
################################

options(warn = -1)

suppressMessages(library("optparse"))

fit_type = 'order'
fit_num = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))

if(is.na(fit_num)){
	fit_type = 'order'
	fit_num = 1
}

print(paste("fit_type is", fit_type))
print(paste("fit_num is", fit_num))

######################
### Body of Script ###
######################


## load libraries
library(ggplot2)
library(ParBayesianOptimization)
library(data.table)
library(plyr)
library(gridExtra)
library(lubridate)
library(parallel)
library(doParallel)
library(grid)
library(plotly)
library(GGally)
library(parallel)
library(doParallel)


num_cores = 20
print(paste("number of cores is", num_cores))

## define paths
library(this.path)
my_path = this.path::here()
setwd(my_path)
savetrainpath <- paste0(my_path, '/../process_data/')
syntheticpath <- paste0(my_path, '/../raw_data/synthetic/output/') 
figsavepath <- paste0(savetrainpath,"figs/")
savemodelpath <- paste0(my_path, '/features2weights/') 

## source in function
source(paste0("../process_data/epi_functions.R"))

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

calibrate_lgbm_wt <- function(input){
  
  ## unpack the training data into the features and class labels
  df       <- na.omit(input)
  df$class <- as.factor(df$class)
  df_train <- df
  lgbm_train <- lightgbm::lgb.Dataset(data = as.matrix(subset(df_train, select=setdiff(names(df_train),c("class","class_wt")))),
                                        label = as.integer(df_train$class)-1,
                                        weight = df_train$class_wt)
  
  # Define a function (based on this data) to calibrate
  fit_model_wrapper = function(num_leaves, #1 < nl <= 131072
                               learning_rate, # lr > 0
                               feature_fraction, #0<= ff <= 1 
                               max_depth,
                               early_stopping_rounds,
                               prop_holdout){ # -1 means no limit
                               #drop_rate){ #0<= dr <= 1

    ## train the LightGBM model
    # Since the prop_holdout is essentially a CV-type criteria, we'll do a sort of hacky version
    # of CV here
    num_cvs = 1
    scores  = c()
    for(cv in 1:num_cvs){

      ## divide df into train and validate
      df_train_ids <- sample(1:nrow(df), prop_holdout*nrow(df), replace=F) 
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
      
      
      ## define parameters for LightGBM model
      params = list(objective = "multiclass",
                    metric = c("multi_error"),
                    num_class = length(unique(df$class)),
                    num_leaves = num_leaves, 
                    learning_rate = learning_rate,
                    feature_fraction = feature_fraction,
                    max_depth = max_depth,
                    #drop_rate = drop_rate,
                    verbose = -1)
      
      ## train the LightGBM model
      lgb_model <- lightgbm::lgb.train(params = params,
                                       data = lgbm_train,
                                       valids = list(train = lgbm_train,
                                                     valid = lgbm_valid),
                                       early_stopping_rounds = early_stopping_rounds,
                                       nrounds = 1000)
      
      min(unlist(lgb_model$record_evals$valid$multi_logloss$eval))
      
      # Find a way to calculation an optimization criterion and store it in scores
      scores = c(scores, min(unlist(lgb_model$record_evals$valid$multi_error$eval)))
    }

    return(list(Score = (1 - mean(scores)),
                Pred  = 0) )
  } 
  
  
  # Perform the Bayesian model calibration.
  param_bounds = list(
    num_leaves = c(1000L, 2921L), 
    learning_rate = c(0.01, 0.3),
    feature_fraction = c(0.7, 1),
    max_depth = c(10L, 250L),
    early_stopping_rounds = c(50L, 200L),
    prop_holdout = c(0.6, 0.95)
    #drop_rate = c(0.01, 0.5)
  )


  print("registered the parallelization, entering the optimization function")
  opt_model = bayesOpt(FUN = fit_model_wrapper,
                                                          bounds = param_bounds,
                                                          initPoints = 30,
							  iters.n = 10,
                                                          iters.k = 5,
							  parallel = F,
							  verbose = 2)

  ## get outta here
  return(opt_model)
}

print("entering the calibrate lgbm protocol")
optObj = calibrate_lgbm_wt(input = subset(train_data, select=setdiff(names(train_data),c("ts_length","ts_id"))))
print(getBestPars(optObj))
write.csv( unlist(getBestPars(optObj)) , file = paste('best_params_for_multierror/best_param_', fit_num, '.csv', sep = ""))


print("Note that no model is saved here!  Just updated the optimal parameters.")


