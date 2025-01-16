######################################################
######################################################
### epiFFORMA Forecasting Minimal Working Example  ###
######################################################
######################################################

### This script provides a minimal working example for how to make new predictions
### using a pre-trained epiFFORMA point forecast model and predictive interval
### widths model.


library(this.path)
setwd(this.path::here()) 

#################################
### Step 0: Define file paths ###
#################################

training_path = './process_data/'
model_path = './fit_model/'
uq_path = './uq/'
source(paste0(training_path,"epi_functions.R"))
source(paste0(training_path,'/../SLURMarray.r'))

######################################
### Step 1: Simulate a time series ###
######################################
x=seq(0.1,10,2*3.14/52)
set.seed(10)
ts = sin(x) + rnorm(length(x),0,0.1)
ts_scaled = pmax(0,round(100*(ts+1)))
info_packet = list(ts = ts_scaled,
                   ts_time_cadence = 'weekly',
                   ts_scale = 'counts')
plot(ts_scaled, xlab = 'Days', ylab = 'Cases')

### Expected value for next 4 observations
offset = 2*3.14/52
x2=seq(x[length(x)]+offset, x[length(x)]+4*offset, offset)
ts_scaled_future = 100*(sin(x2)+1)
plot(x, ts_scaled, xlab = 'Days', ylab = 'Cases', xlim = c(0,10.5))
lines(x2,ts_scaled_future, xlab = 'Days', ylab = 'Cases', col = 'red')


################################################
### Step 2: Estimate features and components ###
################################################

### Read in embedding matrices
embed_mat_X <- data.table::fread(file=paste0(training_path,"embed_mat/embed_mat_X.csv"))
embed_mat_y <- data.table::fread(file=paste0(training_path,"embed_mat/embed_mat_y.csv"))

embed_mat_X_deriv <- data.table::fread(file=paste0(training_path,"embed_mat/embed_mat_X_deriv.csv"))
embed_mat_y_deriv <- data.table::fread(file=paste0(training_path,"embed_mat/embed_mat_y_deriv.csv"))

## handle outliers
info_packet$ts = as.numeric(handle_outliers(info_packet))

## distribute zeros
info_packet$ts = as.numeric(distribute_zeros(info_packet))

## compute the features
features <- make_features(info_packet, h=4)

## compute the components
components <- make_components(info_packet, h=4)

## compute the components
component_intervals <- make_component_intervals(info_packet, h=4)


############################################################
### Step 3: Generate epiFFORMA predictions and intervals ###
############################################################


### Get Point Estimates
feature2wt <- readRDS(paste0(model_path,"/features2weights/fitted_lgb_models_order_multierror.RDS"))
pred_wts_list <- list()
for(j in 1:length(feature2wt)){
  pred_wts_list[[j]] <- predict(feature2wt[[j]], newdata = as.matrix(features))
}
pred_wts <- apply(simplify2array(pred_wts_list), 1:2, mean)
pred_wts_sd <- apply(simplify2array(pred_wts_list), 1:2, sd)
pred_wts[pred_wts<0.01] = 0 
pred_wts = t(apply(pred_wts, 1, FUN = function(x){x/sum(x)}))
components_reordered = NULL
for(j in 1:nrow(components)){
  components_reordered = rbind(components_reordered, sort(unlist(components[j,])))
}
fcst_df <- data.frame(h = 1:4,epifforma = rowSums(pred_wts*components_reordered))


### Get Interval Widths
component_intervals$fcst_mult = 4*fcst_df$epifforma
features2 = features
features2$weight_sd = apply(pred_wts_sd,1,sd)
features2$component_sd = apply(components,1,sd)/fcst_df$epifforma
features2$component_sd[features2$component_sd>1 | is.infinite(features2$component_sd)]=1

mod_interval = readRDS(paste0(uq_path,"/features2weights/interval2.RDS"))
pred_wts_list <- list()
for(j in 1:length(mod_interval)){
  pred_wts_list[[j]] <- predict(mod_interval[[j]], newdata = as.matrix(features2))
}
pred_wts <- apply(simplify2array(pred_wts_list), 1:2, mean)
pred_wts[pred_wts<0.01] = 0 
pred_wts = t(apply(pred_wts, 1, FUN = function(x){x/sum(x)}))
TRUNC = 50
fcst_df$width = rowSums(pred_wts*as.matrix(component_intervals))
fcst_df$width[fcst_df$width>(2*TRUNC*fcst_df$epifforma) & !is.na(fcst_df$width)] = 2*TRUNC*fcst_df$epifforma[fcst_df$width>(2*TRUNC*fcst_df$epifforma) & !is.na(fcst_df$width)]
fcst_df$width = pmax(1e-9,fcst_df$width) 



###########################################
### Plot the Forecast and 95% Intervals ###
###########################################
plot(x, ts_scaled, xlab = 'Days', ylab = 'Cases', xlim = c(0,10.5))
lines(x2,ts_scaled_future, xlab = 'Days', ylab = 'Cases', col = 'red')
lines(x2,fcst_df$epifforma, xlab = 'Days', ylab = 'Cases', col = 'green')
lines(x2,fcst_df$epifforma - 0.5*fcst_df$width, xlab = 'Days', ylab = 'Cases', col = 'green')
lines(x2,fcst_df$epifforma + 0.5*fcst_df$width, xlab = 'Days', ylab = 'Cases', col = 'green')









