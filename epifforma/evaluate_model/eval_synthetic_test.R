########################################################
########################################################
### epiFFORMA Predictions for Synthetic Test Dataset ###
########################################################
########################################################

################################
### Get Global SLURM Options ###
################################

options(warn = -1)

suppressMessages(library("optparse"))

option_list <- list(
  make_option("--eval_key", type="character", default=""),
  make_option("--eval_type", type="character", default=""),
  make_option("--param_type", type="character", default="")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options

eval_key = as.character(opt$eval_key) 
eval_type = as.character(opt$eval_type) 
param_type = as.character(opt$param_type) 
print(paste0(eval_type, ' : ', eval_key, ' : ', param_type))


##############################################################
## STEP 0: PRELIMINARIES
##############################################################

gc() 
## load libraries
library(ggplot2)
library(data.table)
library(plyr)
library(gridExtra)
library(lubridate)
library(parallel)
library(doParallel)
library(grid)
theme_set(theme_bw())


#### set up computer for parallelization
## define number of cores
ncores <- floor(.5*detectCores())

## define socket type
sockettype <- "PSOCK"


## define paths
library(this.path)
my_path = this.path::here()
setwd(my_path)
savetrainpath <- paste0(my_path, '/../process_data/')
savevalidpath <- paste0(my_path, '/../uq/')
syntheticpath <- paste0(my_path, '/../raw_data/synthetic/output/') 
f2wpath <- paste0(my_path,"/../fit_model/features2weights/")
savepath <- paste0(my_path,"/evaluation/")
outputname = 'synthetic_test'

## source in function
source(paste0("../process_data/epi_functions.R"))


## define forecast horizon
h = 4


######################################################
### Read in Test Features and Component Models ###
######################################################

#read in test data
lf <- list.files(paste0(savevalidpath,"validation_data_packets/"))[grep("features_se_",list.files(paste0(savevalidpath,"validation_data_packets/")))]
train_data <- NULL
for(i in 1:length(lf)){
  train_data <- rbind(train_data,fread(paste0(paste0(savevalidpath,"validation_data_packets/"),lf[i])))
}
train_data = train_data[!duplicated(paste0(train_data$ts_length, '_', train_data$ts_id, '_', train_data$h)),]
train_data = subset(train_data, select = -c(class,class_wt)) 
FEATURE_LIST = setdiff(names(train_data),c('ts_length','ts_id'))

#read in component forecasts
comp_ex <- list.files(paste0(savevalidpath,"validation_data_packets/"))[grepl("component",list.files(paste0(savevalidpath,"validation_data_packets/"))) & 
                                                                          !grepl("components_minusk",list.files(paste0(savevalidpath,"validation_data_packets/")))]
component_data <- NULL
for(i in 1:length(comp_ex)){
  component_data <- rbind(component_data,fread(paste0(paste0(savevalidpath,"validation_data_packets/"),comp_ex[i])))
} 
component_data = component_data[!duplicated(paste0(component_data$ts_length, '_', component_data$ts_id, '_', component_data$h)),]
MODEL_LIST = setdiff(names(component_data),c('truth','h','ts_length','ts_id'))

train_data = merge(train_data, component_data, by = c('h', 'ts_length', 'ts_id'), all.x = T, all.y = T) 


##################################################################
## STEP 4: make predictions for a geography
##################################################################

## combination function for foreach when returning a list
## where each item of the list is meant to be combined via rbind()
comb <- function(x, ...) {
  mapply(rbind,x,...,SIMPLIFY=FALSE)
}


## read in fitted model
if(eval_type == 'order'){
  if(param_type == 'multilogloss'){
    mod_wt <- readRDS(paste0(f2wpath,"fitted_lgb_models_order_multilogloss.RDS"))
  }else if(param_type == 'multierror'){
    mod_wt <- readRDS(paste0(f2wpath,"fitted_lgb_models_order_multierror.RDS"))
  }else{
    stop('Invalid param_type')
  }
  components = data.frame(train_data)[,MODEL_LIST]
  components_reordered = t(as.matrix(apply(data.frame(train_data)[,MODEL_LIST],1,FUN = function(x){sort(unlist(x))})))
}else{
  stop('Invalid eval_type')
}


pred_wts_list <- list()
for(j in 1:length(mod_wt)){
  pred_wts_list[[j]] <- predict(mod_wt[[j]], newdata = as.matrix(data.frame(train_data)[,FEATURE_LIST]))
  print(j)
}

## average all the fits
pred_wts <- apply(simplify2array(pred_wts_list), 1:2, mean)
pred_wts_sd <- apply(simplify2array(pred_wts_list), 1:2, sd)

##truncation
pred_wts[pred_wts<0.01] = 0 #changed from 0.05 on 7/7/24
pred_wts = t(apply(pred_wts, 1, FUN = function(x){x/sum(x)}))




## set epifforma equal to the mean forecast
if(eval_type == 'order'){
  fcst_df <- data.frame(epifforma = rowSums(pred_wts*as.matrix(components_reordered)),
                        data.frame(train_data)[,c(MODEL_LIST,'truth',FEATURE_LIST,'ts_length','ts_id')])
  
  pred_wts_reordered = matrix(NA, ncol = ncol(pred_wts), nrow = nrow(pred_wts))
  pred_wts_sd_reordered = matrix(NA, ncol = ncol(pred_wts_sd), nrow = nrow(pred_wts_sd))
  for(j in 1:nrow(components)){
    pred_wts_reordered[j, ] = pred_wts[j,order(order(unlist(components[j,])) )] #double order inverts the ordering operaiton
    pred_wts_sd_reordered[j,] = pred_wts_sd[j,order(order(unlist(components[j,])) )] #double order inverts the ordering operaiton
  }
}else{
  fcst_df <- data.frame(epifforma = pred_wts*scaling_factor,
                        data.frame(train_data)[,c(MODEL_LIST,'truth',FEATURE_LIST,'ts_length','ts_id')])
  for( i in 1:length(MODEL_LIST) ){
    fcst_df[,MODEL_LIST[i]] = fcst_df[,MODEL_LIST[i]]*scaling_factor
  }
}

## melt it
fcst_df_melt <- melt(fcst_df, id.vars = c("h","truth",FEATURE_LIST,'ts_length','ts_id'))
fcst_df_melt$type <- "fcst"
names(fcst_df_melt)[names(fcst_df_melt) == 'variable'] = 'model'
names(fcst_df_melt)[names(fcst_df_melt) == 'value'] = 'fcst'

features = data.frame(train_data)[,c(FEATURE_LIST,'ts_length','ts_id')]


## prepare to leave: pred_wts
pred_wts <- data.frame(pred_wts_reordered)
names(pred_wts) <- paste0(names(data.frame(train_data)[,MODEL_LIST]),"_avg")
names(pred_wts) <- paste0(names(data.frame(train_data)[,MODEL_LIST]),"_avg")
pred_wts$ts_length = train_data$ts_length 
pred_wts$ts_id = train_data$ts_id 
pred_wts$h = train_data$h 

## prepare to leave: pred_wts_sd
pred_wts_sd <- data.frame(pred_wts_sd_reordered)
names(pred_wts_sd) <- paste0(names(data.frame(train_data)[,MODEL_LIST]),"_sd")
pred_wts_sd$ts_length = train_data$ts_length 
pred_wts_sd$ts_id = train_data$ts_id 
pred_wts_sd$h = train_data$h 

## prepare to leave: components
components = data.frame(train_data)[,c(MODEL_LIST,'h','ts_length','ts_id')]


## append results
bigplotoutput <- fcst_df_melt[,-1]#duplicated h
bigfeatures <- features
bigwts <- pred_wts
bigwtssd <- pred_wts_sd
bigcomponents <- components


### add rows for equal_wt
SUB = bigplotoutput[!is.na(bigplotoutput$h) & bigplotoutput$model != 'epifforma',]
SUB = SUB[SUB$model != 'mirror',]
SUB = SUB %>% dplyr::group_by(ts_length, ts_id, h) %>% dplyr::mutate(fcst_mean = mean(fcst, na.rm=T))
SUB = SUB[!duplicated(paste0(SUB$ts_length, '_', SUB$h, '_', SUB$ts_id)),]
SUB$fcst = SUB$fcst_mean
SUB$model = 'equal_wt'
bigplotoutput = rbind(bigplotoutput, SUB[,colnames(bigplotoutput)])



# library(dplyr)
library(dtwclust)
dat_temp = bigplotoutput[bigplotoutput$type == 'fcst',]

## make the accuracy and print it out
dat_temp$fcst <- pmax(1e-10, dat_temp$fcst)
accuracydf <- ddply(dat_temp,.(model),summarise,
                    n = length(fcst),
                    smape_no0 = mean( (abs(fcst - truth)/sum(0.5*(abs(fcst) + abs(truth))))[truth > 0]  ),
                    mae = mean(abs(fcst - truth)),
                    rmse = sqrt(mean(abs(fcst - truth)^2)))
accuracydf$smape_rank <- rank(accuracydf$smape_no0)
accuracydf$mae_rank <- rank(accuracydf$mae)
accuracydf$rmse_rank <- rank(accuracydf$rmse)
accuracydf$avg_rank <- (1/3)*(accuracydf$smape_rank + accuracydf$mae_rank + accuracydf$rmse_rank)
accuracydf$std_smape <- 1-(accuracydf$smape_no0 - min(accuracydf$smape_no0))/(diff(range(accuracydf$smape_no0)))
accuracydf$std_mae <- 1-(accuracydf$mae - min(accuracydf$mae))/(diff(range(accuracydf$mae)))
accuracydf$std_rmse <- 1-(accuracydf$rmse - min(accuracydf$rmse))/(diff(range(accuracydf$rmse)))
accuracydf$std_comb <- accuracydf$std_mae + accuracydf$std_rmse + accuracydf$std_smape
accuracydf <- accuracydf[order(accuracydf$mae, decreasing = T),]
print(accuracydf)
gc()                      


rm(dat_temp)


## combine output into a list
output_list <- list(plot_df = bigplotoutput,
                    accuracy_table = accuracydf,
                    features = bigfeatures,
                    weights = bigwts,
                    weights_sd = bigwtssd,
                    components = bigcomponents)

## save the output
saveRDS(object = output_list, file = paste0(savepath, outputname, "_",eval_type,'_',param_type,".RDS"))



