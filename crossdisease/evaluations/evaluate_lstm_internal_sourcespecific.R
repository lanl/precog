
################################
### Get Global SLURM Options ###
################################

options(warn = -1)

suppressMessages(library("optparse"))

option_list <- list(
  make_option("--alpha", type="character", default=""),
  make_option("--total_run", type="character", default=""),
  make_option("--cur_run", type="character", default=""),
  make_option("--include_disease_source", type="character", default=""),
  make_option("--warm_start", type="character", default=""),
  make_option("--scale_ts", type="character", default="")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options

total_run = as.numeric(opt$total_run) 
cur_run = as.numeric(opt$cur_run) 
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

my_path = '~/GitLab/frodo/evaluations/'
feature_path = '~/GitLab/frodo/features/embeddings_lstm/'
model_path = '~/GitLab/frodo/models/lstm_sourcespecific/'


setwd(my_path)

### Get Files to Evaluate ###
FILES_ALL = list.files(paste0(feature_path))
FILES_ALL = FILES_ALL[!grepl('real_train_mat.csv',FILES_ALL)]
FILES_ALL = FILES_ALL[grepl(include_disease_source, FILES_ALL)]
FILES_CUR = FILES_ALL[seq(cur_run,length(FILES_ALL), by = total_run)]







library(reticulate)
library(keras)
library(tensorflow)
use_virtualenv('~/.virtualenvs/r-reticulate', required=TRUE)
library(keras3)





if(alpha == -1){
  ALPHA_GRID = c(0.025, 0.1, 0.25, 0.75, 0.9, 0.975)
}else{
  ALPHA_GRID = c(alpha)
}

for(alpha in ALPHA_GRID){
  
  
  lut_path = '~/GitLab/frodo/lookup_tables/'
  DISEASE_SOURCES = readxl::read_xlsx(paste0(lut_path, 'Evaluations_ML.xlsx'))
  DISEASE_SOURCES = data.frame(DISEASE_SOURCES)
  DISEASE_SOURCES$FILES = gsub('.RDS','',DISEASE_SOURCES$FILES)
  DISEASE_SOURCES = DISEASE_SOURCES[DISEASE_SOURCES$FILES == include_disease_source, ]
  YEAR_RUNS = c(DISEASE_SOURCES$ML_start:DISEASE_SOURCES$ML_end)
  
  
  
  quantile_loss <- function(q) {
    function(y_true, y_pred) {
      e <- y_true - y_pred
      k_mean(k_maximum(q * e, (q - 1) * e), axis = -1)
    }
  }
  
  if(scale_ts == 'no'){
      mod = replicate(length(YEAR_RUNS), list(NULL))
      for(y in 1:length(YEAR_RUNS)){
        test <- keras_model_sequential() %>%
          layer_masking(mask_value = -1, input_shape = c(52, 1)) %>%
          layer_lstm(units = 64, return_sequences = TRUE) %>%
          layer_dropout(rate = 0.2) %>%
          layer_lstm(units = 32, return_sequences = FALSE) %>%
          layer_dropout(rate = 0.2) %>%
          layer_dense(units = 1)
        test %>% compile(
          loss = quantile_loss(alpha),
          optimizer = 'adam'
        )
        test %>% load_model_weights_hdf5(paste0(model_path,"fit_",YEAR_RUNS[y],'_',include_disease_source,"_",alpha,"_cold_noscale.h5"))
        mod[[y]] <- test
      }
  }else{
      mod = replicate(length(YEAR_RUNS), list(NULL))
      for(y in 1:length(YEAR_RUNS)){
        test <- keras_model_sequential() %>%
          layer_masking(mask_value = -1, input_shape = c(52, 1)) %>%
          layer_lstm(units = 64, return_sequences = TRUE) %>%
          layer_dropout(rate = 0.2) %>%
          layer_lstm(units = 32, return_sequences = FALSE) %>%
          layer_dropout(rate = 0.2) %>%
          layer_dense(units = 1)
        test %>% compile(
          loss = quantile_loss(alpha),
          optimizer = 'adam'
        )
        test %>% load_model_weights_hdf5(paste0(model_path,"fit_",YEAR_RUNS[y],'_',include_disease_source,"_",alpha,"_cold.h5"))
        mod[[y]] <- test
      }
  }
  
  
  
  for(i in 1:length(FILES_CUR)){
    
    eval_key = gsub('.csv','',FILES_CUR[i])
    eval_key = gsub('embed_','',eval_key)
    
    
    if(scale_ts == 'no'){
      if(warm_start == 'yes'){
        FILENAME = paste0(my_path,"lstm_noscale/sourcespecific_warm/real_eval_mat_",eval_key,'_', alpha, ".csv")
      }else{
        FILENAME = paste0(my_path,"lstm_noscale/sourcespecific_cold/real_eval_mat_",eval_key,'_', alpha, ".csv")
      }
    }else{
      if(warm_start == 'yes'){
        FILENAME = paste0(my_path,"lstm/sourcespecific_warm/real_eval_mat_",eval_key,'_', alpha, ".csv")
      }else{
        FILENAME = paste0(my_path,"lstm/sourcespecific_cold/real_eval_mat_",eval_key,'_', alpha, ".csv")
      }
    }
    
    
    if(!file.exists(FILENAME)){
      dat <- read.csv(paste0(feature_path,FILES_CUR[i]))
      
      X_real = as.matrix(dat[,paste0('X_',1:52)])
      y_real = as.matrix(dat[,paste0('y',1:4)])
      
      
      for(j in 1:ncol(X_real)){
        X_real[,j] = pmax(0,X_real[,j])
      }
      
      X_real[is.na(X_real)] = -1
      
      lastobs_real = dat$obs
      dates_real = as.character(dat$date)
      
      
      dat_long = data.frame(X_real)
      colnames(dat_long) = paste0('X_',1:52)
      dat_long$date = dates_real
      dat_long$truth_1 = y_real[,1]
      dat_long$truth_2 = y_real[,2]
      dat_long$truth_3 = y_real[,3]
      dat_long$truth_4 = y_real[,4]
      
      dat_long$year = substr(dat_long$date,1,4)
      dat_long$obs = lastobs_real
      dat_long = dat_long[!is.na(dat_long$y),]
      TO_RUN = sort(intersect(unique(as.numeric(dat_long$year)), YEAR_RUNS))
      
      FEATURE_LIST = paste0('X_',1:52)
      
      PMIN = 0
      if(scale_ts == 'no'){
        PMAX = Inf
      }else{
        PMAX = 3
      }
      if(length(TO_RUN)>0){
        out_data = NULL
        dat_long = dat_long[dat_long$year >= min(TO_RUN),]
        for(y in 1:length(TO_RUN)){
          SUBSET = dat_long[dat_long$year == TO_RUN[y],]
          y2 = which(YEAR_RUNS == TO_RUN[y])
          if(nrow(SUBSET)>1){
            if(scale_ts == 'yes'){
              MAX = 100
            }else{
              MAX = Inf
            }        
            ### h=1
            X_real_sub = as.matrix(data.frame(SUBSET)[,FEATURE_LIST])
            if(scale_ts == 'yes'){
              scale_factor = X_real_sub[,ncol(X_real_sub)]
              scale_factor[scale_factor==0]=1
            }else{
              scale_factor = 1
            }
            for(j in 1:ncol(X_real_sub)){
              X_real_sub[,j] = X_real_sub[,j]/scale_factor
              X_real_sub[,j] = pmin(MAX,X_real_sub[,j])
            }
            pred_list <- list()
            for(j in 1:length(mod[[y2]])){
              pred_list[[j]] <- mod[[y2]] %>% predict(X_real_sub)
              print(j)
            }
            SUBSET$ratio1 = apply(simplify2array(pred_list), 1:2, mean)[,1]
            SUBSET$ratio1 = pmin(pmax(PMIN, SUBSET$ratio1),PMAX)
            SUBSET$fcst1 = SUBSET$ratio1 * scale_factor
            
            ### h=2
            X_real_sub = as.matrix(data.frame(SUBSET)[,FEATURE_LIST])
            X_real_sub = X_real_sub[,-1]
            X_real_sub = cbind(X_real_sub, SUBSET$fcst1)
            if(scale_ts == 'yes'){
              scale_factor = X_real_sub[,ncol(X_real_sub)]
              scale_factor[scale_factor==0]=1
            }else{
              scale_factor = 1
            }
            for(j in 1:ncol(X_real_sub)){
              X_real_sub[,j] = X_real_sub[,j]/scale_factor
              X_real_sub[,j] = pmin(MAX,X_real_sub[,j])
            }
            pred_list <- list()
            for(j in 1:length(mod[[y2]])){
              pred_list[[j]] <- mod[[y2]] %>% predict(X_real_sub)
              print(j)
            }
            SUBSET$ratio2 = apply(simplify2array(pred_list), 1:2, mean)[,1]
            SUBSET$ratio2 = pmin(pmax(PMIN, SUBSET$ratio2),PMAX)
            SUBSET$fcst2 = SUBSET$ratio2 * scale_factor
            
            ### h=3
            X_real_sub = as.matrix(data.frame(SUBSET)[,FEATURE_LIST])
            X_real_sub = X_real_sub[,-(1:2)]
            X_real_sub = cbind(X_real_sub, SUBSET$fcst1, SUBSET$fcst2)
            if(scale_ts == 'yes'){
              scale_factor = X_real_sub[,ncol(X_real_sub)]
              scale_factor[scale_factor==0]=1
            }else{
              scale_factor = 1
            }
            for(j in 1:ncol(X_real_sub)){
              X_real_sub[,j] = X_real_sub[,j]/scale_factor
              X_real_sub[,j] = pmin(MAX,X_real_sub[,j])
            }
            pred_list <- list()
            for(j in 1:length(mod[[y2]])){
              pred_list[[j]] <- mod[[y2]] %>% predict(X_real_sub)
              print(j)
            }
            SUBSET$ratio3 = apply(simplify2array(pred_list), 1:2, mean)[,1]
            SUBSET$ratio3 = pmin(pmax(PMIN, SUBSET$ratio3),PMAX)
            SUBSET$fcst3 = SUBSET$ratio3 * scale_factor
            
            ### h=4
            X_real_sub = as.matrix(data.frame(SUBSET)[,FEATURE_LIST])
            X_real_sub = X_real_sub[,-(1:3)]
            X_real_sub = cbind(X_real_sub, SUBSET$fcst1, SUBSET$fcst2, SUBSET$fcst3)
            if(scale_ts == 'yes'){
              scale_factor = X_real_sub[,ncol(X_real_sub)]
              scale_factor[scale_factor==0]=1
            }else{
              scale_factor = 1
            }
            for(j in 1:ncol(X_real_sub)){
              X_real_sub[,j] = X_real_sub[,j]/scale_factor
              X_real_sub[,j] = pmin(MAX,X_real_sub[,j])
            }
            pred_list <- list()
            for(j in 1:length(mod[[y2]])){
              pred_list[[j]] <- mod[[y2]] %>% predict(X_real_sub)
              print(j)
            }
            SUBSET$ratio4 = apply(simplify2array(pred_list), 1:2, mean)[,1]
            SUBSET$ratio4 = pmin(pmax(PMIN, SUBSET$ratio4),PMAX)
            SUBSET$fcst4 = SUBSET$ratio4 * scale_factor
            
            
            SUBSET$fcst1[SUBSET$obs == 0] = 0
            SUBSET$fcst2[SUBSET$obs == 0] = 0
            SUBSET$fcst3[SUBSET$obs == 0] = 0
            SUBSET$fcst4[SUBSET$obs == 0] = 0
            dat_temp = rbind(data.frame(SUBSET[,c('obs', 'date','year')], truth = SUBSET$truth_1, fcst = SUBSET$fcst1, h=1),
                             data.frame(SUBSET[,c('obs', 'date','year')], truth = SUBSET$truth_2, fcst = SUBSET$fcst2, h=2),
                             data.frame(SUBSET[,c('obs', 'date','year')], truth = SUBSET$truth_3, fcst = SUBSET$fcst3, h=3),
                             data.frame(SUBSET[,c('obs', 'date','year')], truth = SUBSET$truth_4, fcst = SUBSET$fcst4, h=4))
            
            
            out_data = rbind(out_data,dat_temp)
          }
        }
      
  
        
      write.csv(out_data, file = FILENAME, quote = F, row.names = F)
      gc()
      }
    }
    #print(Sys.time())
    print(paste0('Completed: ', FILES_CUR[i]))
  }

}