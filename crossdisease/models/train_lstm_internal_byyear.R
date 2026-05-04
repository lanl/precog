
################################
### Get Global SLURM Options ###
################################

options(warn = -1)

suppressMessages(library("optparse"))

option_list <- list(
  make_option("--year_cutoff", type="character", default=""),
  make_option("--alpha", type="character", default=""),
  make_option("--scale_ts", type="character", default="")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options

# my_path = as.character(opt$my_path) # I don't think we need this anymore.
alpha = as.numeric(opt$alpha) 
year_cutoff = as.numeric(opt$year_cutoff) 
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

################################
### Read in Training Dataset ###
################################

library(data.table)



if(alpha == -1){
  ALPHA_GRID = c(0.025, 0.1, 0.25, 0.75, 0.9, 0.975)
}else{
  ALPHA_GRID = c(alpha)
}



if(year_cutoff == -1){
  YEAR_RUNS = c(2010:2024)
}else{
  YEAR_RUNS = c(year_cutoff)
}



for(year_cutoff in YEAR_RUNS){
  
  for(alpha in ALPHA_GRID){
    

    if(scale_ts == 'no'){
      FILENAME = paste0(my_path,'/lstm/fit_',year_cutoff,'_', alpha, '_noscale.h5')
    }else{
      FILENAME = paste0(my_path,'/lstm/fit_',year_cutoff,'_', alpha,'.h5')
    }
    if(!file.exists(FILENAME)){
      
      ### Read in Embeddings 
      embedding_path = '~/GitLab/frodo/features/embeddings_lstm/'
      dat = data.frame(data.table::fread(paste0(embedding_path,"/real_train_mat.csv")))
      
      
      dates_embed = dat$date
      dates_embed = as.Date(dates_embed,format = '%Y-%m-%d')
      YEAR = lubridate::year(dates_embed)
      dat = dat[YEAR < year_cutoff,]
      
      
      X = as.matrix(dat[,paste0('X_',1:52)])
      y = as.matrix(dat[,'y1'])
      
      for(j in 1:ncol(X)){
        X[,j] = pmax(0,X[,j])
      }
      
      for(j in 1:ncol(y)){
        y[,j] = pmax(0,y[,j])
      }
      # 
      if(scale_ts == 'yes'){
        scale_factor = X[,ncol(X)]
        scale_factor[scale_factor==0]=1
        
        MAX = 100
        for(j in 1:ncol(X)){
          X[,j] = X[,j]/scale_factor
          X[,j] = pmin(MAX,X[,j])
        }
        
        MAX = 3
        for(j in 1:ncol(y)){
          y[,j] = y[,j]/scale_factor
          y[,j] = pmin(MAX,y[,j])
        }
      }
      
      X[is.na(X)] = -1
      y[is.na(y)] = -1
      
      #########################
      ### Fitting Functions ###
      #########################
      library(reticulate)
      library(keras)
      library(tensorflow)
      #virtualenv_create(python = "/projects/opt/centos8/x86_64/miniconda3/py312_24.11.1/bin/python3.12", force = TRUE)
      use_virtualenv('~/.virtualenvs/r-reticulate', required=TRUE)
      #tensorflow::install_tensorflow(virtualenv = 'virtualenv', envname = "r-reticulate")
      #source ~/.virtualenvs/r-reticulate/bin/activate
      library(keras3)
      
      model <- keras_model_sequential() %>%
        layer_masking(mask_value = -1, input_shape = c(ncol(X), 1)) %>% 
        layer_lstm(units = 64, return_sequences = TRUE) %>%
        layer_dropout(rate = 0.2) %>%
        layer_lstm(units = 32, return_sequences = FALSE) %>%
        layer_dropout(rate = 0.2) %>%
        layer_dense(units = 1)
      
      
      quantile_loss <- function(q) {
        function(y_true, y_pred) {
          e <- y_true - y_pred
          k_mean(k_maximum(q * e, (q - 1) * e), axis = -1)
        }
      }
      
      model %>% compile(
        #loss = 'mean_squared_error',
        loss = quantile_loss(alpha),
        optimizer = 'adam'
      )
      
      early_stop <- callback_early_stopping(
        monitor = "val_loss",     # what to monitor
        patience = 1,             # number of epochs with no improvement to wait
        restore_best_weights = TRUE  # revert to best model
      )
      
      # Fit
      model %>% fit(
        X, y,
        epochs = 3,
        batch_size = 32,
        validation_split = 0.1,
        callbacks = list(early_stop)
      )
      
      
      
      #######################
      ### Perform Fitting ###
      #######################
    
      save_model_hdf5(model,FILENAME)
      gc()
      k_clear_session()
    }
        
  }
}
