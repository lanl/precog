
################################
### Get Global SLURM Options ###
################################

options(warn = -1)

suppressMessages(library("optparse"))

option_list <- list(
  make_option("--year_cutoff", type="character", default=""),
  make_option("--alpha", type="character", default=""),
  make_option("--scale_ts", type="character", default=""),
  make_option("--include_mode", type="character", default="")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options

# my_path = as.character(opt$my_path) # I don't think we need this anymore.
alpha = as.numeric(opt$alpha) 
scale_ts = as.character(opt$scale_ts)
include_mode = as.character(opt$include_mode)
year_cutoff = as.numeric(opt$year_cutoff) 

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


if(year_cutoff == -1){
  lut_path = '~/GitLab/frodo/lookup_tables/'
  MODES = readxl::read_xlsx(paste0(lut_path, 'Modes_ML.xlsx'))
  MODES = data.frame(MODES)
  MODES = MODES[MODES$FILES == include_mode, ]
  YEAR_RUNS = c(MODES$ML_start:MODES$ML_end)
}else{
  YEAR_RUNS = c(year_cutoff)
}




if(alpha == -1){
  ALPHA_GRID = c(0.025, 0.1, 0.25, 0.75, 0.9, 0.975)
}else{
  ALPHA_GRID = c(alpha)
}

for(year_cutoff in YEAR_RUNS){
  for(alpha in ALPHA_GRID){
    
    if(scale_ts == 'no'){
      FILENAME = paste0(my_path,'/lstm_modespecific/fit_',year_cutoff,'_', alpha, '_', include_mode, '_noscale.h5')
    }else{
      FILENAME = paste0(my_path,'/lstm_modespecific/fit_',year_cutoff,'_', alpha,'_', include_mode,'.h5')
    }
    if(!file.exists(FILENAME)){
      ### Read in Embeddings 
      embedding_path = '~/GitLab/frodo/features/embeddings_lstm/'
      dat = data.frame(data.table::fread(paste0(embedding_path,"/real_train_mat.csv")))
      
      
      ### Merge in mode of transmission
      FILES_TEMP = unique(dat$FILES)
      SPLITS = strsplit(FILES_TEMP, split = '_')
      DISEASES = rep(NA, length(FILES_TEMP))
      for(i in 1:length(SPLITS)){
        if(length(SPLITS[[i]])==3){
          DISEASES[i] = SPLITS[[i]][1]
        }else if(length(SPLITS[[i]])==4){
          DISEASES[i] = paste0(SPLITS[[i]][1],'_',SPLITS[[i]][2])
        }else if(length(SPLITS[[i]])==5){
          DISEASES[i] = paste0(SPLITS[[i]][1],'_',SPLITS[[i]][2],'_',SPLITS[[i]][3])
        }else{
          DISEASES[i] = paste0(SPLITS[[i]][1],'_',SPLITS[[i]][2],'_',SPLITS[[i]][3],'_',SPLITS[[i]][4])
        }
      }
      lut_path = '~/GitLab/frodo/lookup_tables/'
      disease_lut = readxl::read_xlsx(paste0(lut_path, 'Disease_Mappings.xlsx'))
      TO_MERGE = data.frame(FILES = FILES_TEMP, Disease = DISEASES)
      TO_MERGE = merge(TO_MERGE, disease_lut[,c('Disease', 'Transmission')], by = 'Disease')
      dat = merge(dat, TO_MERGE[,c('FILES', 'Transmission')], all.x = T, all.y = F)
      
      ### Subset by transmission mode
      if(include_mode == 'respiratory'){
        dat = dat[dat$Transmission == 'Respiratory_Secretions',]
      }else if(include_mode == 'fecaloral'){
        dat = dat[dat$Transmission == 'Fecal-Oral' | dat$Transmission == 'Water' | dat$Transmission == 'Fecal-Oral/Bodily Fluids',]
      }else if(include_mode == 'sexual'){
        dat = dat[dat$Transmission == 'Sexual',]
      }else if(include_mode == 'vectorborne'){
        dat = dat[dat$Transmission == 'Vectorborne',]
      }else{
        stop('include_trans not implemented')
      }
      
      
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
    }  
  }
  
}
