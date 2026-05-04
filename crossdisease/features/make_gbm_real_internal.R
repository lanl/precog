
################################
### Get Global SLURM Options ###
################################

options(warn = -1)

suppressMessages(library("optparse"))

option_list <- list(
  make_option("--total_run", type="character", default=""),
  make_option("--cur_run", type="character", default="")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options

total_run = as.numeric(opt$total_run) 
cur_run = as.numeric(opt$cur_run) 


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

my_path = '~/GitLab/frodo/features/'
data_path = '~/GitLab/frodo/raw_data/'
lut_path = '~/GitLab/frodo/lookup_tables/'
setwd(my_path)

## define paths
savetrainpath <- my_path

## define number of cores
ncores <- floor(.5*detectCores())

## define socket type
sockettype <- "PSOCK"

source(paste0(my_path, 'features_functions.R'))

## combination function for foreach when returning a list
## where each item of the list is meant to be combined via rbind()
comb <- function(x, ...) {
  mapply(rbind,x,...,SIMPLIFY=FALSE)
}

FILES_ALL = gsub('.RDS','',list.files(paste0(data_path)))
FILES_CUR = FILES_ALL[seq(cur_run,length(FILES_ALL), by = total_run)]


for(i in 1:length(FILES_CUR)){
  
  eval_key = FILES_CUR[i]
  test <- readRDS(paste0(data_path,"/",eval_key, ".RDS"))
  g_inds = 1:length(test)

  ## subset to only weekly
  CADENCE = unlist(lapply(test, FUN = function(x){x$ts_time_cadence}))
  g_inds = g_inds[CADENCE == 'weekly']
  
  ## set up cluster
  cl <- parallel::makeCluster(spec = ncores,
                              type = sockettype)
  setDefaultCluster(cl)
  registerDoParallel(cl)
  
  print(Sys.time())
  train_data <- foreach(g=g_inds, 
                        .errorhandling = "pass",
                        .verbose = T,
                        .packages = c('dplyr', 'tsfeatures', 'timeDate', 'lubridate','mgcv','zoo','forecast','MMWRweek'))%dopar%{
                          
                          if((!file.exists(paste0(my_path,"results/real_train_mat_",g,'_', eval_key, ".csv")))){
                            library(tsfeatures)
                            library(timeDate)
                            library(lubridate)
                            library(forecast)
                            ts_test = test[[g]]
                            h=4
                              
                            ### correct time series listed as counts but really scaled by population
                            if(ts_test$ts_scale == 'counts' & sum(ts_test$ts != round(ts_test$ts)) > 0){
                              ts_test$ts_scale = 'other'  
                            }
                            
                            ### Subset to valid forecast dates
                            x_inds = 11:unique(pmax(11,(length(ts_test$ts)-h)))
                            nonzero_inds = sort(unique(which(ts_test$ts>0)))
                            nonzero_thresh = ifelse(length(nonzero_inds)>=5, nonzero_inds[5], Inf)
                            x_inds = x_inds[x_inds>=nonzero_thresh]
                            
                            if(length(x_inds)>=1){

                              ################################
                              ### Get Time Series Features ###
                              ################################

                              h=4
                              res = apply(cbind(x_inds),1,
                                          FUN = get_features,
                                          ts = as.vector(ts_test$ts), 
                                          ts_dates = as.character(ts_test$ts_dates),
                                          ts_time_cadence = ts_test$ts_time_cadence,
                                          scaling_type = 'lastobs',
                                          smoothing_type = 'none', 
                                          smooth_truth = T, #just adds one more variable, truthsm
                                          info_packet = ts_test)
                              res_long = do.call('rbind',res)

                              
                              ############################
                              ### Get Disease Features ###
                              ############################
                              #disease characteristics
                              DISEASE = strsplit(eval_key, split = '_')[[1]]
                              DISEASE = DISEASE[1:(length(DISEASE)-1)]
                              DISEASE = paste(DISEASE, collapse = '_')
                              disease_lut = readxl::read_xlsx(paste0(lut_path, 'Disease_Mappings.xlsx'))
                              disease_lut = disease_lut[disease_lut$Disease == DISEASE, ]
                              disease_lut = disease_lut[1,]
                              res_long$transmission = disease_lut$Transmission
                              res_long$viral_bacterial_fungal = as.character(disease_lut[,"Viral/Bacterial/Fungal"])
                              res_long$disease = DISEASE
                              
                              ##########################
                              ### Get Other Features ###
                              ##########################
                              res_long$ts_time_cadence = ts_test$ts_time_cadence
                              res_long$ts_scale = ts_test$ts_scale
                              res_long$ts_measurement_type = ts_test$ts_measurement_type
                              res_long$ts_measurement_type[res_long$ts_measurement_type == 'incidence'] = 'incident_cases'
                              
                              ########################################
                              ### Record-keeping, Not For Training ###
                              ########################################
                              if(is.null(ts_test$ts_geography)){
                                res_long$geography = NA #geography is missing, need to fix
                              }else{
                                res_long$geography = gsub(',','',ts_test$ts_geography)
                              }
                              res_long$disease_source = eval_key
                              A = strsplit(eval_key, split = '_')[[1]]
                              res_long$source = A[length(A)]
                              
                              write.csv(res_long, file = paste0(my_path,"results/real_train_mat_",g,'_', eval_key, ".csv"), quote = F, row.names = F)
                            }
                          }
                          xx <- 1
                          xx
                        }
  stopCluster(cl)
  print(train_data)
  print(Sys.time())
  print(paste0('Completed: ', FILES_CUR[i]))
}