# Create the real-data embeddings used to do the experiments in the supplement
## Author: LJ Beesley
## Date: January 2025

# Note to researcher:
# The raw data to create these embeddings are from public repositories that
# are cited in the main paper.  These data are prohibited from being provided
# on the Github repository.

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

data_path = '/path/to/rds/datafiles'
output_path = '/embeddings_gam_real/'


FILES_ALL = gsub('.RDS','',list.files(paste0(data_path)))
FILES_ALL = FILES_ALL[!grepl('cimmid',FILES_ALL)]
FILES_ALL = FILES_ALL[!grepl('ginkgo',FILES_ALL)]
FILES_ALL = FILES_ALL[!grepl('COVID',FILES_ALL)]
FILES_ALL = FILES_ALL[grepl('Chikungunya_deSouza',FILES_ALL)|
                        grepl('Dengue_opendengue',FILES_ALL)|
                        grepl('Influenza_usflunet',FILES_ALL)|
                        grepl('Influenza_ushhs',FILES_ALL)|
                        grepl('Mpox_who',FILES_ALL)]
FILES_CUR = FILES_ALL[seq(cur_run,length(FILES_ALL), by = total_run)]

####################
### Prepare sMOA ###
####################
k <- 5
h <- 4
library(mgcv)
library(collapse)


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
                      .packages = c('dplyr', 'tsfeatures', 'timeDate', 'lubridate','mgcv','zoo','forecast','collapse'))%dopar%{
                        library(mgcv)
                        library(collapse)
                        eval_key = FILES_CUR[i]
                        list_of_lists <- readRDS(paste0(data_path,"/",eval_key, ".RDS"))
                        NUM_TO_RUN = length(list_of_lists)
                        for(j in 1:NUM_TO_RUN){
                          ts_test = list_of_lists[[j]]
                          fcst_indices <- 7:(length(ts_test$ts)-h)
                          X = NULL
                          y = NULL
                          lastobs = NULL
                          if(!file.exists(paste0(output_path,"embed_",eval_key, '_',j,"_X.csv"))){
                            if((max(ts_test$ts) > 10 & ts_test$ts_scale == 'counts') | (max(ts_test$ts)>1e-3 & ts_test$ts_scale == 'proportion')){
                              for(l in fcst_indices){
                                
                                ## get info_packet and trim
                                info_packet <- ts_test
                                info_packet$ts <- info_packet$ts[1:l]
                                data_till_now = data.frame(value = info_packet$ts, t = 1:l)
                                if(nrow(data_till_now)>50){
                                  data_till_now = data_till_now[(nrow(data_till_now)-50):nrow(data_till_now),]
                                }
                                
                                #### light smoothing and differencing and get last k 
                                data_till_now_smoothed      <- gam(value~ s(t,k=pmax(round(nrow(data_till_now)/2))),data=data_till_now)$fitted.values
                                
                                #### could use gam smoother 
                                to_match_in_moa             <- tail((data_till_now_smoothed),k+1)
                                X_temp = to_match_in_moa
                                y_temp = ts_test$ts[(l+1):(l+h)]

                                X <- rbind(X, X_temp)
                                y <- rbind(y, y_temp)
                                lastobs = c(lastobs, ts_test$ts[l])
                                print(l)
                              }
                              write.csv(lastobs, file = paste0(output_path,"embed_",eval_key, '_',j,"_lastobs.csv"), quote = F, row.names = F)
                              write.csv(X, file = paste0(output_path,"embed_",eval_key, '_',j,"_X.csv"), quote = F, row.names = F)
                              write.csv(y, file = paste0(output_path,"embed_",eval_key, '_',j,"_y.csv"), quote = F, row.names = F)
                            }
                          }
                        }
                        xx <- 1
                        xx
                      }
stopCluster(cl)
