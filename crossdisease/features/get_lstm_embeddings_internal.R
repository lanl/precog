
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

data_path = '~/GitLab/frodo/raw_data/'
my_path = '~/GitLab/frodo/features/'


FILES_ALL = gsub('.RDS','',list.files(paste0(data_path)))
FILES_CUR = FILES_ALL[seq(cur_run,length(FILES_ALL), by = total_run)]


####################
### Prepare sMOA ###
####################
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
                        g_inds = 1:length(list_of_lists)
                        
                        ## subset to only weekly
                        CADENCE = unlist(lapply(list_of_lists, FUN = function(x){x$ts_time_cadence}))
                        g_inds = g_inds[CADENCE == 'weekly']
                        
                        for(j in g_inds){
                          ts_test = list_of_lists[[j]]
                          fcst_indices <- 11:(length(ts_test$ts)-h)
                          X = NULL
                          # X_plus1 = NULL
                          y = NULL
                          lastobs = NULL
                          dates = NULL
                          if(ts_test$ts_scale == 'counts' & sum(ts_test$ts != round(ts_test$ts)) > 0){
                            ts_test$ts_scale = 'other'  
                          }
                          
                          if(!file.exists(paste0(my_path,'embeddings_lstm/',"embed_",eval_key, '_',j,".csv")) & ts_test$ts_time_cadence == 'weekly'){
                            #if((max(ts_test$ts) > 10 & ts_test$ts_scale == 'counts') | (max(ts_test$ts)>1e-3 & ts_test$ts_scale == 'proportion') | ts_test$ts_scale == 'other'){
                            if(TRUE){
                                for(l in fcst_indices){
                                
                                ## get info_packet and trim
                                info_packet <- ts_test
                                info_packet$ts <- info_packet$ts[1:l]
                                
                                WINDOW = 52
                                data_till_now = data.frame(value = info_packet$ts, t = 1:l)
                                if(nrow(data_till_now)>WINDOW){
                                  data_till_now = data_till_now[(nrow(data_till_now)-WINDOW):nrow(data_till_now),]
                                }
                                
                                # data_till_now_plus1 = data.frame(value = ts_test$ts[2:(l+1)], t = 2:(l+1))
                                # if(nrow(data_till_now_plus1)>WINDOW){
                                #   data_till_now_plus1 = data_till_now_plus1[(nrow(data_till_now_plus1)-WINDOW):nrow(data_till_now_plus1),]
                                # }

                                #### light smoothing and differencing and get last k 
                                # data_till_now_smoothed      <- gam(value~ s(t,k=pmax(round(nrow(data_till_now)/2))),data=data_till_now)$fitted.values
                                # data_till_now_smoothed_plus1      <- gam(value~ s(t,k=pmax(round(nrow(data_till_now_plus1)/2))),data=data_till_now_plus1)$fitted.values
                                
                                
                                data_till_now_smoothed = data_till_now$value
                                # data_till_now_smoothed = data_till_now_plus1$value
                                # 
                                #### could use gam smoother 
                                X_temp = tail((data_till_now_smoothed),WINDOW)
                                y_temp = ts_test$ts[(l+1):(l+h)]

                                if(length(X_temp) < WINDOW){
                                  X_temp = c(rep(NA, WINDOW-length(X_temp)), X_temp)
                                  # X_temp_plus1 = c(rep(NA, WINDOW-length(X_temp_plus1)), X_temp_plus1)
                                }
                                # if(length(X_temp) < WINDOW & l > (length(ts_test$ts)-WINDOW+1)){
                                #   X_temp = c(X_temp,rep(NA, WINDOW-length(X_temp)))
                                #   X_temp_plus1 = c(X_temp_plus1,rep(NA, WINDOW-length(X_temp_plus1)))
                                # }
                                
                                X <- rbind(X, X_temp)
                                # X_plus1 <- rbind(X_plus1, X_temp_plus1)
                                y <- rbind(y, y_temp)
                                lastobs = c(lastobs, ts_test$ts[l])
                                ### Assign monthly data date of 1st of month
                                if(ts_test$ts_time_cadence == 'monthly'){
                                  DATE = as.character(gsub('_','-',paste0(ts_test$ts_dates,'_1')))
                                }else{
                                  DATE = as.character(ts_test$ts_dates[l])
                                }
                                dates = c(dates, as.character(ts_test$ts_dates[l]))
                                print(l)
                              }
                              
                              
                              dat = data.frame(obs = unlist(lastobs), date = as.character(dates),
                                               y1 = y[,1], y2 = y[,2], y3 = y[,3], y4 = y[,4],
                                               FILES = paste0(eval_key,'_', j),
                                               ts_scale = ts_test$ts_scale)
                              dat$date = as.character(dat$date)
                              
                              X = data.frame(X)
                              colnames(X) = paste0('X_', 1:WINDOW)
                              dat = data.frame(dat, X)
                              
                              # X_plus1 = data.frame(X_plus1)
                              # colnames(X_plus1) = paste0('Xplus1_', 1:WINDOW)
                              # dat = data.frame(dat, X_plus1)
                              # 
                              write.csv(dat, file = paste0(my_path,'embeddings_lstm/',"embed_",eval_key, '_',j,".csv"), quote = F, row.names = F)
                             
                            }
                          }
                        }
                        xx <- 1
                        xx
                      }
stopCluster(cl)
