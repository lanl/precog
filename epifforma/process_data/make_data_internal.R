########################################
########################################
### Generate epiFFORMA Training Data ###
########################################
########################################


################################
### Get Global SLURM Options ###
################################

options(warn = -1)

suppressMessages(library("optparse"))

option_list <- list(
  make_option("--replicate_num", type="character", default=""),
  make_option("--num_to_run", type="character", default="")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options

# my_path = as.character(opt$my_path) # I don't think we need this anymore.
replicate_num = as.character(opt$replicate_num) 
num_to_run = as.character(opt$num_to_run) 


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
savetrainpath <- here::here("epifforma", "process_data")
syntheticpath <- here::here("epifforma", "raw_data", "synthetic", "output") #paste0(my_path, '../raw_data/synthetic/output/') 
figsavepath <- here::here("epifforma", "process_data", "figs") #paste0(savetrainpath,"figs/")
savemodelpath <- here::here("epifforma", "process_data", "features2weights") #paste0(my_path, 'features2weights/')

## source in function
source(here::here("epifforma", "process_data", "epi_functions.R")) #source("epi_functions.R")

#### set up computer for parallelization
## define number of cores
ncores <- floor(.5*detectCores())

## define socket type
sockettype <- "PSOCK"




##############################################################
## STEP 0: READ IN THE TRAINING DATA (REPLACE WITH JSON(?))
##############################################################

## synthetic time series data for training
synthetic <- readRDS(here::here("epifforma", "raw_data", "synthetic", "output", "synthetic.RDS")) #readRDS(paste0(syntheticpath,"synthetic.RDS"))

## forecast horizon
h = 4

### Read in embedding matrices
embed_mat_X <- data.table::fread(file=paste0(savetrainpath,"/embed_mat/embed_mat_X.csv"))
embed_mat_y <- data.table::fread(file=paste0(savetrainpath,"/embed_mat/embed_mat_y.csv"))

embed_mat_X_deriv <- data.table::fread(file=paste0(savetrainpath,"/embed_mat/embed_mat_X_deriv.csv"))
embed_mat_y_deriv <- data.table::fread(file=paste0(savetrainpath,"/embed_mat/embed_mat_y_deriv.csv"))

gc()

##############################################################
## STEP 1A: MAKE TRAINING DATA PACKETS
##############################################################

## set up cluster
cl <- parallel::makeCluster(spec = ncores,
                            type = sockettype)
setDefaultCluster(cl)
registerDoParallel(cl)

print(Sys.time())
train_data <- foreach(i=1:num_to_run, #added extra 300 to compensate for extra seasonality sims
                      # .errorhandling = "pass",
                      .verbose = T)%dopar%{
                        
                        source(here::here("epifforma", "process_data", "epi_functions.R"))
                        
                        foutput_se <- NULL
                        coutput <- NULL
                        for(j in 1:100){
                          ts_id <- sample(1:length(synthetic),1)
                          #print(ts_id)
                          info_packet <- synthetic[[ts_id]]
                          
                          if(max(info_packet$ts)>0){ #protection against extremely rare all-zero edge case from seasonal sims
                            results <- build_training_data(info_packet, h)
                            gc()
                            features_output_se = results$features_se #includes all models with wt > equal using squared error-based weights
                            features_output_se$ts_id <- ts_id
                            components_output = results$components
                            components_output$ts_id <- ts_id
                            
                            foutput_se <- rbind(foutput_se, features_output_se)
                            coutput <- rbind(coutput, components_output)
                          }
                          print(j)
                        }
                        fwrite(foutput_se, file = paste0(savetrainpath,"/training_data_packets/features_se_",replicate_num,'_',i,"_",Sys.time(),".csv"))
                        fwrite(coutput, file = paste0(savetrainpath,"/training_data_packets/components_",replicate_num,'_',i,"_",Sys.time(),".csv"))
                        
                        x <- 1
                        x
                        
                      }
stopCluster(cl)
print(Sys.time())

print(train_data)
