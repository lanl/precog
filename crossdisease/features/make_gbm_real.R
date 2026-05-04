##################################
##################################
### Generate GBM Training Data ###
##################################
##################################



my_path = '~/GitLab/frodo/features/'
setwd(my_path)
source(paste0(my_path,'/../SLURMarray.r'))

## define paths
savetrainpath <- my_path



################################################
### Generate Features Data Packets via SLURM ###
################################################


### Training Data for Primary Analysis
NUM_RUN = 10
for(i in 1:NUM_RUN){
 cmdLines <- paste0("Rscript --vanilla ",my_path,"/make_gbm_real_internal.R",
                    " --total_run=",NUM_RUN,
                    " --cur_run=",i)
 runSlurm <- slurmarray(cmdLines,sname="PP",stime="540",smem="5G",
                        soutdir=paste0('./logs/'), sparallel = 0)
 system(runSlurm)
}


### Additional Training Data using M4 Time Series
NUM_RUN = 3
for(i in 1:NUM_RUN){
  cmdLines <- paste0("Rscript --vanilla ",my_path,"/make_gbm_real_internal_m4.R",
                       " --total_run=",NUM_RUN,
                       " --cur_run=",i)
  runSlurm <- slurmarray(cmdLines,sname="PP",stime="540",smem="5G",
                           soutdir=paste0('./logs/'), sparallel = 0)
  system(runSlurm)
}


### Additional Training Data using M4 Time Series
NUM_RUN = 3
for(i in 1:NUM_RUN){
  cmdLines <- paste0("Rscript --vanilla ",my_path,"/make_gbm_real_internal_m4daily.R",
                     " --total_run=",NUM_RUN,
                     " --cur_run=",i)
  runSlurm <- slurmarray(cmdLines,sname="PP",stime="540",smem="5G",
                         soutdir=paste0('./logs/'), sparallel = 0)
  system(runSlurm)
}




##########################
##########################
### All Available Data ###
##########################
##########################

library(data.table)
library(parallel)
library(doParallel)

ncores <- floor(.9*detectCores())
sockettype <- "PSOCK"
lf <- list.files(paste0(savetrainpath,"results/"))[grep("real_train_mat_",list.files(paste0(savetrainpath,"results/")))]
lf = lf[!is.na(lf)]
lf = lf[!grepl('real_train_mat.csv',lf)]

lut_path = '~/GitLab/frodo/lookup_tables/'
DISEASES = readxl::read_xlsx(paste0(lut_path, 'Disease_ML.xlsx'))
DISEASES = data.frame(DISEASES)
DISEASES = DISEASES$FILES
lf_all = c()
for(j in 1:length(DISEASES)){
  lf_all = c(lf_all, lf[grepl(DISEASES[j], lf)])
}
lf_all = unique(lf_all)

## set up cluster
cl <- parallel::makeCluster(spec = ncores,
                            type = sockettype)
setDefaultCluster(cl)
registerDoParallel(cl)
print(Sys.time())
train_data <- foreach(i=1:length(lf), 
                      #.errorhandling = "pass",
                      #.verbose = T,
                      .packages = c('dplyr', 'tsfeatures', 'data.table'))%dopar%{
                        dat = data.table::fread(paste0(paste0(savetrainpath,"results/"),lf[i]))
                        dat = dat[dat$obs_smoothed >= 1e-9,] #don't care about the obs_smoothed = 0 situation. Using obs_smoothed instead of obs for outlier detection reasons
                        dat$date = as.character(dat$date)
                        dat = data.frame(dat)
                        if(lf[i] %in% lf_all){
                          scale_factor = 1
                        }else{
                          scale_factor = 2
                        }
                        dat = dat[sample(1:nrow(dat),size = floor(nrow(dat)/scale_factor), replace = F ),]
                        dat
                      }
stopCluster(cl)
print(Sys.time())

train_final = do.call('rbind',train_data) #big big data merge
train_final = subset(train_final, select = -c(last_obs_time,geography,disease_source))
train_final$ts_measurement_type[train_final$ts_measurement_type == 'incidence'] = 'incident_cases'

### training data subsetting ###
train_final = train_final[train_final$ts_time_cadence == 'weekly',]
train_final = subset(train_final, select = -c(day_of_week, day_of_year, day_of_week_forecast, day_of_year_forecast, days_until_xmas, days_until_easter,
                                              days_until_newyears, days_until_usthanksgiving))

data.table::fwrite(x=train_final,file =paste0(savetrainpath,"/results/real_train_mat.csv"),quote = FALSE, row.names = FALSE)



########################
########################
### Disease Specific ###
########################
########################


lut_path = '~/GitLab/frodo/lookup_tables/'
DISEASES = readxl::read_xlsx(paste0(lut_path, 'Disease_ML.xlsx'))
DISEASES = data.frame(DISEASES)
DISEASES = DISEASES$FILES


for(j in 1:length(DISEASES)){
  
  library(data.table)
  library(parallel)
  library(doParallel)
  
  ncores <- floor(.9*detectCores())
  
  ## define socket type
  sockettype <- "PSOCK"
  
  lf <- list.files(paste0(savetrainpath,"results/"))[grep("real_train_mat_",list.files(paste0(savetrainpath,"results/")))]
  lf = lf[!is.na(lf)]
  lf = lf[!grepl('real_train_mat.csv',lf)]
  
  lf_short = lf[grepl(DISEASES[j], lf)]
  if(DISEASES[j] == 'Influenza'){
    lf_short = lf_short[!grepl('Haemophilus_Influenzae', lf_short)]
  }
  
  
  ## set up cluster
  cl <- parallel::makeCluster(spec = ncores,
                              type = sockettype)
  setDefaultCluster(cl)
  registerDoParallel(cl)
  print(Sys.time())
  train_data <- foreach(i=1:length(lf_short), 
                        #.errorhandling = "pass",
                        #.verbose = T,
                        .packages = c('dplyr', 'tsfeatures', 'data.table'))%dopar%{
                          dat = data.table::fread(paste0(paste0(savetrainpath,"results/"),lf_short[i]))
                          dat = dat[dat$obs_smoothed >= 1e-9,] #don't care about the obs_smoothed = 0 situation. Using obs_smoothed instead of obs for outlier detection reasons
                          dat$date = as.character(dat$date)
                          dat = data.frame(dat)
                          dat
                        }
  stopCluster(cl)
  print(Sys.time())
  train_final = do.call('rbind',train_data) #big big data merge
  train_final = subset(train_final, select = -c(last_obs_time,geography,disease_source))
  train_final$ts_measurement_type[train_final$ts_measurement_type == 'incidence'] = 'incident_cases'
  ### training data subsetting ###
  train_final = train_final[train_final$ts_time_cadence == 'weekly',]
  train_final = subset(train_final, select = -c(day_of_week, day_of_year, day_of_week_forecast, day_of_year_forecast, days_until_xmas, days_until_easter,
                                                days_until_newyears, days_until_usthanksgiving))
  data.table::fwrite(x=train_final,file =paste0(savetrainpath,"/results_diseasespecific/",DISEASES[j],"_train_mat.csv"),quote = FALSE, row.names = FALSE)
  
}







########################
########################
### Disease x Source ###
########################
########################

lut_path = '~/GitLab/frodo/lookup_tables/'
DISEASE_SOURCES = readxl::read_xlsx(paste0(lut_path, 'Evaluations_ML.xlsx'))
DISEASE_SOURCES = data.frame(DISEASE_SOURCES)
DISEASE_SOURCES = DISEASE_SOURCES$FILES
DISEASE_SOURCES = gsub('.RDS','',DISEASE_SOURCES)

library(data.table)
library(parallel)
library(doParallel)

ncores <- floor(.9*detectCores())

## define socket type
sockettype <- "PSOCK"

lf <- list.files(paste0(savetrainpath,"results/"))[grep("real_train_mat_",list.files(paste0(savetrainpath,"results/")))]
lf = lf[!is.na(lf)]
lf = lf[!grepl('real_train_mat.csv',lf)]

for(j in 1:length(DISEASE_SOURCES)){

  lf_short = lf[grepl(DISEASE_SOURCES[j], lf)]

  ## set up cluster
  cl <- parallel::makeCluster(spec = ncores,
                              type = sockettype)
  setDefaultCluster(cl)
  registerDoParallel(cl)
  print(Sys.time())
  train_data <- foreach(i=1:length(lf_short), 
                        #.errorhandling = "pass",
                        #.verbose = T,
                        .packages = c('dplyr', 'tsfeatures', 'data.table'))%dopar%{
                          dat = data.table::fread(paste0(paste0(savetrainpath,"results/"),lf_short[i]))
                          dat = dat[dat$obs_smoothed >= 1e-9,] #don't care about the obs_smoothed = 0 situation. Using obs_smoothed instead of obs for outlier detection reasons
                          dat$date = as.character(dat$date)
                          dat = data.frame(dat)
                          dat
                        }
  stopCluster(cl)
  print(Sys.time())
  train_final = do.call('rbind',train_data) #big big data merge
  train_final = subset(train_final, select = -c(last_obs_time,geography,disease_source))
  train_final$ts_measurement_type[train_final$ts_measurement_type == 'incidence'] = 'incident_cases'
  ### training data subsetting ###
  train_final = train_final[train_final$ts_time_cadence == 'weekly',]
  train_final = subset(train_final, select = -c(day_of_week, day_of_year, day_of_week_forecast, day_of_year_forecast, days_until_xmas, days_until_easter,
                                                days_until_newyears, days_until_usthanksgiving))
  data.table::fwrite(x=train_final,file =paste0(savetrainpath,"/results_sourcespecific/",DISEASE_SOURCES[j],"_train_mat.csv"),quote = FALSE, row.names = FALSE)

}








###############
###############
### M4 Data ###
###############
###############

library(data.table)
library(parallel)
library(doParallel)

ncores <- floor(.9*detectCores())
sockettype <- "PSOCK"
lf <- list.files(paste0(savetrainpath,"results_m4/"))[grep("real_train_mat_",list.files(paste0(savetrainpath,"results_m4/")))]
lf = lf[!is.na(lf)]
lf = lf[!grepl('real_train_mat.csv',lf)]

## set up cluster
cl <- parallel::makeCluster(spec = ncores,
                            type = sockettype)
setDefaultCluster(cl)
registerDoParallel(cl)
print(Sys.time())
train_data <- foreach(i=1:length(lf), 
                      #.errorhandling = "pass",
                      #.verbose = T,
                      .packages = c('dplyr', 'tsfeatures', 'data.table'))%dopar%{
                        dat = data.table::fread(paste0(paste0(savetrainpath,"results_m4/"),lf[i]))
                        dat = dat[dat$obs_smoothed >= 1e-9,] #don't care about the obs_smoothed = 0 situation. Using obs_smoothed instead of obs for outlier detection reasons
                        dat$date = as.character(dat$date)
                        dat = data.frame(dat)
                        scale_factor = 1
                        dat = dat[sample(1:nrow(dat),size = floor(nrow(dat)/scale_factor), replace = F ),]
                        dat
                      }
stopCluster(cl)
print(Sys.time())

train_final = do.call('rbind',train_data) #big big data merge
train_final = subset(train_final, select = -c(last_obs_time,geography,disease_source))
train_final$ts_measurement_type[train_final$ts_measurement_type == 'incidence'] = 'incident_cases'

### training data subsetting ###
train_final = train_final[train_final$ts_time_cadence == 'weekly',]
train_final = subset(train_final, select = -c(day_of_week, day_of_year, day_of_week_forecast, day_of_year_forecast, days_until_xmas, days_until_easter,
                                              days_until_newyears, days_until_usthanksgiving))

data.table::fwrite(x=train_final,file =paste0(savetrainpath,"/results_m4/real_train_mat.csv"),quote = FALSE, row.names = FALSE)



