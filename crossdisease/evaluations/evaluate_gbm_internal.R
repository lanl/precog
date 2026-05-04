
################################
### Get Global SLURM Options ###
################################

options(warn = -1)

suppressMessages(library("optparse"))

option_list <- list(
  make_option("--alpha", type="character", default=""),
  make_option("--total_run", type="character", default=""),
  make_option("--cur_run", type="character", default=""),
  make_option("--scale_ts", type="character", default="")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options

alpha = as.numeric(opt$alpha) 
total_run = as.numeric(opt$total_run) 
cur_run = as.numeric(opt$cur_run) 
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


if(scale_ts == 'no'){
  feature_path = '~/GitLab/frodo/features/results_unscaled/'
}else{
  feature_path = '~/GitLab/frodo/features/results/'
}

model_path = '~/GitLab/frodo/models/results/'


setwd(my_path)

FILES_LONG = gsub('.csv','',list.files(paste0(feature_path)))
FILES_LONG = FILES_LONG[FILES_LONG != 'real_train_mat']
lut_path = '~/GitLab/frodo/lookup_tables/'
DISEASE_SOURCES = readxl::read_xlsx(paste0(lut_path, 'Evaluations_ML.xlsx'))
DISEASE_SOURCES = data.frame(DISEASE_SOURCES)
DISEASE_SOURCES = DISEASE_SOURCES$FILES
DISEASE_SOURCES = gsub('.RDS','',DISEASE_SOURCES)
FILES_ALL = c()
for(i in 1:length(DISEASE_SOURCES)){
  FILES_ALL = c(FILES_ALL, FILES_LONG[grepl(DISEASE_SOURCES[i],FILES_LONG)])
}
FILES_ALL = unique(FILES_ALL)
FILES_CUR = FILES_ALL[seq(cur_run,length(FILES_ALL), by = total_run)]



FILENAME = paste0(my_path,'../models/results/cat_to_numeric.RDS')
LEVELS_STORAGE <- readRDS(FILENAME)



if(alpha == -1){
  ALPHA_GRID = c(0.025, 0.1, 0.25, 0.75, 0.9, 0.975)
}else{
  ALPHA_GRID = c(alpha)
}
for(alpha in ALPHA_GRID){




  
  YEAR_RUNS = c(2010:2024)
  
  if(scale_ts == 'no'){
    end_key = '_noscale.RDS'
  }else{
    end_key = '.RDS'
  }
  
  mod = replicate(length(YEAR_RUNS), list(NULL))
  for(y in 1:length(YEAR_RUNS)){
    res = try(readRDS(paste0(model_path,"fit_",YEAR_RUNS[y],"_",alpha, end_key)))
    if(class(res)[1] != 'try-error'){
      mod[[y]] <- readRDS(paste0(model_path,"fit_",YEAR_RUNS[y],"_",alpha, end_key))
    }else{
      mod[[y]] <- list(lgb.load(paste0(model_path,"fit_",YEAR_RUNS[y],"_",alpha, gsub('.RDS','.txt',end_key))))
    }
  }
  
  
  
  for(i in 1:length(FILES_CUR)){
    
    eval_key = FILES_CUR[i]
    eval_key_short = gsub('real_train_mat_','',eval_key)
    real_data <- read.csv(paste0(feature_path,"/",eval_key, ".csv"))
    real_data$scaling_factor = as.numeric(real_data$scaling_factor)
  
    PMIN = 0
    if(scale_ts == 'no'){
      PMAX = Inf
      real_data$scaling_factor = 1
    }else{
      PMAX = 3
    }

    OBS_STORE = real_data$obs
    OBS_STORE_SMOOTHED = real_data$obs_smoothed
    real_data$obs = pmax(pmin(real_data$obs / real_data$scaling_factor,PMAX),PMIN)
    real_data$obs_smoothed = pmax(pmin(real_data$obs_smoothed / real_data$scaling_factor,PMAX),PMIN)
    real_data$obs_smoothed_minus1 = pmax(pmin(real_data$obs_smoothed_minus1 / real_data$scaling_factor,PMAX),PMIN)
    real_data$sloperoll02 = as.numeric(real_data$sloperoll02)
    real_data$sloperoll12 = as.numeric(real_data$sloperoll12)
    
    IS_WEEKLY = as.numeric(real_data$ts_time_cadence[1] == 'weekly')
    real_data$transmission = as.numeric(factor(as.character(real_data$transmission), levels = c(LEVELS_STORAGE$transmission)))
    real_data$viral_bacterial_fungal = as.numeric(factor(as.character(real_data$viral_bacterial_fungal), levels = c(LEVELS_STORAGE$viral_bacterial_fungal)))
    real_data$ts_time_cadence = as.numeric(factor(as.character(real_data$ts_time_cadence), levels = c(LEVELS_STORAGE$ts_time_cadence)))
    real_data$ts_scale = as.numeric(factor(as.character(real_data$ts_scale), levels = c(LEVELS_STORAGE$ts_scale)))
    real_data$ts_measurement_type = as.numeric(factor(as.character(real_data$ts_measurement_type), levels = c(LEVELS_STORAGE$ts_measurement_type)))
    real_data$disease = as.numeric(factor(as.character(real_data$disease), levels = c(LEVELS_STORAGE$disease)))
    real_data$source = as.numeric(factor(as.character(real_data$source), levels = c(LEVELS_STORAGE$source)))
    

    FEATURE_LIST = colnames(real_data)
    FEATURE_LIST = FEATURE_LIST[!(FEATURE_LIST %in% c('truth', 'last_obs_time','geography','disease_source', 'y','date','obs','truth_sm',
                                                      'day_of_week', 'day_of_year', 'day_of_week_forecast', 'day_of_year_forecast', 'days_until_xmas', 'days_until_easter',
                                                      'days_until_newyears', 'days_until_usthanksgiving'))]
    FEATURE_LIST = sort(FEATURE_LIST)
    
  
    real_data$year = substr(real_data$date,1,4)
    TO_RUN = sort(intersect(unique(as.numeric(real_data$year)), YEAR_RUNS))
    
    if(length(TO_RUN)>0 & IS_WEEKLY==1){
      out_data = NULL
      OBS_STORE_SUBSET = OBS_STORE[real_data$year >= min(TO_RUN)]
      OBS_STORE_SMOOTHED_SUBSET = OBS_STORE_SMOOTHED[real_data$year >= min(TO_RUN)]
      real_data = real_data[real_data$year >= min(TO_RUN),]
      for(y in 1:length(TO_RUN)){
        SUBSET = real_data[real_data$year == TO_RUN[y],]
        y2 = which(YEAR_RUNS == TO_RUN[y])
        if(nrow(SUBSET)>1){
          for(h in 1:4){
            SUBSET_h = SUBSET[SUBSET$h == h,]
            
            pred_list<- predict(mod[[y2]][[h]], newdata = as.matrix(data.frame(SUBSET_h)[,FEATURE_LIST]))
            SUBSET_h$fcst = pmin(pmax(PMIN, pred_list),PMAX)*SUBSET_h$scaling_factor
            SUBSET_h$fcst[SUBSET_h$obs_smoothed == 0] = 0
            SUBSET_h$ratio <- pmin(pmax(PMIN, pred_list),PMAX)
            SUBSET_h$ratio[SUBSET_h$obs_smoothed == 0] = 0
    
    
            SUBSET_h$obs = OBS_STORE_SUBSET[real_data$year == TO_RUN[y] & real_data$h == h]
            SUBSET_h$obs_smoothed = OBS_STORE_SMOOTHED_SUBSET[real_data$year == TO_RUN[y] &  real_data$h == h]
            FCSTS = colnames(SUBSET_h)
            FCSTS = FCSTS[grepl('fcst',FCSTS)]
            FCSTS = c(FCSTS, 'truth','obs','obs_smoothed', 'h', 'disease', 'geography','disease_source','last_obs_time','date','ratio')
            FCSTS = unique(FCSTS)
            out_data = rbind(out_data,SUBSET_h[,FCSTS])
          }
        }
      }
    if(scale_ts == 'no'){
        FILENAME = paste0(my_path,"frodo_noscale/all/real_eval_mat_",eval_key_short,"_",alpha,".csv")
    }else{
        FILENAME = paste0(my_path,"frodo/all/real_eval_mat_",eval_key_short, "_",alpha,".csv")
    }
  
    write.csv(out_data, file = FILENAME, quote = F, row.names = F)
    gc()
    }
  
    
    
    print(paste0('Completed: ', FILES_CUR[i]))
  }


}