# Visualizations for the sMOA Paper
## Author: DA Osthus and AC Murph
## Date: August 2024
library(ggplot2)
library(ggpubr)
library(reshape2)
library(tidyverse)
library(dplyr)
library(ggridges)
library(ggrepel)
library(patchwork)
library(gridExtra)
library(grid)
library(ggExtra)
library(latex2exp)
library(tidyverse)
library(this.path)
library(parallel)
library(doParallel)
library(data.table)
theme_set(theme_bw())
setwd(paste0(this.path::here(), '/../'))
source('R/smoa_helpers.R')
savepath = 'data/'

scores          <- read.csv("data/scores_tot.csv")

names_of_models  <- unique(scores$model)
models_to_label <- c("COVIDhub-baseline", "COVIDhub-4_week_ensemble","COVIDhub-trained_ensemble")

file_with_results         <- "data/k_5_num_curves_18387_closest_4422_dispersion_10000_mlebound_10000_state_records"
# file_with_results <- "data/k_5_num_curves_18387_closest_4422_dispersion_10000_mlebound_1250_state_records"
sockettype <- "PSOCK"
savepath = './'
##
## Begin Murph code to gather the data:
full_graph_data                                                                                    <- NULL
mae_comparison_data                                                                                <- NULL
wis_comparison_data                                                                                <- NULL
parfctn = function(file_name_idx){
  library(dplyr)
  full_graph_data                                                                                    <- NULL
  mae_comparison_data                                                                                <- NULL
  wis_comparison_data                                                                                <- NULL
  file_name = list.files(file_with_results)[file_name_idx]
  
  results_list                                                                                     <- NULL
  load(file = paste(file_with_results,"/", file_name, sep = ""))
  
  model_count                                                                                      <- 1
  for(model_name in names_of_models){
    if(is.null(results_list[[model_name]])) next
    
    temp_comparison_data                                                                           <- results_list[[model_name]]
    temp_comparison_data                                                                           <- temp_comparison_data[which(!is.na(temp_comparison_data$abs_error)),]
    mae_comparison_data                                                                            <- rbind(mae_comparison_data, temp_comparison_data)
    
    temp_comparison_data                                                                           <- results_list[[model_name]]
    temp_comparison_data                                                                           <- temp_comparison_data[which(!is.na(temp_comparison_data$wis)),]
    wis_comparison_data                                                                            <- rbind(wis_comparison_data, temp_comparison_data)
    
    if(model_count == 1){
      # We only need to add in the MOA stuff once per state.
      temp_data                                                                                    <- results_list[[model_name]]
      temp_model_results                                                                           <- temp_data
      
      temp_data$abs_error                                                                          <- NULL
      temp_data$abs_error                                                                          <- temp_data$abs_error_moa
      temp_data$abs_error_moa                                                                      <- NULL
      temp_data$wis_error                                                                          <- temp_data$wis_error_moa
      temp_data$wis_error_moa                                                                      <- NULL
      temp_data$model.x                                                                            <- NULL
      temp_data$model.y                                                                            <- NULL
      temp_data$wis                                                                                <- NULL
      temp_data$model_name                                                                         <- rep('sMOA', times = nrow(temp_data))
      full_graph_data                                                                              <- rbind(full_graph_data, temp_data)
      
      temp_data                                                                                    <- temp_model_results 
      temp_data$abs_error_moa                                                                      <- NULL
      temp_data$wis_error                                                                          <- temp_data$wis
      temp_data$wis_error_moa                                                                      <- NULL
      temp_data$model.x                                                                            <- NULL
      temp_data$model.y                                                                            <- NULL
      temp_data$wis                                                                                <- NULL
      temp_data$model_name                                                                         <- rep(model_name, times = nrow(temp_data))
      full_graph_data                                                                              <- rbind(full_graph_data, temp_data)
    }else{
      temp_data                                                                                    <- results_list[[model_name]]
      temp_model_results                                                                           <- temp_data
      
      temp_data                                                                                    <- temp_model_results 
      temp_data$abs_error_moa                                                                      <- NULL
      temp_data$wis_error                                                                          <- temp_data$wis
      temp_data$wis_error_moa                                                                      <- NULL
      temp_data$model.x                                                                            <- NULL
      temp_data$model.y                                                                            <- NULL
      temp_data$wis                                                                                <- NULL
      temp_data$model_name                                                                         <- rep(model_name, times = nrow(temp_data))
      full_graph_data                                                                              <- rbind(full_graph_data, temp_data)
    }
    model_count                                                                                    <- model_count + 1
  }
  
  mae_comparison_data$weekday = weekdays(mae_comparison_data$forecast_date.y)
  wis_comparison_data$weekday = weekdays(wis_comparison_data$forecast_date.y)
  full_graph_data$weekday = weekdays(full_graph_data$forecast_date.y)
  mae_comparison_data = mae_comparison_data %>% subset(weekday%in%c('Sunday','Monday'))
  wis_comparison_data = wis_comparison_data %>% subset(weekday%in%c('Sunday','Monday'))
  full_graph_data =  full_graph_data %>% subset(is.na(weekday) | weekday%in%c('Sunday','Monday'))
  
  return(list(
    mae_comparison_data = mae_comparison_data,
    wis_comparison_data = wis_comparison_data,
    full_graph_data = full_graph_data
  )
  )
}

ncores <- 1
cl <- parallel::makeCluster(spec = ncores,type = sockettype) 
setDefaultCluster(cl)
registerDoParallel(cl)
sim_ts <- foreach(i=1:length(list.files(file_with_results)),
                  .verbose = T)%dopar%{
                    parfctn(i)
                  }
stopCluster(cl)  

full_graph_data                                                                                    <- NULL
mae_comparison_data                                                                                <- NULL
wis_comparison_data                                                                                <- NULL
for(list_idx in 1:length(sim_ts)){
  full_graph_data     = rbind(full_graph_data, sim_ts[[list_idx]]$full_graph_data)
  mae_comparison_data = rbind(mae_comparison_data, sim_ts[[list_idx]]$mae_comparison_data)
  wis_comparison_data = rbind(wis_comparison_data, sim_ts[[list_idx]]$wis_comparison_data)
}

full_graph_data$forecast_date                                                                    <- as.Date(full_graph_data$forecast_date.x)
mae_comparison_data$forecast_date                                                                    <- as.Date(mae_comparison_data$forecast_date.x)
wis_comparison_data$forecast_date                                                                    <- as.Date(wis_comparison_data$forecast_date.x)
names_other_than_smoa                                                                              <- unique(full_graph_data$model_name)
names_other_than_smoa                                                                              <- names_other_than_smoa[which(names_other_than_smoa!='sMOA')]
full_graph_data$model_name                                                                         <- factor(full_graph_data$model_name, levels = c('sMOA', names_other_than_smoa))

HORIZON                                                                                            <- 1

#####
# Dot plots that lauren suggested
mae_comparison_data$abs_error_model                                                                <- mae_comparison_data$abs_error
mae_comparison_data$abs_error_smoa                                                                 <- mae_comparison_data$abs_error_moa
mae_comparison_data$abs_error_moa <- NULL
# mae_comparison_data$abs_error <- NULL
mae_comparison_data$wis_error_moa <- NULL
mae_comparison_data$model.x <- NULL
mae_comparison_data$forecast_date.x <- NULL
mae_comparison_data$forecast_date.y <- NULL
mae_comparison_data$wis <- NULL
wis_comparison_data$wis_error_model                                                                <- wis_comparison_data$wis
wis_comparison_data$wis_error_smoa                                                                 <- wis_comparison_data$wis_error_moa
wis_comparison_data$abs_error_moa <- NULL
wis_comparison_data$wis_error_moa <- NULL
wis_comparison_data$wis <- NULL
# wis_comparison_data$abs_error <- NULL
wis_comparison_data$model.x <- NULL
wis_comparison_data$forecast_date.x <- NULL
wis_comparison_data$forecast_date.y <- NULL
wis_comparison_data$model.x <- NULL
names_other_than_smoa                                                                              <- unique(full_graph_data$model_name)
names_other_than_smoa                                                                              <- names_other_than_smoa[which(names_other_than_smoa!='sMOA')]
full_graph_data$model_name                                                                         <- factor(full_graph_data$model_name, levels = c('sMOA', names_other_than_smoa))


validdf <- read.csv('data/validationdata_2023-12-31.csv')
validdf <- subset(validdf, target_variable == "inc case")
validdf$location <- ifelse(validdf$location <= 9, paste0("0", as.character(validdf$location)), as.character(validdf$location))
validdf$location_name <- validdf$location
validdf$target_end_date <- as.Date(validdf$target_end_date)


# tmp_valid = read.csv('data/validationdata_2023-12-31.csv')

models_to_label <- c("COVIDhub-baseline", "COVIDhub-4_week_ensemble","COVIDhub-trained_ensemble")

## best in class models
bestinclass <- c("LNQ-ens1","COVIDhub-4_week_ensemble","USC-SI_kJalpha",
                 "LANL-GrowthRate","Microsoft-DeepSTIA","COVIDhub-trained_ensemble",
                 "CU-select","BPagano-RtDriven","COVIDhub-baseline","JHUAPL-Bucky")

## use this to add in the validation data
findClosestDate <- function(dates, refDate) {
  # Convert input to Date class if necessary
  dates <- as.Date(dates)
  refDate <- as.Date(refDate)
  
  # Filter for dates that are on or before the reference date
  valid_dates <- dates[dates <= refDate]
  
  # If there are no valid dates, return NA
  if(length(valid_dates) == 0) {
    return(NA)
  }
  
  # Return the maximum valid date (i.e., the closest to refDate)
  return(max(valid_dates))
}


## how many models did sMOA beat if we stopped forecasting at each date?
mae_comparison_data$bestinclass <- "no"
mae_comparison_data[mae_comparison_data$model %in% bestinclass,]$bestinclass <- "yes"
wis_comparison_data$bestinclass <- "no"
wis_comparison_data[wis_comparison_data$model %in% bestinclass,]$bestinclass <- "yes"
unqdates <- sort(unique(mae_comparison_data$forecast_date))
dd1 <- NULL
dd1big <- NULL
parfctn = function(i) {
  # for(i in 1:length(unqdates)){
  library(data.table)
  print(i)
  ## aggregate up to the model level
  tempdfm <- data.table(mae_comparison_data)[forecast_date <= unqdates[i],
                                             list(n_total_submissions = length(abs_error),
                                                  bestinclass = bestinclass[1],
                                                  ## MAE
                                                  competitor = mean(abs_error_model[!is.na(abs_error_model)], na.rm=T),
                                                  smoa = mean(abs_error_smoa[!is.na(abs_error_model)],na.rm=T),
                                                  prop_smoa_better = mean(abs_error_smoa <= abs_error_model)),by=c("model")]
  tempdfm$winner <- "sMOA"
  tempdfm[tempdfm$competitor < tempdfm$smoa,]$winner <- "competitor"
  tempdfm$metric <- "MAE"
  tempdfm$fcst_date_on_or_before <- unqdates[i]
  tempdfm$last_date = unqdates[i]
  tempdfm$us_cases <- subset(validdf, location_name == "US" & target_end_date == findClosestDate(unique(validdf$target_end_date), unqdates[i]))$value
  
  dd1m <- data.frame(metric = "MAE",
                     last_date = unqdates[i],
                     us_cases = subset(validdf, location_name == "US" & target_end_date == findClosestDate(unique(validdf$target_end_date), unqdates[i]))$value,
                     n_models = nrow(tempdfm),
                     n_bestinclass_models = nrow(tempdfm[tempdfm$bestinclass == "yes",]),
                     total_submissions = sum(tempdfm$n_total_submissions),
                     prop_of_models_smoa_beats = mean(tempdfm$winner == "sMOA"),
                     prop_of_bestinclass_models_smoa_beats = mean(tempdfm[tempdfm$bestinclass == "yes",]$winner == "sMOA"))
  
  ## aggregate up to the model level
  tempdfw <- data.table(wis_comparison_data)[forecast_date <= unqdates[i],
                                             list(n_total_submissions = length(abs_error),
                                                  bestinclass = bestinclass[1],
                                                  ## MAE
                                                  competitor = mean(wis_error_model[!is.na(wis_error_model)], na.rm=T),
                                                  smoa = mean(wis_error_smoa[!is.na(wis_error_model)],na.rm=T),
                                                  prop_smoa_better = mean(wis_error_smoa <= wis_error_model)),by=c("model")]
  tempdfw$winner <- "sMOA"
  tempdfw[tempdfw$competitor < tempdfw$smoa,]$winner <- "competitor"
  tempdfw$metric <- "WIS"
  tempdfw$fcst_date_on_or_before <- unqdates[i]
  tempdfw$last_date = unqdates[i]
  tempdfw$us_cases <- subset(validdf, location_name == "US" & target_end_date == findClosestDate(unique(validdf$target_end_date), unqdates[i]))$value
  
  dd1w <- data.frame(metric = "WIS",
                     last_date = unqdates[i],
                     us_cases = subset(validdf, location_name == "US" & target_end_date == findClosestDate(unique(validdf$target_end_date), unqdates[i]))$value,
                     n_models = nrow(tempdfw),
                     n_bestinclass_models = nrow(tempdfw[tempdfw$bestinclass == "yes",]),
                     total_submissions = sum(tempdfw$n_total_submissions),
                     prop_of_models_smoa_beats = mean(tempdfw$winner == "sMOA"),
                     prop_of_bestinclass_models_smoa_beats = mean(tempdfw[tempdfw$bestinclass == "yes",]$winner == "sMOA"))

  
  # dd1 <- rbind(dd1,dd1m, dd1w)
  # dd1big <- rbind(dd1big, tempdfw)
  return(list(dd1m = dd1m, 
              dd1w = dd1w,
              dd1big = tempdfw)
         )
}

sockettype <- "PSOCK"
ncores <- 10

cl <- parallel::makeCluster(spec = ncores,type = sockettype) #, outfile=""
setDefaultCluster(cl)
registerDoParallel(cl)
sim_ts <- foreach(i=1:length(unqdates),
               .verbose = T)%dopar%{
                 print(i)
                 parfctn(i)
               }
stopCluster(cl)   

dd1 <- NULL
dd1big <- NULL
for(list_idx in 1:length(sim_ts)){
  dd1   = rbind(dd1, sim_ts[[list_idx]]$dd1m)
  dd1    = rbind(dd1, sim_ts[[list_idx]]$dd1w)
  dd1big = rbind(dd1big, sim_ts[[list_idx]]$dd1big)
}


####################################################
### Figure 6
namelabeldf <- data.frame(x = ymd("2022-11-01"),
                          metric = c("MAE","MAE","WIS","WIS"),
                          y = c(.525, .85, .725, .95),
                          mylabel = c("Best-in-class Models","All Models","Best-in-class Models","All Models"))

ggplot(data=dd1, aes(x=last_date))+
  geom_hline(aes(yintercept = .5), linetype=I(1),color=I("grey"))+
  geom_line(aes(y=prop_of_models_smoa_beats))+
  geom_point(aes(y=prop_of_models_smoa_beats),size=I(2), fill = I("black"), shape=I(21))+
  geom_line(aes(y=prop_of_bestinclass_models_smoa_beats), color=I("red"))+
  geom_point(aes(y=prop_of_bestinclass_models_smoa_beats),size=I(2), fill = I("red"), color=I("black"), shape=I(21))+
  theme(legend.position = "bottom",
        plot.title = element_text(hjust = 0.5))+
  xlab("")+
  geom_text(aes(x = x, y=y, label=mylabel), data=subset(namelabeldf,mylabel == "Best-in-class Models"), color=I("red"))+
  geom_text(aes(x = x, y=y, label=mylabel), data=subset(namelabeldf,mylabel == "All Models"), color=I("black"))+
  facet_wrap(~metric,ncol=1)+
  scale_y_continuous(limits = c(0,1), breaks = seq(0,1,.2))+
  ylab("Proportion of Models sMOA outperforms")+
  ggtitle("The proportion of all models (black) and best-in-class models (red)\nsMOA outperforms if validation is stopped on the x-axis date")+
  scale_x_date(date_breaks = "3 months", date_labels = "%b %Y") 

