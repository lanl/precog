# Visualizations comparing the online 'adding history' version
## Author: Alexander C. Murph
## Date: October 2024
library(ggplot2)
library(tidyverse)
library(ggridges)
library(ggrepel)
library(patchwork)
library(gridExtra)
library(ggExtra)
library(tidyverse)
library(parallel)
library(doParallel)
theme_set(theme_bw())
setwd("~/GitLab/smoa")

scores                          <- read.csv("data/scores_tot.csv")
names_of_models  <- unique(scores$model)
models_to_label <- c("COVIDhub-baseline", "COVIDhub-4_week_ensemble","COVIDhub-trained_ensemble")

file_with_results                                                                                  <- "data/k_5_num_curves_18387_closest_4422_dispersion_10000_mlebound_10000_state_records"
file_with_results_wOnlineHistories                                                                 <- "data/wHistories_k_5_num_curves_18387_closest_4422_dispersion_10000_mlebound_10000_state_records"
file_with_results_wOnlyHistory                                                                     <- "data/wOnlyHistories_k_5_num_curves_18387_closest_4422_dispersion_10000_mlebound_10000_state_records"

h <- 4
max_y_lim <- 13000

###########
## Gather the data that do not have the online histories:
full_graph_data                                                                                    <- NULL
mae_comparison_data                                                                                <- NULL
wis_comparison_data                                                                                <- NULL
for(file_name in list.files(file_with_results)){
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
      
      temp_data                                                                                    <- temp_model_results #[which(!(is.na(temp_model_results$abs_error)&is.na(temp_model_results$wis))),]
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
      
      temp_data                                                                                    <- temp_model_results #[which(!(is.na(temp_model_results$abs_error)&is.na(temp_model_results$wis))),]
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
}
full_graph_data$target_end_date                                                                    <- as.Date(full_graph_data$target_end_date)
names_other_than_smoa                                                                              <- unique(full_graph_data$model_name)
names_other_than_smoa                                                                              <- names_other_than_smoa[which(names_other_than_smoa!='sMOA')]
full_graph_data$model_name                                                                         <- factor(full_graph_data$model_name, levels = c('sMOA', names_other_than_smoa))

mae_comparison_data$abs_error_model                                                                <- mae_comparison_data$abs_error
mae_comparison_data$abs_error_smoa                                                                 <- mae_comparison_data$abs_error_moa
wis_comparison_data$wis_error_model                                                                <- wis_comparison_data$wis
wis_comparison_data$wis_error_smoa                                                                 <- wis_comparison_data$wis_error_moa



###########
## Gather the data that DOES have the online histories:
full_graph_data_wOnline                                                                                    <- NULL
mae_comparison_data_wOnline                                                                                <- NULL
wis_comparison_data_wOnline                                                                                <- NULL
for(file_name in list.files(file_with_results_wOnlineHistories)){
  results_list                                                                                     <- NULL
  load(file = paste(file_with_results_wOnlineHistories,"/", file_name, sep = ""))
  
  model_count                                                                                      <- 1
  for(model_name in names_of_models){
    if(is.null(results_list[[model_name]])) next
    
    temp_comparison_data                                                                           <- results_list[[model_name]]
    temp_comparison_data                                                                           <- temp_comparison_data[which(!is.na(temp_comparison_data$abs_error)),]
    mae_comparison_data_wOnline                                                                    <- rbind(mae_comparison_data_wOnline, temp_comparison_data)
    
    temp_comparison_data                                                                           <- results_list[[model_name]]
    temp_comparison_data                                                                           <- temp_comparison_data[which(!is.na(temp_comparison_data$wis)),]
    wis_comparison_data_wOnline                                                                            <- rbind(wis_comparison_data_wOnline, temp_comparison_data)
    
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
      full_graph_data_wOnline                                                                      <- rbind(full_graph_data_wOnline, temp_data)
      
      temp_data                                                                                    <- temp_model_results #[which(!(is.na(temp_model_results$abs_error)&is.na(temp_model_results$wis))),]
      temp_data$abs_error_moa                                                                      <- NULL
      temp_data$wis_error                                                                          <- temp_data$wis
      temp_data$wis_error_moa                                                                      <- NULL
      temp_data$model.x                                                                            <- NULL
      temp_data$model.y                                                                            <- NULL
      temp_data$wis                                                                                <- NULL
      temp_data$model_name                                                                         <- rep(model_name, times = nrow(temp_data))
      full_graph_data_wOnline                                                                              <- rbind(full_graph_data_wOnline, temp_data)
    }else{
      temp_data                                                                                    <- results_list[[model_name]]
      temp_model_results                                                                           <- temp_data
      
      temp_data                                                                                    <- temp_model_results #[which(!(is.na(temp_model_results$abs_error)&is.na(temp_model_results$wis))),]
      temp_data$abs_error_moa                                                                      <- NULL
      temp_data$wis_error                                                                          <- temp_data$wis
      temp_data$wis_error_moa                                                                      <- NULL
      temp_data$model.x                                                                            <- NULL
      temp_data$model.y                                                                            <- NULL
      temp_data$wis                                                                                <- NULL
      temp_data$model_name                                                                         <- rep(model_name, times = nrow(temp_data))
      full_graph_data_wOnline                                                                              <- rbind(full_graph_data_wOnline, temp_data)
    }
    model_count                                                                                    <- model_count + 1
  }
}
full_graph_data_wOnline$target_end_date                                                            <- as.Date(full_graph_data_wOnline$target_end_date)
names_other_than_smoa                                                                              <- unique(full_graph_data_wOnline$model_name)
names_other_than_smoa                                                                              <- names_other_than_smoa[which(names_other_than_smoa!='sMOA')]
full_graph_data_wOnline$model_name                                                                 <- factor(full_graph_data_wOnline$model_name, levels = c('sMOA', names_other_than_smoa))


#####
mae_comparison_data_wOnline$abs_error_model                                                                <- mae_comparison_data_wOnline$abs_error
mae_comparison_data_wOnline$abs_error_smoa                                                                 <- mae_comparison_data_wOnline$abs_error_moa
wis_comparison_data_wOnline$wis_error_model                                                                <- wis_comparison_data_wOnline$wis
wis_comparison_data_wOnline$wis_error_smoa                                                                 <- wis_comparison_data_wOnline$wis_error_moa


###########
## Gather the data that ONLY HAS the online history:
full_graph_data_wOnlyHistory                                                                                    <- NULL
mae_comparison_data_wOnlyHistory                                                                                <- NULL
wis_comparison_data_wOnlyHistory                                                                                <- NULL
for(file_name in list.files(file_with_results_wOnlyHistory)){
  results_list                                                                                     <- NULL
  load(file = paste(file_with_results_wOnlyHistory,"/", file_name, sep = ""))
  
  model_count                                                                                      <- 1
  for(model_name in names_of_models){
    if(is.null(results_list[[model_name]])) next
    
    temp_comparison_data                                                                           <- results_list[[model_name]]
    temp_comparison_data                                                                           <- temp_comparison_data[which(!is.na(temp_comparison_data$abs_error)),]
    mae_comparison_data_wOnlyHistory                                                               <- rbind(mae_comparison_data_wOnlyHistory, temp_comparison_data)
    
    temp_comparison_data                                                                           <- results_list[[model_name]]
    temp_comparison_data                                                                           <- temp_comparison_data[which(!is.na(temp_comparison_data$wis)),]
    wis_comparison_data_wOnlyHistory                                                               <- rbind(wis_comparison_data_wOnlyHistory, temp_comparison_data)
    
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
      full_graph_data_wOnlyHistory                                                                 <- rbind(full_graph_data_wOnlyHistory, temp_data)
      
      temp_data                                                                                    <- temp_model_results #[which(!(is.na(temp_model_results$abs_error)&is.na(temp_model_results$wis))),]
      temp_data$abs_error_moa                                                                      <- NULL
      temp_data$wis_error                                                                          <- temp_data$wis
      temp_data$wis_error_moa                                                                      <- NULL
      temp_data$model.x                                                                            <- NULL
      temp_data$model.y                                                                            <- NULL
      temp_data$wis                                                                                <- NULL
      temp_data$model_name                                                                         <- rep(model_name, times = nrow(temp_data))
      full_graph_data_wOnlyHistory                                                                 <- rbind(full_graph_data_wOnlyHistory, temp_data)
    }else{
      temp_data                                                                                    <- results_list[[model_name]]
      temp_model_results                                                                           <- temp_data
      
      temp_data                                                                                    <- temp_model_results #[which(!(is.na(temp_model_results$abs_error)&is.na(temp_model_results$wis))),]
      temp_data$abs_error_moa                                                                      <- NULL
      temp_data$wis_error                                                                          <- temp_data$wis
      temp_data$wis_error_moa                                                                      <- NULL
      temp_data$model.x                                                                            <- NULL
      temp_data$model.y                                                                            <- NULL
      temp_data$wis                                                                                <- NULL
      temp_data$model_name                                                                         <- rep(model_name, times = nrow(temp_data))
      full_graph_data_wOnlyHistory                                                                 <- rbind(full_graph_data_wOnlyHistory, temp_data)
    }
    model_count                                                                                    <- model_count + 1
  }
}
full_graph_data_wOnlyHistory$target_end_date                                                       <- as.Date(full_graph_data_wOnlyHistory$target_end_date)
names_other_than_smoa                                                                              <- unique(full_graph_data_wOnlyHistory$model_name)
names_other_than_smoa                                                                              <- names_other_than_smoa[which(names_other_than_smoa!='sMOA')]
full_graph_data_wOnlyHistory$model_name                                                            <- factor(full_graph_data_wOnlyHistory$model_name, 
                                                                                                             levels = c('sMOA', names_other_than_smoa))
#####
mae_comparison_data_wOnlyHistory$abs_error_model                                                                <- mae_comparison_data_wOnlyHistory$abs_error
mae_comparison_data_wOnlyHistory$abs_error_smoa                                                                 <- mae_comparison_data_wOnlyHistory$abs_error_moa
wis_comparison_data_wOnlyHistory$wis_error_model                                                                <- wis_comparison_data_wOnlyHistory$wis
wis_comparison_data_wOnlyHistory$wis_error_smoa                                                                 <- wis_comparison_data_wOnlyHistory$wis_error_moa


max_y_lim <- 11000
max_x_lim <- 16000
label_size <- 3
force_pull_val <- 2
point_plot_1 <- ggplot(mae_comparison_data %>% dplyr::group_by(model.y) %>% dplyr::summarize(abs_error_model= mean(abs_error_model), 
                                                                                             abs_error_smoa = mean(abs_error_smoa)),
                       aes(x=abs_error_model,y=abs_error_smoa, color=model.y)) + geom_point()  + 
  ggtitle(paste(""))+
  geom_abline() + 
  geom_label_repel(aes(label = ifelse(model.y%in%models_to_label,as.character(model.y),'')),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   max.overlaps = 20,
                   min.segment.length = 0.01,
                   # force_pull = -0.01,
                   size = label_size)+
  theme(legend.position = 'none') + xlab("Mean MAE for Model") + ylab("Mean MAE for sMOA") + 
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 18),
        title = element_text(size=15)) +
  theme(panel.spacing = unit(0, "cm"),
        # plot.margin = margin(0, 0, 0, 0, "cm"), 
        plot.caption = element_blank()) +
  xlim(0,max_x_lim) + ylim(0, max_y_lim)


point_plot_2 <- ggplot(wis_comparison_data %>% dplyr::group_by(model.y) %>% dplyr::summarize(wis_error_model=mean(wis_error_model), 
                                                                                             wis_error_smoa = mean(wis_error_smoa)),
                       aes(x=wis_error_model,y=wis_error_smoa,color=model.y)) + geom_point()  + 
  ggtitle(paste(""))+
  geom_abline() + 
  geom_label_repel(aes(label = ifelse(model.y%in%models_to_label,as.character(model.y),'')),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   max.overlaps = 20,
                   min.segment.length = 0.001,
                   # nudge_y = -100,
                   force_pull = force_pull_val,
                   size = label_size)+
  theme(legend.position = 'none') + xlab("Mean WIS for Model") + ylab("Mean WIS") + 
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 18),
        title = element_text(size=15)) +
  theme(panel.spacing = unit(0, "cm"),
        # plot.margin = margin(0, 0, 0, 0, "cm"), 
        plot.caption = element_blank()) +
  xlim(0,max_x_lim) + ylim(0, max_y_lim)


point_plot_3 <- ggplot(mae_comparison_data_wOnline %>% dplyr::group_by(model.y) %>% dplyr::summarize(abs_error_model= mean(abs_error_model), 
                                                                                                     abs_error_smoa = mean(abs_error_smoa)),
                       aes(x=abs_error_model,y=abs_error_smoa, color=model.y)) + geom_point()  + 
  ggtitle(paste(""))+
  geom_abline() + 
  geom_label_repel(aes(label = ifelse(model.y%in%models_to_label,as.character(model.y),'')),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   max.overlaps = 20,
                   min.segment.length = 0.01,
                   # force_pull = -0.1,
                   size = label_size)+
  theme(legend.position = 'none') + xlab("Mean MAE for Model") + ylab("Mean MAE") + 
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 18),
        title = element_text(size=15)) +
  theme(panel.spacing = unit(0, "cm"),
        # plot.margin = margin(0, 0, 0, 0, "cm"), 
        plot.caption = element_blank()) +
  xlim(0,max_x_lim) + ylim(0, max_y_lim)


point_plot_4 <- ggplot(wis_comparison_data_wOnline %>% dplyr::group_by(model.y) %>% dplyr::summarize(wis_error_model=mean(wis_error_model), 
                                                                                                     wis_error_smoa = mean(wis_error_smoa)),
                       aes(x=wis_error_model,y=wis_error_smoa,color=model.y)) + geom_point()  + 
  ggtitle(paste(""))+
  geom_abline() + 
  geom_label_repel(aes(label = ifelse(model.y%in%models_to_label,as.character(model.y),'')),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   max.overlaps = 20,
                   min.segment.length = 0.01,
                   # nudge_y = -100,
                   force_pull = force_pull_val,
                   size = label_size)+
  theme(legend.position = 'none') + xlab("Mean WIS for Model") + ylab("Mean WIS") + 
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 18),
        title = element_text(size=15)) +
  theme(panel.spacing = unit(0, "cm"),
        # plot.margin = margin(0, 0, 0, 0, "cm"), 
        plot.caption = element_blank()) +
  xlim(0,max_x_lim) + ylim(0, max_y_lim)


point_plot_5 <- ggplot(mae_comparison_data_wOnlyHistory %>% dplyr::group_by(model.y) %>% dplyr::summarize(abs_error_model= mean(abs_error_model), 
                                                                                                     abs_error_smoa = mean(abs_error_smoa)),
                       aes(x=abs_error_model,y=abs_error_smoa, color=model.y)) + geom_point()  + 
  ggtitle(paste(""))+
  geom_abline()  + 
  geom_label_repel(aes(label = ifelse(model.y%in%models_to_label,as.character(model.y),'')),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   max.overlaps = 20,
                   min.segment.length = 0.01,
                   nudge_y = -200,
                   # force_pull = -0.1,
                   size = label_size)+
  theme(legend.position = 'none') + xlab("Mean MAE for Model") + ylab("Mean MAE") + 
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 18),
        title = element_text(size=15)) +
  theme(panel.spacing = unit(0, "cm"),
        # plot.margin = margin(0, 0, 0, 0, "cm"), 
        plot.caption = element_blank()) +
  xlim(0,max_x_lim) + ylim(0, max_y_lim)

point_plot_6 <- ggplot(wis_comparison_data_wOnlyHistory %>% dplyr::group_by(model.y) %>% dplyr::summarize(wis_error_model=mean(wis_error_model), 
                                                                                                     wis_error_smoa = mean(wis_error_smoa)),
                       aes(x=wis_error_model,y=wis_error_smoa,color=model.y)) + geom_point()  + 
  ggtitle(paste(""))+
  geom_abline()  + 
  geom_label_repel(aes(label = ifelse(model.y%in%models_to_label,as.character(model.y),'')),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   max.overlaps = 20,
                   min.segment.length = 0.01,
                   # nudge_y = -200,
                   force_pull = force_pull_val,
                   size = label_size)+
  theme(legend.position = 'none') + xlab("Mean MAE for Model") + ylab("Mean MAE") + 
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 18),
        title = element_text(size=15)) +
  theme(panel.spacing = unit(0, "cm"),
        # plot.margin = margin(0, 0, 0, 0, "cm"), 
        plot.caption = element_blank()) +
  xlim(0,max_x_lim) + ylim(0, max_y_lim)

# w/o online history
wis_comparisons <- wis_comparison_data %>% dplyr::group_by(model.y) %>% dplyr::summarize(wis_error_model=mean(wis_error_model), 
                                                                                         wis_error_smoa = mean(wis_error_smoa))
mean(wis_comparisons$wis_error_model >= wis_comparisons$wis_error_smoa)

mae_comparisons <- mae_comparison_data %>% dplyr::group_by(model.y) %>% dplyr::summarize(abs_error_model= mean(abs_error_model), 
                                                                                         abs_error_smoa = mean(abs_error_smoa))
mean(mae_comparisons$abs_error_model >= mae_comparisons$abs_error_smoa)

# w/ online history
wis_comparisons <- wis_comparison_data_wOnline %>% dplyr::group_by(model.y) %>% dplyr::summarize(wis_error_model=mean(wis_error_model), 
                                                                                         wis_error_smoa = mean(wis_error_smoa))
mean(wis_comparisons$wis_error_model >= wis_comparisons$wis_error_smoa)

mae_comparisons <- mae_comparison_data_wOnline %>% dplyr::group_by(model.y) %>% dplyr::summarize(abs_error_model= mean(abs_error_model), 
                                                                                         abs_error_smoa = mean(abs_error_smoa))
mean(mae_comparisons$abs_error_model >= mae_comparisons$abs_error_smoa)

# w/ only history
wis_comparisons <- wis_comparison_data_wOnlyHistory %>% dplyr::group_by(model.y) %>% dplyr::summarize(wis_error_model=mean(wis_error_model), 
                                                                                                 wis_error_smoa = mean(wis_error_smoa))
mean(wis_comparisons$wis_error_model >= wis_comparisons$wis_error_smoa)

mae_comparisons <- mae_comparison_data_wOnlyHistory %>% dplyr::group_by(model.y) %>% dplyr::summarize(abs_error_model= mean(abs_error_model), 
                                                                                                 abs_error_smoa = mean(abs_error_smoa))
mean(mae_comparisons$abs_error_model >= mae_comparisons$abs_error_smoa)


row1 <- ggplot() + annotate(geom = 'text', x=1, y=1, label='bold("sMOA")', angle = 90, size = 5.8, parse = TRUE) + theme_void() 
row2 <- ggplot() + annotate(geom = 'text', x=1, y=1, label='bold("sMOA with Online History")', angle = 90, size = 5.8, parse = TRUE) + theme_void() 
row3 <- ggplot() + annotate(geom = 'text', x=1, y=1, label='bold("Classic MOA")', angle = 90, size = 5.8, parse = TRUE) + theme_void() 


layoutplot <- "
aeeeeeeeeeeeggggggggggg
aeeeeeeeeeeeggggggggggg
bdddddddddddfffffffffff
bdddddddddddfffffffffff
chhhhhhhhhhhiiiiiiiiiii
chhhhhhhhhhhiiiiiiiiiii
"

plotlist <- list(a = row1, b = row2, c = row3, 
                 e= point_plot_1, g=point_plot_2,
                 d = point_plot_3, f = point_plot_4,
                 h = point_plot_5, i = point_plot_6)
wrap_plots(plotlist, guides = 'collect', design = layoutplot)

