# Visualizations for the sMOA Paper
## Author: Alexander C. Murph
## Date: August 2024
library(ggplot2)
library(tidyverse)
library(ggridges)
library(ggrepel)
library(patchwork)
library(gridExtra)
library(ggExtra)
library(tidyverse)
theme_set(theme_bw())
setwd("~/GitLab/smoa")

scores                          <- read.csv("data/scores_tot.csv")
names_of_models  <- unique(scores$model)
models_to_label <- c("COVIDhub-baseline", "COVIDhub-4_week_ensemble","COVIDhub-trained_ensemble")

file_with_results                                                                                  <- "data/k_5_num_curves_18387_closest_4422_dispersion_10000_mlebound_10000_state_records"

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

HORIZON                                                                                            <- 1


#####
# Dot plots that lauren suggested
mae_comparison_data$abs_error_model                                                                <- mae_comparison_data$abs_error
mae_comparison_data$abs_error_smoa                                                                 <- mae_comparison_data$abs_error_moa

point_plot_1 <- ggplot(mae_comparison_data %>% dplyr::group_by(model.y) %>% dplyr::summarize(abs_error_model= mean(abs_error_model), 
                                                                               abs_error_smoa = mean(abs_error_smoa)),
       aes(x=abs_error_model,y=abs_error_smoa, color=model.y)) + geom_point()  + 
  ggtitle(paste("Individual MAE Comparisons between\nsMOA and ForecastHub Models"))+
  geom_abline() + 
  geom_label_repel(aes(label = ifelse(model.y%in%models_to_label,as.character(model.y),'')),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   max.overlaps = 20,
                   min.segment.length = 0.01,
                   force_pull = -0.1,
                   size = 5)+
  theme(legend.position = 'none') + xlab("Mean MAE for Model") + ylab("Mean MAE for sMOA") + 
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 18),
        title = element_text(size=15)) +
  theme(panel.spacing = unit(0, "cm"),
        plot.margin = margin(0, 0, 0, 0, "cm"), 
        plot.caption = element_blank()) +
   xlim(0,16500) + ylim(0, 10000)


wis_comparison_data$wis_error_model                                                                <- wis_comparison_data$wis
wis_comparison_data$wis_error_smoa                                                                 <- wis_comparison_data$wis_error_moa
point_plot_2 <- ggplot(wis_comparison_data %>% dplyr::group_by(model.y) %>% dplyr::summarize(wis_error_model=mean(wis_error_model), 
                                                                               wis_error_smoa = mean(wis_error_smoa)),
       aes(x=wis_error_model,y=wis_error_smoa,color=model.y)) + geom_point()  + 
  ggtitle(paste("Individual WIS Comparisons between\nsMOA and ForecastHub Models"))+
  geom_abline() + 
  geom_label_repel(aes(label = ifelse(model.y%in%models_to_label,as.character(model.y),'')),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   max.overlaps = 20,
                   min.segment.length = 0.001,
                   force_pull = -0.04,
                   size = 5)+
  theme(legend.position = 'none') + xlab("Mean WIS for Model") + ylab("Mean WIS for sMOA") + 
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 18),
        title = element_text(size=15)) +
  theme(panel.spacing = unit(0, "cm"),
        plot.margin = margin(0, 0, 0, 0, "cm"), 
        plot.caption = element_blank()) +
  xlim(0,16500) + ylim(0, 10000)

layoutplot                                                                                         <- "
eeeeeeeeeeegggggggggggg
"
plotlist                                                                                           <- list(e= point_plot_1, g=point_plot_2)
wrap_plots(plotlist, guides = 'collect', design = layoutplot)

wis_comparisons <- wis_comparison_data %>% dplyr::group_by(model.y) %>% dplyr::summarize(wis_error_model=mean(wis_error_model), 
                                                                                  wis_error_smoa = mean(wis_error_smoa))
mean(wis_comparisons$wis_error_model >= wis_comparisons$wis_error_smoa)

mae_comparisons <- mae_comparison_data %>% dplyr::group_by(model.y) %>% dplyr::summarize(abs_error_model= mean(abs_error_model), 
                                                                                  abs_error_smoa = mean(abs_error_smoa))
mean(mae_comparisons$abs_error_model >= mae_comparisons$abs_error_smoa)


#####
# California for ABS error:
graph_data                                                                                         <- subset(full_graph_data, subset = (horizon==HORIZON))
graph_data                                                                                         <- subset(graph_data, subset = (location_name=='California'))
# ggplot(graph_data, aes(x = target_end_date, y = abs_error, color = model_name, size = model_name)) + geom_line() +
#   scale_size_manual("model_name", values = c(1.5, rep(0.5, times = (length(unique(graph_data$model_name))-1) ))) +
#   ggtitle(paste("MAE by Date for each Model in California for Forecast Horizon", HORIZON))


#####
# Median across states for ABS error:
temp_graph_data                                                                                    <- subset(full_graph_data, subset = (horizon==HORIZON))
graph_data                                                                                         <- NULL
for(curr_model_name in unique(temp_graph_data$model_name)){
  for(curr_date in unique(temp_graph_data$target_end_date)){
    # for(curr_location in state.name){
    subset_data                                                                                    <- subset(temp_graph_data, subset = (target_end_date==as.Date(curr_date))&(model_name==curr_model_name) ) #(location_name==curr_location)&
    temp_row                                                                                       <- data.frame(target_end_date = as.Date(curr_date), model_name = curr_model_name, abs_error = median(subset_data$abs_error), wis_error = median(subset_data$wis_error))
    graph_data                                                                                     <- rbind(graph_data, temp_row)
    # }
  }
}

graph_data$target_end_date                                                                         <- as.Date(graph_data$target_end_date)
names_other_than_smoa                                                                              <- unique(graph_data$model_name)
names_other_than_smoa                                                                              <- names_other_than_smoa[which(names_other_than_smoa!='sMOA')]
graph_data$model_name                                                                              <- factor(graph_data$model_name, levels = c('sMOA', names_other_than_smoa))


#####
# Median across states RANKED according to ABS error:
# (across all horizons)
temp_graph_data                                                                                    <- full_graph_data
graph_data                                                                                         <- NULL
for(curr_model_name in unique(temp_graph_data$model_name)){
  for(curr_date in unique(temp_graph_data$target_end_date)){
    subset_data                                                                                    <- subset(temp_graph_data, subset = (target_end_date==as.Date(curr_date))&(model_name==curr_model_name) ) #(location_name==curr_location)&
    temp_row                                                                                       <- data.frame(target_end_date = as.Date(curr_date), model_name = curr_model_name, abs_error = mean(subset_data$abs_error), wis_error = mean(subset_data$wis_error))
    graph_data                                                                                     <- rbind(graph_data, temp_row)
  }
}

rank_graph_data                                                                                    <- NULL
for(curr_date in unique(graph_data$target_end_date)){
  subset_data                                                                                      <- subset(graph_data, subset = (target_end_date==as.Date(curr_date)) ) #(location_name==curr_location)&
  subset_data                                                                                      <- subset_data[which(!is.na(subset_data$abs_error)),]
  subset_data$rank                                                                                 <- order(subset_data$abs_error)
  subset_data$rank                                                                                 <- subset_data$rank / max(subset_data$rank)
  temp_row                                                                                         <- data.frame(target_end_date = as.Date(subset_data$target_end_date), model_name = subset_data$model_name, abs_error_rank = subset_data$rank)
  rank_graph_data                                                                                  <- rbind(rank_graph_data, temp_row)
}
rank_graph_data$model_name                                                                         <- as.character(rank_graph_data$model_name)
rank_graph_data$model_name[which(!(rank_graph_data$model_name%in%c('sMOA', 'COVIDhub-baseline', 'COVIDhub-4_week_ensemble', 'COVIDhub-trained_ensemble')))] <- 'Other Model'
rank_graph_data$model_name                                                                         <- factor(rank_graph_data$model_name, levels = c('sMOA', 'COVIDhub-baseline', 'COVIDhub-4_week_ensemble', 'COVIDhub-trained_ensemble', 'Other Model'))
rank_graph_data$Model                                                                              <- factor(rank_graph_data$model_name, levels = c('sMOA', 'COVIDhub-baseline', 'COVIDhub-4_week_ensemble', 'COVIDhub-trained_ensemble', 'Other Model'))
rank_graph_data = subset(rank_graph_data, subset = (model_name != 'Other Model'))
rank_graph_data = subset(rank_graph_data, subset = (target_end_date <= '2023-02-24')) # Added so the 4 week ensemble would have the same number of forecasts as the sMOA and the baseline

mean_data <- rank_graph_data %>% 
  dplyr::group_by(Model) %>% 
  dplyr::summarize(mean=mean(abs_error_rank)) 
p1 <- ggplot(rank_graph_data, aes(x=abs_error_rank, color=Model)) + 
  geom_density(alpha=0, size = 2)+ 
  geom_vline(data = mean_data, aes(xintercept = mean,  
                                color = Model), size=0.5)+ 
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 15),
        title = element_text(size=15),strip.text.x = element_text(size = 15),strip.text.y = element_text(size = 15)) +
  theme(panel.spacing = unit(0, "cm"),
        plot.margin = margin(0, 0, 0, 0, "cm"), 
        plot.caption = element_blank()) + theme(legend.text=element_text(size=15)) + ylab('Density') + xlab('MAE Rank at Date / Maximum MAE Rank at that Date')

rank_graph_data_mae <- rank_graph_data


rank_graph_data                                                                                    <- NULL
for(curr_date in unique(graph_data$target_end_date)){
  subset_data                                                                                      <- subset(graph_data, subset = (target_end_date==as.Date(curr_date)) ) #(location_name==curr_location)&
  subset_data                                                                                      <- subset_data[which(!is.na(subset_data$wis_error)),]
  subset_data$rank                                                                                 <- order(subset_data$wis_error)
  subset_data$rank                                                                                 <- subset_data$rank / max(subset_data$rank)
  temp_row                                                                                         <- data.frame(target_end_date = as.Date(subset_data$target_end_date), 
                                                                                                                 model_name = subset_data$model_name, 
                                                                                                                 wis_error_rank = subset_data$rank)
  rank_graph_data                                                                                  <- rbind(rank_graph_data, temp_row)
}
rank_graph_data$model_name                                                                         <- as.character(rank_graph_data$model_name)
rank_graph_data$model_name[which(!(rank_graph_data$model_name%in%c('sMOA', 'COVIDhub-baseline', 'COVIDhub-4_week_ensemble', 'COVIDhub-trained_ensemble')))] <- 'Other Model'
rank_graph_data$model_name                                                                         <- factor(rank_graph_data$model_name, levels = c('sMOA', 'COVIDhub-baseline', 'COVIDhub-4_week_ensemble', 'COVIDhub-trained_ensemble', 'Other Model'))
rank_graph_data$Model                                                                              <- factor(rank_graph_data$model_name, levels = c('sMOA', 'COVIDhub-baseline', 'COVIDhub-4_week_ensemble', 'COVIDhub-trained_ensemble', 'Other Model'))
rank_graph_data <- subset(rank_graph_data, subset = (model_name != 'Other Model'))
rank_graph_data = subset(rank_graph_data, subset = (target_end_date <= '2023-02-24')) # Added so the 4 week ensemble would have the same number of forecasts as the sMOA and the baseline

mean_data <- rank_graph_data %>% 
  dplyr::group_by(Model) %>% 
  dplyr::summarize(mean=mean(wis_error_rank)) 
p2 <- ggplot(rank_graph_data, aes(x=wis_error_rank, color=Model)) + 
  geom_density(alpha=0, size = 2)+ 
  geom_vline(data = mean_data, aes(xintercept = mean,  
                                   color = Model), size=0.5) + 
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 15),
        title = element_text(size=15),strip.text.x = element_text(size = 15),strip.text.y = element_text(size = 15)) +
  theme(panel.spacing = unit(0, "cm"),
        plot.margin = margin(0, 0, 0, 0, "cm"), 
        plot.caption = element_blank()) + theme(legend.text=element_text(size=15)) + ylab('Density') + xlab('WIS Rank at Date / Maximum WIS Rank at that Date')

lst_p = list(p1, p2)
grid.arrange(lst_p[[1]], lst_p[[2]],
             layout_matrix = matrix(c(1, 2),
                                    byrow = TRUE, nrow = 2, ncol = 1))

rank_graph_data                                                                                    <- NULL
for(curr_date in unique(graph_data$target_end_date)){
  subset_data                                                                                      <- subset(graph_data, subset = (target_end_date==as.Date(curr_date)) )
  subset_data                                                                                      <- subset_data[which(!is.na(subset_data$wis_error)),]
  subset_data$rank                                                                                 <- order(subset_data$wis_error)
  
  temp_row                                                                                         <- data.frame(target_end_date = as.Date(subset_data$target_end_date), model_name = subset_data$model_name, wis_error_rank = subset_data$rank)
  rank_graph_data                                                                                  <- rbind(rank_graph_data, temp_row)
}
rank_graph_data$model_name                                                                         <- as.character(rank_graph_data$model_name)
rank_graph_data$model_name[which(!(rank_graph_data$model_name%in%c('sMOA', 'COVIDhub-baseline')))] <- 'Other Model'
rank_graph_data$model_name                                                                         <- factor(rank_graph_data$model_name, levels = c('sMOA', 'COVIDhub-baseline', 'Other Model'))
rank_graph_data$Model                                                                              <- factor(rank_graph_data$model_name, levels = c('sMOA', 'COVIDhub-baseline', 'Other Model'))

temp_rank_graph_data <- subset(rank_graph_data, subset = (Model%in%c('sMOA', 'COVIDhub-baseline')))
wis_graph_data <- NULL
for(temp_date in unique(temp_rank_graph_data$target_end_date)){
  temp_data <- subset(temp_rank_graph_data, subset = (target_end_date==temp_date))
  comparison_sMOA <- temp_data[which(temp_data$Model=='sMOA'),]$wis_error_rank
  comparison_covid <- temp_data[which(temp_data$Model=='COVIDhub-baseline'),]$wis_error_rank
  sMOA_wins <- comparison_sMOA<comparison_covid
  wis_graph_data <- rbind(wis_graph_data, data.frame(target_end_date = as.Date(temp_date),
                                             wins = ifelse(sMOA_wins, 'sMOA', 'COVIDhub-baseline'),
                                             value = 1))
}
wis_graph_data$wins <- factor(wis_graph_data$wins, levels = c('sMOA', 'COVIDhub-baseline'))

rank_graph_data <- rank_graph_data_mae      
temp_rank_graph_data <- subset(rank_graph_data, subset = (Model%in%c('sMOA', 'COVIDhub-baseline')))
mae_graph_data <- NULL
for(temp_date in unique(temp_rank_graph_data$target_end_date)){
  temp_data <- subset(temp_rank_graph_data, subset = (target_end_date==temp_date))
  comparison_sMOA <- temp_data[which(temp_data$Model=='sMOA'),]$abs_error_rank
  comparison_covid <- temp_data[which(temp_data$Model=='COVIDhub-baseline'),]$abs_error_rank
  sMOA_wins <- comparison_sMOA<comparison_covid
  mae_graph_data <- rbind(mae_graph_data, data.frame(target_end_date = as.Date(temp_date),
                                             wins = ifelse(sMOA_wins, 'sMOA', 'COVIDhub-baseline'),
                                             value = 1))
}
mae_graph_data$wins <- factor(mae_graph_data$wins, levels = c('sMOA', 'COVIDhub-baseline'))

mean(ifelse(mae_graph_data$wins=='sMOA', 1, 0))
mean(ifelse(wis_graph_data$wins=='sMOA', 1, 0))





# Plot examples of several synthetic data
load(file ="data/synthetic_logs/synthetic_simidx_1_num_curves_22142.RData")

data_to_plot <- c(1, 499, 1055, 20000)

temp_data <- sim_ts[[data_to_plot[1]]]
graph_data <- data.frame(time = 1:length(temp_data$ts), value = temp_data$ts)
p1 <- ggplot(graph_data, aes(x = time, y = value)) + geom_line(linewidth = 1)+ 
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 15),
        title = element_text(size=15)) +
  theme(panel.spacing = unit(0, "cm"),
        plot.margin = margin(0, 0, 0, 0, "cm"), 
        plot.caption = element_blank()) 

temp_data <- sim_ts[[data_to_plot[2]]]
graph_data <- data.frame(time = 1:length(temp_data$ts), value = temp_data$ts)
p2 <- ggplot(graph_data, aes(x = time, y = value)) + geom_line(linewidth = 1)+ 
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 15),
        title = element_text(size=15)) +
  theme(panel.spacing = unit(0, "cm"),
        plot.margin = margin(0, 0, 0, 0, "cm"), 
        plot.caption = element_blank()) 

temp_data <- sim_ts[[data_to_plot[3]]]
graph_data <- data.frame(time = 1:length(temp_data$ts), value = temp_data$ts)
p3 <- ggplot(graph_data, aes(x = time, y = value)) + geom_line(linewidth = 1)+ 
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 15),
        title = element_text(size=15)) +
  theme(panel.spacing = unit(0, "cm"),
        plot.margin = margin(0, 0, 0, 0, "cm"), 
        plot.caption = element_blank()) 

temp_data <- sim_ts[[data_to_plot[4]]]
graph_data <- data.frame(time = 1:length(temp_data$ts), value = temp_data$ts)
p4 <- ggplot(graph_data, aes(x = time, y = value)) + geom_line(linewidth = 1)+ 
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 15),
        title = element_text(size=15)) +
  theme(panel.spacing = unit(0, "cm"),
        plot.margin = margin(0, 0, 0, 0, "cm"), 
        plot.caption = element_blank()) 


lst_p = list(p1, p2, p3, p4)
grid.arrange(lst_p[[1]], lst_p[[2]], lst_p[[3]], lst_p[[4]],
             layout_matrix = matrix(c(1, 2, 3, 4),
                                    byrow = TRUE, nrow = 2, ncol = 2),
             top = textGrob("Examples of Synthetic Data Time Series",gp=gpar(fontsize=20,font=3)))


