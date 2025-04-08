# Visualizations for the sMOA Paper
## Author: Alexander C. Murph
## Date: August 2024
library(ggplot2)
library(ggpubr)
library(tidyverse)
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
theme_set(theme_bw())
setwd(paste0(this.path::here(), '/../'))

scores                          <- read.csv("data/scores_tot.csv")
names_of_models  <- unique(scores$model)
models_to_label <- c("COVIDhub-baseline", "COVIDhub-4_week_ensemble","COVIDhub-trained_ensemble")

file_with_results                                                                                  <- "data/k_5_num_curves_18387_closest_4422_dispersion_10000_mlebound_10000_state_records"
sockettype <- "PSOCK"

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
names_other_than_smoa                                                                              <- unique(full_graph_data$model_name)
names_other_than_smoa                                                                              <- names_other_than_smoa[which(names_other_than_smoa!='sMOA')]
full_graph_data$model_name                                                                         <- factor(full_graph_data$model_name, levels = c('sMOA', names_other_than_smoa))

HORIZON                                                                                            <- 1


#####
# Dot plots that lauren suggested
mae_comparison_data$abs_error_model                                                                <- mae_comparison_data$abs_error
mae_comparison_data$abs_error_smoa                                                                 <- mae_comparison_data$abs_error_moa
mae_comparison_data$abs_error_moa <- NULL
mae_comparison_data$abs_error <- NULL
mae_comparison_data$wis_error_moa <- NULL
mae_comparison_data$model.x <- NULL
mae_comparison_data$forecast_date.x <- NULL
mae_comparison_data$forecast_date.y <- NULL
mae_comparison_data$wis <- NULL
# Triple check that we are making a fair comparison:
anyNA(mae_comparison_data)
save(mae_comparison_data, file = 'data/mae_comparison_data.RData')

ggplot_data = mae_comparison_data %>% dplyr::group_by(model) %>% dplyr::summarize(abs_error_model= mean(abs_error_model), 
                                                                                  abs_error_smoa = mean(abs_error_smoa))
ggplot_data = subset(ggplot_data, subset=((abs_error_model<16500)&(abs_error_smoa<13000)) )

triangle_data = data.frame(abs_error_model = c(-Inf, -Inf, 0, 13000),
                           abs_error_smoa = c(Inf, -Inf, -Inf, Inf),
                           model = rep(NA, times = 4))
point_plot_1 <- ggplot(ggplot_data,
                       aes(x=abs_error_model,y=abs_error_smoa, color=model)) + geom_point(size = 4)+
  theme(panel.grid = element_blank(),
        panel.border = element_blank()) +  # Add abline
  geom_polygon(data = triangle_data, aes(x=abs_error_model,y=abs_error_smoa), fill = "lightblue", alpha = 0.2) + 
  geom_point(
    data = subset(ggplot_data, (model%in%models_to_label)), # Highlighted points
    aes(x = abs_error_model, y = abs_error_smoa),
    shape = 21,
    color = "black",
    stroke = 1.5,
    size = 3
  )  + 
  ggtitle(paste("Individual MAE Comparisons between\nsMOA and ForecastHub Models"))+
  geom_abline(slope = 1, intercept = 0) + 
  geom_label_repel(data = ggplot_data, aes(label = ifelse(model%in%models_to_label,as.character(model),'')),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   max.overlaps = 20,
                   min.segment.length = 0.01,
                   # force_pull = -0.1,
                   size = 5) +
  theme(legend.position = 'none') + xlab("Mean MAE for Model") + ylab("Mean MAE for sMOA") + 
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 18),
        title = element_text(size=15)) +
  theme(panel.spacing = unit(0, "cm"),
        plot.margin = margin(0, 0, 0, 0, "cm"), 
        plot.caption = element_blank()) +
  annotate("text", x = min(ggplot_data$abs_error_model), y = max(ggplot_data$abs_error_smoa) - 1000, 
           label = "Model Outperforms sMOA", hjust = 0, vjust = 1, size=6) +  
  annotate("text", x = max(ggplot_data$abs_error_model)-1000, y = min(ggplot_data$abs_error_smoa)+1000, 
           label = "Model Outperformed by sMOA", hjust = 1, vjust = 0, size=6)+
  # xlim(0,16500) + ylim(0,10000)
  scale_x_continuous(limits=c(0,16500),expand = c(0,0),
                     oob = scales::oob_keep) +
  scale_y_continuous(limits=c(0, 13000),expand = c(0,0),
                     oob = scales::oob_keep)


wis_comparison_data$wis_error_model                                                                <- wis_comparison_data$wis
wis_comparison_data$wis_error_smoa                                                                 <- wis_comparison_data$wis_error_moa
wis_comparison_data$abs_error_moa <- NULL
wis_comparison_data$wis_error_moa <- NULL
wis_comparison_data$wis <- NULL
wis_comparison_data$abs_error <- NULL
wis_comparison_data$model.x <- NULL
wis_comparison_data$forecast_date.x <- NULL
wis_comparison_data$forecast_date.y <- NULL
wis_comparison_data$model.x <- NULL
# Confirm that we are making the appropriate, fair comparison:
anyNA(wis_comparison_data)

ggplot_data = wis_comparison_data %>% dplyr::group_by(model) %>% dplyr::summarize(wis_error_model= mean(wis_error_model), 
                                                                                  wis_error_smoa = mean(wis_error_smoa))
ggplot_data = subset(ggplot_data, subset=((wis_error_model<16500)&(wis_error_smoa<13000)) )
triangle_data = data.frame(wis_error_model = c(-Inf, -Inf, 0, 13000),
                           wis_error_smoa = c(Inf, -Inf, -Inf, Inf),
                           model = rep(NA, times = 4))

point_plot_2 = ggplot(ggplot_data,
                      aes(x=wis_error_model,y=wis_error_smoa, color=model)) + geom_point(size = 4)+
  theme(panel.grid = element_blank(),
        panel.border = element_blank()) +  # Add abline
  geom_polygon(data = triangle_data, aes(x=wis_error_model,y=wis_error_smoa), fill = "lightblue", alpha = 0.2) + 
  geom_point(
    data = subset(ggplot_data, (model%in%models_to_label)), # Highlighted points
    aes(x = wis_error_model, y = wis_error_smoa),
    shape = 21,
    color = "black",
    stroke = 1.5,
    size = 3
  )  + 
  ggtitle(paste("Individual WIS Comparisons between\nsMOA and ForecastHub Models"))+
  geom_abline(slope = 1, intercept = 0) + 
  geom_label_repel(data = ggplot_data, aes(label = ifelse(model%in%models_to_label,as.character(model),'')),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   max.overlaps = 20,
                   min.segment.length = 0.01,
                   # force_pull = -0.1,
                   size = 5) +
  theme(legend.position = 'none') + xlab("Mean WIS for Model") + ylab("Mean WIS for sMOA") +
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 18),
        title = element_text(size=15)) +
  theme(panel.spacing = unit(0, "cm"),
        plot.margin = margin(0, 0, 0, 0, "cm"), 
        plot.caption = element_blank()) +
  annotate("text", x = min(ggplot_data$wis_error_model), y = max(ggplot_data$wis_error_smoa)+1000, 
           label = "Model Outperforms sMOA", hjust = 0, vjust = 1, size=6) +  
  annotate("text", x = max(ggplot_data$wis_error_model)-500, y = min(ggplot_data$wis_error_smoa)+1000, 
           label = "Model Outperformed by sMOA", hjust = 1, vjust = 0, size=6) +
  scale_x_continuous(limits=c(0,16500),expand = c(0,0),
                     oob = scales::oob_keep) +
  scale_y_continuous(limits=c(0, 13000),expand = c(0,0),
                     oob = scales::oob_keep)

layoutplot                                                                                         <- "
eeeeeeeeeeegggggggggggg
"
plotlist                                                                                           <- list(e= point_plot_1, g=point_plot_2)
wrap_plots(plotlist, guides = 'collect', design = layoutplot)
save(wis_comparison_data, file = 'wis_comparison_data.rds')

## best in class models
bestinclass <- c("LNQ-ens1","COVIDhub-4_week_ensemble","USC-SI_kJalpha",
                 "LANL-GrowthRate","Microsoft-DeepSTIA","COVIDhub-trained_ensemble",
                 "CU-select","BPagano-RtDriven","COVIDhub-baseline","JHUAPL-Bucky")


wis_comparisons <- wis_comparison_data %>% dplyr::group_by(model) %>% dplyr::summarize(wis_error_model=mean(wis_error_model), 
                                                                                       wis_error_smoa = mean(wis_error_smoa))
mean(wis_comparisons$wis_error_model >= wis_comparisons$wis_error_smoa)
models_beat_us = wis_comparisons[wis_comparisons$wis_error_model < wis_comparisons$wis_error_smoa,]$model
print("best in class models that beat us in terms of wis")
intersect(bestinclass, models_beat_us)
print("best in class models that sMOA out-performed in terms of wis")
print(bestinclass[!bestinclass%in%models_beat_us])

save(mae_comparison_data, file = 'mae_comparison_data.rds')
mae_comparisons <- mae_comparison_data %>% dplyr::group_by(model) %>% dplyr::summarize(abs_error_model= mean(abs_error_model), 
                                                                                       abs_error_smoa = mean(abs_error_smoa))
mean(mae_comparisons$abs_error_model >= mae_comparisons$abs_error_smoa)
models_beat_us = mae_comparisons[(mae_comparisons$abs_error_model < mae_comparisons$abs_error_smoa),]$model
print("best in class models that beat us in terms of mae")
intersect(bestinclass, models_beat_us)
print("best in class models that sMOA out-performed in terms of mae")
print(bestinclass[!bestinclass%in%models_beat_us])
