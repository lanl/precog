# Script for the Figure 3 in the epiFFORMA paper.
## Authors: Lauren Beesley and Alexander C. Murph
## Date: Novemember 2024
library(ggrepel)
library(patchwork)
library(gridExtra)
library(ggExtra)
library(tidyverse)
library(reshape2)
library(ggplot2)
library(data.table)
library(GGally)
library(viridis)
library(dplyr)
library(plyr)
theme_set(theme_classic())

# Start by making the boxplots.

library(this.path)
my_path = this.path::here()
setwd(my_path)

## define paths
savepath <- paste0(my_path,"/../evaluate_model/evaluation/")

## get .RDS file names
FILES = list.files(savepath)
FILES = FILES[grepl('.RDS',FILES)]
FILES = FILES[grepl('_order',FILES)] 
FILES = FILES[!grepl('synthetic',FILES)] ###added 1/16
 

### TAKES A LONG TIME!!!! ###
RESULTS = replicate(length(FILES),list(NULL))
RANKS = replicate(length(FILES),list(NULL))
RANKS_EQUAL_WT = replicate(length(FILES),list(NULL))
for(i in 1:length(FILES)){

  output_list = readRDS(paste0(savepath, FILES[i]))
  
  ### Get Accuracy Table
  RESULTS[[i]] = data.frame(output_list['accuracy_table']$accuracy_table,
                            fit_type = ifelse(grepl('multilogloss',FILES[i]), 'multi_logloss','multi_error'))
  
  
  ### Get Accuracy Table by Last Obs Time
  accuracydf = output_list['plot_df']$plot_df
  accuracydf = accuracydf[accuracydf$h != 0 & !is.na(accuracydf$h) & accuracydf$truth > 0 & accuracydf$type == 'fcst',]
  accuracydf$smape = abs(accuracydf$fcst - accuracydf$truth)/(0.5*(abs(accuracydf$fcst) + abs(accuracydf$truth)))
  accuracydf$mae = abs(accuracydf$fcst - accuracydf$truth)
  if('last_obs_time' %in% names(accuracydf)){
    accuracydf = accuracydf %>% dplyr::group_by(geography,last_obs_time,h) %>% dplyr::mutate(mae_rank = rank(mae, ties = 'min'))
  }else{
    accuracydf = accuracydf[,-1]
    accuracydf = accuracydf[!duplicated(paste0(accuracydf$ts_id, '_', accuracydf$ts_length, '_', accuracydf$h, '_', accuracydf$model)),]
    accuracydf = accuracydf %>% dplyr::group_by(ts_id, ts_length, h) %>% dplyr::mutate(mae_rank = rank(mae, ties = 'min'))
  }
  accuracydf_epifforma = accuracydf[accuracydf$model == 'epifforma',]
  AGG_MAE = prop.table(table(factor(accuracydf_epifforma$mae_rank, levels = c(1:max(accuracydf_epifforma$mae_rank)))))
  RANKS[[i]] = data.frame(disease = gsub('_order_multierror.RDS','',gsub('_order_multilogloss.RDS','',FILES[i])),
                          fit_type = ifelse(grepl('multilogloss',FILES[i]), 'multi_logloss','multi_error'),
                          prop = as.numeric(unlist(AGG_MAE)),
                          rank = 1:length(AGG_MAE),
                          type = 'MAE')
  
  accuracydf_equalwt = accuracydf[accuracydf$model == 'equal_wt',]
  AGG_MAE = prop.table(table(factor(accuracydf_equalwt$mae_rank, levels = c(1:max(accuracydf_equalwt$mae_rank)))))
  RANKS_EQUAL_WT[[i]] = data.frame(disease = gsub('_order_multierror.RDS','',gsub('_order_multilogloss.RDS','',FILES[i])),
                                   fit_type = ifelse(grepl('multilogloss',FILES[i]), 'multi_logloss','multi_error'),
                                   prop = as.numeric(unlist(AGG_MAE)),
                                   rank = 1:length(AGG_MAE),
                                   type = 'MAE')
  
  
  rm(output_list)
  rm(accuracydf)
  gc()
  print(paste0('Finished: ', i, ' of ', length(FILES)))
}

DISEASES = c('synthetic_training','synthetic_test',
             'us_covid_rollercoaster','global_covid_rollercoaster','us_ili_rollercoaster','dengue_rollercoaster',
             'us_diphtheria','us_measles','us_mumps','us_polio',
             'us_rubella','us_smallpox',
             'chikungunya')
DISEASE_LABELS = c('Synthetic (Training)', 'Synthetic (Test)',
                   'COVID-19 (US)', 'COVID-19 (Global)', 'ILI (US)','Dengue Fever (Global)', 
                   'Diphtheria (US)', 'Measles (US)', 'Mumps (US)',
                   'Polio (US)', 'Rubella (US)', 'Smallpox (US)', 
                   'Chikungunya (Brazil)')


RESULTS_STACKED = NULL
for(i in 1:length(FILES)){
  A = data.frame(RESULTS[[i]], key = gsub('.RDS','',FILES[i]))
  if(!('sbd_median' %in% names(A))){
    A$sbd_median = NA
    A$sbd_rank = NA
  }
  if(!('spearman' %in% names(A))){
    A$spearman = NA
    A$spearman_rank = NA
  }
  if(i > 1){
    RESULTS_STACKED = rbind(RESULTS_STACKED, A[,colnames(RESULTS_STACKED)])
  }else{
    RESULTS_STACKED = rbind(RESULTS_STACKED, A)
    
  }
}
RESULTS_STACKED$disease = gsub('_order_multierror','',gsub('_order_multilogloss','',RESULTS_STACKED$key))
RESULTS_STACKED = merge(RESULTS_STACKED, data.frame(disease = DISEASES, disease_pretty = DISEASE_LABELS), by = 'disease', all.x = T,
                        all.y = F)


RESULTS_STACKED_EXTRA = NULL
for(i in 1:length(FILES)){
  A = data.frame(RESULTS[[i]], key = gsub('.RDS','',FILES[i]))
  if('sbd_median' %in% names(A) & 'spearman' %in% names(A)){
    RESULTS_STACKED_EXTRA = rbind(RESULTS_STACKED_EXTRA, A)
  }
}
RESULTS_STACKED_EXTRA$disease = gsub('_order_multierror','',gsub('_order_multilogloss','',RESULTS_STACKED_EXTRA$key))
RESULTS_STACKED_EXTRA = merge(RESULTS_STACKED_EXTRA, data.frame(disease = DISEASES, disease_pretty = DISEASE_LABELS), by = 'disease', all.x = T,
                              all.y = F)


SUBSET = RESULTS_STACKED[RESULTS_STACKED$fit_type == 'multi_error',]
SUBSET = SUBSET %>% dplyr::group_by(disease) %>% dplyr::mutate(mae_ratio = mae/mae[which(model=='equal_wt')])
SUBSET = SUBSET %>% dplyr::group_by(disease) %>% dplyr::mutate(rmse_ratio = rmse/rmse[which(model=='equal_wt')])
SUBSET = SUBSET[,c(2,21,22,23)] 
SUBSET_melted = melt(SUBSET, id = c(1,2))
SUBSET_melted = SUBSET_melted[complete.cases(SUBSET_melted),]
SUBSET_melted = SUBSET_melted[which(SUBSET_melted$model!='equal_wt'),]
SUBSET_melted$model = factor(SUBSET_melted$model, levels = c('epifforma', 'rw', 'theta', 'arima', 'gam', 'gam2mirror', 'meanfcst', 'mirror', 'moa', 'moa_deriv'))
p2 = ggplot(SUBSET_melted)+
  geom_boxplot(aes(x=model, y = value, fill = variable))+
  # coord_cartesian(ylim = c(1,3))+
  theme_classic() + 
  #geom_boxplot(aes(y=value, x = model, fill = variable), alpha=.5, data = SUBSET_melted[SUBSET_melted$model == 'equal_wt',], color = 'black', size = 1.1)+
  geom_boxplot(aes(y=value, x = model, fill = variable), alpha=.5, data = SUBSET_melted[SUBSET_melted$model == 'epifforma',], color = 'black', size = 1.1)+ 
  theme(legend.position = 'top',axis.text=element_text(size=15),
        axis.title=element_text(size=18,face="bold"), plot.title = element_text(size=22),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=15))
boxplot_figure = p2 + geom_hline(yintercept = 1) + scale_y_continuous(trans='log10')
boxplot_figure = boxplot_figure +  scale_fill_discrete(type = c('#A3D5D3', '#ECD5B3'),
                                                       name = 'Metric',
                                                       labels = c('mae_ratio' = "Relative MAE", 
                                                                  'rmse_ratio' = 'Relative RMSE'),
                                                       breaks = c('mae_ratio', 'rmse_ratio'))




## define paths
savepath <- paste0(my_path,"/../evaluate_model/evaluation/")

## get .RDS file names
FILES = list.files(savepath)
FILES = FILES[grepl('.RDS',FILES)]
FILES = FILES[grepl('_order',FILES)]
FILES = FILES[!grepl('synthetic',FILES)]
FILES = FILES[grepl('multierror',FILES)]


### TAKES A LONG TIME!!!! ###
RESULTS_MAE = replicate(length(FILES),list(NULL))
RESULTS_WIS = replicate(length(FILES),list(NULL))
for(i in 1:length(FILES)){
  
  
  ### MAE Table
  output_list = readRDS(paste0(savepath, FILES[i]))
  RESULTS_MAE[[i]] = data.frame(output_list['accuracy_table']$accuracy_table,
                                key = FILES[i])
  rm(output_list)
  
  
  ### WIS Table
  output_list = readRDS(paste0(savepath,'../../evaluate_probabilistic/evaluation/',gsub('_multierror','' ,FILES[i])))
  accuracydf = ddply(output_list['results']$results,.(type),summarise,
                     median_truth = median(truth))
  RESULTS_WIS[[i]] = data.frame(output_list['accuracydf']$accuracydf,
                                median_truth = accuracydf$median_truth[1],
                                key = FILES[i])
  rm(output_list)
  
  
  gc()
  print(paste0('Finished: ', i, ' of ', length(FILES)))
}

DISEASES = c('synthetic_training','synthetic_test',
             'us_covid_rollercoaster','global_covid_rollercoaster','us_ili_rollercoaster','dengue_rollercoaster',
             'us_diphtheria','us_measles','us_mumps','us_polio',
             'us_rubella','us_smallpox',
             'chikungunya')
DISEASE_LABELS = c('Synthetic (Training)', 'Synthetic (Test)',
                   'COVID-19 (US)', 'COVID-19 (Global)', 'ILI (US)','Dengue Fever (Global)', 
                   'Diphtheria (US)', 'Measles (US)', 'Mumps (US)',
                   'Polio (US)', 'Rubella (US)', 'Smallpox (US)', 
                   'Chikungunya (Brazil)')


RESULTS_STACKED = NULL
for(i in 1:length(FILES)){
  A = RESULTS_MAE[[i]]
  B = RESULTS_WIS[[i]]
  A = A[A$model %in% c('equal_wt', 'epifforma'),]
  B = B[B$type %in% c('Epifforma (Interval)', 'Equalwt (Interval)'),]
  B$model = c('epifforma', 'equal_wt')
  
  RES = merge(A[,c('model', 'mae', 'rmse','key')], B[,c('model', 'median_width_scaled','coverage','median_msis','median_wis','median_truth')],
              by= c('model'), all.x = T, all.y = T)
  RESULTS_STACKED = rbind(RESULTS_STACKED, RES)
}
RESULTS_STACKED$disease = gsub('_order_multierror.RDS','',gsub('_order_multilogloss.RDS','',RESULTS_STACKED$key))
RESULTS_STACKED = merge(RESULTS_STACKED, data.frame(disease = DISEASES, disease_pretty = DISEASE_LABELS), by = 'disease', all.x = T,
                        all.y = F)

RESULTS_STACKED$disease = factor(RESULTS_STACKED$disease, levels = rev(DISEASES))
SUBSET = RESULTS_STACKED
SUBSET = SUBSET %>% dplyr::group_by(disease) %>% dplyr::mutate(mae_scaled = mae/median_truth)
SUBSET = SUBSET %>% dplyr::group_by(disease) %>% dplyr::mutate(wis_scaled = median_wis/median_truth)
SUBSET = SUBSET %>% dplyr::group_by(disease) %>% dplyr::mutate(rmse_scaled = rmse/median_truth)

SUBSET_all = SUBSET[,c(2,9,11,12,13,14)]
SUBSET_1 = SUBSET_all[which(SUBSET_all$model=='epifforma'),]
SUBSET_2 = SUBSET_all[which(SUBSET_all$model=='equal_wt'),]
graph_data = merge(SUBSET_1, SUBSET_2, by.x = 'disease_pretty', by.y = 'disease_pretty')
# graph_data = graph_data[graph_data$disease_pretty!='ILI (US)',]

models_to_label <- c("COVID-19 (US)", "ILI (US)", 'Dengue Fever (Global)', 'Mumps (US)')

graph_data$side_of_line <- factor(as.numeric(graph_data$mae_scaled.x<=graph_data$mae_scaled.y), levels = c('0','1'))
p1 <- ggplot(graph_data,aes(x = mae_scaled.x, y = mae_scaled.y,fill=side_of_line)) + 
  geom_point(color=I('black'), shape=I(21), size = I(3)) + geom_abline()+ 
  scale_fill_discrete(labels = c('0' = "0", 
                                 '1' = "1"),
                      breaks = c('0','1'),
                      type = c("lavender", "white"))+
  geom_label_repel(aes(label = ifelse(disease_pretty%in%models_to_label,as.character(disease_pretty),'')),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   size = 5)+ scale_x_continuous(trans='log10')+ scale_y_continuous(trans='log10')+
  theme(legend.position = 'none') + xlab("Average MAE for epiFFORMA") + ylab("Average MAE for Equal Weights") + 
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 18),
        title = element_text(size=15)) + coord_flip() +
  theme(panel.spacing = unit(0, "cm"),
        plot.margin = margin(0, 0, 0, 0, "cm"), 
        plot.caption = element_blank())

graph_data$side_of_line <- factor(as.numeric(graph_data$rmse_scaled.x<=graph_data$rmse_scaled.y), levels = c('0','1'))
p2 <- ggplot(graph_data,aes(x = rmse_scaled.x, y = rmse_scaled.y,fill=side_of_line)) + 
  geom_point(color=I('black'), shape=I(21), size = I(3)) + geom_abline()+ 
  scale_fill_discrete(labels = c('0' = "0", 
                                 '1' = "1"),
                      breaks = c('0','1'),
                      type = c("lavender", "white"))+
  geom_label_repel(aes(label = ifelse(disease_pretty%in%models_to_label,as.character(disease_pretty),'')),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   size = 5)+ scale_x_continuous(trans='log10')+ scale_y_continuous(trans='log10')+
  theme(legend.position = 'none') + xlab("Average RMSE for epiFFORMA") + ylab("Average RMSE for Equal Weights") + 
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 18),
        title = element_text(size=15)) + coord_flip() +
  theme(panel.spacing = unit(0, "cm"),
        plot.margin = margin(0, 0, 0, 0, "cm"), 
        plot.caption = element_blank()) 

graph_data$side_of_line <- factor(as.numeric(graph_data$wis_scaled.x<=graph_data$wis_scaled.y), levels = c('0','1'))
p3 <- ggplot(graph_data,aes(x = wis_scaled.x, y = wis_scaled.y,fill=side_of_line)) + 
  geom_point(color=I('black'), shape=I(21), size = I(3)) + geom_abline()+ 
  scale_fill_discrete(labels = c('0' = "0", 
                                 '1' = "1"),
                      breaks = c('0','1'),
                      type = c("white","lavender"))+
  geom_label_repel(aes(label = ifelse(disease_pretty%in%models_to_label,as.character(disease_pretty),'')),
                   box.padding   = 0.35, 
                   point.padding = 0.5,
                   segment.color = 'grey50',
                   size = 5)+ scale_x_continuous(trans='log10')+ scale_y_continuous(trans='log10')+
  theme(legend.position = 'none') + xlab("Median WIS for epiFFORMA") + ylab("Median WIS for Equal Weights") + 
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 18),
        title = element_text(size=15)) + coord_flip() +
  theme(panel.spacing = unit(0, "cm"),
        plot.margin = margin(0, 0, 0, 0, "cm"), 
        plot.caption = element_blank()) 


# Now visualize all plots together. 

layoutplot   <- "
aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa
eeeeeeeeeeegggggggggggghhhhhhhhhhhh
eeeeeeeeeeegggggggggggghhhhhhhhhhhh
"

plotlist  <- list(a=boxplot_figure, e = p1, g = p2, h = p3)

pdf(file = paste0(my_path,'/figure3.pdf'), width = 16, height = 10)
wrap_plots(plotlist, design = layoutplot)
dev.off()