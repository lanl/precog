# Script for the Figure 3 in the epiFFORMA paper.
## Authors: Lauren Beesley and Alexander C. Murph
## Date: Novemember 2024
library(ggplot2)
library(data.table)
library(GGally)
library(viridis)
library(ggrepel)
library(dplyr)
library(gridGraphics)
library(plyr)
theme_set(theme_classic())

library(this.path)
my_path = this.path::here()
setwd(my_path)

## define paths
savepath <- paste0(my_path,"/../evaluate_model/evaluation/")

## get .RDS file names
FILES = list.files(savepath)
FILES = FILES[grepl('.RDS',FILES)]
FILES = FILES[grepl('_order',FILES)]
#FILES = FILES[!grepl('synthetic',FILES)]

### TAKES A LONG TIME!!!! ###
RESULTS = replicate(length(FILES),list(NULL))
RANKS = replicate(length(FILES),list(NULL))
RANKS_EQUAL_WT = replicate(length(FILES),list(NULL))
for(i in 1:length(FILES)){
  #for(i in 4:5){
  
  output_list = readRDS(paste0(savepath, FILES[i]))
  
  ### Get Accuracy Table
  RESULTS[[i]] = data.frame(output_list['accuracy_table']$accuracy_table,
                            fit_type = ifelse(grepl('multilogloss',FILES[i]), 'multi_logloss','multi_error'))
  
  
  ### Get Accuracy Table by Last Obs Time
  library(plyr)
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



RESULTS_STACKED$disease = factor(RESULTS_STACKED$disease, levels = rev(DISEASES))
SUBSET = RESULTS_STACKED[RESULTS_STACKED$fit_type == 'multi_error',]
p2 = ggplot(SUBSET)+
  geom_tile(aes(y=disease, x = mae_rank, color = model, fill = model))+
  geom_tile(aes(y=disease, x = mae_rank, fill = model), color = 'black', data = SUBSET[SUBSET$model == 'epifforma',], linewidth = 2)+
  theme_classic()+
  xlab('Rank')+
  ylab('Disease')+
  labs(title = 'MAE Rank')+
  #facet_grid(paste0('LightGBM: ', method)~., scales = 'free_x')+
  scale_fill_brewer('Model',palette = 'Spectral')+
  scale_color_brewer('Model',palette = 'Spectral')+
  scale_x_continuous(expand=c(0,0), breaks=c(1:length(unique(SUBSET$model))))+
  scale_y_discrete(breaks = DISEASES, labels = DISEASE_LABELS)+
  theme(legend.position = 'none')+
  guides(fill = guide_legend(ncol = 6))
p1 = ggplot(SUBSET)+
  geom_tile(aes(y=disease, x = rmse_rank, color = model, fill = model))+
  geom_tile(aes(y=disease, x = rmse_rank, fill = model), color = 'black', data = SUBSET[SUBSET$model == 'epifforma',], linewidth = 2)+
  theme_classic()+
  xlab('Rank')+
  ylab('')+
  labs(title = 'RMSE Rank')+
  #facet_grid(paste0('LightGBM: ', method)~., scales = 'free_x')+
  scale_fill_brewer('Model',palette = 'Spectral')+
  scale_color_brewer('Model',palette = 'Spectral')+
  scale_x_continuous(expand=c(0,0), breaks=c(1:length(unique(SUBSET$model))))+
  scale_y_discrete(breaks = DISEASES, labels = DISEASE_LABELS)+
  theme(legend.position = 'top')+
  guides(fill = guide_legend(ncol = 6))
legend = cowplot::get_plot_component(p1, 'guide-box-top', return_all = TRUE)
p1 = p1+theme(legend.position = "none")
p1 = p1 + theme(axis.text.y = element_blank())
p0 = cowplot::plot_grid(p2,p1,ncol = 2, nrow = 1, rel_widths = c(1.4,1))
p00 = cowplot::plot_grid(legend,p0,ncol = 1, nrow = 2, rel_heights = c(0.15,1))


pdf(file = paste0(my_path,'/figure4.pdf'), width = 8, height = 6)
p00
dev.off()


