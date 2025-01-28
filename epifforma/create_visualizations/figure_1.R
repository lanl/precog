# Script for the Figure 1 in the epiFFORMA paper.
## Authors: Lauren Beesley and Alexander C. Murph
## Date: Novemember 2024


library(ggplot2)
library(data.table)
library(GGally)
library(viridis)
library(ggrepel)
library(gridExtra)
theme_set(theme_classic())

library(this.path)
my_path = this.path::here()
setwd(my_path)

## define paths
savepath <- paste0(my_path,"/../evaluate_model/evaluation/")
savetrainpath <- paste0(my_path,"/../process_data/")
syntheticpath <- paste0(my_path,"/../raw_data/synthetic/output/") # Murph note 4/9/24: This used to be synthetic_laurentest.  I changed this!!!!
codepath <- paste0(my_path,"/")
outputname <- "us_covid_rollercoaster"
plotpath <- paste0(my_path,"/../evaluate_model/evaluation/",outputname,'/')

## read in forecasting results
output_list_long = readRDS(file = paste0(savepath, outputname, "_order_multierror.RDS"))

FEATURE_LIST = setdiff(names(output_list_long['features']$features),c('geography','disease','last_obs_time','h'))
results = data.frame(output_list_long['plot_df']$plot_df)
results = merge(results, data.frame(output_list_long['features']$features), by = c('geography','disease','last_obs_time','h'), all.x = T, all.y = T)
results$error = (results$fcst - results$truth)^2
results = results %>% dplyr::group_by(geography,disease,last_obs_time,h) %>% dplyr::mutate(error_epifforma = error[model == 'epifforma'][1]) #epifforma error
results = results %>% dplyr::group_by(geography,disease,last_obs_time) %>% dplyr::mutate(error_epifforma_max = error[model == 'epifforma'][1]) #max epifforma error across h


FEATURE_LIST = setdiff(names(output_list_long['features']$features),c('geography','disease','last_obs_time','h'))
results = data.frame(output_list_long['plot_df']$plot_df)
results = merge(results, data.frame(output_list_long['features']$features), by = c('geography','disease','last_obs_time','h'), all.x = T, all.y = T)
results$error = (results$fcst - results$truth)^2
results = results %>% dplyr::group_by(geography,disease,last_obs_time,h) %>% dplyr::mutate(error_epifforma = error[model == 'epifforma'][1]) #epifforma error
results = results %>% dplyr::group_by(geography,disease,last_obs_time) %>% dplyr::mutate(error_epifforma_max = error[model == 'epifforma'][1]) #max epifforma error across h


MODEL_LIST = setdiff(names(output_list_long['components']$components),c('geography','disease','last_obs_time','h'))
components_wide = merge(output_list_long['weights']$weights, output_list_long['weights_sd']$weights_sd, by = c('geography','disease','last_obs_time','h'), all.x = T, all.y = T)
components_wide = merge(components_wide, output_list_long['components']$components, by = c('geography','disease','last_obs_time','h'), all.x = T, all.y = T)
names(components_wide) = ifelse(!grepl('_sd$',names(components_wide)) & !grepl('_avg$',names(components_wide)) & !(names(components_wide) %in% c('geography','disease','last_obs_time','h')), paste0('fcst_',names(components_wide)),names(components_wide))
names(components_wide) = ifelse(grepl('_avg$',names(components_wide)), gsub('_avg$','',paste0('weights_',names(components_wide))),names(components_wide))
names(components_wide) = ifelse(grepl('_sd$',names(components_wide)), gsub('_sd$','',paste0('weights_sd_',names(components_wide))),names(components_wide))
components_long = NULL
for(i in 1:length(MODEL_LIST)){
  dat_temp = data.frame(components_wide[,c('geography','disease','last_obs_time','h')], 
                        weights = components_wide[,paste0('weights_',MODEL_LIST[i])],
                        weights_sd = components_wide[,paste0('weights_sd_',MODEL_LIST[i])],
                        fcst = components_wide[,paste0('fcst_',MODEL_LIST[i])],
                        model = MODEL_LIST[i])
  components_long = rbind(components_long, dat_temp)
}
components_long = merge(components_long, results[results$model != 'obs',c('geography','disease','last_obs_time','h','truth','model','error','error_epifforma')], by = c('geography','disease','last_obs_time','h','model'), all.x = T, all.y = F)





# Stuff that gives me the visualization for the paper.
K = 50
k = 25
results = subset(results, subset = (geography=='United States California'))
worst_errors = sort(unique(results$error_epifforma_max),decreasing = F)[1:K]

covid_CA_dates = c("2020-01-05", "2020-01-12", "2020-01-19", "2020-01-26", "2020-02-02", 
                   "2020-02-09", "2020-02-16", "2020-02-23", "2020-03-01", "2020-03-08", 
                   "2020-03-15", "2020-03-22", "2020-03-29", "2020-04-05",
                   "2020-04-12", "2020-04-19", "2020-04-26", "2020-05-03", "2020-05-10", "2020-05-17", 
                   "2020-05-24", "2020-05-31", "2020-06-07", "2020-06-14", "2020-06-21", "2020-06-28", "2020-07-05", "2020-07-12",
                   "2020-07-19", "2020-07-26", "2020-08-02", "2020-08-09", "2020-08-16", "2020-08-23", 
                   "2020-08-30", "2020-09-06", "2020-09-13", "2020-09-20", "2020-09-27", "2020-10-04", "2020-10-11", "2020-10-18",
                   "2020-10-25", "2020-11-01", "2020-11-08", "2020-11-15", "2020-11-22", "2020-11-29", 
                   "2020-12-06", "2020-12-13", "2020-12-20", "2021-01-03", "2021-01-10", "2021-01-17", "2021-01-24", "2021-01-31",
                   "2021-02-07", "2021-02-14", "2021-02-21", "2021-02-28", "2021-03-07", "2021-03-14", 
                   "2021-03-21", "2021-03-28", "2021-04-04", "2021-04-11", "2021-04-18", "2021-04-25", "2021-05-02", "2021-05-09",
                   "2021-05-16", "2021-05-23", "2021-05-30", "2021-06-06", "2021-06-13", "2021-06-20", 
                   "2021-06-27", "2021-07-04", "2021-07-11", "2021-07-18", "2021-07-25", "2021-08-01", "2021-08-08", "2021-08-15",
                   "2021-08-22", "2021-08-29", "2021-09-05", "2021-09-12", "2021-09-19", "2021-09-26", 
                   "2021-10-03", "2021-10-10", "2021-10-17", "2021-10-24", "2021-10-31", "2021-11-07", "2021-11-14", "2021-11-21",
                   "2021-11-28", "2021-12-05", "2021-12-12", "2021-12-19", "2022-01-02", "2022-01-09", 
                   "2022-01-16", "2022-01-23", "2022-01-30", "2022-02-06", "2022-02-13", "2022-02-20", "2022-02-27", "2022-03-06")



TO_RUN = data.frame(results[results$error_epifforma_max == worst_errors[k],])

geography = TO_RUN$geography[1]
disease = TO_RUN$disease[1]
last_obs_time = TO_RUN$last_obs_time[1]

results$dates = covid_CA_dates[results$x]
results$last_obs_time_dates = covid_CA_dates[results$last_obs_time]

### Subset Data 
results_sub = data.frame(results)
results_sub = results_sub[results_sub$geography == geography & results_sub$disease == disease & results_sub$last_obs_time == last_obs_time,]
results_sub = results_sub[!(results_sub$model %in% c('rw_external','equal_wt')),]
components_sub = components_long[components_long$geography == geography & 
                                   components_long$disease == disease & 
                                   components_long$last_obs_time == last_obs_time,]

p1 = ggplot(results_sub[results_sub$model != 'epifforma',])+
  geom_vline(aes(xintercept = as.Date(dates)),data = results_sub[results_sub$x==(max(results_sub$x)- max(results_sub$h)) & results_sub$model != 'epifforma',], linetype = 2, color = 'gray')+
  geom_line(aes(x=as.Date(dates),y=truth))+
  geom_point(aes(x=as.Date(dates),y=truth))+
  geom_line(aes(x=as.Date(dates), y = fcst, group = model),data = results_sub[results_sub$h>0 & results_sub$model != 'epifforma',],color = 'black', linewidth = 1.1)+
  geom_line(aes(x=as.Date(dates), y = fcst, group = model, color = model),data = results_sub[results_sub$h>0 & results_sub$model != 'epifforma',], linewidth = 2)+
  geom_line(aes(x=as.Date(dates), y = fcst, group = model),data = results_sub[results_sub$h>0 & results_sub$model == 'epifforma',], linewidth =2, linetype = 3, color = 'black')+
  theme_classic()+
  theme(legend.position = c(0.1,0.7))+
  xlab('Time')+
  ylab('Cases')+
  labs(title = paste0('Location: ',results_sub$geography[1]))+
  scale_color_brewer('Models',palette = 'Spectral', limits = sort(unique(results_sub$model[!(results_sub$model %in% c('obs','epifforma'))]))) + 
  xlim(as.Date("2020-06-15"), NA) +
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 15),
        title = element_text(size=15),legend.text = element_text(size = 15)) +
  theme(panel.spacing = unit(0, "cm"),
        plot.margin = margin(0, 0, 0, 0, "cm"), 
        plot.caption = element_blank()) 

### Light GBM Forecast Weights
p3 = ggplot(components_sub)+
  geom_bar(aes(x=factor(h),y=weights,fill = model), position = 'stack', stat = 'identity')+
  theme_classic()+
  scale_y_continuous(expand=c(0,0))+
  xlab('Forecast Horizon')+
  ylab('Average Weight')+
  labs(title = 'LightGBM Forecast Weights')+
  scale_fill_brewer('Models',palette = 'Spectral')+
  theme(legend.position = 'top')+ 
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 15),
        title = element_text(size=15),legend.text = element_text(size = 15)) +
  theme(panel.spacing = unit(0, "cm"),
        plot.margin = margin(0, 0, 0, 0, "cm"), 
        plot.caption = element_blank()) 

lst_p = list(p1, p3)
pdf(file = paste0(my_path,'/figure1.pdf'), width = 16, height = 10)
grid.arrange(lst_p[[1]], lst_p[[2]],
             layout_matrix = matrix(c(1, 2),
                                    byrow = TRUE, nrow = 2, ncol = 1))
dev.off()
