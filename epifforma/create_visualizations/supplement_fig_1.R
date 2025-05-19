# Script for the Figure 3 in the epiFFORMA paper.
## Authors: Lauren Beesley and Alexander C. Murph
## Date: Novemember 2024
library(ggridges)
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

## get .RDS file names
FILES = list.files(savepath)
FILES = FILES[grepl('.RDS',FILES)]
FILES = FILES[grepl('_order',FILES)]
FILES = FILES[!grepl('synthetic',FILES)]
FILES = FILES[grepl('multierror',FILES)]


### TAKES A LONG TIME!!!! ###
RESULTS = replicate(length(FILES),list(NULL))
for(u in 1:length(FILES)){
  
  output_list = readRDS(paste0(savepath, FILES[u]))
  
  MODEL_LIST = setdiff(names(output_list['components']$components),c('geography','disease','last_obs_time','h'))
  components_wide = merge(output_list['weights']$weights, output_list['weights_sd']$weights_sd, by = c('geography','disease','last_obs_time','h'), all.x = T, all.y = T)
  components_wide = merge(components_wide, output_list['components']$components, by = c('geography','disease','last_obs_time','h'), all.x = T, all.y = T)
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
  components_long$x = components_long$last_obs_time + components_long$h
  DAT = data.frame(output_list['plot_df']$plot_df)
  DAT = DAT[!is.na(DAT$truth),]
  DAT = DAT[!duplicated(paste0(DAT$disease, '_', DAT$geography, '_', DAT$x)),]
  components_long = merge(components_long ,data.frame(DAT)[,c('geography', 'disease','x', 'truth')],
                          by = c('geography', 'disease','x'), all.x = T, all.y =F)
  
  ### Get Accuracy Table
  RESULTS[[u]] = data.frame(components_long, key = FILES[u])
  
  
  rm(output_list)
  gc()
  print(paste0('Finished: ', u, ' of ', length(FILES)))
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

RESULTS_STACKED = do.call('rbind', RESULTS)
RESULTS_STACKED$disease = gsub('_order_multierror.RDS','',gsub('_order_multilogloss.RDS','',RESULTS_STACKED$key))
RESULTS_STACKED = merge(RESULTS_STACKED, data.frame(disease = DISEASES, disease_pretty = DISEASE_LABELS), by = 'disease', all.x = T,
                        all.y = F)


RESULTS_STACKED$error = abs(RESULTS_STACKED$fcst - RESULTS_STACKED$truth)
RESULTS_STACKED = data.frame(RESULTS_STACKED)
RESULTS_STACKED = RESULTS_STACKED %>% dplyr::group_by(disease, geography, last_obs_time, h) %>% dplyr::mutate(error_scaled = error/min(error), 
                                                                                                              error_minmax = (error - min(error))/(max(error)-min(error)),
                                                                                                              mae_rank = rank(error, ties = 'min'),
                                                                                                              weights_rank = rank(1-weights, ties = 'min'))


A = RESULTS_STACKED[which.max(as.numeric(as.character(unlist(RESULTS_STACKED$mae_rank)))),]
RESULTS_STACKED = data.frame(RESULTS_STACKED)
B = RESULTS_STACKED[RESULTS_STACKED$disease == A$disease & RESULTS_STACKED$geography == A$geography & RESULTS_STACKED$last_obs_time == A$last_obs_time & RESULTS_STACKED$h == A$h,]
B = B[B$model == 'rw',]

RESULTS_STACKED$mae_rank = factor(RESULTS_STACKED$mae_rank, levels = rev(sort(unique(RESULTS_STACKED$mae_rank))))
p1 = ggplot(data = RESULTS_STACKED, aes(y=mae_rank, x = weights, group = mae_rank, fill = factor(stat(quantile))))+
  stat_density_ridges(
    geom = "density_ridges_gradient", calc_ecdf = TRUE,
    quantiles = 4, quantile_lines = TRUE,
  ) +
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=18,face="bold"), plot.title = element_text(size=22),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=15))+
  scale_fill_viridis_d(name = "Quartiles")+
  coord_cartesian(xlim = c(0,0.3))+
  xlab('Weights')+
  ylab('MAE Rank')+
  labs(title = 'Distribution of Epifforma Weights by Component Rank')

ridges_plot <- p1

table(RESULTS_STACKED$mae_rank, RESULTS_STACKED$weights_rank)

TO_PLOT = data.frame(nums = as.vector(unlist(table(RESULTS_STACKED$mae_rank, RESULTS_STACKED$weights_rank))),
                     mae_rank = rep(9:1, 9),
                     weights_rank = rep(1:9, each = 9))
TO_PLOT = TO_PLOT %>% dplyr::group_by(mae_rank) %>% dplyr::mutate(nums = nums/sum(nums))
TO_PLOT$mae_rank = factor(TO_PLOT$mae_rank, levels = sort(unique(TO_PLOT$mae_rank)))
TO_PLOT$weights_rank = factor(TO_PLOT$weights_rank, levels = rev(sort(unique(TO_PLOT$mae_rank))))
p1=ggplot(TO_PLOT)+
  geom_tile(aes(x=mae_rank, y = weights_rank, fill = nums, color = nums))+
  geom_text(aes(x=mae_rank, y = weights_rank, label = round(nums,2)), color = 'white', size = 10)+
  scale_fill_viridis(name = "Proportion\nof Column")+
  scale_color_viridis(name = "Proportion\nof Column")+
  theme_classic()+
  xlab('MAE Rank (1=Lowest)')+
  ylab('Weight Rank (1=Highest)')+
  theme(axis.text=element_text(size=15),
        axis.title=element_text(size=18,face="bold"), plot.title = element_text(size=22),
        legend.title=element_text(size=18), 
        legend.text=element_text(size=15))+
  labs(title = 'Correspondence between Epifforma Weight and MAE Rankings across Components')+
  scale_x_discrete(expand=c(0,0))+
  scale_y_discrete(expand=c(0,0))

rank_grid_plot <- p1

lst_p = list(rank_grid_plot, ridges_plot)

# pdf(file = paste0(my_path,'/figure2.pdf'), width = 16, height = 10)
grid.arrange(lst_p[[1]], lst_p[[2]], layout_matrix = matrix(c(1,2),ncol = 1))
# dev.off()

