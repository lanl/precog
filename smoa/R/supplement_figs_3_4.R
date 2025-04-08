# Create Figures 3 & 4 in the supplement
## Author: LJ Beesley
## Date: January 2025

## Note to researcher
# Data files are compiled using the get_real_embeddings_gam_internal.R file.

#############################
#############################
### Results Visualization ###
#############################
#############################
library(ggplot2)
library(data.table)
library(GGally)
library(viridis)
library(ggrepel)
library(dplyr)
library(this.path)
library(patchwork)
library(gridExtra)
setwd(paste0(this.path::here(),"/../"))
theme_set(theme_classic())

output_path = 'data/evaluations/'

FILES = list.files(output_path)

#######################
#######################
### FIGURE 3 ##########
#######################
#######################
RESULTS = NULL
for(i in 1:length(FILES)){
  output = read.csv(paste0(output_path,FILES[i]))
  output$fcst = pmax(0,output$fcst)
  
  output = output %>% dplyr::group_by(row_num) %>% dplyr::mutate(mae_smoa = mean(abs(truth-fcst),na.rm=T),
                                                                 mae_rw = mean(abs(truth-obs),na.rm=T),
                                                                 mse_smoa = mean((truth-fcst)^2,na.rm=T),
                                                                 mse_rw = mean((truth-obs)^2,na.rm=T))
  output = output[!duplicated(output$row_num),]
  
  SPLIT = strsplit(gsub('.csv','',gsub('real_eval_mat_','',FILES[i])), split = '_')[[1]]
  output$disease_source = paste0(SPLIT[1],'_',SPLIT[2])
  output$disease = SPLIT[1]
  output$N = nrow(output)
  output$FILES = FILES[i]
  RESULTS = rbind(RESULTS, output[,c('disease_source','disease','N','FILES','mae_smoa','mae_rw','mse_smoa','mse_rw')])
  print(paste0('Finished: ', i, ' of ', length(FILES)))
}

### Read in COVID Results from Murph
load(paste0(output_path,"../mae_comparison_data.RData"))
mae_comparison_data = mae_comparison_data[mae_comparison_data$model == 'COVIDhub-baseline',]
mae_comparison_data$obs = mae_comparison_data$abs_error_model + mae_comparison_data$true_values
mae_comparison_data$fcst_smoa = mae_comparison_data$abs_error_smoa + mae_comparison_data$true_values
mae_comparison_data$forecasts = pmax(0,mae_comparison_data$forecasts)
mae_comparison_data = mae_comparison_data %>% dplyr::group_by(fcast_date, location) %>% dplyr::mutate(mae_smoa = mean(abs(true_values-fcst_smoa),na.rm=T),
                                                                                                      mae_rw = mean(abs(true_values-obs),na.rm=T),
                                                                                                      mse_smoa = mean((true_values-fcst_smoa)^2,na.rm=T),
                                                                                                      mse_rw = mean((true_values-obs)^2,na.rm=T))
mae_comparison_data = mae_comparison_data[!duplicated(paste0(mae_comparison_data$fcast_date,'_',mae_comparison_data$location)),]
mae_comparison_data$disease = 'COVID'
mae_comparison_data$disease_source = 'COVID'
mae_comparison_data$N = NA
mae_comparison_data$FILES = NA
RESULTS_LONG = rbind(RESULTS, mae_comparison_data[,c('disease_source','disease','N','FILES','mae_smoa','mae_rw','mse_smoa','mse_rw')])

RESULTS_LONG$rmse_smoa = sqrt(RESULTS_LONG$mse_smoa)
RESULTS_LONG$rmse_rw= sqrt(RESULTS_LONG$mse_rw)

RESULTS_LONG$disease_source[RESULTS_LONG$disease_source == 'Chikungunya_deSouza'] = 'Chikungunya'
RESULTS_LONG$disease_source[RESULTS_LONG$disease_source == 'Influenza_ushhs'] = 'Influenza Hospitalizations'
RESULTS_LONG$disease_source[RESULTS_LONG$disease_source == 'Influenza_usflunet'] = 'ILI Incidence'
RESULTS_LONG$disease_source[RESULTS_LONG$disease_source == 'Dengue_opendengue'] = 'Dengue Fever'
RESULTS_LONG$disease_source[RESULTS_LONG$disease_source == 'COVID'] = 'COVID-19'

RESULTS_AGG = RESULTS_LONG %>% dplyr::group_by(disease_source) %>% dplyr::mutate(mae_smoa_agg = mean(mae_smoa),
                                                                                 rmse_smoa_agg = sqrt(mean(rmse_smoa^2)),
                                                                                 mae_rw_agg = mean(mae_rw),
                                                                                 rmse_rw_agg = sqrt(mean(rmse_rw^2)),
                                                                                 N_agg = sum(N))
RESULTS_AGG = RESULTS_AGG[!duplicated(RESULTS_AGG$disease_source),]



RESULTS_AGG = RESULTS_AGG[order(RESULTS_AGG$mae_smoa_agg/(RESULTS_AGG$mae_rw_agg+1e-6)),]
DISEASE_ORDER = rev(RESULTS_AGG$disease_source)
p1 = ggplot(RESULTS_LONG)+
  #geom_point(aes(x=factor(disease_source, levels = DISEASE_ORDER), y=mae_smoa/(mae_rw+1e-6)), alpha = 0.01)+
  geom_boxplot(aes(x=factor(disease_source, levels = DISEASE_ORDER), y=mae_smoa/(mae_rw+1e-6)), outlier.size = 0.1, outlier.color = 'darkgray')+
  geom_point(aes(x=factor(disease_source, levels = DISEASE_ORDER), y=mae_smoa_agg/(mae_rw_agg+1e-6)), color = 'red', data = RESULTS_AGG, size = 2)+
  xlab('')+
  ylab('MAE(sMOA) / MAE(baseline)')+
  geom_hline(yintercept = 1, color = 'gray', linetype = 2)+
  coord_flip(ylim=c(0,3))+ 
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 15)) + 
  theme(plot.title = element_blank(),
        # axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1),
        panel.spacing = unit(0, "cm"),
        plot.margin = margin(0, 0, 0, 0, "cm"), 
        plot.caption = element_blank()) + labs(fill=TeX("Nondim. Time [$y/\\tau$]"))

RESULTS_AGG = RESULTS_AGG[order(RESULTS_AGG$rmse_smoa_agg/RESULTS_AGG$rmse_rw_agg),]
DISEASE_ORDER = rev(RESULTS_AGG$disease_source)
p2 = ggplot(RESULTS_LONG)+
  geom_boxplot(aes(x=factor(disease_source, levels = DISEASE_ORDER), y=rmse_smoa/(rmse_rw+1e-6)), outlier.size = 0.1, outlier.color = 'darkgray')+
  geom_point(aes(x=factor(disease_source, levels = DISEASE_ORDER), y=rmse_smoa_agg/rmse_rw_agg), color = 'red', data = RESULTS_AGG, size = 2)+
  xlab('')+
  ylab('rMSE(sMOA) / rMSE(baseline)')+
  geom_hline(yintercept = 1, color = 'gray', linetype = 2)+
  coord_flip(ylim=c(0,3))+ 
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 15)) + 
  theme(plot.title = element_blank(),
        # axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1),
        panel.spacing = unit(0, "cm"),
        plot.margin = margin(0, 0, 0, 0, "cm"), 
        plot.caption = element_blank()) + labs(fill=TeX("Nondim. Time [$y/\\tau$]"))

p1 / p2

#####################
### FIGURE 4 ########
#####################
RESULTS = NULL
for(i in 1:length(FILES)){
  output = read.csv(paste0(output_path,FILES[i]))
  output$fcst = pmax(0,output$fcst)
  DAT = output
  DAT = DAT %>% dplyr::group_by(row_num) %>% dplyr::mutate(mae_smoa = mean(abs(truth - fcst), na.rm=T),
                                                           mae_rw= mean(abs(truth - obs), na.rm=T))
  SPLIT = strsplit(gsub('.csv','',gsub('real_eval_mat_','',FILES[i])), split = '_')[[1]]
  DAT$disease_source = paste0(SPLIT[1],'_',SPLIT[2])
  DAT$disease = SPLIT[1]
  DAT = DAT[DAT$h == 1,]
  RESULTS = rbind(RESULTS,DAT)
  print(paste0('Finished: ', i, ' of ', length(FILES)))
}

### Read in COVID Results from Murph
load(paste0("data/mae_comparison_data.RData"))
mae_comparison_data = mae_comparison_data[mae_comparison_data$model == 'COVIDhub-baseline',]
mae_comparison_data$obs = mae_comparison_data$abs_error_model + mae_comparison_data$true_values
mae_comparison_data$fcst_smoa = mae_comparison_data$abs_error_smoa + mae_comparison_data$true_values
mae_comparison_data$forecasts = pmax(0,mae_comparison_data$forecasts)
mae_comparison_data = mae_comparison_data %>% dplyr::group_by(fcast_date, location) %>% dplyr::mutate(mae_smoa = mean(abs(true_values-fcst_smoa),na.rm=T),
                                                                                                      mae_rw = mean(abs(true_values-obs),na.rm=T),
                                                                                                      mse_smoa = mean((true_values-fcst_smoa)^2,na.rm=T),
                                                                                                      mse_rw = mean((true_values-obs)^2,na.rm=T))
mae_comparison_data = mae_comparison_data[!duplicated(paste0(mae_comparison_data$fcast_date,'_',mae_comparison_data$location)),]
mae_comparison_data$disease = 'COVID'
mae_comparison_data$disease_source = 'COVID'
mae_comparison_data$N = NA
mae_comparison_data$FILES = NA
RESULTS_LONG = rbind(RESULTS, mae_comparison_data[,c('disease_source','disease','N','FILES','mae_smoa','mae_rw','mse_smoa','mse_rw', 'obs', 'min_dist')])

RESULTS_LONG$rmse_smoa = sqrt(RESULTS_LONG$mse_smoa)
RESULTS_LONG$rmse_rw= sqrt(RESULTS_LONG$mse_rw)

RESULTS_LONG$disease_source[RESULTS_LONG$disease_source == 'Chikungunya_deSouza'] = 'Chikungunya'
RESULTS_LONG$disease_source[RESULTS_LONG$disease_source == 'Influenza_ushhs'] = 'Influenza Hospitalizations'
RESULTS_LONG$disease_source[RESULTS_LONG$disease_source == 'Influenza_usflunet'] = 'ILI Incidence'
RESULTS_LONG$disease_source[RESULTS_LONG$disease_source == 'Dengue_opendengue'] = 'Dengue Fever'
RESULTS_LONG$disease_source[RESULTS_LONG$disease_source == 'COVID'] = 'COVID-19'


RESULTS_LONG$ratio = RESULTS_LONG$mae_smoa/RESULTS_LONG$mae_rw#(RESULTS_LONG$mae_rw+1e-6)
RESULTS_LONG$dist_ratio = pmin(1,RESULTS_LONG$min_dist/(RESULTS_LONG$obs+1e-6))
RESULTS_LONG$mae_ratio2 = RESULTS_LONG$mae_smoa/(RESULTS_LONG$obs+1e-6)
RESULTS_LONG$mae_rw_ratio2 = RESULTS_LONG$mae_rw/(RESULTS_LONG$obs+1e-6)

#####################
# Create Supplement Figure 4:
RESULTS_LONG_SUB <- RESULTS_LONG[RESULTS_LONG$obs > 0, ]

RESULTS_AGG = RESULTS_LONG_SUB[order(RESULTS_LONG_SUB$dist_ratio),]
DISEASE_ORDER = rev(c("COVID-19", "ILI Incidence", "Dengue Fever", "Influenza Hospitalizations", "Chikungunya"))
  
RESULTS_AGG$disease_source <- factor(
  RESULTS_AGG$disease_source,
  levels = DISEASE_ORDER
)

# Plot
ggplot(RESULTS_AGG) +
  geom_boxplot(aes(x = disease_source, y = dist_ratio), fill = 'white') +
  xlab("") +
  ylab("Minimum Distance / Last Observed Value") +
  guides(fill = 'none') +
  scale_x_discrete(labels = unique(DISEASE_ORDER))+
  scale_y_continuous(expand = c(0, 0)) +
  coord_flip()+ 
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 15)) 


