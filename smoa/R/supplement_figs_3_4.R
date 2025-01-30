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
theme_set(theme_classic())

output_path = '~/Documents/evaluations/'

FILES = list.files(output_path)

#######################
#######################
### Aggregate Evals ###
#######################
#######################
RESULTS = data.frame(FILES = FILES)
RESULTS$disease = NA
RESULTS$disease_source = NA
RESULTS$N = NA
RESULTS$mae_smoa = NA
RESULTS$rmse_smoa = NA
RESULTS$mae_rw = NA
RESULTS$rmse_rw = NA
RESULTS$min_dist = NA
RESULTS$max_dist = NA
RESULTS$obs = NA
for(i in 1:length(FILES)){
  output = read.csv(paste0(output_path,FILES[i]))
  output$fcst = pmax(0,output$fcst)
  #output = output[output$obs > 0,]
  # if(!grepl('usflunet',FILES[i])){
  #   output$fcst[output$fcst<0.5]=0
  # }
  # output$fcst[!is.na(output$obs) & output$obs == 0] = 0
  SPLIT = strsplit(gsub('.csv','',gsub('real_eval_mat_','',FILES[i])), split = '_')[[1]]
  RESULTS$disease_source[i] = paste0(SPLIT[1],'_',SPLIT[2])
  RESULTS$disease[i] = SPLIT[1]
  RESULTS$N[i] = nrow(output)
  RESULTS$mae_smoa[i] = mean(abs(output$truth - output$fcst), na.rm=T)
  RESULTS$rmse_smoa[i] = sqrt(mean((output$truth - output$fcst)^2, na.rm=T))
  RESULTS$mae_rw[i] = mean(abs(output$truth - output$obs), na.rm=T)
  RESULTS$rmae_rw[i] = sqrt(mean((output$truth - output$obs)^2, na.rm=T))
  RESULTS$truth_avg[i] = mean(output$truth, na.rm=T)
  RESULTS$min_dist[i] = output$min_dist[1]
  RESULTS$max_dist[i] = output$max_dist[1]
  RESULTS$obs[i] = output$obs[1]
  print(paste0('Finished: ', i, ' of ', length(FILES)))
}

RESULTS_AGG = RESULTS %>% dplyr::group_by(disease_source) %>% dplyr::mutate(mae_smoa_agg = N*mae_smoa,
                                                                          rmse_smoa_agg = N*(rmse_smoa^2),
                                                                          mae_rw_agg = N*mae_rw,
                                                                          rmse_rw_agg = N*(rmse_rw^2),
                                                                          N_agg = sum(N))
RESULTS_AGG$mae_smoa_agg = RESULTS_AGG$mae_smoa_agg/RESULTS_AGG$N_agg
RESULTS_AGG$rmse_smoa_agg = sqrt(RESULTS_AGG$rmse_smoa_agg/RESULTS_AGG$N_agg)
RESULTS_AGG$mae_rw_agg = RESULTS_AGG$mae_rw_agg/RESULTS_AGG$N_agg
RESULTS_AGG$rmse_rw_agg = sqrt(RESULTS_AGG$rmse_rw_agg/RESULTS_AGG$N_agg)
RESULTS_AGG = RESULTS_AGG[!duplicated(RESULTS_AGG$disease_source),]


A = RESULTS[grepl('Chik',RESULTS$disease_source),]

RESULTS_AGG = RESULTS_AGG[order(RESULTS_AGG$mae_smoa_agg/(RESULTS_AGG$mae_rw_agg+1e-6)),]
DISEASE_ORDER = rev(RESULTS_AGG$disease_source)

#####################
# Create Supplement Figure 3 (first part)
ggplot(RESULTS)+
  geom_point(aes(x=factor(disease_source, levels = DISEASE_ORDER), y=mae_smoa/(mae_rw+1e-6)))+
  geom_point(aes(x=factor(disease_source, levels = DISEASE_ORDER), y=mae_smoa_agg/(mae_rw_agg+1e-6)), color = 'red', data = RESULTS_AGG, size = 2)+
  xlab('Disease_Source')+
  ylab('MAE(smoa) / MAE(random walk)')+
  geom_hline(yintercept = 1, color = 'gray', linetype = 2)+
  #ylim(0,1.5)+
  coord_flip()

RESULTS_AGG = RESULTS_AGG[order(RESULTS_AGG$rmse_smoa_agg/RESULTS_AGG$rmse_rw_agg),]
DISEASE_ORDER = rev(RESULTS_AGG$disease_source)

#####################
# Create Supplement Figure 3 (second part)
ggplot(RESULTS)+
  geom_point(aes(x=factor(disease_source, levels = DISEASE_ORDER), y=rmse_smoa/rmse_rw))+
  geom_point(aes(x=factor(disease_source, levels = DISEASE_ORDER), y=rmse_smoa_agg/rmse_rw_agg), color = 'red', data = RESULTS_AGG, size = 2)+
  xlab('Disease_Source')+
  ylab('rMSE(smoa) / rMSE(random walk)')+
  geom_hline(yintercept = 1, color = 'gray', linetype = 2)+
  #ylim(0,1.5)+
  coord_flip()


#####################
### Dist vs Error ###
#####################
RESULTS_LONG = NULL
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
  RESULTS_LONG = rbind(RESULTS_LONG,DAT)
  print(paste0('Finished: ', i, ' of ', length(FILES)))
}


RESULTS_LONG$ratio = RESULTS_LONG$mae_smoa/RESULTS_LONG$mae_rw#(RESULTS_LONG$mae_rw+1e-6)
RESULTS_LONG$dist_ratio = pmin(1,RESULTS_LONG$min_dist/(RESULTS_LONG$obs+1e-6))
RESULTS_LONG$mae_ratio2 = RESULTS_LONG$mae_smoa/(RESULTS_LONG$obs+1e-6)
RESULTS_LONG$mae_rw_ratio2 = RESULTS_LONG$mae_rw/(RESULTS_LONG$obs+1e-6)

#####################
# Create Supplement Figure 4:
ggplot(RESULTS_LONG[RESULTS_LONG$obs > 0,])+
  geom_boxplot(aes(x=disease_source, y = dist_ratio, fill = disease_source))+
  xlab('Disease')+
  ylab('Minimum Distance / Last Observed Value')+
  guides(fill = 'none')+
  scale_y_continuous(expand = c(0,0))+
  coord_flip()

