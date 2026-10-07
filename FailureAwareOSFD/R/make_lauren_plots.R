## Originally Figures 3 & 4 in the supplement of sMOA.  Working to make comparison with Failure-aware OSFD
## Author: LJ Beesley + AC Murph
## Date: August 2026

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
library(tidyr)
library(this.path)
library(patchwork)
library(gridExtra)
library(latex2exp)
library(parallel)
library(doParallel)
theme_set(theme_classic())


#######################
### Read in Results ###
#######################

output_path = here::here("data", "embeddings_gam_real")
# output_path = paste0("/Users/lbeesley/Desktop/mutantigen_output_expansion/smoa/data/", "embeddings_gam_real")

FILES = list.files(output_path)


## set up cluster
## define number of cores
ncores <- 90
## define socket type
sockettype <- "PSOCK"

cl <- parallel::makeCluster(spec = ncores,
                            type = sockettype)
setDefaultCluster(cl)
registerDoParallel(cl)

RESULTS <- foreach(i = 1:length(FILES), 
      .errorhandling = "pass",
      .combine = "rbind",
      .verbose = T,
      .packages = c('dplyr', 'tsfeatures', 'timeDate', 'lubridate','mgcv','zoo','forecast','collapse'))%dopar%{
  # for(i in 1:length(FILES)){
    SPLIT = strsplit(gsub('.csv','',gsub('real_eval_mat_','',FILES[i])), split = '_')[[1]]
    if(SPLIT[length(SPLIT)]!='mat') next
    output = read.csv(paste0(output_path,"/",FILES[i]))
    output$fcst = pmax(0,output$moa_OSFD)
    output$fcst_old = pmax(0, output$moa_ISFD)
    output$regime = output$regime_of_real_data
    output = output %>% dplyr::mutate(mae_smoa_failure_aware_osfd = mean(abs(truth-fcst),na.rm=T),
                                      mae_smoa_basic_lhs = mean(abs(truth-fcst_old),na.rm=T),
                                      mae_rw = mean(abs(truth-obs),na.rm=T),
                                      mse_smoa_failure_aware_osfd = mean((truth-fcst)^2,na.rm=T),
                                      mse_smoa_basic_lhs = mean((truth-fcst_old)^2,na.rm=T),
                                      mse_rw = mean((truth-obs)^2,na.rm=T), 
                                      mean_obs = mean(obs, na.rm=T), #obs_actual
                                      mean_truth = mean(truth, na.rm=T),
                                      )
    # output = output[!duplicated(output$row_num),]
    
    SPLIT = strsplit(gsub('.csv','',gsub('real_eval_mat_','',FILES[i])), split = '_')[[1]]
    output$disease_source = paste0(SPLIT[2],'_',SPLIT[3])
    output$disease = SPLIT[2]
    output$N = nrow(output)
    output$FILES = FILES[i]
    output[,c('disease_source','disease','N','FILES','mae_smoa_failure_aware_osfd','mae_smoa_basic_lhs','mae_rw','mse_smoa_failure_aware_osfd','mse_smoa_basic_lhs','mse_rw','min_dist_ISFD',
              'min_dist_OSFD','mean_obs', 'mean_truth','regime')]#,'mean_dist_ISFD','mean_dist_OSFD','regime','sd_to_match')]
    # print(paste0('Finished: ', i, ' of ', length(FILES)))
}
stopCluster(cl)

RESULTS_LONG = RESULTS
RESULTS_LONG$disease_source[RESULTS_LONG$disease_source == 'Chikungunya_deSouza'] = 'Chikungunya'
RESULTS_LONG$disease_source[RESULTS_LONG$disease_source == 'Influenza_ushhs'] = 'Influenza Hospitalizations'
RESULTS_LONG$disease_source[RESULTS_LONG$disease_source == 'Influenza_usflunet'] = 'ILI Incidence'
RESULTS_LONG$disease_source[RESULTS_LONG$disease_source == 'Dengue_opendengue'] = 'Dengue Fever'
RESULTS_LONG$disease_source[RESULTS_LONG$disease_source == 'COVID_jhuowid'] = 'COVID-19'
RESULTS_LONG$disease_source[RESULTS_LONG$disease_source == 'Mpox_who'] = 'Monkeypox'




########################
### Simple Summaries ###
########################
table(RESULTS_LONG$disease_source, RESULTS_LONG$regime)
prop.table(table(RESULTS_LONG$disease_source, RESULTS_LONG$regime), margin = 1)



##########################################
### Plotting Aggregates by Time Series ###
##########################################

RESULTS_SUB = RESULTS_LONG %>% dplyr::group_by(FILES) %>% dplyr::mutate(mean_dist_failure_aware_osfd = mean(min_dist_OSFD),
                                                                        mean_dist_basic_lhs = mean(min_dist_ISFD),
                                                                        mae_smoa_agg_failure_aware_osfd = mean(mae_smoa_failure_aware_osfd),
                                                                        mae_smoa_agg_basic_lhs = mean(mae_smoa_basic_lhs),
                                                                        mae_smoa_agg_rw = mean(mae_rw))
RESULTS_SUB = RESULTS_SUB[!duplicated(RESULTS_SUB$FILES),]

p1 = ggplot(RESULTS_SUB)+
  geom_boxplot(aes(x=disease_source, y = mean_dist_failure_aware_osfd/mean_dist_basic_lhs, group = disease_source, fill = disease_source), alpha = 0.5)+
  geom_hline(yintercept = 1)+
  xlab('')+
  ylab('Distance(Failure-aware OSFD) / Distance(Basic LHS)')+
  scale_fill_viridis(option = 'magma', discrete = T)+
  guides(fill = 'none')+
  coord_flip(ylim = c(0.7,1.1))+
  labs(title = 'Mean Minimum Distance to Library')

p2 = ggplot(RESULTS_SUB)+
  geom_boxplot(aes(x=disease_source, y = mae_smoa_agg_failure_aware_osfd/mae_smoa_agg_basic_lhs, group = disease_source, fill = disease_source), alpha = 0.5)+
  geom_hline(yintercept = 1)+
  xlab('')+
  ylab('MAE(Failure-aware OSFD) / MAE(Basic LHS)')+
  scale_fill_viridis(option = 'magma', discrete = T)+
  guides(fill = 'none')+
  coord_flip(ylim = c(0.85,1.1))+
  labs(title = 'Mean Absolute Prediction Error')

metric_plot = p1 | p2

scale = 0.5
ggsave(
  filename = here::here("viz", "performance_boxplots.png"),
  plot = metric_plot,
  width = 20*scale,      # adjust width as needed
  height = 8*scale,     # adjust height as needed
  dpi = 300       # high-quality resolution
)




### Note: These versions are uniformly worse than this version of random walk
### But note also that "obs" is actually a smoothed version of obs here.
### Could likely improve performance over random walk by increasing size of snippet library
p3 = ggplot(RESULTS_SUB)+
  geom_boxplot(aes(x=disease_source, y = mae_smoa_agg_failure_aware_osfd/mae_smoa_agg_rw, group = disease_source, fill = disease_source), alpha = 0.5)+
  geom_hline(yintercept = 1)+
  xlab('')+
  ylab('MAE(Failure-aware OSFD) / MAE(Persistence)')+
  scale_fill_viridis(option = 'magma', discrete = T)+
  guides(fill = 'none')+
  coord_flip(ylim = c(0.85,5))+
  labs(title = 'Failure-aware OSFD Performance vs. Persistance Model')

ggsave(
  filename = here::here("viz", "performance_rw.png"),
  plot = p3,
  width = 5,      # adjust width as needed
  height = 4,     # adjust height as needed
  dpi = 300       # high-quality resolution
)



###################################################
### Plotting Aggregates by Time Series x Regime ###
###################################################

RESULTS_SUB = RESULTS_LONG %>% dplyr::group_by(FILES, regime) %>% dplyr::mutate(mean_dist_failure_aware_osfd = mean(min_dist_OSFD),
                                                                        mean_dist_basic_lhs = mean(min_dist_ISFD),
                                                                        mae_smoa_agg_failure_aware_osfd = mean(mae_smoa_failure_aware_osfd),
                                                                        mae_smoa_agg_basic_lhs = mean(mae_smoa_basic_lhs),
                                                                        mae_smoa_agg_rw = mean(mae_rw))
RESULTS_SUB = RESULTS_SUB[!duplicated(paste0(RESULTS_SUB$FILES,'_',RESULTS_SUB$regime)),]


RESULTS_SUB$regime[RESULTS_SUB$regime == 'Dec'] = 'Decreasing'
RESULTS_SUB$regime[RESULTS_SUB$regime == 'Inc'] = 'Increasing'


p1 = ggplot(RESULTS_SUB)+
  geom_boxplot(aes(x=disease_source, y = mean_dist_failure_aware_osfd/mean_dist_basic_lhs, group = paste0(disease_source,'_',regime), fill = regime), alpha = 0.5)+
  geom_hline(yintercept = 1)+
  xlab('')+
  ylab('Distance(Failure-aware OSFD) / Distance(Basic LHS)')+
  scale_fill_viridis('',option = 'magma', discrete = T)+
  guides(fill = guide_legend(nrow = 2))+
  #guides(fill = 'none')+
  theme(legend.position = 'top')+
  coord_flip(ylim = c(0.5,1.5))+
  labs(title = 'Mean Minimum Distance to Library')

p2 = ggplot(RESULTS_SUB)+
  geom_boxplot(aes(x=disease_source, y = mae_smoa_agg_failure_aware_osfd/mae_smoa_agg_basic_lhs, group = paste0(disease_source,'_',regime), fill = regime), alpha = 0.5)+
  geom_hline(yintercept = 1)+
  xlab('')+
  ylab('MAE(Failure-aware OSFD) / MAE(Basic LHS)')+
  scale_fill_viridis('',option = 'magma', discrete = T)+
  guides(fill = guide_legend(nrow = 2))+
  theme(legend.position = 'top')+
  coord_flip(ylim = c(0.85,1.1))+
  labs(title = 'Mean Absolute Prediction Error')

metric_plot = p1 | p2

scale = 0.5
ggsave(
  filename = here::here("viz", "performance_boxplots_byregime.png"),
  plot = metric_plot,
  width = 21*scale,      # adjust width as needed
  height = 12*scale,     # adjust height as needed
  dpi = 300       # high-quality resolution
)




p3 = ggplot(RESULTS_SUB)+
  geom_boxplot(aes(x=disease_source, y = mae_smoa_agg_failure_aware_osfd/mean_obs, group = paste0(disease_source,'_',regime), fill = regime), alpha = 0.5, outliers = F)+
  xlab('')+
  ylab('MAE(Failure-aware OSFD) / Last Observed')+
  scale_fill_viridis('',option = 'magma', discrete = T)+
  theme(legend.position = 'top')+
  coord_flip()+
  labs(title = 'Mean Absolute Prediction Error')

#####################################
### Plotting Aggregates by Regime ###
#####################################

RESULTS_SUB = RESULTS_LONG %>% dplyr::group_by(disease_source, regime) %>% dplyr::mutate(mean_dist_failure_aware_osfd = mean(min_dist_OSFD),
                                                                                mean_dist_basic_lhs = mean(min_dist_ISFD),
                                                                                mae_smoa_agg_failure_aware_osfd = mean(N*mae_smoa_failure_aware_osfd)/sum(N),
                                                                                mae_smoa_agg_basic_lhs = mean(N*mae_smoa_basic_lhs)/sum(N),
                                                                                mae_smoa_agg_rw = mean(N*mae_rw)/sum(N),
                                                                                sum_N = sum(N))
RESULTS_SUB = RESULTS_SUB[!duplicated(paste0(RESULTS_SUB$disease_source,'_',RESULTS_SUB$regime)),]


p1 = ggplot(RESULTS_SUB)+
  geom_bar(aes(x = disease_source, y=100*(1-(mae_smoa_agg_failure_aware_osfd/mae_smoa_agg_basic_lhs)), group = paste0(disease_source,'_',regime), fill = regime), color = 'black',alpha = 0.5, stat = 'identity', position = position_dodge())+
  geom_hline(yintercept = 0)+
  xlab('')+
  ylab('% Reduction in MAE')+
  scale_fill_viridis('',option = 'magma', discrete = T)+
  theme(legend.position = 'top')+
  coord_cartesian(ylim = c(0,5))+
  labs(title = '% Reduction in Overall MAE')+
  scale_y_continuous(labels = function(x) x + 1) 


######################################
### Visualizing Library Embeddings ###
######################################

load(here::here("data", "synthetic_embeddings.RData"))

embed_mat_all = scale(rbind(embed_mat_X_OSFD, embed_mat_X_ISFD))
embed_mat_all2 = scale(rbind(cbind(embed_mat_y_OSFD), cbind(embed_mat_y_ISFD)))

embed_labels = c(rep('Failure-aware OSFD', nrow(embed_mat_X_OSFD)), rep('Basic LHS', nrow(embed_mat_X_ISFD)))
embed_regimes = c(regime_OSFD, regime_ISFD)

table(embed_labels, embed_regimes)
prop.table(table(embed_labels, embed_regimes), margin = 1)
chisq.test(table(embed_labels, embed_regimes))

for(u in 1:ncol(embed_mat_X_OSFD)){
  print(t.test(x=embed_mat_X_OSFD[,u], y = embed_mat_X_ISFD[,u])$p.value)
}
for(u in 1:ncol(embed_mat_y_OSFD)){
  print(t.test(x=embed_mat_y_OSFD[,u], y = embed_mat_y_ISFD[,u])$p.value)
}

p = replicate(length(unique(regime_OSFD)),list(NULL))
for( i in 1:length(unique(regime_OSFD))){
  regime = unique(regime_OSFD)[i]
  X_OSFD = embed_mat_X_OSFD[regime_OSFD == regime,]
  y_OSFD = embed_mat_y_OSFD[regime_OSFD == regime,]
  
  X_ISFD = embed_mat_X_ISFD[regime_ISFD == regime,]
  y_ISFD = embed_mat_y_ISFD[regime_ISFD == regime,]
  
  
  QUANT_X_OSFD = apply(X_OSFD, 2, FUN = function(x){quantile(x, probs = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))})
  QUANT_y_OSFD = apply(y_OSFD, 2, FUN = function(x){quantile(x, probs = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))})
  QUANT_X_ISFD = apply(X_ISFD, 2, FUN = function(x){quantile(x, probs = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))})
  QUANT_y_ISFD = apply(y_ISFD, 2, FUN = function(x){quantile(x, probs = c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95))})
  
  
  RES = rbind(data.frame(quant = rep(c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95), ncol(X_OSFD)),
                         value = as.vector(unlist(QUANT_X_OSFD)), analysis = 'OSFD', input = 'X',
                         x = rep(1:5, each = 7)),
              data.frame(quant = rep(c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95), ncol(y_OSFD)),
                         value = as.vector(unlist(QUANT_y_OSFD)), analysis = 'OSFD', input = 'y',
                         x = rep(5+1:4, each = 7)),
              data.frame(quant = rep(c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95), ncol(X_ISFD)),
                         value = as.vector(unlist(QUANT_X_ISFD)), analysis = 'ISFD', input = 'X',
                         x = rep(1:5, each = 7)),
              data.frame(quant = rep(c(0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95), ncol(y_ISFD)),
                         value = as.vector(unlist(QUANT_y_ISFD)), analysis = 'ISFD', input = 'y',
                         x = rep(5+1:4, each = 7)))
  
  p[[i]] = ggplot(RES)+
    geom_line(aes(x=x,y=value, group = paste0(quant,'_', analysis), color = quant, linetype = analysis))+
    xlab('Time')+
    ylab('Value')+
    coord_cartesian(ylim=c(-10,10))+
    labs(title = regime)+
    guides(color = 'none', linetype = 'none')
}

metric_plot = p[[1]]|p[[2]]|p[[4]]|p[[5]]




SUBSET = sample(1:nrow(embed_mat_all), nrow(embed_mat_all), replace = F)

mypca_diff <- prcomp(embed_mat_all[SUBSET, ],center = TRUE,scale. = TRUE)

# 2D PCA coordinates
mypca_diff_2d <- data.frame(mypca_diff$x[, 1:2])
colnames(mypca_diff_2d) = c('X1','X2')
mypca_diff_2d$source = embed_labels[SUBSET]
mypca_diff_2d$regime = embed_regimes[SUBSET]


mypca_diff_2d$regime[mypca_diff_2d$regime == 'Dec'] = 'Decreasing'
mypca_diff_2d$regime[mypca_diff_2d$regime == 'Inc'] = 'Increasing'

p1 = ggplot(mypca_diff_2d)+
  geom_point(aes(x=X1, y = X2, color = regime, shape = regime), alpha = 0.3)+
  scale_color_discrete('')+
  scale_shape_discrete('')+
  #scale_shape_manual('', values = c(0,1), breaks = c('OSFD', 'ISFD'))+
  theme_classic()+
  theme(legend.position = 'top')+
  xlab('Dimension 1')+
  ylab('Dimension 2')+
  #coord_cartesian(xlim=c(-50,50), ylim=c(-50,50))+
  facet_grid(.~source)+
  guides(color = guide_legend(override.aes = list(alpha = 1, size = 2)))+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  theme(
    panel.grid.major = element_line(color = "grey80"),
    panel.grid.minor = element_line(color = "grey90")
  )




ggsave(
  filename = here::here("viz", "pca.png"),
  plot = p1,
  width = 10,      # adjust width as needed
  height = 6,     # adjust height as needed
  dpi = 300       # high-quality resolution
)




