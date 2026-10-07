# Originally Figures 3 & 4 in the supplement of sMOA.  Working to make comparison with OSFD
## Author: AC Murph + LJ Beesley
## Date: March 2026

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

output_path = here::here("data", "embeddings_gam_real")

FILES = list.files(output_path)

#######################
#######################
### FIGURE 3 ##########
#######################
#######################
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
    
    output = output %>% dplyr::mutate(mae_smoa_OSFD = mean(abs(truth-fcst),na.rm=T),
                                      mae_smoa_ISFD = mean(abs(truth-fcst_old),na.rm=T),
                                      mae_rw = mean(abs(truth-obs),na.rm=T),
                                      mse_smoa_OSFD = mean((truth-fcst)^2,na.rm=T),
                                      mse_smoa_ISFD = mean((truth-fcst_old)^2,na.rm=T),
                                      mse_rw = mean((truth-obs)^2,na.rm=T))
    # output = output[!duplicated(output$row_num),]
    
    SPLIT = strsplit(gsub('.csv','',gsub('real_eval_mat_','',FILES[i])), split = '_')[[1]]
    output$disease_source = paste0(SPLIT[2],'_',SPLIT[3])
    output$disease = SPLIT[2]
    output$N = nrow(output)
    output$FILES = FILES[i]
    output[,c('disease_source','disease','N','FILES','mae_smoa_OSFD','mae_smoa_ISFD','mae_rw','mse_smoa_OSFD','mse_smoa_ISFD','mse_rw')]
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

RESULTS_AGG = RESULTS_LONG %>% dplyr::group_by(disease_source) %>% dplyr::mutate(mae_smoa_agg_OSFD = mean(mae_smoa_OSFD),
                                                                                 rmse_smoa_agg_OSFD = sqrt(mean(mse_smoa_OSFD)),
                                                                                 mae_smoa_agg_ISFD = mean(mae_smoa_ISFD),
                                                                                 rmse_smoa_agg_ISFD = sqrt(mean(mse_smoa_ISFD)),
                                                                                 mae_rw_agg = mean(mae_rw),
                                                                                 rmse_rw_agg = sqrt(mean(mse_rw^2)),
                                                                                 N_agg = sum(N))
RESULTS_AGG = RESULTS_AGG[!duplicated(RESULTS_AGG$disease_source),]


RESULTS_AGG = RESULTS_AGG[order(RESULTS_AGG$mae_smoa_agg_OSFD/(RESULTS_AGG$mae_rw_agg+1e-6)),]
DISEASE_ORDER = rev(RESULTS_AGG$disease_source)
xx = RESULTS_LONG$mae_smoa_OSFD/(RESULTS_LONG$mae_smoa_ISFD+1e-6)
p1 = ggplot(RESULTS_LONG)+
  #geom_point(aes(x=factor(disease_source, levels = DISEASE_ORDER), y=mae_smoa/(mae_rw+1e-6)), alpha = 0.01)+
  geom_boxplot(aes(x=factor(disease_source, levels = DISEASE_ORDER), y=mae_smoa_OSFD/(mae_smoa_ISFD+1e-6)), outlier.size = 0.1, outlier.color = 'darkgray')+
  # geom_point(aes(x=factor(disease_source, levels = DISEASE_ORDER), y=mae_smoa_agg_OSFD/(mae_smoa_agg_ISFD+1e-6)), color = 'red', data = RESULTS_AGG, size = 2)+
  xlab('')+
  ylab('MAE(sMOA OSFD) / MAE(sMOA ISFD)')+
  geom_hline(yintercept = 1, color = 'gray', linetype = 2)+
  coord_flip(ylim=c(quantile(xx, probs = 0.01),quantile(xx, probs = 0.99)))+ 
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 15)) + 
  theme(plot.title = element_blank(),
        # axis.text.x = element_text(angle = 270, vjust = 0.5, hjust=1),
        panel.spacing = unit(0, "cm"),
        plot.margin = margin(0, 0, 0, 0, "cm"), 
        plot.caption = element_blank()) + 
    theme(
    axis.text = element_text(size = 18),
    axis.text.x = element_text(angle = 0, hjust = 1, size = 18),
    axis.title = element_text(size = 20),
    legend.title = element_blank(),
    legend.text = element_text(size = 18),
    plot.title = element_text(size = 22),
    plot.subtitle = element_text(size = 20)
  )
RESULTS_AGG = RESULTS_AGG[order(RESULTS_AGG$rmse_smoa_agg_OSFD/RESULTS_AGG$rmse_rw_agg),]
DISEASE_ORDER = rev(RESULTS_AGG$disease_source)

metric_plot = p1 

scale = 0.5
ggsave(
  filename = paste0(here::here("viz", "smoa_performances.png")),
  plot = metric_plot,
  width = 40*scale,      # adjust width as needed
  height = 20*scale,     # adjust height as needed
  dpi = 300       # high-quality resolution
)

#####################
### FIGURE 4 ########
#####################
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
  output = read.csv(paste0(output_path,"/",FILES[i]))
  # output$fcst = pmax(0,output$fcst)
  DAT = output
  # DAT = DAT  %>% dplyr::mutate(mae_smoa = mean(abs(truth-fcst),na.rm=T),
  #                                   mae_smoa_ISFD = mean(abs(truth-fcst_old),na.rm=T),
  #                                   mae_rw = mean(abs(truth-obs),na.rm=T)
  #                                 )
  SPLIT = strsplit(gsub('.csv','',gsub('real_eval_mat_','',FILES[i])), split = '_')[[1]]
  DAT$disease_source = paste0(SPLIT[2],'_',SPLIT[3])
  DAT$disease = SPLIT[2]
  DAT = DAT[DAT$h == 1,]
  DAT
  # print(paste0('Finished: ', i, ' of ', length(FILES)))
}
stopCluster(cl)

RESULTS_LONG = RESULTS

# RESULTS_LONG$rmse_smoa = sqrt(RESULTS_LONG$mse_smoa)
# RESULTS_LONG$rmse_rw= sqrt(RESULTS_LONG$mse_rw)

RESULTS_LONG$disease_source[RESULTS_LONG$disease_source == 'Chikungunya_deSouza'] = 'Chikungunya'
RESULTS_LONG$disease_source[RESULTS_LONG$disease_source == 'Influenza_ushhs'] = 'Influenza Hospitalizations'
RESULTS_LONG$disease_source[RESULTS_LONG$disease_source == 'Influenza_usflunet'] = 'ILI Incidence'
RESULTS_LONG$disease_source[RESULTS_LONG$disease_source == 'Dengue_opendengue'] = 'Dengue Fever'
RESULTS_LONG$disease_source[RESULTS_LONG$disease_source == 'COVID_jhuowid'] = 'COVID-19'
RESULTS_LONG$disease_source[RESULTS_LONG$disease_source == 'Mpox_who'] = 'Monkeypox'


# RESULTS_LONG$dist_ratio = pmin(1,RESULTS_LONG$min_dist)
# RESULTS_LONG$dist_ratio_ISFD = pmin(1,RESULTS_LONG$min_dist_ISFD)
RESULTS_LONG$dist_ratio_OSFD = RESULTS_LONG$min_dist_OSFD #/ (pmax(RESULTS_LONG$obs,0)+1e-5)
RESULTS_LONG$dist_ratio_ISFD = RESULTS_LONG$min_dist_ISFD #/ (pmax(RESULTS_LONG$obs,0)+1e-5)

#####################
# Create Supplement Figure 4:
RESULTS_LONG_SUB <- RESULTS_LONG[RESULTS_LONG$obs > 0, ]

RESULTS_AGG = RESULTS_LONG_SUB[order(RESULTS_LONG_SUB$dist_ratio_OSFD),]
DISEASE_ORDER = rev(c(
  "COVID-19",
  "ILI Incidence",
  "Dengue Fever",
  "Influenza Hospitalizations",
  "Chikungunya",
  "Monkeypox"
))

RESULTS_AGG <- RESULTS_LONG_SUB %>%
  mutate(
    disease_source = factor(disease_source, levels = DISEASE_ORDER)
  )

RESULTS_PLOT <- RESULTS_AGG %>%
  select(disease_source, dist_ratio_OSFD, dist_ratio_ISFD) %>%
  pivot_longer(
    cols = c(dist_ratio_OSFD, dist_ratio_ISFD),
    names_to = "distance_type",
    values_to = "distance_value"
  ) %>%
  mutate(
    distance_type = recode(
      distance_type,
      dist_ratio_OSFD = "OSFD",
      dist_ratio_ISFD = "ISFD"
    ),
    distance_type = factor(distance_type, levels = c("OSFD", "ISFD"))
  )

p1 <- ggplot(
  RESULTS_PLOT,
  aes(x = disease_source, y = distance_value, fill = distance_type)
) +
  geom_boxplot(
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  xlab("") +
  ylab("Minimum Distance Synthetic (Scaled)") +
  coord_flip() +
  scale_fill_manual(values = c("OSFD" = "steelblue", "ISFD" = "tomato")) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15),
    legend.title = element_blank(),
    legend.text = element_text(size = 13)
  )

p2 <- ggplot(
  RESULTS_PLOT,
  aes(x = disease_source, y = distance_value, fill = distance_type)
) +
  stat_summary(
    fun = mean,
    geom = "bar",
    position = position_dodge(width = 0.8),
    width = 0.7
  ) +
  xlab("") +
  ylab("Average Minimum Distance Synthetic (Scaled)") +
  scale_fill_manual(values = c(
                    "OSFD" = "#3B5F8A",   # muted blue
                    "ISFD" = "#E07A5F"    # warm terracotta
                    )) +
  scale_y_continuous(expand = c(0, 0)) +
  theme(
    axis.text = element_text(size = 15),
    axis.title = element_text(size = 15),
    legend.title = element_blank(),
    legend.text = element_text(size = 13)
  )

RESULTS_DIFF <- RESULTS_PLOT %>%
  group_by(disease_source, distance_type) %>%
  summarize(mean_distance = mean(distance_value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = distance_type,
    values_from = mean_distance
  ) %>%
  mutate(diff = ISFD - OSFD)

p <- ggplot(
  RESULTS_DIFF,
  aes(x = disease_source, y = diff, fill = diff > 0)
) +
  geom_col(width = 0.7) +
  xlab("") +
  ylab("Average Distance ISFD - Average Distance OSFD") +
  scale_fill_manual(
    values = c("TRUE" = "#3B5F8A", "FALSE" = "#E07A5F"),
    guide = "none"
  ) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  theme(
    axis.text = element_text(size = 18),
    axis.text.x = element_text(angle = 45, hjust = 1, size = 18),
    axis.title = element_text(size = 20),
    legend.title = element_blank(),
    legend.text = element_text(size = 18),
    plot.title = element_text(size = 22),
    plot.subtitle = element_text(size = 20)
  )

scale = 0.5
ggsave(
  filename = paste0(here::here("viz", "average_distances_to_synthetic.png")),
  plot = p,
  width = 40*scale,      # adjust width as needed
  height = 20*scale,     # adjust height as needed
  dpi = 300       # high-quality resolution
)
