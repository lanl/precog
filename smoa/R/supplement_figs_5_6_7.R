## Dave Osthus
## 12-17-24
## view embeddings

library(ggplot2)
library(data.table)
library(umap)
library(dplyr)
library(ggExtra)
library(gridExtra)
library(this.path)
setwd(paste0(this.path::here(),"/../"))
theme_set(theme_bw())

savepath <- "data/embed_synthetic_w_data"
filepath_covid <- "data/forecasting_23-26/embed_synthetic_w_data/"
filepath_other <- "data/forecasting_23-26/embed_synthetic_w_data/embeddings_for_murph/"
load("data/forecasting_23-26/embed_synthetic_w_data/synthetic_X.RData")

## pick 100k rows of synthetic data
set.seed(75600)
nsample <- 100000
X <- t(apply(X,1,diff))
smallX <- data.frame(X[sample(1:nrow(X),nsample,replace=F),])
rm(X)

## plot first and second indices
## two groups clearly emerge (proportions and counts)
qplot(X1, X2, data=smallX, color=I("grey"))+
  scale_x_log10()+
  scale_y_log10()


## read in covid embeddings (by state)
statefiles <- list.files("data/forecasting_23-26/embed_synthetic_w_data/embeddings_by_location/")
statefiles <- statefiles[grep("_X.RData",statefiles)]
realX <- NULL
for(i in 1:length(statefiles)){
  print(i)
  load(paste0("data/forecasting_23-26/embed_synthetic_w_data/embeddings_by_location/",statefiles[i]))
  xx = data.frame(t(apply(state_X,1,diff)))
  state_X <- data.frame(state = 'COVID-19', xx, max_val = apply(xx,1,max))
  realX <- rbind(realX, state_X)
}

## read in Chikungunya
state_X <- read.csv("data/forecasting_23-26/embed_synthetic_w_data/embeddings_for_murph/embed_Chikungunya_deSouza_X.csv",
                    col.names = c('X0','X1','X2','X3','X4','X5'))
xx = data.frame(t(apply(state_X,1,diff)))
state_X <- data.frame(state = 'Chikungunya', xx, max_val = apply(xx,1,max))
realX <- rbind(realX, state_X)

## read in Dengue
state_X <- read.csv("data/forecasting_23-26/embed_synthetic_w_data/embeddings_for_murph/embed_Dengue_opendengue_X.csv",
                    col.names = c('X0','X1','X2','X3','X4','X5'))
xx = data.frame(t(apply(state_X,1,diff)))
state_X <- data.frame(state = 'Dengue', xx, max_val = apply(xx,1,max))
realX <- rbind(realX, state_X)

## read in Influenza (usflu)
state_X <- read.csv("data/forecasting_23-26/embed_synthetic_w_data/embeddings_for_murph/embed_Influenza_usflunet_X.csv",
                    col.names = c('X0','X1','X2','X3','X4','X5'))
xx = data.frame(t(apply(state_X,1,diff)))
state_X <- data.frame(state = 'Influenza (usflu)', xx, max_val = apply(xx,1,max))
realX <- rbind(realX, state_X)

## read in Influenza (ushhs)
state_X <- read.csv("data/forecasting_23-26/embed_synthetic_w_data/embeddings_for_murph/embed_Influenza_ushhs_X.csv",
                    col.names = c('X0','X1','X2','X3','X4','X5'))
xx = data.frame(t(apply(state_X,1,diff)))
state_X <- data.frame(state = 'Influenza (ushhs)', xx, max_val = apply(xx,1,max))
realX <- rbind(realX, state_X)

## read in Mpox
state_X <- read.csv("data/forecasting_23-26/embed_synthetic_w_data/embeddings_for_murph/embed_Mpox_who_X.csv",
                    col.names = c('X0','X1','X2','X3','X4','X5'))
xx = data.frame(t(apply(state_X,1,diff)))
state_X <- data.frame(state = 'Monkey Pox', xx, max_val = apply(xx,1,max))
realX <- rbind(realX, state_X)


# ## prepare for umap
umapX <- rbind(cbind(data.frame(state = "Synthetic"),smallX, max_val = apply(smallX,1,max)),
               data.frame(realX))
umapX$index = 1:nrow(umapX)

# fit umap to reduce down to 2 dimensions
set.seed(1128)
myumap <- umap(umapX[,-1], n_components = 2)

# ## prepare for plotting
umap_layout <- data.frame(myumap$layout)
umapdf <- data.frame(index = umapX$index, state=umapX$state, data.frame(umap_layout))
umapdf$type <- "real"
umapdf[umapdf$state == "Synthetic",]$type <- "synthetic"
umapdf$sum_snippet <- rowSums(umapX[,-1])
umapdf <- umapdf[order(umapdf$sum_snippet),]
umapdf$sum_snippet2 <- 1:nrow(umapdf)


umapdf$state = factor(umapdf$state, levels = unique(umapdf$state))
umapdf = as_tibble(umapdf)
save(umapdf, file='UMAP_data_alldiseases_diff.RData')
load(file='data/UMAP_data_alldiseases_diff.RData')

p_all_unbounded = ggplot(umapdf, aes(x=X1, y = X2, color = state)) + geom_point() + theme_bw() + 
  scale_color_manual(values = c('Synthetic'='#1b9e77', 
                                'COVID-19' = '#d95f02',
                                "Influenza (usflu)" = '#7570b3', 
                                "Chikungunya" = '#e7298a',
                                "Dengue" = '#66a61e',
                                "Influenza (ushhs)" = '#e6ab02',
                                "Monkey Pox" = '#a6761d')) 
ggsave(plot = p_all_unbounded, file = 'diff_figs/p_all_unbounded.png',
       width = 10, height = 10)

p_all = ggplot(umapdf, aes(x=X1, y = X2, color = state)) + geom_point() + theme_bw() + 
  ylim(-200,75) + xlim(-45,25) + 
  scale_color_manual(values = c('Synthetic'='#1b9e77', 
                                'COVID-19' = '#d95f02',
                                "Influenza (usflu)" = '#7570b3', 
                                "Chikungunya" = '#e7298a',
                                "Dengue" = '#66a61e',
                                "Influenza (ushhs)" = '#e6ab02',
                                "Monkey Pox" = '#a6761d'))
ggsave(plot = p_all, file = 'diff_figs/p_all.png',
       width = 10, height = 10)

# Investigate what the outlier data look like.
umapdf_outliers = umapdf
umapdf_outliers$outlier = ifelse( !(((umapdf_outliers$X1>= -25)&(umapdf_outliers$X1< 17.5))&
                                      ((umapdf_outliers$X2>= -25)&(umapdf_outliers$X2< 25))),
                                  1,0)
umapdf_outliers$in_tail = ifelse( ((umapdf_outliers$X1>-5)&(umapdf_outliers$X1<=7))&
                                    ((umapdf_outliers$X2> 0)&(umapdf_outliers$X2<=15)),
                                  1,0)
graph_data = umapdf_outliers
graph_data$UMAP1 = umapdf_outliers$X1
graph_data$UMAP2 = umapdf_outliers$X2
graph_data$X1 = NULL
graph_data$X2 = NULL

umapX$state=NULL
umapX = as_tibble(umapX)
umapdf_outliers = left_join(graph_data, umapX, by = 'index')
umapdf_outliers = umapdf_outliers[,-c(1,3,4,5,8, 9)]
umapdf_outliers_melted = melt(umapdf_outliers, id.vars = c('outlier', 'in_tail', 'state'))
umapdf_outliers_melted$X = c(rep(1, times = nrow(umapdf_outliers)),
                             rep(2, times = nrow(umapdf_outliers)),
                             rep(3, times = nrow(umapdf_outliers)),
                             rep(4, times = nrow(umapdf_outliers)),
                             rep(5, times = nrow(umapdf_outliers)))
umapdf_outliers_melted$data_num = rep(1:nrow(umapdf_outliers), times = 5)
for(data_type in unique(umapdf_outliers$state)){
  p1 = subset(graph_data, subset = (state==data_type)&(outlier==0) ) %>% 
    arrange(state) %>%
    ggplot( aes(x=UMAP1, y = UMAP2, color = state)) + geom_point(alpha = 0.7) + theme_bw() + 
    labs(color='Data') + 
    scale_color_manual(values = c('Synthetic'='#1b9e77', 
                                  'COVID-19' = '#d95f02',
                                  "Influenza (usflu)" = '#7570b3', 
                                  "Chikungunya" = '#e7298a',
                                  "Dengue" = '#66a61e',
                                  "Influenza (ushhs)" = '#e6ab02',
                                  "Monkey Pox" = '#a6761d')) +
    theme(legend.position = "none") + ggtitle(paste(data_type, ' Without Grid'))
  
  p2 = subset(graph_data, subset = (state==data_type)&(outlier==1) ) %>% 
    arrange(state) %>%
    ggplot( aes(x=UMAP1, y = UMAP2, color = state)) + geom_point(alpha = 0.7) + theme_bw() + 
    labs(color='Data') + 
    scale_color_manual(values = c('Synthetic'='#1b9e77', 
                                  'COVID-19' = '#d95f02',
                                  "Influenza (usflu)" = '#7570b3', 
                                  "Chikungunya" = '#e7298a',
                                  "Dengue" = '#66a61e',
                                  "Influenza (ushhs)" = '#e6ab02',
                                  "Monkey Pox" = '#a6761d')) +
    theme(legend.position = "none") + ggtitle(paste(data_type, ' Grid Values'))
  
  p3 = subset(graph_data, subset = (state==data_type)&(in_tail==1) ) %>% 
    arrange(state) %>%
    ggplot( aes(x=UMAP1, y = UMAP2, color = state)) + geom_point(alpha = 0.7) + theme_bw() + 
    labs(color='Data') + 
    scale_color_manual(values = c('Synthetic'='#1b9e77', 
                                  'COVID-19' = '#d95f02',
                                  "Influenza (usflu)" = '#7570b3', 
                                  "Chikungunya" = '#e7298a',
                                  "Dengue" = '#66a61e',
                                  "Influenza (ushhs)" = '#e6ab02',
                                  "Monkey Pox" = '#a6761d')) +
    theme(legend.position = "none") + ggtitle(paste(data_type, ' in tail'))
  
  p4 = subset(umapdf_outliers_melted, subset = (state==data_type)&(outlier==0) ) %>% 
    arrange(state) %>%
    ggplot( aes(x=X, y = value, color = data_num, group = data_num)) + geom_line(alpha = 0.7) + theme_bw() + 
    labs(color='Data') + 
    # scale_color_manual(values = c('Synthetic'='#1b9e77', 
    #                               'COVID-19' = '#d95f02',
    #                               "Influenza (usflu)" = '#7570b3', 
    #                               "Chikungunya" = '#e7298a',
    #                               "Dengue" = '#66a61e',
    #                               "Influenza (ushhs)" = '#e6ab02',
    #                               "Monkey Pox" = '#a6761d'))+
    theme(legend.position = "none") + ggtitle(paste(data_type, ' Without Grid'))
  
  p5 = subset(umapdf_outliers_melted, subset = (state==data_type)&(outlier==1) ) %>% 
    arrange(state) %>%
    ggplot( aes(x=X, y = value, color = data_num, group = data_num)) + geom_line(alpha = 0.7) + theme_bw() + 
    labs(color='Data') + 
    # scale_color_manual(values = c('Synthetic'='#1b9e77', 
    #                               'COVID-19' = '#d95f02',
    #                               "Influenza (usflu)" = '#7570b3', 
    #                               "Chikungunya" = '#e7298a',
    #                               "Dengue" = '#66a61e',
    #                               "Influenza (ushhs)" = '#e6ab02',
    #                               "Monkey Pox" = '#a6761d'))+
    theme(legend.position = "none") + ggtitle(paste(data_type, ' Grid Values'))
  
  p6 = subset(umapdf_outliers_melted, subset = (state==data_type)&(in_tail==1) ) %>% 
    arrange(state) %>%
    ggplot( aes(x=X, y = value, color = data_num, group = data_num)) + geom_line(alpha = 0.7) + theme_bw() + 
    labs(color='Data') + 
    # scale_color_manual(values = c('Synthetic'='#1b9e77', 
    #                               'COVID-19' = '#d95f02',
    #                               "Influenza (usflu)" = '#7570b3', 
    #                               "Chikungunya" = '#e7298a',
    #                               "Dengue" = '#66a61e',
    #                               "Influenza (ushhs)" = '#e6ab02',
    #                               "Monkey Pox" = '#a6761d')) +
    theme(legend.position = "none") + ggtitle(paste(data_type, ' in tail'))
  
  png(paste0('diff_figs/dashboard', data_type,'.png'),width = 480*3, height = 480*2,)
  grid.arrange(p1,p2,p3,p4,p5,p6, layout_matrix = matrix(1:6, ncol = 3, byrow = T))
  dev.off()
  
}

p_covid = subset(umapdf, subset = (state%in%c('Synthetic', 'COVID-19'))) %>% 
  arrange(state) %>%
  ggplot( aes(x=X1, y = X2, color = state)) + geom_point(alpha = 0.7) + theme_bw() + 
  labs(color='Data') + 
  scale_color_manual(values = c('Synthetic'='#1b9e77', 
                                'COVID-19' = '#d95f02',
                                "Influenza (usflu)" = '#7570b3', 
                                "Chikungunya" = '#e7298a',
                                "Dengue" = '#66a61e',
                                "Influenza (ushhs)" = '#e6ab02',
                                "Monkey Pox" = '#a6761d'))+ 
  ylim(-200,75) + xlim(-45, 25)
ggsave(plot = p_covid, file = 'diff_figs/p_covid.png',
       width = 10, height = 10)


p_covid_unbounded = subset(umapdf, subset = (state%in%c('Synthetic', 'COVID-19'))) %>% 
  arrange(state) %>%
  ggplot( aes(x=X1, y = X2, color = state)) + geom_point(alpha = 0.7) + theme_bw() + 
  labs(color='Data') + 
  scale_color_manual(values = c('Synthetic'='#1b9e77', 
                                'COVID-19' = '#d95f02',
                                "Influenza (usflu)" = '#7570b3', 
                                "Chikungunya" = '#e7298a',
                                "Dengue" = '#66a61e',
                                "Influenza (ushhs)" = '#e6ab02',
                                "Monkey Pox" = '#a6761d'))
ggsave(plot = p_covid_unbounded, file = 'diff_figs/p_covid_unbounded.png',
       width = 10, height = 10)

p_chik = subset(umapdf, subset = (state%in%c('Synthetic', 'Chikungunya'))) %>% 
  arrange(state) %>%
  ggplot( aes(x=X1, y = X2, color = state)) + geom_point(alpha = 0.7) + theme_bw() + 
  labs(color='Data') + 
  scale_color_manual(values = c('Synthetic'='#1b9e77', 
                                'COVID-19' = '#d95f02',
                                "Influenza (usflu)" = '#7570b3', 
                                "Chikungunya" = '#e7298a',
                                "Dengue" = '#66a61e',
                                "Influenza (ushhs)" = '#e6ab02',
                                "Monkey Pox" = '#a6761d'))+ 
  ylim(-200,75) + xlim(-45, 25)
ggsave(plot = p_chik, file = 'diff_figs/p_chik.png',
       width = 10, height = 10)

p_chik_unbounded = subset(umapdf, subset = (state%in%c('Synthetic', 'Chikungunya'))) %>% 
  arrange(state) %>%
  ggplot( aes(x=X1, y = X2, color = state)) + geom_point(alpha = 0.7) + theme_bw() + 
  labs(color='Data') + 
  scale_color_manual(values = c('Synthetic'='#1b9e77', 
                                'COVID-19' = '#d95f02',
                                "Influenza (usflu)" = '#7570b3', 
                                "Chikungunya" = '#e7298a',
                                "Dengue" = '#66a61e',
                                "Influenza (ushhs)" = '#e6ab02',
                                "Monkey Pox" = '#a6761d'))
ggsave(plot = p_chik_unbounded, file = 'diff_figs/p_chik_unbounded.png',
       width = 10, height = 10)

p_usflu = subset(umapdf, subset = (state%in%c('Synthetic', "Influenza (usflu)"))) %>% 
  arrange(state) %>%
  ggplot( aes(x=X1, y = X2, color = state)) + geom_point(alpha = 0.7) + theme_bw() + 
  labs(color='Data') + 
  scale_color_manual(values = c('Synthetic'='#1b9e77', 
                                'COVID-19' = '#d95f02',
                                "Influenza (usflu)" = '#7570b3', 
                                "Chikungunya" = '#e7298a',
                                "Dengue" = '#66a61e',
                                "Influenza (ushhs)" = '#e6ab02',
                                "Monkey Pox" = '#a6761d'))+ 
  ylim(-200,75) + xlim(-45, 25)
ggsave(plot = p_usflu, file = 'diff_figs/p_usflu.png',
       width = 10, height = 10)

p_usflu_unbounded = subset(umapdf, subset = (state%in%c('Synthetic', "Influenza (usflu)"))) %>% 
  arrange(state) %>%
  ggplot( aes(x=X1, y = X2, color = state)) + geom_point(alpha = 0.7) + theme_bw() + 
  labs(color='Data') + 
  scale_color_manual(values = c('Synthetic'='#1b9e77', 
                                'COVID-19' = '#d95f02',
                                "Influenza (usflu)" = '#7570b3', 
                                "Chikungunya" = '#e7298a',
                                "Dengue" = '#66a61e',
                                "Influenza (ushhs)" = '#e6ab02',
                                "Monkey Pox" = '#a6761d'))
ggsave(plot = p_usflu_unbounded, file = 'diff_figs/p_usflu_unbounded.png',
       width = 10, height = 10)

p_dengue = subset(umapdf, subset = (state%in%c('Synthetic', "Dengue"))) %>% 
  arrange(state) %>%
  ggplot( aes(x=X1, y = X2, color = state)) + geom_point(alpha = 0.7) + theme_bw() + 
  labs(color='Data') + 
  scale_color_manual(values = c('Synthetic'='#1b9e77', 
                                'COVID-19' = '#d95f02',
                                "Influenza (usflu)" = '#7570b3', 
                                "Chikungunya" = '#e7298a',
                                "Dengue" = '#66a61e',
                                "Influenza (ushhs)" = '#e6ab02',
                                "Monkey Pox" = '#a6761d'))+ 
  ylim(-200,75) + xlim(-45, 25)
ggsave(plot = p_dengue, file = 'diff_figs/p_dengue.png',
       width = 10, height = 10)

p_dengue_unbounded = subset(umapdf, subset = (state%in%c('Synthetic', "Dengue"))) %>% 
  arrange(state) %>%
  ggplot( aes(x=X1, y = X2, color = state)) + geom_point(alpha = 0.7) + theme_bw() + 
  labs(color='Data') + 
  scale_color_manual(values = c('Synthetic'='#1b9e77', 
                                'COVID-19' = '#d95f02',
                                "Influenza (usflu)" = '#7570b3', 
                                "Chikungunya" = '#e7298a',
                                "Dengue" = '#66a61e',
                                "Influenza (ushhs)" = '#e6ab02',
                                "Monkey Pox" = '#a6761d'))
ggsave(plot = p_dengue_unbounded, file = 'diff_figs/p_dengue_unbounded.png',
       width = 10, height = 10)

p_ushhs = subset(umapdf, subset = (state%in%c('Synthetic', "Influenza (ushhs)"))) %>% 
  arrange(state) %>%
  ggplot( aes(x=X1, y = X2, color = state)) + geom_point(alpha = 0.7) + theme_bw() + 
  labs(color='Data') + 
  scale_color_manual(values = c('Synthetic'='#1b9e77', 
                                'COVID-19' = '#d95f02',
                                "Influenza (usflu)" = '#7570b3', 
                                "Chikungunya" = '#e7298a',
                                "Dengue" = '#66a61e',
                                "Influenza (ushhs)" = '#e6ab02',
                                "Monkey Pox" = '#a6761d'))+ 
  ylim(-200,75) + xlim(-45, 25)
ggsave(plot = p_ushhs, file = 'diff_figs/p_ushhs.png',
       width = 10, height = 10)

p_ushhs_unbounded = subset(umapdf, subset = (state%in%c('Synthetic', "Influenza (ushhs)"))) %>% 
  arrange(state) %>%
  ggplot( aes(x=X1, y = X2, color = state)) + geom_point(alpha = 0.7) + theme_bw() + 
  labs(color='Data') + 
  scale_color_manual(values = c('Synthetic'='#1b9e77', 
                                'COVID-19' = '#d95f02',
                                "Influenza (usflu)" = '#7570b3', 
                                "Chikungunya" = '#e7298a',
                                "Dengue" = '#66a61e',
                                "Influenza (ushhs)" = '#e6ab02',
                                "Monkey Pox" = '#a6761d'))
ggsave(plot = p_ushhs_unbounded, file = 'diff_figs/p_ushhs_unbounded.png',
       width = 10, height = 10)

p_mpox = subset(umapdf, subset = (state%in%c('Synthetic', "Monkey Pox"))) %>% 
  arrange(state) %>%
  ggplot( aes(x=X1, y = X2, color = state)) + geom_point(alpha = 0.7) + theme_bw() + 
  labs(color='Data') + 
  scale_color_manual(values = c('Synthetic'='#1b9e77', 
                                'COVID-19' = '#d95f02',
                                "Influenza (usflu)" = '#7570b3', 
                                "Chikungunya" = '#e7298a',
                                "Dengue" = '#66a61e',
                                "Influenza (ushhs)" = '#e6ab02',
                                "Monkey Pox" = '#a6761d'))+ 
  ylim(-200,75) + xlim(-45, 25)
ggsave(plot = p_mpox, file = 'diff_figs/p_mpox.png',
       width = 10, height = 10)

p_mpox_unbounded = subset(umapdf, subset = (state%in%c('Synthetic', "Monkey Pox"))) %>% 
  arrange(state) %>%
  ggplot( aes(x=X1, y = X2, color = state)) + geom_point(alpha = 0.7) + theme_bw() + 
  labs(color='Data') + 
  scale_color_manual(values = c('Synthetic'='#1b9e77', 
                                'COVID-19' = '#d95f02',
                                "Influenza (usflu)" = '#7570b3', 
                                "Chikungunya" = '#e7298a',
                                "Dengue" = '#66a61e',
                                "Influenza (ushhs)" = '#e6ab02',
                                "Monkey Pox" = '#a6761d'))
ggsave(plot = p_mpox_unbounded, file = 'diff_figs/p_mpox_unbounded.png',
       width = 10, height = 10)





p_covid_unbounded = subset(umapdf, subset = (state%in%c('Synthetic', 'COVID-19'))) %>% 
  arrange(state) %>%
  ggplot( aes(x=X1, y = X2, color = state)) + geom_point(alpha = 0.7) + theme_bw() + 
  labs(color='Data') + 
  scale_color_manual(values = c('Synthetic'='#1b9e77', 
                                'COVID-19' = '#d95f02',
                                "Influenza (usflu)" = '#7570b3', 
                                "Chikungunya" = '#e7298a',
                                "Dengue" = '#66a61e',
                                "Influenza (ushhs)" = '#e6ab02',
                                "Monkey Pox" = '#a6761d'))
ggsave(plot = p_covid_unbounded, file = 'diff_figs/p_covid_unbounded.png',
       width = 10, height = 10)

p_chik_unbounded = subset(umapdf, subset = (state%in%c('Synthetic', 'Chikungunya'))) %>% 
  arrange(state) %>%
  ggplot( aes(x=X1, y = X2, color = state)) + geom_point(alpha = 0.7) + theme_bw() + 
  labs(color='Data') + 
  scale_color_manual(values = c('Synthetic'='#1b9e77', 
                                'COVID-19' = '#d95f02',
                                "Influenza (usflu)" = '#7570b3', 
                                "Chikungunya" = '#e7298a',
                                "Dengue" = '#66a61e',
                                "Influenza (ushhs)" = '#e6ab02',
                                "Monkey Pox" = '#a6761d'))
ggsave(plot = p_chik_unbounded, file = 'diff_figs/p_chik_unbounded.png',
       width = 10, height = 10)

p_usflu_unbounded = subset(umapdf, subset = (state%in%c('Synthetic', "Influenza (usflu)"))) %>% 
  arrange(state) %>%
  ggplot( aes(x=X1, y = X2, color = state)) + geom_point(alpha = 0.7) + theme_bw() + 
  labs(color='Data') + 
  scale_color_manual(values = c('Synthetic'='#1b9e77', 
                                'COVID-19' = '#d95f02',
                                "Influenza (usflu)" = '#7570b3', 
                                "Chikungunya" = '#e7298a',
                                "Dengue" = '#66a61e',
                                "Influenza (ushhs)" = '#e6ab02',
                                "Monkey Pox" = '#a6761d'))
ggsave(plot = p_usflu_unbounded, file = 'diff_figs/p_usflu_unbounded.png',
       width = 10, height = 10)

p_dengue_unbounded = subset(umapdf, subset = (state%in%c('Synthetic', "Dengue"))) %>% 
  arrange(state) %>%
  ggplot( aes(x=X1, y = X2, color = state)) + geom_point(alpha = 0.7) + theme_bw() + 
  labs(color='Data') + 
  scale_color_manual(values = c('Synthetic'='#1b9e77', 
                                'COVID-19' = '#d95f02',
                                "Influenza (usflu)" = '#7570b3', 
                                "Chikungunya" = '#e7298a',
                                "Dengue" = '#66a61e',
                                "Influenza (ushhs)" = '#e6ab02',
                                "Monkey Pox" = '#a6761d'))
ggsave(plot = p_dengue_unbounded, file = 'diff_figs/p_dengue_unbounded.png',
       width = 10, height = 10)

p_ushhs_unbounded = subset(umapdf, subset = (state%in%c('Synthetic', "Influenza (ushhs)"))) %>% 
  arrange(state) %>%
  ggplot( aes(x=X1, y = X2, color = state)) + geom_point(alpha = 0.7) + theme_bw() + 
  labs(color='Data') + 
  scale_color_manual(values = c('Synthetic'='#1b9e77', 
                                'COVID-19' = '#d95f02',
                                "Influenza (usflu)" = '#7570b3', 
                                "Chikungunya" = '#e7298a',
                                "Dengue" = '#66a61e',
                                "Influenza (ushhs)" = '#e6ab02',
                                "Monkey Pox" = '#a6761d'))
ggsave(plot = p_ushhs_unbounded, file = 'diff_figs/p_ushhs_unbounded.png',
       width = 10, height = 10)

p_mpox_unbounded = subset(umapdf, subset = (state%in%c('Synthetic', "Monkey Pox"))) %>% 
  arrange(state) %>%
  ggplot( aes(x=X1, y = X2, color = state)) + geom_point(alpha = 0.7) + theme_bw() + 
  labs(color='Data') + 
  scale_color_manual(values = c('Synthetic'='#1b9e77', 
                                'COVID-19' = '#d95f02',
                                "Influenza (usflu)" = '#7570b3', 
                                "Chikungunya" = '#e7298a',
                                "Dengue" = '#66a61e',
                                "Influenza (ushhs)" = '#e6ab02',
                                "Monkey Pox" = '#a6761d'))
ggsave(plot = p_mpox_unbounded, file = 'diff_figs/p_mpox_unbounded.png',
       width = 10, height = 10)


