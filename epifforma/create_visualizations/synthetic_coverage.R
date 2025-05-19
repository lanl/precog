# Examine how well the synthetic data represent features in epifforma.
## Author: DA Osthus, LJ Beesley, and AC Murph
## Date: April 2025
######################
### Body of Script ###
######################
library(ggplot2)
library(data.table)
library(GGally)
library(viridis)
library(dplyr)
library(ggrepel)
library(this.path)
library(gridExtra)
library(ggExtra)
library(grid)
library(umap)
library(stats)
library(this.path)
setwd(paste0(this.path::here(), '/../'))
theme_set(theme_classic())

screen_multivariate_outliers <- function(df, threshold = 0.975) {
  df_numeric <- df[sapply(df, is.numeric)]
  
  center <- colMeans(df_numeric, na.rm = TRUE)
  cov_matrix <- cov(df_numeric, use = "complete.obs")
  
  dists <- mahalanobis(df_numeric, center, cov_matrix)
  cutoff <- qchisq(threshold, df = ncol(df_numeric))
  
  return(df_numeric[dists < cutoff,])
}


my_path = paste0(this.path::here(), '/../')
savepath <- paste0(my_path,"evaluate_model/evaluation/")

FILES = list.files(savepath)
outlier_threshold = 0.85

############################
### Read in and Organize ###
############################

DISEASE_KEY = 'Dengue'
### demonstrating using Dengue
output_list = readRDS(paste0(savepath,FILES[grepl(DISEASE_KEY,FILES, ignore.case = T)]))


sim_list = readRDS(paste0(savepath,FILES[grepl('synthetic_training',FILES)]))


### organize features
FEATURE_LIST = setdiff(names(output_list['features']$features),c('geography','disease','last_obs_time','h'))
features_wide = output_list['features']$features


train_data = sim_list['features']$features

rm(output_list)


#########################################
### UMAP Analysis ###
#########################################
train_sim_data = train_data %>%
  # dplyr::filter(h==1) %>% 
  dplyr::select(gr12_div_23, last_div_max, coefvar, 
                gam_with_div_without, avg_recent_div_avg_global,
                diff_zscore, entropy, relative_increases,
                prop_since_peak, seasonality)

# Note I might want to downsample sim data.
set.seed(75600)
nsample <- 10000
set.seed(534300)
# train_sim_data = screen_multivariate_outliers(train_sim_data, threshold = 0.99)
train_sim_data = train_sim_data[sample(1:nrow(train_sim_data),nsample,replace=F),]

######################
## Load Global_Covid
DISEASE_KEY = 'Global_Covid'
### demonstrating using Global_Covid
output_list = readRDS(paste0(savepath,FILES[grepl(DISEASE_KEY,FILES, ignore.case = T)]))

### organize features
FEATURE_LIST = setdiff(names(output_list['features']$features),c('geography','disease','last_obs_time','h'))
features_wide = output_list['features']$features
train_disease_data = features_wide %>%
  # dplyr::filter(h==1) %>%
  dplyr::select(gr12_div_23, last_div_max, coefvar, 
                gam_with_div_without, avg_recent_div_avg_global,
                diff_zscore, entropy, relative_increases,
                prop_since_peak, seasonality)
if(nrow(train_disease_data) > nsample){
  set.seed(9878)
  train_disease_data <- train_disease_data[sample(1:nrow(train_disease_data), nsample, replace = F),]
}
# train_disease_data = screen_multivariate_outliers(train_disease_data, threshold = outlier_threshold)

## combine with synthetic
global_coviddf <- as.matrix(rbind(train_sim_data, train_disease_data))

## fit UMAP (2-dimensions)
set.seed(1128)
global_covid_umap <- umap(global_coviddf, n_components = 2)
global_covid_umap_df <- data.frame(global_covid_umap$layout)
names(global_covid_umap_df) <- paste0("X",1:ncol(global_covid_umap_df))
global_covid_umap_df = screen_multivariate_outliers(global_covid_umap_df, threshold = outlier_threshold)
global_covid_umap_df$type <- DISEASE_KEY
global_covid_umap_df[1:nrow(train_sim_data),]$type <- "Synthetic"

######################
## Load US_Covid
DISEASE_KEY = 'US_Covid'
### demonstrating using US_Covid
output_list = readRDS(paste0(savepath,FILES[grepl(DISEASE_KEY,FILES, ignore.case = T)]))

### organize features
FEATURE_LIST = setdiff(names(output_list['features']$features),c('geography','disease','last_obs_time','h'))
features_wide = output_list['features']$features
train_disease_data = features_wide %>%
  # dplyr::filter(h==1) %>%
  dplyr::select(gr12_div_23, last_div_max, coefvar, 
                gam_with_div_without, avg_recent_div_avg_global,
                diff_zscore, entropy, relative_increases,
                prop_since_peak, seasonality)
if(nrow(train_disease_data) > nsample){
  set.seed(9878)
  train_disease_data <- train_disease_data[sample(1:nrow(train_disease_data), nsample, replace = F),]
}
# train_disease_data = screen_multivariate_outliers(train_disease_data, threshold = outlier_threshold)

## combine with synthetic
us_coviddf <- as.matrix(rbind(train_sim_data, train_disease_data))

## fit UMAP (2-dimensions)
set.seed(1128)
us_covid_umap <- umap(us_coviddf, n_components = 2)
us_covid_umap_df <- data.frame(us_covid_umap$layout)
names(us_covid_umap_df) <- paste0("X",1:ncol(us_covid_umap_df))
us_covid_umap_df = screen_multivariate_outliers(us_covid_umap_df, threshold = outlier_threshold)
us_covid_umap_df$type <- DISEASE_KEY
us_covid_umap_df[1:nrow(train_sim_data),]$type <- "Synthetic"

######################
## Load ILI
DISEASE_KEY = 'ILI'
### demonstrating using ILI
output_list = readRDS(paste0(savepath,FILES[grepl(DISEASE_KEY,FILES, ignore.case = T)]))

### organize features
FEATURE_LIST = setdiff(names(output_list['features']$features),c('geography','disease','last_obs_time','h'))
features_wide = output_list['features']$features
train_disease_data = features_wide %>%
  # dplyr::filter(h==1) %>%
  dplyr::select(gr12_div_23, last_div_max, coefvar, 
                gam_with_div_without, avg_recent_div_avg_global,
                diff_zscore, entropy, relative_increases,
                prop_since_peak, seasonality)
if(nrow(train_disease_data) > nsample){
  set.seed(9878)
  train_disease_data <- train_disease_data[sample(1:nrow(train_disease_data), nsample, replace = F),]
}
# train_disease_data = screen_multivariate_outliers(train_disease_data, threshold = outlier_threshold)

## combine with synthetic
ilidf <- as.matrix(rbind(train_sim_data, train_disease_data))

## fit UMAP (2-dimensions)
set.seed(1128)
ili_umap <- umap(ilidf, n_components = 2)
ili_umap_df <- data.frame(ili_umap$layout)
names(ili_umap_df) <- paste0("X",1:ncol(ili_umap_df))
ili_umap_df = screen_multivariate_outliers(ili_umap_df, threshold = outlier_threshold)
ili_umap_df$type <- DISEASE_KEY
ili_umap_df[1:nrow(train_sim_data),]$type <- "Synthetic"


######################
## Load Rubella
DISEASE_KEY = 'Rubella'
### demonstrating using Rubella
output_list = readRDS(paste0(savepath,FILES[grepl(DISEASE_KEY,FILES, ignore.case = T)]))

### organize features
FEATURE_LIST = setdiff(names(output_list['features']$features),c('geography','disease','last_obs_time','h'))
features_wide = output_list['features']$features
train_disease_data = features_wide %>%
  # dplyr::filter(h==1) %>%
  dplyr::select(gr12_div_23, last_div_max, coefvar, 
                gam_with_div_without, avg_recent_div_avg_global,
                diff_zscore, entropy, relative_increases,
                prop_since_peak, seasonality)
if(nrow(train_disease_data) > nsample){
  set.seed(9878)
  train_disease_data <- train_disease_data[sample(1:nrow(train_disease_data), nsample, replace = F),]
}
# train_disease_data = screen_multivariate_outliers(train_disease_data, threshold = outlier_threshold)


## combine with synthetic
rubelladf <- as.matrix(rbind(train_sim_data, train_disease_data))

## fit UMAP (2-dimensions)
set.seed(1128)
rubella_umap <- umap(rubelladf, n_components = 2)
rubella_umap_df <- data.frame(rubella_umap$layout)
names(rubella_umap_df) <- paste0("X",1:ncol(rubella_umap_df))
rubella_umap_df = screen_multivariate_outliers(rubella_umap_df, threshold = outlier_threshold)
rubella_umap_df$type <- DISEASE_KEY
rubella_umap_df[1:nrow(train_sim_data),]$type <- "Synthetic"

######################
## Load Measles
DISEASE_KEY = 'Measles'
### demonstrating using Measles
output_list = readRDS(paste0(savepath,FILES[grepl(DISEASE_KEY,FILES, ignore.case = T)]))

### organize features
FEATURE_LIST = setdiff(names(output_list['features']$features),c('geography','disease','last_obs_time','h'))
features_wide = output_list['features']$features
train_disease_data = features_wide %>%
  # dplyr::filter(h==1) %>%
  dplyr::select(gr12_div_23, last_div_max, coefvar, 
                gam_with_div_without, avg_recent_div_avg_global,
                diff_zscore, entropy, relative_increases,
                prop_since_peak, seasonality)
if(nrow(train_disease_data) > nsample){
  set.seed(9878)
  train_disease_data <- train_disease_data[sample(1:nrow(train_disease_data), nsample, replace = F),]
}
# train_disease_data = screen_multivariate_outliers(train_disease_data, threshold = outlier_threshold)


## combine with synthetic
measlesdf <- as.matrix(rbind(train_sim_data, train_disease_data))

## fit UMAP (2-dimensions)
set.seed(1128)
measles_umap <- umap(measlesdf, n_components = 2)
measles_umap_df <- data.frame(measles_umap$layout)
names(measles_umap_df) <- paste0("X",1:ncol(measles_umap_df))
measles_umap_df = screen_multivariate_outliers(measles_umap_df, threshold = outlier_threshold)
measles_umap_df$type <- DISEASE_KEY
measles_umap_df[1:nrow(train_sim_data),]$type <- "Synthetic"


######################
## Load Smallpox
DISEASE_KEY = 'Smallpox'
### demonstrating using Smallpox
output_list = readRDS(paste0(savepath,FILES[grepl(DISEASE_KEY,FILES, ignore.case = T)]))

### organize features
FEATURE_LIST = setdiff(names(output_list['features']$features),c('geography','disease','last_obs_time','h'))
features_wide = output_list['features']$features
train_disease_data = features_wide %>%
  # dplyr::filter(h==1) %>%
  dplyr::select(gr12_div_23, last_div_max, coefvar, 
                gam_with_div_without, avg_recent_div_avg_global,
                diff_zscore, entropy, relative_increases,
                prop_since_peak, seasonality)
if(nrow(train_disease_data) > nsample){
  set.seed(9878)
  train_disease_data <- train_disease_data[sample(1:nrow(train_disease_data), nsample, replace = F),]
}
# train_disease_data = screen_multivariate_outliers(train_disease_data, threshold = outlier_threshold)


## combine with synthetic
smallpoxdf <- as.matrix(rbind(train_sim_data, train_disease_data))

## fit UMAP (2-dimensions)
set.seed(1128)
smallpox_umap <- umap(smallpoxdf, n_components = 2)
smallpox_umap_df <- data.frame(smallpox_umap$layout)
names(smallpox_umap_df) <- paste0("X",1:ncol(smallpox_umap_df))
smallpox_umap_df = screen_multivariate_outliers(smallpox_umap_df, threshold = outlier_threshold)
smallpox_umap_df$type <- DISEASE_KEY
smallpox_umap_df[1:nrow(train_sim_data),]$type <- "Synthetic"


######################
## Load Chikungunya
DISEASE_KEY = 'Chikungunya'
### demonstrating using Chikungunya
output_list = readRDS(paste0(savepath,FILES[grepl(DISEASE_KEY,FILES, ignore.case = T)]))

### organize features
FEATURE_LIST = setdiff(names(output_list['features']$features),c('geography','disease','last_obs_time','h'))
features_wide = output_list['features']$features
train_disease_data = features_wide %>%
  # dplyr::filter(h==1) %>%
  dplyr::select(gr12_div_23, last_div_max, coefvar, 
                gam_with_div_without, avg_recent_div_avg_global,
                diff_zscore, entropy, relative_increases,
                prop_since_peak, seasonality)
if(nrow(train_disease_data) > nsample){
  set.seed(9878)
  train_disease_data <- train_disease_data[sample(1:nrow(train_disease_data), nsample, replace = F),]
}
# train_disease_data = screen_multivariate_outliers(train_disease_data, threshold = outlier_threshold)


## combine with synthetic
chikungunyadf <- as.matrix(rbind(train_sim_data, train_disease_data))

## fit UMAP (2-dimensions)
set.seed(1128)
chikungunya_umap <- umap(chikungunyadf, n_components = 2)
chikungunya_umap_df <- data.frame(chikungunya_umap$layout)
names(chikungunya_umap_df) <- paste0("X",1:ncol(chikungunya_umap_df))
chikungunya_umap_df = screen_multivariate_outliers(chikungunya_umap_df, threshold = outlier_threshold)
chikungunya_umap_df$type <- DISEASE_KEY
chikungunya_umap_df[1:nrow(train_sim_data),]$type <- "Synthetic"


######################
## Load Diphtheria
DISEASE_KEY = 'Diphtheria'
### demonstrating using Diphtheria
output_list = readRDS(paste0(savepath,FILES[grepl(DISEASE_KEY,FILES, ignore.case = T)]))

### organize features
FEATURE_LIST = setdiff(names(output_list['features']$features),c('geography','disease','last_obs_time','h'))
features_wide = output_list['features']$features
train_disease_data = features_wide %>%
  # dplyr::filter(h==1) %>%
  dplyr::select(gr12_div_23, last_div_max, coefvar, 
                gam_with_div_without, avg_recent_div_avg_global,
                diff_zscore, entropy, relative_increases,
                prop_since_peak, seasonality)
if(nrow(train_disease_data) > nsample){
  set.seed(9878)
  train_disease_data <- train_disease_data[sample(1:nrow(train_disease_data), nsample, replace = F),]
}
# train_disease_data = screen_multivariate_outliers(train_disease_data, threshold = outlier_threshold)


## combine with synthetic
diphtheriadf <- as.matrix(rbind(train_sim_data, train_disease_data))

## fit UMAP (2-dimensions)
set.seed(1128)
diphtheria_umap <- umap(diphtheriadf, n_components = 2)
diphtheria_umap_df <- data.frame(diphtheria_umap$layout)
names(diphtheria_umap_df) <- paste0("X",1:ncol(diphtheria_umap_df))
diphtheria_umap_df = screen_multivariate_outliers(diphtheria_umap_df, threshold = outlier_threshold)
diphtheria_umap_df$type <- DISEASE_KEY
diphtheria_umap_df[1:nrow(train_sim_data),]$type <- "Synthetic"


######################
## Load Polio
DISEASE_KEY = 'Polio'
### demonstrating using Polio
output_list = readRDS(paste0(savepath,FILES[grepl(DISEASE_KEY,FILES, ignore.case = T)]))

### organize features
FEATURE_LIST = setdiff(names(output_list['features']$features),c('geography','disease','last_obs_time','h'))
features_wide = output_list['features']$features
train_disease_data = features_wide %>%
  # dplyr::filter(h==1) %>%
  dplyr::select(gr12_div_23, last_div_max, coefvar, 
                gam_with_div_without, avg_recent_div_avg_global,
                diff_zscore, entropy, relative_increases,
                prop_since_peak, seasonality)
if(nrow(train_disease_data) > nsample){
  set.seed(9878)
  train_disease_data <- train_disease_data[sample(1:nrow(train_disease_data), nsample, replace = F),]
}
# train_disease_data = screen_multivariate_outliers(train_disease_data, threshold = outlier_threshold)


## combine with synthetic
poliodf <- as.matrix(rbind(train_sim_data, train_disease_data))

## fit UMAP (2-dimensions)
set.seed(1128)
polio_umap <- umap(poliodf, n_components = 2)
polio_umap_df <- data.frame(polio_umap$layout)
names(polio_umap_df) <- paste0("X",1:ncol(polio_umap_df))
polio_umap_df = screen_multivariate_outliers(polio_umap_df, threshold = outlier_threshold)
polio_umap_df$type <- DISEASE_KEY
polio_umap_df[1:nrow(train_sim_data),]$type <- "Synthetic"


######################
## Load Mumps
DISEASE_KEY = 'Mumps'
### demonstrating using Mumps
output_list = readRDS(paste0(savepath,FILES[grepl(DISEASE_KEY,FILES, ignore.case = T)]))

### organize features
FEATURE_LIST = setdiff(names(output_list['features']$features),c('geography','disease','last_obs_time','h'))
features_wide = output_list['features']$features
train_disease_data = features_wide %>%
  # dplyr::filter(h==1) %>%
  dplyr::select(gr12_div_23, last_div_max, coefvar, 
                gam_with_div_without, avg_recent_div_avg_global,
                diff_zscore, entropy, relative_increases,
                prop_since_peak, seasonality)
if(nrow(train_disease_data) > nsample){
  set.seed(9878)
  train_disease_data <- train_disease_data[sample(1:nrow(train_disease_data), nsample, replace = F),]
}
# train_disease_data = screen_multivariate_outliers(train_disease_data, threshold = outlier_threshold)


## combine with synthetic
mumpsdf <- as.matrix(rbind(train_sim_data, train_disease_data))

## fit UMAP (2-dimensions)
set.seed(1128)
mumps_umap <- umap(mumpsdf, n_components = 2)
mumps_umap_df <- data.frame(mumps_umap$layout)
names(mumps_umap_df) <- paste0("X",1:ncol(mumps_umap_df))
mumps_umap_df = screen_multivariate_outliers(mumps_umap_df, threshold = outlier_threshold)
mumps_umap_df$type <- "Mumps"
mumps_umap_df[1:nrow(train_sim_data),]$type <- "Synthetic"


######################
## Load Dengue
DISEASE_KEY = 'Dengue'
### demonstrating using Dengue
output_list = readRDS(paste0(savepath,FILES[grepl(DISEASE_KEY,FILES, ignore.case = T)]))

### organize features
FEATURE_LIST = setdiff(names(output_list['features']$features),c('geography','disease','last_obs_time','h'))
features_wide = output_list['features']$features
train_disease_data = features_wide %>%
  # dplyr::filter(h==1) %>%
  dplyr::select(gr12_div_23, last_div_max, coefvar, 
                gam_with_div_without, avg_recent_div_avg_global,
                diff_zscore, entropy, relative_increases,
                prop_since_peak, seasonality)
if(nrow(train_disease_data) > nsample){
  set.seed(9878)
  train_disease_data <- train_disease_data[sample(1:nrow(train_disease_data), nsample, replace = F),]
}
# train_disease_data = screen_multivariate_outliers(train_disease_data, threshold = outlier_threshold)


## combine with synthetic
denguedf <- as.matrix(rbind(train_sim_data, train_disease_data))

## fit UMAP (2-dimensions)
set.seed(1128)
dengue_umap <- umap(denguedf, n_components = 2)
dengue_umap_df <- data.frame(dengue_umap$layout)
names(dengue_umap_df) <- paste0("X",1:ncol(dengue_umap_df))
dengue_umap_df = screen_multivariate_outliers(dengue_umap_df, threshold = outlier_threshold)
dengue_umap_df$type <- "Dengue"
dengue_umap_df[1:nrow(train_sim_data),]$type <- "Synthetic"


######################
# Let's do the plots:
## plot synthetic on top of real
patall <- grid.arrange(
  ggplot()+
    geom_point(aes(x=X1, y=X2), data=subset(global_covid_umap_df, type != "Synthetic"), color = I('#8dd3c7'))+
    geom_point(aes(x=X1, y=X2), data=subset(global_covid_umap_df, type == "Synthetic"), color = I("black"))+
    ggtitle("Global COVID-19")+
    xlim(-15, 15)+
    ylim(-15,15) +
    theme(plot.title = element_text(hjust = 0.5)),
  ggplot()+
    geom_point(aes(x=X1, y=X2), data=subset(us_covid_umap_df, type != "Synthetic"), color = I('#ffffb3'))+
    geom_point(aes(x=X1, y=X2), data=subset(us_covid_umap_df, type == "Synthetic"), color = I("black"))+
    ggtitle("US COVID-19")+
    xlim(range(dengue_umap_df$X1))+
    ylim(range(dengue_umap_df$X2)) +
    theme(plot.title = element_text(hjust = 0.5)),
  ggplot()+
    geom_point(aes(x=X1, y=X2), data=subset(ili_umap_df, type != "Synthetic"), color = I('#bebada'))+
    geom_point(aes(x=X1, y=X2), data=subset(ili_umap_df, type == "Synthetic"), color = I("black"))+
    ggtitle("ILI")+
    xlim(range(ili_umap_df$X1))+
    ylim(range(ili_umap_df$X2)) +
    theme(plot.title = element_text(hjust = 0.5)),
  ggplot()+
    geom_point(aes(x=X1, y=X2), data=subset(rubella_umap_df, type != "Synthetic"), color = I('#fb8072'))+
    geom_point(aes(x=X1, y=X2), data=subset(rubella_umap_df, type == "Synthetic"), color = I("black"))+
    ggtitle("Rubella")+
    xlim(range(rubella_umap_df$X1))+
    ylim(range(rubella_umap_df$X2)) +
    theme(plot.title = element_text(hjust = 0.5)),
  ggplot()+
    geom_point(aes(x=X1, y=X2), data=subset(diphtheria_umap_df, type != "Synthetic"), color = I('#80b1d3'))+
    geom_point(aes(x=X1, y=X2), data=subset(diphtheria_umap_df, type == "Synthetic"), color = I("black"))+
    ggtitle("Diphtheria")+
    xlim(range(diphtheria_umap_df$X1))+
    ylim(range(diphtheria_umap_df$X2)) +
    theme(plot.title = element_text(hjust = 0.5)),
  ggplot()+
    geom_point(aes(x=X1, y=X2), data=subset(measles_umap_df, type != "Synthetic"), color = I('#fdb462'))+
    geom_point(aes(x=X1, y=X2), data=subset(measles_umap_df, type == "Synthetic"), color = I("black"))+
    ggtitle("Measles")+
    xlim(range(measles_umap_df$X1))+
    ylim(range(measles_umap_df$X2)) +
    theme(plot.title = element_text(hjust = 0.5)),
  ggplot()+
    geom_point(aes(x=X1, y=X2), data=subset(mumps_umap_df, type != "Synthetic"), color = I('#b3de69'))+
    geom_point(aes(x=X1, y=X2), data=subset(mumps_umap_df, type == "Synthetic"), color = I("black"))+
    ggtitle("Mumps")+
    xlim(range(mumps_umap_df$X1))+
    ylim(range(mumps_umap_df$X2)) +
    theme(plot.title = element_text(hjust = 0.5)),
  ggplot()+
    geom_point(aes(x=X1, y=X2), data=subset(polio_umap_df, type != "Synthetic"), color = I('#fccde5'))+
    geom_point(aes(x=X1, y=X2), data=subset(polio_umap_df, type == "Synthetic"), color = I("black"))+
    ggtitle("Polio")+
    xlim(range(polio_umap_df$X1))+
    ylim(range(polio_umap_df$X2)) +
    theme(plot.title = element_text(hjust = 0.5)),
  ggplot()+
    geom_point(aes(x=X1, y=X2), data=subset(smallpox_umap_df, type != "Synthetic"), color = I('#d9d9d9'))+
    geom_point(aes(x=X1, y=X2), data=subset(smallpox_umap_df, type == "Synthetic"), color = I("black"))+
    ggtitle("Smallpox")+
    xlim(range(smallpox_umap_df$X1))+
    ylim(range(smallpox_umap_df$X2)) +
    theme(plot.title = element_text(hjust = 0.5)),
  ggplot()+
    geom_point(aes(x=X1, y=X2), data=subset(chikungunya_umap_df, type != "Synthetic"), color = I('#bc80db'))+
    geom_point(aes(x=X1, y=X2), data=subset(chikungunya_umap_df, type == "Synthetic"), color = I("black"))+
    ggtitle("Chikungunya")+
    xlim(-20,20)+
    ylim(-20,20) +
    theme(plot.title = element_text(hjust = 0.5)),ncol=1, top=textGrob("UMAP: Synthetic over Real"),
  layout_matrix = matrix(1:10, nrow = 5, ncol = 2))


tapall <- grid.arrange(
  ggplot()+
    geom_point(aes(x=X1, y=X2), data=subset(global_covid_umap_df, type == "Synthetic"), color = I("black"))+
    geom_point(aes(x=X1, y=X2), data=subset(global_covid_umap_df, type != "Synthetic"), color = I('#8dd3c7'))+
    ggtitle("Global COVID-19")+
    xlim(-15, 15)+
    ylim(-15,15) +
    theme(plot.title = element_text(hjust = 0.5)),
  ggplot()+
    geom_point(aes(x=X1, y=X2), data=subset(us_covid_umap_df, type == "Synthetic"), color = I("black"))+
    geom_point(aes(x=X1, y=X2), data=subset(us_covid_umap_df, type != "Synthetic"), color = I('#ffffb3'))+
    ggtitle("US COVID-19")+
    xlim(range(dengue_umap_df$X1))+
    ylim(range(dengue_umap_df$X2)) +
    theme(plot.title = element_text(hjust = 0.5)),
  ggplot()+
    geom_point(aes(x=X1, y=X2), data=subset(ili_umap_df, type == "Synthetic"), color = I("black"))+
    geom_point(aes(x=X1, y=X2), data=subset(ili_umap_df, type != "Synthetic"), color = I('#bebada'))+
    ggtitle("ILI")+
    xlim(range(ili_umap_df$X1))+
    ylim(range(ili_umap_df$X2)) +
    theme(plot.title = element_text(hjust = 0.5)),
  ggplot()+
    geom_point(aes(x=X1, y=X2), data=subset(rubella_umap_df, type == "Synthetic"), color = I("black"))+
    geom_point(aes(x=X1, y=X2), data=subset(rubella_umap_df, type != "Synthetic"), color = I('#fb8072'))+
    ggtitle("Rubella")+
    xlim(range(rubella_umap_df$X1))+
    ylim(range(rubella_umap_df$X2)) +
    theme(plot.title = element_text(hjust = 0.5)),
  ggplot()+
    geom_point(aes(x=X1, y=X2), data=subset(diphtheria_umap_df, type == "Synthetic"), color = I("black"))+
    geom_point(aes(x=X1, y=X2), data=subset(diphtheria_umap_df, type != "Synthetic"), color = I('#80b1d3'))+
    ggtitle("Diphtheria")+
    xlim(range(diphtheria_umap_df$X1))+
    ylim(range(diphtheria_umap_df$X2)) +
    theme(plot.title = element_text(hjust = 0.5)),
  ggplot()+
    geom_point(aes(x=X1, y=X2), data=subset(measles_umap_df, type == "Synthetic"), color = I("black"))+
    geom_point(aes(x=X1, y=X2), data=subset(measles_umap_df, type != "Synthetic"), color = I('#fdb462'))+
    ggtitle("Measles")+
    xlim(range(measles_umap_df$X1))+
    ylim(range(measles_umap_df$X2)) +
    theme(plot.title = element_text(hjust = 0.5)),
  ggplot()+
    geom_point(aes(x=X1, y=X2), data=subset(mumps_umap_df, type == "Synthetic"), color = I("black"))+
    geom_point(aes(x=X1, y=X2), data=subset(mumps_umap_df, type != "Synthetic"), color = I('#b3de69'))+
    ggtitle("Mumps")+
    xlim(range(mumps_umap_df$X1))+
    ylim(range(mumps_umap_df$X2)) +
    theme(plot.title = element_text(hjust = 0.5)),
  ggplot()+
    geom_point(aes(x=X1, y=X2), data=subset(polio_umap_df, type == "Synthetic"), color = I("black"))+
    geom_point(aes(x=X1, y=X2), data=subset(polio_umap_df, type != "Synthetic"), color = I('#fccde5'))+
    ggtitle("Polio")+
    xlim(range(polio_umap_df$X1))+
    ylim(range(polio_umap_df$X2)) +
    theme(plot.title = element_text(hjust = 0.5)),
  ggplot()+
    geom_point(aes(x=X1, y=X2), data=subset(smallpox_umap_df, type == "Synthetic"), color = I("black"))+
    geom_point(aes(x=X1, y=X2), data=subset(smallpox_umap_df, type != "Synthetic"), color = I('#d9d9d9'))+
    ggtitle("Smallpox")+
    xlim(range(smallpox_umap_df$X1))+
    ylim(range(smallpox_umap_df$X2)) +
    theme(plot.title = element_text(hjust = 0.5)),
  ggplot()+
    geom_point(aes(x=X1, y=X2), data=subset(chikungunya_umap_df, type == "Synthetic"), color = I("black"))+
    geom_point(aes(x=X1, y=X2), data=subset(chikungunya_umap_df, type != "Synthetic"), color = I('#bc80db'))+
    ggtitle("Chikungunya")+
    xlim(-20,20)+
    ylim(-20,20) +
    theme(plot.title = element_text(hjust = 0.5)),ncol=1, top=textGrob("UMAP: Synthetic over Real"),
  layout_matrix = matrix(1:10, nrow = 5, ncol = 2))




