# Script for the Figure 3 in the epiFFORMA paper.
## Authors: Lauren Beesley and Alexander C. Murph
## Date: December 2025
library(igraph)
library(mclust)
library(FNN)
library(dbscan)
library(Rtsne)
library(ggrepel)
library(patchwork)
library(ggplot2)
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
FULL_WEIGHTS_LIST = NULL
for(i in 1:length(FILES)){
  
  output_list = readRDS(paste0(savepath, FILES[i]))
  FULL_WEIGHTS_LIST = rbind(FULL_WEIGHTS_LIST, output_list$weights)
  
  rm(output_list)
  gc()
  print(paste0('Finished: ', i, ' of ', length(FILES)))
}

# ---- 1) Synthetic 9D data ----
orig = FULL_WEIGHTS_LIST
FULL_WEIGHTS_LIST = rbind(data.frame(geography = "equal", disease = "equal", last_obs_time  = "equal", h = 1, 
                                     gam_avg = 1/9,    moa_avg = 1/9, moa_deriv_avg  = 1/9,   rw_avg = 1/9,  
                                     theta_avg = 1/9, meanfcst_avg  = 1/9,mirror_avg = 1/9, gam2mirror_avg  = 1/9, arima_avg = 1/9),
                          FULL_WEIGHTS_LIST) %>%
  subset(disease == "covid")
X <- FULL_WEIGHTS_LIST[c("gam_avg", "moa_avg", "moa_deriv_avg", "rw_avg", "theta_avg", "meanfcst_avg",
                       "mirror_avg", "gam2mirror_avg", "arima_avg") ] 

k2 = kmeans(X, centers = 5, nstart = 50)

table(k2$cluster) / length(k2$cluster)
head(k2$cluster)


clusters <- k2$cluster

dist_from_eq_wts = sqrt(rowSums((X - rep(1/9, times = 9))^2))

FULL_WEIGHTS_LIST$dist_from_eq_wts = dist_from_eq_wts

ggplot(FULL_WEIGHTS_LIST, aes(x = as.factor(clusters), y = dist_from_eq_wts, group = as.factor(clusters), fill =as.factor(clusters))) + geom_boxplot()



