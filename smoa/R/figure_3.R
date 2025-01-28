# Visualizations for the coverage analysis in the sMOA paper.
## Author: Alexander C. Murph
## Date: October 2024
library(lhs)
library(mgcv)
library(ggplot2)
library(plyr)
library(data.table)
library(LearnBayes)
library(LaplacesDemon)
#library(covidHubUtils)
library(dplyr)
library(spatstat)
library(BASS)
#library(covidcast)
library(GPfit)
library(nnet)
library(dplyr)
library(KernelKnn)
library(tidybayes)
library(gridExtra)
library(ggExtra)
library(dplyr)
setwd("~/GitLab/smoa")
source("R/smoa_helpers.R")
source('R/vecchia_scaled.R')

data_path = 'data/coverage_data'

coverage_data <- NULL
for(file_name in list.files(data_path)){
  temp_coverage_data = read.csv(paste(data_path, '/', file_name, sep = ''))
  temp_f = function(x){as.numeric(unlist(strsplit(x, split="%"))[1])}
  temp_coverage_data$Nominal_Coverage = (100 - 2*sapply(temp_coverage_data$X,temp_f))/100
  
  
  coverage_data <- rbind(coverage_data, temp_coverage_data)
}

graph_data = coverage_data
graph_data %>% dplyr::group_by(Nominal_Coverage) %>% dplyr::summarise(Empirical_Coverage = mean(Coverage),
                                                        Nominal_Coverage = Nominal_Coverage[1]) %>%
  ggplot(aes(x = Nominal_Coverage, y = Empirical_Coverage)) + geom_line() + geom_abline(linetype='dotted') +
  theme_bw() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_text(hjust=1, vjust=1)
  )+ 
  theme(axis.text=element_text(size=15), axis.title = element_text(size = 18),
        title = element_text(size=15),strip.text.x = element_text(size = 15),strip.text.y = element_text(size = 15)) +
  theme(panel.spacing = unit(0, "cm"),
        plot.caption = element_blank()) + theme(legend.text=element_text(size=15)) + 
  ggtitle(paste0('Coverage Plot for sMOA')) +
  xlim(0,1) + ylim(0,1) + xlab("Nominal Coverage") + ylab("Empirical Coverage")
  
