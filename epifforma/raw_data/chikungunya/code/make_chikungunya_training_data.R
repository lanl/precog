## Lauren Beesley 
## 6-12-24
## Download and organize chikungunya data for training in epiFFORMA

### this time series for locations in brazil was downloaded from
#https://github.com/wmarciel/Chikungunya-in-Brazil-2013-2022/blob/main/Fig1/Fig1ab_CHIKV_Brazil_line_chart.csv
#made available as part of this paper
#https://www.ncbi.nlm.nih.gov/pmc/articles/PMC10281060/

### Another data sources if needed:
# https://bmcmedicine.biomedcentral.com/articles/10.1186/s12916-018-1127-2
# https://github.com/mooresea/chikv-colombia-emod/tree/master/inputs 
# https://www.sciencedirect.com/science/article/pii/S1755436517300014
# https://ars.els-cdn.com/content/image/1-s2.0-S1755436517300014-mmc2.csv (includes zika!!!!)
# 

###############################################################

###########################################s###################
## load libraries
library(MMWRweek)
library(lubridate)
library(ggplot2)
library(plyr)

theme_set(theme_bw())

###############################################################
## set path
setwd("~/GitLab/epifforma/")
datapath <- "./raw_data/chikungunya/raw_data/"
savepath <- "./raw_data/chikungunya/output/"
plotpath <- "./raw_data/chikungunya/figs/"

###############################################################
## set proxies
Sys.setenv('https_proxy'='http://proxyout.lanl.gov:8080')## not sure why needed, but see here: https://github.com/curl/curl/issues/1015

###############################################################
#### Get Dengue Data

### Peru
chikv_brazil <- read.csv(paste0(datapath,"Fig1ab_CHIKV_Brazil_line_chart.csv"),header = T)
chikv_brazil$time =52*chikv_brazil$year + as.numeric(substr(chikv_brazil$epi_week, 6, nchar(chikv_brazil$epi_week)))
chikv_brazil$epi_date   <- MMWRweek2Date(MMWRyear = chikv_brazil$year, MMWRweek = as.numeric(substr(chikv_brazil$epi_week, 6, nchar(chikv_brazil$epi_week))))
chikv_brazil$week <- as.numeric(substr(chikv_brazil$epi_week, 6, nchar(chikv_brazil$epi_week)))
chikv_brazil$region = substr(chikv_brazil$uf,3,nchar(chikv_brazil$uf))

# ggplot(chikv_brazil)+
#   geom_line(aes(x=time, y=CHIKV_cases, group = factor(uf), color = factor(uf)))+
#   guides(color = 'none')
# 

# 
# ## subsetting to active locations
library(dplyr)
chikv_brazil = chikv_brazil %>% dplyr::group_by(region) %>% dplyr::mutate(total_cases = sum(CHIKV_cases,na.rm=T),
                                                        total_weeks = length(CHIKV_cases[CHIKV_cases>0]))

chikv_brazil[which.min(chikv_brazil$epi_date),]
chikv_brazil[which.max(chikv_brazil$epi_date),]




## final formatting
create_weekly_dataframe <- function(x, y) {
  seq_dates <- seq(x, y + 7, by = "week")
  date <- seq_dates - as.integer(format(x, '%u'))
  subset(data.frame(date = date), date <= y)
}
chikv_brazil       <- chikv_brazil[order(chikv_brazil$epi_date),]
chikv_brazil = data.frame(chikv_brazil)
temp = create_weekly_dataframe(x=min(chikv_brazil$epi_date,na.rm=T), y=max(chikv_brazil$epi_date, na.rm=T))


chikv_brazil$date_char = as.character(chikv_brazil$epi_date)
chikv_brazil = merge(chikv_brazil, data.frame(date_char = as.character(temp$date), week_index = c(1:length(temp$date))),
                     by = 'date_char', all.x = T, all.y = F)



# 
# 
# A = dat[dat$adm_0_name == 'BRAZIL',]
# A = A[order(A$date),]
# data.frame(A[100:108,])
#################################
### Visualize the Time Series ###
#################################
library(dplyr)
chikv_brazil = chikv_brazil %>% dplyr::group_by(region) %>% dplyr::mutate(cases_scaled = CHIKV_cases/max(CHIKV_cases,na.rm=T))
  p1 = ggplot(chikv_brazil)+
    geom_tile(aes(x=week_index,y=region,fill=cases_scaled,color=cases_scaled))+
    theme_classic()+
    scale_color_viridis_c('Scaled Cases', na.value = 'lightgray')+
    scale_fill_viridis_c('Scaled Cases', na.value = 'lightgray')+
    xlab('Week')+
    ylab('Location')+
    labs(title = 'Scaled Cases')+
    scale_x_continuous(expand=c(0,0))+
    theme(legend.position = 'top')
  p2 = ggplot(chikv_brazil)+
    geom_tile(aes(x=week_index,y=region,fill=cut(CHIKV_cases, breaks = c(0,5,10,100,1000,10000,10^10)),color=cut(CHIKV_cases, breaks = c(0,5,10,100,1000,10000,10^10))))+
    theme_classic()+
    scale_color_viridis_d('Cases', na.value = 'lightgray')+
    scale_fill_viridis_d('Cases', na.value = 'lightgray')+
    xlab('Week')+
    ylab('Location')+
    labs(title = 'Cases')+
    scale_x_continuous(expand=c(0,0))+
    theme(legend.position = 'top')
  p0= cowplot::plot_grid(p1,p2,ncol=2)
  ggsave(p0,filename = paste0(plotpath,'/','chikv_cases.png'), width = 10, height = 8)
  

  

####################################################  

unqstate <- sort(unique(chikv_brazil$region))

## chikv rollercoaster
chikv_rollercoaster_training_data <- list()
cnt <- 0
for(i in 1:length(unqstate)){
  print(i)
  tempdf <- subset(chikv_brazil, region == unqstate[i])
  tempdf <- tempdf[order(tempdf$epi_date),]
  if(nrow(tempdf)>= 10  & var(tempdf$CHIKV_cases, na.rm=T) > 0){
    cnt <- cnt + 1
    ts = rep(0,max(tempdf$week_index) - min(tempdf$week_index) + 1)
    ts[tempdf$week_index - min(tempdf$week_index)+1] = tempdf$CHIKV_cases
    
    ts_dates = rep(0,max(tempdf$week_index) - min(tempdf$week_index) + 1)
    ts_dates[tempdf$week_index - min(tempdf$week_index)+1] = tempdf$epi_date
    templist <- list(ts = ts,
                     ts_dates = ts_dates,
                     ts_exogenous = NULL,
                     ts_real_data = T,
                     ts_isolated_strain = F,
                     ts_multiwave = T,
                     ts_disease = "chikv",
                     ts_measurement_type = "incidence",
                     ts_geography = tempdf$region[1],
                     ts_first_time = min(tempdf$epi_date),
                     ts_last_time = max(tempdf$epi_date),
                     ts_time_cadence = "weekly",
                     ts_scale = "counts",
                     ts_exogenous_scale = NA, 
                     ts_seasonal = F) #could list as seasonal, but it is not reliably seasonal for some locations
    
    chikv_rollercoaster_training_data[[cnt]] <- templist
  }
}


## concatenate all the lists
chikv_training_data <- c(chikv_rollercoaster_training_data)

## save the us ili training data
saveRDS(chikv_training_data, paste0(savepath,"chikv.RDS"))







