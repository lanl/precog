## Lauren Beesley and Casey Gibson
## 1-31-24
## Download and organize dengue data for training in epiFFORMA

###############################################################

###############################################################
## load libraries
library(MMWRweek)
library(lubridate)
library(ggplot2)
library(plyr)

theme_set(theme_bw())

###############################################################
## set path

library(this.path)
setwd(this.path::here()) 
datapath <- "./raw_data/"
savepath <- "./output/"
plotpath <- "./figs/"


###############################################################
#### Get Dengue Data

### Peru
iquitos <- read.csv(paste0(datapath,"Iquitos_Training_Data.csv"),header = T)
iquitos$date <- as.Date(iquitos$week_start_date)
iquitos$adm_0_name = 'IQUITOS_PERU'
iquitos$dengue_total = iquitos$total_cases
iquitos$T_res = 'Week'
iquitos = iquitos[,c('adm_0_name','date', 'dengue_total', 'T_res')]

### San Juan
sanjuan <- read.csv(paste0(datapath,"San_Juan_Training_Data.csv"),header = T)
sanjuan$date <- as.Date(sanjuan$week_start_date)
sanjuan$adm_0_name = 'SANJUAN_PUERTORICO'
sanjuan$dengue_total = sanjuan$total_cases
sanjuan$T_res = 'Week'
sanjuan = sanjuan[,c('adm_0_name','date', 'dengue_total', 'T_res')]

### OpenDengue Data
#https://opendengue.org/data.html
opendengue <- read.csv(paste0(datapath,"National-level data_AFGHANISTAN_19240120_20230930.csv"),header = T)
opendengue = opendengue[opendengue$S_res == 'Admin0',] #only value in dataset
opendengue$date = as.Date(opendengue$calendar_start_date)
opendengue = opendengue[,c('adm_0_name','date', 'dengue_total', 'T_res')]
opendengue = opendengue[opendengue$T_res == 'Week',]

## combine data
dat        <- rbind(iquitos, sanjuan, opendengue)
names(dat) <- tolower(names(dat))


### find missing weeks
create_weekly_dataframe <- function(x, y) {
  seq_dates <- seq(x, y + 7, by = "week")
  date <- seq_dates - as.integer(format(x, '%u'))
  subset(data.frame(date = date), date <= y)
}
df_long = NULL
for(l in 1:length(unique(dat$adm_0_name))){
  SUBSET = dat[dat$adm_0_name == unique(dat$adm_0_name)[l],]
  temp = create_weekly_dataframe(x=min(SUBSET$date), y=max(SUBSET$date))
  temp = merge(data.frame(temp, adm_0_name = unique(dat$adm_0_name)[l]), SUBSET, by = c('date','adm_0_name'), all.x = T, all.y = T)
  temp = temp[!is.na(temp$t_res),]
  df_long = rbind(df_long, temp)
}

# df_long$week <- lubridate::isoweek(df_long$date)
df_long$year <- lubridate::isoyear(df_long$date)



A = df_long[df_long$adm_0_name == 'IQUITOS_PERU',]

### only use data since 1990
df_long = df_long[df_long$year >= 1990,]

dat = df_long

## replace spaces with underscores
dat$region <- gsub("\\ ","_",dat$adm_0_name)

## subsetting to active locations
library(dplyr)
dat = dat %>% dplyr::group_by(region) %>% dplyr::mutate(total_cases = sum(dengue_total,na.rm=T), 
                                                        total_weeks = length(dengue_total[dengue_total>0]))
dat = dat[dat$total_weeks >= 52,]

LOCS_TO_INCLUDE = c('VIET_NAM','UNITED_STATES_OF_AMERICA', 'SRI_LANKA', 
                    'SINGAPORE',
                    'SANJUAN_PUERTORICO', 'PARAGUAY', 'PERU',
                    'NICARAGUA','MEXICO','MALAYSIA', 'LAO_PEOPLE\'S_DEMOCRATIC_REPUBLIC',
                    'IQUITOS_PERU', 'HONDURAS', 'EL_SALVADOR','ECUADOR', 'DOMINICAN_REPUBLIC',
                    'COSTA_RICA','COLOMBIA','CAMBODIA','BRAZIL','BOLIVIA', 'BELIZE',
                    'BARBADOS', 'ARGENTINA')
dat = dat[dat$region %in% LOCS_TO_INCLUDE,]


## final formatting
## final formatting
dat       <- dat[order(dat$date),]
dat = data.frame(dat)


unqstate <- sort(unique(dat$region))

### Spot Cleaning

dat$date[as.character(dat$date) == '2016-02-29'] = dat$date[as.character(dat$date) == '2016-02-29'] - 1 #change date to make regular weekly spacing


df_long = NULL
for(i in 1:length(unqstate)){
  tempdf <- subset(dat, region == unqstate[i])
  tempdf <- tempdf[order(tempdf$date),]
  tempdf$date_char = as.character(tempdf$date)
  #tempdf$date_next = tempdf$date + 7
  DIFF = difftime(tempdf$date[-1], tempdf$date[-length(tempdf$date)])
  if(length(DIFF)>0){
    print(unqstate[i])
    print(tempdf[which(DIFF/7 != round(DIFF/7)),])
  }
  temp = create_weekly_dataframe(x=min(tempdf$date,na.rm=T), y=max(tempdf$date, na.rm=T))
  
  # print(tempdf[which(DIFF/7 != round(DIFF/7))-1,])
  # print(tempdf[which(DIFF/7 != round(DIFF/7))+1,])
  # 
  if(length(DIFF) == 0){
    tempdf$date_char = as.character(tempdf$date)
    tempdf = merge(tempdf, data.frame(date_char = as.character(temp$date+1), week_index = c(1:length(temp$date))),
                by = 'date_char', all.x = T, all.y = F)
  }else{ #only Iquitos, Peru and San Juan. Checked that weekly works. They restart reporting on 01/01 each year 
    tempdf$week_index = 1:length(tempdf[,1])
  }
  df_long = rbind(df_long, tempdf)
}


dat = df_long


# 
# 
# A = dat[dat$adm_0_name == 'BRAZIL',]
# A = A[order(A$date),]
# data.frame(A[100:108,])
#################################
### Visualize the Time Series ###
#################################
dat = dat %>% dplyr::group_by(region) %>% dplyr::mutate(cases_scaled = dengue_total/max(dengue_total,na.rm=T))
  p1 = ggplot(dat)+
    geom_tile(aes(x=week_index,y=region,fill=cases_scaled,color=cases_scaled))+
    theme_classic()+
    scale_color_viridis_c('Scaled Cases', na.value = 'lightgray')+
    scale_fill_viridis_c('Scaled Cases', na.value = 'lightgray')+
    xlab('Week')+
    ylab('Location')+
    labs(title = 'Scaled Cases')+
    scale_x_continuous(expand=c(0,0))+
    theme(legend.position = 'top')
  p2 = ggplot(dat)+
    geom_tile(aes(x=week_index,y=region,fill=cut(dengue_total, breaks = c(0,5,10,100,1000,10000,10^10)),color=cut(dengue_total, breaks = c(0,5,10,100,1000,10000,10^10))))+
    theme_classic()+
    scale_color_viridis_d('Cases', na.value = 'lightgray')+
    scale_fill_viridis_d('Cases', na.value = 'lightgray')+
    xlab('Week')+
    ylab('Location')+
    labs(title = 'Cases')+
    scale_x_continuous(expand=c(0,0))+
    theme(legend.position = 'top')
  p0= cowplot::plot_grid(p1,p2,ncol=2)
  ggsave(p0,filename = paste0(plotpath,'/','dengue_cases.png'), width = 10, height = 8)
  

  

####################################################  

unqstate <- sort(unique(dat$region))

## DENGUE rollercoaster
dengue_rollercoaster_training_data <- list()
cnt <- 0
for(i in 1:length(unqstate)){
  print(i)
  tempdf <- subset(dat, region == unqstate[i])
  tempdf <- tempdf[order(tempdf$date),]
  if(nrow(tempdf)>= 10  & var(tempdf$dengue_total, na.rm=T) > 0){
    cnt <- cnt + 1
    ts = rep(0,max(tempdf$week_index) - min(tempdf$week_index) + 1)
    ts[tempdf$week_index - min(tempdf$week_index)+1] = tempdf$dengue_total
    
    
    
    ts_dates = rep(NA,max(tempdf$week_index) - min(tempdf$week_index) + 1)
    ts_dates[tempdf$week_index - min(tempdf$week_index)+1] = as.character(tempdf$date)
    
    
    
    templist <- list(ts = ts,
                     ts_dates = ts_dates,
                     ts_exogenous = NULL,
                     ts_real_data = T,
                     ts_isolated_strain = F,
                     ts_multiwave = T,
                     ts_disease = "dengue",
                     ts_measurement_type = "incidence",
                     ts_geography = tempdf$region[1],
                     ts_first_time = min(tempdf$date),
                     ts_last_time = max(tempdf$date),
                     ts_time_cadence = "weekly",
                     ts_scale = "counts",
                     ts_exogenous_scale = NA, 
                     ts_seasonal = F)
    
    dengue_rollercoaster_training_data[[cnt]] <- templist
  }
}


## concatenate all the lists
dengue_training_data <- c(dengue_rollercoaster_training_data)

## save the us ili training data
saveRDS(dengue_training_data, paste0(savepath,"dengue.RDS"))







