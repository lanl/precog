#########################
#########################
### Clean Dengue Data ###
#########################
#########################


## load libraries
library(MMWRweek)
library(lubridate)
library(ggplot2)
library(plyr)
library(lubridate)
theme_set(theme_bw())

## set path
datapath <- "./raw_data/"
savepath <- "../_Organized_Lists/"

### find missing weeks per location and augment with zero rows
create_weekly_dataframe <- function(x, y) {
  seq_dates <- seq(x, y + 7, by = "week")
  date <- seq_dates - as.integer(format(x, '%u'))
  subset(data.frame(date = date), date <= y)
}

####################
### Read in Data ###
####################

### OpenDengue Data
#https://opendengue.org/data.html
opendengue <- read.csv(paste0(datapath,"National-level data_AFGHANISTAN_19240120_20230930.csv"),header = T)
opendengue = opendengue[opendengue$S_res == 'Admin0',] #only value in dataset
opendengue$date = as.Date(opendengue$calendar_start_date)
opendengue = opendengue[,c('adm_0_name','date', 'dengue_total', 'T_res')]
opendengue = opendengue[opendengue$T_res == 'Week',]
opendengue$denv1_cases = NA
opendengue$denv2_cases = NA
opendengue$denv3_cases = NA
opendengue$denv4_cases = NA
opendengue$date_actual = opendengue$date
opendengue_long = NULL
for(l in 1:length(unique(opendengue$adm_0_name))){
  SUBSET = opendengue[opendengue$adm_0_name == unique(opendengue$adm_0_name)[l],]
  temp = create_weekly_dataframe(x=min(SUBSET$date), y=max(SUBSET$date))
  temp = merge(data.frame(temp, adm_0_name = unique(opendengue$adm_0_name)[l]), SUBSET, by = c('date','adm_0_name'), all.x = T, all.y = T)
  temp$dengue_total[is.na(temp$dengue_total)]=0
  opendengue = rbind(opendengue, temp)
}
### make adjustment for one oddball leap year reporting
opendengue = opendengue[as.character(opendengue$date) != '2016-02-28',]
opendengue$date[as.character(opendengue$date) == '2016-02-29'] = opendengue$date[as.character(opendengue$date) == '2016-02-29'] - 1 #change date to make regular weekly spacing
opendengue = opendengue[year(opendengue$date) >= 1990]

### Peru
iquitos <- read.csv(paste0(datapath,"Iquitos_Training_Data.csv"),header = T)
iquitos$date <- as.Date(iquitos$week_start_date)
iquitos$adm_0_name = 'IQUITOS_PERU'
iquitos$dengue_total = iquitos$total_cases
iquitos$T_res = 'Week'
iquitos = iquitos[,c('adm_0_name','date', 'dengue_total', 'T_res', 'denv1_cases', 'denv2_cases','denv3_cases','denv4_cases' )]
iquitos$date_actual = iquitos$date
DIFF = difftime(iquitos$date[-1], iquitos$date[-length(iquitos$date)], units = 'days')
#reporting restarted on the first of each year and avoided xmas. Inconvenient for modeling, pretend the data had been reported on regular weekly basis
#date and date_actual fields will be different, especially for later times
temp = create_weekly_dataframe(x=min(iquitos$date), y=max(iquitos$date))
temp = temp[-1,]
iquitos$date = temp[1:length(iquitos$date)]

### San Juan
sanjuan <- read.csv(paste0(datapath,"San_Juan_Training_Data.csv"),header = T)
sanjuan$date <- as.Date(sanjuan$week_start_date)
sanjuan$adm_0_name = 'SANJUAN_PUERTORICO'
sanjuan$dengue_total = sanjuan$total_cases
sanjuan$T_res = 'Week'
sanjuan = sanjuan[,c('adm_0_name','date', 'dengue_total', 'T_res', 'denv1_cases', 'denv2_cases','denv3_cases','denv4_cases')]
sanjuan$date_actual = sanjuan$date
sanjuan = sanjuan[order(sanjuan$date),]
DIFF = difftime(sanjuan$date[-1], sanjuan$date[-length(sanjuan$date)], units = 'days')
#reporting restarted on the first of each year and avoided xmas. Inconvenient for modeling, pretend the data had been reported on regular weekly basis
#date and date_actual fields will be different, especially for later times
temp = create_weekly_dataframe(x=min(sanjuan$date), y=max(sanjuan$date))
temp = temp[-1,]
sanjuan$date = temp[1:length(sanjuan$date)]





## combine data
dat        <- rbind(iquitos, sanjuan, opendengue)
names(dat) <- tolower(names(dat))

#####################
### Clean up Data ###
#####################


df_long$week <- lubridate::isoweek(df_long$date)
df_long$year <- lubridate::isoyear(df_long$date)


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

dat       <- dat[order(dat$date),]
dat = data.frame(dat)


### add week indices per time series
unqstate <- sort(unique(dat$region))
df_long = NULL
for(i in 1:length(unqstate)){
  tempdf <- subset(dat, region == unqstate[i])
  tempdf <- tempdf[order(tempdf$date),]
  tempdf$date_char = as.character(tempdf$date)
  tempdf$week_index = 1:length(tempdf[,1])
  DIFF = difftime(tempdf$date[-1], tempdf$date[-length(tempdf$date)], units = 'days')
  df_long = rbind(df_long, tempdf)
}



dat = df_long

######################
### Organize Lists ###
######################


unqstate <- sort(unique(dat$region))

## Dengue incidence
dengue_training_data <- list()
cnt <- 0
for(i in 1:length(unqstate)){
  print(i)
  tempdf <- subset(dat, region == unqstate[i])
  tempdf <- tempdf[order(tempdf$date),]
  if(nrow(tempdf)>= 10  & var(tempdf$dengue_total, na.rm=T) > 0){
    cnt <- cnt + 1
    ts = rep(0,max(tempdf$week_index) - min(tempdf$week_index) + 1)
    ts[tempdf$week_index - min(tempdf$week_index)+1] = tempdf$dengue_total
    ts[is.na(ts)]=0
    
    ts_dates = rep(NA,max(tempdf$week_index) - min(tempdf$week_index) + 1)
    ts_dates[tempdf$week_index - min(tempdf$week_index)+1] = as.character(tempdf$date)
    
    ts_dates_actual = tempdf$date_actual
    
    templist <- list(ts = ts,
                     ts_dates = ts_dates,
                     ts_dates_actual = ts_dates_actual, 
                     ts_disease = "dengue",
                     ts_measurement_type = "incidence",
                     ts_geography = tempdf$region[1],
                     ts_first_time = min(tempdf$date),
                     ts_last_time = max(tempdf$date),
                     ts_time_cadence = "weekly",
                     ts_scale = "counts")
    
    dengue_training_data[[cnt]] <- templist
  }
}

## save the us ili training data
saveRDS(dengue_training_data, paste0(savepath,"Dengue_opendengue.RDS"))












unqstate <- sort(unique(dat$region[!is.na(dat$denv1_cases)]))

## DENGUE rollercoaster
dengue_rollercoaster_training_data <- list()
cnt <- 0
for(i in 1:length(unqstate)){
  print(i)
  tempdf <- subset(dat, region == unqstate[i])
  tempdf <- tempdf[order(tempdf$date),]
  if(nrow(tempdf)>= 10  & var(tempdf$denv1_cases, na.rm=T) > 0){
    cnt <- cnt + 1
    ts = rep(0,max(tempdf$week_index) - min(tempdf$week_index) + 1)
    ts[tempdf$week_index - min(tempdf$week_index)+1] = tempdf$denv1_cases
    ts[is.na(ts)]=0
    
    
    ts_dates = rep(NA,max(tempdf$week_index) - min(tempdf$week_index) + 1)
    ts_dates[tempdf$week_index - min(tempdf$week_index)+1] = as.character(tempdf$date)
    ts_dates_actual = tempdf$date_actual
    
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
saveRDS(dengue_rollercoaster_training_data, paste0(savepath,"DengueSero1_opendengue.RDS"))






unqstate <- sort(unique(dat$region[!is.na(dat$denv2_cases)]))

## DENGUE rollercoaster
dengue_rollercoaster_training_data <- list()
cnt <- 0
for(i in 1:length(unqstate)){
  print(i)
  tempdf <- subset(dat, region == unqstate[i])
  tempdf <- tempdf[order(tempdf$date),]
  if(nrow(tempdf)>= 10  & var(tempdf$denv2_cases, na.rm=T) > 0){
    cnt <- cnt + 1
    ts = rep(0,max(tempdf$week_index) - min(tempdf$week_index) + 1)
    ts[tempdf$week_index - min(tempdf$week_index)+1] = tempdf$denv2_cases
    ts[is.na(ts)]=0
    
    
    ts_dates = rep(NA,max(tempdf$week_index) - min(tempdf$week_index) + 1)
    ts_dates[tempdf$week_index - min(tempdf$week_index)+1] = as.character(tempdf$date)
    ts_dates_actual = tempdf$date_actual
    
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
saveRDS(dengue_rollercoaster_training_data, paste0(savepath,"DengueSero2_opendengue.RDS"))




unqstate <- sort(unique(dat$region[!is.na(dat$denv3_cases)]))

## DENGUE rollercoaster
dengue_rollercoaster_training_data <- list()
cnt <- 0
for(i in 1:length(unqstate)){
  print(i)
  tempdf <- subset(dat, region == unqstate[i])
  tempdf <- tempdf[order(tempdf$date),]
  if(nrow(tempdf)>= 10  & var(tempdf$denv3_cases, na.rm=T) > 0){
    cnt <- cnt + 1
    ts = rep(0,max(tempdf$week_index) - min(tempdf$week_index) + 1)
    ts[tempdf$week_index - min(tempdf$week_index)+1] = tempdf$denv3_cases
    ts[is.na(ts)]=0
    
    
    ts_dates = rep(NA,max(tempdf$week_index) - min(tempdf$week_index) + 1)
    ts_dates[tempdf$week_index - min(tempdf$week_index)+1] = as.character(tempdf$date)
    ts_dates_actual = tempdf$date_actual
    
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
saveRDS(dengue_rollercoaster_training_data, paste0(savepath,"DengueSero3_opendengue.RDS"))






unqstate <- sort(unique(dat$region[!is.na(dat$denv4_cases)]))

## DENGUE rollercoaster
dengue_rollercoaster_training_data <- list()
cnt <- 0
for(i in 1:length(unqstate)){
  print(i)
  tempdf <- subset(dat, region == unqstate[i])
  tempdf <- tempdf[order(tempdf$date),]
  if(nrow(tempdf)>= 10  & var(tempdf$denv4_cases, na.rm=T) > 0){
    cnt <- cnt + 1
    ts = rep(0,max(tempdf$week_index) - min(tempdf$week_index) + 1)
    ts[tempdf$week_index - min(tempdf$week_index)+1] = tempdf$denv4_cases
    ts[is.na(ts)]=0
    
    
    ts_dates = rep(NA,max(tempdf$week_index) - min(tempdf$week_index) + 1)
    ts_dates[tempdf$week_index - min(tempdf$week_index)+1] = as.character(tempdf$date)
    ts_dates_actual = tempdf$date_actual
    
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
saveRDS(dengue_rollercoaster_training_data, paste0(savepath,"DengueSero4_opendengue.RDS"))









