#############################
#############################
### Clean WHO FluNet Data ###
#############################
#############################

library(MMWRweek)
library(lubridate)
library(ggplot2)
library(plyr)
library(lubridate)
theme_set(theme_bw())
library(dplyr)

datapath <- "./raw_data/" 
savepath <- "../_Organized_Lists/"

create_weekly_dataframe <- function(x, y) {
  seq_dates <- seq(x, y + 7, by = "week")
  subset(data.frame(date = seq_dates), date <= y)
}

####################
### Read in Data ###
####################
dat = data.frame(data.table::fread(paste0(datapath,'VIW_FNT.csv')))
dat$MMWR_WEEKSTARTDATE = as.Date(dat$MMWR_WEEKSTARTDATE, format = '%Y-%m-%d')
dat$epi_date = dat$MMWR_WEEKSTARTDATE
dat       <- dat[order(dat$epi_date),]
dat = data.frame(dat)







##########################
### Harmonize Diseases ###
##########################

### H1N1
dat$Cases = dat$AH1N12009
dat = dat %>% dplyr::group_by(COUNTRY_AREA_TERRITORY) %>% dplyr::mutate(num_nonzero = sum(Cases[!is.na(Cases) & Cases>0]))
dat_sub = dat[dat$num_nonzero > 50,]
LOC_TO_EXCLUDE = c('Albania', 'Algeria' ,'Armenia', 'Kuwait', 'Krygyzstan','Lithuania','Mauritius',
                   'Malta','North Macedonia')
SUB = dat_sub[!(dat_sub$COUNTRY_AREA_TERRITORY %in% LOC_TO_EXCLUDE),]
data_list = list()
num = 0
LOCATIONS = unique(SUB$COUNTRY_AREA_TERRITORY)
for(l in 1:length(LOCATIONS)){
  SUBSUB = SUB[SUB$COUNTRY_AREA_TERRITORY == LOCATIONS[l],]
  SUBSUB = SUBSUB[order(SUBSUB$epi_date),]
  SUBSUB$Cases[is.na(SUBSUB$Cases)] = 0
  temp = create_weekly_dataframe(x=min(SUBSUB$epi_date,na.rm=T), y=max(SUBSUB$epi_date, na.rm=T))
  SUBSUB$date_char = as.character(SUBSUB$epi_date)
  SUBSUB = merge(SUBSUB, data.frame(date_char = as.character(temp$date), week_index = c(1:length(temp$date))),
              by = 'date_char', all.x = T, all.y = T)
  ts = rep(0,max(SUBSUB$week_index) - min(SUBSUB$week_index) + 1)
  ts[SUBSUB$week_index - min(SUBSUB$week_index)+1] = SUBSUB$Cases
  ts[is.na(ts)]=0
  ts_dates = rep(NA,max(SUBSUB$week_index) - min(SUBSUB$week_index) + 1)
  ts_dates[SUBSUB$week_index - min(SUBSUB$week_index)+1] = as.character(SUBSUB$date_char)
  ts_dates_orig = rep(NA,max(SUBSUB$week_index) - min(SUBSUB$week_index) + 1)
  ts_dates_orig[SUBSUB$week_index - min(SUBSUB$week_index)+1] = as.character(SUBSUB$epi_date)
  
  inds_to_include = min(which(ts>0)):max(which(ts>0))
  if(nrow(SUBSUB)>= 10 & !is.na(var(SUBSUB$Cases, na.rm=T)) & var(SUBSUB$Cases, na.rm=T) > 0){
    num <- num + 1
    templist <- list(ts = pmax(0,ts)[inds_to_include],
                     ts_dates = ts_dates[inds_to_include],
                     ts_dates_actual = ts_dates_orig[inds_to_include],
                     ts_disease = 'H1N1',
                     ts_measurement_type = "incidence",
                     ts_geography = trimws(LOCATIONS[l]),
                     ts_first_time = min(ts_dates[inds_to_include], na.rm=T),
                     ts_last_time = max(ts_dates[inds_to_include], na.rm=T),
                     ts_time_cadence = 'weekly',
                     ts_scale = "counts")
    data_list[[num]] <- templist
  }
}
saveRDS(data_list, paste0(savepath,"H1N1_whoflunet.RDS"))




### Subtype A
dat$Cases = dat$INF_A
dat = dat %>% dplyr::group_by(COUNTRY_AREA_TERRITORY) %>% dplyr::mutate(num_nonzero = sum(Cases[!is.na(Cases) & Cases>0]))
dat_sub = dat[dat$num_nonzero > 50,]
LOC_TO_EXCLUDE = c('Kuwait', 'Malta', 'North Macedonia', 'Mauritius', 'Saudi Arabia', 'Zimbabwe')
SUB = dat_sub[!(dat_sub$COUNTRY_AREA_TERRITORY %in% LOC_TO_EXCLUDE),]
data_list = list()
num = 0
LOCATIONS = unique(SUB$COUNTRY_AREA_TERRITORY)
for(l in 1:length(LOCATIONS)){
  SUBSUB = SUB[SUB$COUNTRY_AREA_TERRITORY == LOCATIONS[l],]
  SUBSUB = SUBSUB[order(SUBSUB$epi_date),]
  SUBSUB$Cases[is.na(SUBSUB$Cases)] = 0
  temp = create_weekly_dataframe(x=min(SUBSUB$epi_date,na.rm=T), y=max(SUBSUB$epi_date, na.rm=T))
  SUBSUB$date_char = as.character(SUBSUB$epi_date)
  SUBSUB = merge(SUBSUB, data.frame(date_char = as.character(temp$date), week_index = c(1:length(temp$date))),
                 by = 'date_char', all.x = T, all.y = T)
  ts = rep(0,max(SUBSUB$week_index) - min(SUBSUB$week_index) + 1)
  ts[SUBSUB$week_index - min(SUBSUB$week_index)+1] = SUBSUB$Cases
  ts[is.na(ts)]=0
  ts_dates = rep(NA,max(SUBSUB$week_index) - min(SUBSUB$week_index) + 1)
  ts_dates[SUBSUB$week_index - min(SUBSUB$week_index)+1] = as.character(SUBSUB$date_char)
  ts_dates_orig = rep(NA,max(SUBSUB$week_index) - min(SUBSUB$week_index) + 1)
  ts_dates_orig[SUBSUB$week_index - min(SUBSUB$week_index)+1] = as.character(SUBSUB$epi_date)
  inds_to_include = min(which(ts>0)):max(which(ts>0))
  if(nrow(SUBSUB)>= 10 & !is.na(var(SUBSUB$Cases, na.rm=T)) & var(SUBSUB$Cases, na.rm=T) > 0){
    num <- num + 1
    templist <- list(ts = pmax(0,ts)[inds_to_include],
                     ts_dates = ts_dates[inds_to_include],
                     ts_dates_actual = ts_dates_orig[inds_to_include],
                     ts_disease = 'InfluenzaA',
                     ts_measurement_type = "incidence",
                     ts_geography = trimws(LOCATIONS[l]),
                     ts_first_time = min(ts_dates[inds_to_include], na.rm=T),
                     ts_last_time = max(ts_dates[inds_to_include], na.rm=T),
                     ts_time_cadence = 'weekly',
                     ts_scale = "counts")
    data_list[[num]] <- templist
  }
}
saveRDS(data_list, paste0(savepath,"InfluenzaA_whoflunet.RDS"))



### Subtype B
dat$Cases = dat$INF_B
dat = dat %>% dplyr::group_by(COUNTRY_AREA_TERRITORY) %>% dplyr::mutate(num_nonzero = sum(Cases[!is.na(Cases) & Cases>0]))
dat_sub = dat[dat$num_nonzero > 50,]
LOC_TO_EXCLUDE = c('Kuwait', 'Malta', 'North Macedonia', 'Mauritius', 'Saudi Arabia', 'Zimbabwe')
SUB = dat_sub[!(dat_sub$COUNTRY_AREA_TERRITORY %in% LOC_TO_EXCLUDE),]
data_list = list()
num = 0
LOCATIONS = unique(SUB$COUNTRY_AREA_TERRITORY)
for(l in 1:length(LOCATIONS)){
  SUBSUB = SUB[SUB$COUNTRY_AREA_TERRITORY == LOCATIONS[l],]
  SUBSUB = SUBSUB[order(SUBSUB$epi_date),]
  SUBSUB$Cases[is.na(SUBSUB$Cases)] = 0
  temp = create_weekly_dataframe(x=min(SUBSUB$epi_date,na.rm=T), y=max(SUBSUB$epi_date, na.rm=T))
  SUBSUB$date_char = as.character(SUBSUB$epi_date)
  SUBSUB = merge(SUBSUB, data.frame(date_char = as.character(temp$date), week_index = c(1:length(temp$date))),
                 by = 'date_char', all.x = T, all.y = T)
  ts = rep(0,max(SUBSUB$week_index) - min(SUBSUB$week_index) + 1)
  ts[SUBSUB$week_index - min(SUBSUB$week_index)+1] = SUBSUB$Cases
  ts[is.na(ts)]=0
  ts_dates = rep(NA,max(SUBSUB$week_index) - min(SUBSUB$week_index) + 1)
  ts_dates[SUBSUB$week_index - min(SUBSUB$week_index)+1] = as.character(SUBSUB$date_char)
  ts_dates_orig = rep(NA,max(SUBSUB$week_index) - min(SUBSUB$week_index) + 1)
  ts_dates_orig[SUBSUB$week_index - min(SUBSUB$week_index)+1] = as.character(SUBSUB$epi_date)
  inds_to_include = min(which(ts>0)):max(which(ts>0))
  if(nrow(SUBSUB)>= 10 & !is.na(var(SUBSUB$Cases, na.rm=T)) & var(SUBSUB$Cases, na.rm=T) > 0){
    num <- num + 1
    templist <- list(ts = pmax(0,ts)[inds_to_include],
                     ts_dates = ts_dates[inds_to_include],
                     ts_dates_actual = ts_dates_orig[inds_to_include],
                     ts_disease = 'InfluenzaB',
                     ts_measurement_type = "incidence",
                     ts_geography = trimws(LOCATIONS[l]),
                     ts_first_time = min(ts_dates[inds_to_include], na.rm=T),
                     ts_last_time = max(ts_dates[inds_to_include], na.rm=T),
                     ts_time_cadence = 'weekly',
                     ts_scale = "counts")
    data_list[[num]] <- templist
  }
}
saveRDS(data_list, paste0(savepath,"InfluenzaB_whoflunet.RDS"))






### Subtype A
dat$Cases = dat$INF_A
dat = dat %>% dplyr::group_by(COUNTRY_AREA_TERRITORY) %>% dplyr::mutate(num_nonzero = sum(Cases[!is.na(Cases) & Cases>0]))
dat_sub = dat[dat$num_nonzero > 50,]
LOC_TO_EXCLUDE = c('Kuwait', 'Malta', 'North Macedonia', 'Mauritius', 'Saudi Arabia', 'Zimbabwe')
SUB = dat_sub[!(dat_sub$COUNTRY_AREA_TERRITORY %in% LOC_TO_EXCLUDE),]
data_list = list()
num = 0
LOCATIONS = unique(SUB$COUNTRY_AREA_TERRITORY)
for(l in 1:length(LOCATIONS)){
  SUBSUB = SUB[SUB$COUNTRY_AREA_TERRITORY == LOCATIONS[l],]
  SUBSUB = SUBSUB[order(SUBSUB$epi_date),]
  SUBSUB$Cases[is.na(SUBSUB$Cases)] = 0
  temp = create_weekly_dataframe(x=min(SUBSUB$epi_date,na.rm=T), y=max(SUBSUB$epi_date, na.rm=T))
  SUBSUB$date_char = as.character(SUBSUB$epi_date)
  SUBSUB = merge(SUBSUB, data.frame(date_char = as.character(temp$date), week_index = c(1:length(temp$date))),
                 by = 'date_char', all.x = T, all.y = T)
  ts = rep(0,max(SUBSUB$week_index) - min(SUBSUB$week_index) + 1)
  ts[SUBSUB$week_index - min(SUBSUB$week_index)+1] = SUBSUB$Cases
  ts[is.na(ts)]=0
  ts_dates = rep(NA,max(SUBSUB$week_index) - min(SUBSUB$week_index) + 1)
  ts_dates[SUBSUB$week_index - min(SUBSUB$week_index)+1] = as.character(SUBSUB$date_char)
  ts_dates_orig = rep(NA,max(SUBSUB$week_index) - min(SUBSUB$week_index) + 1)
  ts_dates_orig[SUBSUB$week_index - min(SUBSUB$week_index)+1] = as.character(SUBSUB$epi_date)
  inds_to_include = min(which(ts>0)):max(which(ts>0))
  if(nrow(SUBSUB)>= 10 & !is.na(var(SUBSUB$Cases, na.rm=T)) & var(SUBSUB$Cases, na.rm=T) > 0){
    num <- num + 1
    templist <- list(ts = pmax(0,ts)[inds_to_include],
                     ts_dates = ts_dates[inds_to_include],
                     ts_dates_actual = ts_dates_orig[inds_to_include],
                     ts_disease = 'InfluenzaA',
                     ts_measurement_type = "incidence",
                     ts_geography = trimws(LOCATIONS[l]),
                     ts_first_time = min(ts_dates[inds_to_include], na.rm=T),
                     ts_last_time = max(ts_dates[inds_to_include], na.rm=T),
                     ts_time_cadence = 'weekly',
                     ts_scale = "counts")
    data_list[[num]] <- templist
  }
}
saveRDS(data_list, paste0(savepath,"InfluenzaA_whoflunet.RDS"))


### Influenza Any
dat$Cases = dat$INF_ALL
dat = dat %>% dplyr::group_by(COUNTRY_AREA_TERRITORY) %>% dplyr::mutate(num_nonzero = sum(Cases[!is.na(Cases) & Cases>0]))
dat_sub = dat[dat$num_nonzero > 50,]
LOC_TO_EXCLUDE = c('Albania', 'Algeria' ,'Armenia','Croatia','Dominican Republic','Ecuador', 'Kuwait', 'Krygyzstan','Lithuania','Mauritius',
                   'Malta','North Macedonia')
SUB = dat_sub[!(dat_sub$COUNTRY_AREA_TERRITORY %in% LOC_TO_EXCLUDE),]
data_list = list()
num = 0
LOCATIONS = unique(SUB$COUNTRY_AREA_TERRITORY)
for(l in 1:length(LOCATIONS)){
  SUBSUB = SUB[SUB$COUNTRY_AREA_TERRITORY == LOCATIONS[l],]
  SUBSUB = SUBSUB[order(SUBSUB$epi_date),]
  SUBSUB$Cases[is.na(SUBSUB$Cases)] = 0
  temp = create_weekly_dataframe(x=min(SUBSUB$epi_date,na.rm=T), y=max(SUBSUB$epi_date, na.rm=T))
  SUBSUB$date_char = as.character(SUBSUB$epi_date)
  SUBSUB = merge(SUBSUB, data.frame(date_char = as.character(temp$date), week_index = c(1:length(temp$date))),
                 by = 'date_char', all.x = T, all.y = T)
  ts = rep(0,max(SUBSUB$week_index) - min(SUBSUB$week_index) + 1)
  ts[SUBSUB$week_index - min(SUBSUB$week_index)+1] = SUBSUB$Cases
  ts[is.na(ts)]=0
  ts_dates = rep(NA,max(SUBSUB$week_index) - min(SUBSUB$week_index) + 1)
  ts_dates[SUBSUB$week_index - min(SUBSUB$week_index)+1] = as.character(SUBSUB$date_char)
  ts_dates_orig = rep(NA,max(SUBSUB$week_index) - min(SUBSUB$week_index) + 1)
  ts_dates_orig[SUBSUB$week_index - min(SUBSUB$week_index)+1] = as.character(SUBSUB$epi_date)
  inds_to_include = min(which(ts>0)):max(which(ts>0))
  if(nrow(SUBSUB)>= 10 & !is.na(var(SUBSUB$Cases, na.rm=T)) & var(SUBSUB$Cases, na.rm=T) > 0){
    num <- num + 1
    templist <- list(ts = pmax(0,ts)[inds_to_include],
                     ts_dates = ts_dates[inds_to_include],
                     ts_dates_actual = ts_dates_orig[inds_to_include],
                     ts_disease = 'Influenza',
                     ts_measurement_type = "incidence",
                     ts_geography = trimws(LOCATIONS[l]),
                     ts_first_time = min(ts_dates[inds_to_include], na.rm=T),
                     ts_last_time = max(ts_dates[inds_to_include], na.rm=T),
                     ts_time_cadence = 'weekly',
                     ts_scale = "counts")
    data_list[[num]] <- templist
  }
}
saveRDS(data_list, paste0(savepath,"Influenza_whoflunet.RDS"))


  



### Bocavirus
dat$Cases = dat$BOCA
dat = dat %>% dplyr::group_by(COUNTRY_AREA_TERRITORY) %>% dplyr::mutate(num_nonzero = sum(Cases[!is.na(Cases) & Cases>0]))
dat_sub = dat[dat$num_nonzero > 50,]
LOC_TO_INCLUDE = c('Japan', 'Qatar')
SUB = dat_sub[(dat_sub$COUNTRY_AREA_TERRITORY %in% LOC_TO_EXCLUDE),]
data_list = list()
num = 0
LOCATIONS = unique(SUB$COUNTRY_AREA_TERRITORY)
for(l in 1:length(LOCATIONS)){
  SUBSUB = SUB[SUB$COUNTRY_AREA_TERRITORY == LOCATIONS[l],]
  SUBSUB = SUBSUB[order(SUBSUB$epi_date),]
  SUBSUB$Cases[is.na(SUBSUB$Cases)] = 0
  temp = create_weekly_dataframe(x=min(SUBSUB$epi_date,na.rm=T), y=max(SUBSUB$epi_date, na.rm=T))
  SUBSUB$date_char = as.character(SUBSUB$epi_date)
  SUBSUB = merge(SUBSUB, data.frame(date_char = as.character(temp$date), week_index = c(1:length(temp$date))),
                 by = 'date_char', all.x = T, all.y = T)
  ts = rep(0,max(SUBSUB$week_index) - min(SUBSUB$week_index) + 1)
  ts[SUBSUB$week_index - min(SUBSUB$week_index)+1] = SUBSUB$Cases
  ts[is.na(ts)]=0
  ts_dates = rep(NA,max(SUBSUB$week_index) - min(SUBSUB$week_index) + 1)
  ts_dates[SUBSUB$week_index - min(SUBSUB$week_index)+1] = as.character(SUBSUB$date_char)
  ts_dates_orig = rep(NA,max(SUBSUB$week_index) - min(SUBSUB$week_index) + 1)
  ts_dates_orig[SUBSUB$week_index - min(SUBSUB$week_index)+1] = as.character(SUBSUB$epi_date)
  inds_to_include = min(which(ts>0)):max(which(ts>0))
  if(nrow(SUBSUB)>= 10 & !is.na(var(SUBSUB$Cases, na.rm=T)) & var(SUBSUB$Cases, na.rm=T) > 0){
    num <- num + 1
    templist <- list(ts = pmax(0,ts)[inds_to_include],
                     ts_dates = ts_dates[inds_to_include],
                     ts_dates_actual = ts_dates_orig[inds_to_include],
                     ts_disease = 'Bocavirus',
                     ts_measurement_type = "incidence",
                     ts_geography = trimws(LOCATIONS[l]),
                     ts_first_time = min(ts_dates[inds_to_include], na.rm=T),
                     ts_last_time = max(ts_dates[inds_to_include], na.rm=T),
                     ts_time_cadence = 'weekly',
                     ts_scale = "counts")
    data_list[[num]] <- templist
  }
}
saveRDS(data_list, paste0(savepath,"Bocavirus_whoflunet.RDS"))






### Adenovirus
dat$Cases = dat$ADENO
dat = dat %>% dplyr::group_by(COUNTRY_AREA_TERRITORY) %>% dplyr::mutate(num_nonzero = sum(Cases[!is.na(Cases) & Cases>0]))
dat_sub = dat[dat$num_nonzero > 50,]
LOC_TO_INCLUDE = c('Argentina', 'Australia', 'Brazil', 'Canada', 'Chile',
                   'Colombia', 'Costa Rica', 'Japan', 'Qatar','Paraguay', 'Oman')
SUB = dat_sub[(dat_sub$COUNTRY_AREA_TERRITORY %in% LOC_TO_EXCLUDE),]
data_list = list()
num = 0
LOCATIONS = unique(SUB$COUNTRY_AREA_TERRITORY)
for(l in 1:length(LOCATIONS)){
  SUBSUB = SUB[SUB$COUNTRY_AREA_TERRITORY == LOCATIONS[l],]
  SUBSUB = SUBSUB[order(SUBSUB$epi_date),]
  SUBSUB$Cases[is.na(SUBSUB$Cases)] = 0
  temp = create_weekly_dataframe(x=min(SUBSUB$epi_date,na.rm=T), y=max(SUBSUB$epi_date, na.rm=T))
  SUBSUB$date_char = as.character(SUBSUB$epi_date)
  SUBSUB = merge(SUBSUB, data.frame(date_char = as.character(temp$date), week_index = c(1:length(temp$date))),
                 by = 'date_char', all.x = T, all.y = T)
  ts = rep(0,max(SUBSUB$week_index) - min(SUBSUB$week_index) + 1)
  ts[SUBSUB$week_index - min(SUBSUB$week_index)+1] = SUBSUB$Cases
  ts[is.na(ts)]=0
  ts_dates = rep(NA,max(SUBSUB$week_index) - min(SUBSUB$week_index) + 1)
  ts_dates[SUBSUB$week_index - min(SUBSUB$week_index)+1] = as.character(SUBSUB$date_char)
  ts_dates_orig = rep(NA,max(SUBSUB$week_index) - min(SUBSUB$week_index) + 1)
  ts_dates_orig[SUBSUB$week_index - min(SUBSUB$week_index)+1] = as.character(SUBSUB$epi_date)
  inds_to_include = min(which(ts>0)):max(which(ts>0))
  if(nrow(SUBSUB)>= 10 & !is.na(var(SUBSUB$Cases, na.rm=T)) & var(SUBSUB$Cases, na.rm=T) > 0){
    num <- num + 1
    templist <- list(ts = pmax(0,ts)[inds_to_include],
                     ts_dates = ts_dates[inds_to_include],
                     ts_dates_actual = ts_dates_orig[inds_to_include],
                     ts_disease = 'Adenovirus',
                     ts_measurement_type = "incidence",
                     ts_geography = trimws(LOCATIONS[l]),
                     ts_first_time = min(ts_dates[inds_to_include], na.rm=T),
                     ts_last_time = max(ts_dates[inds_to_include], na.rm=T),
                     ts_time_cadence = 'weekly',
                     ts_scale = "counts")
    data_list[[num]] <- templist
  }
}
saveRDS(data_list, paste0(savepath,"Adenovirus_whoflunet.RDS"))








### Metapneumovirus
dat$Cases = dat$METAPNEUMO
dat = dat %>% dplyr::group_by(COUNTRY_AREA_TERRITORY) %>% dplyr::mutate(num_nonzero = sum(Cases[!is.na(Cases) & Cases>0]))
dat_sub = dat[dat$num_nonzero > 50,]
LOC_TO_INCLUDE = c('Argentina', 'Australia', 'Brazil', 'Canada', 'Chile',
                   'Colombia', 'Costa Rica', 'Japan', 'Qatar','Paraguay')
SUB = dat_sub[(dat_sub$COUNTRY_AREA_TERRITORY %in% LOC_TO_EXCLUDE),]
data_list = list()
num = 0
LOCATIONS = unique(SUB$COUNTRY_AREA_TERRITORY)
for(l in 1:length(LOCATIONS)){
  SUBSUB = SUB[SUB$COUNTRY_AREA_TERRITORY == LOCATIONS[l],]
  SUBSUB = SUBSUB[order(SUBSUB$epi_date),]
  SUBSUB$Cases[is.na(SUBSUB$Cases)] = 0
  temp = create_weekly_dataframe(x=min(SUBSUB$epi_date,na.rm=T), y=max(SUBSUB$epi_date, na.rm=T))
  SUBSUB$date_char = as.character(SUBSUB$epi_date)
  SUBSUB = merge(SUBSUB, data.frame(date_char = as.character(temp$date), week_index = c(1:length(temp$date))),
                 by = 'date_char', all.x = T, all.y = T)
  ts = rep(0,max(SUBSUB$week_index) - min(SUBSUB$week_index) + 1)
  ts[SUBSUB$week_index - min(SUBSUB$week_index)+1] = SUBSUB$Cases
  ts[is.na(ts)]=0
  ts_dates = rep(NA,max(SUBSUB$week_index) - min(SUBSUB$week_index) + 1)
  ts_dates[SUBSUB$week_index - min(SUBSUB$week_index)+1] = SUBSUB$date_char
  ts_dates_orig = rep(NA,max(SUBSUB$week_index) - min(SUBSUB$week_index) + 1)
  ts_dates_orig[SUBSUB$week_index - min(SUBSUB$week_index)+1] = as.character(SUBSUB$epi_date)
  inds_to_include = min(which(ts>0)):max(which(ts>0))
  if(nrow(SUBSUB)>= 10 & !is.na(var(SUBSUB$Cases, na.rm=T)) & var(SUBSUB$Cases, na.rm=T) > 0){
    num <- num + 1
    templist <- list(ts = pmax(0,ts)[inds_to_include],
                     ts_dates = ts_dates[inds_to_include],
                     ts_dates_actual = ts_dates_orig[inds_to_include],
                     ts_disease = 'Metapneumovirus',
                     ts_measurement_type = "incidence",
                     ts_geography = trimws(LOCATIONS[l]),
                     ts_first_time = min(ts_dates[inds_to_include], na.rm=T),
                     ts_last_time = max(ts_dates[inds_to_include], na.rm=T),
                     ts_time_cadence = 'weekly',
                     ts_scale = "counts")
    data_list[[num]] <- templist
  }
}
saveRDS(data_list, paste0(savepath,"Metapneumovirus_whoflunet.RDS"))






### Parainfluenza
dat$Cases = dat$PARAINFLUENZA
dat = dat %>% dplyr::group_by(COUNTRY_AREA_TERRITORY) %>% dplyr::mutate(num_nonzero = sum(Cases[!is.na(Cases) & Cases>0]))
dat_sub = dat[dat$num_nonzero > 50,]
LOC_TO_INCLUDE = c('Argentina', 'Australia', 'Brazil', 'Canada', 'Chile',
                   'Colombia', 'Costa Rica', 'Japan','Mongolia', 'Qatar','Paraguay',
                   'Mexico')
SUB = dat_sub[(dat_sub$COUNTRY_AREA_TERRITORY %in% LOC_TO_EXCLUDE),]
data_list = list()
num = 0
LOCATIONS = unique(SUB$COUNTRY_AREA_TERRITORY)
for(l in 1:length(LOCATIONS)){
  SUBSUB = SUB[SUB$COUNTRY_AREA_TERRITORY == LOCATIONS[l],]
  SUBSUB = SUBSUB[order(SUBSUB$epi_date),]
  SUBSUB$Cases[is.na(SUBSUB$Cases)] = 0
  temp = create_weekly_dataframe(x=min(SUBSUB$epi_date,na.rm=T), y=max(SUBSUB$epi_date, na.rm=T))
  SUBSUB$date_char = as.character(SUBSUB$epi_date)
  SUBSUB = merge(SUBSUB, data.frame(date_char = as.character(temp$date), week_index = c(1:length(temp$date))),
                 by = 'date_char', all.x = T, all.y = T)
  ts = rep(0,max(SUBSUB$week_index) - min(SUBSUB$week_index) + 1)
  ts[SUBSUB$week_index - min(SUBSUB$week_index)+1] = SUBSUB$Cases
  ts[is.na(ts)]=0
  ts_dates = rep(NA,max(SUBSUB$week_index) - min(SUBSUB$week_index) + 1)
  ts_dates[SUBSUB$week_index - min(SUBSUB$week_index)+1] = SUBSUB$date_char
  ts_dates_orig = rep(NA,max(SUBSUB$week_index) - min(SUBSUB$week_index) + 1)
  ts_dates_orig[SUBSUB$week_index - min(SUBSUB$week_index)+1] = as.character(SUBSUB$epi_date)
  inds_to_include = min(which(ts>0)):max(which(ts>0))
  if(nrow(SUBSUB)>= 10 & !is.na(var(SUBSUB$Cases, na.rm=T)) & var(SUBSUB$Cases, na.rm=T) > 0){
    num <- num + 1
    templist <- list(ts = pmax(0,ts)[inds_to_include],
                     ts_dates = ts_dates[inds_to_include],
                     ts_dates_actual = ts_dates_orig[inds_to_include],
                     ts_disease = 'Parainfluenza',
                     ts_measurement_type = "incidence",
                     ts_geography = trimws(LOCATIONS[l]),
                     ts_first_time = min(ts_dates[inds_to_include], na.rm=T),
                     ts_last_time = max(ts_dates[inds_to_include], na.rm=T),
                     ts_time_cadence = 'weekly',
                     ts_scale = "counts")
    data_list[[num]] <- templist
  }
}
saveRDS(data_list, paste0(savepath,"Parainfluenza_whoflunet.RDS"))




### Rhinovirus
dat$Cases = dat$RHINO
dat = dat %>% dplyr::group_by(COUNTRY_AREA_TERRITORY) %>% dplyr::mutate(num_nonzero = sum(Cases[!is.na(Cases) & Cases>0]))
dat_sub = dat[dat$num_nonzero > 50,]
LOC_TO_INCLUDE = c('Australia', 'Brazil', 'Canada', 'Chile',
                   'Colombia', 'Costa Rica', 'Japan','Mongolia', 'Oman','Qatar',
                   'Mexico')
SUB = dat_sub[(dat_sub$COUNTRY_AREA_TERRITORY %in% LOC_TO_EXCLUDE),]
data_list = list()
num = 0
LOCATIONS = unique(SUB$COUNTRY_AREA_TERRITORY)
for(l in 1:length(LOCATIONS)){
  SUBSUB = SUB[SUB$COUNTRY_AREA_TERRITORY == LOCATIONS[l],]
  SUBSUB = SUBSUB[order(SUBSUB$epi_date),]
  SUBSUB$Cases[is.na(SUBSUB$Cases)] = 0
  temp = create_weekly_dataframe(x=min(SUBSUB$epi_date,na.rm=T), y=max(SUBSUB$epi_date, na.rm=T))
  SUBSUB$date_char = as.character(SUBSUB$epi_date)
  SUBSUB = merge(SUBSUB, data.frame(date_char = as.character(temp$date), week_index = c(1:length(temp$date))),
                 by = 'date_char', all.x = T, all.y = T)
  ts = rep(0,max(SUBSUB$week_index) - min(SUBSUB$week_index) + 1)
  ts[SUBSUB$week_index - min(SUBSUB$week_index)+1] = SUBSUB$Cases
  ts[is.na(ts)]=0
  ts_dates = rep(NA,max(SUBSUB$week_index) - min(SUBSUB$week_index) + 1)
  ts_dates[SUBSUB$week_index - min(SUBSUB$week_index)+1] = as.character(SUBSUB$date_char)
  ts_dates_orig = rep(NA,max(SUBSUB$week_index) - min(SUBSUB$week_index) + 1)
  ts_dates_orig[SUBSUB$week_index - min(SUBSUB$week_index)+1] = as.character(SUBSUB$epi_date)
  inds_to_include = min(which(ts>0)):max(which(ts>0))
  if(nrow(SUBSUB)>= 10 & !is.na(var(SUBSUB$Cases, na.rm=T)) & var(SUBSUB$Cases, na.rm=T) > 0){
    num <- num + 1
    templist <- list(ts = pmax(0,ts)[inds_to_include],
                     ts_dates = ts_dates[inds_to_include],
                     ts_dates_actual = ts_dates_orig[inds_to_include],
                     ts_disease = 'Rhinovirus',
                     ts_measurement_type = "incidence",
                     ts_geography = trimws(LOCATIONS[l]),
                     ts_first_time = min(ts_dates[inds_to_include], na.rm=T),
                     ts_last_time = max(ts_dates[inds_to_include], na.rm=T),
                     ts_time_cadence = 'weekly',
                     ts_scale = "counts")
    data_list[[num]] <- templist
  }
}
saveRDS(data_list, paste0(savepath,"Rhinovirus_whoflunet.RDS"))



### RSV
dat$Cases = dat$RSV
dat = dat %>% dplyr::group_by(COUNTRY_AREA_TERRITORY) %>% dplyr::mutate(num_nonzero = sum(Cases[!is.na(Cases) & Cases>0]))
dat_sub = dat[dat$num_nonzero > 50,]
SUB = dat_sub
data_list = list()
num = 0
LOCATIONS = unique(SUB$COUNTRY_AREA_TERRITORY)
for(l in 1:length(LOCATIONS)){
  SUBSUB = SUB[SUB$COUNTRY_AREA_TERRITORY == LOCATIONS[l],]
  SUBSUB = SUBSUB[order(SUBSUB$epi_date),]
  SUBSUB$Cases[is.na(SUBSUB$Cases)] = 0
  temp = create_weekly_dataframe(x=min(SUBSUB$epi_date,na.rm=T), y=max(SUBSUB$epi_date, na.rm=T))
  SUBSUB$date_char = as.character(SUBSUB$epi_date)
  SUBSUB = merge(SUBSUB, data.frame(date_char = as.character(temp$date), week_index = c(1:length(temp$date))),
                 by = 'date_char', all.x = T, all.y = T)
  ts = rep(0,max(SUBSUB$week_index) - min(SUBSUB$week_index) + 1)
  ts[SUBSUB$week_index - min(SUBSUB$week_index)+1] = SUBSUB$Cases
  ts[is.na(ts)]=0
  ts_dates = rep(NA,max(SUBSUB$week_index) - min(SUBSUB$week_index) + 1)
  ts_dates[SUBSUB$week_index - min(SUBSUB$week_index)+1] = as.character(SUBSUB$date_char)
  ts_dates_orig = rep(NA,max(SUBSUB$week_index) - min(SUBSUB$week_index) + 1)
  ts_dates_orig[SUBSUB$week_index - min(SUBSUB$week_index)+1] = as.character(SUBSUB$epi_date)
  inds_to_include = min(which(ts>0)):max(which(ts>0))
  if(nrow(SUBSUB)>= 10 & !is.na(var(SUBSUB$Cases, na.rm=T)) & var(SUBSUB$Cases, na.rm=T) > 0){
    num <- num + 1
    templist <- list(ts = pmax(0,ts)[inds_to_include],
                     ts_dates = ts_dates[inds_to_include],
                     ts_dates_actual = ts_dates_orig[inds_to_include],
                     ts_disease = 'RSV',
                     ts_measurement_type = "incidence",
                     ts_geography = trimws(LOCATIONS[l]),
                     ts_first_time = min(ts_dates[inds_to_include], na.rm=T),
                     ts_last_time = max(ts_dates[inds_to_include], na.rm=T),
                     ts_time_cadence = 'weekly',
                     ts_scale = "counts")
    data_list[[num]] <- templist
  }
}
saveRDS(data_list, paste0(savepath,"RSV_whoflunet.RDS"))


