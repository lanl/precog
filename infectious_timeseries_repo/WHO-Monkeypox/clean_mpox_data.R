################################
################################
### Clean WHO Monkeypox Data ###
################################
################################


library(MMWRweek)
library(lubridate)
library(ggplot2)
library(plyr)
library(lubridate)
theme_set(theme_bw())
library(dplyr)

datapath <- "./raw_data/" 
savepath <- "../_Organized_Lists/"


####################
### Read in Data ###
####################


dat = data.frame(data.table::fread(paste0(datapath,'Country data by date of symptom onset.csv')))
dat$reference_date = as.Date(dat$reference_date, format = '%Y-%m-%d')
dat = dat[dat$date_type == 'Onset',]

dat = dat %>% dplyr::group_by(country) %>% dplyr::mutate(num_nonzero = sum(cases[!is.na(cases) & cases>0]))
dat_sub = dat[dat$num_nonzero > 50,]

dat_sub = dat_sub[order(dat_sub$reference_date),]
#diff(unique(dat_sub$reference_date)), no missing dates overall (could be missing for individual locations)


dat_sub       <- dat_sub[order(dat_sub$reference_date),]
dat_sub = data.frame(dat_sub)

unqstate <- sort(unique(dat_sub$country))



 
######################
### Organize Lists ###
######################


create_weekly_dataframe <- function(x, y) {
  seq_dates <- seq(x, y + 7, by = "week")
  subset(data.frame(date = seq_dates), date <= y)
}


unqstate <- sort(unique(dat_sub$country))

mpox_training_data <- list()
cnt <- 0
for(i in 1:length(unqstate)){
  print(i)
  tempdf <- subset(dat_sub, country == unqstate[i])
  tempdf <- tempdf[order(tempdf$reference_date),]
  
  temp = create_weekly_dataframe(x=min(tempdf$reference_date,na.rm=T), y=max(tempdf$reference_date, na.rm=T))
  
  
  tempdf$date_char = as.character(tempdf$reference_date)
  tempdf = merge(tempdf, data.frame(date_char = as.character(unique(temp$date)), week_index = c(1:length(unique(temp$date)))),
                  by = 'date_char', all.x = T, all.y = T)
  tempdf <- tempdf[order(tempdf$week_index),]
  
  
  if(nrow(tempdf)>= 10  & var(tempdf$cases, na.rm=T) > 0){
    cnt <- cnt + 1
    ts = rep(0,max(tempdf$week_index) - min(tempdf$week_index) + 1)
    ts[tempdf$week_index - min(tempdf$week_index)+1] = tempdf$cases
    ts[is.na(ts)]=0
    
    ts_dates = rep(NA,max(tempdf$week_index) - min(tempdf$week_index) + 1)
    ts_dates[tempdf$week_index - min(tempdf$week_index)+1] = as.character(tempdf$date_char)
    
    ts_dates_actual = rep(NA,max(tempdf$week_index) - min(tempdf$week_index) + 1)
    ts_dates_actual[tempdf$week_index - min(tempdf$week_index)+1] = as.character(tempdf$reference_date)
    
    templist <- list(ts = ts,
                     ts_dates = ts_dates,
                     ts_dates_actual = ts_dates_actual,
                     ts_disease = "mpox",
                     ts_measurement_type = "incidence",
                     ts_geography = tempdf$region[1],
                     ts_first_time = min(tempdf$reference_date),
                     ts_last_time = max(tempdf$reference_date),
                     ts_time_cadence = "weekly",
                     ts_scale = "counts")
    
    mpox_training_data[[cnt]] <- templist
  }
}


## save the us ili training data
saveRDS(mpox_training_data, paste0(savepath,"Mpox_who.RDS"))







