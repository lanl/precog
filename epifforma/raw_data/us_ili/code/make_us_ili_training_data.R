## Dave Osthus
## 1-31-24
## Download and organize ILI data for training in epiFFORMA

###############################################################
## NOTE: I downloaded the raw data read in within this script from here: https://gis.cdc.gov/grasp/fluview/fluportaldashboard.html
## I would have used the cdcfluview R package, but it does not appear to be operational
## Also note that MMWR weeks run from Sunday to the following Saturday, as defined here: https://ndc.services.cdc.gov/wp-content/uploads/MMWR_Week_overview.pdf

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

###############################################################
#### ILINet DATA
## get national data
natdat <- read.csv(paste0(datapath,"nat_ili/ILINet.csv"),skip=1)

## get hhs data
hhsdat <- read.csv(paste0(datapath,"hhs_ili/ILINet.csv"),skip=1)

## get state data
statedat <- read.csv(paste0(datapath,"state_ili/ILINet.csv"),skip=1)

## get census data
censusdat <- read.csv(paste0(datapath,"census_ili/ILINet.csv"),skip=1)

## combine data
dat        <- rbind(natdat,hhsdat,statedat,censusdat)
names(dat) <- tolower(names(dat))
dat        <- dat[order(dat$region.type,dat$region,dat$year,dat$week),]
  
#### WHO_NREVSS data
## get national data
## tempnatwhodat <- who_nrevss(region = "national",years = NULL)
####

## make epi week timing
dat$epi_date   <- MMWRweek2Date(MMWRyear = dat$year, MMWRweek = dat$week)
dat$epi_time   <- 0
dat[dat$week==30,]$epi_time <- 1
dat$epi_season <- dat$year
dat            <- dat[order(dat$epi_date),]
unqdates <- sort(unique(dat$epi_date))
for(i in 2:length(unqdates)){
  print(i)
  tempdat <- subset(dat,epi_date==unqdates[i])
  if(tempdat$epi_time[1] != 1){
    dat[dat$epi_date == unqdates[i],]$epi_time   <- dat[dat$epi_date == unqdates[i-1],]$epi_time[1] + 1
    dat[dat$epi_date == unqdates[i],]$epi_season <- dat[dat$epi_date == unqdates[i-1],]$epi_season[1]
  }
}

### rename dat
names(dat)[names(dat)=="year"]            <- "epi_year"
names(dat)[names(dat)=="week"]            <- "epi_week"
names(dat)[names(dat)=="x..weighted.ili"]    <- "wili"
dat$wili <- as.numeric(dat$wili)/100
names(dat)[names(dat)=="x.unweighted.ili"]  <- "ili"
dat$ili <- as.numeric(dat$ili)/100

## replace spaces with underscores
dat$region <- gsub("\\ ","_",dat$region)
  
## final formatting
dat       <- dat[order(dat$epi_date),]
regionvec <- c("nat",paste("hhs",1:10,sep=""))
tempew    <- formatC(dat[nrow(dat),]$epi_week,width=2,flag=0)
tempyr    <- dat[nrow(dat),]$epi_year

dat[which.min(dat$epi_date),]
dat[which.max(dat$epi_date),]

A = dat[dat$region.type != 'Census Regions' & dat$region.type != 'HHS Regions' & dat$region.type != 'National',]
A[which.min(A$epi_date),]
A[which.max(A$epi_date),]



A = dat[dat$region.type == 'Census Regions',]
A[which.min(A$epi_date),]
A[which.max(A$epi_date),]


A = dat[dat$region.type == 'HHS Regions',]
A[which.min(A$epi_date),]
A[which.max(A$epi_date),]

## prepare ili data for training data
unq_state_season <- ddply(dat, .(region,epi_season), summarise, max_epi_time = max(epi_time))
us_ili_training_data <- list()
cnt <- 0
for(i in 1:nrow(unq_state_season)){
  print(i)
  tempdf <- subset(dat, region == unq_state_season$region[i] & epi_season == unq_state_season$epi_season[i] & !is.na(ili))
  if(nrow(tempdf)>=10 & var(tempdf$ili, na.rm=T) > 0){
    cnt <- cnt + 1
    templist <- list(ts = tempdf$ili,
                     ts_dates = tempdf$epi_date,
                     ts_exogenous = NULL,
                     ts_real_data = T,
                     ts_isolated_strain = F,
                     ts_multiwave = F,
                     ts_disease = "influenza-like illness",
                     ts_measurement_type = "incidence",
                     ts_geography = tempdf$region[1],
                     ts_first_time = min(tempdf$epi_date),
                     ts_last_time = max(tempdf$epi_date),
                     ts_time_cadence = "weekly",
                     ts_scale = "proportion",
                     ts_exogenous_scale = NA, 
                     ts_seasonal = F)
    us_ili_training_data[[cnt]] <- templist
  }
}

## weighted ILI (not available for states)
us_wili_training_data <- list()
cnt <- 0
for(i in 1:nrow(unq_state_season)){
  print(i)
  tempdf <- subset(dat, region == unq_state_season$region[i] & epi_season == unq_state_season$epi_season[i] & !is.na(wili))
  if(nrow(tempdf)>= 10  & var(tempdf$wili, na.rm=T) > 0){
    cnt <- cnt + 1
    templist <- list(ts = tempdf$wili,
                     ts_dates = tempdf$epi_date,
                     ts_exogenous = NULL,
                     ts_real_data = T,
                     ts_isolated_strain = F,
                     ts_multiwave = F,
                     ts_disease = "weighted influenza-like illness",
                     ts_measurement_type = "incidence",
                     ts_geography = tempdf$region[1],
                     ts_first_time = min(tempdf$epi_date),
                     ts_last_time = max(tempdf$epi_date),
                     ts_time_cadence = "weekly",
                     ts_scale = "proportion",
                     ts_exogenous_scale = NA,
                     ts_seasonal = F)
    
    us_wili_training_data[[cnt]] <- templist
  }
}


####################################################  
## make wili and ili by region (not season), mimicing a roller coaster
## weighted ILI (not available for states)

unqstate <- setdiff(sort(unique(unq_state_season$region)),"X")

## ILI rollercoaster
us_ili_rollercoaster_training_data <- list()
cnt <- 0
for(i in 1:length(unqstate)){
  print(i)
  tempdf <- subset(dat, region == unqstate[i] & !is.na(ili))
  tempdf <- tempdf[order(tempdf$epi_date),]
  if(nrow(tempdf)>= 10  & var(tempdf$ili, na.rm=T) > 0){
    cnt <- cnt + 1
    templist <- list(ts = tempdf$ili,
                     ts_dates = tempdf$epi_date,
                     ts_exogenous = NULL,
                     ts_real_data = T,
                     ts_isolated_strain = F,
                     ts_multiwave = T,
                     ts_disease = "influenza-like illness",
                     ts_measurement_type = "incidence",
                     ts_geography = tempdf$region[1],
                     ts_first_time = min(tempdf$epi_date),
                     ts_last_time = max(tempdf$epi_date),
                     ts_time_cadence = "weekly",
                     ts_scale = "proportion",
                     ts_exogenous_scale = NA,
                     ts_seasonal = F)
    
    us_ili_rollercoaster_training_data[[cnt]] <- templist
  }
}

## wILI rollercoaster
us_wili_rollercoaster_training_data <- list()
cnt <- 0
for(i in 1:length(unqstate)){
  print(i)
  tempdf <- subset(dat, region == unqstate[i] & !is.na(wili))
  tempdf <- tempdf[order(tempdf$epi_date),]
  if(nrow(tempdf)>= 10  & var(tempdf$wili, na.rm=T) > 0){
    cnt <- cnt + 1
    templist <- list(ts = tempdf$wili,
                     ts_dates = tempdf$epi_date,
                     ts_exogenous = NULL,
                     ts_real_data = T,
                     ts_isolated_strain = F,
                     ts_multiwave = T,
                     ts_disease = "weighted influenza-like illness",
                     ts_measurement_type = "incidence",
                     ts_geography = tempdf$region[1],
                     ts_first_time = min(tempdf$epi_date),
                     ts_last_time = max(tempdf$epi_date),
                     ts_time_cadence = "weekly",
                     ts_scale = "proportion",
                     ts_exogenous_scale = NA, 
                     ts_seasonal = F)
    
    us_wili_rollercoaster_training_data[[cnt]] <- templist
  }
}


####################################################  
#### WHO NREVSS DATA
## get national data
natdatwhoclinic <- read.csv(paste0(datapath,"nat_ili/WHO_NREVSS_Clinical_Labs.csv"),skip=1)
hhsdatwhoclinic <- read.csv(paste0(datapath,"hhs_ili/WHO_NREVSS_Clinical_Labs.csv"),skip=1)
statedatwhoclinic <- read.csv(paste0(datapath,"state_ili/WHO_NREVSS_Clinical_Labs.csv"),skip=1)
censusdatwhoclinic <- read.csv(paste0(datapath,"census_ili/WHO_NREVSS_Clinical_Labs.csv"),skip=1)

## combine the WHO clinical lab files
datwhoclinic <- rbind(natdatwhoclinic,
                      hhsdatwhoclinic,
                      statedatwhoclinic,
                      censusdatwhoclinic)
names(datwhoclinic) <- tolower(names(datwhoclinic))
datwhoclinic <- subset(datwhoclinic, region != "X")

## merge dat with datwhoclinic
dat_big <- merge(dat, 
                 subset(datwhoclinic,select=c("region.type","region","year","week","percent.positive","percent.a","percent.b")),
                 by.x=c("region.type","region","epi_year","epi_week"),
                 by.y=c("region.type","region","year","week"),
                 all.x=T)
dat_big$percent.a <- as.numeric(dat_big$percent.a)
dat_big$percent.b <- as.numeric(dat_big$percent.b)
dat_big$percent.positive <- as.numeric(dat_big$percent.positive)
dat_big <- dat_big[!is.na(dat_big$percent.positive),]
dat_big <- dat_big[!is.na(dat_big$percent.a),]
dat_big <- dat_big[!is.na(dat_big$percent.b),]

dat_big$ili_a <- dat_big$ili*(dat_big$percent.a/dat_big$percent.positive)
dat_big$ili_b <- dat_big$ili*(dat_big$percent.b/dat_big$percent.positive)
dat_big[dat_big$percent.positive == 0,]$ili_a <- 0
dat_big[dat_big$percent.positive == 0,]$ili_b <- 0

## prepare ili strain data for training data
unq_state_season_ab <- ddply(dat_big, .(region,epi_season), summarise, max_epi_time = max(epi_time))

## Strain A
us_ilia_training_data <- list()
cnt <- 0
for(i in 1:nrow(unq_state_season_ab)){
  print(i)
  tempdf <- subset(dat_big, region == unq_state_season_ab$region[i] & epi_season == unq_state_season_ab$epi_season[i] & !is.na(ili_a))
  if(nrow(tempdf) >= 10 & var(tempdf$ili_a) > 0){
    cnt <- cnt + 1
    templist <- list(ts = tempdf$ili_a,
                     ts_dates = tempdf$epi_date,
                     ts_exogenous = tempdf$percent.a/100,
                     ts_real_data = T,
                     ts_isolated_strain = T,
                     ts_multiwave = F,
                     ts_disease = "influenza-like illness strain A",
                     ts_measurement_type = "incidence",
                     ts_geography = tempdf$region[1],
                     ts_first_time = min(tempdf$epi_date),
                     ts_last_time = max(tempdf$epi_date),
                     ts_time_cadence = "weekly",
                     ts_scale = "proportion",
                     ts_exogenous_scale = "proportion",
                     ts_seasonal = F)
    us_ilia_training_data[[cnt]] <- templist
  }
}

## Strain B
us_ilib_training_data <- list()
cnt <- 0
for(i in 1:nrow(unq_state_season_ab)){
  print(i)
  tempdf <- subset(dat_big, region == unq_state_season_ab$region[i] & epi_season == unq_state_season_ab$epi_season[i] & !is.na(ili_b))
  if(nrow(tempdf)>= 10 & var(tempdf$ili_b) > 0){
    cnt <- cnt + 1
    templist <- list(ts = tempdf$ili_b,
                     ts_dates = tempdf$epi_date,
                     ts_exogenous = tempdf$percent.b/100,
                     ts_real_data = T,
                     ts_isolated_strain = T,
                     ts_multiwave = F,
                     ts_disease = "influenza-like illness strain B",
                     ts_measurement_type = "incidence",
                     ts_geography = tempdf$region[1],
                     ts_first_time = min(tempdf$epi_date),
                     ts_last_time = max(tempdf$epi_date),
                     ts_time_cadence = "weekly",
                     ts_scale = "proportion",
                     ts_exogenous_scale = "proportion",
                     ts_seasonal = F)
    us_ilib_training_data[[cnt]] <- templist
  }
}


## concatenate all the lists
us_flu_training_data <- c(us_ili_training_data,
                          us_wili_training_data,
                          us_ili_rollercoaster_training_data,
                          us_wili_rollercoaster_training_data,
                          us_ilia_training_data,
                          us_ilib_training_data)

## save the us ili training data
saveRDS(us_flu_training_data, paste0(savepath,"us_ili.RDS"))


# ## save all the lists as real ILI data time series
# tempid <- sample(1:length(us_flu_training_data),1)
# ttt <- us_flu_training_data[[tempid]]
# qplot(ttt$ts_dates, ttt$ts, geom=c("point","line"))+
#   ggtitle(paste0(ttt$description$ts_disease,": ",ttt$description$ts_geography,", ",ttt$description$ts_first_time," through ",ttt$description$ts_last_time))
# tempid





