############################
############################
### Clean US FluNet Data ###
############################
############################



library(MMWRweek)
library(lubridate)
library(ggplot2)
library(plyr)
theme_set(theme_bw())

datapath <- "./raw_data/"
savepath <- "../_Organized_Lists/"

####################
### Read in Data ###
####################

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


######################
### Organize Lists ###
######################


unqstate <- setdiff(sort(unique(dat$region)),"X")

us_ili_training_data <- list()
cnt <- 0
for(i in 1:length(unqstate)){
  print(i)
  tempdf <- subset(dat, region == unqstate[i] & !is.na(ili))
  tempdf <- tempdf[order(tempdf$epi_date),]
  if(nrow(tempdf)>= 10  & var(tempdf$ili, na.rm=T) > 0){
    cnt <- cnt + 1
    templist <- list(ts = tempdf$ili,
                     ts_dates = tempdf$epi_date,
                     ts_dates_actual = tempdf$epi_date,
                     ts_disease = "influenza-like illness",
                     ts_measurement_type = "incidence",
                     ts_geography = tempdf$region[1],
                     ts_first_time = min(tempdf$epi_date),
                     ts_last_time = max(tempdf$epi_date),
                     ts_time_cadence = "weekly",
                     ts_scale = "proportion")
    
    us_ili_training_data[[cnt]] <- templist
  }
}


us_wili_training_data <- list()
cnt <- 0
for(i in 1:length(unqstate)){
  print(i)
  tempdf <- subset(dat, region == unqstate[i] & !is.na(wili))
  tempdf <- tempdf[order(tempdf$epi_date),]
  if(nrow(tempdf)>= 10  & var(tempdf$wili, na.rm=T) > 0){
    cnt <- cnt + 1
    templist <- list(ts = tempdf$wili,
                     ts_dates = tempdf$epi_date,
                     ts_dates_actual = tempdf$epi_date,
                     ts_disease = "weighted influenza-like illness",
                     ts_measurement_type = "incidence",
                     ts_geography = tempdf$region[1],
                     ts_first_time = min(tempdf$epi_date),
                     ts_last_time = max(tempdf$epi_date),
                     ts_time_cadence = "weekly",
                     ts_scale = "proportion")
    
    us_wili_training_data[[cnt]] <- templist
  }
}


## concatenate all the lists
us_flu_training_data <- c(us_ili_training_data,
                          us_wili_training_data)

## save the us ili training data
saveRDS(us_flu_training_data, paste0(savepath,"Influenza_usflunet.RDS"))







###############################
### Read in WHO NREVSS Data ###
###############################

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
unq_state_ab <- ddply(dat_big, .(region), summarise, max_epi_time = max(epi_time))
dat_big = dat_big[order(dat_big$epi_date, dat_big$region),]

## Strain A
us_ilia_training_data <- list()
cnt <- 0
for(i in 1:nrow(unq_state_ab)){
  print(i)
  tempdf <- subset(dat_big, region == unq_state_ab$region[i] & !is.na(ili_a))
  if(nrow(tempdf) >= 10 & var(tempdf$ili_a) > 0){
    cnt <- cnt + 1
    templist <- list(ts = tempdf$ili_a,
                     ts_dates = tempdf$epi_date,
                     ts_dates_actual = tempdf$epi_date,
                     ts_disease = "influenza-like illness strain A",
                     ts_measurement_type = "incidence",
                     ts_geography = tempdf$region[1],
                     ts_first_time = min(tempdf$epi_date),
                     ts_last_time = max(tempdf$epi_date),
                     ts_time_cadence = "weekly",
                     ts_scale = "proportion")
    us_ilia_training_data[[cnt]] <- templist
  }
}

## Strain B
us_ilib_training_data <- list()
cnt <- 0
for(i in 1:nrow(unq_state_ab)){
  print(i)
  tempdf <- subset(dat_big, region == unq_state_ab$region[i] & !is.na(ili_b))
  if(nrow(tempdf)>= 10 & var(tempdf$ili_b) > 0){
    cnt <- cnt + 1
    templist <- list(ts = tempdf$ili_b,
                     ts_dates = tempdf$epi_date,
                     ts_dates_actual = tempdf$epi_date,
                     ts_disease = "influenza-like illness strain B",
                     ts_measurement_type = "incidence",
                     ts_geography = tempdf$region[1],
                     ts_first_time = min(tempdf$epi_date),
                     ts_last_time = max(tempdf$epi_date),
                     ts_time_cadence = "weekly",
                     ts_scale = "proportion")
    us_ilib_training_data[[cnt]] <- templist
  }
}


saveRDS(us_ilia_training_data, paste0(savepath,"InfluenzaA_usflunet.RDS"))
saveRDS(us_ilib_training_data, paste0(savepath,"InfluenzaB_usflunet.RDS"))



