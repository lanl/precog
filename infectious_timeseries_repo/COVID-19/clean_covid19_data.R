###########################
###########################
### Clean COVID-19 Data ###
###########################
###########################

library(MMWRweek)
library(lubridate)
library(ggplot2)
library(plyr)
library(lubridate)
theme_set(theme_bw())
datapath <- "./raw_data/"
savepath <- "../_Organized_Lists/"


########################
### Read in JHU Data ###
########################

Hopkins_LUT_All = read.csv(paste0(datapath,'COVID-19_LUT.csv'))
Hopkins_Static_All = readRDS(paste0(datapath,'COVID-19_Static.rds'))
Hopkins_Cases_All = readRDS(paste0(datapath,'COVID-19.rds'))

#Subset to admin1 and admin0 except for US
Hopkins_LUT = Hopkins_LUT_All[(Hopkins_LUT_All$Admin2 == '' & Hopkins_LUT_All$Admin3 == '') | (Hopkins_LUT_All$Admin0 == 'United States' & Hopkins_LUT_All$Admin2 != '' & Hopkins_LUT_All$Admin3 == ''),]
Hopkins_LUT$Admin01 = paste0(Hopkins_LUT$Admin0,'_',Hopkins_LUT$Admin1) 
Hopkins_LUT$Admin012 = paste0(Hopkins_LUT$Admin0,'_',Hopkins_LUT$Admin1, '_',Hopkins_LUT$Admin2) 

#Remove unassigned areas (note: this will not exclude their cases and deaths from being counted at the higher level Admins)
Hopkins_LUT = Hopkins_LUT[!(Hopkins_LUT$Admin1 %in% c('Unassigned')),]
Hopkins_LUT = Hopkins_LUT[!(Hopkins_LUT$Admin2 %in% c('Unassigned')),]

#Remove "out of" state listings  (note: this will not exclude their cases and deaths from being counted at the higher level Admins)
Hopkins_LUT = Hopkins_LUT[!grepl('Out of',Hopkins_LUT$Admin2),]

#Remove cruise ships
Hopkins_LUT = Hopkins_LUT[!(Hopkins_LUT$Level %in% c('Cruise Ship')),]

#Remove Namibia, not used elsewhere in JHU dataset
Hopkins_LUT = Hopkins_LUT[!is.na(Hopkins_LUT$ID),] 

#Remove Olympics 
Hopkins_LUT = Hopkins_LUT[!(Hopkins_LUT$ID %in% c('OLYMPICS')),]

#Remove Antarctica 
Hopkins_LUT = Hopkins_LUT[!(Hopkins_LUT$ID %in% c('AQ')),]

#Denmark, Netherlands, and France all have duplicate province listings. Remove
Hopkins_LUT = Hopkins_LUT[Hopkins_LUT$ID != 'FRFR',]
Hopkins_LUT = Hopkins_LUT[Hopkins_LUT$ID != 'DKDK',]
Hopkins_LUT = Hopkins_LUT[Hopkins_LUT$ID != 'NLNL',]

### Merge Static Data
Hopkins_LUT = Hopkins_LUT[,c('ID', 'Level', 'Admin0', 'Admin1', 'Admin2', 'Admin01', 'Admin012', 'Population', 'Longitude', 'Latitude',
                             'ISO1_3N','ISO1_3C','ISO1_2C' ,'ISO2','ISO2_UID','FIPS','NUTS','NameID' )]
Hopkins_LUT = merge(Hopkins_LUT,Hopkins_Static_All[,c('ID', 'Diabetes', 'Obesity', 'Smoking', 'COPD', 'CVD','HIV', 'Hypertension', 'WorldPop_Density', 'Access_City', 'Access_Motor', 'Access_Walk', 'WorldPop_65', 'WorldPop')], by = 'ID', all.x = T, all.y = F )

### Use JHU's Listed Population Unless NA
Hopkins_LUT$Population[is.na(Hopkins_LUT$Population)] = Hopkins_LUT$WorldPop[is.na(Hopkins_LUT$Population)]

#Hopkins_LUT[duplicated(Hopkins_LUT$ID),]
#confirmed there are no duplicate rows

### Format Cases and Deaths
Hopkins_Cases = Hopkins_Cases_All[Hopkins_Cases_All$ID %in% Hopkins_LUT$ID,]
Hopkins_Cases$Date = as.Date(Hopkins_Cases$Date, format = '%Y-%m-%d')

CASES = Hopkins_Cases[Hopkins_Cases$Type == 'Confirmed',]
CASES$ID_Date = paste0(CASES$ID, '_', CASES$Date)
CASES = CASES[!is.na(CASES$Cases),]
### Remove Duplicates from Multiple Sources, Prioritizing Better Sources
CASES$Source = factor(CASES$Source, levels = c('JHU', 'JRC', 'NYT',"SES", "RKI", "DPC" ,"CTP", "NYC" ))
CASES = CASES[order(CASES$ID, CASES$Date, CASES$Source),]
CASES = CASES[!duplicated(CASES$ID_Date),]

DEATHS = Hopkins_Cases[Hopkins_Cases$Type == 'Deaths',]
DEATHS$ID_Date = paste0(DEATHS$ID, '_', DEATHS$Date)
DEATHS = DEATHS[!is.na(DEATHS$Cases),]
DEATHS$Source = factor(DEATHS$Source, levels = c('JHU', 'JRC', 'NYT',"SES", "RKI", "DPC" ,"CTP", "NYC" ))
DEATHS = DEATHS[order(DEATHS$ID, DEATHS$Date, DEATHS$Source),]
DEATHS = DEATHS[!duplicated(DEATHS$ID_Date),]

CLINICAL = CASES[,c('ID', 'Date', 'Cases', 'Cases_New')]
CLINICAL = merge(CLINICAL, data.frame(DEATHS[,c('ID', 'Date')], Deaths = DEATHS$Cases, Deaths_New = DEATHS$Cases_New), all.x = T, all.y = T)
#Note: Cases_New and Deaths_New may be negative if the number of cumulative cases decreases between days. Raw data value.
CLINICAL = merge(CLINICAL, Hopkins_LUT, by = c('ID'),all.x = T, all.y = T)

### Additional Cleaning
CLINICAL$ID_Date = paste0(CLINICAL$ID, '_', CLINICAL$Date)
CLINICAL$days_into_2021 = as.numeric(as.Date(CLINICAL$Date, format = '%Y-%m-%d') - as.Date("2021-01-01", format = '%Y-%m-%d'))
CLINICAL = CLINICAL[!is.na(CLINICAL$Date),]

### FINAL ERROR CHECKING ###
DUPLICATES = CLINICAL[duplicated(CLINICAL$ID_Date),] #should not have duplicates
CLINICAL[CLINICAL$Admin0 =='',] #should be empty
CLINICAL[is.na(CLINICAL$Admin0),] #should be empty

CLINICAL = CLINICAL[order(CLINICAL$Admin012, CLINICAL$days_into_2021),]




#########################
### Read in OWID Data ###
#########################
#https://github.com/owid/covid-19-data/tree/master/public/data

OWID = read.csv(paste0(datapath,'owid-covid-data.csv'))
OWID = OWID[!(OWID$location %in% c('Africa', 'High income', 'European Union', 'Upper middle income', 'Asia', 'Oceania', 'World', 
                                   'Europe', 'International', 'Lower middle income', 'Low income', 'North America', 'South America')),]

COUNTRIES = CLINICAL[!duplicated(CLINICAL$ID),c('Admin0', 'Admin1','Admin2','ISO1_3N','ISO1_3C','ISO1_2C' ,'ISO2','ISO2_UID','FIPS','NUTS')]
COUNTRIES = COUNTRIES[COUNTRIES$Admin0 !='',]
COUNTRIES = COUNTRIES[COUNTRIES$Admin2 == '',]

### Merge by ISO1_3C
OWID$Admin0 = ifelse(OWID$iso_code %in% unique(COUNTRIES$ISO1_3C[COUNTRIES$Admin1 == '']), OWID$location, rep(NA,length(OWID$location)))
OWID$Admin1 = ifelse(!is.na(OWID$Admin0), '', NA)

### Remove locations not in JHU
OWID = OWID[!(OWID$location %in% c('Namibia','Macao','Nauru','Niue','Pitcairn','Tokelau','Turkmenistan','Tuvalu','Vatican','Northern Cyprus')),]

### Hard code some locations
OWID$Admin0[OWID$location == 'Bahamas'] = 'The Bahamas'
OWID$Admin1[OWID$location == 'Bahamas'] = ''
OWID$Admin0[OWID$location == 'Bonaire Sint Eustatius and Saba'] = 'Netherlands'
OWID$Admin1[OWID$location == 'Bonaire Sint Eustatius and Saba'] = 'Bonaire, Sint Eustatius and Saba'
OWID$Admin0[OWID$location == 'British Virgin Islands'] = 'United Kingdom'
OWID$Admin1[OWID$location == 'British Virgin Islands'] = 'Virgin Islands'
OWID$Admin0[OWID$location == 'Congo'] = 'Congo (Brazzaville)'
OWID$Admin1[OWID$location == 'Congo'] = ''
OWID$Admin0[OWID$location == 'Democratic Republic of Congo'] = 'Congo (Kinshasa)'
OWID$Admin1[OWID$location == 'Democratic Republic of Congo'] = ''
OWID$Admin0[OWID$location == 'Faeroe Islands'] = 'Denmark'
OWID$Admin1[OWID$location == 'Faeroe Islands'] = 'Faroe'
OWID$Admin0[OWID$location == 'Gambia'] = 'The Gambia'
OWID$Admin1[OWID$location == 'Gambia'] = ''
OWID$Admin0[OWID$location == 'Micronesia (country)'] = 'Federated States of Micronesia'
OWID$Admin1[OWID$location == 'Micronesia (country)'] = ''
OWID$Admin0[OWID$location == 'Northern Cyprus'] = 'Cyprus'
OWID$Admin1[OWID$location == 'Northern Cyprus'] = ''
OWID$Admin0[OWID$location == 'Saint Helena'] = 'United Kingdom'
OWID$Admin1[OWID$location == 'Saint Helena'] = 'St. Helens'
OWID$Admin0[OWID$location == 'Sint Maarten (Dutch part)'] = 'Netherlands'
OWID$Admin1[OWID$location == 'Sint Maarten (Dutch part)'] = 'Sint Maarten'
OWID$Admin0[OWID$location == 'Timor'] = 'East Timor'
OWID$Admin1[OWID$location == 'Timor'] = ''
OWID$Admin0[OWID$location == 'United States Virgin Islands'] = 'United States'
OWID$Admin1[OWID$location == 'United States Virgin Islands'] = 'Virgin Islands'
OWID$Admin0[OWID$location == 'Kosovo'] = 'Kosovo'
OWID$Admin1[OWID$location == 'Kosovo'] = ''
OWID$Admin0[OWID$location == 'Saint Martin (French part)'] = 'France'
OWID$Admin1[OWID$location == 'Saint Martin (French part)'] = 'Saint Martin'
OWID$Admin0[OWID$location == 'United States'] = 'United States' #1091 entries just as United States
OWID$Admin1[OWID$location == 'United States'] = ''

### Some countries are actually regions in JHU
OWID = merge(OWID, data.frame(location = COUNTRIES$Admin1, country_jhu = COUNTRIES$Admin0), by = 'location', all.x = T, all.y = F)
OWID$Admin0[!is.na(OWID$country_jhu)] = OWID$country_jhu[!is.na(OWID$country_jhu)]
OWID$Admin1[!is.na(OWID$country_jhu)] = OWID$location[!is.na(OWID$country_jhu)]

### Error Checking (should be empty!)
table(OWID[is.na(OWID$Admin0), 'location'])

### Define other variables
OWID$Admin01 = paste0(OWID$Admin0,'_',OWID$Admin1)
OWID$days_into_2021 = as.numeric( as.Date(OWID$date, format = '%Y-%m-%d') - as.Date("2021-01-01", format = '%Y-%m-%d'))
OWID[duplicated(paste0(OWID$Admin01,'_',OWID$days_into_2021)),]

OWID = subset(OWID, select = -c(country_jhu, location))




##################
### Merge Data ###
##################

### JHU
CLINICAL = CLINICAL[CLINICAL$Admin2 == '',]
COUNTRIES = CLINICAL[!duplicated(CLINICAL$Admin01),c('Admin0', 'Admin1', 'Admin01')] #should be no duplicates


### OWID
OWID$Date = as.Date(OWID$date, format = '%Y-%m-%d')
names(OWID)[!(names(OWID) %in% c('iso_code', 'date', 'Admin0', 'Admin1','Admin01', 'days_into_2021','Date'))] = paste0('owid_',names(OWID)[!(names(OWID) %in% c('iso_code', 'date', 'Admin0', 'Admin1','Admin01', 'days_into_2021','Date'))])
OWID = subset(OWID, select =-c(iso_code, date))
CLINICAL = merge(CLINICAL, data.frame(OWID), by = c('Admin0', 'Admin1','Admin01', 'days_into_2021','Date'), all.x = T, all.y = T)

### Make sure numeric
CLINICAL$Cases_New = as.numeric(as.character(CLINICAL$Cases_New))
CLINICAL$Deaths_New = as.numeric(as.character(CLINICAL$Deaths_New))

### Define country vs region indicator
CLINICAL$is_country = as.numeric(CLINICAL$Admin1=='')

### Reorder rows
CLINICAL = CLINICAL[order(CLINICAL$Admin01, CLINICAL$days_into_2021),]
CLINICAL = CLINICAL[CLINICAL$Admin0 != '',]

### Error checking: should be empty
table(CLINICAL[duplicated(paste0(CLINICAL$Admin01, '_', CLINICAL$days_into_2021)),'Admin0'])

#dim(CLINICAL[as.character(CLINICAL$Date) == '',])



######################
### Organize Lists ###
######################
library(data.table)

df = CLINICAL

## get all the unique geographies
unqgeo <- sort(unique(df$Admin01))

## add week and day of week
df$year <- year(df$Date)
df$epiweek <- epiweek(df$Date)

## subset to columns we care about
dfsmall <- subset(df, select=c("Date","year","epiweek","Admin01","Cases_New","Deaths_New","owid_weekly_hosp_admissions"))
rm(df)

## melt it
dfsmallmelt <- reshape2::melt(dfsmall, id.vars = c("Date","year","epiweek","Admin01"))

## remove rows where Date = NA
dfsmallmelt <- dfsmallmelt[!is.na(dfsmallmelt$year),]

## turn all value == NAs to 0
dfsmallmelt$value[is.na(dfsmallmelt$value)] = 0

## turn all negative values to 0
dfsmallmelt$value[dfsmallmelt$value < 0] = 0

## make the ts_disease name
dfsmallmelt$disease <- paste0("covid_",gsub("num_","",dfsmallmelt$variable))
dfsmallmelt[dfsmallmelt$variable %in% c("Cases_New","Deaths_New","owid_weekly_hosp_admissions"),]$disease <- "covid"
dfsmallmelt$variable <- as.character(dfsmallmelt$variable)



dfsmallmelt = dfsmallmelt[dfsmallmelt$Date <= max(dfsmallmelt$Date[dfsmallmelt$value > 0]),]
dfsmallmelt = dfsmallmelt[dfsmallmelt$Date >= min(dfsmallmelt$Date[dfsmallmelt$value > 0]),]
aggregate(Date~variable, FUN = min, data = dfsmallmelt[dfsmallmelt$value > 0,])
aggregate(Date~variable, FUN = max, data = dfsmallmelt[dfsmallmelt$value > 0,])


## get unique Admin01/variable
unqgroups <- dfsmallmelt[!duplicated(paste0(dfsmallmelt$Admin01, dfsmallmelt$variable)), c('Admin01', 'variable')]


## loop over all rows of unqgroups
global_covid_training_data_daily <- list()
global_covid_training_data_weekly <- list()

counts_daily <- 0
counts_weekly <- 0
for(i in 1:nrow(unqgroups)){

  ### Daily Data ###
  tempdf <- dfsmallmelt[dfsmallmelt$Admin01 == unqgroups$Admin01[i] & dfsmallmelt$variable == unqgroups$variable[i],]
  tempdf <- tempdf[order(tempdf$Date),]
  tempdf <- tempdf[-1,]
  tempdf <- tempdf[!is.na(tempdf$epiweek),]
  
  ## turn NAs into 0s
  if(sum(is.na(tempdf$value)) > 0){
    tempdf[is.na(tempdf$value),]$value <- 0
  }
  
  ## make list item
  if(nrow(tempdf)>= 10 & !is.na(var(tempdf$value, na.rm=T)) & var(tempdf$value, na.rm=T) > 0 & grepl('United States',tempdf$Admin01[1]) & tempdf$variable[1] != 'owid_weekly_hosp_admissions'){
    counts_daily <- counts_daily + 1
    templist <- list(ts = pmax(0,tempdf$value),
                     ts_dates = tempdf$Date,
                     ts_dates_actual = tempdf$Date, 
                     ts_disease = tempdf$disease[1],
                     ts_measurement_type = ifelse(tempdf$variable[1] == "Deaths_New","deaths",ifelse(tempdf$variable[1] == 'Cases_New',"incidence","hospitalizations")),
                     ts_geography = trimws(gsub("_"," ",tempdf$Admin01[1])),
                     ts_first_time = min(tempdf$Date, na.rm=T),
                     ts_last_time = max(tempdf$Date, na.rm=T),
                     ts_time_cadence = "daily",
                     ts_scale = "counts")
    global_covid_training_data_daily[[counts_daily]] <- templist
  }
  
  ## Weekly Data
  tempdfw <- data.table(tempdf)[,.(ndays = length(value),
                       Date = min(Date),
                       variable = variable[1],
                       disease = disease[1],
                       Admin01 = Admin01[1],
                       value = sum(value, na.rm=T)),by=c("epiweek","year")]
  tempdfw <- tempdfw[!is.na(tempdfw$epiweek) & tempdfw$ndays == 7,]
  
  ## weekly cases
  if(nrow(tempdfw)>= 10 & !is.na(var(tempdfw$value, na.rm=T)) & var(tempdfw$value, na.rm=T) > 0){
    counts_weekly <- counts_weekly + 1
    templist <- list(ts = pmax(0,tempdfw$value),
                     ts_dates = tempdfw$Date,
                     ts_dates_actual = tempdfw$Date, 
                     ts_disease = tempdfw$disease[1],
                     ts_measurement_type = ifelse(tempdf$variable[1] == "Deaths_New","deaths",ifelse(tempdf$variable[1] == 'Cases_New',"incidence","hospitalizations")),
                     ts_geography = trimws(gsub("_"," ",tempdfw$Admin01[1])),
                     ts_first_time = min(tempdfw$Date, na.rm=T),
                     ts_last_time = max(tempdfw$Date, na.rm=T),
                     ts_time_cadence = "weekly",
                     ts_scale = "counts")
    global_covid_training_data_weekly[[counts_weekly]] <- templist
  }
}

## concatenate all the lists
global_covid_training_data <- c(global_covid_training_data_daily,
                                global_covid_training_data_weekly)
length(global_covid_training_data)

## save the us ili training data
saveRDS(global_covid_training_data, paste0(savepath,"COVID_jhuowid.RDS"))



