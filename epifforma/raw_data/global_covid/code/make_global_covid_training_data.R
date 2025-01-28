
###############################################################
## load libraries
library(MMWRweek)
library(lubridate)
library(ggplot2)
library(plyr)
library(lubridate)
theme_set(theme_bw())

###############################################################
## set path

library(this.path)
setwd(this.path::here()) 
#data path not provided, raw data is big
savepath <- "./output/"
plotpath <- "./figs/"

###############################################################
#### Read in Raw COVID Data
df <- data.table::fread(paste0(datapath,"Merged_Regional.csv"))

## get all the unique geographies
unqgeo <- sort(unique(df$Admin01))

## add week and day of week
df$year <- year(df$Date)
df$epiweek <- epiweek(df$Date)

## subset to columns we care about
dfsmall <- subset(df, select=c("Date","year","epiweek","Admin01","Cases_New","Deaths_New"))
rm(df)

## melt it
dfsmallmelt <- data.table(melt(dfsmall, id.vars = c("Date","year","epiweek","Admin01")))

## remove rows where Date = NA
dfsmallmelt <- dfsmallmelt[!is.na(dfsmallmelt$year),]

## turn all value == NAs to 0
dfsmallmelt[is.na(dfsmallmelt$value),]$value <- 0

## turn all negative values to 0
dfsmallmelt[value < 0,]$value <- 0

## make the ts_disease name
dfsmallmelt$disease <- paste0("covid_",gsub("num_","",dfsmallmelt$variable))
dfsmallmelt[dfsmallmelt$variable %in% c("Cases_New","Deaths_New"),]$disease <- "covid"
dfsmallmelt$variable <- as.character(dfsmallmelt$variable)


## subset to 2020-01-02 to 03/31/2023 
dfsmallmelt = dfsmallmelt[dfsmallmelt$Date <= max(dfsmallmelt$Date[dfsmallmelt$value > 0]),]
dfsmallmelt = dfsmallmelt[dfsmallmelt$Date >= min(dfsmallmelt$Date[dfsmallmelt$value > 0]),]

A = dfsmallmelt[dfsmallmelt$Date == min(dfsmallmelt$Date[dfsmallmelt$value > 0]),]
B = A[A$value > 0,]


## get unique Admin01/variable
unqgroups <- dfsmallmelt[,.(n = length(value)), by = c("Admin01","variable")]


### plot country-level data
dfsmallmelt$week_index = as.numeric(unlist(dfsmallmelt$epiweek))+52*as.numeric(unlist(dfsmallmelt$year))

SUBSET = dfsmallmelt[substr(dfsmallmelt$Admin01,nchar(dfsmallmelt$Admin01),nchar(dfsmallmelt$Admin01)) == '_',]
SUBSET = SUBSET[SUBSET$variable == 'Cases_New',]
SUBSET = SUBSET[SUBSET$epiweek != 53,]
SUBSET = SUBSET %>% dplyr::group_by(Admin01) %>% dplyr::mutate(cases_scaled = value/max(value,na.rm=T))
p1 = ggplot(SUBSET)+
  geom_tile(aes(x=week_index,y=Admin01,fill=cases_scaled,color=cases_scaled))+
  theme_classic()+
  scale_color_viridis_c('Scaled Cases', na.value = 'lightgray')+
  scale_fill_viridis_c('Scaled Cases', na.value = 'lightgray')+
  xlab('Week')+
  ylab('Location')+
  labs(title = 'Scaled Cases')+
  scale_x_continuous(expand=c(0,0))+
  theme(legend.position = 'top')
p2 = ggplot(SUBSET)+
  geom_tile(aes(x=week_index,y=Admin01,fill=cut(value, breaks = c(0,10,10000,100000,10^10)),color=cut(value, breaks = c(0,10,10000,100000,10^10))))+
  theme_classic()+
  scale_color_viridis_d('Cases', na.value = 'lightgray')+
  scale_fill_viridis_d('Cases', na.value = 'lightgray')+
  xlab('Week')+
  ylab('Location')+
  labs(title = 'Cases')+
  scale_x_continuous(expand=c(0,0))+
  theme(legend.position = 'top')
p0= cowplot::plot_grid(p1,p2,ncol=2)
ggsave(p0,filename = paste0(plotpath,'/global_cases.png'), width = 10, height = 8)



SUBSET = dfsmallmelt[grepl('United States', dfsmallmelt$Admin01),]
SUBSET = SUBSET[SUBSET$variable == 'Cases_New',]
SUBSET = SUBSET[SUBSET$epiweek != 53,]
SUBSET = SUBSET %>% dplyr::group_by(Admin01) %>% dplyr::mutate(cases_scaled = value/max(value,na.rm=T))
p1 = ggplot(SUBSET)+
  geom_tile(aes(x=week_index,y=Admin01,fill=cases_scaled,color=cases_scaled))+
  theme_classic()+
  scale_color_viridis_c('Scaled Cases', na.value = 'lightgray')+
  scale_fill_viridis_c('Scaled Cases', na.value = 'lightgray')+
  xlab('Week')+
  ylab('Location')+
  labs(title = 'Scaled Cases')+
  scale_x_continuous(expand=c(0,0))+
  theme(legend.position = 'top')
p2 = ggplot(SUBSET)+
  geom_tile(aes(x=week_index,y=Admin01,fill=cut(value, breaks = c(0,10,10000,100000,10^10)),color=cut(value, breaks = c(0,10,10000,100000,10^10))))+
  theme_classic()+
  scale_color_viridis_d('Cases', na.value = 'lightgray')+
  scale_fill_viridis_d('Cases', na.value = 'lightgray')+
  xlab('Week')+
  ylab('Location')+
  labs(title = 'Cases')+
  scale_x_continuous(expand=c(0,0))+
  theme(legend.position = 'top')
p0= cowplot::plot_grid(p1,p2,ncol=2)
ggsave(p0,filename = paste0(plotpath,'/us_cases.png'), width = 10, height = 8)



AGG = aggregate(value~Admin01, FUN = sum, na.rm=T, data = SUBSET)
AGG = AGG[order(AGG$value),]

A = SUBSET[SUBSET$Admin01 == 'United States_',]

## loop over all rows of unqgroups
global_covid_training_data_daily <- list()
global_covid_training_data_weekly <- list()

counts_daily <- 0
counts_weekly <- 0
for(i in 1:nrow(unqgroups)){
  ################################################
  ## where are we?
  print(counts_daily)
  
  ################################################
  #### make daily 
  
  # A = dfsmallmelt[grepl('United States',dfsmallmelt$Admin01) & dfsmallmelt$Admin01 != 'United States_',]
  # A[which.max(A$Date),]
  ## subset to geography
  tempdf <- dfsmallmelt[Admin01 == unqgroups$Admin01[i] & variable == unqgroups$variable[i],,]
  tempdf <- tempdf[order(tempdf$Date),]
  tempdf <- tempdf[-1,]
  tempdf <- tempdf[!is.na(tempdf$epiweek),]
  
  ## turn NAs into 0s
  if(sum(is.na(tempdf$value)) > 0){
    tempdf[is.na(tempdf$value),]$value <- 0
  }
  
  ## make list item
  if(nrow(tempdf)>= 10 & !is.na(var(tempdf$value, na.rm=T)) & var(tempdf$value, na.rm=T) > 0){
    counts_daily <- counts_daily + 1
    templist <- list(ts = pmax(0,tempdf$value),
                     ts_dates = tempdf$Date,
                     ts_exogenous = NULL,
                     ts_real_data = T,
                     ts_isolated_strain = ifelse(tempdf$variable[1] %in% c("Cases_New","Deaths_New"),F,T),
                     ts_multiwave = ifelse(tempdf$variable[1] %in% c("Cases_New","Deaths_New"),T,F),
                     ts_disease = tempdf$disease[1],
                     ts_measurement_type = ifelse(tempdf$variable[1] == "Deaths_New","incident_deaths","incident_cases"),
                     ts_geography = trimws(gsub("_"," ",tempdf$Admin01[1])),
                     ts_first_time = min(tempdf$Date, na.rm=T),
                     ts_last_time = max(tempdf$Date, na.rm=T),
                     ts_time_cadence = "daily",
                     ts_scale = "counts",
                     ts_exogenous_scale = NA)
    global_covid_training_data_daily[[counts_daily]] <- templist
  }
  
  ################################################
  #### make weekly 
  
  ## make weekly data
  tempdfw <- tempdf[,.(ndays = length(value),
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
                     ts_exogenous = NULL,
                     ts_real_data = T,
                     ts_isolated_strain = ifelse(tempdfw$variable[1] %in% c("Cases_New","Deaths_New"),F,T),
                     ts_multiwave = ifelse(tempdfw$variable[1] %in% c("Cases_New","Deaths_New"),T,F),
                     ts_disease = tempdfw$disease[1],
                     ts_measurement_type = ifelse(tempdfw$variable[1] == "Deaths_New","incident_deaths","incident_cases"),
                     ts_geography = trimws(gsub("_"," ",tempdfw$Admin01[1])),
                     ts_first_time = min(tempdfw$Date, na.rm=T),
                     ts_last_time = max(tempdfw$Date, na.rm=T),
                     ts_time_cadence = "weekly",
                     ts_scale = "counts",
                     ts_exogenous_scale = NA, 
                     ts_seasonal = F)
    global_covid_training_data_weekly[[counts_weekly]] <- templist
  }
}

## concatenate all the lists
global_covid_training_data <- c(global_covid_training_data_daily,
                                global_covid_training_data_weekly)
length(global_covid_training_data)

## save the us ili training data
saveRDS(global_covid_training_data, paste0(savepath,"global_covid.RDS"))


# ## save all the lists as real ILI data time series
# tempid <- sample(1:length(global_covid_training_data),1)
# tempid
# ttt <- global_covid_training_data[[tempid]]
# qplot(ttt$ts_dates, ttt$ts, geom=c("point","line"))+
#   ggtitle(paste0(ttt$description$ts_disease,"     ",ttt$description$ts_isolated_strain,"       ",ttt$description$ts_geography,"       ",ttt$description$ts_time_cadence,"     ",ttt$description$ts_measurement_type))





