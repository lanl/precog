#######################
#######################
### Clean AIDO Data ###
#######################
#######################
 #https://aido.bsvgateway.org/api/outbreaks/?format=api


library(MMWRweek)
library(lubridate)
library(ggplot2)
library(plyr)
library(lubridate)
theme_set(theme_bw())
datapath <- "./raw_data/" 
savepath <- "../_Organized_Lists/"
plotpath <- "./figs/"

####################
### Read in Data ###
####################

FILES = list.files(datapath)
df = NULL
for(i in 1:length(FILES)){
  dat = try(read.csv(paste0(datapath, FILES[i])))
  if(class(dat)[1] != 'try-error'){
    df = rbind(df, dat)
  }
}
df$start_date = as.Date(df$start, tz = "UTC")
df$end_date = as.Date(df$end, tz = "UTC")

disease_list = read.csv(paste0("/Users/lbeesley/Desktop/New_Data_Streams/AIDO/disease_list.csv"))
df$id = as.numeric(do.call('rbind',strsplit(df$disease, split = "/"))[,6])
df = merge(df, disease_list, by = 'id', all.x = T, all.y = F)


##########################
### Subset Time Series ###
##########################

AGG = aggregate(value~name+location_fields.full_name, FUN = length, data = df)

### Require at least 50 observations per time series
AGG_SUB = AGG[AGG$value >= 50,]
df_sub = df[paste0(df$location_fields.full_name, '_', df$name) %in% paste0(AGG_SUB$location_fields.full_name, '_', AGG_SUB$name),]

### This dataset contains several diseases, but we ended up focusing on just 4. Other users can pull out more data in the future. 
### Subset by disease and location
DISEASES = c('Zika', 'Cholera', 'Porcine Epidemic Diarrhea Virus (PEDV)', 'Malaria')
df_sub = df[paste0(df$location_fields.full_name, '_', df$name) %in% paste0(AGG_SUB$location_fields.full_name, '_', AGG_SUB$name) & df$name %in% DISEASES,]
df_sub$key = paste0(df_sub$location_fields.name, '_', df_sub$name)
df_sub = df_sub[!(df_sub$key %in% c('Manitoba_Porcine Epidemic Diarrhea Virus (PEDV)','Ontario_Porcine Epidemic Diarrhea Virus (PEDV)')),]

df_sub$lag = difftime(df_sub$end_date,df_sub$start_date)
plot(df_sub$lag)

table(df_sub$lag, df_sub$name)
# Cholera is weekly, Malaria is weekly or daily depending on location
# PEDV is weekly, Zika is daily



######################
### Generate Lists ###
######################


df_sub$name[df_sub$name == 'Porcine Epidemic Diarrhea Virus (PEDV)'] = 'Porcine_Epidemic_Diarrhea'

DISEASES = unique(df_sub$name)

for(d in 1:length(DISEASES)){
  SUB = df_sub[df_sub$name == DISEASES[d],]
  
  data_list = list()
  num = 0
  LOCATIONS = unique(SUB$location_fields.full_name)
  for(l in 1:length(LOCATIONS)){
    SUBSUB = SUB[SUB$location_fields.full_name == LOCATIONS[l],]
    SUBSUB = SUBSUB[order(SUBSUB$start_date),]
    SUBSUB$value[is.na(SUBSUB$value)]=0
    if(nrow(SUBSUB)>= 10 & !is.na(var(SUBSUB$value, na.rm=T)) & var(SUBSUB$value, na.rm=T) > 0){
      num <- num + 1
      templist <- list(ts = pmax(0,SUBSUB$value),
                       ts_dates = as.character(SUBSUB$start_date),
                       ts_dates_actual = as.character(SUBSUB$start_date),
                       ts_disease = DISEASES[d],
                       ts_measurement_type = "incidence",
                       ts_geography = trimws(LOCATIONS[l]),
                       ts_first_time = min(SUBSUB$start_date, na.rm=T),
                       ts_last_time = max(SUBSUB$start_date, na.rm=T),
                       ts_time_cadence = ifelse(unique(SUBSUB$lag)==7,"weekly","daily"),
                       ts_scale = "counts")
      data_list[[num]] <- templist
    }
  }
  saveRDS(data_list, paste0(savepath,DISEASES[d],"_aido.RDS"))
  
}





