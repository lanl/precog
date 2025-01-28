## Dave Osthus
## 2-15-24
## Make US project tycho level 1 data

###############################################################
## NOTE: I downloaded the raw data read in within this script from here: https://www.tycho.pitt.edu/data/#datasets
## This is the entire "level 1" dataset, including 8 diseases by US state and some cities (total 172 locations)


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
setwd("~/GitLab/epifforma")
datapath <- "./raw_data/project_tycho/raw_data/" # Murph issue!
savepath <- "./raw_data/project_tycho/output/"
plotpath <- "./raw_data/project_tycho/figs/"

###############################################################
## set proxies
Sys.setenv('https_proxy'='http://proxyout.lanl.gov:8080')## not sure why needed, but see here: https://github.com/curl/curl/issues/1015

###############################################################
#### Read in Raw COVID Data
df <- data.table::fread(paste0(datapath,"ProjectTycho_Level1_v1.0.0.csv"))

## get all the unique geographies
unqgeo <- sort(unique(df$loc))


## subset to columns we care about
dfsmall <- subset(df, select=c("loc","disease","cases","epi_week"))
rm(df)

#get dataset with zeros
df_long = expand.grid(disease = unique(dfsmall$disease),
                      loc = unique(dfsmall$loc),
                      epi_week = unique(dfsmall$epi_week))
df_long = merge(df_long, dfsmall, by = c("epi_week","loc","disease"),
                all.x = T, all.y = T)
df_long$cases[is.na(df_long$cases)]=0
df_long$cases[df_long$cases == '\\N']=0
df_long$cases = as.numeric(df_long$cases)
df_long$epi_week = as.numeric(as.character(df_long$epi_week))




### collapse data into 4 week cadence

df_long$year = as.numeric(substr(df_long$epi_week,1,4))
df_long$week_in_year = as.numeric(substr(df_long$epi_week,5,6))

### project tycho data doesn't allow for leap year week 53. Only 52 weeks per year. So, we can also ignore the leap weeks.  
df_long$week_index = 52*as.numeric(substr(df_long$epi_week,1,4))+as.numeric(substr(df_long$epi_week,5,6))
df_long$month_index = floor(df_long$week_index/4)


df_long = df_long %>% dplyr::group_by(loc, disease, month_index) %>% dplyr::mutate(epi_week_month = epi_week[which.max(week_index)],
                                                                                   cases_month = sum(cases))

df_long = df_long[!duplicated(paste0(df_long$loc, '_', df_long$disease, '_' ,df_long$month_index)),]
df_long$epi_week = df_long$epi_week_month
df_long$cases = df_long$cases_month




library(dplyr)

K = 5
df_long = df_long %>% dplyr::group_by(loc,disease) %>% dplyr::mutate(min_nonzero = min(week_index[cases > K],na.rm=T)) #actually min greater than K
df_long = df_long %>% dplyr::group_by(loc,disease) %>% dplyr::mutate(max_nonzero = max(week_index[cases > K],na.rm=T)) #actually min greater than K
df_long = df_long %>% dplyr::group_by(loc,disease) %>% dplyr::mutate(total_cases = sum(cases,na.rm=T))
df_long = df_long[!is.infinite(df_long$min_nonzero),]
df_long = df_long[df_long$week_index >= (df_long$min_nonzero - 10),]
df_long = df_long[df_long$week_index <= (df_long$max_nonzero + 10),]
df_long$year <- as.numeric(substr(df_long$epi_week,1,4))
df_long$epiweek <- as.numeric(substr(df_long$epi_week,5,6))
df_long$variable = 'monthly cases'

df_long$disease = gsub(' ','',df_long$disease)

## look at number of cases per locaiton
# DAT = df_long[!duplicated(paste0(df_long$loc,'_',df_long$disease)), c('loc','disease','total_cases')]
# ggplot(DAT)+
#   geom_boxplot(aes(x=disease, y = log10(total_cases+1)))+
#   geom_hline(yintercept = log10(5000+1),color = 'red')

##subset to at least 5000 cases per location
df_long = df_long[df_long$total_cases >= 5000,]

## get unique location/disease
unqgroups = df_long[!duplicated(paste0(df_long$loc,'_',df_long$disease)), c('loc','disease')]


table(unqgroups$disease)
AGG = aggregate(month_index~disease, FUN = function(x){c(min(x,na.rm=T),max(x,na.rm=T))}, data = df_long[!is.na(df_long$cases) & df_long$cases > 0,])
for(i in 1:length(AGG[,1])){
  print(data.frame(  df_long[df_long$month_index ==AGG[i,2][1], c('epi_week')][1,], vector = AGG[i,1]))
}

for(i in 1:length(AGG[,1])){
  print(data.frame(  df_long[df_long$month_index ==AGG[i,2][2], c('epi_week')][1,], vector = AGG[i,1]))
}
#################################
### Visualize the Time Series ###
#################################
DISEASES = unique(df_long$disease)
for(i in 1:length(DISEASES)){
  SUBSET = df_long[df_long$disease == DISEASES[i],]
  SUBSET = SUBSET %>% dplyr::group_by(loc) %>% dplyr::mutate(cases_scaled = cases/max(cases,na.rm=T))
  p1 = ggplot(SUBSET)+
    geom_tile(aes(x=month_index,y=loc,fill=cases_scaled,color=cases_scaled))+
    theme_classic()+
    scale_color_viridis_c('Scaled Cases', na.value = 'lightgray')+
    scale_fill_viridis_c('Scaled Cases', na.value = 'lightgray')+
    xlab('Week')+
    ylab('Location')+
    labs(title = 'Scaled Cases')+
    scale_x_continuous(expand=c(0,0))+
    theme(legend.position = 'top')
  p2 = ggplot(SUBSET)+
    geom_tile(aes(x=week_index,y=loc,fill=cut(cases, breaks = c(0,1,2,5,10,100,10^10)),color=cut(cases, breaks = c(0,1,2,5,10,100,10^10))))+
    theme_classic()+
    scale_color_viridis_d('Cases', na.value = 'lightgray')+
    scale_fill_viridis_d('Cases', na.value = 'lightgray')+
    xlab('Week')+
    ylab('Location')+
    labs(title = 'Cases')+
    scale_x_continuous(expand=c(0,0))+
    theme(legend.position = 'top')
  p0= cowplot::plot_grid(p1,p2,ncol=2)
  ggsave(p0,filename = paste0(plotpath,'/',DISEASES[i],'_cases.png'), width = 10, height = 8)
  
  AGG = aggregate(cases~loc, FUN = sum, na.rm=T, data = SUBSET)
}




# 
# dim(AGG)
# dim(AGG[AGG$cases>5000,])
# 



# 
# 
# ggplot(df_long[df_long$disease == 'MEASLES',])+
#   geom_point(aes(x=week_index,y=cases, group = factor(loc)))



unqgroups$seasonal = T
unqgroups$seasonal[unqgroups$disease == 'HEPATITISA'] = F #not seasonal


## loop over all rows of unqgroups
us_tycho_data <- list()
counts_weekly <- 0

for(i in 1:nrow(unqgroups)){
  
  ################################################
  #### make daily 
  
  ## subset to geography
  tempdf <- df_long[df_long$loc == unqgroups$loc[i] & df_long$disease == unqgroups$disease[i],]
  tempdf <- tempdf[order(tempdf$year, tempdf$epiweek),]
  
  ## make list item
  if(nrow(tempdf)>= 10 & !is.na(var(tempdf$cases, na.rm=T)) & var(tempdf$cases, na.rm=T) > 0){
    counts_weekly <- counts_weekly + 1
    templist <- list(ts = pmax(0,tempdf$cases),
                     ts_dates = tempdf$epi_week,
                     ts_exogenous = NULL,
                     ts_real_data = T,
                     ts_isolated_strain = F,
                     ts_multiwave = T,
                     ts_disease = tempdf$disease[1],
                     ts_measurement_type = 'incident_cases',
                     ts_geography = trimws(gsub("_"," ",tempdf$loc[1])),
                     ts_first_time = min(tempdf$epi_week, na.rm=T),
                     ts_last_time = max(tempdf$epi_week, na.rm=T),
                     ts_time_cadence = "monthly",
                     ts_scale = "counts",
                     ts_exogenous_scale = NA,
                     ts_seasonal = unqgroups$seasonal[i])
    us_tycho_data[[counts_weekly]] <- templist
  }
}

## save the us ili training data
saveRDS(us_tycho_data, paste0(savepath,"us_project_tycho.RDS"))








# ## save all the lists as real ILI data time series
# tempid <- sample(1:length(global_covid_training_data),1)
# tempid
# ttt <- global_covid_training_data[[tempid]]
# qplot(ttt$ts_dates, ttt$ts, geom=c("point","line"))+
#   ggtitle(paste0(ttt$description$ts_disease,"     ",ttt$description$ts_isolated_strain,"       ",ttt$description$ts_geography,"       ",ttt$description$ts_time_cadence,"     ",ttt$description$ts_measurement_type))





