## Show MutAntiGen Stochasticity 

# load libraries
library(data.table)
library(ggplot2)
theme_set(theme_bw())

# define path
rawdatapath <- paste0(here::here("synthetic_and_genetic_forecasting", "data_raw", "mutantigen"), "/") 
figpath <- paste0(here::here("synthetic_and_genetic_forecasting", "figs"), "/") 

## read in these .timeseries ids: 1 153 305 457 609 913 1065 1217 1369 1673
ts_id <- c(1, 153, 305, 457, 609, 913, 1065, 1217, 1369, 1673)

## read in .timeseries files
ts_df <- NULL
for(i in 1:length(ts_id)){
  tempdf <- fread(paste0(rawdatapath,"out_",ts_id[i],".timeseries"))
  ts_df <- rbind(ts_df, data.table(id = as.factor(ts_id[i]),
                                   newid = i,
                                   date = tempdf$date,
                                   totalCases = tempdf$totalCases))
}


# save it
setEPS()
postscript(paste0(figpath,"mutantigen_stochasticity.eps"),width = 12, height = 4)
  qplot(date, totalCases, data = subset(ts_df,date > .1*max(date)), group = id, geom="line",size=I(.5),show.legend=F)+
    facet_wrap(~newid, nrow = 2, scales="free")+
    xlab("")+
    ylab("Total Cases")
dev.off()





