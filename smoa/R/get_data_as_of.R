library(covidHubUtils)

truth <- load_truth(truth_source = "JHU",hub = "US", target_variable = "inc case",locations = "California")
unique_fcast_dates <- unique(truth$target_end_date)
tdat_list <- list()
for (dat_idx in 30:length(unique_fcast_dates)){
  print (dat_idx)
  truth_as_of <- load_truth(temporal_resolution = "daily",truth_source = "JHU",hub = "US", data_location="covidData", target_variable = "inc case" ,locations = "California",as_of=unique_fcast_dates[dat_idx] + 1)
  tdat_list[[dat_idx]] <- truth_as_of
}


for (dat_idx in 30:length(unique_fcast_dates)){
  tdat_list[[dat_idx]]$as_of <- unique_fcast_dates[dat_idx]
}

tdat_list_tot <- do.call(rbind,tdat_list)
tdat_list_tot$location

write.csv(tdat_list_tot,"tdat_list_tot_daily.csv")

  