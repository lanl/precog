# Script to get the model scores.
# This takes a while, so usually we just run it once and 
# save the output to a csv.
## Author: GC Gibson and AC Murph


######
# Get the Model Scores.
truth_data                              <- load_truth(
  truth_source = "JHU",
  target_variable = "inc case")


truth_weekly                            <- load_truth(truth_source = "JHU",hub = "US",
                                                      target_variable = "inc case", 
                                                      as_of= "2023-03-04",locations =state.name,data_location="covidData")


model_list                              <- get_all_models()
model_list <- c('COVIDhub-4_week_ensemble', 'COVIDhub-trained_ensemble', 'COVIDhub-baseline')

inc_case_targets                        <- paste(1:4, "wk ahead inc case")
forecasts_case_model_res                <- list()

model_list_i                            <- 1
problem_list <- c("UMass-MechBayes", "FAIR-NRAR", "FRBSF_Wilson-Econometric")
for (model_idx in 1:length(model_list)){
  # For model_idx == 110, (mech bayes), this gives nothing...check in with Casey?
  print (model_idx/length(model_list))
  
  if(model_list[model_idx]%in%problem_list){
    forecasts_case                        <- load_forecasts(
      models = model_list[model_idx],
      date_window_size = 6,
      types = c("point","quantile"),
      # targets = inc_case_targets,
      # source = "zoltar",
      # verbose = FALSE,
      locations = state.name,
      # as_of = NULL,
      hub = c("US"))
    forecasts_case <- forecasts_case[which((forecasts_case$horizon <=4)&(forecasts_case$target_variable=="inc cases")),]
  }else {
    forecasts_case                        <- load_forecasts(
      models = model_list[model_idx],
      date_window_size = 6,
      types = c("point","quantile"),
      targets = inc_case_targets,
      source = "zoltar",
      verbose = FALSE,
      locations = state.name,
      as_of = NULL,
      hub = c("US"))
  }
  forecasts_case_model_res[[model_idx]] <- forecasts_case
}

forecasts_case                          <- do.call(rbind,forecasts_case_model_res)

saveRDS(forecasts_case, file = 'data/forecasts.rds')

forecasts_case                          <- forecasts_case[nchar(forecasts_case$location) <= 2, ]


