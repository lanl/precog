# Script to get the model scores.
# This takes a while, so usually we just run it once and 
# save the output to a csv.
## Author: GC Gibson and AC Murph

library(covidHubUtils)
library(dplyr)
library(spatstat)
library(BASS)
library(covidcast)

######
# Get the Model Scores.
# truth_data                              <- load_truth(
#   truth_source = "JHU",
#   target_variable = "inc case")
# 
# 
# truth_weekly                            <- load_truth(truth_source = "JHU",hub = "US",
#                                                       target_variable = "inc case", 
#                                                       as_of= "2023-03-04",locations =state.name,data_location="covidData")

load(file = 'truth_data.rds')
load(file = 'truth_weekly.rds')

model_list                              <- get_all_models()


inc_case_targets                        <- paste(1:4, "wk ahead inc case")
forecasts_case_model_res                <- list()

model_list_i                            <- 1
problem_list <- c("CDDEP-ABM")
for (model_idx in 1:length(model_list)){
  # For model_idx == 110, (mech bayes), this gives nothing...check in with Casey?
  print (model_idx/length(model_list))
  
  if(model_list[model_idx]%in%problem_list){
    next
  }else {
    forecasts_case                        <- load_forecasts(
      models = model_list[model_idx],
      date_window_size = 6,
      types = c("point","quantile"),
      targets = inc_case_targets,
      source = "local_hub_repo",
      hub_repo_path = '~/GitHub/covid19-forecast-hub',
      # verbose = FALSE,
      horizon=1:4,
      locations = state.name,
      # as_of = NULL,
      hub = c("US"))
  }
  forecasts_case_model_res[[model_idx]] <- forecasts_case
}

forecasts_case                          <- do.call(rbind,forecasts_case_model_res)

forecasts_case                          <- forecasts_case[nchar(forecasts_case$location) <= 2, ]
scores                                  <- score_forecasts(
  forecasts = forecasts_case,
  metrics = c("abs_error", "wis"),
  return_format = "wide",
  truth = truth_data
)


mean(scores[scores$location == "06",]$abs_error)
write.csv(x=scores,file = "data/scores_tot_w_deaths.csv")

