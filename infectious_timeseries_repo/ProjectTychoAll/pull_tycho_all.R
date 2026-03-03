#######################
#######################
### Pull Tycho Data ###
#######################
#######################

library(utils)
library(stringr)
library(httr)
library(data.table)


#################
### Functions ###
#################

access_token = ### USER WANTING TO PULL SOURCE DATA NEEDS TO GET THEIR OWN ACCESS TOKEN



# Use the authorization token to access the secured resource
# example: get a listing of events
# query_params_list is a list with name value pairs of query parameters to provide to the api end point
get_resource <- function(access_token = NULL, path =
                           NULL, query_params_list = NULL)
{
  headers = c(
    `Authorization` = paste0("Bearer ", access_token)#,
  )
  query <- query_params_list
  action_path_encoded <- httr::modify_url(path, query = query)
  res <-  GET(url = action_path_encoded, 
             httr::timeout(29), VERBOSE = TRUE)
  result_list <- content(res, type = "text/csv")
  return(result_list)
}



#####################
### Retrieve Data ###
#####################

### List of Pathogens
condition_list <- data.frame(get_resource( access_token = access_token, path = "https://www.tycho.pitt.edu/api/condition?", query_params_list = list(apikey = access_token)))
country_list <- data.frame(get_resource( access_token = access_token, path = "https://www.tycho.pitt.edu/api/country?", query_params_list = list(apikey = access_token)))


library(tidyr)
library(jsonlite)

read_tycho = function(condition, country){
  event_list_long = NULL
  STOP_FLAG = F
  OFFSET = 0
  i = 1
  while(!STOP_FLAG){
    event_list <- get_resource(access_token = access_token, path = "https://www.tycho.pitt.edu/api/query?", 
                               query_params_list = list(apikey = access_token, ConditionName = condition, CountryName = country, limit = 5000, offset = OFFSET))
    if(length(event_list) > 1){
      event_list = as.data.frame(event_list)
      if(!is.null(event_list)){
        event_list_long = dplyr::bind_rows(event_list_long, event_list)
        if(nrow(event_list) < 5000){
          STOP_FLAG = T
        }
      }else{
        STOP_FLAG = T
      }
      i = i + 1
      OFFSET = OFFSET + 5000
    }else{
      STOP_FLAG = T
    }
  }
  return(event_list_long)
}


for(k in 1:length(condition_list$ConditionName)){
  data_long = NULL
  for(u in 1:length(country_list$CountryName)){
    data_out = read_tycho(condition = as.character(unlist(condition_list$ConditionName[k])),
                          country = as.character(unlist(country_list$CountryName[u])))
    if(!is.null(data_out)){
      data_long = rbind(data_long, data_out)
    }
  }
  write.csv(data_long, file = paste0("./raw_data/",as.character(unlist(condition_list$ConditionName[k])),".csv"), quote = T, row.names = F)
  print(k)
}


