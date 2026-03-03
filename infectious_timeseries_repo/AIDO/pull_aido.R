#####################
#####################
### Get AIDO Data ###
#####################
#####################


#https://aido.bsvgateway.org/api/outbreaks/?format=api


#########################
### Pull from Website ###
#########################
get_resource <- function(path =NULL, query_params_list = NULL){
  query <- query_params_list
  action_path_encoded <- httr::modify_url(path, query = query)
  res <-  GET(url = action_path_encoded, 
              httr::timeout(29), VERBOSE = TRUE)
  result_list <- content(res, encoding = "application/json")
  return(result_list)
}


library(tidyr)
library(jsonlite)
library(httr)

### List of Diseases
disease_list <- jsonlite::fromJSON(jsonlite::toJSON(get_resource(path = "https://aido.bsvgateway.org/api/diseases/?format=json")), flatten = T)
disease_list2 <- jsonlite::fromJSON(jsonlite::toJSON(get_resource(path = "https://aido.bsvgateway.org/api/diseases/?format=json", query_params_list = list(page = 3, format = 'json'))), flatten = T)
disease_list3 <- jsonlite::fromJSON(jsonlite::toJSON(get_resource(path = "https://aido.bsvgateway.org/api/diseases/?format=json", query_params_list = list(page = 4, format = 'json'))), flatten = T)
disease_list4 <- jsonlite::fromJSON(jsonlite::toJSON(get_resource(path = "https://aido.bsvgateway.org/api/diseases/?format=json", query_params_list = list(page = 5, format = 'json'))), flatten = T)

disease_list = dplyr::bind_rows(disease_list$results, disease_list2$results)
disease_list = dplyr::bind_rows(disease_list, disease_list3$results)
disease_list = dplyr::bind_rows(disease_list, disease_list4$results)
write.csv(data.frame(id = unlist(disease_list$id), name = unlist(disease_list$name)), file = paste0("./disease_list.csv"), quote = T, row.names = F)


read_aido = function(disease_num){
  event_list_long = NULL
  STOP_FLAG = F
  OFFSET = 0
  i = 1
  while(!STOP_FLAG){
    event_list <- get_resource(path = "https://aido.bsvgateway.org/api/outbreaks/?format=json", query_params_list = list(page = i+1, format = 'json', disease = disease_num))
    if(!('results' %in% names(event_list))){
      event_list <- get_resource(path = "https://aido.bsvgateway.org/api/outbreaks/?format=json", query_params_list = list(format = 'json', disease = disease_num))
    }
    if(length(event_list$results) == 0){
      print(paste0('Failed for ', disease_num)) 
      return(NULL)
    }else{
      df= jsonlite::fromJSON(jsonlite::toJSON(event_list$results), flatten = T)
      
      df_meta = subset(df, select = -c(time_series))
      df_meta[df_meta == 'NULL'] = NA
      for(j in 1:length(colnames(df_meta))){
        df_meta[,colnames(df_meta)[j]] = unlist(df_meta[,colnames(df_meta)[j]])
      }
      for(j in 1:length(df$time_series)){
        df$time_series[[j]] = merge(df$time_series[[j]], df_meta[j,])
      }
      df_long = do.call('rbind',df$time_series)
      for(j in 1:length(colnames(df_long))){
        df_long[,colnames(df_long)[j]] = unlist(df_long[,colnames(df_long)[j]])
      }
  
      event_list_long = dplyr::bind_rows(event_list_long, df_long)
      #print(i)
      i = i + 1
      OFFSET = OFFSET + 200
      if(length(event_list) < 200){
        STOP_FLAG = T
      }
    }
  }
  return(event_list_long)
}

#####################
### Save to Files ###
#####################

for(k in 1:nrow(disease_list)){
  data_out = read_aido(disease_num = unlist(disease_list$id[k]))
  write.csv(data_out, file = paste0("./raw_data/",gsub(' ', '_', unlist(disease_list$name[k])),".csv"), quote = T, row.names = F)
  print(k)
}
