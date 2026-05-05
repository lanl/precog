#rm(list=ls())

#############################################
### find emerging variant characteristics ###
#############################################

# basic libraries
library(rstudioapi)
library(geomtextpath)
library(e1071)
library(lubridate)
library(data.table)
library(tidyr)
library(dplyr) 

dir = "./reporting-delay-pub/"

##############################
### Read in data and clean ###
##############################

merged_df <- read.csv(paste0(dir, "reporting_delay_data.csv")) %>%
  select(-X) %>%
  mutate(submission_date = collection_date + delay_days) %>%
  filter(counts > 0, year(collection_date) < 2023)

############################
### find p_val and p_nrt ###
############################

## only run this section once 

### aggregate by collection date intervals ###
get_AGG <- function(df, collect_date_int, delay_bin){
  
  max_submission_date <- collect_date_int + delay_bin
  during_delay <- paste0("[0,", delay_bin, "]")
  after_delay <- paste0("(", delay_bin, ",Inf]")
  
  collect_date_int_range <- seq(collect_date_int - 6, collect_date_int, by = "days")
  df <- df %>% filter(collection_date %in% collect_date_int_range) %>%
    mutate(submission_date = collect_date_int + delay_days,
           delay_days_cat = ifelse(submission_date <= max_submission_date,
                                   during_delay,
                                   after_delay),
           delay_days_cat = factor(delay_days_cat, levels = c(during_delay, after_delay)))
  
  # aggregate sequenced samples by other variables
  AGG = aggregate(counts ~ pango + delay_days_cat, 
                  data = df, FUN = sum, na.rm = T)
  
  return(AGG)
}

### calculate proportions ###
get_props <- function(loc, df, time_int, delay_bin){
  time_int = as.Date(time_int, origin = '1970-01-01')
  print(paste(time_int, delay_bin, sep = ", "))
  
  # aggregate data
  AGG <- get_AGG(df, time_int, delay_bin)
  
  TABLE <- pivot_wider(AGG, 
                       names_from = 'delay_days_cat', 
                       values_from = 'counts') %>% 
    as.data.frame()
  TABLE[is.na(TABLE)] <- 0
  
  # get proportions
  if(ncol(TABLE) > 2){
    ret_prop = round(apply(TABLE[, -1], 1, sum)/sum(TABLE[, -1]), 3)
    pro_prop = round(TABLE[, 2]/sum(TABLE[, 2]), 3)
    prop_dif = ret_prop - pro_prop
    
    new_table <- data.frame("Admin0" = loc,
                            "collection_date" = time_int,
                            "delay_bin" = delay_bin,
                            "pango" = TABLE[, 1],
                            "ret_prop" = ret_prop,
                            "pro_prop" = pro_prop,
                            "prop_dif" = prop_dif
    )
  } else {
    prop = round(TABLE[, 2]/sum(TABLE[, 2]), 3)
    
    # if everything is reported before delay interval
    if(paste0("[0,", delay_bin, "]") %in% names(TABLE)){
      
      new_table <- data.frame("Admin0" = loc,
                              "collection_date" = time_int,
                              "delay_bin" = delay_bin,
                              "pango" = TABLE[, 1],
                              "ret_prop" = prop,
                              "pro_prop" = prop,
                              "prop_dif" = 0)
      
      # if nothing is reported before delay interval
    } else if(paste0("(", delay_bin, ",Inf]") %in% names(TABLE)){
      new_table <- data.frame("Admin0" = loc,
                              "collection_date" = time_int,
                              "delay_bin" = delay_bin,
                              "pango" = TABLE[, 1],
                              "ret_prop" = prop,
                              "pro_prop" = 0,
                              "prop_dif" = prop)
    } else {
      new_table <- data.frame("Admin0" = loc,
                              "collection_date" = time_int,
                              "delay_bin" = delay_bin,
                              "pango" = NA,
                              "ret_prop" = NA,
                              "pro_prop" = NA,
                              "prop_dif" = NA)
    }
  }
  return(new_table)
}

### save results ###
save_results <- function(loc, delay_val = 30){
  
  # filter data by location
  sub_df <- merged_df %>%
    filter(Admin0 == loc)
  
  # vector of unique dates
  date_vec <- sort(unique(sub_df$collection_date))
  
  all_results <- data.frame()
  for(i in date_vec){
    if((sub_df %>% filter(collection_date == i) %>% nrow()) > 0){
      print(loc)
      
      new_result <- get_props(loc, sub_df, i, delay_val)
      all_results <- rbind(all_results, new_result)
    }
  }
  
  return(all_results)
}

loc_vec <- unique(merged_df$Admin0)

full_results <- data.frame()
for(i in c(7, 14, 21, 30)){
  for(j in loc_vec){
    new_results <- save_results(j, delay_val = i)
    full_results <- rbind(full_results, new_results)
  }
}

write.csv(full_results, file = paste0(dir, "emerging_variant_ts.csv"))
