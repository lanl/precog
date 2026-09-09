#rm(list=ls())

################################################
### Exploring the GISAID Reporting Processes ###
################################################

# load libraries
library(rstudioapi)
library(geomtextpath)
library(dplyr)  
library(lubridate)
library(data.table)
library(tidyr)
library(effectsize)

dir = "./precog/reporting-delay/"
results_dir = paste0(dir, "results/")

### Read in data ###
merged_df <- fread(paste0(dir, "reporting_delay_data.csv")) %>%
  select(-V1)

##################################
###### Functions for running ##### 
##### by collection intervals ####
##################################

### aggregate by collection date intervals ###
get_AGG <- function(df, collect_date_int, delay_bin, sliding_window){
  
  max_submission_date <- collect_date_int + delay_bin
  during_delay <- paste0("[0,", delay_bin, "]")
  after_delay <- paste0("(", delay_bin, ",Inf]")
  
  collect_date_int_range <- seq(collect_date_int - sliding_window, collect_date_int, by = "days")
  df <- df %>% filter(as.Date(collection_date_int) %in% collect_date_int_range) %>%
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

### calculate bias metrics ###
calc_bias <- function(loc, df, time_int, delay_bin, sliding_window){
  
  time_int = as.Date(time_int, origin = '1970-01-01')
  print(paste(time_int, delay_bin, sep = ", "))
  
  # initialize values
  w = NA # Cohen's w for chi-square goodness-of-fit
  Fei = NA # Fei for chi-sqare-goodness-of-fit
  n = NA # near-real-time sample size
  n_variant = NA # number of co-circulating variants
  error_code = NA
  L1 = NA # sum of component magnitudes
  L2 = NA # Euclidean norm
  Linf = NA # largest component magnitude
  r = NA # reporting rate (percentage of samples reported within the delay period)
  p_nrt = NA # near-real-time variant proportion
  p_val = NA # validation variant proportion
  
  # aggregate data for sliding window
  AGG <- get_AGG(df, time_int, delay_bin, sliding_window)
  
  # format
  TABLE <- pivot_wider(AGG, names_from = 'delay_days_cat', 
                       values_from = 'counts') %>% 
    as.data.frame()
  TABLE[is.na(TABLE)] <- 0
  
  # if only two columns are part of TABLE, then one time frame had no samples reported
  if(ncol(TABLE) == 2){
    n <- sum(TABLE[, -1])
    n_variant <- nrow(TABLE)
    
    # if everything is reported before delay interval
    if(paste0("[0,", delay_bin, "]") %in% names(TABLE)){
      error_code = "All samples reported during delay period (1)"
      w = 0
      Fei = 0
      r = 100
      L1 = 0
      L2 = 0
      Linf = 0
      
    # if nothing is reported before delay interval
    } else if(paste0("(", delay_bin, ",Inf]") %in% names(TABLE)){
      error_code = "No samples reported during delay period (1)"
      r = 0
      
      p_val <- round(TABLE[, 2]/sum(TABLE[, 2])) 
      
    }
  }
  
  # if no samples were ever reported
  if(is.null(nrow(TABLE))){
    error_code = "No samples reported within time frame of interest"
    
  # if all samples were of the same variant
  } else if(nrow(TABLE) < 2){
    error_code = "All samples are of the same variant"
    w = 0
    Fei = 0
    n <- sum(TABLE[, 2])
    n_variant <- nrow(TABLE)
    L1 = 0
    L2 = 0
    Linf = 0
    
  # if all samples were reported after the delay period
  } else if(sum(TABLE[, 2]) == 0){
    error_code = "No samples reported during delay period (3)"
    n <- 0
    n_variant <- nrow(TABLE)
    r = 0
    
    p_val <- round(TABLE[, 3]/sum(TABLE[, 3])) 
  
  # otherwise, can calculate metrics    
  } else if(ncol(TABLE) > 2){
    
    # if all samples were reported during the delay period
    if(sum(TABLE[, 3]) == 0){
      error_code = "All samples reported during delay period (2)"
      w = 0
      Fei = 0
      n <- sum(TABLE[, 2])
      n_variant <- nrow(TABLE)
      L1 = 0
      L2 = 0
      Linf = 0
      r = 100
    
    # can calculate metrics
    } else {
      N <- sum(TABLE[, -1])
      n <- sum(TABLE[, 2])
      n_variant <- nrow(TABLE)
      r <- (n/N)*100
      
      MAT <- as.matrix(TABLE[, -1])[, apply(as.matrix(TABLE[, -1]), 2, sum) > 0]
      
      # redefine second column of matrix to be sum of all counts of each pango
      MAT[, 2] <- apply(MAT, 1, sum)
      
      # second column is the p-vector for goodness-of-fit test
      p_vec = MAT[, 2]/sum(MAT[, 2])
      
      p_nrt <- round(MAT[, 1]/sum(MAT[, 1]), 3)
      p_val <- p_vec
      p_dif <- p_val - p_nrt
      
      # calculate vector norms
      L1 = round(norm(as.matrix(p_dif), type = "O"), 3)
      L2 = round(norm(as.matrix(p_dif), type = "F"), 3)
      Linf = round(norm(as.matrix(p_dif), type = "I"), 3)
      
      # calculate Cohen's w
      chi_gof <- as.numeric(chisq.test(MAT[, 1], p = p_vec, rescale.p = T)[1])
      w <- sqrt(chi_gof/n)
      
      # calculate fei
      fei_result <- fei(MAT[, 1], p = p_vec, ci = NULL)
      Fei <- as.numeric(fei_result$Fei)
    
    }
  }
  
  new_results <-  data.frame(w = w,
                             Fei = Fei,
                             K = n_variant,
                             delay_days = delay_bin,
                             Date = time_int,
                             n = n,
                             error_code = error_code,
                             L1 = L1,
                             L2 = L2,
                             Linf = Linf,
                             r = r)
  
  if(any(is.na(p_val))){
    new_results$p_nrt = NA
    new_results$p_val = NA
  } else {
    new_results$p_nrt <- list(p_nrt)
    new_results$p_val <- list(p_val)
  }
  
  return(new_results)
}

### save metrics ###
save_metrics <- function(loc, delay_val = 30, sliding_window){
  
  # filter data by location
  MERGED_DAT <- merged_df %>%
    filter(Admin0 == loc) %>%
    ungroup() %>%
    dplyr::select(collection_date, delay_days, pango, counts) %>%
    mutate(collection_date_int = collection_date)
  
  # vector of unique dates
  date_vec <- sort(unique(MERGED_DAT$collection_date))
  
  all_results <- data.frame()
  for(i in date_vec){
    if((MERGED_DAT %>% filter(collection_date_int == i) %>% nrow()) > 0){
      print(loc)
      new_result <- calc_bias(loc, MERGED_DAT, i, delay_val, sliding_window)
      all_results <- rbind(all_results, new_result)
    }
  }
  
  full_results <- all_results %>% mutate(Location = loc)
  
  return(full_results)
}

############################
### run and save results ###
############################

run_and_save <- function(merged_df, sliding_window, csv_name, save_RData = F){

  loc_vec <- unique(merged_df$Admin0)
  full_results <- data.frame()
  
  for(i in c(7, 14, 21, 30)){
    for(j in loc_vec){
      new_results <- save_metrics(j, delay_val = i, sliding_window = sliding_window)
      full_results <- rbind(full_results, new_results)
    }
  }
  
  # for simulations under the null
  if(save_RData == T){
    for_sim <- full_results %>%
      select(Date, Location, delay_days, w, Fei, L1, L2, Linf, n, r, K, p_nrt, p_val) %>%
      arrange(Date, Location) %>%
      filter(K > 1,
             Date >= as.Date('2020-11-01'),
             Date <= as.Date('2022-12-31'),
             !is.na(n),
             !is.na(w))
    row.names(for_sim) = 1:nrow(for_sim)
    
    save(for_sim, file = paste0(results_dir, "all_metric_results.RData"))
  }
  
  # for csv
  save_results <- full_results %>%
    select(Date, Location, delay_days, w, Fei, L1, L2, Linf, n, r, K, error_code)
  
  write.csv(save_results, file = paste0(results_dir, csv_name))
}

### all metric results, 7 day sliding window ###
run_and_save(merged_df, sliding_window = 6, csv_name = "all_metric_results.csv",
             save_RData = T)

### sensitivity analysis for 3 and 14 day sliding windows ###
sub_df <- merged_df %>%
  filter(Admin0 %in% c("Brazil", "Denmark", "United States"))
  
run_and_save(sub_df, sliding_window = 2, csv_name = "all_metric_results_3_window.csv")
run_and_save(sub_df, sliding_window = 13, csv_name = "all_metric_results_14_window.csv")

### sensitivity analysis for omicron categories ###
collapse_omicron <- sub_df %>%
  mutate(pango = case_when(pango %in% c("Omicron BA.1", "Omicron BA.1.1") ~ "Omicron Cat 1",
                           pango %in% c("Omicron BA.2", "Omicron BA.2.12.1", "Omicron BA.2.75") ~ "Omicron Cat 2",
                           pango %in% c("Omicron BA.4", "Omicron BA.5") ~ "Omicron Cat 3",
                           pango %in% c("Omicron BQ.1", "Omicron XBB") ~ "Omicron Cat 4",
                           .default = pango)
  )

run_and_save(collapse_omicron, sliding_window = 6, csv_name = "all_metric_results_collapse_omicron.csv")
