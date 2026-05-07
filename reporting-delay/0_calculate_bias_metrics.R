#rm(list=ls())

################################################
### Exploring the GISAID Reporting Processes ###
################################################

# basic libraries
library(rstudioapi)
library(geomtextpath)
library(dplyr)  
library(lubridate)
library(data.table)
library(tidyr)

dir = "./precog/reporting-delay/"

### Read in data ###
merged_df <- fread(paste0(dir, "reporting_delay_data.csv")) %>%
  select(-V1)

##################################
###### Functions for running ##### 
##### by collection intervals ####
##################################

### aggregate by collection date intervals ###
get_AGG <- function(df, collect_date_int, delay_bin){
  
  max_submission_date <- collect_date_int + delay_bin
  during_delay <- paste0("[0,", delay_bin, "]")
  after_delay <- paste0("(", delay_bin, ",Inf]")
  
  collect_date_int_range <- seq(collect_date_int - 6, collect_date_int, by = "days")
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

### calculate Cramer's V and vector norms ###
run_cramer <- function(loc, df, time_int, delay_bin){
  
  time_int = as.Date(time_int, origin = '1970-01-01')
  print(paste(time_int, delay_bin, sep = ", "))
  
  # initialize values
  cramer_gof_debias = NA # debiased Cramer's V for goodness-of-fit
  n = NA # near-real-time sample size
  n_variant = NA # number of co-circulating variants
  error_code = NA
  L1 = NA # sum of component magnitudes
  L2 = NA # Euclidean norm
  Linf = NA # largest component magnitude
  omega = NA # reporting rate (proportion of samples reported within the delay period)
  pro_prop = NA # prospective (near-real-time) variant proportion
  ret_prop = NA # retrospective (validation) variant proportion
  
  # aggregate data
  AGG <- get_AGG(df, time_int, delay_bin)
  
  # format for Cramer's V test
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
      cramer_gof_debias = 0
      omega = 100
      L1 = 0
      L2 = 0
      Linf = 0
      
      # if nothing is reported before delay interval
    } else if(paste0("(", delay_bin, ",Inf]") %in% names(TABLE)){
      error_code = "No samples reported during delay period (1)"
      omega = 0
      
      ret_prop <- round(TABLE[, 2]/sum(TABLE[, 2])) 
      
    }
  }
  
  # if no samples were ever reported
  if(is.null(nrow(TABLE))){
    error_code = "No samples reported within time frame of interest"
    
  # if all samples were of the same variant
  } else if(nrow(TABLE) < 2){
    error_code = "All samples are of the same variant"
    cramer_gof_debias = 0
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
    omega = 0
    
    ret_prop <- round(TABLE[, 3]/sum(TABLE[, 3])) 
    
  } else if(ncol(TABLE) > 2){
    
    # if all samples were reported during the delay period
    if(sum(TABLE[, 3]) == 0){
      error_code = "All samples reported during delay period (2)"
      cramer_gof_debias = 0
      n <- sum(TABLE[, 2])
      n_variant <- nrow(TABLE)
      L1 = 0
      L2 = 0
      Linf = 0
      omega = 100
    
    # can calculate Cramer's V  
    } else {
      N <- sum(TABLE[, -1])
      n <- sum(TABLE[, 2])
      n_variant <- nrow(TABLE)
      omega <- (n/N)*100
      
      MAT <- as.matrix(TABLE[, -1])[, apply(as.matrix(TABLE[, -1]), 2, sum) > 0]
      
      # redefine second column of contingency matrix to be sum of all counts of each pango
      MAT[, 2] <- apply(MAT, 1, sum)
      
      # second column is the p-vector for goodness-of-fit test
      p_vec = MAT[, 2]/sum(MAT[, 2])
      
      pro_prop <- round(MAT[, 1]/sum(MAT[, 1]), 3)
      ret_prop <- p_vec
      p_dif <- ret_prop - pro_prop
      
      L1 = round(norm(as.matrix(p_dif), type = "O"), 3)
      L2 = round(norm(as.matrix(p_dif), type = "F"), 3)
      Linf = round(norm(as.matrix(p_dif), type = "I"), 3)
      
      # manual debiasing for GOF
      chi_gof <- as.numeric(chisq.test(MAT[, 1], p = p_vec, rescale.p = T)[1])
      chi_gof_debias <- max(c(0, (chi_gof/n) -(((nrow(MAT)-1)/(n-1)))))
      denom_debias <- nrow(MAT) - ((nrow(MAT) - 1)^2/(n-1)) - 1
      cramer_gof_debias <- sqrt(chi_gof_debias/denom_debias)
      cramer_gof_debias <- ifelse(cramer_gof_debias > 1, 1, cramer_gof_debias)
    
    }
  }
  
  new_results <-  data.frame(V = cramer_gof_debias,
                             K = n_variant,
                             delay_days = delay_bin,
                             Date = time_int,
                             n = n,
                             error_code = error_code,
                             L1 = L1,
                             L2 = L2,
                             Linf = Linf,
                             omega = omega)
  
  if(any(is.na(ret_prop))){
    new_results$p_nrt = NA
    new_results$p_val = NA
  } else {
    new_results$p_nrt <- list(pro_prop)
    new_results$p_val <- list(ret_prop)
  }
  
  return(new_results)
}

### save results ###
save_cramer <- function(loc, delay_val = 30){
  
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
      new_result <- run_cramer(loc, MERGED_DAT, i, delay_val)
      all_results <- rbind(all_results, new_result)
    }
  }
  
  full_results <- all_results %>% mutate(Location = loc)
  
  return(full_results)
}

############################
### run and save results ###
############################

### for all countries and delay periods ###
loc_vec <- unique(merged_df$Admin0)
full_results <- data.frame()

for(i in c(7, 14, 21, 30)){
  for(j in loc_vec){
    new_results <- save_cramer(j, delay_val = i)
    full_results <- rbind(full_results, new_results)
  }
}

# for simulations under the null
for_sim <- full_results %>%
  select(Date, Location, delay_days, V, L1, L2, Linf, n, omega, K, p_nrt, p_val) %>%
  arrange(Date, Location) %>%
  filter(K > 1,
         Date >= as.Date('2020-11-01'),
         Date <= as.Date('2022-12-31'),
         !is.na(n),
         !is.na(V))
row.names(for_sim) = 1:nrow(for_sim)

save(for_sim, file = paste0(dir, "all_cramer_results.RData"))

# for csv
save_results <- full_results %>%
  select(Date, Location, delay_days, V, L1, L2, Linf, n, omega, K, error_code)

write.csv(save_results, file = paste0(dir, "all_cramer_results.csv"))
