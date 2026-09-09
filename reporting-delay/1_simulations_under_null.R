#rm(list=ls())

######################################################
### Bias metrics simulations under null hypothesis ###
######################################################

# basic libraries
library(rstudioapi)
library(geomtextpath)
library(dplyr)  
library(lubridate)
library(data.table)
library(tidyr)
library(effectsize)

###########################
### load and clean data ###
###########################

dir = "./precog/reporting-delay/"
results_dir = paste0(dir, "results/")

### Read in data ###

load(paste0(results_dir, "all_metric_results.RData"))

df <- for_sim

#####################################################
### set up functions prior to running simulations ###
#####################################################

# Cohen's w for chi-square goodness-of-fit
calc_w <- function(Y, n, p_val){
  chisq <- as.numeric(chisq.test(Y, p = p_val, rescale.p = T)[1])
  return(sqrt(chisq/n))
}

calc_fei <- function(Y, p_val){
  fei_result <- fei(Y, p = p_val, ci = NULL)
  return(as.numeric(fei_result$Fei))
}

# calculate vector norms
calc_norms <- function(p_sim, p_val, norm_type){
  
  p_dif <- p_val - p_sim  # proportion difference (true bias)
  
  if(norm_type == "L1"){
    type = "O"  # sum of component magnitudes
  } else if(norm_type == "L2"){
    type = "F"  # Euclidean norm
  } else {
    type = "I"  # component with largest magnitude
  }
  
  norm = round(norm(as.matrix(p_dif), type = type), 3)
  return(norm)
}

# run simulations under null (no bias) for observed n and p_val
run_sim <- function(n, p_val, true_w, true_fei, true_L1, true_L2, true_Linf, M){
  set.seed(12)
  Y_sim = rmultinom(n = M, size = n, prob = p_val)
  
  # Cohen's w
  w_sim <- apply(Y_sim, 2, calc_w, n = n, p_val = p_val)
  w_pvalue <- sum(true_w <= w_sim, na.rm = T)/M
  w_95 = quantile(w_sim, probs = 0.95, na.rm = T)[[1]]
  
  # Fei
  fei_sim <- apply(Y_sim, 2, calc_fei, p_val = p_val)
  fei_pvalue <- sum(true_fei <= fei_sim, na.rm = T)/M
  fei_95 = quantile(fei_sim, probs = 0.95, na.rm = T)[[1]]
  
  # Norms
  p_sim = Y_sim/n
  
  L1_sim <- apply(p_sim, 2, calc_norms, p_val = p_val, norm_type = "L1")
  L1_pvalue <- sum(true_L1 <= L1_sim)/M
  L1_95 = quantile(L1_sim, probs = 0.95, na.rm = T)[[1]]
  
  L2_sim <- apply(p_sim, 2, calc_norms, p_val = p_val, norm_type = "L2")
  L2_pvalue <- sum(true_L2 <= L2_sim)/M
  L2_95 = quantile(L2_sim, probs = 0.95, na.rm = T)[[1]]
  
  Linf_sim <- apply(p_sim, 2, calc_norms, p_val = p_val, norm_type = "Linf")
  Linf_pvalue <- sum(true_Linf <= Linf_sim)/M
  Linf_95 = quantile(Linf_sim, probs = 0.95, na.rm = T)[[1]]
  
  # create data frames
  pvalue_df <- data.frame(
    w = w_pvalue,
    Fei = fei_pvalue,
    L1 = L1_pvalue,
    L2 = L2_pvalue,
    Linf = Linf_pvalue
  )
  
  q95_df <- data.frame(
    w = w_95,
    Fei = fei_95,
    L1 = L1_95,
    L2 = L2_95,
    Linf = Linf_95
  )
  
  return(list("pvalues" = pvalue_df, "q95" = q95_df))
}

# run simulations for a specific location and delay period
sim_loc_delay <- function(loc, delay, filter_df, M = 1000){
  print(paste(loc, delay, sep = ", "))

  pvalues <- data.frame()
  q95 <- data.frame()
  
  for(i in 1:nrow(filter_df)){
    sim_params <- filter_df[i, ]
    n <- sim_params$n
    p_val <- sim_params$p_val[[1]]
    K <- sim_params$K
    true_w <- sim_params$w
    true_fei <- sim_params$Fei
    true_L1 <- sim_params$L1
    true_L2 <- sim_params$L2
    true_Linf <- sim_params$Linf
    
    if(any(is.na(p_val))){
      new_pvalues <- data.frame(w = 0, Fei = 0, L1 = 0, L2 = 0, Linf = 0)
      new_q95 <- data.frame(w = 0, Fei = 0, L1 = 0, L2 = 0, Linf = 0)
    } else {
      new_vals <- run_sim(n, p_val, true_w, true_fei, true_L1, true_L2, true_Linf, M)
      new_pvalues <- new_vals$pvalues
      new_q95 <- new_vals$q95
    }
    
    new_pvalues$Date <- sim_params$Date
    new_q95$Date <- sim_params$Date
    new_pvalues$Location <- sim_params$Location
    new_q95$Location <- sim_params$Location
    new_pvalues$delay_days <- sim_params$delay_days
    new_q95$delay_days <- sim_params$delay_days
    pvalues <- rbind(pvalues, new_pvalues)
    q95 <- rbind(q95, new_q95)
  }
  
  return(list("pvalues" = pvalues, "q95" = q95))
}

####################################
### post-processing organization ###
####################################

### organize results ###
organize_results <- function(pvalue_df, q95_df, file_name){
  df_mod <- df %>%
    select(-c(K, p_nrt, p_val)) %>%
    pivot_longer(cols = c(w, Fei, L1, L2, Linf), names_to = "Metric",
                 values_to = "True_Value")
  
  pvalues_mod <- pvalue_df %>%
    pivot_longer(cols = c(w, Fei, L1, L2, Linf), names_to = "Metric",
                 values_to = "pvalue")
  
  q95_mod <- q95_df %>%
    pivot_longer(cols = c(w, Fei, L1, L2, Linf), names_to = "Metric",
                 values_to = "q95")
  
  metric_comb <- inner_join(df_mod, pvalues_mod) %>%
    inner_join(q95_mod) %>%
    mutate(pvalue_cat = cut(pvalue,
                            breaks = c(-Inf, 0.05, Inf),
                            labels = c("< 0.05", "> 0.05"),
                            right = F),
           n_cat = cut(n,
                       breaks = c(2, 10, 100, 1000, 10000, Inf),
                       labels = c("[0, 10]", "(10, 100]", "(100, 1k]", "(1k, 10k]",
                                  "> 10k"),
                       include.lowest = T),
           r_cat = cut(r, 
                       breaks = seq(0, 100, 10), 
                       labels = c("(0, 10]", "(10, 20]", "(20, 30]", "(30, 40]",
                                  "(40, 50]", "(50, 60]", "(60, 70]", "(70, 80]",
                                  "(80, 90]", "(90, 100)"), 
                       include.lowest = T))
  
  save(metric_comb, file = paste0(results_dir, file_name))
}

#######################
### run simulations ###
#######################

find_pvalues_and_q95 <- function(M, file_name){
  loc_vec <- unique(df$Location) %>% sort()
  delay_vec <- unique(df$delay_days)
  
  pvalue_df <- data.frame()
  q95_df <- data.frame()
  
  for(delay in delay_vec){
    for(loc in loc_vec){
      
      filter_df <- df %>%
        filter(Location == loc,
               delay_days == delay)
      
      if(nrow(filter_df) > 0){
        new_sim <- sim_loc_delay(loc, delay, filter_df, M)
        pvalue_df <- rbind(pvalue_df, new_sim$pvalues)
        q95_df <- rbind(q95_df, new_sim$q95)
      }
    }
  }
  
  organize_results(pvalue_df, q95_df, file_name)
  
}

find_pvalues_and_q95(M = 1000, "all_sim_results.RData")
find_pvalues_and_q95(M = 500, "sim_results_M_500.RData")
find_pvalues_and_q95(M = 5000, "sim_results_M_5000.RData")
