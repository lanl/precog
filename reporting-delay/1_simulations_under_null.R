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

###########################
### load and clean data ###
###########################

dir = "./precog/reporting-delay/"

### Read in data ###

load(paste0(dir, "all_cramer_results.RData"))

df <- for_sim

#####################################################
### set up functions prior to running simulations ###
#####################################################

# manual debiasing for Cramer's V goodness-of-fit
calc_CV <- function(Y, n, p_val, K){
  
  chisq <- as.numeric(chisq.test(Y, p = p_val, rescale.p = T)[1])
  chisq_debias <- max(c(0, (chisq/n) -(((K-1)/(n-1)))))
  denom_debias <- K - ((K- 1)^2/(n-1)) - 1
  CV_debias <- sqrt(chisq_debias/denom_debias)
  CV_debias <- ifelse(CV_debias > 1, 1, CV_debias)
  
  return(CV_debias)
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

# run simulations under null (no bias) for observed n, p_val, and K
run_sim <- function(n, p_val, K, true_V, true_L1, true_L2, true_Linf){
  set.seed(12)
  M = 1000
  Y_sim = rmultinom(n = M, size = n, prob = p_val)
  
  # Cramer's V
  V_sim <- apply(Y_sim, 2, calc_CV, n = n, p_val = p_val, K = K)
  V_pvalue <- sum(true_V <= V_sim, na.rm = T)/M
  V_95 = quantile(V_sim, probs = 0.95, na.rm = T)[[1]]
  
  # Norms
  p_sim = Y_sim/n
  
  L1_sim <- apply(p_sim, 2, calc_norms, p_val = p_val, norm_type = "L1")
  L1_sim <- L1_sim/2  # scale between 0 and 1
  L1_pvalue <- sum(true_L1 <= L1_sim)/M
  L1_95 = quantile(L1_sim, probs = 0.95, na.rm = T)[[1]]
  L1_95 = L1_95/2   # scale between 0 and 1
  
  L2_sim <- apply(p_sim, 2, calc_norms, p_val = p_val, norm_type = "L2")
  L2_sim <- L2_sim/sqrt(2)    # scale between 0 and 1
  L2_pvalue <- sum(true_L2 <= L2_sim)/M
  L2_95 = quantile(L2_sim, probs = 0.95, na.rm = T)[[1]]
  L2_95 = L2_95/sqrt(2)   # scale between 0 and 1
  
  Linf_sim <- apply(p_sim, 2, calc_norms, p_val = p_val, norm_type = "Linf")
  Linf_pvalue <- sum(true_Linf <= Linf_sim)/M
  Linf_95 = quantile(Linf_sim, probs = 0.95, na.rm = T)[[1]]
  
  # create data frames
  pvalue_df <- data.frame(
    V = V_pvalue,
    L1 = L1_pvalue,
    L2 = L2_pvalue,
    Linf = Linf_pvalue
  )
  
  q95_df <- data.frame(
    V = V_95,
    L1 = L1_95,
    L2 = L2_95,
    Linf = Linf_95
  )
  
  return(list("pvalues" = pvalue_df, "q95" = q95_df))
}

# run simulations for a specific location and delay period
sim_loc_delay <- function(loc, delay, filter_df){
  print(paste(loc, delay, sep = ", "))

  pvalues <- data.frame()
  q95 <- data.frame()
  
  for(i in 1:nrow(filter_df)){
    sim_params <- filter_df[i, ]
    n <- sim_params$n
    p_val <- sim_params$p_val[[1]]
    K <- sim_params$K
    true_V <- sim_params$V
    true_L1 <- (sim_params$L1)/2
    true_L2 <- (sim_params$L2)/sqrt(2)
    true_Linf <- sim_params$Linf
    
    if(any(is.na(p_val))){
      new_pvalues <- data.frame(V = 0, L1 = 0, L2 = 0, Linf = 0)
      new_q95 <- data.frame(V = 0, L1 = 0, L2 = 0, Linf = 0)
    } else {
      new_vals <- run_sim(n, p_val, K, true_V, true_L1, true_L2, true_Linf)
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

###########################################################
### run simulations for all locations and delay periods ###
###########################################################

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
      new_sim <- sim_loc_delay(loc, delay, filter_df)
      pvalue_df <- rbind(pvalue_df, new_sim$pvalues)
      q95_df <- rbind(q95_df, new_sim$q95)
    }
  }
}

########################
### organize results ###
########################

df_mod <- df %>%
  select(-c(K, p_nrt, p_val)) %>%
  mutate(L1 = L1/2,
         L2 = L2/sqrt(2)) %>%
  pivot_longer(cols = c(V, L1, L2, Linf), names_to = "Metric",
               values_to = "True_Value")

pvalues_mod <- pvalue_df %>%
  pivot_longer(cols = c(V, L1, L2, Linf), names_to = "Metric",
               values_to = "pvalue")

q95_mod <- q95_df %>%
  pivot_longer(cols = c(V, L1, L2, Linf), names_to = "Metric",
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
         omega_cat = cut(omega, 
                         breaks = seq(0, 100, 10), 
                         labels = c("(0, 10]", "(10, 20]", "(20, 30]", "(30, 40]",
                                    "(40, 50]", "(50, 60]", "(60, 70]", "(70, 80]",
                                    "(80, 90]", "(90, 100)"), 
                         include.lowest = T))

save(metric_comb, file = paste0(dir, "all_sim_results.RData"))
