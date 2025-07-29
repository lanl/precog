# Comparingincidence data model to proper use using new methods
## Author: Alexander C. Murph
## Date: July 2025

## Dave Osthus
## 1-18-24
## Recreating figs and code from Osthus et. al. 2010 AoAS paper
## "Forecasting Seasonal Influenza With A State-Space SIR Model"

## load libraries
library(rjags)
library(runjags)
library(ggplot2)
library(reshape2)
library(data.table)
library(plyr)
library(gridExtra)
library(tmvtnorm)
library(latex2exp)
library(this.path)
library(tidyverse)
setwd(paste0(this.path::here(), "/code"))
theme_set(theme_bw())

set.seed(13)

flu_seasons = c(2014) #, 2014, 2014
lengths = c(8:12)
final_epi_time = 35

graphic_idx = 1
# Set Flu Season
flu_season = flu_seasons[graphic_idx]
# Set forecast horizon


PIV_list = list()
PIT_list = list()

for(length_idx in 1:length(lengths)){
  length_of_time_series = lengths[graphic_idx]
  ## define paths
  ## Note: to recreate this experiment, get the data from Osthus et al 2017 and
  ## place this data in a directory located at the path below.
  datapath <- paste0(this.path::here(), "/data/")
  jagspath <- paste0(this.path::here(), "/code/")
  
  ## load ILI and surveillance data
  dfili <- read.csv(paste0(datapath,"EW08-2020_nat.csv"))
  dfili <- subset(dfili, region == "nat", select=c("epi_year","epi_week","wili","epi_time","epi_season"))
  dfsurveillance <- read.csv(paste0(datapath,"EW08-2020_who_combined_prior_to_2015_16.csv"))
  dfsurveillance <- subset(dfsurveillance, select=c("year","week","percent_positive"), region == "National")
  
  ## merge ILI data with surveillance data
  df <- merge(dfili, dfsurveillance, by.x=c("epi_year","epi_week"), by.y=c("year","week"))
  df$iliplus <- (df$wili/100)*(df$percent_positive/100)
  
  ## focus on 2010, national wILI (Fig 6 of Osthus et. al. AoAS)
  df_forecast_year <- subset(df, epi_season == flu_season & epi_time <= final_epi_time)
  df_forecast_year <- df_forecast_year[order(df_forecast_year$epi_time),]
  
  ###################################
  ## Prepare JAGS
  ## Equation 6.7
  params_eq6.7 <- c(1.62, 7084.10)
  
  ## Equation 6.10
  # The mean and Sigma here should be a truncated normal fit on PIT and PIV points prior to the date
  valid_flu_years = c(2002:2008, 2010:2013)
  valid_flu_years = valid_flu_years[valid_flu_years < flu_season]
  # Gather PIT and PIV
  qoi_df = NULL
  for(flu_year in valid_flu_years){
    df_tmp = subset(df, epi_season == flu_year & epi_time <= final_epi_time)
    df_tmp = df_tmp[order(df_tmp$epi_time),]
    qplot(epi_time, iliplus, data=df_tmp, geom=c("line","point"), group=epi_season)+
      ylab("Proportion Infectious")+
      xlab("t")
    idx_of_max = which.max(df_tmp$iliplus)
    qoi_df = rbind(qoi_df, 
                   data.frame(PIV = df_tmp$iliplus[idx_of_max],
                              PIT = df_tmp$epi_time[idx_of_max],
                              year = flu_year)
    )
  }
  data_mat <- as.matrix(qoi_df[, c("PIV", "PIT")])
  lower <- c(0.001, 1)
  upper <- c(1, 35)
  
  # Log-likelihood function
  neg_log_lik <- function(par) {
    mu <- par[1:2]
    L <- matrix(c(par[3], 0, par[4], par[5]), nrow = 2)  # lower-triangular
    Sigma <- L %*% t(L)
    
    # sum of -log densities
    -sum(dtmvnorm(
      x = data_mat,
      mean = mu,
      sigma = Sigma,
      lower = lower,
      upper = upper,
      log = TRUE
    ))
  }
  mu_init <- colMeans(data_mat)
  sd1 <- sd(data_mat[,1])
  sd2 <- sd(data_mat[,2])
  chol_init <- chol(matrix(c(sd1^2, 0, 0, sd2^2), 2, 2))
  
  par_init <- c(mu_init, chol_init[1,1], chol_init[2,1], chol_init[2,2])
  
  # Optimize
  fit <- optim(
    par = par_init,
    fn = neg_log_lik,
    method = "BFGS",
    control = list(maxit = 1000)
  )
  
  # Extract results
  mu_est <- fit$par[1:2]
  L_est <- matrix(c(fit$par[3], 0, fit$par[4], fit$par[5]), nrow = 2)
  Sigma_est <- L_est %*% t(L_est)
  
  mean_eq6.10 <- mu_est
  Sigma_eq6.10 <- Sigma_est
  samples_df = data.frame(rtmvnorm(1000, mean = mean_eq6.10, sigma = Sigma_eq6.10, lower = lower, upper = upper))
  colnames(samples_df) = c("PIV", "PIT")
  ggplot(samples_df, aes(x = PIT, y = PIV)) + geom_point(color = 'grey', size = 3) + 
    geom_point(data = qoi_df, aes(x = PIT, y = PIV), size = 5) + theme_bw()
  
  ## define number of chains to run
  nchains <- 1
  nsamples <- 5000
  
  
  ##########################################################################
  ################### Incidence Stuff ######################################
  #####################
  # Murph notes: It's here that my methods will update things, I believe.
  ## First of all, we no longer need the rho grids.  We will instead use
  ## the grid of dave_beta, dave_gamma values.
  # Compile Grid from the Darwin files.
  inputgrid = NULL
  for(file_name in list.files("grid_logs")){
    data_name   = paste("grid_logs/", file_name, sep = "")
    temp_rows   = read.csv(data_name)
    temp_rows$X = NULL
    inputgrid   = rbind(inputgrid, temp_rows)
  }
  #####################
  
  # Murph: so dave used a look up table, and then linearly interpolates inside of JAGS.
  # Pretty cool.  I'm also using linear interpolation, but over three variables
  # (trilinear interpolation).
  
  #########################################
  ## Fit Model in JAGS
  
  ## the observed time series
  myy <- df_forecast_year[df_forecast_year$epi_time %in% 1:length_of_time_series,]$iliplus
  
  ## construct data to be passed into JAGS
  jags_data <- list(y = myy,
                    tt = length(myy),
                    TT = 35,
                    params_eq6.7 = params_eq6.7,
                    mean_eq6.10 = mean_eq6.10,
                    Sigma_eq6.10 = Sigma_eq6.10,
                    dave_betas = inputgrid$dave_betas,
                    dave_gammas = inputgrid$dave_gamma)
  
  ## run JAGS
  results_incidence <- run.jags(model = paste0(jagspath,"dbssm_jags_incidence.txt"),
                                data = jags_data,
                                monitor = c("ypred","theta","beta","gamma"),
                                n.chains = nchains,
                                burnin = 12500,
                                sample = nsamples,
                                thin = 10)
  
  ## organize the posterior draws from JAGS
  postdraws_incidence <- NULL
  full_data = NULL
  for(i in 1:nchains){
    full_data = rbind(full_data, data.frame(chain = i, id = 1:nrow(data.frame(results_incidence$mcmc[[i]])), 
                                            data.frame(results_incidence$mcmc[[i]])))
    temppostdraws_incidence <- reshape2::melt(data.frame(chain = i, id = 1:nrow(data.frame(results_incidence$mcmc[[i]])), 
                                                         data.frame(results_incidence$mcmc[[i]])),id.vars = c("chain","id"))
    postdraws_incidence <- rbind(postdraws_incidence, temppostdraws_incidence)
  }
  
  ## organize output
  postdraws_incidence$type <- data.table::tstrsplit(postdraws_incidence$variable,"\\.")[[1]]
  postdraws_incidence$time <- as.numeric(data.table::tstrsplit(postdraws_incidence$variable,"\\.")[[2]])
  postdraws_incidence$theta_type <- as.numeric(data.table::tstrsplit(postdraws_incidence$variable,"\\.")[[3]])
  
  ## just keep ypred and theta_type == 2
  incidence_gammas = subset(postdraws_incidence, type == "gamma")
  incidence_betas  = subset(postdraws_incidence, type == "beta")
  
  postdraws_incidence <- subset(postdraws_incidence, (theta_type == 4 & time <= length_of_time_series + 1) | (type == "ypred" & time > length_of_time_series))
  postdraws_incidence <- subset(postdraws_incidence, select=c("chain","id","type","time","value"))
  
  ## theta time is off shifted by 1
  postdraws_incidence$theta_time <- postdraws_incidence$time - 1
  postdraws_incidence$unqid      <- paste0(postdraws_incidence$chain,"_",postdraws_incidence$id)
  
  ## compute summaries
  postdrawsply_incidence            <- ddply(postdraws_incidence,.(time, type),summarise,lower = quantile(value,probs=.025), avg = mean(value), upper = quantile(value,probs=.975))
  postdrawsply_incidence$theta_time <- postdrawsply_incidence$time - 1
  
  # Use posterior draws to get PIT and PIV draws.
  dbssm_pivs = c()
  dbssm_pits = c()
  for(id_idx in unique(postdraws_incidence$id)){
    tmp_df = postdraws_incidence %>% subset(id == id_idx)
    PIV = max(tmp_df$value)
    PIT = tmp_df$theta_time[which.max(tmp_df$value)]
    dbssm_pivs = c(dbssm_pivs, PIV)
    dbssm_pits = c(dbssm_pits, PIT)
  }
  PIV_list[[length_idx]] = dbssm_pivs
  PIT_list[[length_idx]] = dbssm_pits
}
##########################################################################
##########################################################################

# Do the stochastic thing from Hickman 2015
subset_df = df %>% 
  subset((epi_season %in% valid_flu_years) & (epi_time %in% c(1:final_epi_time))) %>%
  as_tibble()
history_df = matrix(0, nrow = final_epi_time, ncol = length(valid_flu_years))
for(season_idx in 1:length(valid_flu_years)){
  tmp_df = subset_df[which(subset_df$epi_season == valid_flu_years[season_idx]),]
  tmp_df = tmp_df[order(tmp_df$epi_time),]
  history_df[,season_idx] = tmp_df$iliplus
}
sds_of_epi_weeks = apply(history_df, 1, sd)
means_of_epi_weeks = apply(history_df, 1, mean)

stochastic_timeseries = rmvnorm(5000, mean = means_of_epi_weeks, sigma = diag(sds_of_epi_weeks))
PIV_of_stochastic_model = apply(stochastic_timeseries, 1, max)
PIT_of_stochastic_model = apply(stochastic_timeseries, 1, which.max)

# True PIT & PIV
tmp_df = df[which(df$epi_season == flu_seasons[graphic_idx]),]
tmp_df = tmp_df[order(tmp_df$epi_time),]
true_PIV = max(tmp_df$iliplus)
true_PIT = which.max(tmp_df$iliplus)

######################################################################
## Show the PIV, PIT distributions
PIV_df = data.frame(PIV = PIV_of_stochastic_model, Model = rep("Hickman", times = length(PIV_of_stochastic_model)) )
for(length_idx in 1:length(lengths)){
  PIV_df = rbind(PIV_df,
                  data.frame(PIV = PIV_list[[length_idx]], Model = rep(paste0("DBSSM, Horizon = ", 
                                                                              lengths[length_idx]), 
                                                                              times = length(PIV_of_stochastic_model) )
                  ))
}
PIT_df = data.frame(PIT = PIT_of_stochastic_model, Model = rep("Hickman", times = length(PIV_of_stochastic_model)) )
for(length_idx in 1:length(lengths)){
  PIT_df = rbind(PIT_df,
                 data.frame(PIT = PIT_list[[length_idx]], Model = rep(paste0("DBSSM, Horizon = ", 
                                                                             lengths[length_idx]), 
                                                                             times = length(PIV_of_stochastic_model) )
                 ))
}

p1 = PIV_df %>%
  mutate(Model = factor(Model, levels = c("Hickman", "DBSSM, Horizon = 4",
                                          "DBSSM, Horizon = 5",
                                          "DBSSM, Horizon = 6",
                                          "DBSSM, Horizon = 7",
                                          "DBSSM, Horizon = 8",
                                          "DBSSM, Horizon = 9",
                                          "DBSSM, Horizon = 10"))) %>%
  ggplot(aes(x = PIV, fill = Model)) +
  geom_density(alpha = 0.7) + 
  geom_vline(xintercept = true_PIV) + 
  xlim(0,0.2)+ 
  theme(axis.text=element_text(size=18), axis.title = element_text(size = 18), plot.title = element_text(size = 18)) + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1))
p1

p2 = PIT_df %>%
  mutate(Model = factor(Model, levels = c("Hickman", "DBSSM, Horizon = 4",
                                          "DBSSM, Horizon = 5",
                                          "DBSSM, Horizon = 6",
                                          "DBSSM, Horizon = 7",
                                          "DBSSM, Horizon = 8",
                                          "DBSSM, Horizon = 9",
                                          "DBSSM, Horizon = 10"))) %>%
  ggplot(aes(x = PIT, fill = Model)) +
  geom_histogram(alpha = 0.7, color = 'black', position = 'stack') + 
  geom_vline(xintercept = true_PIT) + 
  theme(axis.text=element_text(size=18), axis.title = element_text(size = 18), plot.title = element_text(size = 18)) + 
  theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1))
p2
##########################################################################
##########################################################################å

# save the plots
both_plots <- list(p1, p2)
grid.arrange(both_plots[[1]], both_plots[[2]],
             layout_matrix = matrix(c(1,2), byrow = TRUE, ncol = 1))
saveRDS(both_plots, file = paste0("../figures/pivpitplots.rds"))
