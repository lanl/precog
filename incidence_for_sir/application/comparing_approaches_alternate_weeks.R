# Comparing erroneous use of incidence data to proper use using new methods
## Author: Alexander C. Murph
## Date: March 2024

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
setwd(paste0(this.path::here(), "/code"))
theme_set(theme_bw())

set.seed(13)

flu_seasons = c(2010, 2010) #, 2014, 2014
lengths = c(11, 14, 11, 14)

for(graphic_idx in 1:length(flu_seasons)){
  # Set Flu Season
  flu_season = flu_seasons[graphic_idx]
  # Set forecast horizon
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
  final_epi_time = 35
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
                      # monitor = c("beta","gamma", "z1", "z2", "theta_i0"),
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
  
  
  ##########################################################################
  ##########################################################################
  
  
  
  ##########################################################################
  ################### Prevalence Stuff #####################################
  ## Make grid for g^-1(.) in equation 6.12
  inputgrid <- expand.grid(S0 = 0.9,
                           I0 = seq(0,.05,length.out = 100),
                           PI = seq(0,.1,length.out = 100))
  inputgrid <- subset(inputgrid, PI >= I0)
  inputgrid <- subset(inputgrid, select=c("S0","I0","PI"))
  
  ## Equation 6.11
  f_eq6.11 <- function(x,S0,I0,PI){
    rhs = I0 + S0 - x*(log(S0) + 1 - log(x))
    squared_error <- (PI - rhs)^2
    return(squared_error)
  }
  
  ## solve for rho (assuming rho < 1) for all scenarios
  inputgrid$rho <- NA
  for(i in 1:nrow(inputgrid)){
    print(i)
    rhofit <- optimize(f = f_eq6.11,
                       interval = c(0,1), # this assumes an epidemic (i.e., rho < 1), as is done in the paper between equations 6.11 and 6.12
                       S0 = inputgrid$S0[i],
                       I0 = inputgrid$I0[i],
                       PI = inputgrid$PI[i])
    inputgrid$rho[i] <- rhofit$minimum
  }
  
  ## JAGS can't do any kind of solving (like optimize() or uniroot()) within it, so we need a hack. 
  ## JAGS has a linear interpolation function we can make use of
  ## Equation 6.11: PI = I0 - S0 - rho*(log(S0) + 1 - log(rho)), which implies 
  ## PI - I0 - S0 = -rho*(log(S0) + 1 - log(rho))
  ## In Osthus 2017, he assumes S0 = .9, thus the RHS of the above equation is only a function of rho.
  ## Set rhostar = PI - I0 - S0
  inputgrid$rhostar <- inputgrid$PI - inputgrid$I0 - inputgrid$S0
  inputgridsmall <- ddply(inputgrid, .(rhostar, rho), summarise, n=1)
  
  ## make a grid of rhostar and rho
  ## this will be passed into JAGS to compute rho GIVEN PI, I0, and S0
  rhostargrid <- seq(min(inputgridsmall$rhostar),
                     max(inputgridsmall$rhostar),
                     length.out = 200)
  rhogrid <- approx(x=inputgridsmall$rhostar,
                    y=inputgridsmall$rho,
                    xout = rhostargrid)$y
  
  ## Equation 6.13
  tauhat_eq6.13 <- c(-49.754,
                     -0.9577,
                     -0.0065,
                     -9.4896,
                     -0.3761,
                     -590.0001,
                     -2537.6102,
                     -4756.1828,
                     -3265.2458,
                     -102.2665,
                     -4.0162,
                     -430.9596,
                     -16.7104,
                     -798.3443,
                     -30.6638,
                     -543.8857,
                     -20.7459)
  
  sigma2hat_eq6.13 <- 0.0421^2
  
  
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
                    rhostargrid  = rhostargrid,
                    rhogrid = rhogrid,
                    tauhat_eq6.13 = tauhat_eq6.13,
                    sigma2hat_eq6.13 = sigma2hat_eq6.13)
  
  ## run JAGS
  results_prevalence <- run.jags(model = paste0(jagspath,"dbssm_jags.txt"),
                      data = jags_data,
                      monitor = c("ypred","theta","beta","gamma"),
                      n.chains = nchains,
                      burnin = 12500,
                      sample = nsamples,
                      thin = 10)
  
  ## organize the posterior draws from JAGS
  postdraws_prevalence <- NULL
  for(i in 1:nchains){
    temppostdraws_prevalence <- reshape2::melt(data.frame(chain = i, id = 1:nrow(data.frame(results_prevalence$mcmc[[i]])), data.frame(results_prevalence$mcmc[[i]])),id.vars = c("chain","id"))
    postdraws_prevalence <- rbind(postdraws_prevalence, temppostdraws_prevalence)
  }
  
  ## organize output
  postdraws_prevalence$type <- data.table::tstrsplit(postdraws_prevalence$variable,"\\.")[[1]]
  postdraws_prevalence$time <- as.numeric(data.table::tstrsplit(postdraws_prevalence$variable,"\\.")[[2]])
  postdraws_prevalence$theta_type <- as.numeric(data.table::tstrsplit(postdraws_prevalence$variable,"\\.")[[3]])
  
  prevalence_gammas = subset(postdraws_prevalence, type == "gamma")
  prevalence_betas  = subset(postdraws_prevalence, type == "beta")
  
  ## just keep ypred and theta_type == 2
  postdraws_prevalence <- subset(postdraws_prevalence, (theta_type == 2 & time <= length_of_time_series + 1) | (type == "ypred" & time > length_of_time_series))
  postdraws_prevalence <- subset(postdraws_prevalence, select=c("chain","id","type","time","value"))
  ## theta time is off shifted by 1
  postdraws_prevalence$theta_time <- postdraws_prevalence$time - 1
  postdraws_prevalence$unqid <- paste0(postdraws_prevalence$chain,"_",postdraws_prevalence$id)
  
  ## compute summaries
  postdrawsply_prevalence <- ddply(postdraws_prevalence,.(time, type),summarise,lower = quantile(value,probs=.025), avg = mean(value), upper = quantile(value,probs=.975))
  postdrawsply_prevalence$theta_time <- postdrawsply_prevalence$time - 1
  
  
  
  ######################################################################
  ## plot a version of Fig 6
  ## plot a version of Fig 6
  xx1 = subset(postdrawsply_incidence, type == "ypred")$upper
  xx2 = subset(postdrawsply_prevalence, type == "ypred")$upper
  max_y = max(xx1, xx2, df_forecast_year[df_forecast_year$epi_time %in% 1:35,]$iliplus)
  
  incidence_pred_plot = ggplot()+
    geom_ribbon(aes(x=theta_time, ymin=lower, ymax=upper), data=subset(postdrawsply_incidence, type == "theta"),fill=I("darkgrey"))+
    geom_ribbon(aes(x=time, ymin=lower, ymax=upper), data=subset(postdrawsply_incidence, type == "ypred"),fill=I("lightgrey"))+
    geom_line(aes(x=theta_time, y=avg), data=subset(postdrawsply_incidence, type == "theta"))+
    geom_line(aes(x=time, y=avg), data=subset(postdrawsply_incidence, type == "ypred"))+
    geom_point(aes(x=1:35, y=df_forecast_year[df_forecast_year$epi_time %in% 1:35,]$iliplus), color=I("black"))+
    ylab("Proportion Infectious")+
    xlab("Weeks") + ggtitle("DBSSM Model with Incidence Maps")+ 
    ylim(c(0, max_y+0.005))+ 
    theme(axis.text=element_text(size=18), axis.title = element_text(size = 18), plot.title = element_text(size = 18)) + 
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1))
  
  
  prevalence_pred_plot = ggplot()+
    geom_ribbon(aes(x=theta_time, ymin=lower, ymax=upper), data=subset(postdrawsply_prevalence, type == "theta"),fill=I("darkgrey"))+
    geom_ribbon(aes(x=time, ymin=lower, ymax=upper), data=subset(postdrawsply_prevalence, type == "ypred"),fill=I("lightgrey"))+
    geom_line(aes(x=theta_time, y=avg), data=subset(postdrawsply_prevalence, type == "theta"))+
    geom_line(aes(x=time, y=avg), data=subset(postdrawsply_prevalence, type == "ypred"))+
    geom_point(aes(x=1:35, y=df_forecast_year[df_forecast_year$epi_time %in% 1:35,]$iliplus), color=I("black"))+
    ylab("Proportion Infectious") +
    xlab("Weeks") + ggtitle("DBSSM Model with Prevalence Maps") + 
    ylim(c(0, max_y+0.005))+ 
    theme(axis.text=element_text(size=18), axis.title = element_text(size = 18), plot.title = element_text(size = 18)) + 
    theme(axis.text.x = element_text(angle = 0, vjust = 0.5, hjust=1))
  
  ##########################################################################
  ##########################################################################å
  
  # save the plots
  both_plots <- list(prevalence_pred_plot, incidence_pred_plot)
  grid.arrange(both_plots[[1]], both_plots[[2]],
               layout_matrix = matrix(c(1,2), byrow = TRUE, ncol = 1))
  saveRDS(both_plots, file = paste0("../figures/", flu_season, "_", length_of_time_series, ".rds"))
}

###########
# # Let's compare the dave_beta, dave_gamma values drawn via each method.
# rho_data                = data.frame(variable = rep('rho', times = nrow(incidence_betas)), 
#                                      value = incidence_betas[c(4)]/incidence_gammas[c(4)])
# incidence_data          = rbind(incidence_gammas[c(3,4)],incidence_betas[c(3,4)],rho_data)
# incidence_data$approach = rep('incidence', times = nrow(incidence_data))
# 
# rho_data                 = data.frame(variable = rep('rho', times = nrow(prevalence_betas)), 
#                                       value = prevalence_betas[c(4)]/prevalence_gammas[c(4)])
# prevalence_data          = rbind(prevalence_gammas[c(3,4)],prevalence_betas[c(3,4)], rho_data)
# prevalence_data$approach = rep('prevalence', times = nrow(prevalence_data))
# 
# graph_data = rbind(incidence_data, prevalence_data)
# graph_data$variable = as.factor(graph_data$variable)
# 
# ggplot(graph_data, aes(x = variable, y = value, fill = approach)) + 
#   geom_boxplot() + scale_fill_manual(values=c("#999999", "#FFFFFF")) + 
#   scale_x_discrete(labels = TeX(c('beta' = r"($\beta$)", 
#                                   'gamma' = r"($\gamma$)", 
#                                   'rho' = r"($\rho)")) )+ 
#   labs(fill="Approach")+ 
#   theme(axis.text=element_text(size=20), axis.title = element_text(size = 20),
#         legend.title=element_text(size=20), 
#         legend.text=element_text(size=15)) + 
#   ylab("")+
#   xlab("Parameter")
# 
# graph_data %>%
#   subset(variable == 'gamma' & approach == 'incidence') %>%
#   pull(value) %>%
#   median()
# graph_data %>%
#   subset(variable == 'gamma' & approach == 'prevalence') %>%
#   pull(value) %>%
#   median()
# 
# 
# graph_data %>%
#   subset(variable == 'beta' & approach == 'incidence') %>%
#   pull(value) %>%
#   median()
# graph_data %>%
#   subset(variable == 'beta' & approach == 'prevalence') %>%
#   pull(value) %>%
#   median()
# 
# 
# graph_data %>%
#   subset(variable == 'rho' & approach == 'incidence') %>%
#   pull(value) %>%
#   median()
# graph_data %>%
#   subset(variable == 'rho' & approach == 'prevalence') %>%
#   pull(value) %>%
#   median()


#### Compile all graphics into single graphics.
# Iterate through each flu season
for (season in unique(flu_seasons)) {
  plot_list <- list()  # To store all plots for this flu season
  
  # For each length of time series, load the plots and store them
  for (length_ts in 1:2) {
    season_tmp = season
    shift = 0
    if(season == 2014) shift = 2
    length_ts = lengths[length_ts + shift]
    rds_path <- paste0("../figures/", season, "_", length_ts, ".rds")
    
    if (file.exists(rds_path)) {
      both_plots <- readRDS(rds_path)
      
      # Append both plots to the list
      plot_list <- c(plot_list, both_plots)
    } else {
      warning(paste("Missing file:", rds_path))
    }
  }
  
  # Create output PNG
  png_filename <- paste0("../figures/", season, "_alternate_stacked.png")
  png(png_filename, width = 1200, height = 700)
  
  # Number of plots and desired columns
  num_plots <- length(plot_list)
  n_cols <- 2
  n_rows <- ceiling(num_plots / n_cols)
  
  # Generate layout matrix that fills by column
  layout_matrix <- matrix(1:(n_cols * n_rows), nrow = n_rows, ncol = n_cols, byrow = FALSE)
  
  # Trim matrix to actual number of plots (in case it's over)
  layout_matrix[layout_matrix > num_plots] <- NA
  
  # Use layout_matrix to fill by column
  do.call(grid.arrange, c(plot_list, list(layout_matrix = layout_matrix)))
  
  dev.off()
}
