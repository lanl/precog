## Dave Osthus
## 1-18-24
## Recreating figs and code from Osthus et. al. 2010 AoAS paper
## "Forecasting Seasonal Influenza With A State-Space SIR Model"

## load libraries
setwd("~/GitLab/incidence_for_sir/data_application/code")
library(rjags)
library(runjags)
library(ggplot2)
library(reshape2)
library(data.table)
library(plyr)
theme_set(theme_bw())

set.seed(13)

## define paths
datapath <- "~/GitLab/incidence_for_sir/data_application/data/"
jagspath <- "~/GitLab/incidence_for_sir/data_application/code/"

## load ILI and surveillance data
dfili <- read.csv(paste0(datapath,"EW08-2020_nat.csv"))
dfili <- subset(dfili, region == "nat", select=c("epi_year","epi_week","wili","epi_time","epi_season"))
dfsurveillance <- read.csv(paste0(datapath,"EW08-2020_who_combined_prior_to_2015_16.csv"))
dfsurveillance <- subset(dfsurveillance, select=c("year","week","percent_positive"), region == "National")

## merge ILI data with surveillance data
df <- merge(dfili, dfsurveillance, by.x=c("epi_year","epi_week"), by.y=c("year","week"))
df$iliplus <- (df$wili/100)*(df$percent_positive/100)

## Reproduce the Top of Fig 1
qplot(epi_time, iliplus, data=subset(df, epi_time <= 35 & epi_season %in% c(2002:2007,2010:2013)), geom="line", group=epi_season)+
  ylim(c(0,.025))+
  ylab("ILI+")+
  xlab("Week")+
  scale_x_continuous(breaks = seq(0,100,5))

## focus on 2010, national wILI (Fig 6 of Osthus et. al. AoAS)
df2010 <- subset(df, epi_season == 2010 & epi_time <= 35)
df2010 <- df2010[order(df2010$epi_time),]

## Reproduce Fig 6 
qplot(epi_time, iliplus, data=df2010, geom=c("line","point"), group=epi_season)+
  ylab("Proportion Infectious")+
  xlab("t")


###################################
## Prepare JAGS
## Equation 6.7
params_eq6.7 <- c(1.62, 7084.10)

## Equation 6.10
mean_eq6.10 <- c(0.0144, 17.9)
Sigma_eq6.10 <- matrix(c(0.000036, -0.0187, -0.0187, 16.09), nrow=2)
# 
# ## Make grid for g^-1(.) in equation 6.12
# # Murph note: if you update this grid at all, you must also update the scaling
# # factors inside of JAGS.
# inputgrid <- expand.grid(S0  = 0.9,
#                          I0  = seq(0,.05,length.out = 100),
#                          PIV = seq(0,.1,length.out = 100),
#                          PIT = seq(1,35,length.out = 100))
# # inputgrid <- subset(inputgrid, PIV >= I0) # I know this is ugly, but I'm removing this atm.
# inputgrid <- subset(inputgrid, select=c("S0","I0","PIV","PIT"))

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
# pretty cool.  I'm also using linear interpolation, but over three variables
# trilinear interpolation.

#########################################
## Fit Model in JAGS

## the observed time series
length_of_time_series <- 22
myy <- df2010[df2010$epi_time %in% 1:length_of_time_series,]$iliplus

## construct data to be passed into JAGS
jags_data <- list(y = myy,
                  tt = length(myy),
                  TT = 35,
                  params_eq6.7 = params_eq6.7,
                  mean_eq6.10 = mean_eq6.10,
                  Sigma_eq6.10 = Sigma_eq6.10,
                  dave_betas = inputgrid$dave_betas,
                  dave_gammas = inputgrid$dave_gamma)

## define number of chains to run
nchains <- 2

## run JAGS
results <- run.jags(model = paste0(jagspath,"dbssm_jags_incidence.txt"),
                    data = jags_data,
                    monitor = c("ypred","theta","beta","gamma"),
                    n.chains = nchains,
                    burnin = 12500,
                    sample = 5000,
                    thin = 10)

## organize the posterior draws from JAGS
postdraws <- NULL
for(i in 1:nchains){
  temppostdraws <- reshape2::melt(data.frame(chain = i, id = 1:nrow(data.frame(results$mcmc[[i]])), data.frame(results$mcmc[[i]])),id.vars = c("chain","id"))
  postdraws <- rbind(postdraws, temppostdraws)
}

## organize output
postdraws$type <- data.table::tstrsplit(postdraws$variable,"\\.")[[1]]
postdraws$time <- as.numeric(data.table::tstrsplit(postdraws$variable,"\\.")[[2]])
postdraws$theta_type <- as.numeric(data.table::tstrsplit(postdraws$variable,"\\.")[[3]])

## just keep ypred and theta_type == 2
postdraws <- subset(postdraws, (theta_type == 2 & time <= length_of_time_series + 1) | (type == "ypred" & time > length_of_time_series))
postdraws <- subset(postdraws, select=c("chain","id","type","time","value"))

## theta time is off shifted by 1
postdraws$theta_time <- postdraws$time - 1
postdraws$unqid      <- paste0(postdraws$chain,"_",postdraws$id)

## compute summaries
postdrawsply            <- ddply(postdraws,.(time, type),summarise,lower = quantile(value,probs=.025), avg = mean(value), upper = quantile(value,probs=.975))
postdrawsply$theta_time <- postdrawsply$time - 1

## plot a version of Fig 6
ggplot()+
  geom_ribbon(aes(x=theta_time, ymin=lower, ymax=upper), data=subset(postdrawsply, type == "theta"),fill=I("darkgrey"))+
  geom_ribbon(aes(x=time, ymin=lower, ymax=upper), data=subset(postdrawsply, type == "ypred"),fill=I("lightgrey"))+
  geom_line(aes(x=theta_time, y=avg), data=subset(postdrawsply, type == "theta"))+
  geom_line(aes(x=time, y=avg), data=subset(postdrawsply, type == "ypred"))+
  geom_point(aes(x=1:35, y=df2010[df2010$epi_time %in% 1:35,]$iliplus), color=I("black"))+
  ylab("Proportion Infectious")+
  xlab("t")
  

