## This code was originally written by Dave Osthus using the following license:

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

## plot rhostar vs rho
qplot(rhostargrid, rhogrid)+
  ggtitle("Rho vs Rho*")+
  xlab("Rho* = PI - I0 - S0")+
  ylab("Rho")


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
length_of_time_series <- 22
myy <- df2010[df2010$epi_time %in% 1:length_of_time_series,]$iliplus

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

## define number of chains to run
nchains <- 2

## run JAGS
results <- run.jags(model = paste0(jagspath,"dbssm_jags.txt"),
                    data = jags_data,
                    monitor = c("ypred","theta"),
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
postdraws$unqid <- paste0(postdraws$chain,"_",postdraws$id)

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
  

