## Dave Osthus
## 2-1-24
## Make synthetic training data

## load libraries
library(ggplot2)
library(plyr)
library(data.table)
library(LearnBayes)
library(LaplacesDemon)
theme_set(theme_bw())

## set proxies
Sys.setenv('https_proxy'='http://proxyout.lanl.gov:8080')## not sure why needed, but see here: https://github.com/curl/curl/issues/1015

## define SIRfunctions
sir <- function(beta, gamma, S0, I0, R0, times) {
  require(deSolve) # for the "ode" function
  
  # the differential equations:
  sir_equations <- function(time, variables, parameters) {
    with(as.list(c(variables, parameters)), {
      dS <- -beta * I * S
      dI <-  beta * I * S - gamma * I
      dR <-  gamma * I
      return(list(c(dS, dI, dR)))
    })
  }
  
  # the parameters values:
  parameters_values <- c(beta  = beta, gamma = gamma)
  
  # the initial values of variables:
  initial_values <- c(S = S0, I = I0, R = R0)
  
  # solving
  out <- ode(initial_values, times, sir_equations, parameters_values)
  
  # returning the output:
  as.data.frame(cbind(out, beta = beta, gamma = gamma))
}

gen_curve = function(curve_type){
  library(plyr)
  library(data.table)
  library(LearnBayes)
  library(LaplacesDemon)

  #########################
  ### SIR-ROLLERCOASTER ###
  #########################
  if(curve_type == "sir_rollercoaster"){
    seas_bool = F #not seasonal
    
    ## draw number of waves
    nwaves <- sample(3:7,1)
    #nwaves <- sample(3:15,1)
    ## draw reproduction numbers with preference toward lower values
    #basic_repo <- exp(runif(nwaves,log(1),log(35))) 
    basic_repo <- 1.1+rbeta(nwaves, shape1 = 1, shape2 = 3)*20 #but, have a preference for lower values
    # hist(1.1+rbeta(1000, shape1 = 1, shape2 = 3)*20, breaks = 50)
    
    # x = seq(0,10,0.01)
    # plot(x,dgamma(x,shape=1.5,scale=1.5))
    # hist(exp(pmin(rgamma(1000,shape=1,scale=2),log(35))), breaks = 100)
    # 
    ### randomly choose if wave peak heights generally increasing, decreasing, or random
    order_type = sample(c('increasing', 'decreasing', 'random'), size = 1)
    if(order_type == 'increasing'){
      basic_repo_order <- rev(sample(1:nwaves,nwaves,replace=F,prob = basic_repo))
      basic_repo <- basic_repo[basic_repo_order]
    }else if(order_type == 'decreasing'){
      basic_repo_order <- sample(1:nwaves,nwaves,replace=F,prob = basic_repo)
      basic_repo <- basic_repo[basic_repo_order]
    }
    
    ## draw gamma values
    invgamma <- runif(nwaves, 1, 10)
    gamma <- 1/invgamma
    init_cond <- rdirichlet(nwaves,c(1000,0.1,0.1)) #everyone starts as susceptible
    starttime <- 1
    
    ## waves become less frequent
    cadence <- sample(10:52,nwaves) 
    cadence_order <- sample(1:nwaves,nwaves,replace=F,prob = max(cadence) + 1 - cadence)
    cadence <- cadence[cadence_order]
    tslist <- list()
    for(jj in 1:nwaves){
      S0 <- init_cond[jj,1]
      I0 <- init_cond[jj,2]
      R0 <- init_cond[jj,3]
      beta <- (basic_repo[jj]/S0)*gamma[jj]
      times <- 1:500
      ts <- sir(beta, gamma[jj], S0, I0, R0, times)$I
      cnt = 1
      while(sum(is.nan(ts))!=0){
        basic_repo[jj] <- 1.1+rbeta(nwaves, shape1 = 1, shape2 = 3)*20 #but, have a preference for lower values
        invgamma <- runif(1, 1, 10)
        gamma[jj] <- 1/invgamma
        # init_cond <- rdirichlet(1,c(1000,1,1000))
        init_cond <- rdirichlet(1,c(1000,0.1,0.1))
        S0 <- init_cond[1]
        I0 <- init_cond[2]
        R0 <- init_cond[3]
        beta <- (basic_repo/S0)*gamma
        times <- 1:500
        ts <- sir(beta, gamma, S0, I0, R0, times)$I
        cnt = cnt+1
        if(cnt == 5){
          ts = rep(0.01, length(ts))
        }
      }
      max_cut = max(which(ts/max(ts) > .0001))
      if(is.infinite(max_cut)){
        max_cut = length(ts)
      }
      tslist[[jj]] <- c(rep(0,starttime),ts[1:max_cut])
      ## pad the beginning and end with 0s
      tslist[[jj]] <- c(runif(20,0,quantile(tslist[[jj]],prob=.01, na.rm=T)+1e-9),
                        tslist[[jj]],
                        runif(20,0,quantile(tslist[[jj]],prob=.01, na.rm=T)+1e-9))
      tslist[[jj]][is.nan(tslist[[jj]])] = 0
      starttime <- starttime + cadence[jj]
    }
    tslength <- max(unlist(lapply(tslist,length)))
    ts <- rep(0,tslength)
    for(jj in 1:nwaves){
      tslist[[jj]] <- c(tslist[[jj]],rep(0,tslength - length(tslist[[jj]])))  
      ts <- ts + tslist[[jj]]
    }
    ts = ts/max(ts)
    
    ### drop extra leading zeros
    min_cut = min(which(ts>1e-3)) - 10
    if(min_cut > 1 & !is.infinite(min_cut) & min_cut < (length(ts)-10)){
      ts = ts[-(1:min_cut)]
    }
    max_cut = max(which(ts>1e-3)) + 10
    if(max_cut < length(ts)& !is.infinite(max_cut) & max_cut > 10){
      ts = ts[-(max_cut:length(ts))]
    }
  }
  
  ################################
  ### SIR-ROLLERCOASTER-WIGGLE ###
  ################################
  if(curve_type == "sir_rollercoaster_wiggle"){
    seas_bool = F #not seasonal
    
    ## draw number of waves
    nwaves <- sample(3:7,1)
    #nwaves <- sample(3:15,1)
    
    ## draw reproduction numbers with preference toward lower values
    #basic_repo <- exp(runif(nwaves,log(1),log(35))) 
    basic_repo <- 1.1+rbeta(nwaves, shape1 = 1, shape2 = 3)*20 #but, have a preference for lower values
    
    ### randomly choose if wave peak heights generally increasing, decreasing, or random
    order_type = sample(c('increasing', 'decreasing', 'random'), size = 1)
    if(order_type == 'increasing'){
      basic_repo_order <- rev(sample(1:nwaves,nwaves,replace=F,prob = basic_repo))
      basic_repo <- basic_repo[basic_repo_order]
    }else if(order_type == 'decreasing'){
      basic_repo_order <- sample(1:nwaves,nwaves,replace=F,prob = basic_repo)
      basic_repo <- basic_repo[basic_repo_order]
    }

    ## draw gamma values
    invgamma <- runif(nwaves, 1, 10)
    gamma <- 1/invgamma
    init_cond <- rdirichlet(nwaves,c(1000,0.1,0.1)) #everyone starts as susceptible
    starttime <- 1
    
    ## waves become less frequent
    cadence <- sample(10:52,nwaves) 
    cadence_order <- sample(1:nwaves,nwaves,replace=F,prob = max(cadence) + 1 - cadence)
    cadence <- cadence[cadence_order]
    tslist <- list()
    for(jj in 1:nwaves){
      S0 <- init_cond[jj,1]
      I0 <- init_cond[jj,2]
      R0 <- init_cond[jj,3]
      beta <- (basic_repo[jj]/S0)*gamma[jj]
      times <- 1:500
      ts <- sir(beta, gamma[jj], S0, I0, R0, times)$I
      cnt = 1
      while(sum(is.nan(ts))!=0){
        #basic_repo[jj] <- exp(pmax(0.15,pmin(rgamma(1,shape = 1, rate = 2), log(35)))) #but, have a preference for lower values
        basic_repo[jj] <- 1.1+rbeta(nwaves, shape1 = 1, shape2 = 3)*20
        invgamma <- runif(1, 1, 10)
        gamma[jj] <- 1/invgamma
        # init_cond <- rdirichlet(1,c(1000,1,1000))
        init_cond <- rdirichlet(1,c(1000,0.1,0.1))
        S0 <- init_cond[1]
        I0 <- init_cond[2]
        R0 <- init_cond[3]
        beta <- (basic_repo/S0)*gamma
        times <- 1:500
        ts <- sir(beta, gamma, S0, I0, R0, times)$I
        cnt = cnt+1
        if(cnt == 5){
          ts = rep(0.01, length(ts))
        }
      }
      max_cut = max(which(ts/max(ts) > .0001))
      if(is.infinite(max_cut)){
        max_cut = length(ts)
      }
      tslist[[jj]] <- c(rep(0,starttime),ts[1:max_cut])
      
      ## pad the beginning and end with 0s
      tslist[[jj]] <- c(runif(20,0,quantile(tslist[[jj]],prob=.01, na.rm=T)+1e-9),
                        tslist[[jj]],
                        runif(20,0,quantile(tslist[[jj]],prob=.01, na.rm=T)+1e-9))
      tslist[[jj]][is.nan(tslist[[jj]])] = 0
      starttime <- starttime + cadence[jj]
    }
    tslength <- max(unlist(lapply(tslist,length)))
    ts <- rep(0,tslength)
    for(jj in 1:nwaves){
      tslist[[jj]] <- c(tslist[[jj]],rep(0,tslength - length(tslist[[jj]])))  
      ts <- ts + tslist[[jj]]
    }
    ts = ts/max(ts)
    
    ### drop extra leading zeros
    min_cut = min(which(ts>1e-3)) - 10
    if(min_cut > 1 & !is.infinite(min_cut) & min_cut < (length(ts)-10)){
      ts = ts[-(1:min_cut)]
    }
    max_cut = max(which(ts>1e-3)) + 10
    if(max_cut < length(ts)& !is.infinite(max_cut) & max_cut > 10){
      ts = ts[-(max_cut:length(ts))]
    }
    
    mult_squish <- runif(1,1,2)
    mult_period <- length(ts)/runif(1,sqrt(2),sqrt(10))^2
    mult_sin <- 1+mult_squish*(1+sample(c(-1,1),1)*sin((pi*(1:length(ts)))/mult_period))
    ts <- pmax(0,mult_sin*ts)
    ts = ts/max(ts)
  }
  
  ################
  ### SEASONAL ###
  ################
  if(curve_type == "seasonal"){
    seas_bool = T
    ## get parameters to draw SIR curves
    #basic_repo <- exp(runif(1,log(2),log(5)))
    nwaves <- sample(5:15,1)
    
    # basic_repo <- exp(pmax(0.7,pmin(rgamma(1,shape = 1, rate = 2), log(5)))) #but, have a preference for lower values
    basic_repo <- 1.2+rbeta(1, shape1 = 1, shape2 = 3)*10 #but, have a preference for lower values

    

    
    invgamma <- runif(1, 1, 10)
    gamma <- 1/invgamma
    # init_cond <- rdirichlet(1,c(1000,1,1000))
    init_cond <- rdirichlet(1,c(1000,0.1,0.1))
    S0 <- init_cond[1]
    I0 <- init_cond[2]
    R0 <- init_cond[3]
    beta <- (basic_repo/S0)*gamma
    times <- 1:500
    ts <- sir(beta, gamma, S0, I0, R0, times)$I
    cnt = 1
    while(sum(is.nan(ts))!=0){
      # basic_repo <- exp(pmax(0.7,pmin(rgamma(1,shape = 1, rate = 2), log(5)))) #but, have a preference for lower values
      basic_repo <- 1.2+rbeta(1, shape1 = 1, shape2 = 3)*10 #but, have a preference for lower values
      invgamma <- runif(1, 1, 10)
      gamma <- 1/invgamma
      # init_cond <- rdirichlet(1,c(1000,1,1000))
      init_cond <- rdirichlet(1,c(1000,0.1,0.1))
      S0 <- init_cond[1]
      I0 <- init_cond[2]
      R0 <- init_cond[3]
      beta <- (basic_repo/S0)*gamma
      times <- 1:500
      ts <- sir(beta, gamma, S0, I0, R0, times)$I
      cnt = cnt+1
      if(cnt == 5){
        ts = rep(0.01, length(ts))
      }
    }
    # max_cut = max(which(ts/max(ts) > .00001))
    # if(is.infinite(max_cut)){
    #   max_cut = length(ts)
    # }
    # ts <- ts[1:max_cut]
    

    expit = function(x){exp(x)/(1+exp(x))}
    trend = sample(1:2, size = 1, prob = c(0.7, 0.3))
    if(trend == 1){
      prop_scale = expit(rnorm(1,mean=0,sd = 0.5)*1:nwaves)
      prop_scale = prop_scale/prop_scale[1]
    }else{
      prop_scale = rep(1,nwaves)
    }
    
    ts_long = rep(0, nwaves*pmax(52, length(ts)))
    for(jj in 1:nwaves){
      ts_long[(1:length(ts)) + 52*(jj-1)] = ts_long[(1:length(ts)) + 52*(jj-1)] + prop_scale[jj]*ts
    }
    ts_long[is.na(ts_long)] = 0
    ts = ts_long
    ts = ts[1:floor(52*nwaves)]
    if(max(ts)>1){
      ts = ts/max(ts)
    }
    ts = ts[-c(1:52)] #get rid of start-of-ts effects
    ts = ts[-c((length(ts)-26):length(ts))] #get rid of end-of-ts effects
  }
  
  ##################################
  ### Additional Transformations ###
  ##################################
  
  
  ## pad the beginning and end with 0s
  if(curve_type != "seasonal"){
    ts <- c(runif(20,0,quantile(ts,prob=.011,na.rm=T)+1e-9),
            ts,
            runif(20,0,quantile(ts,prob=.011,na.rm=T)+1e-9))
    
    ## randomly add spike
    # spike = sample(c('no','yes'),1)
    # if(spike == 'yes'){
    #   ts[sample(1:length(ts),1)] = runif(1,0,1)
    # }
  }

  
  time_cadence = 'weekly'
  if(curve_type == 'seasonal'){
    time_cadence <- c("monthly","weekly")[sample(1:2,1)]
    if(time_cadence == "monthly"){
      index = floor(1:length(ts)/4)
      ts <- aggregate(ts~index, FUN = sum)$ts
      if(max(ts)>1){
        ts = ts/max(ts)
      }
    }
  }
  # else{
  #   ts <- approx(x=seq(0,floor(52*length(ts)/365),length.out = length(ts)),
  #                y=ts,
  #                xout=seq(0,floor(52*length(ts)/365),1))$y
  # }
  # 
  # 
  ## choose scale
  ts_scale <- c("proportion","counts")[sample(1:2,1)]
  
  
  ## Add Noise and Rescale

    ## add noise
    alpha <- exp(runif(1,log(5e1),log(1e4)))
    obs_ts <- rbeta(1:length(ts), alpha*ts, alpha*(1-ts))
    
    
    ## put the peak on a reasonable scale
    PI <- exp(runif(1,log(.0005),log(.25)))
    obs_ts <- (obs_ts/max(obs_ts))*PI
    
    ## upscale
    if(ts_scale == "counts"){
      Npop <- round(runif(1,sqrt(2e5),sqrt(1e8))^2,0)
      obs_ts <- round(Npop*obs_ts,0)
    }

  
  ## if peak happens too soon, pad the start of the time series with 0s
  if(which.max(obs_ts) < 30 & curve_type != 'seasonal'){
    obs_ts <- c(runif(sample(10:17,1),0,quantile(obs_ts,prob = .05)+1e-9),obs_ts)
  }
    

  ## don't allow anything to be less than 1e-10
  obs_ts <- pmax(1e-10, obs_ts)
  
  
  ## append to list if no error occurred
  if(!is.na(sum(obs_ts))){
    templist <- list(ts = obs_ts,
                     ts_dates = NULL,
                     ts_exogenous = NULL,
                     ts_real_data = F,
                     ts_isolated_strain = ifelse(curve_type == "sir_rollercoaster",F,T),
                     ts_multiwave = ifelse(curve_type == "sir_rollercoaster",T,F),
                     ts_disease = curve_type, 
                     ts_measurement_type = NA,
                     ts_geography = NA,
                     ts_first_time = NA,
                     ts_last_time = NA,
                     ts_time_cadence = time_cadence,
                     ts_scale = ts_scale,
                     ts_exogenous_scale = NA, 
                     ts_seasonal = seas_bool)
  }
  return(templist)
}



