# depends: regime_functions
#' @export
sir                               <- function(beta, gamma, S0, I0, R0, times) {
  # the differential equations:
  sir_equations                   <- function(time, variables, parameters) {
    with(as.list(c(variables, parameters)), {
      dS                          <- -beta * I * S
      dI                          <- beta * I * S - gamma * I
      dR                          <- gamma * I
      return(list(c(dS, dI, dR)))
    })
  }

  # the parameters values:
  parameters_values               <- c(beta  = beta, gamma = gamma)
  
  # the initial values of variables:
  initial_values                  <- c(S = S0, I = I0, R = R0)
  
  # solving
  # out                             <- deSolve::ode(initial_values, times, sir_equations, parameters_values)
  out <- local({
    zz <- file(nullfile(), "wt"); sink(zz); sink(zz, type="message")
    on.exit({sink(type="message"); sink(); close(zz)}, add=TRUE)

    suppressWarnings(
      deSolve::ode(initial_values, times, sir_equations, parameters_values)
    )
  })

  
  # returning the output:
  as.data.frame(cbind(out, incidence = c(NA,-diff(as.data.frame(out)$S)), beta = beta, gamma = gamma))
}

#' @export
gen_curve = function(nwaves,
                        invgamma1,
                        invgamma2,
                        repo1,
                        repo2,
                        alpha,
                        Npop,
                        PI,
                        curve_type){
  
  nwaves = floor(nwaves) + 1
  invgamma                      <- 1 + 9*rbeta(nwaves, shape1 = invgamma1, shape2 = invgamma2)
  basic_repo                    <- 1.1 + rbeta(nwaves, shape1 = repo1, shape2 = repo2)*20
  Npop = round(Npop,0)
  
  #########################
  ### SIR-ROLLERCOASTER ###
  #########################
  PITs                            <- c()
  PIVs                            <- c()
  s0s                             <- c()
  i0s                             <- c()
  r0s                             <- c()
  if(curve_type == "sir_rollercoaster"){
    seas_bool = F #not seasonal

    order_type = sample(c('increasing', 'decreasing', 'random'), size = 1)
    if(order_type == 'increasing'){
      basic_repo_order            <- rev(sample(1:nwaves,nwaves,replace=F,prob = basic_repo))
      basic_repo                  <- basic_repo[basic_repo_order]
    }else if(order_type == 'decreasing'){
      basic_repo_order            <- sample(1:nwaves,nwaves,replace=F,prob = basic_repo)
      basic_repo                  <- basic_repo[basic_repo_order]
    }
    
    gamma                         <- 1/invgamma
    init_cond                     <- MCMCpack::rdirichlet(nwaves,c(1000,0.1,0.1)) #everyone starts as susceptible
    starttime                     <- 1
    
    ## waves become less frequent
    cadence                       <- sample(10:52,nwaves)
    cadence_order                 <- sample(1:nwaves,nwaves,replace=F,prob = max(cadence) + 1 - cadence)
    cadence                       <- cadence[cadence_order]
    tslist                        <- list()
    for(jj in 1:nwaves){
      S0                          <- init_cond[jj,1]
      I0                          <- init_cond[jj,2]
      R0                          <- init_cond[jj,3]
      beta                        <- (basic_repo[jj]/S0)*gamma[jj]
      times                       <- 1:500
      sir_output                  <- sir(beta, gamma[jj], S0, I0, R0, times)
      ts                          <- sir_output$I
      cnt = 1
      while(sum(is.nan(ts))!=0){
        init_cond                 <- MCMCpack::rdirichlet(1,c(1000,0.1,0.1))
        S0                        <- init_cond[1]
        I0                        <- init_cond[2]
        R0                        <- init_cond[3]
        beta                      <- (basic_repo/S0)*gamma
        times                     <- 1:500
	      sir_output                <- sir(beta, gamma[jj], S0, I0, R0, times)
        ts                        <- sir_output$I
        cnt = cnt+1
        if(cnt == 5){
          stop("failed to generate an SIR curve with those parameters.")
        }
      }
      max_cut = max(which(ts/max(ts) > .0001))
      if(is.infinite(max_cut)){
        max_cut = length(ts)
      }
      tslist[[jj]]                <- c(rep(0,starttime),ts[1:max_cut])
      ## pad the beginning and end with 0s
      tslist[[jj]] <- suppressWarnings(
                        c(
                          runif(20, 0, quantile(tslist[[jj]], prob = 0.01, na.rm = TRUE) + 1e-9),
                          tslist[[jj]],
                          runif(20, 0, quantile(tslist[[jj]], prob = 0.01, na.rm = TRUE) + 1e-9)
                        )
                      )
      tslist[[jj]][is.nan(tslist[[jj]])] = 0
      starttime                   <- starttime + cadence[jj]

      temp_incidences             <- sir_output$incidence
      temp_incidences[1]          <- 0
      max_incidence_time          <- which.max(temp_incidences) - 1
      max_incidence_value         <- max(temp_incidences)
      PITs                        <- c(PITs, max_incidence_time)
      PIVs                        <- c(PIVs, max_incidence_value)
      s0s                         <- c(s0s, S0)
      i0s                         <- c(i0s, I0)
      r0s                         <- c(r0s, R0)

    }
    tslength                      <- max(unlist(lapply(tslist,length)))
    ts                            <- rep(0,tslength)
    for(jj in 1:nwaves){
      tslist[[jj]]                <- c(tslist[[jj]],rep(0,tslength - length(tslist[[jj]])))
      ts                          <- ts + tslist[[jj]]
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
    
    ### randomly choose if wave peak heights generally increasing, decreasing, or random
    order_type = sample(c('increasing', 'decreasing', 'random'), size = 1)
    if(order_type == 'increasing'){
      basic_repo_order            <- rev(sample(1:nwaves,nwaves,replace=F,prob = basic_repo))
      basic_repo                  <- basic_repo[basic_repo_order]
    }else if(order_type == 'decreasing'){
      basic_repo_order            <- sample(1:nwaves,nwaves,replace=F,prob = basic_repo)
      basic_repo                  <- basic_repo[basic_repo_order]
    }

    gamma                         <- 1/invgamma
    init_cond                     <- MCMCpack::rdirichlet(nwaves,c(1000,0.1,0.1)) #everyone starts as susceptible
    starttime                     <- 1
    
    ## waves become less frequent
    cadence                       <- sample(10:52,nwaves)
    cadence_order                 <- sample(1:nwaves,nwaves,replace=F,prob = max(cadence) + 1 - cadence)
    cadence                       <- cadence[cadence_order]
    tslist                        <- list()
    for(jj in 1:nwaves){
      S0                          <- init_cond[jj,1]
      I0                          <- init_cond[jj,2]
      R0                          <- init_cond[jj,3]
      beta                        <- (basic_repo[jj]/S0)*gamma[jj]
      times                       <- 1:500
      sir_output                  <- sir(beta, gamma[jj], S0, I0, R0, times)
      ts                          <- sir_output$I
      cnt = 1
      while(sum(is.nan(ts))!=0){
        init_cond                 <- MCMCpack::rdirichlet(1,c(1000,0.1,0.1))
        S0                        <- init_cond[1]
        I0                        <- init_cond[2]
        R0                        <- init_cond[3]
        beta                      <- (basic_repo/S0)*gamma
        times                     <- 1:500
	    sir_output                <- sir(beta, gamma[jj], S0, I0, R0, times)
        ts                        <- sir_output$I
        cnt = cnt+1
        if(cnt == 5){
          stop("failed to generate SIR curve")
        }
      }
      max_cut = max(which(ts/max(ts) > .0001))
      if(is.infinite(max_cut)){
        max_cut = length(ts)
      }
      tslist[[jj]]                <- c(rep(0,starttime),ts[1:max_cut])
      
      ## pad the beginning and end with 0s
      tslist[[jj]] <- suppressWarnings(
        c(
          runif(20, 0, quantile(tslist[[jj]], prob = 0.01, na.rm = TRUE) + 1e-9),
          tslist[[jj]],
          runif(20, 0, quantile(tslist[[jj]], prob = 0.01, na.rm = TRUE) + 1e-9)
        )
      )
      tslist[[jj]][is.nan(tslist[[jj]])] = 0
      starttime                   <- starttime + cadence[jj]
      temp_incidences             <- sir_output$incidence
      temp_incidences[1]          <- 0
      max_incidence_time          <- which.max(temp_incidences) - 1
      max_incidence_value         <- max(temp_incidences)
      PITs                        <- c(PITs, max_incidence_time)
      PIVs                        <- c(PIVs, max_incidence_value)
      s0s                         <- c(s0s, S0)
      i0s                         <- c(i0s, I0)
      r0s                         <- c(r0s, R0)
    }
    tslength                      <- max(unlist(lapply(tslist,length)))
    ts                            <- rep(0,tslength)
    for(jj in 1:nwaves){
      tslist[[jj]]                <- c(tslist[[jj]],rep(0,tslength - length(tslist[[jj]])))
      ts                          <- ts + tslist[[jj]]
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
    
    mult_squish                   <- runif(1,1,2)
    mult_period                   <- length(ts)/runif(1,sqrt(2),sqrt(10))^2
    mult_sin                      <- 1+mult_squish*(1+sample(c(-1,1),1)*sin((pi*(1:length(ts)))/mult_period))
    ts                            <- pmax(0,mult_sin*ts)
    ts = ts/max(ts)
  }
  
  ################
  ### SEASONAL ###
  ################
  if(curve_type == "seasonal"){
    seas_bool = T
    gamma                         <- 1/invgamma[1]
    # init_cond <- MCMCpack::rdirichlet(1,c(1000,1,1000))
    init_cond                     <- MCMCpack::rdirichlet(1,c(1000,0.1,0.1))
    S0                            <- init_cond[1]
    I0                            <- init_cond[2]
    R0                            <- init_cond[3]
    beta                          <- (basic_repo[1]/S0)*gamma
    times                         <- 1:500
    sir_output                    <- sir(beta, gamma, S0, I0, R0, times)
    ts                            <- sir_output$I
    cnt = 1
    while(sum(is.nan(ts))!=0){
      init_cond                   <- MCMCpack::rdirichlet(1,c(1000,0.1,0.1))
      S0                          <- init_cond[1]
      I0                          <- init_cond[2]
      R0                          <- init_cond[3]
      beta                        <- (basic_repo[1]/S0)*gamma
      times                       <- 1:500
      sir_output                  <- sir(beta, gamma, S0, I0, R0, times)
      ts                          <- sir_output$I
      cnt = cnt+1
      if(cnt == 5){
        stop("Failed to generate SIR.")
      }
    }
    
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

    temp_incidences               <- sir_output$incidence
      temp_incidences[1]          <- 0
      max_incidence_time          <- which.max(temp_incidences) - 1
      max_incidence_value         <- max(temp_incidences)
      PITs                        <- c(PITs, max_incidence_time)
      PIVs                        <- c(PIVs, max_incidence_value)
      s0s                         <- c(s0s, S0)
      i0s                         <- c(i0s, I0)
      r0s                         <- c(r0s, R0)
  }
  
  ##################################
  ### Additional Transformations ###
  ##################################
  
  
  ## pad the beginning and end with 0s
  if(curve_type != "seasonal"){
    ts <- suppressWarnings(
                        c(
                          runif(20, 0, quantile(ts, prob = 0.01, na.rm = TRUE) + 1e-9),
                          ts,
                          runif(20, 0, quantile(ts, prob = 0.01, na.rm = TRUE) + 1e-9)
                        )
                      )
  }

  
  time_cadence = 'weekly'
  if(curve_type == 'seasonal'){
    time_cadence                  <- c("monthly","weekly")[sample(1:2,1)]
    if(time_cadence == "monthly"){
      index = floor(1:length(ts)/4)
      ts                          <- aggregate(ts~index, FUN = sum)$ts
      if(max(ts)>1){
        ts = ts/max(ts)
      }
    }
  }

  ## choose scale
  # ts_scale                        <- c("proportion","counts")[sample(1:2,1)]
  ts_scale = "counts" #only doing counts for now.
  
  
  ## Add Noise and Rescale
  if(any(ts<0)) error("runif issue occured.")

  ts[ts<0]                        <- 1e-9
  ## add noise
  obs_ts                        <- rbeta(1:length(ts), alpha*ts, alpha*(1-ts))
  
  
  ## put the peak on a reasonable scale
  
  obs_ts                        <- (obs_ts/max(obs_ts))*PI
  
  ## upscale
  if(ts_scale == "counts"){
    
    obs_ts                      <- round(Npop*obs_ts,0)
  }

  if(any(is.nan(obs_ts))) return(NULL)
  
  ## if peak happens too soon, pad the start of the time series with 0s
  if(which.max(obs_ts) < 30 & curve_type != 'seasonal'){
    obs_ts                        <- c(runif(sample(10:17,1),0,quantile(obs_ts,prob = .05)+1e-9),obs_ts)
  }
  
  ## don't allow anything to be less than 1e-10
  obs_ts                          <- pmax(1e-10, obs_ts)
  
  ## append to list if no error occurred
  if(!is.na(sum(obs_ts))){
    templist                      <- list(ts = obs_ts,
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
                     ts_seasonal = seas_bool, 
		     PITs = PITs,
		     PIVs = PIVs,
		     s0s = s0s,
		     i0s = i0s,
		     r0s = r0s)
  }
  return(templist)
}

#' @export
create_embed_matrix               <- function(synthetic, h, k = 4){
  s_idx                           <- 1
  for (s in synthetic){
    s$ts_id                       <- s_idx
    synthetic[[s_idx]]            <- s
    s_idx                         <- s_idx + 1
  }
  embed_mat                       <- lapply(synthetic,function(x){ embed( pmax(1e-8,x$ts),k+h)})
  embed_mat                       <- do.call(rbind,embed_mat)
  embed_mat                       <- embed_mat[,ncol(embed_mat):1] #casey
  RowVar                          <- function(x, ...) {
    rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
  }
  rows_to_delete = RowVar(embed_mat[,1:k])
  embed_mat                       <- embed_mat[which(rows_to_delete > 0),]
  embed_mat_X                     <- embed_mat[,1:k]
  embed_mat_y                     <- embed_mat[,(k+1):(k+h)]
  ret_list                        <- list()
  ret_list[[1]]                   <- embed_mat_X
  ret_list[[2]]                   <- embed_mat_y
  return (ret_list)
}


### This function generates the epiFFORMA features from a time series.
### input: ts = single time series
### output: matrix of time series
#' @export
make_features <- function(info_packet, h){
  suppressPackageStartupMessages({
    library(zoo)
    library(quantmod)   # if you load it
  })

  ## define ts
  ts <- info_packet$ts
  
  ## make everything at least as big as 1e-10
  minval <- 1e-10
  ts <- pmax(minval, ts)
  
  ## make the output dataframe
  tsf <- data.frame(h = 1:h)
  
  
  ## fit a gam 
  if(var(ts[pmax(length(ts)-15,1):length(ts)]) > minval & length(unique(ts[pmax(length(ts)-15,1):length(ts)])) > 3){
    
    
    ## make data frame for gam. Used outlier-cleaned time series.
    smooth_df <- data.frame(x = 1:length(ts),
                            y = ts,
                            y_minus1 = c(ts[-length(ts)],NA), #needed for low counts rollmean without last value
                            wt = 1:length(ts))
    
    last16id <- pmax(length(ts)-15,1):length(ts)
    last3id <- pmax(length(ts)-2,1):length(ts)
    
    ## Low Counts: apply 3 week rolling median (smoothing to stabilize) and model entire time series
    MAXIT = 200 #default value, can change if need to speed up in future
    gam_family = 'gaussian'
    if(info_packet$ts_scale == 'counts' & min(smooth_df$y[last3id[-length(last3id)]]) <= 20){
      smooth_df$y = rollapply(smooth_df$y, align = 'center', width = 3, FUN = function(x){median(x,na.rm=T)}, partial = T)
      smooth_df$y_minus1 = rollapply(smooth_df$y_minus1, align = 'center', width = 3, FUN = function(x){median(x,na.rm=T)}, partial = T)
      smooth_df$y_minus1[length(smooth_df$y_minus1)]=NA
      smooth_df = smooth_df[pmax(length(ts)-15,1):length(ts),]
      gam_mod <- suppressWarnings(try(mgcv::bam(y ~ s(x, bs = "ps", m = c(2,1)), data = smooth_df, weights = wt, method = "fREML", discrete = T, select = T, control = list(maxit = MAXIT), family = gam_family), silent = T))
      gam_mod_minus1 <- suppressWarnings(try(mgcv::bam(y_minus1 ~ s(x, bs = "ps", m = c(2,1)), data = smooth_df[-length(smooth_df[,1]),], weights = wt, method = "fREML", discrete = T, select = T, control = list(maxit = MAXIT), family = gam_family), silent = T))
    }else{
      smooth_df = smooth_df[pmax(length(ts)-15,1):length(ts),]
      gam_mod <- suppressWarnings(try(mgcv::bam(y ~ s(x, bs = "ps",  m = c(2,1)), data = smooth_df, weights = wt, method = "fREML", discrete = T, select = T, control = list(maxit = MAXIT), family = gam_family), silent = T))
      gam_mod_minus1 <- suppressWarnings(try(mgcv::bam(y_minus1 ~ s(x, bs = "ps", m = c(2,1)), data = smooth_df[-length(smooth_df[,1]),], weights = wt, method = "fREML", discrete = T, select = T, control = list(maxit = MAXIT), family = gam_family), silent = T))
    }
    
    #############################
    ### Get Smoothed Features ###
    #############################
    
    ### With Last Point
    if(ifelse(class(gam_mod)[1] == 'try-error', TRUE, ifelse(unlist(gam_mod$sp) > 10^6, TRUE, FALSE))){
      smooth_df$y = rollapply(smooth_df$y, align = 'center', width = 4, FUN = function(x){median(x,na.rm=T)}, partial = T)
      ts_smooth = ts
      ts_smooth[smooth_df$x] = smooth_df$y
      pred_ts_smooth = pmax(0,rep(smooth_df$y[length(smooth_df$y)],h))
    }else{
      gam_pred <- predict(gam_mod, newdata = data.frame(x = (length(ts)+1):(length(ts)+h)), type = 'response')
      ts_smooth = ts
      ts_smooth[smooth_df$x] = pmax(0,gam_mod$fitted.values)
      pred_ts_smooth = pmax(0,gam_pred)
    }
    
    ### Without Last Point
    if(ifelse(class(gam_mod_minus1)[1] == 'try-error', TRUE, ifelse(unlist(gam_mod_minus1$sp) > 10^6, TRUE, FALSE))){
      smooth_df$y_minus1 = rollapply(smooth_df$y_minus1, align = 'center', width = 4, FUN = function(x){median(x,na.rm=T)}, partial = T)
      smooth_df$y_minus1[length(smooth_df$y_minus1)] = NA
      ts_smooth_minus1 = ts
      ts_smooth_minus1[smooth_df$x] = smooth_df$y_minus1
      pred_ts_smooth_minus1 = pmax(0,rep(smooth_df$y_minus1[length(smooth_df$y)-1],h))
    }else{
      gam_pred_minus1 <- predict(gam_mod_minus1, newdata = data.frame(x = (length(ts)+1):(length(ts)+h)), type = 'response')
      ts_smooth_minus1 = ts
      ts_smooth_minus1[smooth_df$x] = c(pmax(0,gam_mod_minus1$fitted.values),NA)
      pred_ts_smooth_minus1 <- pmax(0,gam_pred_minus1)
    }
    
  }else{
    ts_smooth <- ts
    pred_ts_smooth <- rep(ts[length(ts)],h)
    pred_ts_smooth_minus1 = rep(ts[length(ts)],h)
  }
  
  
  #####################
  ### LOCAL METRICS ###
  #####################
  last10id <- (length(ts_smooth)-9):length(ts_smooth)
  
  #### Goal: Is the last jump positive or negative, and how much, multiplicatively? 
  ## ratio of gam[t] / gam[t-1]
  if(ts_smooth[length(ts_smooth)-1] <= minval | ts[length(ts)-1] <= minval |
     ts_smooth[length(ts_smooth)] <= minval | ts[length(ts)] <= minval){
    tsf$gr12_div_23 = 0 #pretty flat, possibly just getting started or just ending
  }else{
    tsf$gr12_div_23 = tanh(( (ts_smooth[length(ts_smooth)] + ts_smooth[length(ts_smooth)-1])/(ts_smooth[length(ts_smooth)-1]+ts_smooth[length(ts_smooth)-2])) - 1)
    
  }
  
  #### Goal: How big is the value relative to the values that have been observed recently?
  ts_smooth_smooth = rollapply(ts_smooth, align = 'center', width = 3, FUN = function(x){mean(x,na.rm=T)}, partial = T)
  ## last observation divided by max observation [0,1], smooth
  if(min(ts_smooth_smooth[last10id]) < max(ts_smooth_smooth[last10id]) & max(ts_smooth_smooth[last10id]) > minval){
    tsf$last_div_max <- (ts_smooth_smooth[length(ts_smooth_smooth)] - min(ts_smooth_smooth[last10id]))/(max(ts_smooth_smooth[last10id]) - min(ts_smooth_smooth[last10id]))
  }else{
    tsf$last_div_max <- 0 #saying last obs is the minimum
  }
  
  #### Goal: What does the signal to noise ratio look like in recent past? 
  ## coefvar
  if(mean(ts[last10id])>0){
    tsf$coefvar <- tanh(0.1*((sd(ts[last10id] - ts_smooth[last10id])/mean(ts_smooth[last10id]))-1)) 
  }else{
    tsf$coefvar <- tanh(-0.1) #changed on 5/8/24
  }
  
  ### Goal: Changes with h, might help algorithm really separate out "risky" longer-term forecasts when combined with other metrics
  tsf$gam_with_div_without <- tanh(ifelse(pred_ts_smooth_minus1==0, rep(0,h), 0.1*((pred_ts_smooth/pred_ts_smooth_minus1) - 1)))
  
  
  ##########################
  ### Global-ish Metrics ### (Last 2 Years)
  ##########################
  
  if(info_packet$ts_time_cadence == 'weekly'){
    last2ys <- pmax(1,(length(ts_smooth)-(52*2))):length(ts_smooth)
  }else if(info_packet$ts_time_cadence == 'monthly'){
    last2ys <- pmax(1,(length(ts_smooth)-(13*2))):length(ts_smooth) 
  }else if(info_packet$ts_time_cadence == 'monthly_12'){
    last2ys <- pmax(1,(length(ts)-(12*2))):length(ts)
  }else if(info_packet$ts_time_cadence == 'daily'){
    last2ys <- pmax(1,(length(ts_smooth)-(365*2))):length(ts_smooth)
  }else{
    last2ys <- 1:length(ts_smooth)
  }
  
  #### Goal: How big is the value, relative to the mean of past values, subset num to last to to avoid length of ts effect
  ## mean(y_1:t)/mean(y_1:(t-1)), smooth and normalized
  if(mean(ts_smooth[last2ys[-length(last2ys)]])>minval){
    tsf$avg_recent_div_avg_global <- tanh((mean(ts_smooth[last10id])/mean(ts_smooth[last2ys[-length(last2ys)]])) - 1) 
  }else{
    tsf$avg_recent_div_avg_global = 0 
  }
  
  ### Goal: Is the last jump an outlier relative to other jumps that have been ever seen? Like a global riskiness measure
  ts_smooth_diff <- diff(ts_smooth[last2ys]) 
  if(sd(ts_smooth_diff)>minval){
    ts_smooth_diff_z <- (ts_smooth_diff - mean(ts_smooth_diff))/sd(ts_smooth_diff)
    tsf$diff_zscore <- tanh(0.1 * ts_smooth_diff_z[length(ts_smooth_diff_z)])
  }else{
    tsf$diff_zscore <- 0
  }
  
  ### Goal: Capture the forecastability of time series if you don't condition on anything else
  if(var(ts) == 0){
    tsf$entropy = 1 
  }else{
    tsf$entropy = tsfeatures::entropy(ts)
  }
  
  
  ## Consecutive increase divided by max consecutive increase
  ts_smooth_diff <- as.numeric(diff(ts_smooth)>0)
  if(ts_smooth_diff[length(ts_smooth_diff)]==1){
    first = min(which(rev(ts_smooth_diff) == 0)) #first different value from end
  }
  MAT = data.frame(lengths = rle(ts_smooth_diff == 1)$lengths, values = rle(ts_smooth_diff == 1)$values)
  MAT = MAT[MAT$values,]
  tsf$relative_increases = ifelse(ts_smooth_diff[length(ts_smooth_diff)] == 0, 0, mean( MAT$lengths<(first-1)))
  
  ### Proportion of Year Since Average Yearly Max
  if(info_packet$ts_time_cadence == 'weekly'){
    freq_decomp = 52
  }else if(info_packet$ts_time_cadence == 'monthly'){
    freq_decomp = 13 
  }else if(info_packet$ts_time_cadence == 'monthly_12'){
    freq_decomp <- 12
  }else if(info_packet$ts_time_cadence == 'daily'){
    freq_decomp = 365
  }else{
    freq_decomp = 1
  }
  lags = (1:length(ts)) %% freq_decomp
  AGG = aggregate(ts~lags, FUN = mean)
  tsf$prop_since_peak = ((length(ts)-(which.max(AGG$ts)-1)) %% freq_decomp)/freq_decomp #proportion of year since max counts location
  
  ### Seasonality
  if(sum(ts > minval) == 0){ #if all zeros, metric is zero
    tsf$seasonality = 0
  }else{
    if(min(which(ts > minval)) == 1){
      ts_sub = ts
    }else{
      ts_sub = ts[-c(1:(min(which(ts > minval))-1))] #zeros at the beginning of ts confuse this metric
    }
    if(length(ts_sub) >= (freq_decomp+1)){ #require at least freq_decomp times to calculate this metric, otherwise zero
      max_lag = pmin(ceiling(freq_decomp*1.5), length(ts_sub)-1)
      min_lag = floor(freq_decomp*0.5)
      lags = 0:max_lag 
      tsf$seasonality = pmax(0,max(as.numeric(unlist(acf(ts_sub,max_lag,plot=F)$acf[lags >= min_lag]))))
    }else{
      tsf$seasonality = 0
    }
  }
  
  ## get outta here
  return(tsf)
}

#' @export
generate_smoa_synthetic = function(x, curve_type, os_type, save_location){
  nwaves_bounds = c(2,10)
  invgamma1_bounds = c(1,10)
  invgamma2_bounds = c(1,10)
  repo1_bounds = c(1,10)
  repo2_bounds = c(1,10)
  alpha_bounds = c(50, 9999)
  Npop_bounds = c(200000, 99999999)
  PI_bounds = c(0.005, 0.25)
  
  templist <- gen_curve(
                        nwaves     = nwaves_bounds[1]     + x[1] * diff(nwaves_bounds),
                        invgamma1 = invgamma1_bounds[1] + x[2] * diff(invgamma1_bounds),
                        invgamma2 = invgamma2_bounds[1] + x[3] * diff(invgamma2_bounds),
                        repo1     = repo1_bounds[1]     + x[4] * diff(repo1_bounds),
                        repo2     = repo2_bounds[1]     + x[5] * diff(repo2_bounds),
                        alpha     = alpha_bounds[1]     + x[6] * diff(alpha_bounds),
                        Npop      = Npop_bounds[1]      + x[7] * diff(Npop_bounds),
                        PI        = PI_bounds[1]        + x[8] * diff(PI_bounds),
                        curve_type
                      )
  
  # little sloppy, but probabilistically should work out.
  identity_int = sample.int(1e15, 1)
  save(templist, file = paste0(save_location, "/data/", os_type, '/', identity_int, '.RData'))
  
  tsf = make_features(templist, 1)
  tsf$h = NULL
  return(tsf)
}