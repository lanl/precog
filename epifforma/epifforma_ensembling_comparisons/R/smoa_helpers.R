# Helper Functions for the run_smoa.R main file.
## Author: AC Murph

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

sir                               <- function(beta, gamma, S0, I0, R0, times) {
  require(deSolve) # for the "ode" function
  
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
  out                             <- ode(initial_values, times, sir_equations, parameters_values)
  # returning the output:
  as.data.frame(cbind(out, incidence = c(NA,-diff(as.data.frame(out)$S)), beta = beta, gamma = gamma))
}

#' Fetch the latest or earliest issue for each observation
#'
#' The data returned from `covidcast_signal()` or `covidcast_signals()` can, if
#' called with the `issues` argument, contain multiple issues for a single
#' observation in a single location. These functions filter the data frame to
#' contain only the earliest issue or only the latest issue.
#'
#' @param df A `covidcast_signal` or `covidcast_signal_long` data frame, such as
#'   returned from `covidcast_signal()` or the "long" format of
#'   `aggregate_signals()`.
#' @return A data frame in the same form, but with only the earliest or latest
#'   issue of every observation. Note that these functions sort the data frame
#'   as part of their filtering, so the output data frame rows may be in a
#'   different order.
#' @importFrom rlang .data
#' @export
latest_issue                      <- function(df) {
  return(first_or_last_issue(df, TRUE))
}

#' @rdname latest_issue
#' @export
earliest_issue                    <- function(df) {
  return(first_or_last_issue(df, FALSE))
}


#' Get state abbreviations from state names
#'
#' Look up state abbreviations by state names (including District of Columbia
#' and Puerto Rico); this function is based on `grep()`, and hence allows for
#' regular expressions.
#'
#' @param name Vector of state names to look up.
#' @param ignore.case,perl,fixed Arguments to pass to `grep()`, with the same
#'   defaults as in the latter function. Hence, by default, regular expressions
#'   are used; to match against a fixed string (no regular expressions), set
#'   `fixed = TRUE`.
#' @param ties_method If "first", then only the first match for each name is
#'   returned. If "all", then all matches for each name are returned.
#'
#' @return A vector of state abbreviations if `ties_method` equals "first", and
#'   a list of state abbreviations otherwise.
#'
#' @examples
#' name_to_abbr("Penn")
#' name_to_abbr(c("Penn", "New"), ties_method = "all")
#'
#' @seealso [abbr_to_name()]
#' @export
name_to_abbr                      <- function(name, ignore.case = FALSE, perl = FALSE, fixed = FALSE,
                        ties_method = c("first", "all")) {
  # First get rid of United States from state_census
  state_census                    <- read.csv("data/state_census.csv")
  state_census$X                  <- NULL
  df                              <- state_census %>% dplyr::filter(.data$STATE > 0)

  # Now perform the grep-based look up
  grep_lookup(key = name, keys = df$NAME, values = df$ABBR,
              ignore.case = ignore.case, perl = perl, fixed = fixed,
              ties_method = ties_method)
}

#' Get state names from state abbreviations
#'
#' Look up state names by state abbreviations (including District of Columbia
#' and Puerto Rico); this function is based on `grep()`, and hence allows for
#' regular expressions.
#'
#' @param abbr Vector of state abbreviations to look up.
#' @param ignore.case,perl,fixed Arguments to pass to `grep()`, with the same
#'   defaults as in the latter function. Hence, by default, regular expressions
#'   are used; to match against a fixed string (no regular expressions), set
#'   `fixed = TRUE`.
#' @param ties_method If "first", then only the first match for each name is
#'   returned. If "all", then all matches for each name are returned.
#'
#' @return A vector of state names if `ties_method` equals "first", and a list
#'   of state names otherwise.
#'
#' @examples
#' abbr_to_name("PA")
#' abbr_to_name(c("PA", "PR", "DC"))
#'
#' @seealso [name_to_abbr()]
#' @export
abbr_to_name                      <- function(abbr, ignore.case = FALSE, perl = FALSE, fixed = FALSE,
                        ties_method = c("first", "all")) {
  # First get rid of United States from state_census
  state_census                    <- read.csv("data/state_census.csv")
  state_census$X                  <- NULL
  
  # First get rid of United States from state_census
  df                              <- state_census %>% dplyr::filter(.data$STATE > 0)

  # Perform the grep-based look up
  grep_lookup(key = abbr, keys = df$ABBR, values = df$NAME,
              ignore.case = ignore.case, perl = perl, fixed = fixed,
              ties_method = ties_method)
}

#' Get state abbreviations from FIPS codes
#'
#' Look up state abbreviations by FIPS codes (including District of Columbia and
#' Puerto Rico). Will match the first two digits of the input codes, so should
#' work for 5-digit county codes, or even longer tract and census block FIPS
#' codes.
#'
#' @param code Vector of FIPS codes to look up; will match the first two digits
#'             of the code. Note that these are treated as strings; the number
#'             1 will not match "01".
#'
#' @return A vector of state abbreviations.
#'
#' @examples
#' fips_to_abbr("42000")
#' fips_to_abbr(c("42", "72", "11"))
#'
#' @seealso [abbr_to_fips()]
#' @export
fips_to_abbr                      <- function(code)
{
  fips                            <- sprintf("%02d", covidcast::state_census$STATE)
  index                           <- match(substr(code, 1, 2), fips)
  state_census                    <- read.csv('data/state_census.csv')
  state_census$X                  <- NULL
  output                          <- covidcast::state_census$ABBR[index]
  names(output)                   <- fips[index]
  output
}

# This is the core lookup function
grep_lookup                       <- function(key, keys, values, ignore.case = FALSE, perl = FALSE,
                       fixed = FALSE,  ties_method = c("first", "all")) {
  ties_method                     <- match.arg(ties_method)

  # Only do grep lookup for unique keys, to keep the look sort. It's a common
  # use case to, say, call state_fips_to_name on a covidcast_signal data frame
  # over many days of data with many repeat observations of the same locations.
  unique_key                      <- unique(key)

  res                             <- vector("list", length(unique_key))
  for (i in seq_along(unique_key)) {
    ind                           <- grep(unique_key[i], keys, ignore.case = ignore.case, perl = perl,
               fixed = fixed)
    if (length(ind) == 0) {
      res[[i]]                    <- NA
      names(res[[i]])             <- unique_key[i]
    } else {
      res[[i]]                    <- values[ind]
      names(res[[i]])             <- keys[ind]
    }
  }

  # Restore to original length, including duplicate keys.
  res                             <- res[match(key, unique_key)]

  # If they ask for all matches, then return the list
  if (ties_method == "all") {
    return(res)
  }

  # Otherwise, format into a vector, and warn if needed
  if (length(unlist(res)) > length(key)) {
    warning(paste("Some inputs were not uniquely matched; returning only the",
               "first match in each case."),
         res = res, key = key, class = "grep_lookup_nonunique_match")
  }
  return(unlist(lapply(res, `[`, 1)))
}

# Simple convenience functions for FIPS formatting
format_fips                       <- function(fips) { sprintf("%05d", fips) }
format_state_fips                 <- function(fips) { sprintf("%02d000", fips) }

gen_curve = function(curve_type){
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
    
    ## draw number of waves
    nwaves                        <- sample(3:7,1)
    #nwaves <- sample(3:15,1)
    ## draw reproduction numbers with preference toward lower values
    #basic_repo <- exp(runif(nwaves,log(1),log(35))) 
    basic_repo                    <- 1.1+rbeta(nwaves, shape1 = 1, shape2 = 3)*20 #but, have a preference for lower values
    # hist(1.1+rbeta(1000, shape1 = 1, shape2 = 3)*20, breaks = 50)
    
    # x = seq(0,10,0.01)
    # plot(x,dgamma(x,shape=1.5,scale=1.5))
    # hist(exp(pmin(rgamma(1000,shape=1,scale=2),log(35))), breaks = 100)
    # 
    ### randomly choose if wave peak heights generally increasing, decreasing, or random
    order_type = sample(c('increasing', 'decreasing', 'random'), size = 1)
    if(order_type == 'increasing'){
      basic_repo_order            <- rev(sample(1:nwaves,nwaves,replace=F,prob = basic_repo))
      basic_repo                  <- basic_repo[basic_repo_order]
    }else if(order_type == 'decreasing'){
      basic_repo_order            <- sample(1:nwaves,nwaves,replace=F,prob = basic_repo)
      basic_repo                  <- basic_repo[basic_repo_order]
    }
    
    ## draw gamma values
    invgamma                      <- runif(nwaves, 1, 10)
    gamma                         <- 1/invgamma
    init_cond                     <- rdirichlet(nwaves,c(1000,0.1,0.1)) #everyone starts as susceptible
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
        basic_repo[jj]            <- 1.1+rbeta(nwaves, shape1 = 1, shape2 = 3)*20 #but, have a preference for lower values
        invgamma                  <- runif(1, 1, 10)
        gamma[jj]                 <- 1/invgamma
        # init_cond <- rdirichlet(1,c(1000,1,1000))
        init_cond                 <- rdirichlet(1,c(1000,0.1,0.1))
        S0                        <- init_cond[1]
        I0                        <- init_cond[2]
        R0                        <- init_cond[3]
        beta                      <- (basic_repo/S0)*gamma
        times                     <- 1:500
	sir_output                <- sir(beta, gamma[jj], S0, I0, R0, times)
        ts                        <- sir_output$I
        cnt = cnt+1
        if(cnt == 5){
          ts = rep(0.01, length(ts))
        }
      }
      max_cut = max(which(ts/max(ts) > .0001))
      if(is.infinite(max_cut)){
        max_cut = length(ts)
      }
      tslist[[jj]]                <- c(rep(0,starttime),ts[1:max_cut])
      ## pad the beginning and end with 0s
      tslist[[jj]]                <- c(runif(20,0,quantile(tslist[[jj]],prob=.01, na.rm=T)+1e-9),
                        tslist[[jj]],
                        runif(20,0,quantile(tslist[[jj]],prob=.01, na.rm=T)+1e-9))
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
    
    ## draw number of waves
    nwaves                        <- sample(3:7,1)
    #nwaves <- sample(3:15,1)
    
    ## draw reproduction numbers with preference toward lower values
    #basic_repo <- exp(runif(nwaves,log(1),log(35))) 
    basic_repo                    <- 1.1+rbeta(nwaves, shape1 = 1, shape2 = 3)*20 #but, have a preference for lower values
    
    ### randomly choose if wave peak heights generally increasing, decreasing, or random
    order_type = sample(c('increasing', 'decreasing', 'random'), size = 1)
    if(order_type == 'increasing'){
      basic_repo_order            <- rev(sample(1:nwaves,nwaves,replace=F,prob = basic_repo))
      basic_repo                  <- basic_repo[basic_repo_order]
    }else if(order_type == 'decreasing'){
      basic_repo_order            <- sample(1:nwaves,nwaves,replace=F,prob = basic_repo)
      basic_repo                  <- basic_repo[basic_repo_order]
    }

    ## draw gamma values
    invgamma                      <- runif(nwaves, 1, 10)
    gamma                         <- 1/invgamma
    init_cond                     <- rdirichlet(nwaves,c(1000,0.1,0.1)) #everyone starts as susceptible
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
        #basic_repo[jj] <- exp(pmax(0.15,pmin(rgamma(1,shape = 1, rate = 2), log(35)))) #but, have a preference for lower values
        basic_repo[jj]            <- 1.1+rbeta(nwaves, shape1 = 1, shape2 = 3)*20
        invgamma                  <- runif(1, 1, 10)
        gamma[jj]                 <- 1/invgamma
        # init_cond <- rdirichlet(1,c(1000,1,1000))
        init_cond                 <- rdirichlet(1,c(1000,0.1,0.1))
        S0                        <- init_cond[1]
        I0                        <- init_cond[2]
        R0                        <- init_cond[3]
        beta                      <- (basic_repo/S0)*gamma
        times                     <- 1:500
	sir_output                <- sir(beta, gamma[jj], S0, I0, R0, times)
        ts                        <- sir_output$I
        cnt = cnt+1
        if(cnt == 5){
          ts = rep(0.01, length(ts))
        }
      }
      max_cut = max(which(ts/max(ts) > .0001))
      if(is.infinite(max_cut)){
        max_cut = length(ts)
      }
      tslist[[jj]]                <- c(rep(0,starttime),ts[1:max_cut])
      
      ## pad the beginning and end with 0s
      tslist[[jj]]                <- c(runif(20,0,quantile(tslist[[jj]],prob=.01, na.rm=T)+1e-9),
                        tslist[[jj]],
                        runif(20,0,quantile(tslist[[jj]],prob=.01, na.rm=T)+1e-9))
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
    ## get parameters to draw SIR curves
    #basic_repo <- exp(runif(1,log(2),log(5)))
    nwaves                        <- sample(5:15,1)
    
    # basic_repo <- exp(pmax(0.7,pmin(rgamma(1,shape = 1, rate = 2), log(5)))) #but, have a preference for lower values
    basic_repo                    <- 1.2+rbeta(1, shape1 = 1, shape2 = 3)*10 #but, have a preference for lower values

    

    
    invgamma                      <- runif(1, 1, 10)
    gamma                         <- 1/invgamma
    # init_cond <- rdirichlet(1,c(1000,1,1000))
    init_cond                     <- rdirichlet(1,c(1000,0.1,0.1))
    S0                            <- init_cond[1]
    I0                            <- init_cond[2]
    R0                            <- init_cond[3]
    beta                          <- (basic_repo/S0)*gamma
    times                         <- 1:500
    sir_output                    <- sir(beta, gamma, S0, I0, R0, times)
    ts                            <- sir_output$I
    cnt = 1
    while(sum(is.nan(ts))!=0){
      # basic_repo <- exp(pmax(0.7,pmin(rgamma(1,shape = 1, rate = 2), log(5)))) #but, have a preference for lower values
      basic_repo                  <- 1.2+rbeta(1, shape1 = 1, shape2 = 3)*10 #but, have a preference for lower values
      invgamma                    <- runif(1, 1, 10)
      gamma                       <- 1/invgamma
      # init_cond <- rdirichlet(1,c(1000,1,1000))
      init_cond                   <- rdirichlet(1,c(1000,0.1,0.1))
      S0                          <- init_cond[1]
      I0                          <- init_cond[2]
      R0                          <- init_cond[3]
      beta                        <- (basic_repo/S0)*gamma
      times                       <- 1:500
      sir_output                  <- sir(beta, gamma, S0, I0, R0, times)
      ts                          <- sir_output$I
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
    ts                            <- c(runif(20,0,quantile(ts,prob=.011,na.rm=T)+1e-9),
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
    time_cadence                  <- c("monthly","weekly")[sample(1:2,1)]
    if(time_cadence == "monthly"){
      index = floor(1:length(ts)/4)
      ts                          <- aggregate(ts~index, FUN = sum)$ts
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
  ts_scale                        <- c("proportion","counts")[sample(1:2,1)]
  
  
  ## Add Noise and Rescale
    if(any(ts<0)) error("runif issue occured.")

  ts[ts<0]                        <- 1e-9
    ## add noise
    alpha                         <- exp(runif(1,log(5e1),log(1e4)))
    obs_ts                        <- rbeta(1:length(ts), alpha*ts, alpha*(1-ts))
    
    
    ## put the peak on a reasonable scale
    PI                            <- exp(runif(1,log(.0005),log(.25)))
    obs_ts                        <- (obs_ts/max(obs_ts))*PI
    
    ## upscale
    if(ts_scale == "counts"){
      Npop                        <- round(runif(1,sqrt(2e5),sqrt(1e8))^2,0)
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


gen_curve_lhs = function(curve_type, sir_curve_lhs){

  PITs                            <- c()
  PIVs                            <- c()
  s0s                             <- c()
  i0s                             <- c()
  r0s                             <- c()
  #########################
  ### SIR-ROLLERCOASTER ###
  #########################
  if(curve_type == "sir_rollercoaster"){
    seas_bool = F #not seasonal
    
    ## draw number of waves
    nwaves                        <- sample(3:7,1)
    idx_of_waves                  <- sample(1:nrow(sir_curve_lhs), size = nwaves)
    sir_curve_lhs_subset          <- sir_curve_lhs[idx_of_waves,]
    basic_repo                    <- sir_curve_lhs_subset[, 6]

    ### randomly choose if wave peak heights generally increasing, decreasing, or random
    order_type = sample(c('increasing', 'decreasing', 'random'), size = 1)
    if(order_type == 'increasing'){
      basic_repo_order            <- rev(sample(1:nwaves,nwaves,replace=F,prob = basic_repo))
      basic_repo                  <- basic_repo[basic_repo_order]
      sir_curve_lhs_subset        <- sir_curve_lhs_subset[basic_repo_order,]
    }else if(order_type == 'decreasing'){
      basic_repo_order            <- sample(1:nwaves,nwaves,replace=F,prob = basic_repo)
      sir_curve_lhs_subset        <- sir_curve_lhs_subset[basic_repo_order,]
      basic_repo                  <- basic_repo[basic_repo_order]
    }
   
    starttime                     <- 1

    ## waves become less frequent
    cadence                       <- sample(10:52,nwaves)
    cadence_order                 <- sample(1:nwaves,nwaves,replace=F,prob = max(cadence) + 1 - cadence)
    cadence                       <- cadence[cadence_order]
    tslist                        <- list()
    for(jj in 1:nwaves){

      S0                          <- sir_curve_lhs_subset[jj,3]
      I0                          <- sir_curve_lhs_subset[jj,4]
      R0                          <- sir_curve_lhs_subset[jj,5]
      beta                        <- sir_curve_lhs_subset[jj, 1]
      gamma                       <- sir_curve_lhs_subset[jj,2]
      times                       <- 1:500

      sir_output                  <- sir(beta, gamma, S0, I0, R0, times)
      ts                          <- sir_output$I

      cnt = 1
      while(sum(is.nan(ts))!=0){
	temp_idx                  <- sample(1:nrow(sir_curve_lhs), size = nwaves)
        sir_curve_lhs_subset[jj,] <- sir_curve_lhs[temp_idx,]

	S0                        <- sir_curve_lhs_subset[jj,3]
        I0                        <- sir_curve_lhs_subset[jj,4]
        R0                        <- sir_curve_lhs_subset[jj,5]
        beta                      <- sir_curve_lhs_subset[jj, 1]
        gamma                     <- sir_curve_lhs_subset[jj,2]

        times                     <- 1:500
	sir_output                <- sir(beta, gamma, S0, I0, R0, times)
        ts                        <- sir_output$I
        cnt = cnt+1
        if(cnt == 5){
          ts = rep(0.01, length(ts))
        }
      }

      temp_incidences             <- sir_output$incidence
      temp_incidences[1]          <- 0
      max_incidence_time          <- which.max(temp_incidences) - 1
      max_incidence_value         <- max(temp_incidences)
      PITs                        <- c(PITs, max_incidence_time)
      PIVs                        <- c(PIVs, max_incidence_value)
      s0s                         <- c(s0s, S0)
      i0s                         <- c(i0s, I0)
      r0s                         <- c(r0s, R0)

      max_cut = max(which(ts/max(ts) > .0001))
      if(is.infinite(max_cut)){
        max_cut = length(ts)
      }
      tslist[[jj]]                <- c(rep(0,starttime),ts[1:max_cut])
      ## pad the beginning and end with 0s
      tslist[[jj]]                <- c(runif(20,0,quantile(tslist[[jj]],prob=.01, na.rm=T)+1e-9),
                        tslist[[jj]],
                        runif(20,0,quantile(tslist[[jj]],prob=.01, na.rm=T)+1e-9))
      tslist[[jj]][is.nan(tslist[[jj]])] = 0
      starttime                   <- starttime + cadence[jj]
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

    ## draw number of waves
    nwaves                        <- sample(3:7,1)
    idx_of_waves                  <- sample(1:nrow(sir_curve_lhs), size = nwaves)
    sir_curve_lhs_subset          <- sir_curve_lhs[idx_of_waves,]
    basic_repo                    <- sir_curve_lhs_subset[, 6]

    ### randomly choose if wave peak heights generally increasing, decreasing, or random
    order_type = sample(c('increasing', 'decreasing', 'random'), size = 1)
    if(order_type == 'increasing'){
      basic_repo_order            <- rev(sample(1:nwaves,nwaves,replace=F,prob = basic_repo))
      basic_repo                  <- basic_repo[basic_repo_order]
      sir_curve_lhs_subset        <- sir_curve_lhs_subset[basic_repo_order,]
    }else if(order_type == 'decreasing'){
      basic_repo_order            <- sample(1:nwaves,nwaves,replace=F,prob = basic_repo)
      sir_curve_lhs_subset        <- sir_curve_lhs_subset[basic_repo_order,]
      basic_repo                  <- basic_repo[basic_repo_order]
    }  

    ## draw gamma values
    starttime                     <- 1
    
    ## waves become less frequent
    cadence                       <- sample(10:52,nwaves)
    cadence_order                 <- sample(1:nwaves,nwaves,replace=F,prob = max(cadence) + 1 - cadence)
    cadence                       <- cadence[cadence_order]
    tslist                        <- list()
    for(jj in 1:nwaves){
      S0                          <- sir_curve_lhs_subset[jj,3]
      I0                          <- sir_curve_lhs_subset[jj,4]
      R0                          <- sir_curve_lhs_subset[jj,5]
      beta                        <- sir_curve_lhs_subset[jj, 1]
      gamma                       <- sir_curve_lhs_subset[jj,2]
      times                       <- 1:500

      sir_output                  <- sir(beta, gamma, S0, I0, R0, times)
      ts                          <- sir_output$I
      cnt = 1
      while(sum(is.nan(ts))!=0){
        temp_idx                  <- sample(1:nrow(sir_curve_lhs), size = nwaves)
        sir_curve_lhs_subset[jj,] <- sir_curve_lhs[temp_idx,]

        S0                        <- sir_curve_lhs_subset[jj,3]
        I0                        <- sir_curve_lhs_subset[jj,4]
        R0                        <- sir_curve_lhs_subset[jj,5]
        beta                      <- sir_curve_lhs_subset[jj, 1]
        gamma                     <- sir_curve_lhs_subset[jj,2]
	
	times                     <- 1:500

	sir_output                <- sir(beta, gamma, S0, I0, R0, times)
        ts                        <- sir_output$I
        cnt = cnt+1
        if(cnt == 5){
          ts = rep(0.01, length(ts))
        }
      }

      temp_incidences             <- sir_output$incidence
      temp_incidences[1]          <- 0
      max_incidence_time          <- which.max(temp_incidences) - 1
      max_incidence_value         <- max(temp_incidences)
      PITs                        <- c(PITs, max_incidence_time)
      PIVs                        <- c(PIVs, max_incidence_value)
      s0s                         <- c(s0s, S0)
      i0s                         <- c(i0s, I0)
      r0s                         <- c(r0s, R0)

      max_cut = max(which(ts/max(ts) > .0001))
      if(is.infinite(max_cut)){
        max_cut = length(ts)
      }
      tslist[[jj]]                <- c(rep(0,starttime),ts[1:max_cut])
      
      ## pad the beginning and end with 0s
      tslist[[jj]]                <- c(runif(20,0,quantile(tslist[[jj]],prob=.01, na.rm=T)+1e-9),
                        tslist[[jj]],
                        runif(20,0,quantile(tslist[[jj]],prob=.01, na.rm=T)+1e-9))
      tslist[[jj]][is.nan(tslist[[jj]])] = 0
      starttime                   <- starttime + cadence[jj]
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
    ## get parameters to draw SIR curves
    #basic_repo <- exp(runif(1,log(2),log(5)))
    nwaves                        <- sample(5:15,1)
   
    idx_of_waves                  <- sample(1:nrow(sir_curve_lhs), size = nwaves)
    sir_curve_lhs_subset          <- sir_curve_lhs[idx_of_waves,]
    basic_repo                    <- sir_curve_lhs_subset[, 6]
   
    S0                            <- sir_curve_lhs_subset[jj,3]
    I0                            <- sir_curve_lhs_subset[jj,4]
    R0                            <- sir_curve_lhs_subset[jj,5]
    beta                          <- sir_curve_lhs_subset[jj, 1]
    gamma                         <- sir_curve_lhs_subset[jj,2]

    times                         <- 1:500

    sir_output                    <- sir(beta, gamma, S0, I0, R0, times)
    ts                            <- sir_output$I
    cnt = 1
    while(sum(is.nan(ts))!=0){
      temp_idx                    <- sample(1:nrow(sir_curve_lhs), size = nwaves)
      sir_curve_lhs_subset[jj,]   <- sir_curve_lhs[temp_idx,]

      S0                          <- sir_curve_lhs_subset[jj,3]
      I0                          <- sir_curve_lhs_subset[jj,4]
      R0                          <- sir_curve_lhs_subset[jj,5]
      beta                        <- sir_curve_lhs_subset[jj, 1]
      gamma                       <- sir_curve_lhs_subset[jj,2]
	    
      times                       <- 1:500
      sir_output                  <- sir(beta, gamma, S0, I0, R0, times)
      ts                          <- sir_output$I
      cnt = cnt+1
      if(cnt == 5){
        ts = rep(0.01, length(ts))
      }
    }

    temp_incidences               <- sir_output$incidence
    temp_incidences[1]            <- 0
    max_incidence_time            <- which.max(temp_incidences) - 1
    max_incidence_value           <- max(temp_incidences)
    PITs                          <- c(PITs, max_incidence_time)
    PIVs                          <- c(PIVs, max_incidence_value)
    s0s                           <- c(s0s, S0)
    i0s                           <- c(i0s, I0)
    r0s                           <- c(r0s, R0)

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
    ts                            <- c(runif(20,0,quantile(ts,prob=.011,na.rm=T)+1e-9),
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
  ts_scale                        <- c("proportion","counts")[sample(1:2,1)]
  
  
  ## Add Noise and Rescale

    ## add noise
    alpha                         <- exp(runif(1,log(5e1),log(1e4)))
    obs_ts                        <- rbeta(1:length(ts), alpha*ts, alpha*(1-ts))
    
    
    ## put the peak on a reasonable scale
    PI                            <- exp(runif(1,log(.0005),log(.25)))
    obs_ts                        <- (obs_ts/max(obs_ts))*PI
    
    ## upscale
    if(ts_scale == "counts"){
      Npop                        <- round(runif(1,sqrt(2e5),sqrt(1e8))^2,0)
      obs_ts                      <- round(Npop*obs_ts,0)
    }

  
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

#' Compute weighted interval score
#'
#' Computes weighted interval score (WIS), a well-known quantile-based
#' approximation of the commonly-used continuous ranked probability score
#' (CRPS). WIS is a proper score, and can be thought of as a distributional
#' generalization of absolute error. For example, see [Bracher et
#' al. (2020)](https://arxiv.org/abs/2005.12881) for discussion in the context
#' of COVID-19 forecasting.
#' 
#' @param quantile vector of forecasted quantiles
#' @param value vector of forecasted values
#' @param actual_value Actual value.
#' 
#' @export
weighted_interval_score           <- function(quantile, value, actual_value) {
  score_func_param_checker(quantile, value, actual_value, "weighted_interval_score")
  if (all(is.na(actual_value))) return(NA)
  # `score_func_param_checker` above has already checked for uniqueness, so we
  # can save a bit of effort and just take the first actual.
  actual_value                    <- actual_value[[1L]]

  value                           <- value[!is.na(quantile)]
  quantile                        <- quantile[!is.na(quantile)]
  
  # per Ryan: WIS is equivalent to quantile loss modulo an extra 0.5 AE term
  # for the median forecast (counted twice). 
  # 
  # update: WIS is now being redefined to match exactly, still some question
  # about the correct denominator but the formula seems to be  1 / (K + 0.5)
  # 
  # Finally, the multiplication by 2 is because alpha_k = 2*quantile_k
  # 
  med                             <- value[find_quantile_match(quantile, 0.5)]
  
  if (length(med) > 1L) return(NA)
  
  wis                             <- 2 * mean(pmax(
    quantile * (actual_value - value),
    (1 - quantile) * (value - actual_value), 
    na.rm = TRUE))
  
  return(wis)
}

#' Common parameter checks for score functions
#'
#' A set of common checks for score functions, meant to identify common causes
#' of issues. Avoids `assert_that` for speed.
#'
#' @param quantiles vector of forecasted quantiles
#' @param values vector of forecasted values
#' @param actual_value actual_value, either as a scalar or a vector of the same
#'   length as `quantiles` and `values`
#' @param id string to identify the caller of the function and displayed in
#'   error messages (recommended to be the parent function's name)
score_func_param_checker          <- function(quantiles, values, actual_value, id = ""){
  id_str = paste0(id, ": ")
  if (length(actual_value) > 1) {
    if (length(actual_value) != length(values)) {
      stop(paste0(id_str,
                  "actual_value must be a scalar or the same length",
                  " as values"))
    }
    actual_value = unique(actual_value)
  }

  if (length(actual_value) != 1) {
    stop(paste0(id_str,
                "actual_value must have exactly 1 unique value"))
  }
  if (length(quantiles) != length(values)) {
    stop(paste0(id_str,
                "quantiles and values must be of the same length"))
  }

  if (anyDuplicated(quantiles)) {
    stop(paste0(id_str,
                "quantiles must be unique.")
    )
  }
}

find_quantile_match               <- function(quantiles, val_to_match, tol=1e-8){
  return(abs(quantiles - val_to_match) < tol  & !is.na(quantiles))
}

#' Get FIPS or CBSA codes from county or metropolitan area locations
#'
#' @param data a vector of strings for fips code, CBSA codes, location names
#' such as "Hampshire County, MA", "Alabama", "United Kingdom".
#' A US county location names must include state abbreviation.
#' @param hub character vector, where the first element indicates the hub
#' from which to load forecasts. Possible options are `"US"`, `"ECDC"` and `"FluSight"`
#' @return A vector of FIPS or CBSA codes
#' @export
name_to_fips <- function(data, hub = c("US", "ECDC")){
  if (is.null(data)){
    return (NULL)
  }
  
  if (hub[1] == "US") {
    load("data/hub_locations.rda")
    fips_tmp <- hub_locations %>%
      dplyr::mutate(full_location_name = fips)
    locations <- dplyr::bind_rows(hub_locations, fips_tmp)
    if(!all(data %in% locations$full_location_name)){
      stop("Error in name_to_fips: Please provide valid location name.eg: Bullock County,AL")
    }
    return (locations[locations$full_location_name %in% data, ]$fips)
  } else if (hub[1] == "ECDC") {
    load("data/hub_locations_ecdc.rda")
    location_tmp <- hub_locations_ecdc %>%
      dplyr::mutate(location_name = location)
    locations <- dplyr::bind_rows(hub_locations_ecdc, location_tmp)
    if(!all(data %in% locations$location_name)){
      stop("Error in name_to_fips: Please provide valid location name.")
    }
    return(locations[locations$location_name %in% data, ]$location)
  } else if (hub[1] == "FluSight") {
    load("data/hub_locations_flusight.rda")
    fips_tmp <- hub_locations_flusight %>%
      dplyr::mutate(location_name = fips)
    fips_tmp[fips_tmp$location_name == "US",]$location_name <- "United States"
    locations <- dplyr::bind_rows(hub_locations_flusight, fips_tmp)
    if(!all(data %in% locations$location_name)){
      stop("Error in name_to_fips: Please provide valid location name.eg: Bullock County,AL")
    }
    return (locations[locations$location_name %in% data, ]$fips)
  }
}


