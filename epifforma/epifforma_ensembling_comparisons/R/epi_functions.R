library(tsfeatures)
library(forecast)
library(randomForest)
library(data.table)
library(ranger)
library(lightgbm)
library(e1071)
library(deepgp)
library(FNN)
library(reshape2)
library(plyr)
library(collapse) #this is for moa speed up




create_embed_matrix <- function(synthetic, h, k = 4){
  s_idx <- 1
  for (s in synthetic){
    s$ts_id <- s_idx
    synthetic[[s_idx]] <- s
    s_idx <- s_idx + 1
  }
  embed_mat <- lapply(synthetic,function(x){ embed( pmax(1e-8,x$ts),k+h)})
  
  embed_mat <- do.call(rbind,embed_mat)
  embed_mat <- embed_mat[,ncol(embed_mat):1] #casey
  
  #rows_to_delete <- apply(embed_mat[,1:k],1,sd)
  ### below code is WAAAAY faster
  RowVar <- function(x, ...) {
    rowSums((x - rowMeans(x, ...))^2, ...)/(dim(x)[2] - 1)
  }
  rows_to_delete = RowVar(embed_mat[,1:k])

  embed_mat <- embed_mat[which(rows_to_delete > 0),]
  embed_mat_X <- embed_mat[,1:k]
  
  
  embed_mat_y <- embed_mat[,(k+1):(k+h)]
  # embed_mat_y_scaled <- (embed_mat_y - embed_mat_X_mean)/embed_mat_X_sd
  
  ret_list <- list()
  ret_list[[1]] <- embed_mat_X
  ret_list[[2]] <- embed_mat_y
  return (ret_list)
}


## input: ts_mat,h
## output: fitted fforma

build_training_data <- function(info_packet, h, split_method = 'diff'){
  #info_packet = A
  # info_packet <- synthetic[[ts_id]]
  
  ## first split time series
  train_test <- split_ts(info_packet = info_packet, h = h, split_method)
  
  ## replace the time series for validation
  info_packet$ts <- train_test$ts_train
  
  ## get the length of the time series
  TT <- length(info_packet$ts)
  #print(TT)
  ## trim the dates 
  if(sum(is.na(info_packet$ts_dates)) == 0 & sum(is.null(info_packet$ts_dates)) == 0){
    info_packet$ts_dates <- info_packet$ts_dates[1:TT]
  }
  
  ## trim the exogenous vars
  if(sum(is.na(info_packet$ts_exogenous)) == 0 & sum(is.null(info_packet$ts_exogenous)) == 0){
    info_packet$ts_exogenous <- info_packet$ts_exogenous[1:TT]
  }
  

  ## handle outliers
  info_packet$ts = as.numeric(handle_outliers(info_packet))

  ## distribute zeros
  info_packet$ts = as.numeric(distribute_zeros(info_packet))
  
  ## compute the features of the valiation time series
  ts_features <- data.frame(make_features(info_packet = info_packet, h = h))
  
  ## component forecasts
  components <- make_components(info_packet, h = h)

  ## eval components using squared error
  evals_se <- eval_components_se(forecasts = t(components), truth = train_test$ts_test, seasonal = info_packet$ts_seasonal)
  
  ## eval components using order statistics
  evals_twobest <- eval_components_twobest(forecasts = t(components), truth = train_test$ts_test)
  
  ## combine features and evals
  temp_output_se <- merge(ts_features, evals_se, by = "h", all.y=T)
  temp_output_se$ts_length <- length(info_packet$ts)
  
  temp_output_twobest <- merge(ts_features, evals_twobest, by = "h", all.y=T)
  temp_output_twobest$ts_length <- length(info_packet$ts)
  
  ## store component model results
  components$truth = train_test$ts_test
  components$h = c(1:h)
  components$ts_length = length(info_packet$ts)

  ## get outta here
  return(list(features_se = temp_output_se, features_twobest = temp_output_twobest, components = components))
  
}

build_training_data_incluster <- function(info_packet, h, embed_list, split_method = 'diff'){
  #info_packet = A
  # info_packet <- synthetic[[ts_id]]

  ## first split time series
  train_test <- split_ts(info_packet = info_packet, h = h, split_method)

  ## replace the time series for validation
  info_packet$ts <- train_test$ts_train

  ## get the length of the time series
  TT <- length(info_packet$ts)
  #print(TT)
  ## trim the dates
  if(sum(is.na(info_packet$ts_dates)) == 0 & sum(is.null(info_packet$ts_dates)) == 0){
    info_packet$ts_dates <- info_packet$ts_dates[1:TT]
  }

  ## trim the exogenous vars
  if(sum(is.na(info_packet$ts_exogenous)) == 0 & sum(is.null(info_packet$ts_exogenous)) == 0){
    info_packet$ts_exogenous <- info_packet$ts_exogenous[1:TT]
  }


  ## handle outliers
  info_packet$ts = as.numeric(handle_outliers(info_packet))

  ## distribute zeros
  info_packet$ts = as.numeric(distribute_zeros(info_packet))

  ## compute the features of the valiation time series
  ts_features <- data.frame(make_features(info_packet = info_packet, h = h))

  ## component forecasts
  components <- make_components_incluster(info_packet, h = h, embed_list)

  ## eval components using squared error
  evals_se <- eval_components_se(forecasts = t(components), truth = train_test$ts_test, seasonal = info_packet$ts_seasonal)

  ## eval components using order statistics
  evals_twobest <- eval_components_twobest(forecasts = t(components), truth = train_test$ts_test)

  ## combine features and evals
  temp_output_se <- merge(ts_features, evals_se, by = "h", all.y=T)
  temp_output_se$ts_length <- length(info_packet$ts)

  temp_output_twobest <- merge(ts_features, evals_twobest, by = "h", all.y=T)
  temp_output_twobest$ts_length <- length(info_packet$ts)

  ## store component model results
  components$truth = train_test$ts_test
  components$h = c(1:h)
  components$ts_length = length(info_packet$ts)

  ## get outta here
  return(list(features_se = temp_output_se, features_twobest = temp_output_twobest, components = components))

}


build_validation_data <- function(info_packet, h, split_method = 'diff'){
  #info_packet = A
  # info_packet <- synthetic[[ts_id]]
  
  ## first split time series
  train_test <- split_ts(info_packet = info_packet, h = h, split_method)
  
  ## replace the time series for validation
  info_packet$ts <- train_test$ts_train
  
  ## get the length of the time series
  TT <- length(info_packet$ts)
  #print(TT)
  ## trim the dates 
  if(sum(is.na(info_packet$ts_dates)) == 0 & sum(is.null(info_packet$ts_dates)) == 0){
    info_packet$ts_dates <- info_packet$ts_dates[1:TT]
  }
  
  ## trim the exogenous vars
  if(sum(is.na(info_packet$ts_exogenous)) == 0 & sum(is.null(info_packet$ts_exogenous)) == 0){
    info_packet$ts_exogenous <- info_packet$ts_exogenous[1:TT]
  }
  
  
  ## handle outliers
  info_packet$ts = as.numeric(handle_outliers(info_packet))
  
  ## distribute zeros
  info_packet$ts = as.numeric(distribute_zeros(info_packet))
  
  ## compute the features of the valiation time series
  ts_features <- data.frame(make_features(info_packet = info_packet, h = h))
  
  ## component forecasts
  components <- make_components(info_packet, h = h)
  
  ## forecast intervals
  component_intervals <- make_component_intervals(info_packet, h = h)
  
  ## eval components using squared error
  evals_se <- eval_components_se(forecasts = t(components), truth = train_test$ts_test, seasonal = info_packet$ts_seasonal)
  
  ## eval components using order statistics
  evals_twobest <- eval_components_twobest(forecasts = t(components), truth = train_test$ts_test)
  
  ## combine features and evals
  temp_output_se <- merge(ts_features, evals_se, by = "h", all.y=T)
  temp_output_se$ts_length <- length(info_packet$ts)
  
  temp_output_twobest <- merge(ts_features, evals_twobest, by = "h", all.y=T)
  temp_output_twobest$ts_length <- length(info_packet$ts)
  
  ## store component model results
  components$truth = train_test$ts_test
  components$h = c(1:h)
  components$ts_length = length(info_packet$ts)
  
  component_intervals$truth = train_test$ts_test
  component_intervals$h = c(1:h)
  component_intervals$ts_length = length(info_packet$ts)
  
  ## get outta here
  return(list(features_se = temp_output_se, features_twobest = temp_output_twobest, components = components, component_intervals = component_intervals))
  
}




### input: ts = single time series, h = horizon
### output: list[[1]] = time series to train, list[[2]] = time series to test

split_ts <- function(info_packet, h, random=T, split_method = 'diff'){
  
  ## get ts
  ts <- info_packet$ts
  
  ## which is the first indices are larger than the min
  idx_bigger_min <- intersect(which(ts > min(ts)), 11:(length(ts)-h)) 
  
  ## get the ts difference
  ts_diff <- sqrt(c(0,abs(diff(ts))))
  if(max(ts_diff) == 0){
    ts_diff <- 1e-10 + ts_diff
  }
  
  ## get the eligible split points
  if(length(idx_bigger_min) > 0){
    sample_idx <- seq(min(idx_bigger_min), max(idx_bigger_min), 1)
  }else{
    sample_idx <- 11:(length(ts)-h)
  }
  
  ## get the split point
  if(split_method == 'diff'){
    if(length(sample_idx) == 1){ #lauren changed on 7/10/24 to bug fix
      T_to_split = sample_idx
    }else{
      T_to_split <- sample(x = sample_idx, size = 1, prob = ts_diff[sample_idx]) #original
    }
  }else if(split_method == 'diff_before'){
    T_to_split <- sample(x = sample_idx, size = 1, prob = c(ts_diff[-1],1e-10)[sample_idx]) #upweights times before big jump (maybe helps on upticks?)
  }else if(split_method == 'diff_sq'){
    T_to_split <- sample(x = sample_idx, size = 1, prob = ts_diff[sample_idx]^2) #really prefer big jumps
  }else{
    stop('Invalid split_method')
  }
  
  ## return list
  ret_list <- list(ts_train = ts[1:T_to_split],
                   ts_test = ts[(T_to_split + 1):(T_to_split + h)])
  
  ## get outta here
  return(ret_list)
  
}


## handle outliers. Call obs an outlier if
## outlier in at least two metrics AND
## not in interval of multiple outliers unless value is zero AND
## not last value unless last value is zero or > 10 x max
handle_outliers <- function(info_packet){
  ts = info_packet$ts
  if(max(ts[-length(ts)]) > 0){
    ### Initial Outlier Screening
    out = forecast::tsoutliers(ts, lambda = NULL, iterate = 5)
    inds = out$index
    if(length(inds)>1){
      not_outlier = inds[apply(cbind(inds),1,FUN = function(x,inds){min(abs(x-inds[inds!=x]))},inds= inds)==1]
      inds = inds[!(inds %in% not_outlier) | inds == length(ts)]
    }
    #changed to 5 and last two years by Lauren on 5/14/24
    if(info_packet$ts_time_cadence == 'weekly'){
      last2ys <- pmax(1,(length(ts)-(52*2))):length(ts)
    }else if(info_packet$ts_time_cadence == 'daily'){
      last2ys <- pmax(1,(length(ts)-(365*2))):length(ts)
    }else if(info_packet$ts_time_cadence == 'monthly'){
      last2ys <- pmax(1,(length(ts)-(13*2))):length(ts)#monthly is actually every 4 weeks for Project Tycho
    }else if(info_packet$ts_time_cadence == 'monthly_12'){
      last2ys <- pmax(1,(length(ts)-(12*2))):length(ts)#monthly is actually every 4 weeks for Project Tycho
    }else{
      last2ys <- 1:length(ts)
    }
    if(length(inds) >= 1 & ts[length(ts)] != 0 & (max(ts)/max(ts[last2ys[-length(last2ys)]])) < 5){
      inds = inds[inds != length(ts)]
    }
    if(length(inds) >= 1 & var(ts)>0){
      ts_clean = ts
      ts_clean[inds] = NA
      ts_clean = forecast::na.interp(ts_clean)
      ts = as.numeric(ts_clean)
    }
    
    ### Extra Zeros Screening using ziP model (important for diseases like Hep-A)
    ### Lauren changed from ziP on 5/10/24
    if(info_packet$ts_scale == 'counts'){
      smooth_df <- data.frame(x = 1:length(ts),y = round(as.numeric(ts)))
      K = pmax(round(length(as.numeric(ts))/10),5)
      MAXIT = 50
      gam_mod <- try(mgcv::bam(y ~ s(x, k = K, bs = "tp", m = c(2,1)), data = smooth_df, method = "fREML", gamma = 10, discrete = T, select = T, family = 'poisson', control = list(maxit = MAXIT)), silent = T)
      if(class(gam_mod)[1]!='try-error'){
        probs = apply(cbind(predict(gam_mod, type = 'response')), 1, FUN = function(x){ppois(0,lambda = x, lower.tail = TRUE)})
        inds = which(probs < 0.05 & smooth_df$y == 0)
        if(length(inds)>1){
          not_outlier = inds[apply(cbind(inds),1,FUN = function(x,inds){min(abs(x-inds[inds!=x]))},inds= inds)==1]
          inds = inds[!(inds %in% not_outlier)]
        }
        if(length(inds) >= 1 & var(ts)>0){
          ts_clean = ts
          ts_clean[inds] = NA
          ts_clean = forecast::na.interp(ts_clean)
          ts = ts_clean
        }
      }
    }
    
    if(info_packet$ts_scale == 'proportion'){
      smooth_df <- data.frame(x = 1:length(ts),y = as.numeric(ts))
      K = pmax(round(length(as.numeric(ts))/10),5)
      MAXIT = 50
      gam_mod <- try(mgcv::bam(y ~ s(x, k = K, bs = "tp", m = c(2,1)), data = smooth_df, method = "fREML", gamma = 10, discrete = T, select = T, family = 'gaussian', control = list(maxit = MAXIT)), silent = T)
      if(class(gam_mod)[1]!='try-error'){
        probs = apply(cbind(predict(gam_mod, type = 'response')), 1, FUN = function(x, sigma){pnorm(0,mean = x, sd = sigma, lower.tail = TRUE)}, sigma = sqrt(summary(gam_mod)$dispersion))
        inds = which(probs < 0.05 & smooth_df$y == 0)
        if(length(inds)>1){
          not_outlier = inds[apply(cbind(inds),1,FUN = function(x,inds){min(abs(x-inds[inds!=x]))},inds= inds)==1]
          inds = inds[!(inds %in% not_outlier)]
        }
        if(length(inds) >= 1){
          ts_clean = ts
          ts_clean[inds] = NA
          ts_clean = forecast::na.interp(ts_clean)
          ts = ts_clean
        }
      }
    }
  }
  return(ts)
}


distribute_zeros = function(info_packet){
  inds = which(info_packet$ts == 0)
  inds = inds[inds>min(which(info_packet$ts>0))]
  ts = info_packet$ts
  #Lauren changed to 5 point rollmean from 3 point rollmean on 5/14/24
  ts_smooth = zoo::rollapply(ts, align = 'center', width = 5, FUN = function(x){mean(x,na.rm=T)}, partial = T)
  to_update = unique(c(inds, inds - 1, inds + 1, inds - 2, inds + 2))
  to_update = to_update[to_update>0]
  to_update = to_update[to_update<=length(ts)]
  ts[to_update] = ts_smooth[to_update]
  if(info_packet$ts_scale == 'counts'){
    ts = round(ts) #get back zeros in zero/one regions
  }
  return(ts)
}

### input: ts = single time series
### output: matrix of time series
make_features <- function(info_packet, h){
  
  ## define ts
  ts <- info_packet$ts
  
  ## make everything at least as big as 1e-10
  minval <- 1e-10
  ts <- pmax(minval, ts)
  
  ## make the output
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
    
    ## Low Counts: apply 4 week rolling median (smoothing to stabilize) and model entire time series
    MAXIT = 200 #default value, can change if need to speed up in future
    gam_family = 'gaussian'
    if(info_packet$ts_scale == 'counts' & min(smooth_df$y[last3id[-length(last3id)]]) <= 20){
      smooth_df$y = zoo::rollapply(smooth_df$y, align = 'center', width = 3, FUN = function(x){median(x,na.rm=T)}, partial = T)
      smooth_df$y_minus1 = zoo::rollapply(smooth_df$y_minus1, align = 'center', width = 3, FUN = function(x){median(x,na.rm=T)}, partial = T)
      smooth_df$y_minus1[length(smooth_df$y_minus1)]=NA
      smooth_df = smooth_df[pmax(length(ts)-15,1):length(ts),]
      gam_mod <- try(mgcv::bam(y ~ s(x, bs = "ps", m = c(2,1)), data = smooth_df, weights = wt, method = "fREML", discrete = T, select = T, control = list(maxit = MAXIT), family = gam_family), silent = T)
      gam_mod_minus1 <- try(mgcv::bam(y_minus1 ~ s(x, bs = "ps", m = c(2,1)), data = smooth_df[-length(smooth_df[,1]),], weights = wt, method = "fREML", discrete = T, select = T, control = list(maxit = MAXIT), family = gam_family), silent = T)
    }else{
      smooth_df = smooth_df[pmax(length(ts)-15,1):length(ts),]
      gam_mod <- try(mgcv::bam(y ~ s(x, bs = "ps",  m = c(2,1)), data = smooth_df, weights = wt, method = "fREML", discrete = T, select = T, control = list(maxit = MAXIT), family = gam_family), silent = T)
      gam_mod_minus1 <- try(mgcv::bam(y_minus1 ~ s(x, bs = "ps", m = c(2,1)), data = smooth_df[-length(smooth_df[,1]),], weights = wt, method = "fREML", discrete = T, select = T, control = list(maxit = MAXIT), family = gam_family), silent = T)
    }
    
    #############################
    ### Get Smoothed Features ###
    #############################

    
    ### With Last Point
    if(ifelse(class(gam_mod)[1] == 'try-error', TRUE, ifelse(unlist(gam_mod$sp) > 10^6, TRUE, FALSE))){
      smooth_df$y = zoo::rollapply(smooth_df$y, align = 'center', width = 4, FUN = function(x){median(x,na.rm=T)}, partial = T)
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
      smooth_df$y_minus1 = zoo::rollapply(smooth_df$y_minus1, align = 'center', width = 4, FUN = function(x){median(x,na.rm=T)}, partial = T)
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
    
    # # # ## Compare Lines
    # plot(ts[smooth_df$x])
    # lines(ts_smooth[smooth_df$x], col = 'orange')
    # lines(ts_smooth_minus1[smooth_df$x], col = 'red')

  }else{
    ts_smooth <- ts
    pred_ts_smooth <- rep(ts[length(ts)],h)
    pred_ts_smooth_minus1 = rep(ts[length(ts)],h)
  }
  
  
  #####################
  ### LOCAL METRICS ###
  #####################
  last10id <- pmax((length(ts_smooth)-9),1):length(ts_smooth)
  
  #### Goal: Is the last jump positive or negative, and how much, multiplicatively? 
  ## ratio of gam[t] / gam[t-1]
  if(ts_smooth[length(ts_smooth)-1] <= minval | ts[length(ts)-1] <= minval |
     ts_smooth[length(ts_smooth)] <= minval | ts[length(ts)] <= minval){
    tsf$gr12_div_23 = 0 #pretty flat, possibly just getting started or just ending
  }else{
    tsf$gr12_div_23 = tanh(( (ts_smooth[length(ts_smooth)] + ts_smooth[length(ts_smooth)-1])/(ts_smooth[length(ts_smooth)-1]+ts_smooth[length(ts_smooth)-2])) - 1)
    
  }
  
  #### Goal: How big is the value relative to the values that have been observed recently?
  ### Changed to smooth_smooth by Lauren on 5/14/24
  ts_smooth_smooth = zoo::rollapply(ts_smooth, align = 'center', width = 3, FUN = function(x){mean(x,na.rm=T)}, partial = T)
  ## last observation divided by max observation [0,1], smooth
  if(min(ts_smooth_smooth[last10id]) < max(ts_smooth_smooth[last10id]) & max(ts_smooth_smooth[last10id]) > minval){
    tsf$last_div_max <- (ts_smooth_smooth[length(ts_smooth_smooth)] - min(ts_smooth_smooth[last10id]))/(max(ts_smooth_smooth[last10id]) - min(ts_smooth_smooth[last10id]))
  }else{
    tsf$last_div_max <- 0 #saying last obs is the minimum
  }
  #lauren edited 5/8/24 to make more stable, changed interpretation
  # if(min(ts_smooth[last10id]) < max(ts_smooth[last10id]) & median(ts_smooth[last10id]) > minval){
  #   tsf$last_div_median <- tanh(0.1*((ts_smooth[length(ts_smooth)]/median(ts_smooth[last10id])) - 1))
  # }else{
  #   tsf$last_div_median <- 0 #saying last obs is the minimum
  # }
  # 
  
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
  ### Lauren changed to last two years on 5/14/24 to deal with waning diseases 
  
  if(info_packet$ts_time_cadence == 'weekly'){
    last2ys <- pmax(1,(length(ts_smooth)-(52*2))):length(ts_smooth)
  }else if(info_packet$ts_time_cadence == 'monthly'){
    last2ys <- pmax(1,(length(ts_smooth)-(13*2))):length(ts_smooth) #monthly is actually every 4 weeks for Project Tycho
  }else if(info_packet$ts_time_cadence == 'monthly_12'){
    last2ys <- pmax(1,(length(ts)-(12*2))):length(ts)#monthly is actually every 4 weeks for Project Tycho
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
    tsf$entropy = 1 #is 1 the right value here? 
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
    freq_decomp = 13 #monthly is actually every 4 weeks for Project Tycho
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
  

  # ts_smooth = sin(seq(0,100.5, 0.1))+1
  # #ts_smooth = 0.1*(1:length(ts_smooth))*ts_smooth
  # ts_smooth = (sin(seq(0,100.5, 0.1))+0.1)*ts_smooth
  # freq_decomp = 63
  # tsf$seasonality = as.numeric(info_packet$ts_seasonal)

  ## get outta here
  return(tsf)
  
}



### input: ts = single time series, h= horizon
### output: list of fitted models
make_components <- function(info_packet, h){
  
  ## define ts
  ts <- info_packet$ts
  
  ## make fcst holder
  ret_mat <- data.frame(h = 1:h)
  minval <- 1e-10

  ## make gam 
  if(var(ts[pmax(length(ts)-15,1):length(ts)]) > minval & length(unique(ts[pmax(length(ts)-15,1):length(ts)])) > 3){
    

    ## get rid of zeros for ts_gam
    ts_gam <- pmax(minval, ts)
    
    ## make data frame for gam
    smooth_df <- data.frame(x = 1:length(ts_gam), y = ts_gam, wt = 1:length(ts_gam))
    
    last16id <- pmax(length(ts)-15,1):length(ts)
    last3id <- pmax(length(ts)-2,1):length(ts)
    
    ## Apply 4 week rolling median (smoothing to stabilize)
   MAXIT = 200 #default value, can change if need to speed up in future
   gam_family = 'gaussian'
   if(info_packet$ts_scale == 'counts' & min(smooth_df$y[last3id[-length(last3id)]]) <= 20){
     smooth_df$y = zoo::rollapply(smooth_df$y, align = 'center', width = 3, FUN = function(x){median(x,na.rm=T)}, partial = T)
     smooth_df = smooth_df[pmax(length(ts)-15,1):length(ts),]
     gam_mod <- try(mgcv::bam(y ~ s(x, bs = "ps", m = c(2,1)), data = smooth_df, weights = wt, method = "fREML", discrete = T, select = T, control = list(maxit = MAXIT), family = gam_family), silent = T)
   }else{
    smooth_df = smooth_df[pmax(length(ts)-15,1):length(ts),]
     gam_mod <- try(mgcv::bam(y ~ s(x, bs = "ps",  m = c(2,1)), data = smooth_df, weights = wt, method = "fREML", discrete = T, select = T, control = list(maxit = MAXIT), family = gam_family), silent = T)
   }

    ### regular gam
    if(ifelse(class(gam_mod)[1] == 'try-error', TRUE, ifelse(unlist(gam_mod$sp) > 10^6, TRUE, FALSE))){
      smooth_df$y = zoo::rollapply(smooth_df$y, align = 'center', width = 4, FUN = function(x){median(x,na.rm=T)}, partial = T)
      ts_smooth = ts
      ts_smooth[smooth_df$x] = smooth_df$y
      ret_mat$gam = pmax(0,rep(smooth_df$y[length(smooth_df$y)],h))
    }else{
      gam_pred <- predict(gam_mod, newdata = data.frame(x = (length(ts_gam)+1):(length(ts_gam)+h)), type = 'response')
      ts_smooth = pmax(0,ts)
      ts_smooth[smooth_df$x] = pmax(0,gam_mod$fitted.values)
      ret_mat$gam <- pmax(0, gam_pred)
    }

  
      
  }else{
    ## ts_smooth
    ts_smooth <- pmax(0, ts) 
    
    ## add prediction
    ret_mat$gam <- pmax(0, rep(ts[length(ts)],h))
    #ret_mat$gam_smoother <- pmax(0, rep(mean(ts),h))
    
  }
  
  #### make MOA 
  k <- 4
  dist_to_test = rowSums(abs(embed_mat_X %r-% tail(ts,k))) 
  min_dist <- sort(dist_to_test,index.return = TRUE)$ix[1:2000]
  ret_mat$moa <- pmax(0,apply(embed_mat_y[min_dist,],2,median))
  rm('dist_to_test')

  
  #### make MOA Derivative 
  k <- 4 #this k = 4 b/c diff is inside tail call
  dist_to_test = rowSums(abs(embed_mat_X_deriv %r-% tail(diff(ts),k)))  
  min_dist <- sort(dist_to_test,index.return = TRUE)$ix[1:2000]
  ret_mat$moa_deriv <- pmax(0,tail(ts,1) + cumsum(apply(embed_mat_y_deriv[min_dist,],2,median)))
  rm('dist_to_test')
  
  ## make persistence model
  ret_mat$rw <- pmax(0,forecast::forecast(forecast::rwf(ts_smooth, drift = F), h = h)$mean)
  
  ## make theta, modified to refer to last 6 weeks by Lauren on 4/26/24
  ret_mat$theta <- pmax(0, forecast::forecast(forecast::thetaf(ts), h = h)$mean)
  
  ## make mean 
  ## Lauren edited to mean from last two years on 5/14/24
  if(info_packet$ts_time_cadence == 'weekly'){
    last2ys <- pmax(1,(length(ts_smooth)-(52*2))):length(ts_smooth)
  }else if(info_packet$ts_time_cadence == 'monthly'){
    last2ys <- pmax(1,(length(ts_smooth)-(13*2))):length(ts_smooth) #monthly is actually every 4 weeks for Project Tycho
  }else if(info_packet$ts_time_cadence == 'daily'){
    last2ys <- pmax(1,(length(ts_smooth)-(365*2))):length(ts_smooth)
  }else{
    last2ys <- 1:length(ts_smooth)
  }
  ret_mat$meanfcst <- pmax(0, rep(mean(ts[last2ys]),h))
  
  ## make mirror
  ret_mat$mirror <- pmax(0, ts_smooth[((length(ts_smooth)-1):(length(ts_smooth)-h))])
  
  ## make gam2mirror
  gam_wt <- seq(.9, .1, length.out=h)
  ret_mat$gam2mirror<- gam_wt*ret_mat$gam + (1-gam_wt)*ret_mat$mirror
  
  # ## arima 
  ret_mat$arima <- pmax(0,forecast::forecast(forecast::auto.arima(ts_smooth),h=h)$mean)

  
  ## make equal wt
  #ret_mat$equal_wt <- pmax(0, rowMeans(subset(ret_mat, select = setdiff(names(ret_mat),c("h")))))

  ## seasonal component
  if(info_packet$ts_time_cadence == 'weekly'){
    freq_decomp = 52
  }else if(info_packet$ts_time_cadence == 'monthly'){
    freq_decomp = 13 #monthly is actually every 4 weeks for Project Tycho
  }else if(info_packet$ts_time_cadence == 'monthly_12'){
    freq_decomp <- 12
  }else if(info_packet$ts_time_cadence == 'daily'){
    freq_decomp = 365
  }else{
    freq_decomp = 1
  }
  
  
  ## reformat ret_mat
  ret_mat <- subset(ret_mat, select=setdiff(names(ret_mat),"h"))
  
  ## get outta here
  return(ret_mat = ret_mat)
}


make_component_intervals <- function(info_packet, h){
  
  ## define ts
  ts <- info_packet$ts
  
  ## make fcst holder
  ret_mat <- data.frame(h = 1:h)
  minval <- 1e-10
  
  
  ## Apply 4 week rolling median (smoothing to stabilize)
  MAXIT = 200 #default value, can change if need to speed up in future
  gam_family = 'gaussian'
  
  # ts = sin(seq(0,100,0.1))
  # ts = ts + rnorm(length(ts),0,0.1)
  smooth_df <- data.frame(x = 1:length(ts), y = ts)
  gam_mod <- try(mgcv::bam(y ~ s(x, bs = "ps", m = c(2,1), k = floor(nrow(smooth_df)/2)), data = smooth_df, method = "fREML", discrete = T, select = T, control = list(maxit = MAXIT), family = gam_family), silent = T)
  theta_mod = try(forecast::thetaf(smooth_df$y, level = 95), silent = T)
  arima_mod = try(forecast::auto.arima(smooth_df$y, max.p = 10, max.q = 10), silent = T)
  
  if(ifelse(class(gam_mod)[1] == 'try-error', TRUE, ifelse(unlist(gam_mod$sp) > 10^6, TRUE, FALSE))){
    ret_mat$gam = pmax(0,2*1.96*sd(smooth_df$y))
  }else{
    preds = predict(gam_mod, newdata = data.frame(x = 1:length(ts)), type = 'response')
    res = predict(gam_mod, newdata = data.frame(x = (length(ts)+1):(length(ts)+h)), type = 'response', se.fit = T)
    preds_forecast = res$fit
    gam_error = try(mgcv::bam(error ~ s(ts, bs = "ps", m = c(2,1)), data = data.frame(ts = ts, error = abs(preds-ts)),method = "fREML", discrete = T, select = T, control = list(maxit = MAXIT), family = gam_family), silent = T)
    if(class(gam_error)[1]!='try-error'){
      preds_se = sqrt((pmax(0,predict(gam_error, newdata = data.frame(ts = preds_forecast)))^2)+(res$se.fit^2))
    }else{
      preds_se = sd(smooth_df$y)
    }
    ret_mat$gam <- pmax(0, 2*1.96*preds_se)
  }
  
  if(class(theta_mod)[1] == 'try-error'){
    ret_mat$theta = pmax(0,2*1.96*sd(smooth_df$y))
  }else{
    preds = theta_mod$fitted
    res = forecast::forecast(theta_mod, h = h, bootstrap = T)
    res$var = pmax(0,(((res$upper - res$lower)/(1.96*2))^2) - (theta_mod$model$sigma^2))
    res$se = pmax(0,sqrt(res$var))
    preds_forecast = res$mean
    gam_error = try(mgcv::bam(error ~ s(ts, bs = "ps", m = c(2,1)), data = data.frame(ts = ts, error = abs(preds-ts)),method = "fREML", discrete = T, select = T, control = list(maxit = MAXIT), family = gam_family), silent = T)
    if(class(gam_error)[1]!='try-error'){
      preds_se = sqrt((pmax(0,predict(gam_error, newdata = data.frame(ts = preds_forecast)))^2) + (res$se^2))
    }else{
      preds_se = sd(smooth_df$y)
    }
    ret_mat$theta <- pmax(0, 2*1.96*preds_se)
  }
  
  if(class(arima_mod)[1] == 'try-error'){
    ret_mat$arima = pmax(0,2*1.96*sd(smooth_df$y))
  }else{
    preds = arima_mod$fitted
    res = forecast::forecast(arima_mod, h = h, level = 95)
    res$var = pmax(0,(((res$upper - res$lower)/(1.96*2))^2) - (arima_mod$sigma2))
    res$se = pmax(0,sqrt(res$var))
    preds_forecast = res$mean
    gam_error = try(mgcv::bam(error ~ s(ts, bs = "ps", m = c(2,1)), data = data.frame(ts = ts, error = abs(preds-ts)),method = "fREML", discrete = T, select = T, control = list(maxit = MAXIT), family = gam_family), silent = T)
    if(class(gam_error)[1]!='try-error'){
      preds_se = sqrt((pmax(0,predict(gam_error, newdata = data.frame(ts = preds_forecast)))^2) + (res$se^2))
    }else{
      preds_se = sd(smooth_df$y)
    }
    ret_mat$arima <- pmax(0, 2*1.96*preds_se)
  }
  
  ## reformat ret_mat
  ret_mat <- subset(ret_mat, select=setdiff(names(ret_mat),"h"))
  
  ## get outta here
  return(ret_mat = ret_mat)
}



### input: forecasts, truth
### output: vector of class probabilities

eval_components_se <- function(forecasts,truth,seasonal = T){
  
  ## compute each component forecasting models MAE with the test time series
  err_by_component_and_horizon <- abs(forecasts - matrix(rep(truth,nrow(forecasts)),nrow=nrow(forecasts),byrow=T))^2
  
  ## return the best component forecasting model per horizon
  ## NOTE: When there is a tie, this returns all models
  labeldf <- NULL
  for(i in 1:ncol(err_by_component_and_horizon)){
    templabeldf <- data.table(h = i, class = 1:nrow(err_by_component_and_horizon))
    
    if(min(err_by_component_and_horizon[,i]) == 0){
      minidx <- which(err_by_component_and_horizon[,i] == min(err_by_component_and_horizon[,i]))
      templabeldf$class_wt <- 0
      templabeldf$class_wt[minidx] <- 1/length(minidx)
    }else{
      # if(!seasonal){
      #   s1 <- sum(err_by_component_and_horizon[rownames(err_by_component_and_horizon) != 'snaive',i])/err_by_component_and_horizon[rownames(err_by_component_and_horizon) != 'snaive',i]
      #   templabeldf$class_wt = 0 #snaive model give weight = 0 and otherwise not impacting other weights
      #   templabeldf$class_wt[rownames(err_by_component_and_horizon) != 'snaive'] <- s1/sum(s1)
      # }else{
        s1 <- sum(err_by_component_and_horizon[,i])/err_by_component_and_horizon[,i]
        templabeldf$class_wt <- s1/sum(s1)
      # }
    }
    labeldf <- rbind(labeldf, templabeldf)
  }
  
  ## only return the larger than equal weights
  labeldf <- subset(labeldf, class_wt > .99/length(unique(labeldf$class)))

  ## get outta here
  return(labeldf)
}


eval_components_twobest <- function(forecasts,truth){
  
  labeldf <- NULL
  for(jj in 1:ncol(forecasts)){
    
    ## forecasts too low
    minidfb <- data.frame(model = rownames(forecasts),
                          mae = abs(forecasts[,jj] - truth[jj]),
                          sign = sign(forecasts[,jj] - truth[jj]))
    
    minidf <- ddply(minidfb, .(sign), summarise, bestmodel = model[which.min(mae)][1])
    
    if(nrow(minidf) == 1){
      if(minidf$sign == 1){
        minidf <- rbind(minidf,
                        data.frame(sign = -1,
                                   bestmodel = "none"))
      }else{
        minidf <- rbind(minidf,
                        data.frame(sign = 1,
                                   bestmodel = "none"))
      }
    }
    labeldf <- rbind(labeldf, data.frame(h = jj, minidf))
  }
  
  ## get outta here
  return(labeldf)
}




## input: features, forecast components, ts scaling
## output: design matrix
make_design_matrix <- function(features, response = NULL){
  
  ## make the design matrix
  if(is.null(response)){
    temp_output <- cbind(data.frame(features))
  }else{
    temp_output <- cbind(data.frame(features),
                         response = response)
  }
  
  ## get outta here
  return(temp_output)
}






## input: output listfrom build_training_data
## output: fitted a lightgbm model object using the lightgbm package
## max_depth = -1 means no limit on depth
##
fit_lgbm_wt <- function(input, nmodels = 1, max_depth = -1, num_leaves = NA){
  
  ## unpack the training data into the features and class labels
  df <- na.omit(input)
  df$class <- as.factor(df$class)
  
  ## fit nmodels
  model_list <- list()
  for(i in 1:nmodels){
    
    ## divide df into train and validate
    df_train_ids <- sample(1:nrow(df), .8*nrow(df), replace=F) 
    df_train <- df[df_train_ids,]
    df_valid_id <- setdiff(1:nrow(df),df_train_ids)
    df_valid <- df[df_valid_id,]
    
    ## convert the data to LightGBM dataset format: training data
    lgbm_train <- lightgbm::lgb.Dataset(data = as.matrix(subset(df_train, select=setdiff(names(df_train),c("class","class_wt")))),
                                        label = as.integer(df_train$class)-1,
                                        weight = df_train$class_wt)
    
    ## convert the data to LightGBM dataset format: validation data
    lgbm_valid <- lightgbm::lgb.Dataset(data = as.matrix(subset(df_valid, select=setdiff(names(df_valid),c("class","class_wt")))),
                                        label = as.integer(df_valid$class)-1,
                                        weight = df_valid$class_wt,
                                        reference = lgbm_train)
    
    
    ## define parameters for LightGBM model
    params = list(objective = "multiclass",
                  metric = c("multi_logloss"),
                  num_class = length(unique(df$class)),
                  num_leaves = ifelse(!is.na(num_leaves),num_leaves,round(nrow(df_train)^(1/2),0)),
                  learning_rate = .1,
                  feature_fraction = 0.9,
                  max_depth = max_depth,
                  verbose = -1)
    
    ## train the LightGBM model
    lgb_model <- lightgbm::lgb.train(params = params,
                                     data = lgbm_train,
                                     valids = list(train = lgbm_train,
                                                   valid = lgbm_valid),
                                     early_stopping_rounds = 100,
                                     nrounds = 10000)
    
    model_list[[i]] <- lgb_model
  }
  
  ## get outta here
  return(model_list)
  
}








## input: output listfrom build_training_data
## output: fit a quantile regression lightgbm model object using the lightgbm package
## max_depth = -1 means no limit on depth
##
fit_lgbm_wt <- function(input, nmodels = 1, max_depth = -1, num_leaves = NA){
  
  ## unpack the training data into the features and class labels
  df <- na.omit(input)
  df$class <- as.factor(df$class)
  
  ## fit nmodels
  model_list <- list()
  for(i in 1:nmodels){
    
    ## divide df into train and validate
    df_train_ids <- sample(1:nrow(df), .8*nrow(df), replace=F) 
    df_train <- df[df_train_ids,]
    df_valid_id <- setdiff(1:nrow(df),df_train_ids)
    df_valid <- df[df_valid_id,]
    
    ## convert the data to LightGBM dataset format: training data
    lgbm_train <- lightgbm::lgb.Dataset(data = as.matrix(subset(df_train, select=setdiff(names(df_train),c("class","class_wt")))),
                                        label = as.integer(df_train$class)-1,
                                        weight = df_train$class_wt)
    
    ## convert the data to LightGBM dataset format: validation data
    lgbm_valid <- lightgbm::lgb.Dataset(data = as.matrix(subset(df_valid, select=setdiff(names(df_valid),c("class","class_wt")))),
                                        label = as.integer(df_valid$class)-1,
                                        weight = df_valid$class_wt,
                                        reference = lgbm_train)
    
    
    ## define parameters for LightGBM model
    params = list(objective = "multiclass",
                  metric = c("multi_logloss"),
                  num_class = length(unique(df$class)),
                  num_leaves = ifelse(!is.na(num_leaves),num_leaves,round(nrow(df_train)^(1/2),0)),
                  learning_rate = .1,
                  feature_fraction = 0.9,
                  max_depth = max_depth,
                  verbose = -1)
    
    ## train the LightGBM model
    lgb_model <- lightgbm::lgb.train(params = params,
                                     data = lgbm_train,
                                     valids = list(train = lgbm_train,
                                                   valid = lgbm_valid),
                                     early_stopping_rounds = 100,
                                     nrounds = 10000)
    
    model_list[[i]] <- lgb_model
  }
  
  ## get outta here
  return(model_list)
  
}




## input: info_packet, fitted lgb model, forecast horizon
## output: forecasts, lgb weights, features, components

predict_epifforma <- function(info_packet, feature2wt, h){
  
  ## handle outliers
  info_packet$ts = handle_outliers(info_packet)
  
  ## distribute zeros
  info_packet$ts = as.numeric(distribute_zeros(info_packet))
  
  ## compute the features
  features <- make_features(info_packet, h)
  
  ## compute the components
  components <- make_components(info_packet, h)

  ## predict the class probabilities
  pred_wts_list <- list()
  for(j in 1:length(feature2wt)){
    pred_wts_list[[j]] <- predict(feature2wt[[j]], newdata = as.matrix(features))
  }
  
  ## average all the fits
  pred_wts <- apply(simplify2array(pred_wts_list), 1:2, mean)
  pred_wts_sd <- apply(simplify2array(pred_wts_list), 1:2, sd)
  
  ## set epifforma equal to the mean forecast
  fcst_df <- data.frame(h = 1:h,
                        epifforma = rowSums(pred_wts*components),
                        components)
  
  
  ## melt it
  fcst_df_melt <- melt(fcst_df, id.vars = "h")
  fcst_df_melt$type <- "fcst"
  fcst_df_melt$x <- fcst_df_melt$h + i
  
  ##fcst truth
  fcst_truth <- data.frame(h = 1:h, truth = ts_test$ts[(i+1):(i+h)])
  fcst_df_melt <- merge(fcst_df_melt, fcst_truth, by="h", all.x=T)
  
  ## observation df
  obs_df <- data.frame(h = NA,
                       variable = "obs",
                       value = ts_test$ts[1:(length(info_packet$ts)+h)],
                       type = "obs",
                       x = 1:(length(info_packet$ts)+h),
                       truth = ts_test$ts[1:(length(info_packet$ts)+h)])
  
  ## combine observations and fcst
  output_df <- rbind(obs_df, fcst_df_melt)
  names(output_df) <- c("h","model","fcst","type","x","truth")
  
  ## data for plotting
  output_df$geography <- info_packet$ts_geography
  output_df$disease <- info_packet$ts_disease
  output_df$last_obs_time <- i
  output_df <- data.table(subset(output_df, select=c("geography","disease","last_obs_time","h",setdiff(names(output_df),c("geography","disease","last_obs_time","h")))))
  
  ## prepare to leave: features
  features$geography <- info_packet$ts_geography
  features$disease <- info_packet$ts_disease
  features$last_obs_time <- i
  features <- subset(features, select=c("geography","disease","last_obs_time","h",setdiff(names(features),c("geography","disease","last_obs_time","h"))))
  
  ## prepare to leave: pred_wts
  pred_wts <- data.frame(pred_wts)
  names(pred_wts) <- paste0(names(components),"_avg")
  pred_wts$geography <- info_packet$ts_geography
  pred_wts$disease <- info_packet$ts_disease
  pred_wts$last_obs_time <- i
  pred_wts$h <- 1:h
  pred_wts <- subset(pred_wts, select=c("geography","disease","last_obs_time","h",setdiff(names(pred_wts),c("geography","disease","last_obs_time","h"))))
  
  ## prepare to leave: pred_wts_sd
  pred_wts_sd <- data.frame(pred_wts_sd)
  names(pred_wts_sd) <- paste0(names(components),"_sd")
  pred_wts_sd$geography <- info_packet$ts_geography
  pred_wts_sd$disease <- info_packet$ts_disease
  pred_wts_sd$last_obs_time <- i
  pred_wts_sd$h <- 1:h
  pred_wts_sd <- subset(pred_wts_sd, select=c("geography","disease","last_obs_time","h",setdiff(names(pred_wts_sd),c("geography","disease","last_obs_time","h"))))
  
  ## prepare to leave: components
  components$geography <- info_packet$ts_geography
  components$disease <- info_packet$ts_disease
  components$last_obs_time <- i
  components$h <- 1:h
  components <- subset(components, select=c("geography","disease","last_obs_time","h",setdiff(names(components),c("geography","disease","last_obs_time","h"))))
  
  
  ## get outta here
  return(list(output_df = output_df,
              pred_wts = pred_wts,
              components = components,
              features = features,
              pred_wts_sd = pred_wts_sd))
  
}




#Lauren removed feature2wt dependence on 5/9/24 to avoid lots of copying of large model file
#will search in other environments
predict_epifforma_ordering <- function(info_packet, h){
  
  ## handle outliers
  info_packet$ts = handle_outliers(info_packet)
  
  ## distribute zeros
  info_packet$ts = as.numeric(distribute_zeros(info_packet))
  
  ## compute the features
  features <- make_features(info_packet, h)
  
  ## compute the components
  components <- make_components(info_packet, h)
  
  ## predict the class probabilities
  pred_wts_list <- list()
  for(j in 1:length(feature2wt)){
    pred_wts_list[[j]] <- predict(feature2wt[[j]], newdata = as.matrix(features))
  }
  rm(feature2wt) #trying to conserve memory here
  
  ## average all the fits
  pred_wts <- apply(simplify2array(pred_wts_list), 1:2, mean)
  pred_wts_sd <- apply(simplify2array(pred_wts_list), 1:2, sd)
  
  components_reordered = NULL
  pred_wts_reordered = NULL
  pred_wts_sd_reordered = NULL
  for(j in 1:nrow(components)){
    components_reordered = rbind(components_reordered, sort(unlist(components[j,])))
    pred_wts_reordered = rbind(pred_wts_reordered, pred_wts[j,order(order(unlist(components[j,])) )]) #double order inverts the ordering operaiton
    pred_wts_sd_reordered = rbind(pred_wts_sd_reordered, pred_wts_sd[j,order(order(unlist(components[j,])) )]) #double order inverts the ordering operaiton
  }
  
  ## set epifforma equal to the mean forecast
  fcst_df <- data.frame(h = 1:h,
                        epifforma = rowSums(pred_wts*components_reordered),
                        components)
  ## melt it
  fcst_df_melt <- melt(fcst_df, id.vars = "h")
  fcst_df_melt$type <- "fcst"
  fcst_df_melt$x <- fcst_df_melt$h + i
  
  ##fcst truth
  fcst_truth <- data.frame(h = 1:h, truth = ts_test$ts[(i+1):(i+h)])
  fcst_df_melt <- merge(fcst_df_melt, fcst_truth, by="h", all.x=T)
  
  ## observation df
  obs_df <- data.frame(h = NA,
                       variable = "obs",
                       value = ts_test$ts[1:(length(info_packet$ts)+h)],
                       type = "obs",
                       x = 1:(length(info_packet$ts)+h),
                       truth = ts_test$ts[1:(length(info_packet$ts)+h)])
  
  ## combine observations and fcst
  output_df <- rbind(obs_df, fcst_df_melt)
  names(output_df) <- c("h","model","fcst","type","x","truth")
  
  ## data for plotting
  output_df$geography <- info_packet$ts_geography
  output_df$disease <- info_packet$ts_disease
  output_df$last_obs_time <- i
  output_df <- data.table(subset(output_df, select=c("geography","disease","last_obs_time","h",setdiff(names(output_df),c("geography","disease","last_obs_time","h")))))
  
  ## prepare to leave: features
  features$geography <- info_packet$ts_geography
  features$disease <- info_packet$ts_disease
  features$last_obs_time <- i
  features <- subset(features, select=c("geography","disease","last_obs_time","h",setdiff(names(features),c("geography","disease","last_obs_time","h"))))
  
  ## prepare to leave: pred_wts
  pred_wts <- data.frame(pred_wts_reordered)
  names(pred_wts) <- paste0(names(components),"_avg")
  pred_wts$geography <- info_packet$ts_geography
  pred_wts$disease <- info_packet$ts_disease
  pred_wts$last_obs_time <- i
  pred_wts$h <- 1:h
  pred_wts <- subset(pred_wts, select=c("geography","disease","last_obs_time","h",setdiff(names(pred_wts),c("geography","disease","last_obs_time","h"))))
  
  ## prepare to leave: pred_wts_sd
  pred_wts_sd <- data.frame(pred_wts_sd_reordered)
  names(pred_wts_sd) <- paste0(names(components),"_sd")
  pred_wts_sd$geography <- info_packet$ts_geography
  pred_wts_sd$disease <- info_packet$ts_disease
  pred_wts_sd$last_obs_time <- i
  pred_wts_sd$h <- 1:h
  pred_wts_sd <- subset(pred_wts_sd, select=c("geography","disease","last_obs_time","h",setdiff(names(pred_wts_sd),c("geography","disease","last_obs_time","h"))))
  
  ## prepare to leave: components
  components$geography <- info_packet$ts_geography
  components$disease <- info_packet$ts_disease
  components$last_obs_time <- i
  components$h <- 1:h
  components <- subset(components, select=c("geography","disease","last_obs_time","h",setdiff(names(components),c("geography","disease","last_obs_time","h"))))
  
  
  ## get outta here
  return(list(output_df = output_df,
              pred_wts = pred_wts,
              components = components,
              features = features,
              pred_wts_sd = pred_wts_sd))
  
}

#will search in other environments
predict_epifforma_orderingtrunc <- function(info_packet, h, additional_features = F){
  
  ## handle outliers
  info_packet$ts = as.numeric(handle_outliers(info_packet))

  ## distribute zeros
  info_packet$ts = as.numeric(distribute_zeros(info_packet))
  
  ## compute the features
  features <- make_features(info_packet, h)
  
  ## compute the components
  components <- make_components(info_packet, h)
  
  if(additional_features){
    EDGEVALS = t(apply(as.matrix(components),1,FUN = function(x){c(min(x),max(x))}))
    EDGEVALS[,1] = EDGEVALS[,1]/components$rw
    EDGEVALS[,2] = EDGEVALS[,2]/components$rw
    EDGEVALS[(is.infinite(EDGEVALS) | EDGEVALS > 3) & !is.na(EDGEVALS)] = 3
    EDGEVALS[is.na(EDGEVALS)] = 0
    features = cbind(features, EDGEVALS)
  }

  ## predict the class probabilities
  pred_wts_list <- list()
  for(j in 1:length(feature2wt)){
    pred_wts_list[[j]] <- predict(feature2wt[[j]], newdata = as.matrix(features))
  }
  #rm(feature2wt) #trying to conserve memory here
  
  ## average all the fits
  pred_wts <- apply(simplify2array(pred_wts_list), 1:2, mean)
  pred_wts_sd <- apply(simplify2array(pred_wts_list), 1:2, sd)
  
  ##truncation
  pred_wts[pred_wts<0.01] = 0 #changed from 0.05 on 7/7/24
  # if(!info_packet$ts_seasonal){ #force snaive to have zero weight if not seasonal
  #   pred_wts[,ncol(pred_wts)] = 0 
  # }
  pred_wts = t(apply(pred_wts, 1, FUN = function(x){x/sum(x)}))
  
  ##if not seasonal, snaive is driven to infinity so will be the largest. Force zero weight on that component here

  
  components_reordered = NULL
  pred_wts_reordered = NULL
  pred_wts_sd_reordered = NULL
  for(j in 1:nrow(components)){
    components_reordered = rbind(components_reordered, sort(unlist(components[j,])))
    pred_wts_reordered = rbind(pred_wts_reordered, pred_wts[j,order(order(unlist(components[j,])) )]) #double order inverts the ordering operaiton
    pred_wts_sd_reordered = rbind(pred_wts_sd_reordered, pred_wts_sd[j,order(order(unlist(components[j,])) )]) #double order inverts the ordering operaiton
  }
  
  ## set epifforma equal to the mean forecast
  fcst_df <- data.frame(h = 1:h,
                        epifforma = rowSums(pred_wts*components_reordered),
                        components)
  ## melt it
  fcst_df_melt <- melt(fcst_df, id.vars = "h")
  fcst_df_melt$type <- "fcst"
  fcst_df_melt$x <- fcst_df_melt$h + i
  
  ##fcst truth
  fcst_truth <- data.frame(h = 1:h, truth = ts_test$ts[(i+1):(i+h)])
  fcst_df_melt <- merge(fcst_df_melt, fcst_truth, by="h", all.x=T)
  
  ## observation df
  obs_df <- data.frame(h = NA,
                       variable = "obs",
                       value = ts_test$ts[1:(length(info_packet$ts)+h)],
                       type = "obs",
                       x = 1:(length(info_packet$ts)+h),
                       truth = ts_test$ts[1:(length(info_packet$ts)+h)])
  
  ## combine observations and fcst
  output_df <- rbind(obs_df, fcst_df_melt)
  names(output_df) <- c("h","model","fcst","type","x","truth")
  
  ## data for plotting
  output_df$geography <- info_packet$ts_geography
  output_df$disease <- info_packet$ts_disease
  output_df$last_obs_time <- i
  output_df <- data.table(subset(output_df, select=c("geography","disease","last_obs_time","h",setdiff(names(output_df),c("geography","disease","last_obs_time","h")))))
  
  ## prepare to leave: features
  features$geography <- info_packet$ts_geography
  features$disease <- info_packet$ts_disease
  features$last_obs_time <- i
  features <- subset(features, select=c("geography","disease","last_obs_time","h",setdiff(names(features),c("geography","disease","last_obs_time","h"))))
  
  ## prepare to leave: pred_wts
  pred_wts <- data.frame(pred_wts_reordered)
  names(pred_wts) <- paste0(names(components),"_avg")
  pred_wts$geography <- info_packet$ts_geography
  pred_wts$disease <- info_packet$ts_disease
  pred_wts$last_obs_time <- i
  pred_wts$h <- 1:h
  pred_wts <- subset(pred_wts, select=c("geography","disease","last_obs_time","h",setdiff(names(pred_wts),c("geography","disease","last_obs_time","h"))))
  
  ## prepare to leave: pred_wts_sd
  pred_wts_sd <- data.frame(pred_wts_sd_reordered)
  names(pred_wts_sd) <- paste0(names(components),"_sd")
  pred_wts_sd$geography <- info_packet$ts_geography
  pred_wts_sd$disease <- info_packet$ts_disease
  pred_wts_sd$last_obs_time <- i
  pred_wts_sd$h <- 1:h
  pred_wts_sd <- subset(pred_wts_sd, select=c("geography","disease","last_obs_time","h",setdiff(names(pred_wts_sd),c("geography","disease","last_obs_time","h"))))
  
  ## prepare to leave: components
  components$geography <- info_packet$ts_geography
  components$disease <- info_packet$ts_disease
  components$last_obs_time <- i
  components$h <- 1:h
  components <- subset(components, select=c("geography","disease","last_obs_time","h",setdiff(names(components),c("geography","disease","last_obs_time","h"))))
  
  
  ## get outta here
  return(list(output_df = output_df,
              pred_wts = pred_wts,
              components = components,
              features = features,
              pred_wts_sd = pred_wts_sd))
  
}




## input: info_packet, fitted lgb model, forecast horizon
## output: forecasts, lgb weights, features, components

predict_epifforma_theseus <- function(info_packet, h, additional_features = F){
  
  ## handle outliers
  info_packet$ts = handle_outliers(info_packet)
  
  ## distribute zeros
  info_packet$ts = as.numeric(distribute_zeros(info_packet))
  
  ## compute the features
  features <- make_features(info_packet, h)
  
  ## compute the components
  components <- make_components(info_packet, h)

  ## get scaling factors 
  scaling_factor = ifelse(info_packet$ts[length(info_packet$ts)]>=1e-9,info_packet$ts[length(info_packet$ts)],ifelse(max(info_packet$ts)>=1e-9,max(info_packet$ts),1))
  components = components/scaling_factor
  components[components>3]=3

  
  ## predict the class probabilities
  pred_wts_list <- list()
  for(j in 1:length(feature2wt)){
    pred_wts_list[[j]] <- as.matrix(predict(feature2wt[[j]], newdata = cbind(as.matrix(features),as.matrix(components))))
  }
  
  ## average all the fits
  pred_wts <- apply(simplify2array(pred_wts_list), 1, mean)
  pred_wts_sd <- apply(simplify2array(pred_wts_list), 1, sd)
  
  pred_wts[pred_wts>3]=3
  
  ## set epifforma equal to the mean forecast
  fcst_df <- data.frame(h = 1:h,
                        epifforma = pred_wts*scaling_factor,
                        components*scaling_factor)
  
  
  ## melt it
  fcst_df_melt <- melt(fcst_df, id.vars = "h")
  fcst_df_melt$type <- "fcst"
  fcst_df_melt$x <- fcst_df_melt$h + i
  
  ##fcst truth
  fcst_truth <- data.frame(h = 1:h, truth = ts_test$ts[(i+1):(i+h)])
  fcst_df_melt <- merge(fcst_df_melt, fcst_truth, by="h", all.x=T)
  
  ## observation df
  obs_df <- data.frame(h = NA,
                       variable = "obs",
                       value = ts_test$ts[1:(length(info_packet$ts)+h)],
                       type = "obs",
                       x = 1:(length(info_packet$ts)+h),
                       truth = ts_test$ts[1:(length(info_packet$ts)+h)])
  
  ## combine observations and fcst
  output_df <- rbind(obs_df, fcst_df_melt)
  names(output_df) <- c("h","model","fcst","type","x","truth")
  
  ## data for plotting
  output_df$geography <- info_packet$ts_geography
  output_df$disease <- info_packet$ts_disease
  output_df$last_obs_time <- i
  output_df <- data.table(subset(output_df, select=c("geography","disease","last_obs_time","h",setdiff(names(output_df),c("geography","disease","last_obs_time","h")))))
  
  ## prepare to leave: features
  features$geography <- info_packet$ts_geography
  features$disease <- info_packet$ts_disease
  features$last_obs_time <- i
  features <- subset(features, select=c("geography","disease","last_obs_time","h",setdiff(names(features),c("geography","disease","last_obs_time","h"))))
  
  ## prepare to leave: pred_wts
  pred_wts <- data.frame(pred_wts)
  #names(pred_wts) <- paste0(names(components),"_avg")
  pred_wts$geography <- info_packet$ts_geography
  pred_wts$disease <- info_packet$ts_disease
  pred_wts$last_obs_time <- i
  pred_wts$h <- 1:h
  pred_wts <- subset(pred_wts, select=c("geography","disease","last_obs_time","h",setdiff(names(pred_wts),c("geography","disease","last_obs_time","h"))))
  
  ## prepare to leave: pred_wts_sd
  pred_wts_sd <- data.frame(pred_wts_sd)
  #names(pred_wts_sd) <- paste0(names(components),"_sd")
  pred_wts_sd$geography <- info_packet$ts_geography
  pred_wts_sd$disease <- info_packet$ts_disease
  pred_wts_sd$last_obs_time <- i
  pred_wts_sd$h <- 1:h
  pred_wts_sd <- subset(pred_wts_sd, select=c("geography","disease","last_obs_time","h",setdiff(names(pred_wts_sd),c("geography","disease","last_obs_time","h"))))
  
  ## prepare to leave: components
  components$geography <- info_packet$ts_geography
  components$disease <- info_packet$ts_disease
  components$last_obs_time <- i
  components$h <- 1:h
  components <- subset(components, select=c("geography","disease","last_obs_time","h",setdiff(names(components),c("geography","disease","last_obs_time","h"))))
  
  
  ## get outta here
  return(list(output_df = output_df,
              pred_wts = pred_wts,
              components = components,
              features = features,
              pred_wts_sd = pred_wts_sd))
  
}




