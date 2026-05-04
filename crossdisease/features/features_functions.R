#########################################
#########################################
### Functions for Generating Features ###
#########################################
#########################################

##############################################
### Main Function for Time Series Features ###
##############################################

get_features   = function(x,ts, ts_dates, ts_time_cadence, scaling_type = 'lastobs', smoothing_type = 'none', info_packet, smooth_truth = F, bootstrap = F){

  ### Prune and Smooth Time Series ###
  info_packet$ts = info_packet$ts[1:x]
  info_packet$ts = as.numeric(handle_outliers(info_packet)) 
  info_packet$ts = as.numeric(distribute_zeros(info_packet))

  ### Optional Time Series Smoothing
  ts_mod = info_packet$ts
  if(smoothing_type == 'wtmean'){
    ts_smoothed = zoo::rollapply(ts_mod, width = 3, FUN = function(x){weighted.mean(x, w = c(0.1, 0.3, 0.6)[(3-length(as.vector(x))+1):3], na.rm=T)}, partial = T, align = 'right', fill = NA)
  }else if(smoothing_type == 'gam' & info_packet$ts_scale == 'counts'){
    to_fit = unique(sort(pmax(1,c((x-50):x))))
    dat = data.frame(y = round(ts_mod[to_fit]), x = 1:length(to_fit))
    ts_smoothed = zoo::rollapply(ts_mod, width = 3, FUN = function(x){weighted.mean(x, w = c(0.1, 0.3, 0.6)[(3-length(as.vector(x))+1):3], na.rm=T)}, partial = T, align = 'right', fill = NA)
    if(var(dat$y)>0){
      try_gam = try(mgcv::bam(y~s(x, bs = 'tp',k = floor(nrow(dat)/5)), data = dat, family = poisson(link = 'identity'))$fitted.values, silent = T)
      if(class(try_gam)[1] != 'try-error'){
        ts_smoothed[to_fit] = try_gam
      }else{
        try_gam = try(mgcv::bam(y~s(x, bs = 'tp',k = floor(nrow(dat)/5)), data = dat, family = gaussian(link = 'identity'))$fitted.values, silent = T)
        if(class(try_gam)[1] != 'try-error'){
          ts_smoothed[to_fit] = try_gam
        }
      }
    }
  }else if(smoothing_type == 'gam' & (info_packet$ts_scale == 'proportion' | info_packet$ts_scale == 'other')){
    to_fit = unique(sort(pmax(1,c((x-50):x))))
    dat = data.frame(y = ts_mod[to_fit], x = 1:length(to_fit))
    ts_smoothed = zoo::rollapply(ts_mod, width = 3, FUN = function(x){weighted.mean(x, w = c(0.1, 0.3, 0.6)[(3-length(as.vector(x))+1):3], na.rm=T)}, partial = T, align = 'right', fill = NA)
    if(var(dat$y)>0){
      try_gam = try(mgcv::bam(y~s(x, bs = 'tp',k = floor(nrow(dat)/1.5)), data = dat)$fitted.values, silent = T)
      if(class(try_gam)[1] != 'try-error'){
        ts_smoothed[to_fit] = try_gam
      }
    }
  }else{
    ts_smoothed = ts_mod
  }
  
  ### Scale Time Series ###
  scale_vec = ts_smoothed
  if(info_packet$ts_scale == 'counts'){
    scale_vec = round(scale_vec)
  }
  if(scaling_type == 'lastobs'){
    scaling_factor = ifelse(scale_vec[x]>=1e-9,scale_vec[x],1)
    ts_rescaled = ts_smoothed/scaling_factor
  }else if(scaling_type == 'q95'){
    scaling_factor = ifelse(quantile(scale_vec,0.95)>=1e-9, quantile(scale_vec,0.95),1) 
    ts_rescaled = ts_smoothed/scaling_factor
  }else{
    scaling_factor = 1
    ts_rescaled = ts_smoothed/scaling_factor
  }
  
  ### Assign monthly data date of 1st of month
  if(ts_time_cadence == 'monthly'){
    ts_dates = gsub('_','-',paste0(ts_dates,'_1'))
  }
  
  ### Get Feature Set 1
  temp = data.frame(get_featureset1(ts_rescaled[1:x], h=h))
  temp$scaling_factor = scaling_factor
  temp$h = c(1:h)
  temp$truth = ts[(x+1):(x+h)]
  temp$obs = ts[x] 
  temp$obs_smoothed_minus1 = ifelse(info_packet$ts_scale == 'counts', round(ts_smoothed[x-1]),ts_smoothed[x-1])
  temp$obs_smoothed = ifelse(info_packet$ts_scale == 'counts', round(ts_smoothed[x]),ts_smoothed[x])
  temp$last_obs_time = x
  temp$date = as.character(ts_dates[x])
  

  ### Optional generation of smoothed truth, does not impact anything else. 
  if(smooth_truth){
    if(info_packet$ts_scale == 'counts'){
      ts_mod = ts[1:(x+h)]
      to_fit = unique(sort(pmax(1,c((x-50):(x+h)))))
      dat = data.frame(y = round(ts_mod[to_fit]), x = 1:length(to_fit))
      ts_smoothed = zoo::rollapply(ts_mod, width = 3, FUN = function(x){weighted.mean(x, w = c(0.1, 0.3, 0.6)[(3-length(as.vector(x))+1):3], na.rm=T)}, partial = T, align = 'right', fill = NA)
      if(var(dat$y)>0){
        try_gam = try(mgcv::bam(y~s(x, bs = 'tp',k = floor(nrow(dat)/5)), data = dat, family = poisson(link = 'identity'))$fitted.values, silent = T)
        if(class(try_gam)[1] != 'try-error'){
          ts_smoothed[to_fit] = pmax(0,try_gam)
        }else{
          try_gam = try(mgcv::bam(y~s(x, bs = 'tp',k = floor(nrow(dat)/5)), data = dat, family = gaussian(link = 'identity'))$fitted.values, silent = T)
          if(class(try_gam)[1] != 'try-error'){
            ts_smoothed[to_fit] = pmax(0,try_gam)
          }
        }
      }
    }else if(info_packet$ts_scale == 'proportion' | info_packet$ts_scale == 'other'){
      ts_mod = ts[1:(x+h)]
      to_fit = unique(sort(pmax(1,c((x-50):x))))
      dat = data.frame(y = ts_mod[to_fit], x = 1:length(to_fit))
      ts_smoothed = zoo::rollapply(ts_mod, width = 3, FUN = function(x){weighted.mean(x, w = c(0.1, 0.3, 0.6)[(3-length(as.vector(x))+1):3], na.rm=T)}, partial = T, align = 'right', fill = NA)
      if(var(dat$y)>0){
        try_gam = try(mgcv::bam(y~s(x, bs = 'tp',k = floor(nrow(dat)/1.5)), data = dat)$fitted.values, silent = T)
        if(class(try_gam)[1] != 'try-error'){
          ts_smoothed[to_fit] = pmax(0,try_gam)
        }
      }
    }
    temp$truth_sm = tail(ts_smoothed,h)
  }
  

  ### Get Features Set 2
  temp2 = data.frame(get_featureset2(ts_rescaled[1:x],h = h, ts_time_cadence = ts_time_cadence))
  temp = data.frame(temp, subset(temp2, select = -c(h)))

  
  ### Get Feature Set 3
  temp3 =data.frame(get_featureset3(ts_dates = ts_dates[x], ts_dates_forecast = ts_dates[(x+1):(x+h)],
                                    h = h, ts_time_cadence = ts_time_cadence))
  temp = data.frame(temp, subset(temp3, select = -c(h)))

  
  ### Add past day/week/month summaries
  if(ts_time_cadence == 'weekly'){
    dat = data.frame(week_forecast = MMWRweek::MMWRweek(as.Date(ts_dates[1:x], format = '%Y-%m-%d'))$MMWRweek,
                     y = ts_rescaled[1:x])
    AGG = aggregate(y~week_forecast, FUN = median, data = dat)
    names(AGG)[names(AGG)=='y']='past_avg_forecast'
    temp = merge(temp, AGG, by = 'week_forecast', all.x = T, all.y = F)
    temp = temp[order(temp$h),]
    
    dat = data.frame(week = MMWRweek::MMWRweek(as.Date(ts_dates[1:x], format = '%Y-%m-%d'))$MMWRweek,
                     y = ts_rescaled[1:x])
    AGG = aggregate(y~week, FUN = median, data = dat)
    names(AGG)[names(AGG)=='y']='past_avg'
    temp = merge(temp, AGG, by = 'week', all.x = T, all.y = F)
    temp = temp[order(temp$h),]
    
  }else if(ts_time_cadence == 'daily'){
    dat = data.frame(day_of_year_forecast = data.table::yday(as.Date(ts_dates[1:x], format = '%Y-%m-%d')),
                     y = ts_rescaled[1:x])
    AGG = aggregate(y~day_of_year_forecast, FUN = median, data = dat)
    names(AGG)[names(AGG)=='y']='past_avg_forecast'
    temp = merge(temp, AGG, by = 'day_of_year_forecast', all.x = T, all.y = F)
    temp = temp[order(temp$h),]
    
    dat = data.frame(day_of_year = data.table::yday(as.Date(ts_dates[1:x], format = '%Y-%m-%d')),
                     y = ts_rescaled[1:x])
    AGG = aggregate(y~day_of_year, FUN = median, data = dat)
    names(AGG)[names(AGG)=='y']='past_avg'
    temp = merge(temp, AGG, by = 'day_of_year', all.x = T, all.y = F)
    temp = temp[order(temp$h),]
    
  }else if(ts_time_cadence == 'monthly'){
    dat = data.frame(month_forecast = month(as.Date(ts_dates[1:x], format = '%Y-%m-%d')),
                     y = ts_rescaled[1:x])
    AGG = aggregate(y~month_forecast, FUN = median, data = dat)
    names(AGG)[names(AGG)=='y']='past_avg_forecast'
    temp = merge(temp, AGG, by = 'month_forecast', all.x = T, all.y = F)
    temp = temp[order(temp$h),]
    
    dat = data.frame(month = month(as.Date(ts_dates[1:x], format = '%Y-%m-%d')),
                     y = ts_rescaled[1:x])
    AGG = aggregate(y~month, FUN = median, data = dat)
    names(AGG)[names(AGG)=='y']='past_avg'
    temp = merge(temp, AGG, by = 'month', all.x = T, all.y = F)
    temp = temp[order(temp$h),]
  }else{
    temp$past_avg_forecast = mean(ts_rescaled[1:x])
    temp$past_avg = mean(ts_rescaled[1:x])
  }
  temp$past_avg_forecast[is.na(temp$past_avg_forecast)]=1
  temp$past_avg[is.na(temp$past_avg)]=1
  
  temp$is_long = as.numeric(x>=50)
  
  
  ### Some component forecasts
  if(info_packet$ts_time_cadence == 'weekly'){
    freq_decomp = 52
  }else if(info_packet$ts_time_cadence == 'monthly'){
    freq_decomp <- 12
  }else if(info_packet$ts_time_cadence == 'daily'){
    freq_decomp = 365
  }else{
    freq_decomp = 1
  }
  last2ys <- pmax(1,(x-(freq_decomp*2))):x
  temp$meanfcst <- pmin(3,pmax(0, rep(mean(ts_rescaled[last2ys]),h)))
  temp$theta <- pmin(3,pmax(0, forecast(thetaf(ts_rescaled[unique(sort(pmax(1,x-50))):x], h = h), h = h)$mean))
  temp$arima <- pmin(3,pmax(0, forecast(auto.arima(ts_rescaled[unique(sort(pmax(1,x-50))):x]),h = h)$mean))
  
  
  ### Wiggly gam
  to_fit = unique(sort(pmax(1,c((x-50):x))))
  dat = data.frame(y = ts_rescaled[to_fit], x = 1:length(to_fit))
  if(var(dat$y)>0){
    try_gam =  try(predict(mgcv::bam(y~s(x, bs = 'tp',k = floor(nrow(dat)/1.5)), data = dat),newdata = data.frame(x=(x+1):(x+h))), silent = T)
    if(class(try_gam)[1] != 'try-error'){
      pred_gam = try_gam
    }else{
      pred_gam = rep(ts_rescaled[length(ts_rescaled)],h)
    }
  }else{
    pred_gam = rep(ts_rescaled[length(ts_rescaled)],h)
  }
  temp$gam <- pmin(3,pmax(0, pred_gam))
  
  out = tsoutliers(ts_rescaled[1:x], lambda = NULL, iterate = 5)
  temp$is_outlier = as.numeric(x %in% out$index)
  temp$last_outlier_remaining = x - max(out$index)
  

  return(temp)
}





get_features_unscaled  = function(x,ts, ts_dates, ts_time_cadence, scaling_type = 'lastobs', smoothing_type = 'none', info_packet, smooth_truth = F, bootstrap = F){
  
  ### Prune and Smooth Time Series ###
  info_packet$ts = info_packet$ts[1:x]
  info_packet$ts = as.numeric(handle_outliers(info_packet)) 
  info_packet$ts = as.numeric(distribute_zeros(info_packet))
  

  ### Optional Time Series Smoothing
  ts_mod = info_packet$ts
  if(smoothing_type == 'wtmean'){
    ts_smoothed = zoo::rollapply(ts_mod, width = 3, FUN = function(x){weighted.mean(x, w = c(0.1, 0.3, 0.6)[(3-length(as.vector(x))+1):3], na.rm=T)}, partial = T, align = 'right', fill = NA)
  }else if(smoothing_type == 'gam' & info_packet$ts_scale == 'counts'){
    to_fit = unique(sort(pmax(1,c((x-50):x))))
    dat = data.frame(y = round(ts_mod[to_fit]), x = 1:length(to_fit))
    ts_smoothed = zoo::rollapply(ts_mod, width = 3, FUN = function(x){weighted.mean(x, w = c(0.1, 0.3, 0.6)[(3-length(as.vector(x))+1):3], na.rm=T)}, partial = T, align = 'right', fill = NA)
    if(var(dat$y)>0){
      try_gam = try(mgcv::bam(y~s(x, bs = 'tp',k = floor(nrow(dat)/5)), data = dat, family = poisson(link = 'identity'))$fitted.values, silent = T)
      if(class(try_gam)[1] != 'try-error'){
        ts_smoothed[to_fit] = try_gam
      }else{
        try_gam = try(mgcv::bam(y~s(x, bs = 'tp',k = floor(nrow(dat)/5)), data = dat, family = gaussian(link = 'identity'))$fitted.values, silent = T)
        if(class(try_gam)[1] != 'try-error'){
          ts_smoothed[to_fit] = try_gam
        }
      }
    }
  }else if(smoothing_type == 'gam' & (info_packet$ts_scale == 'proportion' | info_packet$ts_scale == 'other')){
    to_fit = unique(sort(pmax(1,c((x-50):x))))
    dat = data.frame(y = ts_mod[to_fit], x = 1:length(to_fit))
    ts_smoothed = zoo::rollapply(ts_mod, width = 3, FUN = function(x){weighted.mean(x, w = c(0.1, 0.3, 0.6)[(3-length(as.vector(x))+1):3], na.rm=T)}, partial = T, align = 'right', fill = NA)
    if(var(dat$y)>0){
      try_gam = try(mgcv::bam(y~s(x, bs = 'tp',k = floor(nrow(dat)/1.5)), data = dat)$fitted.values, silent = T)
      if(class(try_gam)[1] != 'try-error'){
        ts_smoothed[to_fit] = try_gam
      }
    }
  }else{
    ts_smoothed = ts_mod
  }
  
  ### Scale Time Series ###
  scale_vec = ts_smoothed
  if(info_packet$ts_scale == 'counts'){
    scale_vec = round(scale_vec)
  }
  if(scaling_type == 'lastobs'){
    scaling_factor = ifelse(scale_vec[x]>=1e-9,scale_vec[x],1)
    ts_rescaled = ts_smoothed/scaling_factor
  }else if(scaling_type == 'q95'){
    scaling_factor = ifelse(quantile(scale_vec,0.95)>=1e-9, quantile(scale_vec,0.95),1) 
    ts_rescaled = ts_smoothed/scaling_factor
  }else{
    scaling_factor = 1
    ts_rescaled = ts_smoothed/scaling_factor
  }
  
  ### Assign monthly data date of 1st of month
  if(ts_time_cadence == 'monthly'){
    ts_dates = gsub('_','-',paste0(ts_dates,'_1'))
  }
  
  ### Get Feature Set 1
  temp = data.frame(get_featureset1(ts_rescaled[1:x], h=h))
  temp$scaling_factor = scaling_factor
  temp$h = c(1:h)
  temp$truth = ts[(x+1):(x+h)]
  temp$obs = ts[x] 
  temp$obs_smoothed_minus1 = ifelse(info_packet$ts_scale == 'counts', round(ts_smoothed[x-1]),ts_smoothed[x-1])
  temp$obs_smoothed = ifelse(info_packet$ts_scale == 'counts', round(ts_smoothed[x]),ts_smoothed[x])
  temp$last_obs_time = x
  temp$date = as.character(ts_dates[x])
  

  ### Optional generation of smoothed truth, does not impact anything else. 
  if(smooth_truth){
    if(info_packet$ts_scale == 'counts'){
      ts_mod = ts[1:(x+h)]
      to_fit = unique(sort(pmax(1,c((x-50):(x+h)))))
      dat = data.frame(y = round(ts_mod[to_fit]), x = 1:length(to_fit))
      ts_smoothed = zoo::rollapply(ts_mod, width = 3, FUN = function(x){weighted.mean(x, w = c(0.1, 0.3, 0.6)[(3-length(as.vector(x))+1):3], na.rm=T)}, partial = T, align = 'right', fill = NA)
      if(var(dat$y)>0){
        try_gam = try(mgcv::bam(y~s(x, bs = 'tp',k = floor(nrow(dat)/5)), data = dat, family = poisson(link = 'identity'))$fitted.values, silent = T)
        if(class(try_gam)[1] != 'try-error'){
          ts_smoothed[to_fit] = pmax(0,try_gam)
        }else{
          try_gam = try(mgcv::bam(y~s(x, bs = 'tp',k = floor(nrow(dat)/5)), data = dat, family = gaussian(link = 'identity'))$fitted.values, silent = T)
          if(class(try_gam)[1] != 'try-error'){
            ts_smoothed[to_fit] = pmax(0,try_gam)
          }
        }
      }
    }else if(info_packet$ts_scale == 'proportion' | info_packet$ts_scale == 'other'){
      ts_mod = ts[1:(x+h)]
      to_fit = unique(sort(pmax(1,c((x-50):x))))
      dat = data.frame(y = ts_mod[to_fit], x = 1:length(to_fit))
      ts_smoothed = zoo::rollapply(ts_mod, width = 3, FUN = function(x){weighted.mean(x, w = c(0.1, 0.3, 0.6)[(3-length(as.vector(x))+1):3], na.rm=T)}, partial = T, align = 'right', fill = NA)
      if(var(dat$y)>0){
        try_gam = try(mgcv::bam(y~s(x, bs = 'tp',k = floor(nrow(dat)/1.5)), data = dat)$fitted.values, silent = T)
        if(class(try_gam)[1] != 'try-error'){
          ts_smoothed[to_fit] = pmax(0,try_gam)
        }
      }
    }
    temp$truth_sm = tail(ts_smoothed,h)
  }
  

  ### Get Feature Set 2
  temp2 = data.frame(get_featureset2(ts_rescaled[1:x],h = h, ts_time_cadence = ts_time_cadence))
  temp = data.frame(temp, subset(temp2, select = -c(h)))

  ### Get Feature Set 3
  temp3 =data.frame(get_featureset3(ts_dates = ts_dates[x], ts_dates_forecast = ts_dates[(x+1):(x+h)],
                                    h = h,
                                    ts_time_cadence = ts_time_cadence))
  temp = data.frame(temp, subset(temp3, select = -c(h)))

  
  ### Add past day/week/month summaries
  if(ts_time_cadence == 'weekly'){
    dat = data.frame(week_forecast = MMWRweek::MMWRweek(as.Date(ts_dates[1:x], format = '%Y-%m-%d'))$MMWRweek,
                     y = ts_rescaled[1:x])
    AGG = aggregate(y~week_forecast, FUN = median, data = dat)
    names(AGG)[names(AGG)=='y']='past_avg_forecast'
    temp = merge(temp, AGG, by = 'week_forecast', all.x = T, all.y = F)
    temp = temp[order(temp$h),]
    
    dat = data.frame(week = MMWRweek::MMWRweek(as.Date(ts_dates[1:x], format = '%Y-%m-%d'))$MMWRweek,
                     y = ts_rescaled[1:x])
    AGG = aggregate(y~week, FUN = median, data = dat)
    names(AGG)[names(AGG)=='y']='past_avg'
    temp = merge(temp, AGG, by = 'week', all.x = T, all.y = F)
    temp = temp[order(temp$h),]
    
  }else if(ts_time_cadence == 'daily'){
    dat = data.frame(day_of_year_forecast = data.table::yday(as.Date(ts_dates[1:x], format = '%Y-%m-%d')),
                     y = ts_rescaled[1:x])
    AGG = aggregate(y~day_of_year_forecast, FUN = median, data = dat)
    names(AGG)[names(AGG)=='y']='past_avg_forecast'
    temp = merge(temp, AGG, by = 'day_of_year_forecast', all.x = T, all.y = F)
    temp = temp[order(temp$h),]
    
    dat = data.frame(day_of_year = data.table::yday(as.Date(ts_dates[1:x], format = '%Y-%m-%d')),
                     y = ts_rescaled[1:x])
    AGG = aggregate(y~day_of_year, FUN = median, data = dat)
    names(AGG)[names(AGG)=='y']='past_avg'
    temp = merge(temp, AGG, by = 'day_of_year', all.x = T, all.y = F)
    temp = temp[order(temp$h),]
    
  }else if(ts_time_cadence == 'monthly'){
    dat = data.frame(month_forecast = month(as.Date(ts_dates[1:x], format = '%Y-%m-%d')),
                     y = ts_rescaled[1:x])
    AGG = aggregate(y~month_forecast, FUN = median, data = dat)
    names(AGG)[names(AGG)=='y']='past_avg_forecast'
    temp = merge(temp, AGG, by = 'month_forecast', all.x = T, all.y = F)
    temp = temp[order(temp$h),]
    
    dat = data.frame(month = month(as.Date(ts_dates[1:x], format = '%Y-%m-%d')),
                     y = ts_rescaled[1:x])
    AGG = aggregate(y~month, FUN = median, data = dat)
    names(AGG)[names(AGG)=='y']='past_avg'
    temp = merge(temp, AGG, by = 'month', all.x = T, all.y = F)
    temp = temp[order(temp$h),]
  }else{
    temp$past_avg_forecast = mean(ts_rescaled[1:x])
    temp$past_avg = mean(ts_rescaled[1:x])
  }
  temp$past_avg_forecast[is.na(temp$past_avg_forecast)]=1
  temp$past_avg[is.na(temp$past_avg)]=1
  
  temp$is_long = as.numeric(x>=50)
  
  
  ### Some component forecasts
  if(info_packet$ts_time_cadence == 'weekly'){
    freq_decomp = 52
  }else if(info_packet$ts_time_cadence == 'monthly'){
    freq_decomp <- 12
  }else if(info_packet$ts_time_cadence == 'daily'){
    freq_decomp = 365
  }else{
    freq_decomp = 1
  }
  PMAX = Inf
  last2ys <- pmax(1,(x-(freq_decomp*2))):x
  temp$meanfcst <- pmin(PMAX,pmax(0, rep(mean(ts_rescaled[last2ys]),h)))
  temp$theta <- pmin(PMAX,pmax(0, forecast(thetaf(ts_rescaled[unique(sort(pmax(1,x-50))):x], h = h), h = h)$mean))
  temp$arima <- pmin(PMAX,pmax(0, forecast(auto.arima(ts_rescaled[unique(sort(pmax(1,x-50))):x]),h = h)$mean))
  
  
  ### Wiggly gam
  to_fit = unique(sort(pmax(1,c((x-50):x))))
  dat = data.frame(y = ts_rescaled[to_fit], x = 1:length(to_fit))
  if(var(dat$y)>0){
    try_gam =  try(predict(mgcv::bam(y~s(x, bs = 'tp',k = floor(nrow(dat)/1.5)), data = dat),newdata = data.frame(x=(x+1):(x+h))), silent = T)
    if(class(try_gam)[1] != 'try-error'){
      pred_gam = try_gam
    }else{
      pred_gam = rep(ts_rescaled[length(ts_rescaled)],h)
    }
  }else{
    pred_gam = rep(ts_rescaled[length(ts_rescaled)],h)
  }
  temp$gam <- pmin(PMAX,pmax(0, pred_gam))
  
  out = tsoutliers(ts_rescaled[1:x], lambda = NULL, iterate = 5)
  temp$is_outlier = as.numeric(x %in% out$index)
  temp$last_outlier_remaining = x - max(out$index)
  
  return(temp)
}








##############################
### Define Feature Subsets ###
##############################

get_featureset1 <- function(ts, h=4, is.cases = T){
  
  ## length of time series
  n <- length(ts)
  
  ## get 4 week rolling average ts
  ts_roll <- c(rep(NA,3),zoo::rollmean(ts,4,align = "right"))
  
  ## compute difference of the rolling average
  tsdiff_roll <- c(NA,diff(ts_roll))
  
  ## ts difference
  tsdiff <- diff(ts)
  ndiff <- length(tsdiff)
  
  ## get running means
  tdf <- data.frame(## total summaries
    roll_min  = min(ts_roll, na.rm=T),
    roll_q.10 = as.numeric(quantile(ts_roll, probs = .10, na.rm=T)),
    roll_q.25 = as.numeric(quantile(ts_roll, probs = .25, na.rm=T)),
    roll_q.50 = as.numeric(quantile(ts_roll, probs = .50, na.rm=T)),
    roll_q.75 = as.numeric(quantile(ts_roll, probs = .75, na.rm=T)),
    roll_q.90 = as.numeric(quantile(ts_roll, probs = .90, na.rm=T)),
    roll_max  = max(ts_roll, na.rm=T),
    ## means
    m0 = ts[n],
    m1 = mean(ts[pmax(1,(n-1)):n], na.rm=T),
    m3 = mean(ts[pmax(1,(n-3)):n], na.rm=T),
    m9 = mean(ts[pmax(1,(n-9)):n], na.rm=T),
    m19 = mean(ts[pmax(1,n-19):n], na.rm=T),
    mtotal = mean(ts, na.rm=T),
    #
    rollm0 = ts_roll[n],
    rollm1 = mean(ts_roll[pmax(1,(n-1)):n], na.rm=T),
    rollm3 = mean(ts_roll[pmax(1,(n-3)):n], na.rm=T),
    rollm9 = mean(ts_roll[pmax(1,(n-9)):n], na.rm=T),
    rollm19 = mean(ts_roll[pmax(1,n-19):n], na.rm=T),
    rollmtotal = mean(ts_roll, na.rm=T),
    #
    diffrollm0 = tsdiff_roll[n],
    diffrollm1 = mean(tsdiff_roll[pmax(1,(n-1)):n], na.rm=T),
    diffrollm3 = mean(tsdiff_roll[pmax(1,(n-3)):n], na.rm=T),
    diffrollm9 = mean(tsdiff_roll[pmax(1,(n-9)):n], na.rm=T),
    diffrollm19 = mean(tsdiff_roll[pmax(1,n-19):n], na.rm=T),
    diffrollmtotal = mean(tsdiff_roll, na.rm=T),
    ## slopes (LAUREN CHANGED THE BELOW FROM == 0 on 12/30/24)
    sloperoll01 = as.numeric(ifelse(ts_roll[n-1] < 0.001, ts_roll[n], ts_roll[n]/ts_roll[n-1])),
    sloperoll02 = as.numeric(ifelse(ts_roll[n-2] < 0.001, ts_roll[n], ts_roll[n]/ts_roll[n-2])),
    sloperoll12 = as.numeric(ifelse(ts_roll[n-2] < 0.001, ts_roll[n-1], ts_roll[n-1]/ts_roll[n-2])),
    ## local modes
    loc2   = ifelse(max(ts[pmax(1,n-2):n], na.rm=T) < 0.001,  .5, ts[n]/max(ts[pmax(1,n-2):n], na.rm=T)),
    loc5   = ifelse(max(ts[pmax(1,n-5):n], na.rm=T) < 0.001,  .5, ts[n]/max(ts[pmax(1,n-5):n], na.rm=T)),
    loc10  = ifelse(max(ts[pmax(1,n-10):n], na.rm=T) < 0.001, .5, ts[n]/max(ts[pmax(1,n-10):n], na.rm=T)),
    loc15  = ifelse(max(ts[pmax(1,n-15):n], na.rm=T) < 0.001, .5, ts[n]/max(ts[pmax(1,n-15):n], na.rm=T)),
    loctotalmax = ifelse(max(ts, na.rm=T) < 0.001, .5, ts[n]/max(ts, na.rm=T)),
    loctotalmin = ifelse(ts[n] < 0.001, .5, min(ts, na.rm=T)/ts[n]),
    ## prop zero
    propzero5 = mean(ts[pmax(1,n-5):n] == 0, na.rm=T),
    propzero10 = mean(ts[pmax(1,n-10):n] == 0, na.rm=T),
    propzero15 = mean(ts[pmax(1,n-15):n] == 0, na.rm=T),
    propzerototal = mean(ts == 0, na.rm=T),
    ## runs
    runs1 = sum(sign(rev(tsdiff))[1:pmin(ndiff,1)]),
    runs2 = sum(sign(rev(tsdiff))[1:pmin(ndiff,2)]),
    runs3 = sum(sign(rev(tsdiff))[1:pmin(ndiff,3)]),
    runs4 = sum(sign(rev(tsdiff))[1:pmin(ndiff,4)]),
    runs5 = sum(sign(rev(tsdiff))[1:pmin(ndiff,5)]),
    runs8 = sum(sign(rev(tsdiff))[1:pmin(ndiff,8)]),
    runs10 = sum(sign(rev(tsdiff))[1:pmin(ndiff,10)]))
  
  
  ttdf <- NULL
  for(i in 1:h){
    ttdf <- rbind(ttdf, tdf)
  }
  ## slope signs
  ttdf$slopesign01 <- sign(ttdf$sloperoll01-1)
  ttdf$slopesign02 <- sign(ttdf$sloperoll02-1)
  ttdf$slopesign12 <- sign(ttdf$sloperoll12-1)
  
  ## mean differences
  ttdf$m0_minus_m1 <- ttdf$m0 - ttdf$m1
  ttdf$m0_minus_m3 <- ttdf$m0 - ttdf$m3
  ttdf$m0_minus_m9 <- ttdf$m0 - ttdf$m9
  ttdf$m0_minus_m19 <- ttdf$m0 - ttdf$m19
  
  ## get outta here
  return(ttdf)
}




get_featureset2 <- function(ts, h=4, ts_time_cadence){
  
  ## make everything at least as big as 1e-10
  minval <- 1e-10

  ## make the output
  tsf <- data.frame(h = 1:h)

  ## get 4 week rolling average ts
  ts_smooth <- c(zoo::rollapply(ts,4,align = "right", partial = T, FUN = mean)) 
  
  
  #####################
  ### LOCAL METRICS ###
  #####################
  last10id <- (length(ts_smooth)-9):length(ts_smooth)
  
  #### Goal: Is the last jump positive or negative, and how much, multiplicatively? 
  ## ratio of gam[t] / gam[t-1]
  if(abs(ts_smooth[length(ts_smooth)-1]) <= minval | abs(ts[length(ts)-1]) <= minval |
     abs(ts_smooth[length(ts_smooth)]) <= minval | abs(ts[length(ts)]) <= minval){
    tsf$gr12_div_23 = 0 #pretty flat, possibly just getting started or just ending
  }else{
    tsf$gr12_div_23 = tanh(( abs(ts_smooth[length(ts_smooth)] + ts_smooth[length(ts_smooth)-1])/abs(ts_smooth[length(ts_smooth)-1]+ts_smooth[length(ts_smooth)-2])) - 1)
    
  }
  
  #### Goal: How big is the value relative to the values that have been observed recently?
  ts_smooth_smooth = zoo::rollapply(ts_smooth, align = 'center', width = 3, FUN = function(x){mean(x,na.rm=T)}, partial = T)
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
    tsf$coefvar <- tanh(-0.1) 
  }

  ##########################
  ### Global-ish Metrics ### (Last 2 Years)
  ##########################

  if(ts_time_cadence == 'weekly'){
    last2ys <- pmax(1,(length(ts_smooth)-(52*2))):length(ts_smooth)
  }else if(ts_time_cadence == 'monthly'){
    last2ys <- pmax(1,(length(ts)-(12*2))):length(ts)
  }else if(ts_time_cadence == 'daily'){
    last2ys <- pmax(1,(length(ts_smooth)-(365*2))):length(ts_smooth)
  }else{
    last2ys <- 1:length(ts_smooth)
  }
  
  #### Goal: How big is the value, relative to the mean of past values, subset num to last to to avoid length of ts effect
  ## mean(y_1:t)/mean(y_1:(t-1)), smooth and normalized
  if(mean(abs(ts_smooth[last2ys[-length(last2ys)]]))>minval){
    tsf$avg_recent_div_avg_global <- tanh((mean(abs(ts_smooth[last10id]))/mean(abs(ts_smooth[last2ys[-length(last2ys)]]))) - 1) 
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
  
  if(var(ts_smooth[-length(ts_smooth)]) == 0){
    tsf$cur_relative_max = 1 
  }else{
    tsf$cur_relative_max = tail(ts_smooth,1)/max(ts_smooth[-length(ts_smooth)])
  }
  
  
  ## Consecutive increase divided by max consecutive increase
  ts_smooth_diff <- as.numeric(diff(ts_smooth)>0)
  if(ts_smooth_diff[length(ts_smooth_diff)]==1){
    first = min(which(rev(ts_smooth_diff) == 0)) #first different value from end
  }
  MAT = data.frame(lengths = rle(ts_smooth_diff == 1)$lengths, values = rle(ts_smooth_diff == 1)$values)
  MAT = MAT[MAT$values,]
  tsf$relative_increases = ifelse(ts_smooth_diff[length(ts_smooth_diff)] == 0, 0, median( MAT$lengths<(first-1))) 
  

  ### Proportion of Year Since Average Yearly Max
  if(ts_time_cadence == 'weekly'){
    freq_decomp = 52
  }else if(ts_time_cadence == 'monthly'){
    freq_decomp = 12 
  }else if(ts_time_cadence == 'daily'){
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
      max_lag = pmin(ceiling(freq_decomp*1.25), length(ts_sub)-1) 
      min_lag = floor(freq_decomp*0.75) 
      lags = 0:max_lag 
      tsf$seasonality = pmax(0,max(as.numeric(unlist(acf(ts_sub,max_lag,plot=F)$acf[lags >= min_lag]))))
    }else{
      tsf$seasonality = 0
    }
  }

  
  ### Seasonality with spearman
  spearman_acf <- function(x, max_lag) {
    acf_values <- numeric(max_lag + 1)  # Include lag 0
    
    for (lag in 0:max_lag) {
      # Compute Spearman rank correlation for the given lag
      if (lag == 0) {
        acf_values[lag + 1] <- 1  # The correlation at lag 0 is always 1 (perfect correlation with itself)
      } else {
        acf_values[lag + 1] <- cor(x[(1 + lag):length(x)], x[1:(length(x) - lag)], method = "spearman")
      }
    }
    
    return(acf_values)
  }
  if(sum(ts > minval) == 0){ #if all zeros, metric is zero
    tsf$seasonality_spearman = 0
  }else{
    if(min(which(ts > minval)) == 1){
      ts_sub = ts
    }else{
      ts_sub = ts[-c(1:(min(which(ts > minval))-1))] #zeros at the beginning of ts confuse this metric
    }
    if(length(ts_sub) >= (freq_decomp+1)){ #require at least freq_decomp times to calculate this metric, otherwise zero
      max_lag = pmin(ceiling(freq_decomp*1.25), length(ts_sub)-1) 
      min_lag = floor(freq_decomp*0.75) 
      lags = 0:max_lag 
      tsf$seasonality_spearman = pmax(0,max(as.numeric(unlist(spearman_acf(ts_sub,max_lag)[lags >= min_lag]))))
    }else{
      tsf$seasonality_spearman = 0
    }
  }
  ## get outta here
  return(tsf)
  
}



get_holiday_lag = function(ts_dates, Holiday){
  date_cur = as.Date(holiday(year = year(as.Date(ts_dates, format = '%Y-%m-%d')), Holiday = Holiday))
  date_next = as.Date(holiday(year = year(as.Date(ts_dates, format = '%Y-%m-%d'))+1, Holiday = Holiday))
  diff_cur = round(as.numeric(difftime(date_cur, ts_dates, units = 'days')))
  diff_next = round(as.numeric(difftime(date_next, ts_dates, units = 'days')))
  diff_return = ifelse(diff_cur>0,diff_cur,diff_next)
  return(diff_return)
}

get_featureset3 <- function(ts_dates, ts_dates_forecast, h=4, ts_time_cadence){
  
  
  ## make the output
  tsf <- data.frame(h = 1:h)
  
  
  if(ts_time_cadence == 'weekly'){
    tsf$month = month(as.Date(ts_dates, format = '%Y-%m-%d'))
    tsf$week = MMWRweek::MMWRweek(as.Date(ts_dates, format = '%Y-%m-%d'))$MMWRweek
    tsf$day_of_week = NA
    tsf$day_of_year = NA
    tsf$month_forecast = month(as.Date(ts_dates_forecast, format = '%Y-%m-%d'))
    tsf$week_forecast = MMWRweek::MMWRweek(as.Date(ts_dates_forecast, format = '%Y-%m-%d'))$MMWRweek
    tsf$day_of_week_forecast = NA
    tsf$day_of_year_forecast = NA
    tsf$days_until_xmas = NA
    tsf$days_until_easter = NA
    tsf$days_until_newyears = NA
    tsf$days_until_usthanksgiving = NA
    tsf$weeks_until_xmas = floor(get_holiday_lag(ts_dates_forecast, 'ChristmasDay')/7)
    tsf$weeks_until_easter = floor(get_holiday_lag(ts_dates_forecast, 'Easter')/7)
    tsf$weeks_until_newyears = floor(get_holiday_lag(ts_dates_forecast, 'NewYearsDay')/7)
    tsf$weeks_until_usthanksgiving = floor(get_holiday_lag(ts_dates_forecast, 'USThanksgivingDay')/7)
  }else if(ts_time_cadence == 'monthly'){
    tsf$month = month(as.Date(ts_dates, format = '%Y-%m-%d'))
    tsf$week = NA
    tsf$day_of_week = NA
    tsf$day_of_year = NA
    tsf$month_forecast = month(as.Date(ts_dates_forecast, format = '%Y-%m-%d'))
    tsf$week_forecast = NA
    tsf$day_of_week_forecast = NA
    tsf$day_of_year_forecast = NA
    tsf$days_until_xmas = NA
    tsf$days_until_easter = NA
    tsf$days_until_newyears = NA
    tsf$days_until_usthanksgiving = NA
    tsf$weeks_until_xmas = NA
    tsf$weeks_until_easter = NA
    tsf$weeks_until_newyears = NA
    tsf$weeks_until_usthanksgiving = NA
  }else if(ts_time_cadence == 'daily'){
    tsf$month = month(as.Date(ts_dates, format = '%Y-%m-%d'))
    tsf$week = MMWRweek::MMWRweek(as.Date(ts_dates, format = '%Y-%m-%d'))$MMWRweek
    tsf$day_of_week = weekdays(as.Date(ts_dates, format = '%Y-%m-%d'))
    tsf$day_of_year = data.table::yday(as.Date(ts_dates, format = '%Y-%m-%d'))
    tsf$month_forecast = month(as.Date(ts_dates_forecast, format = '%Y-%m-%d'))
    tsf$week_forecast = MMWRweek::MMWRweek(as.Date(ts_dates_forecast, format = '%Y-%m-%d'))$MMWRweek
    tsf$day_of_week_forecast = weekdays(as.Date(ts_dates_forecast, format = '%Y-%m-%d'))
    tsf$day_of_year_forecast = data.table::yday(as.Date(ts_dates_forecast, format = '%Y-%m-%d'))
    tsf$days_until_xmas = get_holiday_lag(ts_dates_forecast, 'ChristmasDay')
    tsf$days_until_easter = get_holiday_lag(ts_dates_forecast, 'Easter')
    tsf$days_until_newyears = get_holiday_lag(ts_dates_forecast, 'NewYearsDay')
    tsf$days_until_usthanksgiving = get_holiday_lag(ts_dates_forecast, 'USThanksgivingDay')
    tsf$weeks_until_xmas = NA
    tsf$weeks_until_easter = NA
    tsf$weeks_until_newyears = NA
    tsf$weeks_until_usthanksgiving = NA
  }else{
    tsf$month = NA
    tsf$week =NA
    tsf$day_of_week = NA
    tsf$day_of_year = NA
    tsf$month_forecast = NA
    tsf$week_forecast = NA
    tsf$day_of_week_forecast = NA
    tsf$day_of_year_forecast = NA
    tsf$day_of_week_forecast = NA
    tsf$days_until_xmas = NA
    tsf$days_until_easter = NA
    tsf$days_until_newyears = NA
    tsf$days_until_usthanksgiving = NA
    tsf$weeks_until_xmas = NA
    tsf$weeks_until_easter = NA
    tsf$weeks_until_newyears = NA
    tsf$weeks_until_usthanksgiving = NA
  }
  
  ## get outta here
  return(tsf)
  
}






### This function screens for outliers in a time series.
## handle outliers. Call obs an outlier if
## outlier in at least two metrics AND
## not in interval of multiple outliers unless value is zero AND
## not last value unless last value is zero or > 5 x max in last 2 years
handle_outliers <- function(info_packet){
  ts = info_packet$ts
  if(max(ts[-length(ts)]) > 0){
    
    
    if(info_packet$ts_time_cadence == 'weekly'){
      last2ys <- pmax(1,(length(ts)-(52*2))):length(ts)
    }else if(info_packet$ts_time_cadence == 'daily'){
      last2ys <- pmax(1,(length(ts)-(365*2))):length(ts)
    }else if(info_packet$ts_time_cadence == 'monthly'){
      last2ys <- pmax(1,(length(ts)-(12*2))):length(ts)
    }else{
      last2ys <- 1:length(ts)
    }
    
    
    ### Initial Outlier Screening, Excluding last value
    out = tsoutliers(ts, lambda = NULL, iterate = 5)
    inds = out$index
    inds = inds[ts[inds] != 0] #ignore the zeros here 
    if(length(inds)>1){
      not_outlier = inds[apply(cbind(inds),1,FUN = function(x,inds){min(abs(x-inds[inds!=x]))},inds= inds)==1]
      inds = inds[!(inds %in% not_outlier) | inds == length(ts)]
    }
    A = ts[last2ys[-length(last2ys)]]
    ROLLMAX = max(zoo::rollmedian(A,3))
    inds2 = which(ts > (ROLLMAX*3))
    inds = sort(unique(c(inds, inds2)))
    inds = inds[inds != length(ts)] #exclude last value for now
    if(length(inds) >= 1 & var(ts)>0){
      ts_clean = ts
      ts_clean[inds] = NA
      ts_clean = forecast::na.interp(ts_clean)
      ts = as.numeric(ts_clean)
    }
    
    ### Second Pass Outlier Screening, Focusing on last value
    out = tsoutliers(ts, lambda = NULL, iterate = 5)
    inds = out$index
    if(length(inds)>1){
      not_outlier = inds[apply(cbind(inds),1,FUN = function(x,inds){min(abs(x-inds[inds!=x]))},inds= inds)==1]
      inds = inds[!(inds %in% not_outlier) | inds == length(ts)]
    }
    A = ts[last2ys[-length(last2ys)]]
    ROLLMAX = max(zoo::rollmedian(A,3))
    ROLLMAXDIFF = max(zoo::rollmedian(abs(diff(A)),3))
    if(ROLLMAX > 0 & ts[length(ts)-1] > 0 & ts[length(ts)]>0 & ROLLMAXDIFF > 0){
      if((ts[length(ts)]/ROLLMAX) > 3  | (abs(diff(tail(ts,2)))/ROLLMAXDIFF) > 3){
        ts_clean = ts
        ts_clean[length(ts)] = NA
        ts_clean = forecast::na.interp(ts_clean)
        ts = as.numeric(ts_clean)
      }
    }

    
    ### Deal with outlier zeros
    if(info_packet$ts_scale == 'counts'){
      n = length(ts)
      tail_of_ts = unique(sort(pmax(1,(n-50):n)))
      smooth_df <- data.frame(x = tail_of_ts,y = round(as.numeric(ts)[tail_of_ts]))
      K = pmax(round(length(as.numeric(smooth_df$y))/10),5)
      MAXIT = 50
      gam_mod <- try(mgcv::bam(y ~ s(x, k = K, bs = "tp", m = c(2,1)), data = smooth_df, method = "fREML", gamma = 10, discrete = T, select = T, family = 'poisson', control = list(maxit = MAXIT)), silent = T)
      if(class(gam_mod)[1]!='try-error'){
        probs = apply(cbind(predict(gam_mod, type = 'response')), 1, FUN = function(x){ppois(0,lambda = x, lower.tail = TRUE)})
        inds = which(probs < 0.05 & smooth_df$y == 0)
        if(length(inds)>1){
          inds = inds + (min(smooth_df$x)-1) #correct for shortened ts, helps with gam speed
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
    ### Deal with outlier zeros
    if(info_packet$ts_scale == 'proportion'){
      n = length(ts)
      tail_of_ts = unique(sort(pmax(1,(n-50):n)))
      smooth_df <- data.frame(x = tail_of_ts,y = as.numeric(ts)[tail_of_ts])
      K = pmax(round(length(as.numeric(smooth_df$y))/10),5)
      MAXIT = 50
      gam_mod <- try(mgcv::bam(y ~ s(x, k = K, bs = "tp", m = c(2,1)), data = smooth_df, method = "fREML", gamma = 10, discrete = T, select = T, family = 'gaussian', control = list(maxit = MAXIT)), silent = T)
      if(class(gam_mod)[1]!='try-error'){
        probs = apply(cbind(predict(gam_mod, type = 'response')), 1, FUN = function(x, sigma){pnorm(0,mean = x, sd = sigma, lower.tail = TRUE)}, sigma = sqrt(summary(gam_mod)$dispersion))
        inds = which(probs < 0.05 & smooth_df$y == 0)
        if(length(inds)>1){
          inds = inds + (min(smooth_df$x)-1) #correct for shortened ts, helps with gam speed
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


### This function distributes errant zeros in the observed time series. 
distribute_zeros = function(info_packet){
  inds = which(info_packet$ts == 0)
  inds = inds[inds>min(which(info_packet$ts>0))]
  ts = info_packet$ts
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


