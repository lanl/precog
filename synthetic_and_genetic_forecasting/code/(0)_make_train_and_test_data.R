## Dave Osthus
## 2-11-26
## Make the training data 

## define libraries
library(data.table)
library(ggplot2)
library(lubridate)
library(zoo)
library(here)

## -----------------------------------------------------------------------------
## define paths
rawdatapath <- paste0(here::here("synthetic_and_genetic_forecasting", "data_raw"), "/") 
cleandatapath <- paste0(here::here("synthetic_and_genetic_forecasting", "data_clean"), "/") 
figpath <- paste0(here::here("synthetic_and_genetic_forecasting", "figs"), "/") 


## extract the VACs
bin_counts <- function(DT, nbins, time_col = "time", group_col = "group",
                       tmin = NULL, tmax = NULL, right_closed = TRUE,
                       alpha = 0) {
  
  DT <- as.data.table(DT)
  alpha <- as.numeric(alpha)
  if (!is.finite(alpha) || alpha < 0 || alpha > 1) stop("`alpha` must be between 0 and 1.")
  nbins <- as.integer(nbins)
  if (!is.finite(nbins) || nbins < 1L) stop("`nbins` must be an integer >= 1.")
  
  n_total <- nrow(DT)
  thresh  <- alpha * n_total
  
  # --- keep groups by total count threshold (based on ORIGINAL DT) ---
  g_counts <- DT[, .(N = .N), by = group_col]
  setnames(g_counts, group_col, "group_tmp")
  keep_groups <- unique(g_counts[N > thresh, group_tmp])
  
  if (!length(keep_groups)) {
    # still return overall (group=0) if nothing kept
    return(data.table(group = 0, bin = seq_len(nbins), n = 0L)[, .(group, n, bin)])
  }
  
  DTk <- DT[get(group_col) %in% keep_groups]
  
  # --- breaks (range after filtering unless provided) ---
  if (is.null(tmin)) tmin <- min(DTk[[time_col]], na.rm = TRUE)
  if (is.null(tmax)) tmax <- max(DTk[[time_col]], na.rm = TRUE)
  brks <- seq(tmin, tmax, length.out = nbins + 1L)
  
  # --- bin assignment helper ---
  assign_bin <- function(d) {
    d[, bin := findInterval(get(time_col), brks,
                            rightmost.closed = TRUE,
                            left.open = !right_closed)]
    d[bin < 1L, bin := 1L]
    d[bin > nbins, bin := nbins]
    d[]
  }
  
  DTk <- assign_bin(copy(DTk))
  DTall <- assign_bin(copy(DT))  # for overall counts across ALL groups
  
  # --- per-group counts (kept groups only) ---
  counts <- DTk[, .(n = .N), by = .(group = get(group_col), bin)]
  
  grid <- CJ(group = sort(keep_groups), bin = seq_len(nbins), unique = TRUE)
  out_groups <- counts[grid, on = .(group, bin)]
  out_groups[is.na(n), n := 0L]
  
  # --- overall counts across ALL groups as group = 0 ---
  out_total <- DTall[, .(n = .N), by = bin]
  grid0 <- data.table(bin = seq_len(nbins))
  out_total <- out_total[grid0, on = "bin"]
  out_total[is.na(n), n := 0L]
  out_total[, group := 0]
  
  out <- rbindlist(list(out_total[, .(group, n, bin)],
                        out_groups[, .(group, n, bin)]),
                   use.names = TRUE)
  
  setorder(out, group, bin)
  out[]
}

## -----------------------------------------------------------------------------
## read in the timeseries data
lf_tc <- list.files(paste0(rawdatapath,"mutantigen/"), pattern = ".timeseries")

length(lf_tc)

## -----------------------------------------------------------------------------
## make train_syn_tc (training data of total cases)
train_syn_tc <- list()
cnt <- 0
Ncoarse <- 10
df_plot_tc <- NULL
for(i in 1:length(lf_tc)){
  
  ## grab the data frame
  tempdf_orig <- fread(paste0(rawdatapath,"mutantigen/",lf_tc[i]))
  
  ## subset columns and remove first 10% (sometimes there is a weird spike)
  tempdf_orig <- subset(tempdf_orig, select = c("date","totalCases"), date > .1*max(date))
  
  ## order tempdf_orig
  tempdf_orig <- tempdf_orig[order(tempdf_orig$date),]
  
  ## add in series id and time
  tempdf_orig$series_id <- paste0("syn_tc_",i)
  tempdf_orig$time <- 1:nrow(tempdf_orig)
  
  ## subset again and rename cols
  tempdf_orig <- subset(tempdf_orig, select = c("series_id","time","totalCases"))
  names(tempdf_orig) <- c("series_id","time","cases")
  
  ## coarsen the data
  for(j in 1:Ncoarse){
    print(j)
    ## sample the number of time steps
    possible_ids <- intersect(1:nrow(tempdf_orig), 52:(10*52))
    if(length(possible_ids) > 1){
      nsteps <- sample(possible_ids, 1)
    }else{
      nsteps <- nrow(tempdf_orig)
    }
    
    ## make new, coarsened tempdf
    tempdf <- data.table(series_id = paste0(tempdf_orig$series_id[1],"_",j),
                         time = 1:nsteps)
    tempdf$cases <- as.numeric(approx(x = tempdf_orig$time,
                                      y = tempdf_orig$cases,
                                      xout = seq(min(tempdf_orig$time), max(tempdf_orig$time), length.out = nsteps))$y)
    
    ## keep going if nrow(tempdf2) > 20
    if(nrow(tempdf) > 20){
      
      ## ---------------------------------------------------------------------
      ## add in cases1: add in outliers
      tempdf$cases1 <- tempdf$cases
      set.seed(i*j+1)
      if(runif(1) < .25){
        n_outliers <- round(runif(1, 5, 10))
        idx_outliers <- sort(sample(1:nrow(tempdf), n_outliers, replace = F))
        mult_outliers <- runif(n_outliers, runif(1,2,3), runif(1,5,10))
        zero_outlier_idx <- sort(sample(1:length(mult_outliers), round(length(mult_outliers)/2,0), replace = F))
        mult_outliers[zero_outlier_idx] <- runif(1,0,.05)
        for(k in 1:n_outliers){
          tempdf$cases1[idx_outliers[k]] <- tempdf$cases[idx_outliers[k]]*mult_outliers[k]
        }
      }
      
      ## make sure all cases are nonneg and between 0 and 1
      tempdf$cases1 <- pmax(0, tempdf$cases1)
      if(max(tempdf$cases1, na.rm=T) > 0){
        tempdf$cases1 <- tempdf$cases1/max(tempdf$cases1, na.rm=T)
      }
      
      
      ## ---------------------------------------------------------------------
      ## add in cases2: add in noise and outliers
      tempdf$cases2 <- tempdf$cases
      
      ## flip and coin to add in random noise
      set.seed(i*j+19)
      noise_kappa <- runif(1,1.5,3.5)
      tempdf$cases2 <- tempdf$cases2 * runif(nrow(tempdf), 1/noise_kappa, noise_kappa)
      
      
      ## flip a coin to add in outliers
      set.seed(2*i*j+19)
      if(runif(1) < .25){
        n_outliers <- round(runif(1, 5, 10))
        idx_outliers <- sort(sample(1:nrow(tempdf), n_outliers, replace = F)) 
        mult_outliers <- runif(n_outliers, runif(1,2,3), runif(1,5,10))
        zero_outlier_idx <- sort(sample(1:length(mult_outliers), round(length(mult_outliers)/2,0), replace = F))
        mult_outliers[zero_outlier_idx] <- runif(1,0,.05)
        for(k in 1:n_outliers){
          tempdf$cases2[idx_outliers[k]] <- tempdf$cases2[idx_outliers[k]]*mult_outliers[k]
        }
      }
      
      ## make sure all cases are nonneg and between 0 and 1
      tempdf$cases2 <- pmax(0, tempdf$cases2)
      if(max(tempdf$cases2, na.rm=T) > 0){
        tempdf$cases2 <- tempdf$cases2/max(tempdf$cases2, na.rm=T)
      }
      
      ## ------------------------------------------------------------------------
      ## Append to list
      names(tempdf)[which(names(tempdf) == "time_step")] <- "time"
      
      ## cases1
      cnt <- cnt + 1
      tempdf1 <- subset(tempdf, select = c("series_id","time","cases1"))
      names(tempdf1) <- c("series_id","time","cases")
      train_syn_tc[[cnt]] <- tempdf1
      
      ## cases2
      cnt <- cnt + 1
      tempdf2 <- subset(tempdf, select = c("series_id","time","cases2"))
      names(tempdf2) <- c("series_id","time","cases")
      train_syn_tc[[cnt]] <- tempdf2
      
    }
  }
  
  ## save data for plotting
  if(i == 1){
    df_plot_tc <- rbindlist(train_syn_tc, idcol = "id")
    df_plot_tc$type <- "realization"
    df_plot_tc$group <- paste0("Realization ",1+floor((df_plot_tc$id-.01)/2))
    df_plot_tc$noise <- "Scaling +\nOutliers"
    df_plot_tc[df_plot_tc$id %in% seq(2,20,2),]$noise <- "Scaling +\nOutliers +\nNoise"
    
    df_plot_tc_orig <- data.table(id = 0, tempdf_orig)
    df_plot_tc_orig$cases <- df_plot_tc_orig$cases/max(df_plot_tc_orig$cases)
    df_plot_tc_orig$type <- "original"
    df_plot_tc_orig2 <- NULL
    for(ll in 1:length(unique(df_plot_tc$id))){
      if(unique(df_plot_tc$id)[ll] %in% seq(2,20,2)){
        df_plot_tc_orig$group <- df_plot_tc[df_plot_tc$id == unique(df_plot_tc$id)[ll],]$group[1]
        df_plot_tc_orig$noise <- "Scaling +\nOutliers +\nNoise"
      }else{
        df_plot_tc_orig$group <- df_plot_tc[df_plot_tc$id == unique(df_plot_tc$id)[ll],]$group[1]
        df_plot_tc_orig$noise <- "Scaling +\nOutliers"
      }
      df_plot_tc_orig2 <- rbind(df_plot_tc_orig2, df_plot_tc_orig)
    }
    
    df_plot_tc <- rbind(df_plot_tc, df_plot_tc_orig2)
    df_plot_tc$group <- factor(as.factor(df_plot_tc$group), levels = paste0("Realization ",1:10))
    
    ## plot it
      setEPS()     
      postscript(paste0(figpath,"example_synthetic_observation_model.eps"),width = 10, height = 4)
      print(
        ggplot(data = subset(df_plot_tc, type == "realization" & group %in% paste0("Realization ",1:5)))+
          geom_line(aes(time, cases))+
          facet_grid(noise~group, scales="free_x")+
          theme_bw() +
          ylab("Normalized Cases")+
          xlab("")+
          guides(color = "none")+
          theme(
            strip.text = element_text(size = 10)
          )
      )
    dev.off()
    
  }
}

## save train_syn_tc
saveRDS(object = train_syn_tc,
        file = paste0(cleandatapath,"train_syn_tc.RDS"))

## remove the object
rm(train_syn_tc)



## -----------------------------------------------------------------------------
## read in the mutationSeries data data
lf_vac <- list.files(paste0(rawdatapath,"mutantigen/"), pattern = "antitypes")

## the number of MutAntiGen runs
length(lf_vac)

## plot some of the runs for the paper
set.seed(90256)
vac_ids <- sort(sample(1:length(lf_vac),4,replace = F))
ex_vac_df <- NULL
for(i in 1:length(vac_ids)){
  tempdf <- fread(paste0(rawdatapath,"mutantigen/",lf_vac[vac_ids[i]]))
  tempdf_tc <- fread(paste0(rawdatapath,"mutantigen/",lf_tc[vac_ids[i]]))
  
  
  res <- bin_counts(DT = subset(tempdf,year > .1*max(tempdf$year)), 
                    nbins = round(.9*nrow(tempdf_tc),0),
                    time_col = "year",
                    group_col = "antigen",
                    alpha = 0)
  res$std_value <- res$n/max(res$n)
  res$id <- paste0("RunID = ",vac_ids[i])
  res$type <- "VAC"
  res[res$group == 0,]$type <- "TC"
  ex_vac_df <- rbind(ex_vac_df, res)
}
ex_vac_df$id <- factor(as.factor(ex_vac_df$id), levels = paste0("RunID = ",vac_ids))

## plot it
setEPS()     
postscript(paste0(figpath,"example_synthetic_data.eps"),width = 10, height = 4)
print(
  ggplot(data = ex_vac_df)+
    geom_line(aes(x = bin, y=std_value), data = subset(ex_vac_df, type == "TC"))+
    geom_line(aes(x = bin, y=std_value, group = group, color = as.factor(group)), data = subset(ex_vac_df, type == "VAC"), show.legend=F)+
    facet_grid(type ~ id, scales="free")+
    ylab("Normalized Cases")+
    xlab("")+
    ggtitle("Synthetic Runs")+
    theme_bw()+
    theme(plot.title = element_text(hjust = 0.5))
)
dev.off()



## -----------------------------------------------------------------------------
## make train_syn_vac (training data of variant attributable cases)
train_syn_vac <- list()
cnt <- 0
for(i in 1:length(lf_vac)){
  print(i)
  
  ## grab the data frame
  tempdf_orig <- fread(paste0(rawdatapath,"mutantigen/",lf_vac[i]))
  tempdf_orig <- data.frame(tempdf_orig)
  
  ## reformat the and scale
  set.seed(i*64)
  tempdf_reformat <- bin_counts(DT = subset(tempdf_orig,year > .1*max(tempdf_orig$year)), # throw out the first 10% as "burn-in"
                                nbins = sample(52:520,1),  ## pick the time steps 
                                time_col = "year",
                                group_col = "antigen",
                                alpha = 0)
  
  ## pick up to 10 variants
  prob_df <- tempdf_reformat[,list(prob = sum(n)),by="group"]
  pick <- sort(sample(prob_df$group, size = min(c(nrow(prob_df), 10)), prob = sqrt(prob_df$prob), replace = F))
  
  ## continue if there are at least 2 variants with cases
  if(length(pick) > 1){
    
    ## coarsen the data
    for(j in 1:length(pick)){
      
      ## get the column id
      groupid <- pick[j]
      
      ## subset columns and remove first 10% (sometimes there is a weird spike)
      tempdf_small <- subset(tempdf_reformat, group == groupid)
      tempdf_small <- tempdf_small[order(tempdf_small$bin),]
      tempdf <- data.table(series_id = paste0("syn_vac_",i,"_",j),
                            time = tempdf_small$bin,
                            cases = tempdf_small$n)
      
      ## keep going if nrow(tempdf2) > 20
      if(nrow(tempdf) > 20){
        
        ## ---------------------------------------------------------------------
        ## add in cases1: add in outliers
        tempdf$cases1 <- tempdf$cases
        set.seed(i*j+2)
        if(runif(1) < .25){
          n_outliers <- round(runif(1, 5, 10))
          idx_outliers <- sort(sample(1:nrow(tempdf), n_outliers, replace = F))
          mult_outliers <- runif(n_outliers, runif(1,2,3), runif(1,5,10))
          zero_outlier_idx <- sort(sample(1:length(mult_outliers), round(length(mult_outliers)/2,0), replace = F))
          mult_outliers[zero_outlier_idx] <- runif(1,0,.05)
          for(k in 1:n_outliers){
            tempdf$cases1[idx_outliers[k]] <- tempdf$cases[idx_outliers[k]]*mult_outliers[k]
          }
        }
        
        ## make sure all cases are nonneg and between 0 and 1
        tempdf$cases1 <- pmax(0, tempdf$cases1)
        if(max(tempdf$cases1, na.rm=T) > 0){
          tempdf$cases1 <- tempdf$cases1/max(tempdf$cases1, na.rm=T)
        }
        
        
        ## ---------------------------------------------------------------------
        ## add in cases2: add in noise and outliers
        tempdf$cases2 <- tempdf$cases
        
        ## flip and coin to add in random noise
        set.seed(i*j+21)
        noise_kappa <- runif(1,1.5,3.5)
        if(sum(tempdf$cases2) == 0){
          tempdf$cases2 <- rpois(nrow(tempdf), runif(nrow(tempdf),0,1))
        }else{
          tempdf$cases2 <- tempdf$cases2 * runif(nrow(tempdf), 1/noise_kappa, noise_kappa)
        }
        
        ## flip a coin to add in outliers
        set.seed(2*i*j+21)
        if(runif(1) < .25){
          n_outliers <- round(runif(1, 5, 10))
          idx_outliers <- sort(sample(1:nrow(tempdf), n_outliers, replace = F))
          mult_outliers <- runif(n_outliers, runif(1,2,3), runif(1,5,10))
          zero_outlier_idx <- sort(sample(1:length(mult_outliers), round(length(mult_outliers)/2,0), replace = F))
          mult_outliers[zero_outlier_idx] <- runif(1,0,.05)
          for(k in 1:n_outliers){
            tempdf$cases2[idx_outliers[k]] <- tempdf$cases2[idx_outliers[k]]*mult_outliers[k]
          }
        }
        
        ## make sure all cases are nonneg and between 0 and 1
        tempdf$cases2 <- pmax(0, tempdf$cases2)
        if(max(tempdf$cases2, na.rm=T) > 0){
          tempdf$cases2 <- tempdf$cases2/max(tempdf$cases2, na.rm=T)
        }
        
        ## ------------------------------------------------------------------------
        ## Append to list
        names(tempdf)[which(names(tempdf) == "time_step")] <- "time"
        
        ## cases1
        cnt <- cnt + 1
        tempdf1 <- subset(tempdf, select = c("series_id","time","cases1"))
        names(tempdf1) <- c("series_id","time","cases")
        train_syn_vac[[cnt]] <- tempdf1
        
        ## cases2
        cnt <- cnt + 1
        tempdf2 <- subset(tempdf, select = c("series_id","time","cases2"))
        names(tempdf2) <- c("series_id","time","cases")
        train_syn_vac[[cnt]] <- tempdf2
      
      }
    }
  }
}

## save train_syn_vac
saveRDS(object = train_syn_vac,
        file = paste0(cleandatapath,"train_syn_vac.RDS"))

## remove the object
rm(train_syn_vac)


## -----------------------------------------------------------------------------
## make train_real (training data of Non-COVID-19, real respiratory data)
df_resp <- readRDS(paste0(rawdatapath,"real_respiratory_data_complete.RDS"))

## get the disease list
disease_list <- NULL
for(i in 1:length(df_resp)){
  
  ## grab the data frame
  tempdf <- df_resp[[i]]
  
  ## append disease list
  disease_list <- c(disease_list, gsub(",","",tolower(tempdf$ts_disease[1])))
}

## get the df_resp indices that are not covid
noncovid_indices <- grep("covid",disease_list,invert=T)

## get the disease list
train_real <- list()
cnt <- 0
real_summary <- NULL
for(i in noncovid_indices){
  
  ## grab the data frame
  tempdf <- data.table(df_resp[[i]])
  
  ## subset to data before January 1st, 2020
  tempdf$ts_date <- ymd(tempdf$ts_dates)
  tempdf <- subset(tempdf, ts_dates < ymd("2020-01-01"))
  
  ## keep going if there are rows
  if(nrow(tempdf)>0){
    ## update the count
    cnt <- cnt + 1
    
    ## update the summary
    real_summary <- rbind(real_summary,
                          data.table(train_real_id = cnt,
                                     mx_ts = max(tempdf$cases),
                                     geography = tempdf$ts_geography[1],
                                     disease = tolower(tempdf$ts_disease[1]),
                                     first_date = min(tempdf$ts_dates),
                                     last_date = max(tempdf$ts_dates),
                                     nobs = nrow(tempdf)))
    ## rename time_step
    names(tempdf)[which(names(tempdf) == "time_step")] <- "time"
    
    ## reorder cols
    tempdf <- subset(tempdf, select = c("series_id","time","cases",setdiff(names(tempdf), c("series_id","time","cases"))))
    
    ## make sure all cases are nonneg and between 0 and 1
    tempdf$cases <- pmax(0, tempdf$cases)
    if(max(tempdf$cases, na.rm=T) > 0){
      tempdf$cases <- tempdf$cases/max(tempdf$cases, na.rm=T)
    }
    
    ## append
    train_real[[cnt]] <- tempdf
  }
  
}

# --------------------------------------
## some summaries of the real, non-COVID-19 respiratory data

# number of time series
length(train_real)

# summary of time series lengths
summary(real_summary$nobs)

# total number of observations
sum(real_summary$nobs)

# what proportion have a length >= 52?
nrow(subset(real_summary, nobs >= 52))/nrow(real_summary)


# --------------------------------------
# plot a few of the time series
ids <- c(806, 405, 67, 1724)
pr1 <- qplot(ts_date, cases*real_summary[real_summary$train_real_id == ids[1],]$mx_ts, data = train_real[[ids[1]]], geom="line")+
  ggtitle("Influenza\nAustralia")+
  ylab("Cases")+
  xlab("")+
  theme(plot.title = element_text(hjust = 0.5))
pr2 <- qplot(ts_date, cases*real_summary[real_summary$train_real_id == ids[2],]$mx_ts, data = train_real[[ids[2]]], geom="line")+
  ggtitle("Influenza\nSouth Africa")+
  ylab("Cases")+
  xlab("")+
  theme(plot.title = element_text(hjust = 0.5))
pr3 <- qplot(ts_date, cases*real_summary[real_summary$train_real_id == ids[3],]$mx_ts, data = train_real[[ids[3]]], geom="line")+
  ggtitle("Diphtheria\nNorth Dakota, USA")+
  ylab("Cases")+
  xlab("")+
  theme(plot.title = element_text(hjust = 0.5))
pr4 <- qplot(ts_date, cases*real_summary[real_summary$train_real_id == ids[4],]$mx_ts, data = train_real[[ids[4]]], geom="line")+
  ggtitle("Pneumonia\nEssex County, Massachusetts, USA")+
  ylab("Cases")+
  xlab("")+
  theme(plot.title = element_text(hjust = 0.5))

## save the example noncovid data
setEPS()     
postscript(paste0(figpath,"example_real_noncovid_data.eps"),width = 14, height = 4)
print(pr1 | pr2 | pr3 | pr4)
dev.off()

## save train_syn_tc
saveRDS(object = train_real,
        file = paste0(cleandatapath,"train_real.RDS"))

## remove the object
rm(train_real)



## -----------------------------------------------------------------------------
## make test_covid (training data of Non-COVID-19, real respiratory data)
df_covid <- readRDS(paste0(rawdatapath,"dfreal_var_attr.RDS"))

## number of sequences being analyzed between June 1st, 202 through December 31st, 2022.
nseqdf <- rbindlist(lapply(df_covid, function(x){data.frame(series_id = x$series_id[1], n_seq = sum(subset(x, date >= ymd("2020-06-01") & date <= ymd("2022-12-31"))$num_seq))}))
nseqdf <- nseqdf[order(nseqdf$n_seq),]
nseqdf$series_id <- factor(as.factor(nseqdf$series_id), levels = nseqdf$series_id)
sum(nseqdf$n_seq)

# plot it by state
qplot(series_id, n_seq, data = nseqdf) + 
  scale_y_log10()+
  coord_flip()


## get the disease list
test_covid <- list()
covid_summary <- NULL
for(i in 1:length(df_covid)){
  
  ## grab the data frame
  tempdf <- data.table(df_covid[[i]])
  
  ## rename time_step
  names(tempdf)[which(names(tempdf) == "time_step")] <- "time"
  
  ## remove "num_seq"
  tempdf <- subset(tempdf, select = setdiff(names(tempdf),"num_seq"))
  
  ## reorder cols
  tempdf <- subset(tempdf, select = c("series_id","time","cases",setdiff(names(tempdf), c("series_id","time","cases"))))
  
  ## make sure all cases are nonneg and between 0 and 1
  tempdf$cases <- pmax(0, tempdf$cases)
  
  
  ## order by time
  tempdf <- tempdf[order(tempdf$date),]
  
  ## append
  test_covid[[i]] <- tempdf
  
  ## update summary
  covid_summary <- rbind(covid_summary, 
                         data.table(test_covid_id = i,
                                    state = gsub("unitedstates_","",tempdf$series_id[1])))
}

## save train_syn_tc
saveRDS(object = test_covid,
        file = paste0(cleandatapath,"test_covid.RDS"))

## remove the object
rm(test_covid)

