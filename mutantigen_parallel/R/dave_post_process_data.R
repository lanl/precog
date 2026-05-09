## Dave Osthus
## 6-11-25
## Format synthetic MutAntiGen Data
## the input .timeseries files can be found at /projects/epifforma/mutantigen_parallel/outfiles/*.timeseries*
## NOTE: I have not confirmed that this is the correct location 

## set paths
library(this.path)
datapath <- paste0(this.path::here(), "/../data_post_processing/")
syndatapath1 <- paste0(this.path::here(), "/../emma_features/")
syndatapath2 <- paste0(this.path::here(), "/../outfiles/")
codepath <- paste0(this.path::here(), "/")

################################################################################
## load libraries, config, and functions
source(paste0(codepath,"load_libraries.R"))
source(paste0(codepath,"config.R"))
source(paste0(codepath,"functions.R"))


################################################################################
## read in real data
real_df <- fread(paste0(datapath,"Real_Pretty.csv"))

## get the sequence names
real_seq_names <- grep("SEQ",names(real_df),value = T)

## get populations
real_pops <- real_df[grep("unitedstates",series_id),list(population = mean(population)), by = c("series_id")]

## read in, format, and assemble data
lf1 <- list.files(syndatapath1, pattern = ".csv")
lf2 <- list.files(syndatapath2, pattern = ".timeseries")
lf <- c(lf1, lf2)
common <- intersect(gsub("_features.csv","",lf1), gsub(".timeseries","",lf2))
tabulate_df <- data.frame(filename = lf)
tabulate_df$on_both_lists <- NA
tabulate_df$in_final_data <- "no"
tabulate_df$csv <- NA
tabulate_df$nrow <- NA
tabulate_df$prop_seq_NA <- NA
tabulate_df$prop_non_zero_cases <- NA
tabulate_df$common_name <- NA

syn_list <- list()
cnt <- 0
for(i in 1:length(lf)){
  print(paste0(i," of ",length(lf)))
  
  ## get name
  temp_file_name <- lf[i]
  tabulate_df$filename[i] <- temp_file_name
  
  ## classify type
  temp_common <- sub("^([^_]+_[^_.]+).*", "\\1", temp_file_name) %in% common
  tabulate_df$on_both_lists[i] <- temp_common
  tabulate_df$common_name[i] <- sub("^([^_]+_[^_.]+).*", "\\1", temp_file_name)
  
  temp_csv <- ifelse(length(grep(".csv",temp_file_name)) > 0, T, F)
  tabulate_df$csv[i] <- temp_csv
  
  ## if temp_common is on the common list and temp_file_name is a .csv, process it as .csv
  ## if temp_common is on the common list and temp_file_name is a .timeseries, skip it
  ## if temp_common is not on the common list, then process it as either .csv or .timeseries
  temp_process <- ifelse((temp_common & temp_csv) | (!temp_common & !temp_csv), T, F)
  
  ## read in data
  if(temp_process){
    
    if(temp_csv){
      
      ## read in data
      tempdf <- fread(paste0(syndatapath1,lf[i]))
      tabulate_df$nrow[i] <- nrow(tempdf)
      tabulate_df$prop_seq_NA[i] <- mean(is.na(tempdf$SEQ1_nucdiv))
      tabulate_df$prop_non_zero_cases[i] <- mean(tempdf$Cases > 0, na.rm=T)
      
      names(tempdf)[which(names(tempdf) == "Cases")] <- "totalCases"
      names(tempdf)[which(names(tempdf) == "Time")] <- "time_step"
      syn_seq_names <- grep("SEQ",names(tempdf),value=T)
      tempdf <- subset(tempdf, select = c(setdiff(names(tempdf),syn_seq_names), real_seq_names))
      
      parts <- strsplit(lf1[i], "_")[[1]]
      keep <- paste0(paste(parts[1:2], collapse = "_"),".timeseries")
      temp_pop_name <- grep(keep, lf2, value = T)
      temppop <- fread(paste0(syndatapath2,temp_pop_name))
      temppop <- temppop[order(temppop$date),]
      temppop$time_step <- 0:(nrow(temppop)-1)
      
      tempdf <- merge(tempdf, subset(temppop, select = c("time_step","totalN")), by = "time_step", all.x=T)
      names(tempdf)[which(names(tempdf) == "totalN")] <- "population"
      
      if(nrow(tempdf[tempdf$totalCases > 0,])/nrow(tempdf) > .5 & nrow(tempdf) > 104){
        cnt <- cnt + 1
        tabulate_df$in_final_data[i] <- "yes"
        
        # change number of time steps
        set.seed(cnt)
        new_n <- sample(round(nrow(tempdf)/2,0):nrow(tempdf), 1)
        tempdf <- data.table(series_id = paste0("synthetic_",i),
                             population = tempdf$population[1],
                             time_step = 1:new_n,
                             totalCases = as.numeric(interpolate_with_na_check(x = 1:nrow(tempdf),
                                                                               y = tempdf$totalCases,
                                                                               xnew = seq(1,nrow(tempdf),length.out = new_n))),
                             SEQ1_nucdiv = as.numeric(interpolate_with_na_check(x = 1:nrow(tempdf),
                                                                                y = tempdf$SEQ1_nucdiv,
                                                                                xnew = seq(1,nrow(tempdf),length.out = new_n))))
        
        ## add in outliers
        n_outliers <- sample(1:10,1)
        idx_outliers <- sample(1:nrow(tempdf), n_outliers, replace = F)
        mult_outliers <- runif(n_outliers,0,10)
        mult_vec <- rep(1, nrow(tempdf))
        mult_vec[idx_outliers] <- mult_outliers
        tempdf$totalCases <- tempdf$totalCases * mult_vec
        noise_level <- truncnorm::rtruncnorm(n = 1, a = 0, b = 0.5, mean = 0, sd = .2)
        tempdf$totalCases <- round(runif(nrow(tempdf), (1-noise_level)*tempdf$totalCases, (1+noise_level)*tempdf$totalCases),0)
        tempdf$synthetic <- 1
        names(tempdf)[which(names(tempdf) == "totalCases")] <- "cases"
        tempdf <- subset(tempdf, select = c("series_id","time_step","cases","population","synthetic",real_seq_names))
        
        print(mean(is.na(tempdf$SEQ1_nucdiv)))
        
        ## make sure tempdf has the same column order as real_df
        missing_cols <- setdiff(names(real_df), names(tempdf))
        for (col in missing_cols) {
          class_type <- class(real_df[[col]])[1]
          
          if (class_type == "IDate") {
            tempdf[, (col) := as.IDate(NA)]
          } else if (class_type == "factor") {
            tempdf[, (col) := factor(NA, levels = levels(real_df[[col]]))]
          } else {
            tempdf[, (col) := vector(class_type, .N)]
            set(tempdf, i = NULL, j = col, value = NA)
          }
        }
        setcolorder(tempdf, names(real_df))
        
        tempdf$sim_number = as.numeric(sub(".*_(\\d+)", "\\1", tabulate_df$common_name[i]))
        syn_list[[cnt]] <- tempdf
      }
        
    }else{
      tempdf <- fread(paste0(syndatapath2,lf[i]))
      tabulate_df$nrow[i] <- nrow(tempdf)
      tabulate_df$prop_seq_NA[i] <- 1
      tabulate_df$prop_non_zero_cases[i] <- mean(tempdf$totalCases > 0, na.rm=T)
      
      tempdf <- tempdf[order(tempdf$date),]
      tempdf$time_step <- 1:nrow(tempdf)
      names(tempdf)[which(names(tempdf) == "totalN")] <- "population"
      # tempdf <- subset(tempdf, time_step > .1*max(tempdf$time_step))
      # tempdf$time_step <- 1:nrow(tempdf)
      
      if(nrow(tempdf[tempdf$totalCases > 0,])/nrow(tempdf) > .5 & nrow(tempdf) > 104){
        cnt <- cnt + 1
        tabulate_df$in_final_data[i] <- "yes"
        
        # change number of time steps
        set.seed(cnt)
        new_n <- sample(round(nrow(tempdf)/2,0):nrow(tempdf), 1)
        tempdf <- data.table(series_id = paste0("synthetic_",i),
                             population = tempdf$population[1],
                             time_step = 1:new_n,
                             totalCases = as.numeric(interpolate_with_na_check(x = 1:nrow(tempdf),
                                                                               y = tempdf$totalCases,
                                                                               xnew = seq(1,nrow(tempdf),length.out = new_n))),
                             SEQ1_nucdiv = NA)
        
        ## add in outliers
        n_outliers <- sample(1:10,1)
        idx_outliers <- sample(1:nrow(tempdf), n_outliers, replace = F)
        mult_outliers <- runif(n_outliers,0,10)
        mult_vec <- rep(1, nrow(tempdf))
        mult_vec[idx_outliers] <- mult_outliers
        tempdf$totalCases <- tempdf$totalCases * mult_vec
        noise_level <- truncnorm::rtruncnorm(n = 1, a = 0, b = 0.5, mean = 0, sd = .2)
        tempdf$totalCases <- round(runif(nrow(tempdf), (1-noise_level)*tempdf$totalCases, (1+noise_level)*tempdf$totalCases),0)
        tempdf$synthetic <- 1
        names(tempdf)[which(names(tempdf) == "totalCases")] <- "cases"
        tempdf <- subset(tempdf, select = c("series_id","time_step","cases","population","synthetic",real_seq_names))
        
        print(mean(is.na(tempdf$SEQ1_nucdiv)))
        
        ## make sure tempdf has the same column order as real_df
        missing_cols <- setdiff(names(real_df), names(tempdf))
        for (col in missing_cols) {
          class_type <- class(real_df[[col]])[1]
          
          if (class_type == "IDate") {
            tempdf[, (col) := as.IDate(NA)]
          } else if (class_type == "factor") {
            tempdf[, (col) := factor(NA, levels = levels(real_df[[col]]))]
          } else {
            tempdf[, (col) := vector(class_type, .N)]
            set(tempdf, i = NULL, j = col, value = NA)
          }
        }
        setcolorder(tempdf, names(real_df))

        tempdf$sim_number = as.numeric(sub(".*_(\\d+)", "\\1", tabulate_df$common_name[i]))
        syn_list[[cnt]] <- tempdf
      }
    }
  }
}

## turn list into data.table
syn_df <- rbindlist(syn_list)

## change any Inf to NA
for (col in names(syn_df)) {
  set(syn_df, i = which(is.infinite(syn_df[[col]])), j = col, value = NA)
}

## print summary
print(summary(real_pops$population))
print(summary(syn_df$population))
print(head(syn_df))
print(sum(sapply(syn_df, function(col) sum(is.infinite(col)))))

# write out the cleaned synthetic data
fwrite(x = syn_df,
       file = paste0(datapath,"Synthetic_Pretty.csv"),
       row.names=F)

fwrite(x = tabulate_df,
       file = paste0(datapath,"tabulate_df.csv"),
       row.names=F)

################################################################################
## some analysis one the synthetic files

# ## subset to runs that did not make trees
# df2 <- subset(tabulate_df, on_both_lists == F)
# 
# ## subset to runs where at least half of the cases are greater than 0
# df3 <- subset(df2, prop_non_zero_cases > .5)
# 
# ## make sure it's long enough
# df4 <- subset(df3, nrow > 104)
# df4$id <- as.numeric(gsub("out_","",df4$common_name))
# df4 <- df4[order(df4$id),]
# fwrite(x = data.frame(filename = df4$filename), file = "/Users/dosthus/Desktop/good_runs_no_trees.csv")









