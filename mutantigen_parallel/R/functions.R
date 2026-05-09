## Dave Osthus
## 6-4-25
## 4 quadrants paper functions

################################################################################
## list the functions
function_names <- c(
  "make_train_X",
  "prep4fcst",
  "fit_gbm",
  "ts2features",
  "next_weeks",
  "prep_data",
  "safe_prcomp",
  "safe_plsr",
 "weighted_interval_score",
 "shrink_int",
 "generate_balanced_pairs",
 "lagged_diff"
)

###### FUNCTION TEMPLATE ################################################
#' Title of the Function
#'
#' One-line description of what the function does.
#'
#' @param param1 Description of first parameter.
#' @param param2 Description of second parameter (if any).
#' @return Description of the return value.
#' @export
function_name <- function(param1, param2) {
  # function body
}
###### END OF FUNCTION TEMPLATE #########################################


################################################################################
#' prepare all the data for regression
#'
#' takes the data frame, max horizon and max_lag as input
#'
#' @param df a data.frame/data.table with columns: series_id, time, and cases
#' @param max_horizon is the integer max forecast horizon
#' @param max_lag is the integer max lag for features creation
#' @return a list: a single data frame with case response and a fitted pls object
#' @export
make_train_X <- function(df, 
                         epi_cols = NULL,
                         seq_cols = NULL,
                         max_horizon = 4, 
                         max_lag = 12,
                         nrow_synthetic = 50000){
  
  ##############################################################################
  ## prep training data
  ##############################################################################
  
  ## save df, then subset df to time < last_obs_date for real data, keeping all the synthetic
  full_df <- df
  
  
  ##############################################################################
  ## start making training data for real data
  ##############################################################################
  print("Prepare real training data")
  start_time <- Sys.time()
  
  ## subset to real data
  real_df <- subset(full_df, synthetic == 0)
  
  
  ## get unique regions and dates
  unqregions <- sort(unique(real_df[grep("synthetic",real_df$series_id, invert = T),]$series_id))
  unqdates <- sort(unique(real_df$date))
  
  ## prepare for loop
  real_train_list <- list()
  cnt <- 0
  
  ## here we go
  print("start making real data features")
  for(i in 1:length(unqregions)){
    ## print progress
    print(unqregions[i])
    
    ## loop over dates
    for(j in 1:length(unqdates)){
      
      ## subset data
      tempdf <- subset(real_df, series_id == unqregions[i] & date <= unqdates[j])
      temptestdf <- subset(real_df, series_id == unqregions[i] & date > unqdates[j])
      
      ## if time series is long enough, do some more work
      if((nrow(tempdf) > max_lag) & (nrow(temptestdf) > 0)){
        
        ## unique id counter
        cnt <- cnt + 1

        ########################################################################
        ## set up epi feature vector
        epi_feature_vec <- NULL
        
        ## loop over epi_cols
        if(!is.null(epi_cols)){
          
          ## loop over epi_cols
          for(ee in 1:length(epi_cols)){
            ## order the data frame by date
            tempdf <- tempdf[order(tempdf$time_step),]
            
            ## get the features
            tempfeatures <- ts2features(ts = tempdf[[epi_cols[ee]]], max_lag = max_lag, cases = T)
            
            ## extract the feature vector
            epi_feature_vec <- c(epi_feature_vec, tempfeatures$features)
          }
          
          ## turn vector to matrix and rename
          epi_feature_vec <- data.frame(matrix(epi_feature_vec, nrow=1))
          names(epi_feature_vec) <- gsub("X","EPI",names(epi_feature_vec))
          
        }

        ########################################################################
        ## set up seq feature vector
        seq_feature_vec <- NULL
        
        ## loop over epi_cols
        if(!is.null(seq_cols)){
          
          ## loop over seq_cols
          for(ss in 1:length(seq_cols)){
            ## order the data frame by date
            tempdf <- tempdf[order(tempdf$time_step),]
            
            ## get the features
            tempfeatures <- ts2features(ts = tempdf[[seq_cols[ss]]], max_lag = max_lag, cases = F)
            
            ## extract the feature vector
            seq_feature_vec <- c(seq_feature_vec, tempfeatures$features)
          }
          
          ## turn vector to matrix and rename
          seq_feature_vec <- data.frame(matrix(seq_feature_vec, nrow=1))
          names(seq_feature_vec) <- gsub("X","SEQ",names(seq_feature_vec))
          
        }
        
        ## concatenate epi_feature_vec and seq_feature_vec
        if(!is.null(epi_feature_vec) & !is.null(seq_feature_vec)){
          feature_vec <- cbind(epi_feature_vec, seq_feature_vec)
        }
        if(!is.null(epi_feature_vec) & is.null(seq_feature_vec)){
          feature_vec <- epi_feature_vec
        }
        if(is.null(epi_feature_vec) & !is.null(seq_feature_vec)){
          feature_vec <- seq_feature_vec
        }
        

        ## unpack the features
        temprow <- data.table(cbind(data.frame(unq_id = cnt,
                                               unq_id_long = paste0(unqregions[i],"_",unqdates[j]),
                                               series_id = unqregions[i],
                                               last_date = unqdates[j],
                                               last_time_step = max(tempdf$time_step),
                                               synthetic = tempdf$synthetic[1],
                                               max_lag = max_lag,
                                               std_shift = tempfeatures$shift_std,
                                               std_scale = tempfeatures$scale_std,
                                               std_last_obs = rev(tempfeatures$std_ts_segment)[1]),
                                    feature_vec))

        ## make a matrix of max_horizon by ncol(temprow)
        tempmat <- do.call(rbind, replicate(max_horizon, temprow, simplify = FALSE))
        
        ## order the data frame by date
        temptestdf <- temptestdf[order(temptestdf$time_step),]
        
        ## get the next max_horizon dates
        temptestdates <- next_weeks(max(tempdf$date), n_weeks = max_horizon)
        temptestdates <- ymd(temptestdates)
        
        ## subset test df to only the future dates
        temptestdf$date <- ymd(temptestdf$date)
        temptestdf <- data.frame(subset(temptestdf, date %in% temptestdates))
        
        ## continue if nrow(temptestdf) > 0
        if(nrow(temptestdf) > 0){
          
          ## add in a step ahead label
          temptestdatesdf <- data.frame(date = temptestdates,
                                        step_ahead = 1:length(temptestdates))
          
          ## merge step ahead label with temptestdf
          temptestdf <- merge(temptestdf, temptestdatesdf, by = "date")
          
          ## reformat
          if(tempfeatures$scale_std > 0){
            std_y <- (temptestdf$cases - tempfeatures$shift_std)/tempfeatures$scale_std
          }else{
            std_y <- (temptestdf$cases - tempfeatures$shift_std)
          }
          
          ## add in standardized cases
          temptestdf$std_cases <- std_y
          
          ## rename 
          names(temptestdf)[which(names(temptestdf) == "date")] <- "step_ahead_date"
          
          ## expand training
          tempmat <- tempmat[1:nrow(temptestdf), , drop = FALSE]
          
          ## combine
          tempdf <- cbind(tempmat, 
                          subset(temptestdf, select = setdiff(names(temptestdf),"series_id")))
          
          ## change X to EPI
          names(tempdf) <- gsub("X","EPI",names(tempdf))
          
          ## reorder cols
          tempdf <- subset(tempdf, select = c("unq_id","unq_id_long","series_id","last_date","synthetic","population",
                                              "step_ahead","step_ahead_date",
                                              "std_shift","std_scale","max_lag","std_last_obs",
                                              "cases","std_cases",
                                              names(feature_vec)))
          
          ## concatenate
          real_train_list[[cnt]] <- tempdf
          # train_df <- rbind(train_df, tempdf)
        }
      }
    }
  }
  
  ## make real data frame
  if(length(real_train_list)>0){
    real_train_df <- rbindlist(real_train_list)
  }else{
    real_train_df <- NULL
  }
  
  print(head(real_train_df))
  end_time <- Sys.time()
  print(end_time - start_time)
  
  
  
  
  ##############################################################################
  ## start making training data for synthetic data
  ##############################################################################
  print("Prepare synthetic training data")
  start_time <- Sys.time()
  
  ## continue of there is synthetic data
  syn_train_list <- list()
  syn_cnt <- 0
  
  if(nrow(subset(full_df, synthetic == 1))>0){
    
    ## subset to synthetic data
    syn_df <- subset(full_df, synthetic == 1)
    
    ## get unique regions and dates
    unqregions <- sort(unique(syn_df$series_id))
    
    ## n cuts per synthetic time series
    nrow_synthetic <- max(c(nrow_synthetic, max_horizon*length(unqregions)))
    cuts_per_region <- round((nrow_synthetic/length(unqregions))/max_horizon,0)
    
    print("start making synthetic data features")
    ## loop over synthetic time series'
    for(i in 1:length(unqregions)){
      print(unqregions[i])
      
      ## subset data
      tempdf_series_id <- subset(syn_df, series_id == unqregions[i])
      
      ## determine candidate cuts
      sample_idx <- (max_lag + 1):(nrow(tempdf_series_id) - max_horizon - 1)
      
      ## sample cuts
      temp_cuts <- sort(sample(sample_idx, 
                               min(c(cuts_per_region, length(sample_idx))), 
                               replace = F,
                               prob = sqrt(tempdf_series_id$cases[sample_idx] + 1)))
      
      ## loop over cuts
      for(j in temp_cuts){
        
        ## update counter
        syn_cnt <- syn_cnt + 1
        
        ## subset data
        tempdf <- subset(tempdf_series_id, time_step <= j)
        temptestdf <- subset(tempdf_series_id,  (time_step > j) & (time_step <= (j + max_horizon)))
        
        ## if time series is long enough, do some more work
        if((nrow(tempdf) > max_lag) & (nrow(temptestdf) > 0)){
          
          ## unique id counter
          cnt <- cnt + 1
          
          ########################################################################
          ## set up epi feature vector
          epi_feature_vec <- NULL
          
          ## loop over epi_cols
          if(!is.null(epi_cols)){
            
            ## loop over epi_cols
            for(ee in 1:length(epi_cols)){
              ## order the data frame by date
              tempdf <- tempdf[order(tempdf$time_step),]
              
              ## get the features
              tempfeatures <- ts2features(ts = tempdf[[epi_cols[ee]]], max_lag = max_lag, cases = T)
              
              ## extract the feature vector
              epi_feature_vec <- c(epi_feature_vec, tempfeatures$features)
            }
            
            ## turn vector to matrix and rename
            epi_feature_vec <- data.frame(matrix(epi_feature_vec, nrow=1))
            names(epi_feature_vec) <- gsub("X","EPI",names(epi_feature_vec))
            
          }
          
          
          ########################################################################
          ## set up seq feature vector
          seq_feature_vec <- NULL
          
          ## loop over epi_cols
          if(!is.null(seq_cols)){
            
            ## loop over seq_cols
            for(ss in 1:length(seq_cols)){
              ## order the data frame by date
              tempdf <- tempdf[order(tempdf$time_step),]
              
              ## get the features
              tempfeatures <- ts2features(ts = tempdf[[seq_cols[ss]]], max_lag = max_lag, cases = F)
              
              ## extract the feature vector
              seq_feature_vec <- c(seq_feature_vec, tempfeatures$features)
            }
            
            ## turn vector to matrix and rename
            seq_feature_vec <- data.frame(matrix(seq_feature_vec, nrow=1))
            names(seq_feature_vec) <- gsub("X","SEQ",names(seq_feature_vec))
            
          }
          
          ## concatenate epi_feature_vec and seq_feature_vec
          if(!is.null(epi_feature_vec) & !is.null(seq_feature_vec)){
            feature_vec <- cbind(epi_feature_vec, seq_feature_vec)
          }
          if(!is.null(epi_feature_vec) & is.null(seq_feature_vec)){
            feature_vec <- epi_feature_vec
          }
          if(is.null(epi_feature_vec) & !is.null(seq_feature_vec)){
            feature_vec <- seq_feature_vec
          }
          

          ## unpack the features
          temprow <- data.table(cbind(data.frame(unq_id = cnt,
                                                 unq_id_long = paste0(unqregions[i],"_",j),
                                                 series_id = unqregions[i],
                                                 last_date = ymd(NA),
                                                 last_time_step = max(tempdf$time_step),
                                                 synthetic = tempdf$synthetic[1],
                                                 max_lag = max_lag,
                                                 std_shift = tempfeatures$shift_std,
                                                 std_scale = tempfeatures$scale_std,
                                                 std_last_obs = rev(tempfeatures$std_ts_segment)[1]),
                                      feature_vec))
          
          ## make a matrix of max_horizon by ncol(temprow)
          tempmat <- do.call(rbind, replicate(max_horizon, temprow, simplify = FALSE))
          
          ## order the data frame by date
          temptestdf <- temptestdf[order(temptestdf$time_step),]
          
          ## get the next max_horizon dates
          temptestdates <- temptestdf$time_step
          
          ## subset test df to only the future dates
          temptestdf <- data.frame(subset(temptestdf, time_step %in% temptestdates))
          
          ## continue if nrow(temptestdf) > 0
          if(nrow(temptestdf) > 0){
            
            ## add in a step ahead label
            temptestdatesdf <- data.frame(time_step = temptestdates,
                                          step_ahead = 1:length(temptestdates))
            
            ## merge step ahead label with temptestdf
            temptestdf <- merge(temptestdf, temptestdatesdf, by = "time_step")
            
            ## reformat
            if(tempfeatures$scale_std > 0){
              std_y <- (temptestdf$cases - tempfeatures$shift_std)/tempfeatures$scale_std
            }else{
              std_y <- (temptestdf$cases - tempfeatures$shift_std)
            }
            
            ## add in standardized cases
            temptestdf$std_cases <- std_y
            
            ## rename 
            names(temptestdf)[which(names(temptestdf) == "date")] <- "step_ahead_date"
            
            ## expand training
            tempmat <- tempmat[1:nrow(temptestdf), , drop = FALSE]
            
            ## combine
            tempdf <- cbind(tempmat, 
                            subset(temptestdf, select = setdiff(names(temptestdf),"series_id")))
            
            ## reorder cols
            tempdf <- subset(tempdf, select = c("unq_id","unq_id_long","series_id","last_date","synthetic","population",
                                                "step_ahead","step_ahead_date",
                                                "std_shift","std_scale","max_lag","std_last_obs",
                                                "cases","std_cases",
                                                names(feature_vec)))
            
            ## concatenate
            syn_train_list[[syn_cnt]] <- tempdf

          }
        }
      }
    }
  }
  
  ## make synthetic data frame
  if(length(syn_train_list)>0){
    syn_train_df <- rbindlist(syn_train_list)
  }else{
    syn_train_df <- NULL
  }
  
  end_time <- Sys.time()
  print(end_time - start_time)
  
  ## rbind train_df with syn_train_df
  train_df <- rbind(real_train_df, syn_train_df)
  
  ## get outta here
  return(train_df)
  
}






################################################################################
#' prepare all the data for regression
#'
#' takes the data frame, max horizon and max_lag as input
#'
#' @param df a data.frame/data.table with columns: series_id, time, and cases
#' @param max_horizon is the integer max forecast horizon
#' @param max_lag is the integer max lag for features creation
#' @return a list: a single data frame with case response and a fitted pls object
#' @export
prep4fcst <- function(train_df, test_df){
  
  
  ##############################################################################
  ## perform PCA on training data
  ##############################################################################
  print("PCA")
  start_time <- Sys.time()
  
  ## get covariate columns
  covariate_cols <- c(names(train_df)[grep("EPI",names(train_df))],
                      names(train_df)[grep("SEQ",names(train_df))])
  
  ## perform partial least squares on train_df
  train_X <- data.frame(subset(train_df, select = c("step_ahead",covariate_cols)))
  print(dim(train_X))
  
  ## fit PCA
  pca <- safe_prcomp(x = train_X)
  
  ## Variance explained by each PC
  var_explained <- pca$sdev^2 / sum(pca$sdev^2)
  
  ## Cumulative variance
  cum_var <- cumsum(var_explained)
  
  ## Minimum number of PCs to reach at least 95% cumulative variance
  num_components <- which(cum_var >= 0.95)[1]
  
  ## get principal component features
  train_Z <- data.table(pca$x[,1:num_components])
  names(train_Z) <- paste0("PC",1:ncol(train_Z))
  
  ## append to train_df 
  train_df <- cbind(train_df, train_Z)
  
  ## finish the time
  end_time <- Sys.time()
  print(end_time - start_time)
  

  ##############################################################################
  ## make PCA for the test data
  test_X <- data.frame(subset(test_df, select = c("step_ahead",covariate_cols)))
  
  ## transform into PCA
  test_Z <- data.table(predict(pca, test_X)[,1:num_components])
  
  ## append to test_df 
  test_df <- cbind(test_df, test_Z)
  

  ##############################################################################
  ## add a United States label
  train_df$us <- 0
  if(length(grep("unitedstates",train_df$series_id)) > 0){
    train_df[grep("unitedstates",train_df$series_id),]$us <- 1
  }
  test_df$us <- 0
  if(length(grep("unitedstates",test_df$series_id)) > 0){
    test_df[grep("unitedstates",test_df$series_id),]$us <- 1
  }
  
  
  ##############################################################################
  ## perform Partial least-squares (PLS) on the PC inputs
  start_time <- Sys.time()
  print("Start PLS")
  
  ## make inputs for PLS
  train_X_pls <- data.frame(train_df)[,c(which(names(train_df) == "step_ahead"),
                                         grep("EPI", names(train_df)))]
  test_X_pls <- data.frame(test_df)[,c(which(names(test_df) == "step_ahead"),
                                       grep("EPI", names(test_df)))]
  
  ## fit partial least squares
  ncomp_pls <- max(c(10, num_components))
  mod_pls <- safe_plsr(x = train_X_pls, y = train_df$std_cases, ncomp = ncomp_pls)

  ## get pls train 
  Xtrain_pls <- as.matrix(predict(mod_pls, type = "scores"))
  attr(Xtrain_pls, "dimnames") <- NULL
  Xtrain_pls <- data.table(Xtrain_pls)
  colnames(Xtrain_pls) <- paste0("PLS",1:ncol(Xtrain_pls))
  
  ## get pls for test
  Xtest_pls <- as.matrix(predict(mod_pls, newdata = test_X_pls, type = "scores"))
  attr(Xtest_pls, "dimnames") <- NULL
  Xtest_pls <- data.table(Xtest_pls)
  colnames(Xtest_pls) <- paste0("PLS",1:ncol(Xtest_pls))
  
  ## concatenate
  train_df <- cbind(train_df, Xtrain_pls)
  test_df <- cbind(test_df, Xtest_pls)
  
  ## finish the time
  end_time <- Sys.time()
  print(end_time - start_time)
  
  ## get outta here
  return(list(train_df = train_df,
              test_df = test_df,
              pca = pca))
}









################################################################################
#' Train LightGBM
#'
#' takes a training data set and fits a GBM
#'
#' @param df is a training data set
#' @param df is a training data set
#' @return a fitted GBM object
#' @export
fit_gbm <- function(train_df, 
                    test_df,
                    nrounds = 50000,
                    learning_rate = 0.05,
                    num_leaves = 15,
                    early_stopping_rounds = 20,
                    feature_fraction = 1,
                    bagging_fraction = 1,
                    bagging_freq = 1,
                    seed = 42,
                    nfits = 2,
                    alphavec = c(.01,.025,.1,.25,.5,.75,.9,.975,.99),
                    num_threads = 1,
                    n_lhs = 10,
                    model = temp_model){
  
  # Check data.table format
  if (!is.data.table(train_df) || !is.data.table(test_df)) {
    stop("Both train_df and test_df must be data.table objects.")
  }
  if (!"std_cases" %in% colnames(train_df) || !"std_cases" %in% colnames(test_df)) {
    stop("Both datasets must contain a 'std_cases' column.")
  }
  
  ## order test_df
  test_df <- test_df[order(test_df$step_ahead, test_df$last_date, test_df$series_id),]
  
  # Predictor columns
  predictor_cols <- NULL
  predictor_cols <- c(predictor_cols, 
                      grep("UMAP", colnames(train_df), value = TRUE))
  if (length(predictor_cols) == 0) stop("No valid predictor columns found.")
  
  
  ####################################################################
  ## set up for validation for hyper parameter selection
  ####################################################################
  
  scenario_grid_df <- expand.grid(alpha = alphavec,
                                  horizons = sort(unique(test_df$step_ahead)))
  scenario_grid_df_small1 <- subset(scenario_grid_df, alpha %in% c(max(alphavec[alphavec<.5]), min(alphavec[alphavec>.5])) & horizons %in% setdiff(unique(scenario_grid_df$horizons), range(scenario_grid_df$horizons)))
  scenario_grid_df_small2 <- subset(scenario_grid_df, alpha %in% c(0.5, range(scenario_grid_df$alpha)) & horizons %in% range(scenario_grid_df$horizons))
  scenario_grid_df_small <- rbind(scenario_grid_df_small1,
                                  scenario_grid_df_small2)
  
  ## param limits
  param_df <- data.frame(param = c("num_leaves","learning_rate","min_data_in_leaf","lambda_l2"),
                         lower = c(15, .01, 5, 0),
                         upper = c(63, .1, 100, 10))
  
  ## draw an LHS design
  param_grid <- data.frame(maximinLHS(n = n_lhs, k = nrow(param_df)))
  for(cc in 1:ncol(param_grid)){
    param_grid[,cc] <- param_grid[,cc]*(param_df$upper[cc] - param_df$lower[cc]) + param_df$lower[cc]
  }
  names(param_grid) <- param_df$param
  param_grid$num_leaves <- round(param_grid$num_leaves,0)
  param_grid$min_data_in_leaf <- round(param_grid$min_data_in_leaf,0)
  
  
  
  ###############################################################################
  ## split train into train and validate
  if(model %in% setdiff(models, synthetic_only_models)){
    
    ###############################################################################
    ## pick the candidate dates for removal
    us_dates_df <- train_df[grepl("unitedstates",series_id),
                            list(n = length(step_ahead),
                                 n_step_ahead = length(unique(step_ahead))),
                            by = "last_date"]
    
    ## get candidate dates
    candidate_valid_dates <- sort(unique(subset(us_dates_df, n_step_ahead == max(n_step_ahead) & n > (0.5* max(n)))$last_date))
    
    ## remove the first valid date
    candidate_valid_dates <- candidate_valid_dates[-which.min(candidate_valid_dates)]
    
    ## pick the hold out dates
    if(length(candidate_valid_dates) < nfits){
      nfits <- length(candidate_valid_dates)
    }
    
    ## pick the validation dates
    candidate_valid_dates <- sort(candidate_valid_dates)
    set.seed(seed)
    seed <- seed + 1
    valid_dates <- sort(sample(candidate_valid_dates, nfits, replace = F))


    ####################################################################
    ## 3) parameter optimization
    optim_dates <- sort(ymd(setdiff(as.character(candidate_valid_dates), as.character(valid_dates))))
    if(length(optim_dates) == 0){
      set.seed(seed)
      seed <- seed + 1
      optim_dates <- sample(candidate_valid_dates, 1, replace = F)
    }else{
      set.seed(seed)
      seed <- seed + 1
      optim_dates <- sort(sample(optim_dates, min(nfits, length(optim_dates)), replace=F))
    }
  }else{
      unq_ts <- sort(unique(train_df$series_id))
      set.seed(seed)
      seed <- seed + 1
      valid_out_folds <- sort(sample(unq_ts, size = min(c(50, round(.1*length(unq_ts),0))), replace = F))
  }
  
  
  ## loop over scenarios
  hyperparam_start_time <- Sys.time()
  param_search_df <- NULL
  for(gg in 1:nrow(scenario_grid_df_small)){
    print(gg)
    temp_horizon <- scenario_grid_df_small$horizons[gg]
    temp_alpha <- scenario_grid_df_small$alpha[gg]
    
    ## subset to forecast horizon
    if(model %in% setdiff(models, synthetic_only_models)){
      valid_idx <- intersect(which(grepl("unitedstates",train_df$series_id)),
                             which(grepl(paste(optim_dates, collapse = "|"), train_df$last_date)))
      train_idx <- setdiff(1:nrow(train_df), valid_idx)
    }else{
      valid_idx <- which(train_df$series_id %in% valid_out_folds)
      train_idx <- setdiff(1:nrow(train_df), valid_idx)
    }

    
    ## separate train and valid by hold out date
    valid_sub <- train_df[valid_idx,]
    train_sub <- train_df[train_idx,]
    
    ## further subset by horizon
    valid_sub <- subset(valid_sub, step_ahead == temp_horizon)
    train_sub <- subset(train_sub, step_ahead == temp_horizon)
    
    ## prepare training data
    X_train <- as.matrix(train_sub[, ..predictor_cols])
    resid_train <- train_sub[["std_resid"]]
    
    ## prepare validation data
    X_valid <- as.matrix(valid_sub[, ..predictor_cols])
    resid_valid <- valid_sub[["std_resid"]]
    
    ## set up train and validation
    dtrain <- lightgbm::lgb.Dataset(data = Matrix(X_train, sparse = TRUE),
                                    label = resid_train)
    dvalid <- lightgbm::lgb.Dataset(data = Matrix(X_valid, sparse = TRUE),
                                    label = resid_valid)
    
    
    ## loop over grid choices
    for(pp in 1:nrow(param_grid)){
      temp_lr <- param_grid$learning_rate[pp]
      temp_num_leaves <- param_grid$num_leaves[pp]
      temp_mdil <- param_grid$min_data_in_leaf[pp]
      temp_l2 <- param_grid$lambda_l2[pp]
      
      # ---- LightGBM parameters ----
      params <- list(
        objective = "quantile",
        alpha = temp_alpha,
        metric = "quantile",
        learning_rate = temp_lr,
        num_leaves = temp_num_leaves,
        min_data_in_leaf = temp_mdil,
        feature_pre_filter = FALSE,
        lambda_l2 = temp_l2,
        feature_fraction = feature_fraction,
        bagging_fraction = bagging_fraction,
        bagging_freq = bagging_freq,
        verbosity = 1,
        num_threads = num_threads,
        deterministic = TRUE,
        seed = seed+10,
        bagging_seed = seed+20,
        feature_fraction_seed = seed+30,
        force_col_wise = T
      )
      
      ## LightGBM model fit
      mod_gbm <- lightgbm::lgb.train(
        params = params,
        data = dtrain,
        valids = list(validation = dvalid),
        nrounds = nrounds,
        early_stopping_rounds = early_stopping_rounds
      )
      
      param_search_df <- rbind(param_search_df,
                               data.table(horizon = temp_horizon,
                                          alpha = temp_alpha,
                                          learning_rate = temp_lr,
                                          num_leaves = temp_num_leaves,
                                          min_data_in_leaf = temp_mdil,
                                          lambda_l2 = temp_l2,
                                          best_score = mod_gbm$best_score,
                                          best_iter = mod_gbm$best_iter))
      
      # print(tail(param_search_df))
    }
  }
  
  ## order param_search_df
  param_search_df <- param_search_df[order(param_search_df$horizon,
                                           param_search_df$alpha),]
  
  ## get weights
  best_param_search_df <- NULL
  for(hh in sort(unique(param_search_df$horizon))){
    for(aa in sort(unique(param_search_df$alpha))){
      eee <- subset(param_search_df, alpha == aa & horizon == hh)
      if(nrow(eee)>0){
        if(diff(range(eee$best_score)) == 0){
          eee$wt <- 1
        }else{
          eee$wt <- (1 - (eee$best_score - min(eee$best_score))/diff(range(eee$best_score)))^2
        }
        eee$wt <- eee$wt/sum(eee$wt)
        best_lr <- weighted.mean(x = eee$learning_rate, w = eee$wt)
        best_num_leaves <- weighted.mean(x = eee$num_leaves, w = eee$wt)
        best_min_data_in_leaf <- weighted.mean(x = eee$min_data_in_leaf, w = eee$wt)
        best_l2 <- weighted.mean(x = eee$lambda_l2, w = eee$wt)
        best_param_search_df <- rbind(best_param_search_df,
                                      data.table(horizon = hh,
                                                 alpha = aa,
                                                 best_lr = best_lr,
                                                 best_num_leaves = best_num_leaves,
                                                 best_min_data_in_leaf = best_min_data_in_leaf,
                                                 best_l2 = best_l2))
      }
    }
  }

  ## fill in scenario_grid_df for learning_rate
  scenario_grid_df <- unique(scenario_grid_df[, c("alpha", "horizons")])
  
  ## merge with best param search df
  scenario_grid_df2 <- merge(scenario_grid_df,
                             best_param_search_df,
                             by.x=c("horizons","alpha"),
                             by.y=c("horizon","alpha"),
                             all.x=T)
  
  best_param_search_df$horizon_std <- (best_param_search_df$horizon - min(best_param_search_df$horizon))/diff(range(best_param_search_df$horizon))
  scenario_grid_df2$horizon_std <- (scenario_grid_df2$horizons - min(best_param_search_df$horizon))/diff(range(best_param_search_df$horizon))
  
  for(ii in 1:nrow(scenario_grid_df2)){
    if(is.na(scenario_grid_df2$best_lr[ii])){
      temp_idx <- as.numeric(FNN::knnx.index(data = subset(best_param_search_df, select=c("horizon_std","alpha")),
                                             query = scenario_grid_df2[ii,c("horizon_std","alpha")],
                                             k = 3))
      temp_dist <- as.numeric(FNN::knnx.dist(data = subset(best_param_search_df, select=c("horizon_std","alpha")),
                                             query = scenario_grid_df2[ii,c("horizon_std","alpha")],
                                             k = 3))
      temp_wt <- (1/temp_dist)/sum(1/temp_dist)
      temp_lr <- weighted.mean(x = best_param_search_df[temp_idx,]$best_lr, w = temp_wt)
      temp_num_leaves <- weighted.mean(x = best_param_search_df[temp_idx,]$best_num_leaves, w = temp_wt)
      temp_min_data_in_leaf <- weighted.mean(x = best_param_search_df[temp_idx,]$best_min_data_in_leaf, w = temp_wt)
      temp_l2 <- weighted.mean(x = best_param_search_df[temp_idx,]$best_l2, w = temp_wt)
      
      scenario_grid_df2[ii,]$best_lr <- temp_lr
      scenario_grid_df2[ii,]$best_num_leaves <- temp_num_leaves
      scenario_grid_df2[ii,]$best_min_data_in_leaf <- temp_min_data_in_leaf
      scenario_grid_df2[ii,]$best_l2 <- temp_l2
      
    }
  }
  
  scenario_grid_df <- scenario_grid_df2
  scenario_grid_df$best_num_leaves <- round(scenario_grid_df$best_num_leaves,0)
  scenario_grid_df$best_min_data_in_leaf <- round(scenario_grid_df$best_min_data_in_leaf,0)
  
  # ## plot it
  # grid.arrange(
  #   qplot(alpha, horizons, color = best_lr, data = scenario_grid_df, size = I(10))+scale_color_viridis_c(direction = -1),
  #   qplot(alpha, horizons, color = best_num_leaves, data = scenario_grid_df, size = I(10))+scale_color_viridis_c(direction = 1),
  #   qplot(alpha, horizons, color = best_min_data_in_leaf, data = scenario_grid_df, size = I(10))+scale_color_viridis_c(direction = 1),
  #   qplot(alpha, horizons, color = best_l2, data = scenario_grid_df, size = I(10))+scale_color_viridis_c(direction = 1))
  hyperparam_end_time <- Sys.time()
  hyperparam_runtime <- difftime(hyperparam_end_time, hyperparam_start_time, units = "secs")
  

  
  
  ##############################################################################
  ## now do the actual fits, with the CV-selected params
  
  ## set up for storage
  test_df_final <- list()
  
  ## set up for different fits
  lgb_start_time <- Sys.time()
  cnt <- 0
  for(j in 1:nfits){
    
    print(paste0("Fit LightGBM Model ",j," of ",nfits))
    
    ## subset to forecast horizon
    if(model %in% setdiff(models, synthetic_only_models)){
      valid_idx <- intersect(which(grepl("unitedstates",train_df$series_id)),
                             which(grepl(valid_dates[j], train_df$last_date)))
      train_idx <- setdiff(1:nrow(train_df), valid_idx)
    }else{
      set.seed(seed)
      seed <- seed + 1
      pick50_ts <- sort(sample(unique(train_df$series_id), 50, replace=F))
      
      ids <- do.call(rbind, lapply(pick50_ts, function(g) {
        subset_rows <- train_df[train_df$series_id == g, ]
        set.seed(seed)
        seed <- seed + 1
        subset_rows[sample(nrow(subset_rows), 1), ]$unq_id_long
      }))
      
      valid_idx <- sort(which(train_df$unq_id_long %in% ids))
      train_idx <- setdiff(1:nrow(train_df), valid_idx)
    }
    
    ## separate train and valid by hold out date
    main_valid_sub <- train_df[valid_idx,]
    main_train_sub <- train_df[train_idx,]
    
    ## make storage
    test_preds <- NULL
    
    ## different model per forecast horizon
    unq_step_ahead <- sort(unique(test_df$step_ahead))
    for(ss in 1:length(unq_step_ahead)){
      print(paste0("Steps-ahead ",ss," of ",length(unq_step_ahead)))
      
      ## subset to forecast horizon
      train_sub <- subset(main_train_sub, step_ahead == unq_step_ahead[ss])
      valid_sub <- subset(main_valid_sub, step_ahead == unq_step_ahead[ss])
      test_sub  <- subset(test_df, step_ahead == unq_step_ahead[ss])
      test_sub <- test_sub[order(test_sub$step_ahead, test_sub$last_date, test_sub$series_id),]
      
      ## prepare training data
      X_train <- as.matrix(train_sub[, ..predictor_cols])
      # y_train <- train_sub[["std_cases"]]
      # resid_train <- y_train - train_sub[["std_last_obs"]]
      resid_train <- train_sub[["std_resid"]]
      
      ## prepare validation data
      X_valid <- as.matrix(valid_sub[, ..predictor_cols])
      # y_valid <- valid_sub[["std_cases"]]
      # resid_valid <- y_valid - valid_sub[["std_last_obs"]]
      resid_valid <- valid_sub[["std_resid"]]
      
      ## prepare testing data
      X_test  <- as.matrix(test_sub[, ..predictor_cols])
      
      
      ## set up train and validation
      dtrain <- lightgbm::lgb.Dataset(data = Matrix(X_train, sparse = TRUE), 
                                      label = resid_train)
      dvalid <- lightgbm::lgb.Dataset(data = Matrix(X_valid, sparse = TRUE), 
                                      label = resid_valid)
      
      ## prepare for fit
      test_preds_step_ahead <- matrix(NA, ncol = length(alphavec), nrow = nrow(test_sub))
      
      ## loop over quantiles
      for(a in 1:length(alphavec)){
        print(paste0("quantile ",alphavec[a]))
        
        cnt <- cnt + 1
        
        ## unpack the params
        temp_alpha <- alphavec[a]
        temp_lr <- subset(scenario_grid_df, alpha == temp_alpha & horizons == ss)$best_lr[1]    
        temp_nl <- subset(scenario_grid_df, alpha == temp_alpha & horizons == ss)$best_num_leaves[1]    
        temp_mdil <- subset(scenario_grid_df, alpha == temp_alpha & horizons == ss)$best_min_data_in_leaf[1]    
        temp_l2 <- subset(scenario_grid_df, alpha == temp_alpha & horizons == ss)$best_l2[1]    
        
        # ---- LightGBM parameters ----
        params <- list(
          boosting = "gbdt",
          objective = "quantile",
          alpha = temp_alpha,
          metric = "quantile",
          learning_rate = temp_lr,
          num_leaves = temp_nl,
          min_data_in_leaf = temp_mdil,
          feature_pre_filter = FALSE,
          lambda_l2 = temp_l2,
          feature_fraction = feature_fraction,
          bagging_fraction = bagging_fraction,
          bagging_freq = bagging_freq,
          verbosity = -1,
          num_threads = num_threads,
          deterministic = TRUE,
          seed = seed+1,
          bagging_seed = seed+2,
          feature_fraction_seed = seed+3,
          force_col_wise = T
        )
        
        ## LightGBM model fit
        mod_gbm <- lightgbm::lgb.train(
          params = params,
          data = dtrain,
          valids = list(validation = dvalid),
          nrounds = nrounds,
          early_stopping_rounds = early_stopping_rounds
        )
        
        # ---- add predictions to test_sub ----
        test_std_pred  <- test_sub$std_last_obs + as.numeric(predict(mod_gbm, X_test))
        test_preds_step_ahead[,a] <- test_std_pred*test_sub$std_scale + test_sub$std_shift
      }
      
      ## append to test_preds
      test_preds <- rbind(test_preds, test_preds_step_ahead)
    }
    
    ## sort each row
    test_preds <- t(apply(test_preds, 1, sort))
    
    ## ensure no quantile estimates are larger than the population
    test_preds2 <- test_preds
    
    # Find rows where population is NOT NA
    valid_rows <- !is.na(test_df$population)
    
    # Apply row-wise clipping using matrix recycling
    test_preds2[valid_rows, ] <- pmin(test_preds[valid_rows, , drop = FALSE], test_df$population[valid_rows])
    
    ## add to list
    test_df_final[[j]] <- test_preds2
  }
  lgb_end_time <- Sys.time()
  lgb_runtime <- difftime(lgb_end_time, lgb_start_time, units = "secs")
  
  
  ## compute the Tukey's biweight function (essentially take averages but remove outliers)
  test_df_quantile_mean <- apply(simplify2array(test_df_final), c(1, 2), function(x) mean(x, na.rm = TRUE))
  test_df_quantile_median <- apply(simplify2array(test_df_final), c(1, 2), function(x) median(x, na.rm = TRUE))
  test_df_quantile <- 0.5*(test_df_quantile_mean + test_df_quantile_median)
  # test_df_quantile <- apply(simplify2array(test_df_final), c(1, 2), robust_mean_mcd)
  
  ## double check no negatives
  test_df_quantile[test_df_quantile < 0] <- 0
  
  ## double check sorted
  test_df_quantile <- t(apply(test_df_quantile, 1, sort))
  
  
  ## for every step ahead and series id, make sure the quantiles respect an ordering
  test_df_quantile_list <- list()
  for(iii in 1:nrow(test_df)){
    tempdf <- cbind(test_df[iii,][rep(1, ncol(test_df_quantile)), , drop = FALSE],
                    data.frame(alpha = alphavec,
                               model_preds = as.numeric(test_df_quantile[iii,])))
    tempwis <- weighted_interval_score(y = test_df$cases[iii],
                                       alpha = alphavec,
                                       quantiles = as.numeric(test_df_quantile[iii,]))
    tempdf$wis <- tempwis$wis
    tempdf$wis_prop_width <- tempwis$contributions$width
    tempdf$wis_prop_lower_penalty <- tempwis$contributions$lower
    tempdf$wis_prop_upper_penalty <- tempwis$contributions$upper
    tempdf$wis_prop_median <- tempwis$contributions$med_part
    
    tempdf <- cbind(tempdf, tempwis$coverage[rep(1, nrow(tempdf)), , drop = FALSE])
    tempdf <- cbind(tempdf, tempwis$widths[rep(1, nrow(tempdf)), , drop = FALSE])
    
    test_df_quantile_list[[iii]] <- tempdf
  }
  
  
  ## stack the matrices
  test_df_quantile <- rbindlist(test_df_quantile_list)
  test_df_quantile$quantile_type <- "GBM"
  
  ## get outta here
  return(list(train_df = train_df, 
              test_df = test_df_quantile,
              n_features = length(predictor_cols),
              n_instances = nrow(train_df),
              n_gbm_fits = cnt,
              nfits = nfits,
              param_grid = scenario_grid_df,
              hyperparam_runtime = hyperparam_runtime,
              lgb_runtime = lgb_runtime))
}








################################################################################
#' Time series to features
#'
#' Takes a time series as input and returns a row vector of time series features
#'
#' @param ts_list is the output of ts2stdts()
#' @param ts is the original ts
#' @param max_lag the maximum lag to consider
#' @return a list of features, standardization mean, and standardization standard deviation
#' @export
ts2features <- function(ts, max_lag, jitter = T, cases = T) {
  
  ## ensure all non-negative
  if(cases){
    ts <- pmax(0, ts)
  }
  
  ## length of time series
  TT <- length(ts)
  
  ## truncate time series to max_lag. 
  ## if TT < max_lag, pad with -1s
  if(TT >= max_lag){
    ts_small <- rev(rev(ts)[1:max_lag])
  }else{
    ts_small <- rev(rev(c(rep(-1,max_lag), ts))[1:max_lag])
  }
  
  ## make sure everything is non-negative
  if(cases){
    ts_small <- pmax(0, ts_small)
  }
  
  ## compute shift and scale of ts_small
  median_small <- median(ts_small, na.rm=T)
  mean_small <- median(ts_small, na.rm=T)
  shift_small <- 0.5*(median_small + mean_small)
  
  mad_small <- mad(ts_small, na.rm=T)
  sd_small  <- sd(ts_small, na.rm=T)
  scale_small <- 0.5*(mad_small + sd_small)
  
  ## standardize ts_small
  if(is.na(scale_small) | is.na(shift_small)){
    ts_small_std <- ts_small
  }else{
    if(scale_small > 0){
      ts_small_std <- (ts_small - shift_small)/scale_small
    }else{
      ts_small_std <- (ts_small - shift_small)
    }
  }

  ts_small_std <-as.numeric(ts_small_std)
  
  ##############################################################################
  ## make the features vector
  feature_vec <- NULL
  

  ##############################################################################
  ## 2) relative time series (scaled by whole time series)
  if(sum(!is.na(ts)) <= 0){
    ts_relative <- ts_small
  }else{
    if(max(ts, na.rm=T) > 0){
      ts_relative <- (ts_small - min(ts, na.rm=T))/(max(ts, na.rm=T) - min(ts, na.rm=T))
    }else{
      ts_relative <- ts_small
    }
  }
  ts_relative_diff <- diff(ts_relative)
  ts_diff <- diff(ts_small_std)

  ## compute the z-score of each entry of ts_small if that entry was removed
  ts_z <- holdout_z_score(ts_small)
  
  
  ## make feature vec
  feature_vec <- c(ts_small, # the last max_lags 
                   ts_relative, # the last max_lags / max(ts)
                   ts_small_std, # the last max_lags standardized
                   ts_z, # hold out z-score
                   ## first-order diffs of the above
                   diff(ts_small, lag = 1), 
                   diff(ts_relative, lag = 1),
                   diff(ts_small_std, lag = 1),
                   diff(ts_z, lag = 1),
                   mean(is.na(ts_small)))
  
  ## get outta here
  return(list(features = feature_vec,
              std_ts_segment = as.numeric(ts_small_std),
              shift_std = shift_small,
              scale_std = scale_small))
}




################################################################################
#' Get the next dates by week
#'
#' Takes a date and integer as input and returns the next integer weeks from the date
#'
#' @param start_date a yyyy-mm-dd start date
#' @param n_weeks the number of weeks from the start date to be returned
#' @return a vector of dates
#' @export
next_weeks <- function(start_date, n_weeks) {
  
  # Convert the input string to a Date object
  start_date <- as.Date(start_date)
  
  # Generate a vector of weekly dates (1 to n_weeks weeks ahead)
  dates <- start_date + 7 * seq_len(n_weeks)
  
  # Return as character vector in yyyy-mm-dd format
  return(as.character(dates))
}



################################################################################
#' do the proper subsetting of full_train_df to prepare for test
#'
#' Takes full_train_df, a last_obs_date, and a model (M1, M2, M3, M4, M5, M6) as input and returns train_df and test_df for the experiment
#'
#' @param full_train_df full training data
#' @param last_obs_date is the date yyyy_mm_dd of the last observation
#' @param model is either model M1, M2, M3, or M4, where
#'  M1: real data, epi features only
#'  M2: real and synthetic data, epi features only
#'  M3: real data, epi and sequence features
#'  M4: real and synthetic data, epi and sequence features
#' @return a vector of dates
#' @export
prep_data <- function(full_train_df, last_obs_date, model, eval_range) {
  
  ##############################################################################
  ## define training data column types
  epi_cols <- grep("EPI", names(full_train_df), value = T)
  seq_cols <- grep("SEQ", names(full_train_df), value = T)
  id_cols <- setdiff(names(full_train_df),c(epi_cols, seq_cols))
  
  ##############################################################################
  ## define columns for each model
  if(model == "M1" | model == "M2" | model == "M3"){
    mycols <- c(id_cols, epi_cols)
  }
  if(model == "M4" | model == "M5" | model == "M6"){
    mycols <- c(id_cols, epi_cols, seq_cols)
  }

  ##############################################################################
  ## only real data
  if(model %in% c("M1","M4")){
    ## training data
    train_df <- subset(full_train_df, 
                       select = mycols, 
                       synthetic == 0 & 
                         last_date <= last_obs_date & 
                         step_ahead_date <= last_obs_date)
    
    ## test data
    test_df <- subset(full_train_df, 
                      select = mycols,
                      synthetic == 0 & 
                        last_date == last_obs_date & 
                        grepl("unitedstates", series_id))
  }
  
  ##############################################################################
  ## real and synthetic data
  if(model %in% c("M3","M6")){
    ## training data
    train_df <- subset(full_train_df, 
                       select = mycols, 
                       (synthetic == 1) | (synthetic == 0 & last_date <= last_obs_date & step_ahead_date <= last_obs_date))
    
    ## test data
    test_df <- subset(full_train_df, 
                      select = mycols,
                      synthetic == 0 & 
                        last_date == last_obs_date & 
                        grepl("unitedstates", series_id))
  }

  ##############################################################################
  ## only synthetic data
  if(model %in% synthetic_only_models){
    ## training data
    train_df <- subset(full_train_df, 
                       select = mycols, 
                       synthetic == 1)
    
    ## test data
    test_df <- subset(full_train_df, 
                      select = mycols,
                      synthetic == 0 & 
                      last_date >= min(eval_range) &
                      last_date <= max(eval_range) &
                      grepl("unitedstates", series_id))
    }
  

  ## get outta here
  return(list(train_df = train_df,
              test_df = test_df))
}





################################################################################
#' do principal components safely
#'
#' do principal components 
#'
#' @param full_train_df full training data
#' @param last_obs_date is the date yyyy_mm_dd of the last observation
#' @param model is either model M1, M2, M3, or M4, where
#'  M1: real data, epi features only
#'  M2: real and synthetic data, epi features only
#'  M3: real data, epi and sequence features
#'  M4: real and synthetic data, epi and sequence features
#' @return a vector of dates
#' @export
safe_prcomp <- function(x, center = TRUE) {
  tryCatch(
    mod <- prcomp(x, center = center, scale. = TRUE),
    error = function(e) {
      message("prcomp() with scale. = TRUE failed, switching to scale. = FALSE")
      mod <- prcomp(x, center = center, scale. = FALSE)
    },
    warning = function(w) {
      # Handle warnings silently or print message if you prefer
      invokeRestart("muffleWarning")  # suppress warning
    }
  )
  
  ## get outta here
  return(pca = mod)
  
}



################################################################################
#' do PLS safely
#'
#' do PLS components 
#'
#' @param full_train_df full training data
#' @param last_obs_date is the date yyyy_mm_dd of the last observation
#' @param model is either model M1, M2, M3, or M4, where
#'  M1: real data, epi features only
#'  M2: real and synthetic data, epi features only
#'  M3: real data, epi and sequence features
#'  M4: real and synthetic data, epi and sequence features
#' @return a vector of dates
#' @export
safe_plsr <- function(x, y, ncomp = 2, ...) {
  tryCatch(
    mod <- plsr(y ~ . -1, data = x, scale=T, validation = "CV", ncomp = ncomp),
    error = function(e) {
      message("plsr() with scale = TRUE failed, switching to scale = FALSE")
      mod <- plsr(y ~ . -1, data = x, scale=F, validation = "CV", ncomp = ncomp)
    },
    warning = function(w) {
      invokeRestart("muffleWarning")  # suppress warnings
    }
  )
  
  ## get outta here
  return(pls = mod)
  
}





################################################################################
#' Compute the Weighted Interval Score (WIS) with Coverage and Width Details
#'
#' This function evaluates a set of quantile-based forecasts using the Weighted Interval Score (WIS).
#' It also returns the interval-wise coverage indicators and interval widths for all symmetric alpha pairs.
#'
#' @param y Numeric scalar. The observed true value.
#' @param alpha Numeric vector of quantile levels (e.g., c(0.01, 0.5, 0.99)). Should include symmetric pairs.
#' @param quantiles Numeric vector of predicted quantiles corresponding to each alpha level.
#'
#' @return A named list with:
weighted_interval_score <- function(y, alpha, quantiles) {
  if (length(alpha) != length(quantiles)) {
    stop("alpha and quantiles must be the same length.")
  }
  if (any(alpha <= 0 | alpha >= 1)) {
    stop("All alpha values must be strictly between 0 and 1.")
  }
  
  # Order by alpha
  ord <- order(alpha)
  alpha <- alpha[ord]
  quantiles <- quantiles[ord]
  
  # Identify lower and upper quantiles for symmetric intervals
  lower_idx <- which(alpha < 0.5)
  upper_idx <- rev(which(alpha > 0.5))
  
  if (length(lower_idx) != length(upper_idx)) {
    stop("alpha must contain symmetric central interval pairs (e.g., 0.1 and 0.9).")
  }
  
  # Initialize tracking variables
  width_total <- 0
  lower_total <- 0
  upper_total <- 0
  n_intervals <- length(lower_idx)
  
  coverage_vals <- numeric(n_intervals)
  width_vals <- numeric(n_intervals)
  coverage_names <- character(n_intervals)
  width_names <- character(n_intervals)
  
  ## compute wis
  wis <- 0.5*abs(y - quantiles[which(alpha == .5)])
  med_part <- 0.5*abs(y - quantiles[which(alpha == .5)])
  width_vec <- 0
  lower_penalty_vec <- 0
  upper_penalty_vec <- 0
  for (i in seq_len(n_intervals)) {
    a_lower <- alpha[lower_idx[i]]
    a_upper <- alpha[upper_idx[i]]
    lower <- quantiles[lower_idx[i]]
    upper <- quantiles[upper_idx[i]]
    
    alpha_level <- (1-a_upper) + a_lower
    wt <- alpha_level/2
    coverage_level <- 1 - alpha_level
    label <- round(100 * coverage_level)
    
    width <- upper - lower
    penalty_lower <- ifelse(y < lower, (2 / alpha_level) * (lower - y), 0)
    penalty_upper <- ifelse(y > upper, (2 / alpha_level) * (y - upper), 0)
    int_score <- width + penalty_lower + penalty_upper
    
    ## update
    wis <- wis + wt*int_score
    width_vec <- width_vec + wt*width
    lower_penalty_vec <- lower_penalty_vec + wt*penalty_lower
    upper_penalty_vec <- upper_penalty_vec + wt*penalty_upper
    
    coverage_vals[i] <- as.integer(y >= lower & y <= upper)
    width_vals[i] <- width
    coverage_names[i] <- paste0("cov", label)
    width_names[i] <- paste0("width", label)
  }
  
  ## scale wis
  wis <- (1/(n_intervals + .5))*wis
  width_vec <-(1/(n_intervals + .5))*width_vec
  lower_penalty_vec <- (1/(n_intervals + .5))*lower_penalty_vec
  upper_penalty_vec <- (1/(n_intervals + .5))*upper_penalty_vec
  med_part <- (1/(n_intervals + .5))*med_part
  
  # Contributions as data.frame
  contributions <- data.frame(
    width = width_vec / wis,
    lower = lower_penalty_vec / wis,
    upper = upper_penalty_vec / wis,
    med_part = med_part / wis)

  
  coverage_df <- setNames(as.data.frame(t(coverage_vals)), coverage_names)
  widths_df <- setNames(as.data.frame(t(width_vals)), width_names)
  
  return(list(
    wis = wis,
    contributions = contributions,
    coverage = coverage_df,
    widths = widths_df
  ))
}






############################################################################
#' Compute the Weighted Interval Score (WIS) with Coverage and Width Details
#'
#' This function finds an adjustment to the GBM quantile regression that tries to minimize WIS which keeping good coverage
#' @param mat matrix if quantiles
#' @param med the median of the quantiles
#' @param y the observations from the validation set
#' @param alphavec the vectors of alpha levels
#'
#' @return the coverage weighted WIS
shrink_int <- function(x, mat, med, y, alphavec){
  
  alpha1 <- 10*(exp(x[1])/(1+exp(x[1])))
  alpha2 <- 10*(exp(x[2])/(1+exp(x[2])))
  
  ## center on the median
  lower_cols <- which(1:ncol(mat) < which(alphavec == .5))
  upper_cols <- which(1:ncol(mat) > which(alphavec == .5))
  
  ## make new mat
  matnew <- cbind(alpha1*(mat[,lower_cols] - med) + med,
                  med,
                  alpha2*(mat[,upper_cols] - med) + med)
  
  ## set up storage of WIS and coverage
  wis_vec <- NULL
  cov_df <- NULL
  
  ## loop over instances
  for(iii in 1:nrow(matnew)){
    
    ## compute WIS
    tempwis <- weighted_interval_score(y = y[iii],
                                       alpha = alphavec,
                                       quantiles = as.numeric(matnew[iii,]))
    
    ## unpack WIS
    wis_vec <- c(wis_vec, tempwis$wis)
    
    ## unpack coverage
    cov_df <- rbind(cov_df, tempwis$coverage)
    
  }
  
  ## get wis
  total_wis <- sum(wis_vec)
  
  ## get empirical coverage
  empirical_coverage <- as.numeric(colMeans(cov_df))
  
  ## get nominal coverage
  nominal_coverage <- as.numeric(gsub("cov","",names(cov_df)))/100
  
  ## get average absolute difference between empirical and nominal coverage
  coverage_diff <- mean(abs(empirical_coverage - nominal_coverage))
  
  ## compute objective function
  objective <- total_wis*(1+sqrt(coverage_diff))
  
  ## get outta here
  return(objective)
}




################################################################################
#' Generate N diverse, random pairs from 1:P with balanced roles
#'
#' @param P Integer. Size of the set to sample from (1:P).
#' @param N Integer. Number of pairs to generate.
#'
#' @return A matrix of size N x 2, where each row is a unique pair (i ≠ j),
#'         using up to 2N unique values when possible. Numbers are as evenly distributed
#'         as possible across columns 1 and 2.
#'
#' @examples
#' set.seed(42)
#' generate_balanced_pairs(P = 20, N = 10)
#'
#' @export
generate_balanced_pairs <- function(P, N) {
  # Step 1: Case where 2N unique values are possible
  if (2 * N <= P) {
    elems <- sample(1:P, 2 * N)
    pair_matrix <- matrix(elems, ncol = 2, byrow = TRUE)
  } else {
    # Step 2: Fall back to greedy selection + filler
    all_pairs <- as.data.frame(t(combn(1:P, 2)))
    all_pairs <- all_pairs[sample(nrow(all_pairs)), ]
    
    selected <- list()
    used <- integer(0)
    
    for (i in 1:nrow(all_pairs)) {
      pair <- as.integer(all_pairs[i, ])
      if (length(selected) >= N) break
      
      if (!all(pair %in% used)) {
        selected[[length(selected) + 1]] <- pair
        used <- union(used, pair)
      }
    }
    
    if (length(selected) < N) {
      remaining <- setdiff(1:nrow(all_pairs), 
                           sapply(selected, function(p) which(all_pairs[,1]==p[1] & all_pairs[,2]==p[2])))
      fill_indices <- sample(remaining, N - length(selected))
      selected <- c(selected, lapply(fill_indices, function(i) as.integer(all_pairs[i, ])))
    }
    
    pair_matrix <- do.call(rbind, selected)
  }
  
  # Step 3: Balance appearances across columns
  balance_columns <- function(mat, max_iter = 1000) {
    for (iter in seq_len(max_iter)) {
      col1_tab <- table(mat[,1])
      col2_tab <- table(mat[,2])
      
      all_vals <- union(names(col1_tab), names(col2_tab))
      col1_cts <- as.integer(col1_tab[match(all_vals, names(col1_tab))])
      col2_cts <- as.integer(col2_tab[match(all_vals, names(col2_tab))])
      col1_cts[is.na(col1_cts)] <- 0
      col2_cts[is.na(col2_cts)] <- 0
      
      imbalance <- col1_cts - col2_cts
      
      if (all(abs(imbalance) <= 1)) break  # good enough
      
      # Try flipping one pair to reduce imbalance
      flipped <- FALSE
      for (i in seq_len(nrow(mat))) {
        a <- mat[i, 1]
        b <- mat[i, 2]
        a_idx <- match(a, all_vals)
        b_idx <- match(b, all_vals)
        
        before <- abs(imbalance[a_idx]) + abs(imbalance[b_idx])
        imbalance_tmp <- imbalance
        imbalance_tmp[a_idx] <- imbalance_tmp[a_idx] - 2
        imbalance_tmp[b_idx] <- imbalance_tmp[b_idx] + 2
        after <- abs(imbalance_tmp[a_idx]) + abs(imbalance_tmp[b_idx])
        
        if (after < before) {
          mat[i, ] <- c(b, a)
          flipped <- TRUE
          break
        }
      }
      
      if (!flipped) break  # no improvement possible — exit
    }
    
    mat
  }
  
  
  pair_matrix <- balance_columns(pair_matrix)
  return(pair_matrix)
}

################################################################################
#' compute the lagged difference for a time series
#'
#' @param x is a time series as a vector
#' @param lag is a positive integer for the lag size
#' @return returns a vector of the lags
#' @export
lagged_diff <- function(x, lag) {
  if (length(x) <= lag) return(numeric(0))
  x[(lag + 1):length(x)] - x[1:(length(x) - lag)]
}


################################################################################
#' Compute Mean Differences Between Front and Back Segments of a Vector (NA-Tolerant, NA-Filled with 0)
#'
#' For a numeric vector `x`, this function computes the difference between the mean of the first `k` elements and the mean of the last `k` elements, for all `k` from 1 up to \code{floor(length(x) / 2)}.
#' Missing values (\code{NA}) are ignored in the mean calculations. If a segment consists entirely of \code{NA}s, the resulting difference is set to \code{0}.
#'
#' @param x A numeric vector, potentially containing \code{NA} values.
#'
#' @return A numeric vector of length \code{floor(length(x) / 2)} containing differences \code{mean(x[1:k]) - mean(x[(n - k + 1):n])}, with \code{NA}s replaced by \code{0}.
#'
#' @examples
#' front_back_diff(c(1, 2, 3, 4, 5, NA, 100))
#' front_back_diff(rep(NA, 10))
#'
#' @export
front_back_diff <- function(x) {
  n <- length(x)
  k_max <- floor(n / 2)
  diffs <- sapply(1:k_max, function(k) {
    mean_front <- mean(x[1:k], na.rm = TRUE)
    mean_back <- mean(x[(n - k + 1):n], na.rm = TRUE)
    mean_front - mean_back
  })
  
  # Replace any resulting NA with 0
  diffs[is.na(diffs)] <- 0
  return(diffs)
}


#' Compute Leave-One-Out Z-Scores for a Time Series
#'
#' Calculates a z-score time series where each value is standardized based on
#' the mean and standard deviation of the rest of the series (excluding the current value).
#' Handles edge cases including NA, Inf, constant values, and zero variance.
#'
#' @param x A numeric vector (time series) possibly containing NA, NaN, Inf, or -Inf.
#'
#' @return A numeric vector of z-scores with the same length as `x`. 
#' If `x[i]` is NA, NaN, or not finite, `z[i]` is NA. 
#' If the standard deviation of the hold-out set is zero or undefined, `z[i]` is set to 0.
#'
#' @examples
#' x <- c(1, 2, 3, NA, 2, 1)
#' holdout_z_score(x)
#'
#' x <- c(0, 0, 0)
#' holdout_z_score(x)
#'
#' x <- c(Inf, -Inf, 1, 2)
#' holdout_z_score(x)
#'
#' @export
holdout_z_score <- function(x) {
  n <- length(x)
  z <- rep(NA_real_, n)
  
  if (all(!is.finite(x))) {
    return(z)  # All values are NA, NaN, Inf, or -Inf
  }
  
  for (i in seq_len(n)) {
    if (!is.finite(x[i])) {
      z[i] <- NA  # Skip if x[i] is NA, NaN, Inf, or -Inf
      next
    }
    
    xi <- x[-i]
    xi_valid <- xi[is.finite(xi)]  # Remove NA, NaN, Inf, -Inf
    
    if (length(xi_valid) < 2) {
      z[i] <- NA
      next
    }
    
    mu <- mean(xi_valid)
    sigma <- sd(xi_valid)
    
    if (is.na(sigma) || sigma < .Machine$double.eps) {
      z[i] <- 0
    } else {
      z[i] <- (x[i] - mu) / sigma
    }
  }
  
  return(z)
}




#' Safely interpolate a time series using nearest finite neighbors
#'
#' @param x Numeric vector (may include NA, NaN, Inf)
#' @return A numeric vector z, interpolated from neighbors
#'
#' @export
safe_linear_interp <- function(x) {
  T_len <- length(x)
  z <- rep(NA_real_, T_len)
  
  is_valid <- function(val) is.finite(val)
  
  find_next_valid <- function(x, start_idx, direction) {
    i <- start_idx
    while (i >= 1 && i <= length(x)) {
      if (is_valid(x[i])) return(x[i])
      i <- i + direction
    }
    return(NA_real_)
  }
  
  for (i in seq_len(T_len)) {
    if (i > 1 && i < T_len) {
      left  <- find_next_valid(x, i - 1, -1)
      right <- find_next_valid(x, i + 1, 1)
    } else if (i == 1) {
      left  <- find_next_valid(x, 2, 1)
      right <- find_next_valid(x, 3, 1)
    } else if (i == T_len) {
      left  <- find_next_valid(x, T_len - 1, -1)
      right <- find_next_valid(x, T_len - 2, -1)
    }
    
    if (is_valid(left) && is_valid(right)) {
      z[i] <- 0.5 * (left + right)
    } else if (is_valid(left)) {
      z[i] <- left
    } else if (is_valid(right)) {
      z[i] <- right
    } else {
      z[i] <- x[i]  # No valid values found
    }
  }
  
  return(z)
}


################################################################################
#' Compute a Robust Mean with Non-Negativity Constraint via MCD
#'
#' This function estimates a robust mean using the Minimum Covariance Determinant (MCD),
#' then enforces a non-negativity constraint on the final result.
#' It uses `MASS::cov.rob` under the hood. If the input has too few values or no spread,
#' it falls back to the arithmetic mean.
#'
#' @param x A numeric vector (may contain NAs or negative values).
#'
#' @return A single numeric value representing a non-negative robust estimate of central tendency.
#'
#' @examples
#' robust_mean_mcd(c(100, 20, rnorm(8)))     # robust to large outliers
#' robust_mean_mcd(c(-3, -1, 0.5, 1.2))      # allows small negatives, clamps result
#'
#' @importFrom MASS cov.rob
#' @export
robust_mean_mcd <- function(x) {
  x <- x[!is.na(x)]
  if (length(x) < 2 || IQR(x) == 0) {
    return(max(mean(x), 0))  # fallback + clamp
  }
  est <- as.numeric(MASS::cov.rob(as.matrix(x), method = "mcd")$center)
  return(max(est, 0))  # enforce non-negativity on the output
}



################################################################################
#' NA-aware linear interpolation with safe extrapolation
#'
#' Interpolates at new x-values (`xnew`) using linear interpolation, but only
#' if the surrounding y-values are not `NA`. Extrapolates flat using boundary
#' y-values (`y[min(x)]`, `y[max(x)]`) unless those are `NA`.
#'
#' @param x Numeric vector of x-values (must be increasing).
#' @param y Numeric vector of y-values (same length as x), may contain NAs.
#' @param xnew Numeric vector of new x-values at which to interpolate.
#'
#' @return A numeric vector of interpolated values or NAs.
#'
#' @examples
#' x <- 1:5
#' y <- c(1, NA, NA, 10, 20)
#' xnew <- seq(0, 6, by = 0.5)
#' interpolate_with_na_check(x, y, xnew)
#'
#' @export
interpolate_with_na_check <- function(x, y, xnew) {
  sapply(xnew, function(xi) {
    # Below range: return y[min(x)], even if NA
    if (xi < min(x)) return(y[which.min(x)])
    
    # Above range: return y[max(x)], even if NA
    if (xi > max(x)) return(y[which.max(x)])
    
    # Exactly at an x value
    if (xi %in% x) {
      match_idx <- which(x == xi)[1]
      return(y[match_idx])
    }
    
    # In range: do NA-aware linear interpolation
    left_idx <- max(which(x < xi))
    right_idx <- min(which(x > xi))
    
    yl <- y[left_idx]
    yr <- y[right_idx]
    xl <- x[left_idx]
    xr <- x[right_idx]
    
    if (is.na(yl) || is.na(yr)) {
      return(NA)
    } else {
      slope <- (yr - yl) / (xr - xl)
      return(yl + slope * (xi - xl))
    }
  })
}


