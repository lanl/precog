############################################################
## Two-Dataset Stacking Simulation:
## - Small "original" dataset (what we care about)
## - Larger "representative" dataset (for learning weights)
## We:
##   1. Learn stacking weights on representative data
##   2. Fit base learners on original data
##   3. Apply equal-weights vs transferred-stacking to
##      original test set and compare MSE
############################################################

set.seed(123)

## Packages ----
## Install if needed:
## install.packages(c("randomForest", "glmnet"))
library(randomForest)
library(glmnet)
library(parallel)
library(doParallel)

############################################################
## 1. Data-generating processes
############################################################

# Representative DGP: Y = f_rep(X) + noise
simulate_rep_data <- function(n, p = 5, noise_sd = 1) {
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("x", 1:p)
  
  f_rep <- function(x) {
    x1 <- x[1]; x2 <- x[2]; x3 <- x[3]; x4 <- x[4]
    # Nonlinear structure
    2 * x1^2 + sin(pi * x2) + x3 * x4 - 0.5 * x4^2
  }
  
  f_x <- apply(X, 1, f_rep)
  y <- f_x + rnorm(n, sd = noise_sd)
  
  data.frame(y = y, X)
}

# Original DGP: slightly different but "representative"
# (e.g., slightly shifted coefficients / nonlinearity)
simulate_orig_data <- function(n, p = 5, noise_sd = 1.2) {
  X <- matrix(rnorm(n * p), nrow = n, ncol = p)
  colnames(X) <- paste0("x", 1:p)
  
  f_orig <- function(x) {
    x1 <- x[1]; x2 <- x[2]; x3 <- x[3]; x4 <- x[4]
    # Slightly perturbed structure vs representative
    1.8 * x1^2 + 1.1 * sin(pi * x2) + 0.9 * x3 * x4 - 0.4 * x4^2
  }
  
  f_x <- apply(X, 1, f_orig)
  y <- f_x + rnorm(n, sd = noise_sd)
  
  data.frame(y = y, X)
}

############################################################
## 2. Base learners (generic utilities)
############################################################

fit_base_learners <- function(train_dat) {
  # train_dat: data.frame with y and x1...xp
  
  # Linear regression
  mod_lm <- lm(y ~ ., data = train_dat)
  
  # Random forest
  mod_rf <- randomForest(y ~ ., data = train_dat)
  
  # Penalized linear model (LASSO or ridge etc.)
  x_mat <- as.matrix(train_dat[ , grepl("^x", names(train_dat))])
  y_vec <- train_dat$y
  mod_glmnet <- cv.glmnet(x_mat, y_vec, alpha = 1)  # LASSO
  
  list(
    lm     = mod_lm,
    rf     = mod_rf,
    glmnet = mod_glmnet
  )
}

predict_base_learners <- function(models, new_dat) {
  x_mat_new <- as.matrix(new_dat[ , grepl("^x", names(new_dat))])
  
  p_lm     <- predict(models$lm,     newdata = new_dat)
  p_rf     <- predict(models$rf,     newdata = new_dat)
  p_glmnet <- as.numeric(predict(models$glmnet, newx = x_mat_new, s = "lambda.min"))
  
  cbind(
    lm     = p_lm,
    rf     = p_rf,
    glmnet = p_glmnet
  )
}

############################################################
## 3. Stacking utilities
############################################################

# Get K-fold out-of-fold predictions on *some* dataset
# (we'll use this for the representative data)
get_oof_predictions <- function(train_dat, K = 5) {
  n <- nrow(train_dat)
  folds <- sample(rep(1:K, length.out = n))
  
  oof_preds <- matrix(NA_real_, nrow = n, ncol = 3)
  colnames(oof_preds) <- c("lm", "rf", "glmnet")
  
  for (k in 1:K) {
    idx_train <- which(folds != k)
    idx_valid <- which(folds == k)
    
    fold_train <- train_dat[idx_train, , drop = FALSE]
    fold_valid <- train_dat[idx_valid, , drop = FALSE]
    
    models_k <- fit_base_learners(fold_train)
    preds_k  <- predict_base_learners(models_k, fold_valid)
    
    oof_preds[idx_valid, ] <- preds_k
  }
  
  list(
    preds = oof_preds,
    y     = train_dat$y
  )
}

# Meta-learner: ridge regression on OOF preds
fit_meta_learner <- function(oof_preds, y) {
  x_meta <- as.matrix(oof_preds)
  cv.glmnet(x_meta, y, alpha = 0)  # ridge regression
}

# Use meta-learner to combine base predictions
predict_stacked <- function(meta_model, base_preds) {
  x_meta_test <- as.matrix(base_preds)
  as.numeric(predict(meta_model, newx = x_meta_test, s = "lambda.min"))
}

############################################################
## 4. One simulation replicate
##    - Generate representative data
##    - Learn stacking weights on representative data
##    - Generate original (train + test) data
##    - Train base learners on original TRAIN
##    - Evaluate equal-weights vs transferred stacking on
##      original TEST only
############################################################

one_simulation <- function(
    n_rep       = 2000,
    n_orig_train = 100,
    n_orig_test  = 200,
    p           = 5,
    K_rep       = 5
) {
  ## 4.1 Representative data + stacking weights ----
  rep_dat <- simulate_rep_data(n = n_rep, p = p)
  
  # Out-of-fold predictions for stacking on representative data
  oof_rep <- get_oof_predictions(rep_dat, K = K_rep)
  
  # Fit meta-learner on rep data
  meta_rep <- fit_meta_learner(oof_rep$preds, oof_rep$y)
  
  ## 4.2 Original data (train + test) ----
  orig_train <- simulate_orig_data(n = n_orig_train, p = p)
  orig_test  <- simulate_orig_data(n = n_orig_test,  p = p)
  
  # Fit base learners on original training data
  base_models_orig <- fit_base_learners(orig_train)
  
  # Predictions of base models on original test set
  base_preds_orig_test <- predict_base_learners(base_models_orig, orig_test)
  
  # Equal-weights ensemble on original test
  eq_ensemble_pred <- rowMeans(base_preds_orig_test)
  
  # Transferred stacking ensemble on original test:
  # reuse meta_learner learned on representative data
  stacked_pred_orig <- predict_stacked(meta_rep, base_preds_orig_test)
  
  # Compute MSEs on original test set
  y_test <- orig_test$y
  
  mse_lm        <- mean((y_test - base_preds_orig_test[, "lm"])^2)
  mse_rf        <- mean((y_test - base_preds_orig_test[, "rf"])^2)
  mse_glmnet    <- mean((y_test - base_preds_orig_test[, "glmnet"])^2)
  mse_equal     <- mean((y_test - eq_ensemble_pred)^2)
  mse_stacked   <- mean((y_test - stacked_pred_orig)^2)
  
  c(
    lm      = mse_lm,
    rf      = mse_rf,
    glmnet  = mse_glmnet,
    equal   = mse_equal,
    weighted = mse_stacked
  )
}

############################################################
## 5. Run the full simulation
############################################################

parfctn = function(n_sim){
  one_simulation(
    n_rep        = 2000,
    n_orig_train = 100,
    n_orig_test  = 200,
    p            = 5,
    K_rep        = 5
  )
}

sockettype <- "PSOCK"

## Uncomment this to work with a simple example (one run).
# parfctn(3)
ncores = 11
n_sims <- 200
cl <- parallel::makeCluster(spec = ncores,type = sockettype) #, outfile=""
setDefaultCluster(cl)
registerDoParallel(cl)
results <- foreach(i=1:n_sims,
                  .verbose = T,
                  .combine = 'rbind',
                  .packages = c("randomForest",
                                "glmnet"))%dopar%{
                    print(i)
                    parfctn(i)
                  }
stopCluster(cl)

# results <- t(results)  # n_sims x methods

head(results)

############################################################
## 6. Summaries and comparison: original test set only
############################################################

cat("Average MSE over simulations (original test set only):\n")
print(colMeans(results))

cat("\nStandard deviation of MSEs:\n")
print(apply(results, 2, sd))

cat("\nProportion of simulations where transferred stacking beats equal-weights (MSE):\n")
prop_better <- mean(results[, "weighted"] < results[, "equal"])
print(prop_better)

cat("\nPaired t-test (weighted vs equal) on original test-set MSEs:\n")
print(t.test(results[, "weighted"], results[, "equal"], paired = TRUE))

############################################################
## 7. Optional: visualization
############################################################

boxplot(
  results[, c("lm", "glmnet", "equal", "rf", "weighted")],
  main = "MSE on 'real' timeseries across 200 simulations",
  ylab = "MSE", xlab = "Model",
  ylim = c(0,20)
)
