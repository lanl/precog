############################################################
# Adaptive Point-Forecast Ensembles: 5 Methods + Comparison
############################################################

set.seed(123)

##############################
# 0. Simulation Setup
##############################

T_max <- 200          # number of time steps
M     <- 3            # number of base models

# True data-generating process (DGP)
# y_t = 0.5 * y_{t-1} + epsilon_t
y_true <- numeric(T_max)
epsilon_sd <- 1

# Base model forecasts:
# Model 1: correctly specified AR(1) with some noise
# Model 2: biased (too high)
# Model 3: underestimates y (shrunk toward 0)
f_mat <- matrix(NA_real_, nrow = T_max, ncol = M)

y_true[1] <- rnorm(1, 0, epsilon_sd)
for (t in 2:T_max) {
  y_true[t] <- 0.5 * y_true[t - 1] + rnorm(1, 0, epsilon_sd)
}

for (t in 2:T_max) {
  # model 1: good AR(1) with noise
  f_mat[t, 1] <- 0.5 * y_true[t - 1] + rnorm(1, 0, 0.5)
  # model 2: biased upward
  f_mat[t, 2] <- 0.7 * y_true[t - 1] + 1 + rnorm(1, 0, 0.7)
  # model 3: shrunk toward zero
  f_mat[t, 3] <- 0.3 * y_true[t - 1] + rnorm(1, 0, 0.4)
}
# For t=1, just initialize forecasts to y_true[1]
f_mat[1, ] <- y_true[1]

################################
# Helper: MSE function
################################
mse <- function(x, y) mean((x - y)^2, na.rm = TRUE)

##################################################
# 1. Exponentially Weighted Average (Hedge / EWA)
##################################################

make_ewa_ensemble <- function(num_models, eta = 0.05) {
  list(
    M       = num_models,
    eta     = eta,
    weights = rep(1 / num_models, num_models)
  )
}

update_ewa_ensemble <- function(ens, forecasts, y_t) {
  M   <- ens$M
  w   <- ens$weights
  eta <- ens$eta
  
  # Squared-error loss
  loss <- (y_t - forecasts)^2
  
  # Hedge/EWA update: w_m <- w_m * exp(-eta * loss_m), then renormalize
  new_w <- w * exp(-eta * loss)
  new_w <- new_w / sum(new_w)
  
  ens$weights <- new_w
  ens
}

########################################
# 2. Online Linear Regression (RLS)
########################################

make_rls_ensemble <- function(num_models,
                              lambda = 0.99,  # forgetting factor (~1 = slow forgetting)
                              delta  = 1e3) { # initial covariance scale (large = weak prior)
  list(
    M       = num_models,
    lambda  = lambda,
    # Initialize weights at zero (can also use equal weights)
    w       = rep(0, num_models),
    # P is the covariance of the weight estimate (M x M)
    P       = diag(delta, num_models)
  )
}

update_rls_ensemble <- function(ens, forecasts, y_t) {
  # RLS for y_t ~ w^T x_t with forgetting factor lambda
  x <- as.numeric(forecasts)   # predictor vector
  w <- ens$w
  P <- ens$P
  lambda <- ens$lambda
  
  # Gain vector K_t
  Px <- P %*% x
  denom <- lambda + as.numeric(t(x) %*% Px)
  K <- Px / denom
  
  # Prediction error
  y_pred <- sum(w * x)
  e_t    <- y_t - y_pred
  
  # Update weights
  w_new <- w + as.numeric(K * e_t)
  
  # Update covariance
  P_new <- (P - K %*% t(x) %*% P) / lambda
  
  ens$w <- w_new
  ens$P <- P_new
  ens
}

###################################################
# 3. Online Ridge Regression (Regularized RLS-ish)
###################################################
# Simple approach: RLS with a stronger prior on ||w||^2
# by choosing smaller delta and no forgetting (lambda=1),
# plus a prior precision tau controlling initial P.

make_ridge_rls_ensemble <- function(num_models,
                                    tau   = 1.0,   # prior precision on weights
                                    delta = 1.0) { # initial scale for P
  list(
    M      = num_models,
    # Equivalent to lambda = 1 (no forgetting) but with ridge prior
    lambda = 1.0,
    w      = rep(0, num_models),
    # Prior covariance: (1/tau) * I  (smaller = stronger shrinkage)
    P      = (1 / tau) * diag(delta, num_models)
  )
}

update_ridge_rls_ensemble <- function(ens, forecasts, y_t) {
  # Same as RLS but with lambda fixed = 1 (no forgetting)
  x <- as.numeric(forecasts)
  w <- ens$w
  P <- ens$P
  lambda <- ens$lambda  # = 1
  
  Px <- P %*% x
  denom <- lambda + as.numeric(t(x) %*% Px)
  K <- Px / denom
  
  y_pred <- sum(w * x)
  e_t    <- y_t - y_pred
  
  w_new <- w + as.numeric(K * e_t)
  P_new <- (P - K %*% t(x) %*% P) / lambda
  
  ens$w <- w_new
  ens$P <- P_new
  ens
}

#############################################
# 5. Kalman Filter on Time-Varying Weights
#############################################
# State-space model:
#   w_t = w_{t-1} + q_t,      q_t ~ N(0, Q)
#   y_t = x_t^T w_t + r_t,    r_t ~ N(0, R)
#
# Q controls how quickly weights can change; R is observation noise.

make_kalman_ensemble <- function(num_models,
                                 Q_scale = 1e-4,
                                 R       = 1.0) {
  list(
    M       = num_models,
    w       = rep(0, num_models),      # initial weights
    P       = diag(1, num_models),     # initial state covariance
    Q       = Q_scale * diag(1, num_models), # process noise
    R       = R                        # observation noise variance
  )
}

update_kalman_ensemble <- function(ens, forecasts, y_t) {
  x <- as.numeric(forecasts)
  w <- ens$w
  P <- ens$P
  Q <- ens$Q
  R <- ens$R
  
  # ---- Predict step ----
  w_pred <- w           # random walk: w_t = w_{t-1} + noise
  P_pred <- P + Q       # add process noise
  
  # ---- Update step ----
  # Predicted observation
  y_pred <- sum(w_pred * x)
  
  # Innovation
  v_t <- y_t - y_pred
  
  # Innovation covariance
  S_t <- as.numeric(t(x) %*% P_pred %*% x + R)
  
  # Kalman gain
  K_t <- (P_pred %*% x) / S_t
  
  # Updated weights
  w_new <- w_pred + as.numeric(K_t * v_t)
  
  # Updated covariance
  P_new <- P_pred - K_t %*% t(x) %*% P_pred
  
  ens$w <- w_new
  ens$P <- P_new
  ens
}

###########################################
# 6. Rolling-Window Regression (ridge)
###########################################
# At each t, regress y on forecasts from last "window" points:
#   y ~ X w,  w = argmin ||y - Xw||^2 + lambda * ||w||^2

make_rolling_reg_ensemble <- function(num_models,
                                      window = 40,
                                      lambda = 0.1) {
  list(
    M       = num_models,
    window  = window,
    lambda  = lambda,
    # We'll store past X and y and recompute w each step
    X_hist  = NULL,  # matrix of forecasts
    y_hist  = NULL,  # vector of y
    w       = rep(1 / num_models, num_models) # start equal
  )
}

update_rolling_reg_ensemble <- function(ens, forecasts, y_t) {
  x <- matrix(as.numeric(forecasts), nrow = 1)
  # Update history
  if (is.null(ens$X_hist)) {
    ens$X_hist <- x
    ens$y_hist <- y_t
  } else {
    ens$X_hist <- rbind(ens$X_hist, x)
    ens$y_hist <- c(ens$y_hist, y_t)
  }
  
  # Truncate to most recent "window" samples
  n_hist <- nrow(ens$X_hist)
  if (n_hist > ens$window) {
    keep_idx <- (n_hist - ens$window + 1):n_hist
    ens$X_hist <- ens$X_hist[keep_idx, , drop = FALSE]
    ens$y_hist <- ens$y_hist[keep_idx]
  }
  
  # Only fit if we have enough history (say at least M points)
  if (nrow(ens$X_hist) >= ens$M) {
    X <- ens$X_hist
    y <- ens$y_hist
    
    # Ridge solution: (X^T X + lambda I)^(-1) X^T y
    lambda <- ens$lambda
    XtX <- crossprod(X)                 # t(X) %*% X
    XtY <- crossprod(X, y)              # t(X) %*% y
    w_ridge <- solve(XtX + lambda * diag(ens$M), XtY)
    ens$w <- as.numeric(w_ridge)
  }
  
  ens
}

##############################################
# Run all methods on the simulated data
##############################################

# Initialize ensembles
ewa    <- make_ewa_ensemble(M, eta = 0.05)
rls    <- make_rls_ensemble(M, lambda = 0.99, delta = 1e3)
ridge  <- make_ridge_rls_ensemble(M, tau = 1.0, delta = 1.0)
kalman <- make_kalman_ensemble(M, Q_scale = 1e-4, R = 1.0)
roll   <- make_rolling_reg_ensemble(M, window = 40, lambda = 0.1)

# Storage for weights and predictions
ewa_w    <- matrix(NA_real_, nrow = T_max, ncol = M)
rls_w    <- matrix(NA_real_, nrow = T_max, ncol = M)
ridge_w  <- matrix(NA_real_, nrow = T_max, ncol = M)
kalman_w <- matrix(NA_real_, nrow = T_max, ncol = M)
roll_w   <- matrix(NA_real_, nrow = T_max, ncol = M)

ewa_pred    <- rep(NA_real_, T_max)
rls_pred    <- rep(NA_real_, T_max)
ridge_pred  <- rep(NA_real_, T_max)
kalman_pred <- rep(NA_real_, T_max)
roll_pred   <- rep(NA_real_, T_max)

for (t in 1:T_max) {
  x_t <- f_mat[t, ]
  y_t <- y_true[t]
  
  # 1. EWA
  ewa_pred[t] <- sum(ewa$weights * x_t)
  ewa <- update_ewa_ensemble(ewa, x_t, y_t)
  ewa_w[t, ] <- ewa$weights
  
  # 2. RLS
  rls_pred[t] <- sum(rls$w * x_t)
  rls <- update_rls_ensemble(rls, x_t, y_t)
  rls_w[t, ] <- rls$w
  
  # 3. Ridge RLS
  ridge_pred[t] <- sum(ridge$w * x_t)
  ridge <- update_ridge_rls_ensemble(ridge, x_t, y_t)
  ridge_w[t, ] <- ridge$w
  
  # 5. Kalman
  kalman_pred[t] <- sum(kalman$w * x_t)
  kalman <- update_kalman_ensemble(kalman, x_t, y_t)
  kalman_w[t, ] <- kalman$w
  
  # 6. Rolling window regression
  roll_pred[t] <- sum(roll$w * x_t)
  roll <- update_rolling_reg_ensemble(roll, x_t, y_t)
  roll_w[t, ] <- roll$w
}

##########################################
# Compare MSEs (after a burn-in period)
##########################################

burn_in <- 20  # ignore first few points while ensembles stabilize

ewa_mse    <- mse(ewa_pred[burn_in:T_max],    y_true[burn_in:T_max])
rls_mse    <- mse(rls_pred[burn_in:T_max],    y_true[burn_in:T_max])
ridge_mse  <- mse(ridge_pred[burn_in:T_max],  y_true[burn_in:T_max])
kalman_mse <- mse(kalman_pred[burn_in:T_max], y_true[burn_in:T_max])
roll_mse   <- mse(roll_pred[burn_in:T_max],   y_true[burn_in:T_max])

cat("MSE (burn-in =", burn_in, "):\n")
print(c(
  EWA    = ewa_mse,
  RLS    = rls_mse,
  Ridge  = ridge_mse,
  Kalman = kalman_mse,
  Roll   = roll_mse
))

##########################################
# Plot weight evolution for each method
##########################################

par(mfrow = c(3, 2))

matplot(ewa_w, type = "l", lty = 1,
        xlab = "Time", ylab = "Weight", main = "EWA Weights")
legend("topright", legend = paste("Model", 1:M), lty = 1, col = 1:M, bty = "n")

matplot(rls_w, type = "l", lty = 1,
        xlab = "Time", ylab = "Weight", main = "RLS Weights")
legend("topright", legend = paste("Model", 1:M), lty = 1, col = 1:M, bty = "n")

matplot(ridge_w, type = "l", lty = 1,
        xlab = "Time", ylab = "Weight", main = "Ridge-RLS Weights")
legend("topright", legend = paste("Model", 1:M), lty = 1, col = 1:M, bty = "n")

matplot(kalman_w, type = "l", lty = 1,
        xlab = "Time", ylab = "Weight", main = "Kalman Weights")
legend("topright", legend = paste("Model", 1:M), lty = 1, col = 1:M, bty = "n")

matplot(roll_w, type = "l", lty = 1,
        xlab = "Time", ylab = "Weight", main = "Rolling Regression Weights")
legend("topright", legend = paste("Model", 1:M), lty = 1, col = 1:M, bty = "n")

# Final panel: true series and a couple of ensemble predictions
plot(y_true, type = "l", lwd = 2, main = "True vs Ensemble Prediction",
     xlab = "Time", ylab = "Value")
lines(ewa_pred, col = 2)
lines(rls_pred, col = 3)
lines(kalman_pred, col = 4)
legend("topleft",
       legend = c("Truth", "EWA", "RLS", "Kalman"),
       col    = c(1, 2, 3, 4),
       lwd    = c(2, 1, 1, 1),
       bty    = "n")
