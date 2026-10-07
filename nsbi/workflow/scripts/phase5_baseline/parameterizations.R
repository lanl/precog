#!/usr/bin/env Rscript

#' Parameterization Module for SIR MLE Fitting
#' 
#' Provides modular parameterization schemes for MLE optimization.
#' Each parameterization defines how to:
#' 1. Construct the negative log-likelihood function
#' 2. Compute parameter bounds
#' 3. Compute starting values
#' 4. Transform estimates to standard output format

# ============================================================================
# Parameterization Registry
# ============================================================================

get_parameterization <- function(name) {
  #' Get parameterization configuration by name
  #' 
  #' @param name Character: "betagamma" or "r0rt"
  #' @return List with parameterization functions
  
  if (name == "betagamma") {
    return(parameterization_betagamma())
  } else if (name == "r0rt") {
    return(parameterization_r0rt())
  } else {
    stop(sprintf("Unknown parameterization: %s", name))
  }
}

# ============================================================================
# Beta-Gamma Parameterization (Original)
# ============================================================================

parameterization_betagamma <- function() {
  #' Original parameterization: optimize beta and gamma directly
  #' 
  #' Optimization parameters: (beta, gamma)
  #' Derived parameters: R0 = beta/gamma, recovery_time = 1/gamma
  
  list(
    name = "betagamma",
    param_names = c("beta", "gamma"),
    
    create_nll = function(sir_ode, init_state, times, observed_cases) {
      #' Create negative log-likelihood function
      #' @return Function with signature nll(beta, gamma)
      
      function(beta, gamma) {
        if (beta <= 0 || gamma <= 0) return(1e10)
        
        # Extend times by 1 to get proper diff(R) alignment
        # Day 0 corresponds to recoveries from t=0 to t=1, so we need R(0) and R(1)
        # Day n corresponds to recoveries from t=n to t=n+1
        times_extended <- c(times, max(times) + 1)
        
        out <- tryCatch({
          ode(y = init_state, times = times_extended, func = sir_ode,
              parms = c(beta = beta, gamma = gamma))
        }, error = function(e) NULL)
        
        if (is.null(out)) return(1e10)
        
        # Expected new recoveries per day: diff(R) gives recoveries from t[i] to t[i+1]
        # diff(R)[1] = R(1) - R(0) = expected cases for day 0 (t in [0,1))
        # diff(R)[2] = R(2) - R(1) = expected cases for day 1 (t in [1,2))
        R_t <- out[, "R"]
        expected_cases <- diff(R_t)  # FIXED: No prepended 0!
        expected_cases <- pmax(expected_cases, 1e-6)
        
        nll_val <- -sum(dpois(observed_cases, lambda = expected_cases, log = TRUE))
        
        if (!is.finite(nll_val)) return(1e10)
        return(nll_val)
      }
    },
    
    compute_bounds = function(R0_min, R0_max, rt_min, rt_max, buffer = 0.1) {
      #' Compute parameter bounds with buffer
      
      beta_min <- (R0_min / rt_max) * (1 - buffer)
      beta_max <- (R0_max / rt_min) * (1 + buffer)
      gamma_min <- (1 / rt_max) * (1 - buffer)
      gamma_max <- (1 / rt_min) * (1 + buffer)
      
      list(
        lower = c(beta = beta_min, gamma = gamma_min),
        upper = c(beta = beta_max, gamma = gamma_max)
      )
    },
    
    compute_starting_values = function(strategy, true_R0, true_recovery_time, 
                                      R0_min, R0_max, rt_min, rt_max) {
      #' Compute starting values based on strategy
      
      if (strategy == "prior_mean") {
        R0_start <- exp((log(R0_min) + log(R0_max)) / 2)
        rt_start <- (rt_min + rt_max) / 2
      } else if (strategy == "true_values") {
        R0_start <- true_R0
        rt_start <- true_recovery_time
      } else {
        stop(sprintf("Unknown starting value strategy: %s", strategy))
      }
      
      gamma_start <- 1 / rt_start
      beta_start <- R0_start * gamma_start
      
      list(beta = beta_start, gamma = gamma_start)
    },
    
    transform_estimates = function(estimates) {
      #' Transform optimization results to standard output format
      
      beta <- estimates["beta"]
      gamma <- estimates["gamma"]
      R0 <- beta / gamma
      recovery_time <- 1 / gamma
      
      list(
        beta = as.numeric(beta),
        gamma = as.numeric(gamma),
        R0 = as.numeric(R0),
        recovery_time = as.numeric(recovery_time)
      )
    }
  )
}

# ============================================================================
# R0-Recovery Time Parameterization (New)
# ============================================================================

parameterization_r0rt <- function() {
  #' New parameterization: optimize R0 and recovery_time directly
  #' 
  #' Optimization parameters: (R0, recovery_time)
  #' Derived parameters: gamma = 1/recovery_time, beta = R0 * gamma
  #' 
  #' Advantage: R0 is directly constrained within prior bounds (with buffer)
  
  list(
    name = "r0rt",
    param_names = c("R0", "recovery_time"),
    
    create_nll = function(sir_ode, init_state, times, observed_cases) {
      #' Create negative log-likelihood function
      #' @return Function with signature nll(R0, recovery_time)
      
      function(R0, recovery_time) {
        if (R0 <= 0 || recovery_time <= 0) return(1e10)
        
        # Convert from R0 and recovery_time to beta and gamma
        gamma <- 1 / recovery_time
        beta <- R0 * gamma
        
        if (beta <= 0 || gamma <= 0) return(1e10)
        
        # Extend times by 1 to get proper diff(R) alignment
        # Day 0 corresponds to recoveries from t=0 to t=1, so we need R(0) and R(1)
        # Day n corresponds to recoveries from t=n to t=n+1
        times_extended <- c(times, max(times) + 1)
        
        out <- tryCatch({
          ode(y = init_state, times = times_extended, func = sir_ode,
              parms = c(beta = beta, gamma = gamma))
        }, error = function(e) NULL)
        
        if (is.null(out)) return(1e10)
        
        # Expected new recoveries per day: diff(R) gives recoveries from t[i] to t[i+1]
        # diff(R)[1] = R(1) - R(0) = expected cases for day 0 (t in [0,1))
        # diff(R)[2] = R(2) - R(1) = expected cases for day 1 (t in [1,2))
        R_t <- out[, "R"]
        expected_cases <- diff(R_t)  # FIXED: No prepended 0!
        expected_cases <- pmax(expected_cases, 1e-6)
        
        nll_val <- -sum(dpois(observed_cases, lambda = expected_cases, log = TRUE))
        
        if (!is.finite(nll_val)) return(1e10)
        return(nll_val)
      }
    },
    
    compute_bounds = function(R0_min, R0_max, rt_min, rt_max, buffer = 0.1) {
      #' Compute parameter bounds with buffer
      #' 
      #' Advantage: R0 bounded directly [R0_min*(1-buffer), R0_max*(1+buffer)]
      #' E.g., for R0 ∈ [1, 10]: bounds are [0.9, 11.0] instead of [0.06, 171]
      
      R0_lower <- R0_min * (1 - buffer)
      R0_upper <- R0_max * (1 + buffer)
      rt_lower <- rt_min * (1 - buffer)
      rt_upper <- rt_max * (1 + buffer)
      
      list(
        lower = c(R0 = R0_lower, recovery_time = rt_lower),
        upper = c(R0 = R0_upper, recovery_time = rt_upper)
      )
    },
    
    compute_starting_values = function(strategy, true_R0, true_recovery_time,
                                      R0_min, R0_max, rt_min, rt_max) {
      #' Compute starting values based on strategy
      
      if (strategy == "prior_mean") {
        R0_start <- exp((log(R0_min) + log(R0_max)) / 2)
        rt_start <- (rt_min + rt_max) / 2
      } else if (strategy == "true_values") {
        R0_start <- true_R0
        rt_start <- true_recovery_time
      } else {
        stop(sprintf("Unknown starting value strategy: %s", strategy))
      }
      
      list(R0 = R0_start, recovery_time = rt_start)
    },
    
    transform_estimates = function(estimates) {
      #' Transform optimization results to standard output format
      
      R0 <- estimates["R0"]
      recovery_time <- estimates["recovery_time"]
      gamma <- 1 / recovery_time
      beta <- R0 * gamma
      
      list(
        beta = as.numeric(beta),
        gamma = as.numeric(gamma),
        R0 = as.numeric(R0),
        recovery_time = as.numeric(recovery_time)
      )
    }
  )
}
