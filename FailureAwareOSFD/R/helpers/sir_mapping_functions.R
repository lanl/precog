# depends: 
fn1                    = function(alpha, beta, peak_incidence_time, peak_incidence_value, s0, i0, max_reproduction_number){
  tao                  = approx_proper_time_2(peak_incidence_time, alpha, beta, s0, i0)
  tao[is.nan(tao)]     = 0
  return_vector        = ifelse((alpha/beta > max_reproduction_number), 1e3, 
                                (alpha * ( s0 * exp(-alpha * tao) ) * ( (s0+i0) - s0*exp(-alpha*tao) - beta*tao ) - peak_incidence_value)**2)
  return_vector        = ifelse( (((s0+i0) - s0*exp(-alpha*tao) - beta*tao)<=0) , 1e3, 
                                 return_vector)
  return(return_vector)
}

fn2                    = function(alpha, beta, peak_incidence_time, peak_incidence_value, s0, i0, max_reproduction_number){
  tao                  = approx_proper_time_2(peak_incidence_time, alpha, beta, s0, i0)
  tao[is.nan(tao)]     = 0
  return_vector        = ifelse((alpha/beta > max_reproduction_number), 1e3, 
                                (-(s0 + i0) + beta * tao + 2 * s0 * exp(-alpha*tao) - beta / alpha)**2)
  return_vector        = ifelse( (((s0+i0) - s0*exp(-alpha*tao) - beta*tao)<=0) , 1e3, 
                                 return_vector)
  return(return_vector)
}

stochastic_sir_discrete <- function(beta, gamma, S0, I0, R0 = 0,
                                    t_max = 400, dt = 1) {
  times <- seq(0, t_max, by = dt)
  n <- length(times)

  S <- numeric(n)
  I <- numeric(n)
  R <- numeric(n)

  S[1] <- S0
  I[1] <- I0
  R[1] <- R0

  N <- S0 + I0 + R0

  for (k in 2:n) {
    if (I[k - 1] <= 0) {
      S[k] <- S[k - 1]
      I[k] <- I[k - 1]
      R[k] <- R[k - 1]
      next
    }

    p_inf <- 1 - exp(-beta * I[k - 1] / N * dt)
    p_rec <- 1 - exp(-gamma * dt)

    new_inf <- rbinom(1, size = S[k - 1], prob = p_inf)
    new_rec <- rbinom(1, size = I[k - 1], prob = p_rec)

    S[k] <- S[k - 1] - new_inf
    I[k] <- I[k - 1] + new_inf - new_rec
    R[k] <- R[k - 1] + new_rec
  }

  # Calculate incidence as -diff(S), matching deterministic sir_sirstuff
  incidence <- c(NA, -diff(S))

  data.frame(time = times, S = S, I = I, R = R, incidence = incidence)
}

approx_proper_time_2   = function(t, alpha, beta, s0, i0){
  # It is likely that the integrand will have an asymptote somewhere in the range of tao.
  # In the following, I should only consider tao such that i0+s0-s0e... is above zero.
  temp_integrand     = function(x){time_integrand(x, s0, i0, alpha, beta)}
  taos               = 1:1000/1000
  inverse_integrands = 1/(temp_integrand(taos))
  
  # In the following, I should only consider tao such that i0+s0-s0e... is above zero.
  taos               = 1:1000/1000
  inverse_integrands = 1/(sapply(taos, temp_integrand))
  # These should stay above zero.  If they don't, we have to upper bound our taos.
  indices_below_zero = which( inverse_integrands <= 0 )
  if(length(indices_below_zero)==0){
    upper_tao = 1
  } else if((c(1)%in%indices_below_zero)){
    return(0)
  } else {
    upper_tao = taos[min(indices_below_zero)-1]
  }
  
  gradient            = function(x){ time_integrand(x, s0, i0, alpha, beta) }
  gradient_diff       = function(x){ sapply(x, integrator) }
  
  integrator          = function(x){ abs(integrate(temp_integrand, 0, x)$value - t) }
  integrator_diff     = function(x){ sapply(x, integrator) }
  
  tao                 = optim(par = upper_tao/2, fn = integrator_diff,
                              gr = gradient_diff, method = "Brent",
                              lower = 0, upper = upper_tao)$par
  
  # tao                 = optimize(f = integrator_diff,
  #                             lower = 0, upper = upper_tao)$minimum
  
  
  return(tao)
}

time_integrand         = function(tao, s0, i0, alpha, beta){
  return( 1/( (i0+s0) - beta * tao - s0*exp(-alpha*tao) ) )
}

sir_sirstuff                    = function(beta, gamma, S0, I0, R0, times) {
  # the differential equations:
  sir_equations <- function(time, variables, parameters) {
    with(as.list(c(variables, parameters)), {
      dS <- -beta * I * S
      dI <-  beta * I * S - gamma * I
      dR <-  gamma * I
      return(list(c(dS, dI, dR)))
    })
  }
  
  # the parameters values:
  parameters_values <- c(beta  = beta, gamma = gamma)
  
  # the initial values of variables:
  initial_values <- c(S = S0, I = I0, R = R0)
  
  # solving
  out <- data.frame(deSolve::ode(initial_values, times, sir_equations, 
                      parameters_values, method = "lsoda",
                      rtol = 1e-10,
                      atol = 1e-12))
  
  # returning the output:
  as.data.frame(cbind(out, incidence = c(NA,-diff(out$S)), beta = beta, gamma = gamma))
}

#' @export
calculate_alphabeta_from_PIVPIT = function(observed_incidence_peak_time, observed_peak_incidence, s0, 
                                      pit_bounds = c(5, 50),
                                      piv_bounds = c(0.005, 0.6),
                                      s0_bounds = c(0.95, 0.999), 
                                      reproduction_number_bounds = c(1.001, 100)) {

  max_reproduction_number = max(reproduction_number_bounds)
  observed_incidence_peak_time = pit_bounds[1] + observed_incidence_peak_time * diff(pit_bounds)
  observed_peak_incidence = piv_bounds[1] + observed_peak_incidence * diff(piv_bounds)
  s0 = s0_bounds[1] + s0 * diff(s0_bounds)

  alpha_initial                = 1.1788322 / 2
  beta_initial                 = 0.8851870 / 2
  optim_method = "Nelder-Mead"
  N          = 10000
  # s0         = 0.95
  i0         = 0.001
  r0         = 1 - i0 - s0
  temp_fn          = function(x){(fn1(x[1], x[2], observed_incidence_peak_time, observed_peak_incidence, s0, i0, max_reproduction_number) + 
                                  fn2(x[1], x[2], observed_incidence_peak_time, observed_peak_incidence, s0, i0, max_reproduction_number)) }
  xx                    = optim(c(alpha_initial, beta_initial), temp_fn, method = optim_method)
  calculated_alpha = xx$par[1]
  calculated_beta  = xx$par[2]
  return(data.frame(alpha = calculated_alpha, beta = calculated_beta, s0 = s0))
}

#' @export
calculate_sir_from_PIVPIT = function(observed_incidence_peak_time, observed_peak_incidence, s0,
                                      pit_bounds = c(5, 50),
                                      piv_bounds = c(0.005, 0.6),
                                      s0_bounds = c(0.95, 0.999),
                                      reproduction_number_bounds = c(1.001, 100)) {

  max_reproduction_number = max(reproduction_number_bounds)
  observed_incidence_peak_time = pit_bounds[1] + observed_incidence_peak_time * diff(pit_bounds)
  observed_peak_incidence = piv_bounds[1] + observed_peak_incidence * diff(piv_bounds)
  s0 = s0_bounds[1] + s0 * diff(s0_bounds)

  # Method of moments estimator from Amaro (2022) using Gumbel approximation
  # Based on equations (41-44) and (70) - proper time fit without time integral!

  # Step 1: Solve for rho from i_peak using equation (21): i_peak = 1 - (1 + ln(s0*rho))/rho
  estimate_rho_from_ipeak <- function(ipeak, s0_val) {
    objective <- function(rho) {
      (rho * (1 - ipeak) - 1 - log(s0_val * rho))^2
    }
    result <- optimize(objective, interval = c(1.01, max_reproduction_number))
    return(result$minimum)
  }

  # Get rho from peak incidence
  rho_estimated <- estimate_rho_from_ipeak(observed_peak_incidence, s0)

  # Step 2: Use equations (43-44) for Gumbel parameters
  # From equation (43): a = e*ln(s0*rho)/λ
  # From equation (44): a/b = e*(1 - (1 + ln(s0*rho))/rho) = e*i_peak
  # From equation (70): β*b = ln(rho)/(rho - 1 - ln(rho))

  # From (44): a/b = e * i_peak
  a_over_b <- exp(1) * observed_peak_incidence

  # From (70): β*b = ln(rho)/(rho - 1 - ln(rho))
  # Therefore: β = ln(rho)/(b*(rho - 1 - ln(rho)))
  # And: b ≈ βt_peak (approximate relationship for the Gumbel width parameter)

  # Using equation (70) and the relationship between b and peak time:
  # β*b ≈ ln(rho)/(rho - 1 - ln(rho))
  # We need to estimate b from the peak time. From the Gumbel fit,
  # b is related to the width/time scale of the epidemic.

  # Empirical relationship from paper: b is roughly proportional to the timescale
  # From equation (70): β = ln(rho)/(b*(rho - 1 - ln(rho)))
  # Rearranging: b = ln(rho)/(β*(rho - 1 - ln(rho)))

  # For initial estimate, we can use the fact that the peak time scales with 1/β
  # From the paper's Table 2 and Figure 9, b ≈ t_peak/C where C ≈ 1.5-2
  # Let's use a more principled approach from equation (70):

  # β*b = ln(rho)/(rho - 1 - ln(rho))
  # If we assume b is proportional to t_peak: b ≈ k*t_peak for some constant k
  # From the paper's numerical results, k ≈ 0.6-0.8 for typical ρ values
  # Let's use k ≈ 0.65 as a reasonable estimate

  b_estimate <- 0.65 * observed_incidence_peak_time
  beta_initial <- log(rho_estimated) / (b_estimate * (rho_estimated - 1 - log(rho_estimated)))

  # From equation (43) and ρ = λ/β:
  alpha_initial <- rho_estimated * beta_initial

  # Set up initial conditions
  i0 <- 0.001
  r0 <- 1 - i0 - s0

  # Box constraints for (alpha, rho)
  # alpha: [0.001, 100]
  # rho: [reproduction_number_bounds[1]/s0, reproduction_number_bounds[2]]
  #      This ensures s0*rho >= reproduction_number_bounds[1]
  alpha_lower <- 0.001
  alpha_upper <- 100
  rho_lower <- reproduction_number_bounds[1] / s0
  rho_upper <- reproduction_number_bounds[2]

  # Reparameterized optimization: optimize (alpha, rho) instead of (alpha, beta)
  # This makes constraints easier: rho directly constrained, beta = alpha/rho
  # Define objective in terms of (alpha, rho)
  obj_fn_reparameterized <- function(x) {
    alpha <- x[1]
    rho <- x[2]
    beta <- alpha / rho

    # Penalize negative or near-zero values
    if (alpha <= 0 || beta <= 0 || rho <= 0) {
      return(1e6)
    }

    # Simple objective: just sum the two components
    # L-BFGS-B will respect the box constraints automatically
    fn1(alpha, beta, observed_incidence_peak_time, observed_peak_incidence, s0, i0, max_reproduction_number) +
    fn2(alpha, beta, observed_incidence_peak_time, observed_peak_incidence, s0, i0, max_reproduction_number)
  }

  # Use L-BFGS-B for box-constrained optimization
  xx <- optim(
    par = c(alpha_initial, rho_estimated),
    fn = obj_fn_reparameterized,
    method = "L-BFGS-B",
    lower = c(alpha_lower, rho_lower),
    upper = c(alpha_upper, rho_upper)
  )
  calculated_alpha <- xx$par[1]
  calculated_rho <- xx$par[2]
  calculated_beta <- calculated_alpha / calculated_rho
  ### Now let's graph the SIR curve from these parameters and compare it to the initial time and height of the
  ### incidence peak.
  dave_beta  = calculated_alpha
  dave_gamma = calculated_beta
  df2         = sir_sirstuff(beta = dave_beta, gamma = dave_gamma, S0 = s0, I0 = i0, R0 = r0, times = 0:50)
  # Grab the actual incidence qois
  temp_incidences          = df2$incidence
  temp_incidences[1]       = 0
  max_incidence_time_computational  = which.max(temp_incidences) - 2
  max_incidence_value_computational = max(temp_incidences)

  # Check if computational PIV/PIT match the observed values closely
  # Relative tolerance of 10% for PIV and absolute tolerance of 2 time units for PIT
  piv_rel_error <- abs(max_incidence_value_computational - observed_peak_incidence) / observed_peak_incidence
  pit_abs_error <- abs(max_incidence_time_computational - observed_incidence_peak_time)

  # if (piv_rel_error > 0.5) {
  #   stop(sprintf("Computational PIV (%.6f) differs from observed PIV (%.6f) by %.1f%% (tolerance: 10%%)",
  #                max_incidence_value_computational, observed_peak_incidence, piv_rel_error * 100))
  # }

  # if (pit_abs_error > 5) {
  #   stop(sprintf("Computational PIT (%.1f) differs from observed PIT (%.1f) by %.1f time units (tolerance: 2)",
  #                max_incidence_time_computational, observed_incidence_peak_time, pit_abs_error))
  # }

  # # With the reparameterized optimization, s0*rho should always be feasible due to box constraints
  # # But check anyway as a safety measure
  # tmp_rho = calculated_rho
  # if(s0*tmp_rho < reproduction_number_bounds[1]) {
  #   return(c(piv = NA_real_, pit = NA_real_, alpha = NA_real_, beta = NA_real_,
  #            calculate_alpha = NA_real_, calculate_beta = NA_real_, s0 = NA_real_))
  # }

  return(c(piv = max_incidence_value_computational,
           pit = max_incidence_time_computational,
           alpha = dave_beta,
           beta = dave_gamma,
           calculate_alpha = calculated_alpha,
           calculate_beta = calculated_beta,
           s0 = s0))
}

#' @export
calculate_sir_for_outputSpaceFilling = function(alpha_initial, rho_initial, s0_initial,
                                        alpha_bounds = c(0.01, 5),
                                        reproduction_number_bounds = c(1.01, 30),
                                        s0_bounds = c(0.95, 0.999)){
  alpha = alpha_bounds[1] + alpha_initial * diff(alpha_bounds)
  s0 = s0_bounds[1] + s0_initial * diff(s0_bounds)
  rho = reproduction_number_bounds[1] + rho_initial * diff(reproduction_number_bounds)
  beta = alpha / rho

  N          = 10000
  # s0         = 0.95
  i0         = 0.001
  r0         = 1 - i0 - s0

  dave_beta  = alpha
  dave_gamma = beta
  df2         = sir_sirstuff(beta = dave_beta, gamma = dave_gamma, S0 = s0, I0 = i0, R0 = r0, times = 0:50)
  # Grab the actual incidence qois
  temp_incidences          = df2$incidence
  temp_incidences[1]       = 0
  max_incidence_time_computational  = which.max(temp_incidences) - 2
  max_incidence_value_computational = max(temp_incidences)

  tmp_rho = dave_beta / dave_gamma
  # if(s0*tmp_rho <= reproduction_number_bounds[1]) return(NA)

  return(list(piv = max_incidence_value_computational, pit = max_incidence_time_computational))
}

#' @export
calculate_stochastic_sir_for_outputSpaceFilling = function(alpha_initial, rho_initial, s0_initial,
                                        alpha_bounds = c(0.01, 5),
                                        reproduction_number_bounds = c(1.01, 30),
                                        s0_bounds = c(0.95, 0.999)){
  alpha = alpha_bounds[1] + alpha_initial * diff(alpha_bounds)
  s0 = s0_bounds[1] + s0_initial * diff(s0_bounds)
  rho = reproduction_number_bounds[1] + rho_initial * diff(reproduction_number_bounds)
  beta = alpha / rho

  N          = 10000
  # s0         = 0.95
  i0         = 0.001
  r0         = 1 - i0 - s0

  dave_beta  = alpha
  dave_gamma = beta
  S0 = floor(N*s0)
  I0 = floor(N*i0)
  R0 = N - S0 - I0

  df2 = stochastic_sir_discrete(dave_beta, dave_gamma, S0, I0, R0,
                                    t_max = 400, dt = 1)

  # Grab the actual incidence qois
  temp_incidences          = df2$incidence
  temp_incidences[1]       = 0
  tmp_idx = which.max(temp_incidences)
  max_incidence_time_computational  = df2$time[tmp_idx]
  max_incidence_value_computational = max(temp_incidences)

  tmp_rho = dave_beta / dave_gamma
  if(s0*tmp_rho <= reproduction_number_bounds[1]) return(NA)

  return(list(piv = max_incidence_value_computational, pit = max_incidence_time_computational))
}




#' @export
calculate_sir_from_PIVPIT_by_regime = function(observed_incidence_peak_time, observed_peak_incidence, s0, epsilon,
                                      regime_name = "Surge",
                                      snippet_length = 5,
                                      pit_bounds = c(5, 50),
                                      piv_bounds = c(0.005, 0.6),
                                      s0_bounds = c(0.95, 0.999),
                                      epsilon_bounds = c(0.01,2),
                                      max_reproduction_number = 50,
                                      match_length = 5) {

  observed_incidence_peak_time = pit_bounds[1] + observed_incidence_peak_time * diff(pit_bounds)
  observed_peak_incidence = piv_bounds[1] + observed_peak_incidence * diff(piv_bounds)
  s0 = s0_bounds[1] + s0 * diff(s0_bounds)
  epsilon = epsilon_bounds[1] + epsilon * diff(epsilon_bounds)

  alpha_initial                = 1.1788322 / 2
  beta_initial                 = 0.8851870 / 2
  optim_method = "Nelder-Mead"
  N          = 10000
  i0         = 0.001
  r0         = 1 - i0 - s0
  temp_fn          = function(x){(fn1(x[1], x[2], observed_incidence_peak_time, observed_peak_incidence, s0, i0, max_reproduction_number) + 
                                  fn2(x[1], x[2], observed_incidence_peak_time, observed_peak_incidence, s0, i0, max_reproduction_number)) }
  xx                    = optim(c(alpha_initial, beta_initial), temp_fn, method = optim_method)
  calculated_alpha = xx$par[1]
  calculated_beta  = xx$par[2]

  ### Now let's graph the SIR curve from these parameters and compare it to the initial time and height of the
  ### incidence peak.
  dave_beta  = calculated_alpha
  dave_gamma = calculated_beta
  df2         = sir_sirstuff(beta = dave_beta, gamma = dave_gamma, S0 = s0, I0 = i0, R0 = r0, times = 0:50)        

  I = df2$I
  xx = classify_entire_timeseries(I, snippet_length, match_length)
  first_index <- which(xx$regime == regime_name)[1]
  v = xx$value[first_index:(first_index+snippet_length-1)]

  I_noisy <- pmax(v, 1e-8) * (1 + rnorm(length(v), mean = 0, sd = epsilon))
  I_noisy <- pmax(I_noisy, 0) 

  result <- as.list(I_noisy)
  names(result) <- as.character(seq_along(I_noisy))
  result[["calculated_alpha"]] = calculated_alpha
  result[["calculated_beta"]] = calculated_beta
  return(result)
}

#' @export
calculate_sir_for_outputSpaceFilling_by_regime = function(alpha_initial, rho_initial, s0_initial, epsilon_initial, 
                                        regime_name = "Surge",
                                        snippet_length = 5,
                                        alpha_bounds = c(0.01, 5),
                                        reproduction_number_bounds = c(1.01, 30),
                                        s0_bounds = c(0.95, 0.999), 
                                        # epsilon_bounds = c(0.01,2),
                                        match_length = 4,
                                        method = "ISFD"){
  if(match_length > snippet_length) stop("match length must be leq snippet_length.")
  alpha = alpha_bounds[1] + alpha_initial * diff(alpha_bounds)
  s0 = s0_bounds[1] + s0_initial * diff(s0_bounds)
  rho = reproduction_number_bounds[1] + rho_initial * diff(reproduction_number_bounds)
  # epsilon = epsilon_bounds[1] + epsilon_initial * diff(epsilon_bounds)

  beta = alpha / rho

  N          = 10000
  i0         = 0.001
  r0         = 1 - i0 - s0
  # dave_beta  = alpha
  # dave_gamma = beta
  # df2         = sir_sirstuff(beta = dave_beta, gamma = dave_gamma, S0 = s0, I0 = i0, R0 = r0, times = 0:50)

  dave_beta  = alpha
  dave_gamma = beta
  S0 = floor(N*s0)
  I0 = floor(N*i0)
  R0 = N - S0 - I0

  df2 = stochastic_sir_discrete(dave_beta, dave_gamma, S0, I0, R0,
                                    t_max = 400, dt = 1)

  # Grab the actual incidence qois
  I = df2$incidence
  I[1] = 0  # Set first value to 0 (was NA)
  # Remove continuous leading and trailing zeros, then add padding of 9 zeros on each end

  # Find first non-zero index
  first_nonzero <- which(I != 0)[1]
  # Find last non-zero index
  last_nonzero <- which(I != 0)[length(which(I != 0))]

  # Handle edge cases where all values are zero or no zeros exist
  if (is.na(first_nonzero) || is.na(last_nonzero)) {
    # If all zeros or something went wrong, keep original I
    I_trimmed <- I
  } else {
    # Trim to non-zero range
    I_trimmed <- I[first_nonzero:last_nonzero]
    # Add padding of 9 zeros on each end
    I_trimmed <- c(rep(0, 9), I_trimmed, rep(0, 9))
  }

  I <- I_trimmed

  xx = classify_entire_timeseries(I, snippet_length, match_length)

  # Generate random 10 digit number for unique filename
  rand <- sample.int(9999999999, 1)

  # Create directory if it doesn't exist
  output_dir <- here::here('data', 'stochastic_sir')
  if (!dir.exists(output_dir)) {
    dir.create(output_dir, recursive = TRUE)
  }

  # Build a data frame with each row being a complete snippet
  snippet_data <- NULL
  for (i in 1:(length(I) - snippet_length + 1)) {
    snippet_values <- I[i:(i + snippet_length - 1)]
    snippet_regime <- xx$regime[i]  # Regime assignment for first value of snippet

    # Create a row with the regime and all snippet values
    row_data <- c(regime = snippet_regime, snippet_values)
    snippet_data <- rbind(snippet_data, row_data)
  }

  # Convert to data frame with proper column names
  snippet_df <- as.data.frame(snippet_data)
  colnames(snippet_df) <- c("regime", paste0("t", 1:snippet_length))

  # Write all snippets to file
  output_file <- here::here('data', 'stochastic_sir', paste0(rand, '_', method, '.csv'))
  write.csv(snippet_df, output_file, row.names = FALSE)

  # Still return the snippet matching regime_name for backward compatibility
  # IMPORTANT: Randomly sample from all occurrences to get diversity, not just the first!
  regime_indices <- which(xx$regime == regime_name)

  if(length(regime_indices) == 0) {
    # Regime not found in this SIR curve - this is expected sometimes
    # Return NA to signal failure so OSFD can retry with different parameters
    stop(sprintf("Regime '%s' not found in generated SIR curve. Found regimes: %s",
                 regime_name, paste(unique(xx$regime), collapse=", ")))
  }

  # Randomly sample from valid starting indices (ensure we have room for full snippet)
  valid_indices <- regime_indices[regime_indices <= (nrow(xx) - snippet_length + 1)]
  if(length(valid_indices) == 0) {
    stop(sprintf("Regime '%s' found but no valid snippet locations", regime_name))
  }

  first_index <- sample(valid_indices, 1)
  v = xx$value[first_index:(first_index+snippet_length-1)]
  I_noisy <- pmax(v, 1e-8) #* (1 + rnorm(length(v), mean = 0, sd = epsilon))

  result <- as.list(I_noisy)
  names(result) <- as.character(seq_along(I_noisy))
  return(result)
}
