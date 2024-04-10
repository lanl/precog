# Incidence to sir params helper functions
## Author: Alexander C. Murph
## Date: March 2024
library(deSolve) # for the "ode" function

time_integrand         = function(tao, s0, i0, alpha, beta){
  return( 1/( (i0+s0) - beta * tao - s0*exp(-alpha*tao) ) )
}

sir                    = function(beta, gamma, S0, I0, R0, times) {
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
  out <- data.frame(ode(initial_values, times, sir_equations, parameters_values, method = "rk4"))
  
  # returning the output:
  as.data.frame(cbind(out, incidence = c(NA,-diff(out$S)), beta = beta, gamma = gamma))
}

solve_tao_ode          = function(alpha, beta, S0, I0, R0, PIT, times) {
  # the differential equations:
  sir_equations = function(time, variables, parameters) {
    with(as.list(c(variables, parameters)), {
      dR =  beta * (1 - s0*exp(-(R-r0)*(alpha/beta)) - R)
      return(list(c(dR)))
    })
  }
  
  # the parameters values:
  parameters_values <- c(alpha  = alpha, beta = beta,
                         s0    = S0,   r0    = R0)
  
  # the initial values of variables:
  initial_values <- c(R = r0)
  
  # solving
  out = data.frame(ode(initial_values, 0:PIT, sir_equations, parameters_values, method = "rk4"))
  # out = out[complete.cases(out),]
  
  # returning the output:
  return(out)
}

sir_PIT                = function(beta, gamma, S0, I0, R0, times) {
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
  out <- data.frame(ode(initial_values, times, sir_equations, parameters_values, method = "rk4"))
  
  # returning the output:
  df2                      = as.data.frame(cbind(out, incidence = c(NA,-diff(out$S)), beta = beta, gamma = gamma))
  # Grab the actual incidence qois
  temp_incidences          = df2$incidence
  temp_incidences[1]       = 0
  max_incidence_time_df2   = which.max(temp_incidences) - 1
  max_incidence_value_df2  = max(temp_incidences)
  
  return(max_incidence_time_df2)
}

sir_qois               = function(beta, gamma, S0, I0, R0, times) {
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
  out <- data.frame(ode(initial_values, times, sir_equations, parameters_values, method = "rk4"))
  
  # returning the output:
  df2                      = as.data.frame(cbind(out, incidence = c(NA,-diff(out$S)), beta = beta, gamma = gamma))
  df2                      = df2[complete.cases(df2),]
  # Grab the actual incidence qois
  temp_incidences          = df2$incidence
  temp_incidences[1]       = 0
  max_incidence_time_df2   = which.max(temp_incidences) - 1
  max_incidence_value_df2  = max(temp_incidences)
  
  if(is.nan(max_incidence_time_df2)){
    return(list(max_time = 1, max_value = 0.5))
  }
  
  if(is.nan(max_incidence_value_df2)){
    return(list(max_time = 1, max_value = 0.5))
  }
  
  return(list(max_time = max_incidence_time_df2, max_value = max_incidence_value_df2))
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
  
  # browser()
  tao                 = optim(par = upper_tao/2, fn = integrator_diff,
                              gr = gradient_diff, method = "Brent",
                              lower = 0, upper = upper_tao)$par
  
  # tao                 = optimize(f = integrator_diff,
  #                             lower = 0, upper = upper_tao)$minimum
  
  
  return(tao)
}

approx_proper_time_ode = function(t, alpha, beta, s0, i0){
  tao = solve_tao_ode(alpha, beta, s0, i0, (1-s0-i0), t, times = t)
  # print("----------------")
  # print(paste('the number of rows is', nrow(tao), "while the value of t is", t))
  # print(paste('alpha is', alpha, "while beta is", beta))
  # print(tao)
  tao = tao$R[t+1]
  # print(paste("the tao approx was:", tao))
  return(tao)
}

fn1                    = function(alpha, beta, peak_incidence_time, peak_incidence_value, s0, i0){
  tao                  = approx_proper_time_2(peak_incidence_time, alpha, beta, s0, i0)
  max_reproduction_num = 20
  tao[is.nan(tao)]     = 0
  return_vector        = ifelse((alpha/beta > max_reproduction_num), 1e3, 
                                (alpha * ( s0 * exp(-alpha * tao) ) * ( (s0+i0) - s0*exp(-alpha*tao) - beta*tao ) - peak_incidence_value)**2)
  return_vector        = ifelse( (((s0+i0) - s0*exp(-alpha*tao) - beta*tao)<=0) , 1e3, 
                                 return_vector)
  return_vector        = ifelse( (alpha > 5) , 1e3, 
                                 return_vector)
  return(return_vector)
}

fn2                    = function(alpha, beta, peak_incidence_time, peak_incidence_value, s0, i0){
  tao                  = approx_proper_time_2(peak_incidence_time, alpha, beta, s0, i0)
  max_reproduction_num = 20
  tao[is.nan(tao)]     = 0
  return_vector        = ifelse((alpha/beta > max_reproduction_num), 1e3, 
                                (-(s0 + i0) + beta * tao + 2 * s0 * exp(-alpha*tao) - beta / alpha)**2)
  return_vector        = ifelse( (((s0+i0) - s0*exp(-alpha*tao) - beta*tao)<=0) , 1e3, 
                                 return_vector)
  return_vector        = ifelse( (alpha > 5) , 1e3, 
                                 return_vector)
  return(return_vector)
}

fn1_ode                = function(alpha, beta, peak_incidence_time, peak_incidence_value, s0, i0){
  tao                  = approx_proper_time_ode(peak_incidence_time, alpha, beta, s0, i0)
  # print(paste("tao for 1 is:", tao))
  max_reproduction_num = 20
  tao[is.nan(tao)]     = 0
  return_vector        = ifelse(alpha/beta > max_reproduction_num, 1e3, 
                                (alpha * ( s0 * exp(-alpha * tao) ) * ( (s0+i0) - s0*exp(-alpha*tao) - beta*tao ) - peak_incidence_value)**2)
  return_vector        = ifelse( ((s0+i0) - s0*exp(-alpha*tao) - beta*tao)<=0 , 1e3, 
                                 return_vector)
  return(return_vector)
}

fn2_ode                = function(alpha, beta, peak_incidence_time, peak_incidence_value, s0, i0){
  tao                  = approx_proper_time_ode(peak_incidence_time, alpha, beta, s0, i0)
  # print(paste("tao for 2 is:", tao))
  tao[is.nan(tao)]     = 0
  max_reproduction_num = 20
  return_vector        = ifelse(alpha/beta > max_reproduction_num, 1e3, 
                                (-(s0 + i0) + beta * tao + 2 * s0 * exp(-alpha*tao) - beta / alpha)**2)
  return_vector        = ifelse( ((s0+i0) - s0*exp(-alpha*tao) - beta*tao)<=0 , 1e3, 
                                 return_vector)
  return(return_vector)
}

reparmeterized_time_integrand_sans_alpha = function(tao, s0, i0, rho){
  return( 1/( (i0 + s0) - tao/rho - s0*exp(-tao)) )
}

taylor_approx          = function(alpha, beta, peak_incidence_time, s0, i0){
  rho         = alpha / beta
  rho_temp    = ifelse( (s0*rho - 1)>=0.99, 0.99, (s0*rho - 1) )
  rho_temp    = ifelse(rho_temp<=-0.99, -0.99, rho_temp)
  if(rho_temp%in%c(-0.99, 0.99)){
    return(0)
  }
  rho         = (rho_temp + 1)/s0
  
  K           = sqrt((s0*rho - 1)^2 + 2*s0*i0*rho^2)
  phi         = atanh((1/K)*(s0*rho - 1))
  term1       = rho * s0 - 1
  term2       = (beta * K * peak_incidence_time) / 2 - phi
  proper_time = (alpha^2 / s0) * (term1 + K * tanh(term2))
  
  return( (proper_time + (1-s0-i0)) )
}

fn1_taylor             = function(alpha, beta, peak_incidence_time, peak_incidence_value, s0, i0){
  tao                  = taylor_approx(alpha, beta, peak_incidence_time, s0, i0)
  max_reproduction_num = 20
  return_vector        = ifelse( (alpha/beta > max_reproduction_num)|(tao==0) , 1e3, 
                                (alpha * ( s0 * exp(-alpha * tao) ) * ( (s0+i0) - s0*exp(-alpha*tao) - beta*tao ) - peak_incidence_value)**2)
  return_vector        = ifelse( ((s0+i0) - s0*exp(-alpha*tao) - beta*tao)<=0 , 1e3, 
                                 return_vector)
  return(return_vector)
}

fn2_taylor             = function(alpha, beta, peak_incidence_time, peak_incidence_value, s0, i0){
  tao                  = taylor_approx(alpha, beta, peak_incidence_time, s0, i0)
  max_reproduction_num = 20
  return_vector        = ifelse(alpha/beta > max_reproduction_num, 1e3, 
                                (-(s0 + i0) + beta * tao + 2 * s0 * exp(-alpha*tao) - beta / alpha)**2)
  return_vector        = ifelse( ((s0+i0) - s0*exp(-alpha*tao) - beta*tao)<=0 , 1e3, 
                                 return_vector)
  return(return_vector)
}

ode_approx             = function(alpha_initial, beta_initial, peak_incidence_time, peak_incidence_value, s0, i0){
  dave_beta   = alpha
  dave_gamma  = beta
  
  # It is likely that the integrand will have an asymptote somewhere in the range of tao.
  # In the following, I should only consider tao such that i0+s0-s0e... is above zero.
  full_ode_solve_fctn = function(x){
    max_reproduction_num = 20
    # full_ode             = ifelse( (x[1]/x[2] > max_reproduction_num), 0.5, 
    #                               sir_qois(x[1], x[2], s0, i0, (1-s0-i0), times = 0:40) )
    
    full_ode             = sir_qois(x[1], x[2], s0, i0, (1-s0-i0), times = 0:40)
    
    # return_vector        = ifelse( (((s0+i0) - s0*exp(-alpha*tao) - beta*tao)<=0) , 1e3, 
    #                                return_vector)
    # full_ode = sir_qois(x[1], x[2], s0, i0, (1-s0-i0), times = 0:40)
    return((full_ode$max_time - peak_incidence_time)**2 + (full_ode$max_value- peak_incidence_value)**2)
  }
  xx               = optim(par = c(alpha_initial, beta_initial), fn = full_ode_solve_fctn) # , lb = c(0,0), ub = c(Inf,Inf))
  # xx               = fmincon(x0 = c(alpha_initial, beta_initial), fn = full_ode_solve_fctn, lb = c(0,0), ub = c(Inf,Inf))
  calculated_alpha = xx$par[1]
  calculated_beta  = xx$par[2]

  return(list(alpha = calculated_alpha, beta = calculated_beta))
}



