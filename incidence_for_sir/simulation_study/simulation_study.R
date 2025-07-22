# Comparing the brute-force computational approach to the taylor approximation via simulation.
## Author: Alexander C. Murph
## Date: March 2024

library(nleqslv)
library(ggplot2)
library(deSolve)
library(pracma)
library(reshape)
library(parallel)
library(doParallel)
theme_set(theme_bw())
setwd(this.path::here())
source('incidence_to_sir_helpers.R')

set.seed(13)
# Starting alpha,beta for optim calls:
alpha_initial                = 1.1788322 / 2
beta_initial                 = 0.8851870 / 2
# Note that i0 can't be greater than the above peak height...
N          = 10000
s0         = 0.95
i0         = 0.003
r0         = 1 - i0 - s0

optim_method = "Nelder-Mead"

num_of_sims     = 1000
simulation_data = NULL
parfctn = function(iter){
# for(iter in 1:num_of_sims){
  library(nleqslv)
  library(ggplot2)
  library(deSolve)
  library(pracma)
  library(reshape)
  print(paste("Starting sim number", iter))
  # Begin a given iteration of this simulation by chosing a random peak height and
  # time.
  # Chose the initial parameters based on dave's paper:
  peak_height = runif(1, min = 0.005, max = 0.035)
  peak_time   = sample(5:30, size = 1)
  
  ####################################################################################
  ####################################################################################
  ## In the following, I'm getting a reasonable peak incidence value + time.
  ## None of these parameter values are actually used (since they're the 'true values')
  # The following will give the value of rho directly:
  minimize_function          = function(x){ ((s0+i0) - peak_height - 1/x - log(s0*x)/x )**2 }
  rho_of_max_prevalence      = optim(1.2, minimize_function, method = "Brent", lower = 1, upper = 100)$par
  tao_peak_times_alpha       = log(rho_of_max_prevalence*s0)
  temp_integrand             = function(x){reparmeterized_time_integrand_sans_alpha(x, s0, i0, rho_of_max_prevalence)}
  peak_time_times_alpha      = integrate(temp_integrand, 0, tao_peak_times_alpha)$value
  lit_beta                   = peak_time_times_alpha / (peak_time*rho_of_max_prevalence)
  lit_lambda                 = lit_beta * rho_of_max_prevalence
  
  # There are so many different naming conventions across different papers, I 
  # am just going to include this generally.
  lit_beta                   = lit_beta
  lit_lambda                 = lit_lambda
  lit_alpha                  = lit_lambda
  dave_beta                  = lit_lambda
  dave_gamma                 = lit_beta 
  
  original_alphabetarho = c(lit_alpha, lit_beta, lit_alpha/lit_beta)
  
  ## Run the SIR model and get the values for peak incidence / peak incidence time.
  df1                      = sir(beta = original_alphabetarho[1], gamma = original_alphabetarho[2], S0 = s0, I0 = i0, R0 = r0, times = 0:50)
  temp_incidences          = df1$incidence
  temp_incidences[1]       = 0
  max_incidence_time       = which.max(temp_incidences) - 1
  max_incidence_value      = max(temp_incidences)
  
  true_max_incidence_time  = max_incidence_time
  true_max_incidence_value = max_incidence_value
  
  observed_peak_incidence      = max_incidence_value
  observed_incidence_peak_time = max_incidence_time
  ####################################################################################
  ####################################################################################
  
  ####################################################################################
  ####################################################################################
  ## In the following, I get the computational (brute-force) solution.
  temp_fn          = function(x){(fn1(x[1], x[2], observed_incidence_peak_time, observed_peak_incidence, s0, i0) + 
                                    fn2(x[1], x[2], observed_incidence_peak_time, observed_peak_incidence, s0, i0)) }
  
  start_time            = proc.time()
  # xx                    = fminunc(c(alpha_initial, beta_initial), temp_fn)
  xx                    = optim(c(alpha_initial, beta_initial), temp_fn, method = optim_method)
  end_time              = proc.time() - start_time
  computational_runtime = as.numeric(end_time[1])
  
  calculated_alpha = xx$par[1]
  calculated_beta  = xx$par[2]
  
  ### Now let's graph the SIR curve from these parameters and compare it to the initial time and height of the
  ### incidence peak.
  dave_beta  = calculated_alpha
  dave_gamma = calculated_beta
  lit_lambda = dave_beta
  lit_beta   = dave_gamma
  lit_rho    = lit_lambda / lit_beta
  
  ## run an SIR model
  df2         = sir(beta = dave_beta, gamma = dave_gamma, S0 = s0, I0 = i0, R0 = r0, times = 0:50)
  # Grab the actual incidence qois
  temp_incidences          = df2$incidence
  temp_incidences[1]       = 0
  
  max_incidence_time_computational  = which.max(temp_incidences) - 2
  max_incidence_value_computational = max(temp_incidences)
  ####################################################################################
  ####################################################################################
  
  
  ####################################################################################
  ####################################################################################
  ## In the following, I get the taylor-approximation solution.
  temp_fn          = function(x){(fn1_taylor(x[1], x[2], observed_incidence_peak_time, observed_peak_incidence, s0, i0) + 
                                    fn2_taylor(x[1], x[2], observed_incidence_peak_time, observed_peak_incidence, s0, i0)) }
  start_time       = proc.time()
  # xx               = fminunc(c(alpha_initial, beta_initial), temp_fn)
  xx               = optim(c(alpha_initial, beta_initial), temp_fn, method = optim_method)
  end_time         = proc.time() - start_time
  taylor_runtime   = as.numeric(end_time[1])
  
  calculated_alpha = xx$par[1]
  calculated_beta  = xx$par[2]
  
  ### Now let's graph the SIR curve from these parameters and compare it to the initial time and height of the
  ### incidence peak.
  dave_beta  = calculated_alpha
  dave_gamma = calculated_beta
  lit_lambda = dave_beta
  lit_beta   = dave_gamma
  lit_rho    = lit_lambda / lit_beta
  
  ## run an SIR model
  df2         = sir(beta = dave_beta, gamma = dave_gamma, S0 = s0, I0 = i0, R0 = r0, times = 0:50)
  # Grab the actual incidence qois
  temp_incidences          = df2$incidence
  temp_incidences[1]       = 0
  max_incidence_time_taylor   = which.max(temp_incidences) - 2
  max_incidence_value_taylor  = max(temp_incidences)
  ####################################################################################
  ####################################################################################
  
  
  ####################################################################################
  ####################################################################################
  ## In the following, I get the full ode solution
  start_time       = proc.time()
  xx               = ode_approx(alpha_initial, beta_initial, observed_incidence_peak_time, observed_peak_incidence, s0, i0)
  end_time         = proc.time() - start_time
  fullode_runtime  = as.numeric(end_time[1])
  
  calculated_alpha = xx$alpha
  calculated_beta  = xx$beta
  
  ### Now let's graph the SIR curve from these parameters and compare it to the initial time and height of the
  ### incidence peak.
  dave_beta  = calculated_alpha
  dave_gamma = calculated_beta
  lit_lambda = dave_beta
  lit_beta   = dave_gamma
  lit_rho    = lit_lambda / lit_beta
  
  ## run an SIR model
  df2         = sir(beta = dave_beta, gamma = dave_gamma, S0 = s0, I0 = i0, R0 = r0, times = 0:50)
  # Grab the actual incidence qois
  temp_incidences             = df2$incidence
  temp_incidences[1]          = 0
  max_incidence_time_fullode  = which.max(temp_incidences) - 2
  max_incidence_value_fullode = max(temp_incidences)
  if(is.nan(max_incidence_value_fullode)) browser()
  ####################################################################################
  ####################################################################################
  
  ####################################################################################
  ####################################################################################
  ## In the following, I get brute force + ode approx solution
  temp_fn            = function(x){(fn1_ode(x[1], x[2], observed_incidence_peak_time, observed_peak_incidence, s0, i0) + 
                                      fn2_ode(x[1], x[2], observed_incidence_peak_time, observed_peak_incidence, s0, i0)) }
  
  start_time         = proc.time()
  # xx                 = fminunc(c(alpha_initial, beta_initial), temp_fn)
  xx                 = optim(c(alpha_initial, beta_initial), temp_fn, method = optim_method)
  end_time           = proc.time() - start_time
  partialode_runtime = as.numeric(end_time[1])
  calculated_alpha   = xx$par[1]
  calculated_beta    = xx$par[2]
  
  ### Now let's graph the SIR curve from these parameters and compare it to the initial time and height of the
  ### incidence peak.
  dave_beta  = calculated_alpha
  dave_gamma = calculated_beta
  lit_lambda = dave_beta
  lit_beta   = dave_gamma
  lit_rho    = lit_lambda / lit_beta
  
  ## run an SIR model
  df2         = sir(beta = dave_beta, gamma = dave_gamma, S0 = s0, I0 = i0, R0 = r0, times = 0:50)
  # Grab the actual incidence qois
  temp_incidences          = df2$incidence
  temp_incidences[1]       = 0
  max_incidence_time_partialode   = which.max(temp_incidences) - 2
  max_incidence_value_partialode  = max(temp_incidences)
  ####################################################################################
  ####################################################################################
  
  
  # Record the results.
  temp_row = data.frame(computational_time_error  = abs(true_max_incidence_time - max_incidence_time_computational),
                        taylor_time_error         = abs(true_max_incidence_time - max_incidence_time_taylor),
                        fullode_time_error        = abs(true_max_incidence_time - max_incidence_time_fullode),
                        partialode_time_error     = abs(true_max_incidence_time - max_incidence_time_partialode),
                        computational_value_error = abs(true_max_incidence_value - max_incidence_value_computational),
                        taylor_value_error        = abs(true_max_incidence_value - max_incidence_value_taylor),
                        fullode_value_error       = abs(true_max_incidence_value - max_incidence_value_fullode),
                        partialode_value_error    = abs(true_max_incidence_value - max_incidence_value_partialode),
                        computational_time        = computational_runtime,
                        taylor_time               = taylor_runtime,
                        fullode_time              = fullode_runtime,
                        partialode_time           = partialode_runtime
  )
  return(temp_row)
  # simulation_data = rbind(simulation_data, temp_row)
}



sockettype <- "PSOCK"
ncores = 10
cl <- parallel::makeCluster(spec = ncores,type = sockettype) #, outfile=""
setDefaultCluster(cl)
registerDoParallel(cl)
sim_ts <- foreach(i=1:num_of_sims,
                  .verbose = T,
                  .combine = rbind)%dopar%{
                    print(i)
                    parfctn(i)
                  }
stopCluster(cl)  


write.csv(sim_ts, file = 'simulation_data.csv')
# sim_data = read.csv("simulation_data.csv")

