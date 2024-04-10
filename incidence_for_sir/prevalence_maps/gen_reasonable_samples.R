## Testing the process of going from prevalence peak+time to SIR parameters.
#### Author: Alexander C. Murph
library(ggplot2)
library(deSolve)
theme_set(theme_bw())
setwd("~/GitLab/incidence_for_sir")

# Note a change in parameters from how Dave does it to how
# I've found it down in the literature:
# dave <-> literature
# beta <-> lambda/N
# gamma <-> beta
# in lit: rho = lambda / beta

## define SIRfunctions
sir <- function(beta, gamma, S0, I0, R0, times) {
  require(deSolve) # for the "ode" function
  
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

time_integrand = function(tao, s0, rho){
  return( ( 1/(1 - tao - s0*exp(-rho*tao)) ) )
}

### Murph update 2023:
# I want to add in the peak time and height to this simulation to see how well
# it works.
N           = 10000
s0          = .999
# peak_height = 0.065
# peak_time   = 17
peak_height = runif(1, min = 0.005, max = 0.035)
peak_time   = sample(5:30, size = 1)

minimize_function     = function(x){ (1 - peak_height - (1+log(s0*x))/x)**2 }
rho_of_max_prevalence = optim(0.5, minimize_function, method = "Brent", lower = 1, upper = 100)$par
tao_peak              = log(rho_of_max_prevalence*s0) / rho_of_max_prevalence
temp_integrand        = function(x){time_integrand(x, s0, rho_of_max_prevalence)}
peak_time_sans_beta   = integrate(temp_integrand, 0, tao_peak)$value

lit_beta   = peak_time_sans_beta / peak_time
lit_lambda = lit_beta * rho_of_max_prevalence

dave_beta  = lit_lambda
dave_gamma = lit_beta

reasonable_samples = NULL
for(samp in 1:500){
  peak_height = runif(1, min = 0.005, max = 0.035)
  peak_time   = sample(5:30, size = 1)
  
  minimize_function     = function(x){ (1 - peak_height - (1+log(s0*x))/x)**2 }
  rho_of_max_prevalence = optim(0.5, minimize_function, method = "Brent", lower = 1, upper = 100)$par
  tao_peak              = log(rho_of_max_prevalence*s0) / rho_of_max_prevalence
  temp_integrand        = function(x){time_integrand(x, s0, rho_of_max_prevalence)}
  peak_time_sans_beta   = integrate(temp_integrand, 0, tao_peak)$value
  
  lit_beta   = peak_time_sans_beta / peak_time
  lit_lambda = lit_beta * rho_of_max_prevalence
  
  dave_beta  = lit_lambda
  dave_gamma = lit_beta
  temp = data.frame(samp_num = samp, beta = dave_beta, gamma = dave_gamma)
  reasonable_samples = rbind(reasonable_samples, temp)
}

