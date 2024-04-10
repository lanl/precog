## Experiments comparing incidence time/peaks with prevalence time/peaks
## over a reasonable space for SIR parameters.
# Author: Alexander C. Murph
# Experiment suggestions by Dave Osthus.

library(lhs)
library(ggplot2)
library(gridExtra)
# Comment out the following line after the first run of this script.
source("gen_reasonable_samples.R")
theme_set(theme_bw())
set.seed(13)

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

time_integrand = function(tao, s0, rho, beta){
  return( ( 1/(1 - tao - s0*exp(-rho*tao)) )/beta )
}

### Murph update 2023:
# I want to add in the peak time and height to this simulation to see how well
# it works.
N          = 10000
s0         = .999
# dave_beta  = runif(1, min=0.25, max = 2)
# dave_gamma = runif(1, min=0.25, max = dave_beta)
# count      = count + 1
# dave_lower_beta  = min(reasonable_samples['beta'])
# dave_upper_beta  = max(reasonable_samples['beta'])
# dave_upper_gamma = max(reasonable_samples['gamma'])
# dave_lower_gamma = min(reasonable_samples['gamma'])

dave_lower_beta  = 0.25
dave_upper_beta  = 8
dave_upper_gamma = 8
dave_lower_gamma = 0.25


num_to_sample = 8000
lhs_sample    = randomLHS(num_to_sample, 2) 

prev_inci_data_frame = NULL
for(samp_idx in 1:num_to_sample){
  print(samp_idx)
  # Get dave beta and gamma:
  lhs_param_1 = lhs_sample[samp_idx, 1]
  lhs_param_2 = lhs_sample[samp_idx, 2]
  dave_beta   = dave_lower_beta + (dave_upper_beta - dave_lower_beta)*lhs_param_1
  dave_gamma  = dave_lower_gamma + (dave_upper_gamma - dave_lower_gamma)*lhs_param_2
  
  if( (dave_beta/dave_gamma <1)){
    next
  }
  
  # Assign these to literature parameters:
  lit_lambda = dave_beta
  lit_beta   = dave_gamma
  lit_rho    = lit_lambda / lit_beta
  
  ## run an SIR model
  df1         = sir(beta = dave_beta, gamma = dave_gamma, S0 = s0, I0 = 1-s0, R0 = 0, times = 0:50)
  peak_height = 1 - (1 + log(lit_rho*s0))/lit_rho
  
  temp_integrand      = function(x){time_integrand(x, s0, lit_rho, lit_beta)}
  tao_peak            = log(lit_rho*s0) / lit_rho
  
  peak_time           = integrate(temp_integrand, 0, tao_peak)$value
  
  # Get the incidence time + height using convex optimization problem.
  minimize_function        = function(tao){ (-1 + tao + 2*s0*exp(-lit_rho*tao) - 1 / lit_rho)**2 }
  tao_before_max_incidence = optim(0.5, minimize_function, method = "Brent", lower = 0, upper = 1)$par
  prevalence_at_this_time  = 1 - tao_before_max_incidence - s0*exp(-lit_rho*tao_before_max_incidence)
  
  # So, I think the tao above is the tao directly before the tao of max incidence.
  # Can I use the differential equation to move it a single step forwards?
  temp_incidence    = df1$incidence
  temp_incidence[1] = 0
  
  max_incidence_time   = (which.max(temp_incidence)-1)
  max_incidence_value  = max(temp_incidence)
  
  # Collect the relevant info:
  prev_gt_inci_peak = 0
  prev_gt_inci_peak = ifelse(prevalence_at_this_time>max_incidence_value, "prevalence", prev_gt_inci_peak)
  prev_gt_inci_peak = ifelse(prevalence_at_this_time<max_incidence_value, "incidence", prev_gt_inci_peak)
  prev_gt_inci_peak = ifelse(prevalence_at_this_time==max_incidence_value, "equal", prev_gt_inci_peak)
  
  prev_gt_inci_time = 0
  prev_gt_inci_time = ifelse(peak_time>max_incidence_time, "prevalence", prev_gt_inci_time)
  prev_gt_inci_time = ifelse(peak_time<max_incidence_time, "incidence", prev_gt_inci_time)
  prev_gt_inci_time = ifelse(peak_time==max_incidence_time, "equal", prev_gt_inci_time)
  
  temp_row = data.frame(beta = dave_beta, gamma = dave_gamma, 
                        prev_peak = prevalence_at_this_time, prev_peak_time = peak_time,
                        inci_peak = max_incidence_value, inci_peak_time = max_incidence_time,
                        prev_gt_inci_peak = (prev_gt_inci_peak),
                        prev_gt_inci_times = (prev_gt_inci_time) )
  prev_inci_data_frame = rbind(prev_inci_data_frame, temp_row)
}

g1 = ggplot(prev_inci_data_frame, aes(x = beta, y = gamma, color = prev_gt_inci_peak)) + geom_point() + ggtitle("Which is Greater: Incidence Peak or Prevalence Peak?")
g2 = ggplot(prev_inci_data_frame, aes(x = beta, y = gamma, color = prev_gt_inci_times)) + geom_point() + ggtitle("Which is Greater: Incidence Time or Prevalence Time?")
list_gg = list(g1,g2)
grid.arrange(list_gg[[1]], list_gg[[2]], ncol = 1)
