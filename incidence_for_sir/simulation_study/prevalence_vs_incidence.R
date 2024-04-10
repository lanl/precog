# Code to compare the prevalence map to the incidence map.
## Author: Alexander C. Murph
## Date: April 2024

## load libraries
library(ggplot2)
library(deSolve)
library(reshape2)
library(ggpubr)
setwd("~/GitLab/incidence_for_sir/simulation_study")
# Comment out the following two lines after the first run of this script.
# source("../prevalence_maps/gen_reasonable_samples.R")
count = 0

theme_set(theme_bw())

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

### Murph update 2023:
# I want to add in the peak time and height to this simulation to see how well
# it works.
N          = 10000
s0         = .9
i0         = 0.05
r0         = 1 - i0 - s0
# dave_beta  = runif(1, min=0.25, max = 2)
# dave_gamma = runif(1, min=0.25, max = dave_beta)
count      = count + 1
dave_beta  = reasonable_samples[count,'beta']
dave_gamma = reasonable_samples[count,'gamma']/2
lit_lambda = dave_beta
lit_alpha  = dave_beta
lit_beta   = dave_gamma
lit_rho    = lit_lambda / lit_beta
lit_gamma  = 1/lit_rho

time_integrand = function(tao, s0, alpha, beta){
  return( 1/(i0 + s0 - beta*tao - s0*exp(-alpha*tao)) )
}

## run an SIR model
df1         = sir(beta = dave_beta, gamma = dave_gamma, S0 = s0, I0 = i0, R0 = r0, times = 0:50)
tao_peak    = log(lit_rho*s0) / lit_alpha


peak_height = i0 + (s0 - lit_gamma) - lit_gamma * log(s0 / lit_gamma)

temp_integrand      = function(x){time_integrand(x, s0, lit_alpha, lit_beta)}
peak_time           = integrate(temp_integrand, 0, tao_peak)$value

# In the following, I should only consider tao such that i0+s0-s0e... is above zero.
taos               = 1:1000/1000
inverse_integrands = 1/(temp_integrand(taos))
# These should stay above zero.  If they don't, we have to upper bound our taos.
indices_below_zero = which( inverse_integrands <= 0 )
if(length(indices_below_zero)==0){
  upper_tao = 1
} else {
  upper_tao = taos[min(indices_below_zero)-1]
}

minimize_function        = function(tao){ (-(i0 + s0) + lit_beta*tao + 2*s0*exp(-lit_alpha*tao) - 1 / lit_rho)**2 }
tao_before_max_incidence = optim(0.5, minimize_function, method = "Brent", lower = 0, upper = upper_tao)$par
prevalence_at_this_time  = (i0 + s0) - lit_beta*tao_before_max_incidence - s0*exp(-lit_alpha*tao_before_max_incidence)

# So, I think the tao above is the tao directly before the tao of max incidence.
# Can I use the differential equation to move it a single step forwards?
tao_of_max_incidence = tao_before_max_incidence + prevalence_at_this_time * lit_beta
max_incidence_time   = integrate(temp_integrand, 0, tao_of_max_incidence)$value
max_incidence_value  = ((i0+s0) - lit_beta*tao_before_max_incidence - s0*exp(-lit_alpha*tao_before_max_incidence))*s0*exp(-lit_alpha*tao_before_max_incidence)*lit_lambda

temp_incidence    = df1$incidence
temp_incidence[1] = 0

## plot prevalence and incidence for the sir model
melted_df = melt(df1[,1:5],id=1)
ggplot(data = melted_df, aes(x=time,y = value, color = variable)) + geom_line() +
  geom_hline(yintercept = i0, color = 'purple') +
  geom_hline(yintercept = s0, color = 'orange') +
  geom_hline(yintercept = r0, color = 'blue') +
  geom_point(aes(x = max_incidence_time, y = max_incidence_value), color = 'purple', size = 3)+
  geom_point(aes(x = peak_time, y = peak_height), color = 'green', size = 3)


# p1 = ggplot(data = melted_df, aes(x=time,y = value, linetype = variable)) + 
#   geom_line()+
#   scale_linetype_manual(values=c('solid', "twodash", "dotted", 'dashed'))+ 
#   guides(linetype=guide_legend(title="SIR Quantity")) + xlim(0,20)
# 
# write.csv(melted_df, file = "sir_curve_example_1.csv")
# # 
# # p1 = p1 + theme(legend.position = 'none')
# # p2 = p2 + theme(legend.position = 'none')
# # 
# ggarrange(p1, p2, ncol=2, nrow=1, common.legend = TRUE, legend="right")


# ggplot(data = df1, aes(x=time))+
#   geom_line(aes(y=I))+  # prevalence
#   geom_line(aes(y=incidence), color=I("red")) + # incidence
#   geom_hline(yintercept = peak_height) +
#   geom_vline(xintercept = peak_time) +
#   geom_vline(xintercept = (max_incidence_time), color = 'red')+
#   geom_hline(yintercept = max_incidence_value, color = 'red') #+
#   # geom_vline(xintercept = (which.max(temp_incidence)-1), color = 'green') +
#   # geom_hline(yintercept = max(temp_incidence), color = 'green')





