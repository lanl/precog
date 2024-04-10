## Testing the process of going from prevalence peak+time to SIR parameters.
#### Author: Alexander C. Murph
library(ggplot2)
library(deSolve)
library(reshape2)
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

### Murph update 2023:
# I want to add in the peak time and height to this simulation to see how well
# it works.
# peak_height = 0.065
# peak_time   = 17
peak_height = runif(1, min = 0.005, max = 0.04)
peak_time   = sample(5:30, size = 1)

# Note that i0 can't be greater than the above peak height...
N          = 10000
s0         = .99
i0         = peak_height / 200
r0         = 1 - i0 - s0

# The following will give the value of rho directly:
minimize_function     = function(x){ ((s0+i0) - peak_height - 1/x - log(s0*x)/x )**2 }
rho_of_max_prevalence = optim(1.2, minimize_function, method = "Brent", lower = 1, upper = 100)$par

# My previous method for doing this uses tao rather than beta*tao.  
# This requires me to reparameterize the time integral:
reparmeterized_time_integrand_sans_alpha = function(tao, s0, i0, rho){
  return( 1/( (i0 + s0) - tao/rho - s0*exp(-tao)) )
}

tao_peak_times_alpha       = log(rho_of_max_prevalence*s0)
temp_integrand             = function(x){reparmeterized_time_integrand_sans_alpha(x, s0, i0, rho_of_max_prevalence)}
peak_time_times_alpha      = integrate(temp_integrand, 0, tao_peak_times_alpha)$value
lit_beta                   = peak_time_times_alpha / (peak_time*rho_of_max_prevalence)
lit_lambda                 = lit_beta * rho_of_max_prevalence
lit_alpha                  = lit_lambda
dave_beta                  = lit_lambda
dave_gamma                 = lit_beta

## run an SIR model
df1         = sir(beta = dave_beta, gamma = dave_gamma, S0 = s0, I0 = i0, R0 = r0, times = 0:50)

## plot prevalence and incidence for the sir model
ggplot(data = df1, aes(x=time))+
  geom_line(aes(y=I))+  # prevalence
  geom_line(aes(y=incidence), color=I("red")) +
  geom_point(aes(x=peak_time, y=peak_height), colour="blue")

melted_df = melt(df1[,1:5],id=1)
ggplot(data = melted_df, aes(x=time,y = value, color = variable)) + geom_line() +
  geom_hline(yintercept = i0, color = 'purple') +
  geom_hline(yintercept = s0, color = 'orange') +
  geom_hline(yintercept = r0, color = 'blue') +
  # geom_point(aes(x = max_incidence_time, y = max_incidence_value), color = 'purple', size = 3)+
  geom_point(aes(x = peak_time, y = peak_height), color = 'green', size = 3)


