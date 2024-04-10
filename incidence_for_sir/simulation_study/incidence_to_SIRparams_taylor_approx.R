## Testing the process of going from incidence peak+time to SIR parameters.
#### Author: Alexander C. Murph
library(nleqslv)
library(ggplot2)
library(deSolve)
library(pracma)
library(reshape)
theme_set(theme_bw())
setwd("~/GitLab/incidence_for_sir/simulation_study")
source('incidence_to_sir_helpers.R')

# Chose the initial parameters based on dave's paper:
peak_height = runif(1, min = 0.005, max = 0.035)
peak_time   = sample(5:30, size = 1)
# Note that i0 can't be greater than the above peak height...
N          = 10000
s0         = 0.95
i0         = 0.003
r0         = 1 - i0 - s0

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
# slow down the recovery rate:
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

observed_peak_incidence      = max_incidence_value
observed_incidence_peak_time = max_incidence_time
alpha_initial                = 1.1788322 / 2
beta_initial                 = 0.8851870 / 2

# ggplot(data = df1, aes(x=time))+
#   geom_line(aes(y=I))+  # prevalence
#   geom_line(aes(y=incidence), color=I("red")) + # incidence
#   geom_point(aes(x = observed_incidence_peak_time, y = observed_peak_incidence))
####################################################################################
####################################################################################


####################################################################################
####################################################################################
## Next, perform the optimization problem of two variables.

temp_fn          = function(x){(fn1_taylor(x[1], x[2], observed_incidence_peak_time, observed_peak_incidence, s0, i0) + 
                                 fn2_taylor(x[1], x[2], observed_incidence_peak_time, observed_peak_incidence, s0, i0)) }
# f0_abv           = function(x){fn0(x[1], x[2], x[3], observed_incidence_peak_time, observed_peak_incidence)}
fminunc(c(alpha_initial, beta_initial), temp_fn)
print(c(lit_lambda, lit_beta))
xx               = fminunc(c(alpha_initial, beta_initial), temp_fn)
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
max_incidence_time_df2   = which.max(temp_incidences) - 1
max_incidence_value_df2  = max(temp_incidences)
####################################################################################
####################################################################################

## plot prevalence and incidence for the sir model
ggplot(data = df2, aes(x=time))+
        geom_line(aes(y=I))+  # prevalence
        geom_line(aes(y=incidence), color=I("red")) + # incidence
        geom_point(aes(x = observed_incidence_peak_time, y = observed_peak_incidence)) +
  xlim(0,35) + ylim(0,0.05)
