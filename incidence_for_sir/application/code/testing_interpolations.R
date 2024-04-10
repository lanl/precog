# Testing script for all the interpolations I'm doing.
## Author: Murph
## Date: April 2024
setwd("~/GitLab/incidence_for_sir/data_application/code")
source("../../incidence_to_sir_helpers.R")
library(nleqslv)
library(ggplot2)
library(deSolve)
library(pracma)
library(reshape)
library(parallel)

set.seed(15)
# Starting alpha,beta for optim calls:
alpha_initial                = 1.1788322 / 2
beta_initial                 = 0.8851870 / 2

# Compile Grid from the Darwin files.
inputgrid = NULL
for(file_idx in 1:100){
  data_name   = paste("grid_logs/grid_data", file_idx, ".csv", sep = "")
  temp_rows   = read.csv(data_name)
  temp_rows$X = NULL
  inputgrid   = rbind(inputgrid, temp_rows)
}


# Inputs from me and from the algorithm.
i0_values   = inputgrid$I0
piv_values  = inputgrid$PIV
pit_values  = inputgrid$PIT
dave_betas  = inputgrid$dave_betas
dave_gammas = inputgrid$dave_gamma

theta_i0 = 2.161298e-04
z1 = 0.0122876146
z2 = 21.38724

seq(0.00001,0.05001,length.out = 100)

# This stuff needs to all go into JAGS:
i0_scaling_factor  <- (99/(0.05001-0.00001))
piv_scaling_factor <- (99/(0.1001-0.00001))
pit_scaling_factor <- (99/(36-1))

i0_integer_lower  = floor( (theta_i0-0.00001) * i0_scaling_factor) + 1 # We add 1 here b/c we start counting at 1, not zero.
x0                = (i0_integer_lower-1) / i0_scaling_factor + 0.00001
i0_integer_upper  = i0_integer_lower + 1
x1                = (i0_integer_upper-1) / i0_scaling_factor + 0.00001

piv_integer_lower = floor((z1-0.00001) * piv_scaling_factor) + 1
y0                = (piv_integer_lower-1) / piv_scaling_factor + 0.00001
piv_integer_upper = piv_integer_lower + 1
y1                = (piv_integer_upper-1) / piv_scaling_factor + 0.00001

pit_integer_lower = floor( (z2-1.0) * pit_scaling_factor ) + 1
z0                = (pit_integer_lower-1) / pit_scaling_factor + 1.0
pit_integer_upper = pit_integer_lower + 1
ze1               = (pit_integer_upper-1) / pit_scaling_factor + 1.0

# The order of indices is (i0)(PIV)(PIT)
c000_idx = (i0_integer_lower) + (piv_integer_lower-1)*100 + 10000 * (pit_integer_lower-1)
c100_idx = (i0_integer_lower + 1) + (piv_integer_lower-1)*100 + 10000 * (pit_integer_lower-1)
c010_idx = (i0_integer_lower) + (piv_integer_lower-1 + 1)*100 + 10000 * (pit_integer_lower-1)
c001_idx = (i0_integer_lower) + (piv_integer_lower-1)*100 + 10000 * (pit_integer_lower+1-1)
c110_idx = (i0_integer_lower + 1) + (piv_integer_lower-1 + 1)*100 + 10000 * (pit_integer_lower-1)
c011_idx = (i0_integer_lower) + (piv_integer_lower-1 + 1)*100 + 10000 * (pit_integer_lower+1-1)
c101_idx = (i0_integer_lower + 1) + (piv_integer_lower-1)*100 + 10000 * (pit_integer_lower + 1-1)
c111_idx = (i0_integer_lower + 1) + (piv_integer_lower-1 + 1)*100 + 10000 * (pit_integer_lower + 1-1)
all_indices = c(c000_idx, c100_idx, c010_idx, c001_idx, c110_idx, c011_idx, c101_idx, c111_idx)

# Get the trilinear interpolation for dave's beta:
c000  = dave_betas[c000_idx]
c100  = dave_betas[c100_idx]
c010  = dave_betas[c010_idx]
c001  = dave_betas[c001_idx]
c110  = dave_betas[c110_idx]
c011  = dave_betas[c011_idx]
c101  = dave_betas[c101_idx]
c111  = dave_betas[c111_idx]
all_beta_values = c(c000,c100,c010, c001, c110, c011, c101, c111)

a0 = (-c000*x1*y1*ze1+c001*x1*y1*z0+c010*x1*y0*ze1-c011*x1*y0*z0)/((x0-x1)*(y0-y1)*(z0-ze1)) + 
      (c100*x0*y1*ze1-c101*x0*y1*z0-c110*x0*y0*ze1+c111*x0*y0*z0)/((x0-x1)*(y0-y1)*(z0-ze1))

a1 = (c000*y1*ze1-c001*y1*z0-c010*y0*ze1+c011*y0*z0)/((x0-x1)*(y0-y1)*(z0-ze1)) +
      (-c100*y1*ze1+c101*y1*z0+c110*y0*ze1-c111*y0*z0)/((x0-x1)*(y0-y1)*(z0-ze1))

a2 = (c000*x1*ze1-c001*x1*z0-c010*x1*ze1+c011*x1*z0)/((x0-x1)*(y0-y1)*(z0-ze1)) +
      (-c100*x0*ze1+c101*x0*z0+c110*x0*ze1-c111*x0*z0)/((x0-x1)*(y0-y1)*(z0-ze1))

a3 = (c000*x1*y1-c001*x1*y1-c010*x1*y0+c011*x1*y0)/((x0-x1)*(y0-y1)*(z0-ze1)) + 
      (-c100*x0*y1+c101*x0*y1+c110*x0*y0-c111*x0*y0)/((x0-x1)*(y0-y1)*(z0-ze1))

a4 = (-c000*ze1+c001*z0+c010*ze1-c011*z0+c100*ze1-c101*z0-c110*ze1+c111*z0)/((x0-x1)*(y0-y1)*(z0-ze1))

a5 = (-c000*y1+c001*y1+c010*y0-c011*y0+c100*y1-c101*y1-c110*y0+c111*y0)/((x0-x1)*(y0-y1)*(z0-ze1))

a6 = (-c000*x1+c001*x1+c010*x1-c011*x1+c100*x0-c101*x0-c110*x0+c111*x0)/((x0-x1)*(y0-y1)*(z0-ze1))

a7 = (c000-c001-c010+c011-c100+c101+c110-c111)/((x0-x1)*(y0-y1)*(z0-ze1))

beta = a0 + a1*theta_i0 + a2*z1 + a3*z2 + a4*theta_i0*z1 + a5*theta_i0*z2 + a6*z1*z2 + a7*theta_i0*z1*z2

# Get the trilinear interpolation for dave's gamma:
c000 = dave_gammas[c000_idx]
c100 = dave_gammas[c100_idx]
c010 = dave_gammas[c010_idx]
c001 = dave_gammas[c001_idx]
c110 = dave_gammas[c110_idx]
c011 = dave_gammas[c011_idx]
c101 = dave_gammas[c101_idx]
c111 = dave_gammas[c111_idx]

a0 = (-c000*x1*y1*ze1+c001*x1*y1*z0+c010*x1*y0*ze1-c011*x1*y0*z0)/((x0-x1)*(y0-y1)*(z0-ze1)) + 
  (c100*x0*y1*ze1-c101*x0*y1*z0-c110*x0*y0*ze1+c111*x0*y0*z0)/((x0-x1)*(y0-y1)*(z0-ze1))
a1 = (c000*y1*ze1-c001*y1*z0-c010*y0*ze1+c011*y0*z0)/((x0-x1)*(y0-y1)*(z0-ze1)) +
  (-c100*y1*ze1+c101*y1*z0+c110*y0*ze1-c111*y0*z0)/((x0-x1)*(y0-y1)*(z0-ze1))
a2 = (c000*x1*ze1-c001*x1*z0-c010*x1*ze1+c011*x1*z0)/((x0-x1)*(y0-y1)*(z0-ze1)) +
  (-c100*x0*ze1+c101*x0*z0+c110*x0*ze1-c111*x0*z0)/((x0-x1)*(y0-y1)*(z0-ze1))
a3 = (c000*x1*y1-c001*x1*y1-c010*x1*y0+c011*x1*y0)/((x0-x1)*(y0-y1)*(z0-ze1)) + 
  (-c100*x0*y1+c101*x0*y1+c110*x0*y0-c111*x0*y0)/((x0-x1)*(y0-y1)*(z0-ze1))
a4 = (-c000*ze1+c001*z0+c010*ze1-c011*z0+c100*ze1-c101*z0-c110*ze1+c111*z0)/((x0-x1)*(y0-y1)*(z0-ze1))
a5 = (-c000*y1+c001*y1+c010*y0-c011*y0+c100*y1-c101*y1-c110*y0+c111*y0)/((x0-x1)*(y0-y1)*(z0-ze1))
a6 = (-c000*x1+c001*x1+c010*x1-c011*x1+c100*x0-c101*x0-c110*x0+c111*x0)/((x0-x1)*(y0-y1)*(z0-ze1))
a7 = (c000-c001-c010+c011-c100+c101+c110-c111)/((x0-x1)*(y0-y1)*(z0-ze1))

gamma = a0 + a1*theta_i0 + a2*z1 + a3*z2 + a4*theta_i0*z1 + a5*theta_i0*z2 + a6*z1*z2 + a7*theta_i0*z1*z2


print(c(beta,gamma))
print(inputgrid[all_indices,])

