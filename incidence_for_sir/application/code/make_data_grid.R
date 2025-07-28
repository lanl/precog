
# Testing script for all the interpolations I'm doing.
## Author: Murph
library(nleqslv)
library(ggplot2)
library(deSolve)
library(pracma)
library(reshape)
library(parallel)
library(doParallel)
library(this.path)
setwd(this.path::here())
source(paste0(this.path::here(), "/../../simulation_study/incidence_to_sir_helpers.R"))

# Starting alpha,beta for optim calls:
alpha_initial                = 1.2588322 / 2
beta_initial                 = 1.0051870 / 2

slurm_idx = as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID"))
print(paste("Running slurm job:", slurm_idx))
if(is.na(slurm_idx)){
	slurm_idx = 2
}

# This is the code to make the initial grid:
inputgrid <- expand.grid(S0  = 0.9,
                         I0  = seq(0.0001,.0025,length.out = 100),
                         PIV = seq(0.001,.051,length.out = 100),
                         PIT = seq(1,35,length.out = 100))
inputgrid <- subset(inputgrid, select=c("S0","I0","PIV","PIT"))

total_num_of_jobs = 8
rows_for_this_job = nrow(inputgrid) / total_num_of_jobs
lower_idx = 1 + rows_for_this_job*(slurm_idx-1)
upper_idx = lower_idx + (rows_for_this_job-1)

parfctn = function(i){
  library(nleqslv)
  library(ggplot2)
  library(deSolve)
  library(pracma)
  library(reshape)
  
  print(i)
  if(inputgrid$PIV[i] > inputgrid$I0[i]){
    # This is where murph's methods will update what dave did.
    observed_incidence_peak_time = inputgrid$PIT[i]
    observed_peak_incidence      = inputgrid$PIV[i]
    s0                           = 0.9
    i0                           = inputgrid$I0[i]

    # temp_fn          = function(x){(fn1(x[1], x[2], observed_incidence_peak_time, observed_peak_incidence, s0, i0) +
    #                                   fn2(x[1], x[2], observed_incidence_peak_time, observed_peak_incidence, s0, i0)) }
    xx               = ode_approx(alpha_initial, beta_initial, observed_incidence_peak_time, observed_peak_incidence, s0, i0)
    # xx                    = fminunc(c(alpha_initial, beta_initial), temp_fn)
    calculated_alpha = xx$alpha
    calculated_beta  = xx$beta

    ### Now let's graph the SIR curve from these parameters and compare it to the initial time and height of the
    ### incidence peak.
    dave_beta  = calculated_alpha
    dave_gamma = calculated_beta

  } else {
    dave_beta  = 0
    dave_gamma = 0
  }
  return(list(dave_beta=dave_beta, dave_gamma = dave_gamma))
}

sockettype <- "PSOCK"
ncores <- 99
cl <- parallel::makeCluster(spec = ncores,type = sockettype) #, outfile=""
setDefaultCluster(cl)
registerDoParallel(cl)
all_returns <- foreach(i=lower_idx:upper_idx,
                  .verbose = T)%dopar%{
                    print(i)
                    parfctn(i)
                  }
stopCluster(cl)


# all_returns = mclapply(lower_idx:upper_idx, parftcn, mc.cores = 98)

full_data   = NULL
dave_betas  = c()
dave_gammas = c()
for(idx in 1:length(all_returns)){
  dave_betas  = c(dave_betas, all_returns[[idx]]$dave_beta)
  dave_gammas = c(dave_gammas, all_returns[[idx]]$dave_gamma)
}
inputgrid$dave_gamma = NA
inputgrid$dave_betas = NA
inputgrid[lower_idx:upper_idx,]$dave_gamma  = dave_gammas
inputgrid[lower_idx:upper_idx,]$dave_betas = dave_betas
write.csv(inputgrid[lower_idx:upper_idx,], file = paste("grid_logs/grid_data", slurm_idx, '.csv', sep = ""))


