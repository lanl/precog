# General Script to Perform Bayesian Optimization (as Murph understands it).
## Author: Alexander C. Murph
## Date: August 2023
library(hetGP)
library(lhs)
library(parallel)
library(dplyr)
library(foreach)
library(doParallel)
library(plyr)
library(LaplacesDemon)
library(spatstat)
library(BASS)
library(utils)
library(collapse)

parallel_bayesian_optimization                                    <- function(optim_function,
                                                                               parameter_bounds,
                                                                               initial_lhs,
                                                                               number_of_iters,
                                                                               is_discrete,
                                                                               num_cores){
  # Get top-level information:
  param_dim                                                       <- length(parameter_bounds)
  
  # Create initial LHS
  initial_sample                                                  <- randomLHS(initial_lhs, param_dim)
  for(var_idx in 1:length(parameter_bounds)){
    var_name                                                      <- names(parameter_bounds)[var_idx]
    var_bounds                                                    <- parameter_bounds[[var_name]]
    initial_sample[,var_idx]                                      <- initial_sample[,var_idx] * (var_bounds[2] - var_bounds[1]) + var_bounds[1]
    if(as.logical(is_discrete[var_idx])) initial_sample[,var_idx] <- ceiling(initial_sample[,var_idx])
  }
  initial_sample                                                  <- data.frame(initial_sample)
  names(initial_sample)                                           <- names(parameter_bounds)
  
  # Go through the optimization function (in parallel) to get the values for every
  # input in the LHS.
  if(!file.exists('initial_lhs.RData')){
  	sockettype                                                <- "PSOCK"
  	cl                                                        <- parallel::makeCluster(spec = num_cores,type = sockettype,outfile="")
  	setDefaultCluster(cl)
  	registerDoParallel(cl)
  	sim_ts                                                    <- foreach(i=1:initial_lhs,
									       .errorhandling='pass',
  	                                                                          .verbose = T, .packages="foreach")%dopar%{
  	                                                                              optim_function(unlist(initial_sample[i,]))
  	                                                                          }
  	outputs                                                   <- unlist(sim_ts)
  	initial_sample$outputs                                    <- outputs
  	stopCluster(cl) 


  	save(initial_sample, file="initial_lhs.RData")
  } else {
	load(file = 'initial_lhs.RData')										       
  }
  # Fit the initial GP.
  X                                                               <- list(X0 = as.matrix(initial_sample[,1:length(parameter_bounds)]),
                                                                           Z0 = initial_sample[,'outputs'], 
                                                                           mult = rep(1,times = nrow(initial_sample)))
  
  model                                                           <- mleHomGP(X = X, Z = initial_sample[,'outputs'])
  
  # Now we need to go through iterations to determine new points to check:
  # The following criterion is for minimization.
  crit                                                            <- "crit_EI"
  print(paste("Starting to do a sequential design."))
  for(epoch in 1:number_of_iters){
    print(paste("Starting iteration", epoch, "in the sequential design."))
    res                                                           <- crit_optim(model, crit = crit,
                      h=4, control = list(multi.start = 50, maxit = 30),
                      ncores = num_cores)
    print("The current best parameters are:")
    print(res)
    print(res$par)
    print("-------------------------")
    
    newX                                                          <- as.vector(res$par)


    
    newX[is_discrete] = ceiling(newX[is_discrete])
    
    newZ                                                          <- optim_function(newX)
    print("newX is:")
    print(newX)
    print(as.matrix(newX))
    print(t(as.matrix(newX)))
    print("newZ is:")
    print(newZ)
    model                                                         <- update(object = model, Xnew = t(as.matrix(newX)), Znew = newZ$Score)
  }
  res                                                             <- crit_optim(model, crit = crit,
                                                                                  h=4, control = list(multi.start = 50, maxit = 30),
                                                                                  ncores = num_cores)
  optimal_parameters                                              <- as.vector(res$par)
  optimal_parameters[is_discrete]                                 <- ceiling(optimal_parameters[is_discrete])
  
  print("The optimal parameters are:")
  print(res)
  print(res$par)
  return(res$par)
}

# # Example:
# optim_function = function(par0, par2, par3){
#   Sys.sleep(1)
#   return(runif(1))
# }
# parameter_bounds = list(par0 = c(0,1),
#                         par1 = c(0,1),
#                         par2 = c(0,1))
# is_discrete = c(0, 0, 0)
# 
# parallel_bayesian_optimization(optim_function,
#                                            parameter_bounds,
#                                            9,
#                                            9,
#                                            is_discrete,
#                                            4)
  
  

