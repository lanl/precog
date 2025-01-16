## Dave Osthus
## 2-1-24
## Make synthetic training data

## load libraries
library(ggplot2)
library(plyr)
library(data.table)
library(LearnBayes)
library(LaplacesDemon)
library(parallel)
library(doParallel)

theme_set(theme_bw())


library(this.path)
my_path = this.path::here()
setwd(my_path)

## set path
savepath <- "./output/"

## set proxies
Sys.setenv('https_proxy'='http://proxyout.lanl.gov:8080')## not sure why needed, but see here: https://github.com/curl/curl/issues/1015

## define number of cores
ncores <- floor(.5*detectCores())

## define socket type
sockettype <- "PSOCK"

## functions for making synthetic data
source(paste0('./code/',"synthetic_functions.R"))


types_of_curves <- c("sir_rollercoaster", "sir_rollercoaster_wiggle", 'seasonal')
N <- 20000*length(types_of_curves)

############################################
### Generate Synthetic List for Training ###
############################################

cl <- parallel::makeCluster(spec = ncores,type = sockettype)
setDefaultCluster(cl)
registerDoParallel(cl)
print(Sys.time())
sim_ts <- foreach(i=1:N, #added extra 300 to compensate for extra seasonality sims
                      .errorhandling = "pass",
                      .verbose = F)%dopar%{
                        curve_type <- types_of_curves[rep(c(1:length(types_of_curves)), each = 20000)][i]
                        templist = gen_curve(curve_type)
                        templist
                      }
stopCluster(cl)
print(Sys.time())
A = lapply(sim_ts,function(x){ length(x)})
unique(A)
sim_ts = sim_ts[A > 2]
saveRDS(sim_ts, file=paste0(savepath,"synthetic.RDS"))

# A = lapply(sim_ts,function(x){ sum(is.nan(x$ts))})
# unique(A)
# 
# A = lapply(sim_ts,function(x){ sum(is.na(x$ts))})
# unique(A)

# A = lapply(sim_ts,function(x){ length(x)})
# unique(A)
# sim_ts[[which(A==2)[1]]]
# 
# ts = sir(beta, gamma, S0, I0, R0, times)$I
# runif(0,0,quantile(c(rep(0,starttime),ts[1:pmin(max(which(ts/max(ts) > .0001)), length(ts))]),prob=.01))
#######################################
### Generate Synthetic List for MOA ###
#######################################

cl <- parallel::makeCluster(spec = ncores,type = sockettype)
setDefaultCluster(cl)
registerDoParallel(cl)
print(Sys.time())
sim_ts <- foreach(i=1:N, #added extra 300 to compensate for extra seasonality sims
                  .errorhandling = "pass",
                  .verbose = F)%dopar%{
                    curve_type <- types_of_curves[rep(c(1:length(types_of_curves)), each = 20000)][i]
                    templist = gen_curve(curve_type)
                    templist
                  }
stopCluster(cl)
print(Sys.time())
A = lapply(sim_ts,function(x){ length(x)})
unique(A)
sim_ts = sim_ts[A > 2]
saveRDS(sim_ts, file=paste0(savepath,"synthetic_moa.RDS"))




######################################
### Generate Synthetic List for UQ ###
######################################


cl <- parallel::makeCluster(spec = ncores,type = sockettype)
setDefaultCluster(cl)
registerDoParallel(cl)
print(Sys.time())
sim_ts <- foreach(i=1:N, #added extra 300 to compensate for extra seasonality sims
                  .errorhandling = "pass",
                  .verbose = F)%dopar%{
                    curve_type <- types_of_curves[rep(c(1:length(types_of_curves)), each = 20000)][i]
                    templist = gen_curve(curve_type)
                    templist
                  }
stopCluster(cl)
print(Sys.time())
A = lapply(sim_ts,function(x){ length(x)})
unique(A)
sim_ts = sim_ts[A > 2]
saveRDS(sim_ts, file=paste0(savepath,"synthetic_uq.RDS"))




#####################
### Plot Examples ###
#####################
sim_ts = readRDS(paste0(savepath,"synthetic.RDS"))

TYPES = unlist(lapply(sim_ts,function(ll){return(ll$ts_disease)})) 

TS1 = sim_ts[[which(TYPES == 'sir_rollercoaster')[1]]]$ts
TS2 = sim_ts[[which(TYPES == 'sir_rollercoaster_wiggle')[1]]]$ts
TS3 = sim_ts[[which(TYPES == 'seasonal')[4]]]$ts

plot(TS1)
lines(TS1)


plot(TS2)
lines(TS2)


plot(TS3)
lines(TS3)




# 
# curve_type <- types_of_curves[sample(1:length(types_of_curves),1)]
# templist = gen_curve(curve_type)
# plot(templist$ts)
