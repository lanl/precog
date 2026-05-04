################################
################################
### Evaluate Forecastability ###
################################
################################



my_path = '~/GitLab/frodo/features/'
setwd(my_path)
source(paste0(my_path,'/../SLURMarray.r'))

## define paths
savetrainpath <- my_path

#########################################
### Calculate Forecastability Metrics ###
#########################################
NUM_RUN = 6
for(i in 1:NUM_RUN){
  cmdLines <- paste0("Rscript --vanilla ",my_path,"/evaluate_forecastability_internal.R",
                     " --total_run=",NUM_RUN,
                     " --cur_run=",i)
  runSlurm <- slurmarray(cmdLines,sname="PP",stime="540",smem="5G",
                         soutdir=paste0('./logs/'), sparallel = 0)
  system(runSlurm)
}

