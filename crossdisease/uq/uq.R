###############################################
###############################################
### Obtain/Organize Probabilistic Forecasts ###
###############################################
###############################################

my_path = '~/GitLab/frodo/uq/'
setwd(my_path)
source(paste0(my_path,'../SLURMarray.r'))


############
### SMOA ###
############

DIRS = c('~/GitLab/frodo/evaluations/smoanodiff/all/',
         '~/GitLab/frodo/evaluations/smoanodiff/sourcespecific/',
         '~/GitLab/frodo/evaluations/smoanodiff/diseasespecific/',
         '~/GitLab/frodo/evaluations/smoanodiff/modespecific/')
for(i in 1:length(DIRS)){
  if(!dir.exists(gsub('evaluations','uq',DIRS[i]))){dir.create(gsub('evaluations','uq',DIRS[i]))}
   cmdLines <- paste0("Rscript --vanilla ",my_path,"/uq_internal_smoa.R",
                      " --total_run=",1,
                      " --cur_run=",1, 
                      " --path_to_dir=",DIRS[i])
   runSlurm <- slurmarray(cmdLines,sname="PP",stime="540",smem="5G",
                          soutdir=paste0('./logs/'), sparallel = 0)
   system(runSlurm)
}


#############
### OTHER ###
#############


DIRS = c('~/GitLab/frodo/evaluations/frodo/all/',
         '~/GitLab/frodo/evaluations/frodo/sourcespecific_cold/',
         '~/GitLab/frodo/evaluations/frodo/diseasespecific_cold/',
         '~/GitLab/frodo/evaluations/frodo/modespecific_cold/',
         '~/GitLab/frodo/evaluations/lstm/all/',
         '~/GitLab/frodo/evaluations/lstm/sourcespecific_cold/',
         '~/GitLab/frodo/evaluations/lstm/diseasespecific_cold/',
         '~/GitLab/frodo/evaluations/lstm/modespecific_cold/')

for(i in 1:length(DIRS)){
  if(!dir.exists(gsub('evaluations','uq',DIRS[i]))){dir.create(gsub('evaluations','uq',DIRS[i]))}
  cmdLines <- paste0("Rscript --vanilla ",my_path,"/uq_internal_gbmlstm.R",
                     " --total_run=",1,
                     " --cur_run=",1, 
                     " --path_to_dir=",DIRS[i])
  runSlurm <- slurmarray(cmdLines,sname="PP",stime="540",smem="5G",
                         soutdir=paste0('./logs/'), sparallel = 0)
  system(runSlurm)
}






