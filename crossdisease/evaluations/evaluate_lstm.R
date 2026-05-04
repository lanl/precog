################################
################################
### Obtaining LSTM Forecasts ###
################################
################################

my_path = '~/GitLab/frodo/evaluations/'
setwd(my_path)
source(paste0(my_path,'/../SLURMarray.r'))

lut_path = '~/GitLab/frodo/lookup_tables/'
DISEASES = readxl::read_xlsx(paste0(lut_path, 'Disease_ML.xlsx'))
DISEASES = data.frame(DISEASES)
DISEASES = DISEASES$FILES

lut_path = '~/GitLab/frodo/lookup_tables/'
DISEASE_SOURCES = readxl::read_xlsx(paste0(lut_path, 'Evaluations_ML.xlsx'))
DISEASE_SOURCES = data.frame(DISEASE_SOURCES)
DISEASE_SOURCES = DISEASE_SOURCES$FILES
DISEASE_SOURCES = gsub('.RDS','',DISEASE_SOURCES)




#####################
### All Available ###
#####################

### SCALED
if(!dir.exists(paste0(my_path,"lstm/all"))){dir.create(paste0(my_path,"lstm/all"))}
NUM_RUN = 1
for(i in 1:NUM_RUN){
  cmdLines <- paste0("Rscript --vanilla ",my_path,"/evaluate_lstm_internal.R",
                     " --alpha=0.5",
                     " --total_run=",NUM_RUN,
                     " --cur_run=",i,
                     " --scale_ts=yes")
  runSlurm <- slurmarray(cmdLines,sname="PP",stime="540",smem="5G",
                         soutdir=paste0('./logs/'), sparallel = 0)
  system(runSlurm)
}


### SCALED, OTHER ALPHAS
NUM_RUN = 3
for(i in 1:NUM_RUN){
    cmdLines <- paste0("Rscript --vanilla ",my_path,"/evaluate_lstm_internal.R",
                       " --alpha=",-1,
                       " --total_run=",NUM_RUN,
                       " --cur_run=",i,
                       " --scale_ts=yes")
    runSlurm <- slurmarray(cmdLines,sname="PP",stime="540",smem="5G",
                           soutdir=paste0('./logs/'), sparallel = 0)
    system(runSlurm)
}









##############################
### Mode-Specific Modeling ###
##############################

MODES = c('respiratory','fecaloral','sexual','vectorborne')

### SCALED
for(i in 1:length(MODES)){
  cmdLines <- paste0("Rscript --vanilla ",my_path,"/evaluate_lstm_internal_modespecific.R",
                     " --alpha=0.5",
                     " --total_run=",1,
                     " --cur_run=",1,
                     " --include_mode=",MODES[i],
                     " --warm_start=no",
                     " --scale_ts=yes")
  runSlurm <- slurmarray(cmdLines,sname="PP",stime="540",smem="5G",
                         soutdir=paste0('./logs/'), sparallel = 0)
  system(runSlurm)
}



### SCALED, OTHER ALPHAS
for(i in 1:length(MODES)){
    cmdLines <- paste0("Rscript --vanilla ",my_path,"/evaluate_lstm_internal_modespecific.R",
                       " --alpha=",-1,
                       " --total_run=",1,
                       " --cur_run=",1,
                       " --include_mode=",MODES[i],
                       " --warm_start=no",
                       " --scale_ts=yes")
    runSlurm <- slurmarray(cmdLines,sname="PP",stime="540",smem="5G",
                           soutdir=paste0('./logs/'), sparallel = 0)
    system(runSlurm)
}



#################################
### Disease-Specific Modeling ###
#################################

### SCALED
for(i in 1:length(DISEASES)){
  cmdLines <- paste0("Rscript --vanilla ",my_path,"/evaluate_lstm_internal_diseasespecific.R",
                     " --alpha=0.5",
                     " --total_run=",1,
                     " --cur_run=",1,
                     " --include_disease=",DISEASES[i],
                     " --warm_start=no",
                     " --scale_ts=yes")
  runSlurm <- slurmarray(cmdLines,sname="PP",stime="540",smem="5G",
                         soutdir=paste0('./logs/'), sparallel = 0)
  system(runSlurm)
}




### SCALED, OTHER ALPHAS
for(i in 1:length(DISEASES)){
  cmdLines <- paste0("Rscript --vanilla ",my_path,"/evaluate_lstm_internal_diseasespecific.R",
                     " --alpha=",-1,
                     " --total_run=",1,
                     " --cur_run=",1,
                     " --include_disease=",DISEASES[i],
                     " --warm_start=no",
                     " --scale_ts=yes")
  runSlurm <- slurmarray(cmdLines,sname="PP",stime="540",smem="5G",
                         soutdir=paste0('./logs/'), sparallel = 0)
  system(runSlurm)
}







##########################################
### Disease x Source Specific Modeling ###
##########################################

### SCALED
for(i in 1:length(DISEASE_SOURCES)){
  cmdLines <- paste0("Rscript --vanilla ",my_path,"/evaluate_lstm_internal_sourcespecific.R",
                     " --alpha=0.5",
                     " --total_run=",1,
                     " --cur_run=",1,
                     " --include_disease_source=",DISEASE_SOURCES[i],
                     " --warm_start=no",
                     " --scale_ts=yes")
  runSlurm <- slurmarray(cmdLines,sname="PP",stime="540",smem="5G",
                         soutdir=paste0('./logs/'), sparallel = 0)
  system(runSlurm)
}


for(i in 1:length(DISEASE_SOURCES)){
  cmdLines <- paste0("Rscript --vanilla ",my_path,"/evaluate_lstm_internal_sourcespecific.R",
                     " --alpha=",-1,
                     " --total_run=",1,
                     " --cur_run=",1,
                     " --include_disease_source=",DISEASE_SOURCES[i],
                     " --warm_start=no",
                     " --scale_ts=yes")
  runSlurm <- slurmarray(cmdLines,sname="PP",stime="540",smem="5G",
                         soutdir=paste0('./logs/'), sparallel = 0)
  system(runSlurm)
}




