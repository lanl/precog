###############################
###############################
### Obtaining MOA Forecasts ###
###############################
###############################

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
NUM_RUN = 5
for(i in 1:NUM_RUN){
  cmdLines <- paste0("Rscript --vanilla ",my_path,"/evaluate_moa_internal.R",
                     " --total_run=",NUM_RUN,
                     " --cur_run=",i,
                     " --scale_ts=yes")
  runSlurm <- slurmarray(cmdLines,sname="PP",stime="540",smem="5G",
                         soutdir=paste0('./logs/'), sparallel = 0)
  system(runSlurm)
}


### NEIGHBORHOOD FREQUENCIES FOR COVID
NUM_RUN = 2
for(i in 1:NUM_RUN){
  cmdLines <- paste0("Rscript --vanilla ",my_path,"/evaluate_moa_internal_covidfrequencies.R",
                     " --total_run=",NUM_RUN,
                     " --cur_run=",i,
                     " --scale_ts=yes")
  runSlurm <- slurmarray(cmdLines,sname="PP",stime="540",smem="5G",
                         soutdir=paste0('./logs/'), sparallel = 0)
  system(runSlurm)
}


### LIBRARY FREQUENCIES FOR COVID
NUM_RUN = 2
for(i in 1:NUM_RUN){
  cmdLines <- paste0("Rscript --vanilla ",my_path,"/evaluate_moa_internal_marginalfrequencies.R",
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
  cmdLines <- paste0("Rscript --vanilla ",my_path,"/evaluate_moa_internal_modespecific.R",
                     " --total_run=",1,
                     " --cur_run=",1,
                     " --include_mode=",MODES[i],
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
  cmdLines <- paste0("Rscript --vanilla ",my_path,"/evaluate_moa_internal_diseasespecific.R",
                     " --total_run=",1,
                     " --cur_run=",1,
                     " --include_disease=",DISEASES[i],
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
  cmdLines <- paste0("Rscript --vanilla ",my_path,"/evaluate_moa_internal_sourcespecific.R",
                     " --total_run=",1,
                     " --cur_run=",1,
                     " --include_disease_source=",DISEASE_SOURCES[i],
                     " --scale_ts=yes")
  runSlurm <- slurmarray(cmdLines,sname="PP",stime="540",smem="5G",
                         soutdir=paste0('./logs/'), sparallel = 0)
  system(runSlurm)
}


