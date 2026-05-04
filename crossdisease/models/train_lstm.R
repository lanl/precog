######################################
### Train Fancier Component Models ###
######################################


my_path = "~/GitLab/frodo/models"
source(paste0(my_path,'/../SLURMarray.r'))
lut_path = '~/GitLab/frodo/lookup_tables/'

setwd(my_path)

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
YEAR_RUNS = c(2010:2024)
for(i in 1:length(YEAR_RUNS)){
  cmdLines <- paste0("Rscript --vanilla ",my_path,"/train_lstm_internal_byyear.R",
                     " --year_cutoff=", YEAR_RUNS[i],
                     " --alpha=0.5",
                     " --scale_ts=yes")
  runSlurm <- slurmarray(cmdLines,sname="PP",stime="540",smem="50G",
                         soutdir=paste0('./logs/'), sparallel = 0)
  system(runSlurm)
}


### SCALED, OTHER ALPHAS
ALPHA_GRID = c(0.025, 0.1, 0.25, 0.75, 0.9, 0.975)
YEAR_RUNS = c(2010:2024)
#for(i in 1:length(YEAR_RUNS)){
  for(j in 1:length(ALPHA_GRID)){
    cmdLines <- paste0("Rscript --vanilla ",my_path,"/train_lstm_internal_byyear.R",
                       " --year_cutoff=", -1,#YEAR_RUNS[i],
                       " --alpha=",ALPHA_GRID[j],
                       " --scale_ts=yes")
    runSlurm <- slurmarray(cmdLines,sname="PP",stime="540",smem="50G",
                           soutdir=paste0('./logs/'), sparallel = 0)
    system(runSlurm)
    
  }
#}

##############################
### Mode-Specific Modeling ###
##############################


MODES = c('respiratory','fecaloral','sexual','vectorborne')

### SCALED
for(j in 1:length(MODES)){
  cmdLines <- paste0("Rscript --vanilla ",my_path,"/train_lstm_internal_byyearmodeonly.R",
                     " --year_cutoff=", -1,
                     " --alpha=0.5",
                     " --scale_ts=yes",
                     " --include_mode=",MODES[j])
  runSlurm <- slurmarray(cmdLines,sname="PP",stime="540",smem="50G",
                         soutdir=paste0('./logs/'), sparallel = 0)
  system(runSlurm)
}





### SCALED, OTHER ALPHAS
ALPHA_GRID = c(-1)
for(j in 1:length(MODES)){
  for(k in 1:length(ALPHA_GRID)){
    cmdLines <- paste0("Rscript --vanilla ",my_path,"/train_lstm_internal_byyearmodeonly.R",
                       " --year_cutoff=", -1,
                       " --alpha=",ALPHA_GRID[k],
                       " --scale_ts=yes",
                       " --include_mode=",MODES[j])
    runSlurm <- slurmarray(cmdLines,sname="PP",stime="540",smem="50G",
                           soutdir=paste0('./logs/'), sparallel = 0)
    system(runSlurm)
  }
}


#################################
### Disease-Specific Modeling ###
#################################


### SCALED
for(i in 1:length(DISEASES)){
  cmdLines <- paste0("Rscript --vanilla ",my_path,"/train_lstm_internal_byyeardiseaseonly.R",
                     " --alpha=0.5",
                     " --include_disease=",DISEASES[i],
                     " --warm_start=no",
                     " --scale_ts=yes")
  runSlurm <- slurmarray(cmdLines,sname="PP",stime="540",smem="50G",
                         soutdir=paste0('./logs/'), sparallel = 0)
  system(runSlurm)
}



### SCALED, OTHER ALPHAS
ALPHA_GRID = c(0.025, 0.1, 0.25, 0.75, 0.9, 0.975)
DISEASES = c('Influenza')
for(i in 1:length(DISEASES)){
  for(j in 1:length(ALPHA_GRID)){
    cmdLines <- paste0("Rscript --vanilla ",my_path,"/train_lstm_internal_byyeardiseaseonly.R",
                       " --alpha=",ALPHA_GRID[j],
                       " --include_disease=",DISEASES[i],
                       " --warm_start=no",
                       " --scale_ts=yes")
    runSlurm <- slurmarray(cmdLines,sname="PP",stime="540",smem="50G",
                           soutdir=paste0('./logs/'), sparallel = 0)
    system(runSlurm)
  }
}


##########################################
### Disease x Source Specific Modeling ###
##########################################

### SCALED
for(i in 1:length(DISEASE_SOURCES)){
  cmdLines <- paste0("Rscript --vanilla ",my_path,"/train_lstm_internal_byyearsourceonly.R",
                     " --alpha=0.5",
                     " --include_disease_source=",DISEASE_SOURCES[i],
                     " --warm_start=no",
                     " --scale_ts=yes")
  runSlurm <- slurmarray(cmdLines,sname="PP",stime="540",smem="50G",
                         soutdir=paste0('./logs/'), sparallel = 0)
  system(runSlurm)
}



### SCALED, OTHER ALPHAS
ALPHA_GRID = c(-1)
for(i in 1:length(DISEASE_SOURCES)){
  for(j in 1:length(ALPHA_GRID)){
    cmdLines <- paste0("Rscript --vanilla ",my_path,"/train_lstm_internal_byyearsourceonly.R",
                     " --alpha=",ALPHA_GRID[j],
                     " --include_disease_source=",DISEASE_SOURCES[i],
                     " --warm_start=no",
                     " --scale_ts=yes")
    runSlurm <- slurmarray(cmdLines,sname="PP",stime="540",smem="50G",
                           soutdir=paste0('./logs/'), sparallel = 0)
    system(runSlurm)
  }
}






