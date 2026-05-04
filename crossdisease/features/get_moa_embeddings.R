######################################
######################################
### Generate MOA Snippet Libraries ###
######################################
######################################

my_path = '~/GitLab/frodo/features/'
setwd(my_path)
source(paste0(my_path,'/../SLURMarray.r'))

## define paths
savetrainpath <- my_path


#####################################
### Generate Embeddings via SLURM ###
#####################################

### Training Data for Primary Analysis
if(!dir.exists(paste0(my_path,"embeddings"))){dir.create(paste0(my_path,"embeddings"))}
NUM_RUN = 10
for(i in 1:NUM_RUN){
  cmdLines <- paste0("Rscript --vanilla ",my_path,"/get_moa_embeddings_gam_internal.R",
                     " --total_run=",NUM_RUN,
                     " --cur_run=",i)
  runSlurm <- slurmarray(cmdLines,sname="PP",stime="540",smem="5G",
                         soutdir=paste0('./logs/'), sparallel = 0)
  system(runSlurm)
}

### Additional Training Data using M4 Time Series
NUM_RUN = 3
for(i in 1:NUM_RUN){
  cmdLines <- paste0("Rscript --vanilla ",my_path,"/get_moa_embeddings_gam_internal_m4.R",
                       " --total_run=",NUM_RUN,
                       " --cur_run=",i)
  runSlurm <- slurmarray(cmdLines,sname="PP",stime="540",smem="5G",
                           soutdir=paste0('./logs/'), sparallel = 0)
  system(runSlurm)
}

### Additional Training Data using M4 Time Series
NUM_RUN = 3
for(i in 1:NUM_RUN){
  cmdLines <- paste0("Rscript --vanilla ",my_path,"/get_moa_embeddings_gam_internal_m4daily.R",
                     " --total_run=",NUM_RUN,
                     " --cur_run=",i)
  runSlurm <- slurmarray(cmdLines,sname="PP",stime="540",smem="5G",
                         soutdir=paste0('./logs/'), sparallel = 0)
  system(runSlurm)
}


####################################
### Merge Giant Training Dataset ###
####################################

# library(data.table)
# library(parallel)
# library(doParallel)
# 
# ncores <- floor(.9*detectCores())
# 
# ## define socket type
# sockettype <- "PSOCK"
# 
# lf <- list.files(paste0(savetrainpath,"embeddings/"))
# lf = lf[lf != 'real_train_mat.csv']
# lf = lf[!is.na(lf)]
# 
# lut_path = '~/GitLab/frodo/lookup_tables/'
# DISEASES = readxl::read_xlsx(paste0(lut_path, 'Disease_ML.xlsx'))
# DISEASES = data.frame(DISEASES)
# DISEASES = DISEASES$FILES
# lf_all = c()
# for(j in 1:length(DISEASES)){
#   lf_all = c(lf_all, lf[grepl(DISEASES[j], lf)])
# }
# lf_all = unique(lf_all)
# 
# 
# ## set up cluster
# cl <- parallel::makeCluster(spec = ncores,
#                             type = sockettype)
# setDefaultCluster(cl)
# registerDoParallel(cl)
# print(Sys.time())
# train_data <- foreach(i=1:length(lf), 
#                       #.errorhandling = "pass",
#                       #.verbose = T,
#                       .packages = c('dplyr', 'tsfeatures', 'data.table'))%dopar%{
#                         dat = data.table::fread(paste0(paste0(savetrainpath,"embeddings/"),lf[i]))
#                         dat$date = as.character(dat$date)
#                         dat = data.frame(dat)
#                         if(lf[i] %in% lf_all){
#                           scale_factor = 1
#                         }else{
#                           scale_factor = 2
#                         }
#                         dat = dat[sample(1:nrow(dat),size = floor(nrow(dat)/scale_factor), replace = F ),]
#                         dat
#                       }
# stopCluster(cl)
# print(Sys.time())
# 
# train_final = do.call('rbind',train_data) #big big data merge
# data.table::fwrite(x=train_final,file =paste0(savetrainpath,"/embeddings/real_train_mat.csv"),quote = FALSE, row.names = FALSE)
# 
# 



# library(data.table)
# library(parallel)
# library(doParallel)
# 
# ncores <- floor(.9*detectCores())
# 
# ## define socket type
# sockettype <- "PSOCK"
# 
# lf <- list.files(paste0(savetrainpath,"embeddings/"))
# lf = lf[lf != 'real_train_mat.csv']
# lf = lf[!is.na(lf)]
# 
# lut_path = '~/GitLab/frodo/lookup_tables/'
# DISEASES = readxl::read_xlsx(paste0(lut_path, 'Disease_ML.xlsx'))
# DISEASES = data.frame(DISEASES)
# DISEASES = DISEASES$FILES
# lf_all = c()
# for(j in 1:length(DISEASES)){
#   lf_all = c(lf_all, lf[grepl(DISEASES[j], lf)])
# }
# lf_all = unique(lf_all)
# 
# 
# ## set up cluster
# cl <- parallel::makeCluster(spec = ncores,
#                             type = sockettype)
# setDefaultCluster(cl)
# registerDoParallel(cl)
# print(Sys.time())
# train_data <- foreach(i=1:length(lf), 
#                       #.errorhandling = "pass",
#                       #.verbose = T,
#                       .packages = c('dplyr', 'tsfeatures', 'data.table'))%dopar%{
#                         dat = data.table::fread(paste0(paste0(savetrainpath,"embeddings/"),lf[i]))
#                         dat = dat[dat$obs > 0,]
#                         dat$date = as.character(dat$date)
#                         dat = data.frame(dat)
#                         if(lf[i] %in% lf_all){
#                           scale_factor = 1
#                         }else{
#                           scale_factor = 1
#                         }
#                         dat = dat[sample(1:nrow(dat),size = floor(nrow(dat)/scale_factor), replace = F ),]
#                         dat
#                       }
# stopCluster(cl)
# print(Sys.time())
# 
# train_final = do.call('rbind',train_data) #big big data merge
# data.table::fwrite(x=train_final,file =paste0(savetrainpath,"/embeddings/real_train_mat.csv"),quote = FALSE, row.names = FALSE)









####################################
### Merge M4 Training Dataset ###
####################################

library(data.table)
library(parallel)
library(doParallel)

ncores <- floor(.9*detectCores())

## define socket type
sockettype <- "PSOCK"

lf <- list.files(paste0(savetrainpath,"embeddings_m4/"))
lf = lf[lf != 'real_train_mat.csv']
lf = lf[!is.na(lf)]
lf = lf[grepl('Weekly',lf)]


## set up cluster
cl <- parallel::makeCluster(spec = ncores,
                            type = sockettype)
setDefaultCluster(cl)
registerDoParallel(cl)
print(Sys.time())
train_data <- foreach(i=1:length(lf), 
                      #.errorhandling = "pass",
                      #.verbose = T,
                      .packages = c('dplyr', 'tsfeatures', 'data.table'))%dopar%{
                        dat = data.table::fread(paste0(paste0(savetrainpath,"embeddings_m4/"),lf[i]))
                        dat$date = as.character(dat$date)
                        dat = data.frame(dat)
                        scale_factor = 1
                        dat = dat[sample(1:nrow(dat),size = floor(nrow(dat)/scale_factor), replace = F ),]
                        dat
                      }
stopCluster(cl)
print(Sys.time())

train_final = do.call('rbind',train_data) #big big data merge
data.table::fwrite(x=train_final,file =paste0(savetrainpath,"/embeddings_m4/real_train_mat.csv"),quote = FALSE, row.names = FALSE)



