
################################
### Get Global SLURM Options ###
################################

options(warn = -1)

suppressMessages(library("optparse"))

option_list <- list(
  make_option("--total_run", type="character", default=""),
  make_option("--cur_run", type="character", default=""),
  make_option("--scale_ts", type="character", default=""),
  make_option("--include_mode", type="character", default="")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options

total_run = as.numeric(opt$total_run) 
cur_run = as.numeric(opt$cur_run) 
scale_ts = as.character(opt$scale_ts) 
include_mode = as.character(opt$include_mode) 

######################
### Body of Script ###
######################


## load libraries
library(ggplot2)
library(data.table)
library(plyr)
library(gridExtra)
library(lubridate)
library(parallel)
library(doParallel)
library(grid)
library(plotly)
library(GGally)
theme_set(theme_bw())

embedding_path = '~/GitLab/frodo/features/embeddings/'
if(scale_ts == 'yes'){
  output_path = '~/GitLab/frodo/evaluations/smoanodiff/modespecific'
}else{
  output_path = '~/GitLab/frodo/evaluations/smoanodiff_noscale/modespecific'
}


FILES_LONG = list.files(paste0(embedding_path))
FILES_LONG = FILES_LONG[!grepl('real_train_mat.csv',FILES_LONG)]
lut_path = '~/GitLab/frodo/lookup_tables/'
DISEASE_SOURCES = readxl::read_xlsx(paste0(lut_path, 'Evaluations_ML.xlsx'))
DISEASE_SOURCES = data.frame(DISEASE_SOURCES)
DISEASE_SOURCES = DISEASE_SOURCES$FILES
DISEASE_SOURCES = gsub('.RDS','',DISEASE_SOURCES)
FILES_ALL = c()
for(i in 1:length(DISEASE_SOURCES)){
  FILES_ALL = c(FILES_ALL, FILES_LONG[grepl(DISEASE_SOURCES[i],FILES_LONG)])
}
FILES_ALL = unique(FILES_ALL)


### Merge in mode of transmission
FILES_TEMP = gsub('.csv','',unique(gsub('embed_','',FILES_ALL)))
SPLITS = strsplit(FILES_TEMP, split = '_')
DISEASES = rep(NA, length(FILES_TEMP))
for(i in 1:length(SPLITS)){
  if(length(SPLITS[[i]])==3){
    DISEASES[i] = SPLITS[[i]][1]
  }else if(length(SPLITS[[i]])==4){
    DISEASES[i] = paste0(SPLITS[[i]][1],'_',SPLITS[[i]][2])
  }else if(length(SPLITS[[i]])==5){
    DISEASES[i] = paste0(SPLITS[[i]][1],'_',SPLITS[[i]][2],'_',SPLITS[[i]][3])
  }else{
    DISEASES[i] = paste0(SPLITS[[i]][1],'_',SPLITS[[i]][2],'_',SPLITS[[i]][3],'_',SPLITS[[i]][4])
  }
}
lut_path = '~/GitLab/frodo/lookup_tables/'
disease_lut = readxl::read_xlsx(paste0(lut_path, 'Disease_Mappings.xlsx'))
TO_MERGE = data.frame(FILES_ALL = unique(FILES_ALL), FILES = FILES_TEMP, Disease = DISEASES)
TO_MERGE = merge(TO_MERGE, disease_lut[,c('Disease', 'Transmission')], by = 'Disease')


### Subset by transmission mode
if(include_mode == 'respiratory'){
  TO_MERGE = TO_MERGE[TO_MERGE$Transmission == 'Respiratory_Secretions',]
}else if(include_mode == 'fecaloral'){
  TO_MERGE = TO_MERGE[TO_MERGE$Transmission == 'Fecal-Oral' | TO_MERGE$Transmission == 'Water' | TO_MERGE$Transmission == 'Fecal-Oral/Bodily Fluids',]
}else if(include_mode == 'sexual'){
  TO_MERGE = TO_MERGE[TO_MERGE$Transmission == 'Sexual',]
}else if(include_mode == 'vectorborne'){
  TO_MERGE = TO_MERGE[TO_MERGE$Transmission == 'Vectorborne',]
}else{
  stop('include_trans not implemented')
}
FILES_ALL = TO_MERGE$FILES_ALL

FILES_CUR = FILES_ALL[seq(cur_run,length(FILES_ALL), by = total_run)]

####################
### Prepare sMOA ###
####################
library(mgcv)
library(collapse)

### Read in Embeddings 
dat = data.frame(data.table::fread(paste0(embedding_path,"/real_train_mat.csv")))

### Merge in mode of transmission
FILES_TEMP = unique(dat$FILES)
SPLITS = strsplit(FILES_TEMP, split = '_')
DISEASES = rep(NA, length(FILES_TEMP))
for(i in 1:length(SPLITS)){
  if(length(SPLITS[[i]])==3){
    DISEASES[i] = SPLITS[[i]][1]
  }else if(length(SPLITS[[i]])==4){
    DISEASES[i] = paste0(SPLITS[[i]][1],'_',SPLITS[[i]][2])
  }else if(length(SPLITS[[i]])==5){
    DISEASES[i] = paste0(SPLITS[[i]][1],'_',SPLITS[[i]][2],'_',SPLITS[[i]][3])
  }else{
    DISEASES[i] = paste0(SPLITS[[i]][1],'_',SPLITS[[i]][2],'_',SPLITS[[i]][3],'_',SPLITS[[i]][4])
  }
}
lut_path = '~/GitLab/frodo/lookup_tables/'
disease_lut = readxl::read_xlsx(paste0(lut_path, 'Disease_Mappings.xlsx'))
TO_MERGE = data.frame(FILES = FILES_TEMP, Disease = DISEASES)
TO_MERGE = merge(TO_MERGE, disease_lut[,c('Disease', 'Transmission')], by = 'Disease')
dat = merge(dat, TO_MERGE[,c('FILES', 'Transmission')], all.x = T, all.y = F)
 
unique(dat$Transmission)

A = dat[is.na(dat$Transmission),]

### Subset by transmission mode
if(include_mode == 'respiratory'){
  dat = dat[dat$Transmission == 'Respiratory_Secretions',]
}else if(include_mode == 'fecaloral'){
  dat = dat[dat$Transmission == 'Fecal-Oral' | dat$Transmission == 'Water' | dat$Transmission == 'Fecal-Oral/Bodily Fluids',]
}else if(include_mode == 'sexual'){
  dat = dat[dat$Transmission == 'Sexual',]
}else if(include_mode == 'vectorborne'){
  dat = dat[dat$Transmission == 'Vectorborne',]
}else{
  stop('include_trans not implemented')
}

dat = dat[!is.na(dat$Transmission),]
A = dat[is.na(dat$Transmission),]
dat$X1[dat$ts_scale == 'counts'] = round(dat$X1[dat$ts_scale == 'counts'], 0)
dat$X2[dat$ts_scale == 'counts'] = round(dat$X2[dat$ts_scale == 'counts'], 0)
dat$X3[dat$ts_scale == 'counts'] = round(dat$X3[dat$ts_scale == 'counts'], 0)
dat$X4[dat$ts_scale == 'counts'] = round(dat$X4[dat$ts_scale == 'counts'], 0)
dat$X5[dat$ts_scale == 'counts'] = round(dat$X5[dat$ts_scale == 'counts'], 0)
dat$X6[dat$ts_scale == 'counts'] = round(dat$X6[dat$ts_scale == 'counts'], 0)

X = as.matrix(dat[,c('X1','X2','X3','X4','X5','X6')])


X[,1] = pmax(0,X[,1])
X[,2] = pmax(0,X[,2])
X[,3] = pmax(0,X[,3])
X[,4] = pmax(0,X[,4])
X[,5] = pmax(0,X[,5])
X[,6] = pmax(0,X[,6])


y = as.matrix(dat[,c('ysmooth1','ysmooth2','ysmooth3','ysmooth4')])


if(scale_ts == 'yes'){
  center_factor = 0
  scale_factor = X[,6]
  scale_factor[scale_factor==0]=1
  y[,1] = (y[,1]-center_factor)/scale_factor
  y[,2] = (y[,2]-center_factor)/scale_factor
  y[,3] = (y[,3]-center_factor)/scale_factor
  y[,4] = (y[,4]-center_factor)/scale_factor
  
  X[,1] = (X[,1]-center_factor)/scale_factor
  X[,2] = (X[,2]-center_factor)/scale_factor
  X[,3] = (X[,3]-center_factor)/scale_factor
  X[,4] = (X[,4]-center_factor)/scale_factor
  X[,5] = (X[,5]-center_factor)/scale_factor
  X[,6] = (X[,6]-center_factor)/scale_factor
}

dates_embed = dat$date
dates_embed = as.Date(dates_embed,format = '%Y-%m-%d')
### convert to differences

rm(list=c('dat'))


# The following values are all from Murph's sMOA github as of 12/20/24
# https://github.com/lanl/precog/blob/main/smoa/R/run_smoa.R
k <- 5
closest <- 4422
#don't let closest be more than 10% of data
closest = pmin(closest, round(0.1*nrow(X)))
h <- 4


library(mgcv)
library(collapse)


## define number of cores
ncores <- pmin(floor(.5*detectCores()),10)

## define socket type
sockettype <- "PSOCK"

## set up cluster
cl <- parallel::makeCluster(spec = ncores,
                            type = sockettype)
setDefaultCluster(cl)
registerDoParallel(cl)

print(Sys.time())
train_data <- foreach(i=1:length(FILES_CUR), 
                      .errorhandling = "pass",
                      .verbose = T,
                      .packages = c('dplyr', 'tsfeatures', 'timeDate', 'lubridate','mgcv','zoo','forecast','collapse'))%dopar%{  
  eval_key = gsub('_X.csv','',FILES_CUR[i])
  eval_key = gsub('embed_','',eval_key)
  
  dat <- read.csv(paste0(embedding_path,FILES_CUR[i]))
  
  dat$X1[dat$ts_scale == 'counts'] = round(dat$X1[dat$ts_scale == 'counts'], 0)
  dat$X2[dat$ts_scale == 'counts'] = round(dat$X2[dat$ts_scale == 'counts'], 0)
  dat$X3[dat$ts_scale == 'counts'] = round(dat$X3[dat$ts_scale == 'counts'], 0)
  dat$X4[dat$ts_scale == 'counts'] = round(dat$X4[dat$ts_scale == 'counts'], 0)
  dat$X5[dat$ts_scale == 'counts'] = round(dat$X5[dat$ts_scale == 'counts'], 0)
  dat$X6[dat$ts_scale == 'counts'] = round(dat$X6[dat$ts_scale == 'counts'], 0)
  
  
  X_real = as.matrix(dat[,c('X1','X2','X3','X4','X5','X6')])
  y_real = as.matrix(dat[,c('y1','y2','y3','y4')])
  lastobs_real = dat$obs
  dates_real = as.character(dat$date)
  
  if(dat$ts_scale[1] == 'counts'){
    X_real = round(X_real, 0)
  }
  times = as.numeric(difftime(as.Date(dates_real,format = '%Y-%m-%d'),as.Date(paste0('2010-01-01'), format = '%Y-%m-%d'), units = 'days'))
  SUBSET_LONG = which(times >= 0)
  
  if(!file.exists(paste0(output_path,"/real_eval_mat_",eval_key, ".csv"))){
    RESULTS = NULL
    for(j in SUBSET_LONG){
      #### get real embedding row
      to_match_in_moa             <- pmax(0,as.vector(as.numeric(X_real[j,])))
      to_match_in_moa = pmax(0,to_match_in_moa)
      to_match_in_moa_orig = to_match_in_moa
      

      if(scale_ts == 'yes'){
        center_factor = 0
        scale_factor = to_match_in_moa[length(to_match_in_moa)]
        if(scale_factor == 0){
          scale_factor = 1
        }
        to_match_in_moa = (to_match_in_moa-center_factor)/scale_factor
      }else{
        scale_factor = 1
      }
      #### This is where the sMOA calculation actually happens! 
      #### take our matrix to match (X_diff) and efficiently compute distances from to_match_in_moa to X_diff
      times = as.numeric(difftime(as.Date(dates_real[j],format = '%Y-%m-%d'),dates_embed, units = 'days'))
      SUBSET = which(times > 30) #at least 30 days after to ensure data are observed
      if(length(SUBSET) > 0){
        closest_temp = pmax(pmin(closest, floor(0.1*nrow(X[SUBSET,]))),2)
        dist_to_test <- rowSums(abs(X[SUBSET,2:6] %r-% tail(to_match_in_moa,k)))
        
        #### get the closest ids in the X_diff mat
        dists_of_closest = DescTools::Small(dist_to_test, k = closest_temp, unique = FALSE, na.last = NA)
        closest_ids = which(dist_to_test<dists_of_closest[length(dists_of_closest)])
        closest_ids2 = c(closest_ids, which(dist_to_test == dists_of_closest[length(dists_of_closest)]))
        closest_ids2 = closest_ids2[1:pmin(closest_temp,length(closest_ids2))]
        point <- pmax(0,scale_factor*head(apply(y[SUBSET,][closest_ids2,],2,median),h))
        if(lastobs_real[j] == 0){
          point = 0
        }
        DAT = data.frame(row_num = j, h = c(1:4), truth = as.vector(as.numeric(y_real[j,])),
                         obs = lastobs_real[j], fcst = point, 
                         date = dates_real[j],
                         min_dist = dists_of_closest[1], 
                         max_dist = dists_of_closest[length(dists_of_closest)])
        RESULTS <- rbind(RESULTS, DAT)
      }
      print(j)
    }
    write.csv(RESULTS, file = paste0(output_path,"/real_eval_mat_",eval_key, ".csv"), quote = F, row.names = F)
  }

  rm('RESULTS')
  xx <- 1
  xx
}
stopCluster(cl)



