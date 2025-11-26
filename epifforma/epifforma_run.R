#########################################
#########################################
### epiFFORMA Training and Evaluation ###
#########################################
#########################################

### This function runs the entire epiFFORMA pipeline starting after the synthetic data library generation. 
### See ./raw_data for implementation of synthetic data generation. 
# Combine all required packages (deduplicated)
required_packages <- unique(c(
  "tsfeatures", "forecast", "randomForest", "data.table",
  "deSolve", "optparse", "here", "dtwclust",
  "ranger", "lightgbm", "e1071", "deepgp", "FNN",
  "reshape2", "plyr", "collapse",
  "ggplot2", "LearnBayes", "LaplacesDemon", "parallel", "doParallel",
  "this.path", "gridExtra", "lubridate", "grid", "plotly", 
  "GGally", "viridis", "ggrepel", "dplyr", "plyr"
))

# Identify missing packages
installed_packages <- rownames(installed.packages())
missing_packages <- setdiff(required_packages, installed_packages)

# Install missing packages
if (length(missing_packages) > 0) {
  message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
  install.packages(missing_packages, dependencies = TRUE)
}

# Load all packages
invisible(lapply(required_packages, library, character.only = TRUE))
setwd(this.path::here()) 

# List of required packages
required_packages <- c(
  "tsfeatures", "forecast", "randomForest", "data.table",
  "ranger", "lightgbm", "e1071", "deepgp", "FNN",
  "reshape2", "plyr", "collapse"
)

# Install any missing packages
installed_packages <- rownames(installed.packages())
missing_packages <- setdiff(required_packages, installed_packages)

if (length(missing_packages) > 0) {
  message("Installing missing packages: ", paste(missing_packages, collapse = ", "))
  install.packages(missing_packages, dependencies = TRUE)
}

#################################
### Step 0: Define File Paths ###
#################################


training_path = './process_data/'
if(!dir.exists(paste0(training_path,"training_data_packets"))){dir.create(paste0(training_path,"training_data_packets"))}
if(!dir.exists(paste0(training_path,"embed_mat"))){dir.create(paste0(training_path,"embed_mat"))}

model_path = './fit_model/'
if(!dir.exists(paste0(model_path,"features2weights"))){dir.create(paste0(model_path,"features2weights"))}
if(!dir.exists(paste0(model_path,"best_params_for_multierror"))){dir.create(paste0(model_path,"best_params_for_multierror"))}

eval_path = './evaluate_model/'
if(!dir.exists(paste0(eval_path,"evaluation"))){dir.create(paste0(eval_path,"evaluation"))}

uq_path = './uq/'
if(!dir.exists(paste0(uq_path,"validation_data_packets"))){dir.create(paste0(uq_path,"validation_data_packets"))}

evalprob_path = './evaluate_probabilistic/'
if(!dir.exists(paste0(eval_path,"evaluation"))){dir.create(paste0(eval_path,"evaluation"))}

synthetic_path <- paste0(training_path, '../raw_data/synthetic/output/') 

source(paste0(training_path,"epi_functions.R"))
source(paste0(training_path,'/../SLURMarray.r'))

######################################################################
### Step 1: Generate Training Data from Existing Synthetic Library ### 
######################################################################


if(!file.exists(paste0(synthetic_path,"synthetic_moa.RDS")) | 
   !file.exists(paste0(synthetic_path,"synthetic_uq.RDS")) | 
   !file.exists(paste0(synthetic_path,"synthetic.RDS"))){
  # setwd(paste0(synthetic_path,"../code/"))
	source(paste0(synthetic_path,"../code/", "make_synthetic_training_data.R"))
  # setwd(this.path::here()) 
}

### Create Embedding Matrices for sMOA ### 
if(!file.exists(paste0(training_path,"embed_mat/embed_mat_y.csv"))){
  h = 4
  synthetic <- readRDS(paste0(synthetic_path,"synthetic_moa.RDS")) #CHANGED TO SEPARATE SYNTHETIC!!!!!!
  embed_mat <- create_embed_matrix(synthetic,h)
  NMAX = 10000000 #don't allow more than 10M snippets
  set.seed('1234')
  inds = sample(1:nrow(embed_mat[[1]]), NMAX, replace = F)
  write.csv(embed_mat[[1]][inds,], file = paste0(training_path,"embed_mat/embed_mat_X.csv"), quote = F, row.names = F)
  write.csv(embed_mat[[2]][inds,], file = paste0(training_path,"embed_mat/embed_mat_y.csv"), quote = F, row.names = F)
  rm(embed_mat)
}
if(!file.exists(paste0(training_path,"embed_mat/embed_mat_y_deriv.csv"))){
  h = 4
  synthetic <- readRDS(paste0(synthetic_path,"synthetic_moa.RDS")) #CHANGED TO SEPARATE SYNTHETIC!!!!!!
  embed_mat_deriv <- create_embed_matrix(synthetic, h = h, k=5)
  embed_mat_X_deriv <- apply(embed_mat_deriv[[1]],1,diff)
  embed_mat_y_deriv <- apply(embed_mat_deriv[[2]],1,diff)
  embed_mat_X_deriv = t(embed_mat_X_deriv)
  embed_mat_y_deriv = t(embed_mat_y_deriv)
  embed_mat_y_deriv = cbind(embed_mat_deriv[[2]][,1] - embed_mat_deriv[[1]][,ncol(embed_mat_deriv[[1]])], embed_mat_y_deriv) #changed on 5/12/25
  NMAX = 10000000 #don't allow more than 10M snippets
  set.seed('4321')
  inds = sample(1:nrow(embed_mat_X_deriv), NMAX, replace = F)
  write.csv(embed_mat_X_deriv[inds,], file = paste0(training_path,"embed_mat/embed_mat_X_deriv.csv"), quote = F, row.names = F)
  write.csv(embed_mat_y_deriv[inds,], file = paste0(training_path,"embed_mat/embed_mat_y_deriv.csv"), quote = F, row.names = F)
  rm(embed_mat_deriv)
}

### Generate Point-Prediction Training Data Packets via SLURM ###
#6 x 300 = 1800 total runs of 100 
for(i in 1:6){
 cmdLines <- paste0("Rscript --vanilla ",training_path,"make_data_internal.R",
                    " --replicate_num=",i,
                    " --num_to_run=",300)
 runSlurm <- slurmarray(cmdLines,sname="PP",stime="540",smem="5G",
                        soutdir=paste0('./logs/'), sparallel = 0)
 system(runSlurm)
}

### Generate UQ Training Data Packets via SLURM ***
#6 x 300 = 1800 total runs of 100 
for(i in 1:6){
  cmdLines <- paste0("Rscript --vanilla ",uq_path,"make_validation_internal.R",
                     " --replicate_num=",i,
                     " --num_to_run=",300)
  runSlurm <- slurmarray(cmdLines,sname="PP",stime="540",smem="5G",
                         soutdir=paste0('./logs/'), sparallel = 0)
  system(runSlurm)
}


### *** Wait for Step 1 to complete before proceeding to Step 2 *** ###



######################################
### Step 2: Train epiFFORMA Models ### 
######################################

### Perform Bayesian optimization for hyperparameters ###
runSlurm = readLines(paste0(model_path,"run_calibration.sh"))
system(runSlurm)

### Train 7 Replicates ###
for(i in 1:7){
  cmdLines <- paste0("Rscript --vanilla ",model_path,"fit_lightgbm_internal.R",
                     " --fit_type=order",
                     " --fit_num=",i,
                     " --param_type=multierror")
  runSlurm <- slurmarray(cmdLines,sname="PP",stime="540",smem="5G",
                         soutdir=paste0('./logs/'), sparallel = 0)
  system(runSlurm)
}

### *** Wait for training to complete before proceeding to aggregation *** ###

### Aggregate Model Fits ###
lf = list.files(paste0(model_path,'features2weights/'))
lf = lf[grepl('multierror', lf) &  !grepl('fitted_lgb_models_order_multierror.RDS', lf)]
mod_wt = replicate(length(lf),list(NULL))
for(i in 1:length(lf)){
  mod_wt[[i]] = readRDS(paste0(model_path,'features2weights/',lf[i]))
}
mod_wt <- unlist(mod_wt)
saveRDS(mod_wt, file = paste0(model_path,'features2weights/',"fitted_lgb_models_order_multierror.RDS"))


### *** Wait for Step 2 to complete before proceeding to Step 3 *** ###




########################################################################
### Step 3: Obtain epiFFORMA Predictions for Real and Synthetic Data ### 
########################################################################

### Get predictions for synthetic training and test data ###
cmdLines <- paste0("Rscript --vanilla ",eval_path,"eval_synthetic_training.R",
                   " --eval_key=synthetic_training", " --eval_type=order", " --param_type=multierror")
runSlurm <- slurmarray(cmdLines,sname="PP",stime="540",smem="5G", soutdir=paste0('./logs/'), sparallel = 0)
system(runSlurm)

cmdLines <- paste0("Rscript --vanilla ",eval_path,"eval_synthetic_test.R",
                   " --eval_key=synthetic_test", " --eval_type=order",  " --param_type=multierror")
runSlurm <- slurmarray(cmdLines,sname="PP",stime="540",smem="5G", soutdir=paste0('./logs/'), sparallel = 0)
system(runSlurm)


### Get Predictions for Real Data ###
EVAL_RUNS = c('us_covid_rollercoaster','global_covid_rollercoaster',
              'us_ili_rollercoaster', 'dengue_rollercoaster',
              'us_measles', 'us_mumps','us_diphtheria',
              'us_rubella','us_polio','us_smallpox',
              'chikungunya')
for(i in 1:length(EVAL_RUNS)){
  if(grepl('measles',EVAL_RUNS[i])|grepl('us_ili_rollercoaster',EVAL_RUNS[i])){
    cmdLines <- paste0("Rscript --vanilla ",eval_path,"eval_real.R",
                       " --eval_key=", EVAL_RUNS[i], " --eval_type=order", " --param_type=multierror")
    extraoptions = c('--qos=long')
    runSlurm <- slurmarray(cmdLines,sname="PP",stime="2500",smem="5G",
                           soutdir=paste0('./logs/'), sparallel = 0,
                           extraOption = extraoptions)
  }else{
    cmdLines <- paste0("Rscript --vanilla ",eval_path,"eval_real.R",
                       " --eval_key=", EVAL_RUNS[i], " --eval_type=order", " --param_type=multierror")
    runSlurm <- slurmarray(cmdLines,sname="PP",stime="540",smem="5G",
                           soutdir=paste0('./logs/'), sparallel = 0)
  }
  system(runSlurm)
}





### *** Wait for Step 3 to complete before proceeding to Step 4 *** ###




#############################################################################
### Step 4: Train Interval-Based epiFFORMA for Uncertainty Quantification ### 
#############################################################################

cmdLines <- paste0("Rscript --vanilla ",uq_path,"make_probabilistic_epifforma.R")
runSlurm <- slurmarray(cmdLines,sname="PP",stime="540",smem="5G", soutdir=paste0('./logs/'), sparallel = 0)
system(runSlurm)

cmdLines <- paste0("Rscript --vanilla ",uq_path,"make_probabilistic_equal_wt.R")
runSlurm <- slurmarray(cmdLines,sname="PP",stime="540",smem="5G", soutdir=paste0('./logs/'), sparallel = 0)
system(runSlurm)

### *** Wait for Step 4 to complete before proceeding to Step 5 *** ###



#########################################################################
### Step 5: Obtain epiFFORMA Predictive Interval Widths for Real Data ### 
#########################################################################

### Get Predictions for Real Data ###
EVAL_RUNS = c('us_covid_rollercoaster','global_covid_rollercoaster',
              'us_ili_rollercoaster', 'dengue_rollercoaster',
              'us_measles', 'us_mumps','us_diphtheria',
              'us_rubella','us_polio','us_smallpox',
              'chikungunya')
for(i in 1:length(EVAL_RUNS)){
  if(grepl('measles',EVAL_RUNS[i])|grepl('us_ili_rollercoaster',EVAL_RUNS[i])){
    cmdLines <- paste0("Rscript --vanilla ",evalprob_path,"eval_probabilistic_internal.R",
                       " --eval_key=", EVAL_RUNS[i], " --eval_type=order")
    extraoptions = c('--qos=long')
    runSlurm <- slurmarray(cmdLines,sname="PP",stime="2500",smem="5G",
                           soutdir=paste0('./logs/'), sparallel = 0,
                           extraOption = extraoptions)
  }else{
    cmdLines <- paste0("Rscript --vanilla ",evalprob_path,"eval_probabilistic_internal.R",
                       " --eval_key=", EVAL_RUNS[i], " --eval_type=order")
    runSlurm <- slurmarray(cmdLines,sname="PP",stime="540",smem="5G",
                           soutdir=paste0('./logs/'), sparallel = 0)
  }
  system(runSlurm)
}


### *** Wait for Step 5 to complete before proceeding to Step 6 *** ###



#############################################################
### Step 6: Summarize Results and Generate Visualizations ### 
#############################################################

source('./create_visualizations/figure_1.R')
source('./create_visualizations/figure_2.R')
source('./create_visualizations/figure_3.R')

