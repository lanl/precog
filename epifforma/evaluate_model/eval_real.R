###############################################
###############################################
### epiFFORMA Predictions for Real Datasets ###
###############################################
###############################################

################################
### Get Global SLURM Options ###
################################

options(warn = -1)

suppressMessages(library("optparse"))

option_list <- list(
  make_option("--eval_key", type="character", default=""),
  make_option("--eval_type", type="character", default=""),
  make_option("--param_type", type="character", default="")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options

eval_key = as.character(opt$eval_key) 
eval_type = as.character(opt$eval_type) 
param_type = as.character(opt$param_type) 
print(paste0(eval_type, ' : ', eval_key, ' : ', param_type))


######################
### Body of Script ###
######################


## load libraries
library(data.table)
library(plyr)
library(gridExtra)
library(lubridate)
library(parallel)
library(doParallel)
library(grid)


## define socket type
sockettype <- "PSOCK"

library(this.path)
my_path = this.path::here()
setwd(my_path)

## define paths
codepath <- paste0(my_path,"/../process_data/")
f2wpath <- paste0(my_path,"/../fit_model/features2weights/")
savepath <- paste0(my_path,"/evaluation/")
savetrainpath =  './process_data/'
outputname <- eval_key

## define forecast horizon
h = 4


## source in function
source(paste0(codepath,"epi_functions.R"))


## read in fitted model
if(eval_type == 'order'){
  if(param_type == 'multilogloss'){
    feature2wt <- readRDS(paste0(f2wpath,"fitted_lgb_models_order_multilogloss.RDS"))
  }else if(param_type == 'multierror'){
    feature2wt <- readRDS(paste0(f2wpath,"fitted_lgb_models_order_multierror.RDS"))
  }else{
    stop('Invalid param_type')
  }
  pred_func <- predict_epifforma_orderingtrunc
}else{
  stop('Invalid eval_type')
}



### Read in embedding matrices
embed_mat_X <- data.table::fread(file=paste0(savetrainpath,"embed_mat/embed_mat_X.csv"))
embed_mat_y <- data.table::fread(file=paste0(savetrainpath,"embed_mat/embed_mat_y.csv"))

embed_mat_X_deriv <- data.table::fread(file=paste0(savetrainpath,"embed_mat/embed_mat_X_deriv.csv"))
embed_mat_y_deriv <- data.table::fread(file=paste0(savetrainpath,"embed_mat/embed_mat_y_deriv.csv"))


##############
### DENGUE ###
##############
if(eval_key == 'dengue_rollercoaster'){
  
  ## load data
  test <- readRDS(paste0(my_path,"/../raw_data/dengue/output/","dengue.RDS"))
  
  ## number of cores, can differ based on dataset size
  ncores <- pmin(floor(.5*detectCores()),5)
  
  ## determine the weekly, incidenct cases, of US states for the "roller coaster"
  criteria_list <- list()
  criteria_list[[1]] <- which(unlist(lapply(test,function(ll){return(ll$ts_time_cadence)})) == "weekly")
  criteria_list[[2]] <- which(unlist(lapply(test,function(ll){return(ll$ts_disease == "dengue")})) == 1)
  criteria_list[[3]] <- which(unlist(lapply(test,function(ll){return(ll$ts_measurement_type)})) == "incidence")
  
  ## get runs
  testingids = rev(sort(Reduce("intersect", criteria_list)))
  mag_df <- data.frame(tsid = testingids,
                       geography = unlist(lapply(test[testingids],function(ll){return(ll$ts_geography)})),
                       disease = unlist(lapply(test[testingids],function(ll){return(ll$ts_disease)})),
                       tot = unlist(lapply(test[testingids],function(ll){return(sum(ll$ts))})))
  mag_df <- mag_df[order(mag_df$tot, decreasing = T),]
  testingids <- mag_df$tsid
}


##############################
### US COVID Rollercoaster ###
##############################
if(eval_key == 'us_covid_rollercoaster'){
  
  ## load data
  test <- readRDS(paste0(my_path,"/../raw_data/global_covid/output/","global_covid.RDS"))
  
  ## number of cores, can differ based on dataset size
  ncores <- pmin(floor(.5*detectCores()),5)
  
  ## determine the weekly, incidenct cases, of US states for the "roller coaster"
  criteria_list <- list()
  criteria_list[[1]] <- which(unlist(lapply(test,function(ll){return(ll$ts_time_cadence)})) == "weekly")
  criteria_list[[2]] <- which(unlist(lapply(test,function(ll){return(length(grep("United States",ll$ts_geography)))})) > 0)
  criteria_list[[3]] <- which(unlist(lapply(test,function(ll){return(ll$ts_disease == "covid")})) == 1)
  criteria_list[[4]] <- which(unlist(lapply(test,function(ll){return(ll$ts_measurement_type)})) == "incident_cases")
  
  ## get runs
  testingids = rev(sort(Reduce("intersect", criteria_list)))
  mag_df <- data.frame(tsid = testingids,
                       geography = unlist(lapply(test[testingids],function(ll){return(ll$ts_geography)})),
                       disease = unlist(lapply(test[testingids],function(ll){return(ll$ts_disease)})),
                       tot = unlist(lapply(test[testingids],function(ll){return(sum(ll$ts))})))
  mag_df <- mag_df[order(mag_df$tot, decreasing = T),]
  mag_df <- subset(mag_df, geography != "United States")
  testingids <- mag_df$tsid
  
  ## process Florida, right roll to avoid contamination with future!!!!
  FLORIDA = test[[mag_df$tsid[grepl('Florida', mag_df$geography)]]]
  FLORIDA$ts[111:length(FLORIDA$ts)] =  zoo::rollapply(FLORIDA$ts[111:length(FLORIDA$ts)], align = 'right', width = 4, FUN = function(x){mean(x,na.rm=T)}, partial = T)
  test[[mag_df$tsid[grepl('Florida', mag_df$geography)]]] = FLORIDA
}



##############################
### Global COVID Rollercoaster ###
##############################
if(eval_key == 'global_covid_rollercoaster'){
  
  ## load data
  test <- readRDS(paste0(my_path,"/../raw_data/global_covid/output/","global_covid.RDS"))
  
  ## number of cores, can differ based on dataset size
  ncores <- pmin(floor(.5*detectCores()),5)
  
  COUNTRIES = c(  "Afghanistan"      ,                "Albania",                          "Algeria"         ,                
                  "Andorra"           ,               "Angola"  ,                         "Antigua and Barbuda",             
                  "Argentina"          ,              "Armenia"  ,                        "Australia"           ,            
                  "Austria"             ,             "Azerbaijan",                       "Bahrain"              ,           
                  "Bangladesh"           ,            "Barbados"   ,                      "Belarus"               ,          
                  "Belgium"               ,           "Belize"      ,                     "Benin"                  ,         
                  "Bhutan"                 ,          "Bolivia"      ,                    "Bosnia and Herzegovina"  ,        
                  "Botswana"                ,         "Brazil"        ,                   "Brunei"                   ,       
                  "Bulgaria"                 ,        "Burkina Faso"   ,                  "Burundi"                   ,      
                  "Cambodia"                  ,       "Cameroon"        ,                 "Canada"                     ,     
                  "Cape Verde"   ,                    "Central African Republic"   ,      "Chad"                        ,    
                  "Chile"         ,                   "China"        ,                    "Colombia"                     ,   
                  "Comoros"        ,                  "Congo (Brazzaville)" ,             "Congo (Kinshasa)" ,               
                  "Costa Rica"      ,                 "Cote d'Ivoire"        ,            "Croatia"           ,              
                  "Cuba"             ,                "Cyprus"                ,           "Czechia"            ,             
                  "Denmark"           ,               "Djibouti"               ,          "Dominica"            ,            
                  "Dominican Republic" ,              "East Timor"              ,         "Ecuador"              ,           
                  "Egypt"               ,             "El Salvador"              ,        "Equatorial Guinea"     ,          
                  "Eritrea"              ,            "Estonia"                   ,       "Eswatini"               ,         
                  "Ethiopia"              ,           "Federated States of Micronesia",   "Fiji"                    ,        
                  "Finland"                ,          "France" ,                          "Gabon"                    ,       
                  "Georgia" ,                         "Germany" ,                         "Ghana"                     ,      
                  "Greece"   ,                        "Grenada"  ,                        "Guatemala"                  ,     
                  "Guinea"    ,                       "Guinea-Bissau" ,                   "Guyana"                      ,    
                  "Haiti"      ,                      "Honduras"       ,                  "Hungary"                      ,   
                  "Iceland"     ,                     "India"           ,                 "Indonesia" ,                      
                  "Iran"         ,                    "Iraq"             ,                "Ireland"    ,                     
                  "Israel"        ,                   "Italy"             ,               "Jamaica"     ,                    
                  "Japan"          ,                  "Jordan"             ,              "Kazakhstan"   ,                   
                  "Kenya"           ,                 "Kiribati"            ,             "Kosovo"        ,                  
                  "Kuwait"           ,                "Kyrgyzstan"           ,            "Laos"           ,                 
                  "Latvia"            ,               "Lebanon"               ,           "Lesotho"         ,                
                  "Liberia"            ,              "Libya"                  ,          "Liechtenstein"    ,               
                  "Lithuania"           ,             "Luxembourg"              ,         "Madagascar"        ,              
                  "Malawi"               ,            "Malaysia"                 ,        "Maldives"           ,             
                  "Mali"                  ,           "Malta"                     ,       "Marshall Islands"    ,            
                  "Mauritania"             ,          "Mauritius"                  ,      "Mexico"               ,           
                  "Moldova"                 ,         "Monaco"                      ,     "Mongolia"              ,          
                  "Montenegro"               ,        "Morocco"                      ,    "Mozambique"             ,         
                  "Myanmar"                   ,       "Nauru"                         ,   "Nepal"                   ,        
                  "Netherlands"                ,      "New Zealand" ,                     "Nicaragua"                ,       
                  "Niger"                       ,     "Nigeria"      ,                   "North Korea"                ,     
                  "North Macedonia"              ,    "Norway"        ,                  "Oman"      ,                      
                  "Pakistan" ,                        "Palau"           ,                 "Palestine" ,                      
                  "Panama"    ,                       "Papua New Guinea" ,                "Paraguay"   ,                     
                  "Peru"       ,                      "Philippines"       ,               "Poland"      ,                    
                  "Portugal"    ,                     "Qatar"              ,              "Romania"      ,                   
                  "Russia"       ,                    "Rwanda"              ,             "Saint Kitts and Nevis" ,          
                  "Saint Lucia"   ,                   "Saint Vincent and the Grenadines", "Samoa"                  ,         
                  "San Marino"     ,                  "Sao Tome and Principe" ,           "Saudi Arabia"            ,        
                  "Senegal"         ,                 "Serbia"                  ,         "Seychelles"               ,       
                  "Sierra Leone"     ,                "Singapore"                ,        "Slovakia"                  ,      
                  "Slovenia"          ,               "Solomon Islands"           ,       "Somalia"                    ,     
                  "South Africa"       ,              "South Korea"                ,      "South Sudan"                 ,    
                  "Spain"               ,             "Sri Lanka"                   ,     "Sudan"                        ,   
                  "Suriname"             ,            "Sweden"                       ,    "Switzerland" ,                    
                  "Syria"                 ,           "Taiwan"                        ,   "Tajikistan"   ,                   
                  "Tanzania"               ,          "Thailand"                       ,  "The Bahamas"   ,                  
                  "The Gambia"              ,         "Togo"       ,                      "Tonga"          ,                 
                  "Trinidad and Tobago"      ,        "Tunisia"     ,                     "Turkey"          ,                
                  "Tuvalu"                    ,       "Uganda"       ,                    "Ukraine"          ,               
                  "United Arab Emirates"       ,      "United Kingdom",                   "United States"     ,              
                  "Uruguay"                     ,     "Uzbekistan"     ,                  "Vanuatu"            ,             
                  "Vatican City"                 ,    "Venezuela"       ,                 "Vietnam"             ,            
                  "Western Sahara"                ,   "Yemen"            ,                "Zambia"               ,           
                  "Zimbabwe"   )
  
  ### drop 3 locations with less than 1000 total cases
  COUNTRIES = COUNTRIES[!(COUNTRIES %in% c('North Korea', 'Vatican City', 'Western Sahara'))]
  
  ## determine the weekly, incident cases, of US states for the "roller coaster"
  criteria_list <- list()
  criteria_list[[1]] <- which(unlist(lapply(test,function(ll){return(ll$ts_time_cadence)})) == "weekly")
  criteria_list[[2]] <- which(unlist(lapply(test,function(ll){return(as.numeric(ll$ts_geography %in% COUNTRIES))})) > 0)
  criteria_list[[3]] <- which(unlist(lapply(test,function(ll){return(ll$ts_disease == "covid")})) == 1)
  criteria_list[[4]] <- which(unlist(lapply(test,function(ll){return(ll$ts_measurement_type)})) == "incident_cases")
  
  ## get runs
  testingids = rev(sort(Reduce("intersect", criteria_list)))
  mag_df <- data.frame(tsid = testingids,
                       geography = unlist(lapply(test[testingids],function(ll){return(ll$ts_geography)})),
                       disease = unlist(lapply(test[testingids],function(ll){return(ll$ts_disease)})),
                       tot = unlist(lapply(test[testingids],function(ll){return(sum(ll$ts))})))
  mag_df <- mag_df[order(mag_df$tot, decreasing = T),]
  testingids <- mag_df$tsid
  
}



#########################
### ILI Rollercoaster ###
#########################
if(eval_key == 'us_ili_rollercoaster'){
  
  ## load data
  test <- readRDS(paste0(my_path,"/../raw_data/us_ili/output/","us_ili.RDS"))
  
  ## number of cores, can differ based on dataset size
  ncores <- pmin(floor(.5*detectCores()),5)
  
  ## determine the weekly, incidenct cases, of US states for the "roller coaster"
  criteria_list <- list()
  criteria_list[[1]] <- which(unlist(lapply(test,function(ll){return(ll$ts_time_cadence)})) == "weekly")
  criteria_list[[2]] <- which(unlist(lapply(test,function(ll){return(ll$ts_disease == "influenza-like illness")})) == 1)
  criteria_list[[3]] <- which(unlist(lapply(test,function(ll){return(ll$ts_measurement_type)})) == "incidence")
  criteria_list[[4]] <- which(unlist(lapply(test,function(ll){return(year(ll$ts_last_time) - year(ll$ts_first_time))})) > 1)
  criteria_list[[5]] <- which(unlist(lapply(test,function(ll){return(length(grep("region",tolower(ll$ts_geography))))})) == 0)
  
  ## get runs
  testingids = rev(sort(Reduce("intersect", criteria_list)))
  mag_df <- data.frame(tsid = testingids,
                       geography = unlist(lapply(test[testingids],function(ll){return(ll$ts_geography)})),
                       disease = unlist(lapply(test[testingids],function(ll){return(ll$ts_disease)})),
                       tot = unlist(lapply(test[testingids],function(ll){return(sum(ll$ts))})))
  mag_df <- mag_df[order(mag_df$tot, decreasing = T),]
  mag_df <- subset(mag_df, geography != "United States")
  testingids <- mag_df$tsid
}


################
### US Polio ###
################
if(eval_key == 'us_polio'){
  
  ## load data
  test <- readRDS(paste0(my_path,"/../raw_data/project_tycho/output/","us_project_tycho.RDS"))
  
  ## number of cores, can differ based on dataset size
  ncores <- pmin(floor(.5*detectCores()),5)
  
  ## determine the weekly, incidenct cases, of US states for the "roller coaster"
  DISEASES = unique(unlist(lapply(test,function(ll){return(ll$ts_disease)})))
  criteria_list <- list()
  criteria_list[[1]] <- which(unlist(lapply(test,function(ll){return(ll$ts_disease == "POLIO")})) == 1)
  
  ## get runs
  testingids = rev(sort(Reduce("intersect", criteria_list)))
  mag_df <- data.frame(tsid = testingids,
                       geography = unlist(lapply(test[testingids],function(ll){return(ll$ts_geography)})),
                       disease = unlist(lapply(test[testingids],function(ll){return(ll$ts_disease)})),
                       tot = unlist(lapply(test[testingids],function(ll){return(sum(ll$ts))})))
  mag_df <- mag_df[order(mag_df$tot, decreasing = T),]
  testingids <- mag_df$tsid
}


#####################
### US Diphtheria ###
#####################
if(eval_key == 'us_diphtheria'){
  
  ## load data
  test <- readRDS(paste0(my_path,"/../raw_data/project_tycho/output/","us_project_tycho.RDS"))
  
  ## number of cores, can differ based on dataset size
  ncores <- pmin(floor(.5*detectCores()),5)
  
  ## determine the weekly, incidenct cases, of US states for the "roller coaster"
  DISEASES = unique(unlist(lapply(test,function(ll){return(ll$ts_disease)})))
  criteria_list <- list()
  criteria_list[[1]] <- which(unlist(lapply(test,function(ll){return(ll$ts_disease == "DIPHTHERIA")})) == 1)
  
  ## get runs
  testingids = rev(sort(Reduce("intersect", criteria_list)))
  mag_df <- data.frame(tsid = testingids,
                       geography = unlist(lapply(test[testingids],function(ll){return(ll$ts_geography)})),
                       disease = unlist(lapply(test[testingids],function(ll){return(ll$ts_disease)})),
                       tot = unlist(lapply(test[testingids],function(ll){return(sum(ll$ts))})))
  mag_df <- mag_df[order(mag_df$tot, decreasing = T),]
  testingids <- mag_df$tsid
}



###################
### US Smallpox ###
###################
if(eval_key == 'us_smallpox'){
  
  ## load data
  test <- readRDS(paste0(my_path,"/../raw_data/project_tycho/output/","us_project_tycho.RDS"))
  
  ## number of cores, can differ based on dataset size
  ncores <- pmin(floor(.5*detectCores()),5)
  
  ## determine the weekly, incidenct cases, of US states for the "roller coaster"
  DISEASES = unique(unlist(lapply(test,function(ll){return(ll$ts_disease)})))
  criteria_list <- list()
  criteria_list[[1]] <- which(unlist(lapply(test,function(ll){return(ll$ts_disease == "SMALLPOX")})) == 1)
  
  ## get runs
  testingids = rev(sort(Reduce("intersect", criteria_list)))
  mag_df <- data.frame(tsid = testingids,
                       geography = unlist(lapply(test[testingids],function(ll){return(ll$ts_geography)})),
                       disease = unlist(lapply(test[testingids],function(ll){return(ll$ts_disease)})),
                       tot = unlist(lapply(test[testingids],function(ll){return(sum(ll$ts))})))
  mag_df <- mag_df[order(mag_df$tot, decreasing = T),]
  testingids <- mag_df$tsid
}




################
### US Mumps ###
################
if(eval_key == 'us_mumps'){
  
  ## load data
  test <- readRDS(paste0(my_path,"/../raw_data/project_tycho/output/","us_project_tycho.RDS"))
  
  ## number of cores, can differ based on dataset size
  ncores <- pmin(floor(.5*detectCores()),5)
  
  ## determine the weekly, incidenct cases, of US states for the "roller coaster"
  DISEASES = unique(unlist(lapply(test,function(ll){return(ll$ts_disease)})))
  criteria_list <- list()
  criteria_list[[1]] <- which(unlist(lapply(test,function(ll){return(ll$ts_disease == "MUMPS")})) == 1)
  
  ## get runs
  testingids = rev(sort(Reduce("intersect", criteria_list)))
  mag_df <- data.frame(tsid = testingids,
                       geography = unlist(lapply(test[testingids],function(ll){return(ll$ts_geography)})),
                       disease = unlist(lapply(test[testingids],function(ll){return(ll$ts_disease)})),
                       tot = unlist(lapply(test[testingids],function(ll){return(sum(ll$ts))})))
  mag_df <- mag_df[order(mag_df$tot, decreasing = T),]
  testingids <- mag_df$tsid
}


##################
### US Measles ###
##################
if(eval_key == 'us_measles'){
  
  ## load data
  test <- readRDS(paste0(my_path,"/../raw_data/project_tycho/output/","us_project_tycho.RDS"))
  
  ## number of cores, can differ based on dataset size
  ncores <- pmin(floor(.5*detectCores()),5)
  
  ## determine the weekly, incidenct cases, of US states for the "roller coaster"
  DISEASES = unique(unlist(lapply(test,function(ll){return(ll$ts_disease)})))
  criteria_list <- list()
  criteria_list[[1]] <- which(unlist(lapply(test,function(ll){return(ll$ts_disease == "MEASLES")})) == 1)
  
  ## get runs
  testingids = rev(sort(Reduce("intersect", criteria_list)))
  mag_df <- data.frame(tsid = testingids,
                       geography = unlist(lapply(test[testingids],function(ll){return(ll$ts_geography)})),
                       disease = unlist(lapply(test[testingids],function(ll){return(ll$ts_disease)})),
                       tot = unlist(lapply(test[testingids],function(ll){return(sum(ll$ts))})))
  mag_df <- mag_df[order(mag_df$tot, decreasing = T),]
  testingids <- mag_df$tsid
}

##################
### US Rubella ###
##################
if(eval_key == 'us_rubella'){
  
  ## load data
  test <- readRDS(paste0(my_path,"/../raw_data/project_tycho/output/","us_project_tycho.RDS"))
  
  ## number of cores, can differ based on dataset size
  ncores <- pmin(floor(.5*detectCores()),5)
  
  ## determine the weekly, incidenct cases, of US states for the "roller coaster"
  DISEASES = unique(unlist(lapply(test,function(ll){return(ll$ts_disease)})))
  criteria_list <- list()
  criteria_list[[1]] <- which(unlist(lapply(test,function(ll){return(ll$ts_disease == "RUBELLA")})) == 1)
  
  ## get runs
  testingids = rev(sort(Reduce("intersect", criteria_list)))
  mag_df <- data.frame(tsid = testingids,
                       geography = unlist(lapply(test[testingids],function(ll){return(ll$ts_geography)})),
                       disease = unlist(lapply(test[testingids],function(ll){return(ll$ts_disease)})),
                       tot = unlist(lapply(test[testingids],function(ll){return(sum(ll$ts))})))
  mag_df <- mag_df[order(mag_df$tot, decreasing = T),]
  testingids <- mag_df$tsid
}









##################
### Chikungunya ###
##################
if(eval_key == 'chikungunya'){
  
  ## load data
  test <- readRDS(paste0(my_path,"/../raw_data/chikungunya/output/","chikv.RDS"))
  
  ## number of cores, can differ based on dataset size
  ncores <- pmin(floor(.5*detectCores()),5)
  
  ## determine the weekly, incidenct cases, of US states for the "roller coaster"
  DISEASES = unique(unlist(lapply(test,function(ll){return(ll$ts_disease)})))
  criteria_list <- list()
  criteria_list[[1]] <- which(unlist(lapply(test,function(ll){return(ll$ts_disease == "chikv")})) == 1)
  
  ## get runs
  testingids = rev(sort(Reduce("intersect", criteria_list)))
  mag_df <- data.frame(tsid = testingids,
                       geography = unlist(lapply(test[testingids],function(ll){return(ll$ts_geography)})),
                       disease = unlist(lapply(test[testingids],function(ll){return(ll$ts_disease)})),
                       tot = unlist(lapply(test[testingids],function(ll){return(sum(ll$ts))})))
  mag_df <- mag_df[order(mag_df$tot, decreasing = T),]
  testingids <- mag_df$tsid
}




#####################################
### Apply Basic Outlier Screening ###
#####################################
outlier_ratios = unlist(lapply(test[testingids],function(ll){return(max(ll$ts)/max(ll$ts[-which.max(ll$ts)]))}))
for(j in testingids[outlier_ratios>10]){
  RATIO = max(test[[j]]$ts)/max(test[[j]]$ts[-which.max(test[[j]]$ts)])
  cnt = 0
  while(RATIO > 10){
    cnt = cnt + 1
    test[[j]]$ts[which.max(test[[j]]$ts)] = round(0.5*(test[[j]]$ts[which.max(test[[j]]$ts)+1] + test[[j]]$ts[which.max(test[[j]]$ts)-1]))
    RATIO = max(test[[j]]$ts)/max(test[[j]]$ts[-which.max(test[[j]]$ts)])
  }
  print(paste0('Changed: ', cnt))
}


##################################################################
## STEP 4: make predictions for a geography
##################################################################

## combination function for foreach when returning a list
## where each item of the list is meant to be combined via rbind()
comb <- function(x, ...) {
  mapply(rbind,x,...,SIMPLIFY=FALSE)
}



## parallelize over time series, loop within time series
cl <- parallel::makeCluster(spec = ncores,type = sockettype)
setDefaultCluster(cl)
registerDoParallel(cl)
print(Sys.time())
myoutput <- foreach(g = 1:length(testingids),
                    .verbose = F)%dopar%{
                      
                      ## load functions
                      source(paste0(codepath,"epi_functions.R"))
                      
                      ## grab the geography to be forecasted
                      ts_test <- test[[testingids[g]]]
                      
                      ## print out where we are
                      print("")
                      print(paste0(g," of ", length(testingids)))
                      print(paste0(ts_test$ts_geography,"__",ts_test$ts_disease))
                      
                      ## get the forecast indices
                      ## each time series needs at least 11 observations and at least h future observations
                      fcst_indices <- 11:(length(ts_test$ts)-h)
                      
                      bigplotoutput = bigfeatures = bigwts = bigwtssd = bigcomponents = NULL
                      for(i in fcst_indices){
                        
                        ## get info_packet
                        info_packet <- ts_test
                        
                        ## trim the time series
                        info_packet$ts <- info_packet$ts[1:i]
                        if(!is.null(info_packet$ts_dates)){
                          info_packet$ts_dates <- info_packet$ts_dates[1:i]
                        }
                        if(!is.null(info_packet$ts_exogenous)){
                          info_packet$ts_exogenous <- info_packet$ts_exogenous[1:i]
                        }
                        
                        ## put info_packet into predict_epifforma()
                        temp_pred <- pred_func(info_packet = info_packet, h = h, additional_features = ifelse(eval_type == 'bounds', T, F))
                        
                        # hard set forecast to zero if last 5 values all 1s or 0s
                        if(info_packet$ts_scale == 'counts'){
                          if(sum(tail(as.numeric(info_packet$ts),n=5) != 0 & tail(as.numeric(info_packet$ts),n=5) != 1) == 0){
                            temp_pred$output_df$fcst[temp_pred$output_df$model == 'epifforma'] = 0
                          }
                        }
                        
                        # add equal_wt benchmark to epifforma
                        SUB = temp_pred$output_df[temp_pred$output_df$model != 'epifforma',]
                        SUB = SUB[!is.na(SUB$h),]
                        SUB = SUB %>% dplyr::group_by(last_obs_time, h, geography) %>% dplyr::mutate(fcst_mean = mean(fcst, na.rm=T))
                        SUB = SUB[!duplicated(paste0(SUB$last_obs_time, '_', SUB$h, '_', SUB$geography)),]
                        SUB$fcst = SUB$fcst_mean
                        SUB$model = 'equal_wt'
                        temp_pred$output_df = rbind(temp_pred$output_df, SUB[,colnames(temp_pred$output_df)])
                        
                        ## pack up. rbinds will get slow when fcst_indices is super long
                        bigplotoutput <- rbind(bigplotoutput, temp_pred$output_df)
                        bigfeatures <- rbind(bigfeatures, temp_pred$features)
                        bigwts <- rbind(bigwts, temp_pred$pred_wts)
                        bigwtssd <- rbind(bigwtssd, temp_pred$pred_wts_sd)
                        bigcomponents <- rbind(bigcomponents, temp_pred$components)
                        
                        if(i %% 250 == 0){
                          gc() #will slow things down, but trying to conserve memory
                        }
                        print(i)
                      }
                      
                      ## combine output into a list
                      results <- list(plot_df = bigplotoutput,
                                      features = bigfeatures,
                                      weights = bigwts,
                                      weights_sd = bigwtssd,
                                      components = bigcomponents)                          
                      
                      results
                    }                     
stopCluster(cl)
print(Sys.time())

rm(feature2wt)


############################
### Organize the Results ###
############################

## prepare for the run
bigplotoutput <- NULL
bigfeatures <- NULL
bigwts <- NULL
bigwtssd <- NULL
bigcomponents <- NULL
for(g in 1:length(myoutput)){
  bigplotoutput = rbind(bigplotoutput, myoutput[[g]]['plot_df']$plot_df)
  bigfeatures = rbind(bigfeatures, myoutput[[g]]['features']$features)
  bigwts = rbind(bigwts, myoutput[[g]]['weights']$weights)
  bigwtssd = rbind(bigwtssd, myoutput[[g]]['weights_sd']$weights_sd)
  bigcomponents = rbind(bigcomponents, myoutput[[g]]['components']$components)
}
rm(myoutput)



if(nrow(bigplotoutput[bigplotoutput$type == 'fcst',]) < Inf){ #changed so always eval this one
  library(dtwclust)
  dat_agg = bigplotoutput[bigplotoutput$type == 'fcst',] %>% dplyr::group_by(last_obs_time, disease, geography, model) %>% dplyr::mutate(sbd = dtwclust::SBD(fcst, truth, znorm = FALSE, error.check = FALSE, return.shifted = FALSE))
  gc()
  
  ## make the accuracy and print it out
  dat_agg$fcst <- pmax(1e-10, dat_agg$fcst)
  accuracydf <- ddply(dat_agg,.(model),summarise,
                      n = length(fcst),
                      smape_no0 = mean( (abs(fcst - truth)/sum(0.5*(abs(fcst) + abs(truth))))[truth>0]   ),#zeros messing this metric up and being unduly influential. Exclude.
                      mae = mean(abs(fcst - truth)),
                      rmse = sqrt(mean(abs(fcst - truth)^2)),
                      sbd_median = median(sbd[!is.infinite(sbd)]), #shape-based distance
                      spearman = cor(fcst, truth, method = 'spearman'))
  accuracydf$smape_rank <- rank(accuracydf$smape_no0)
  accuracydf$mae_rank <- rank(accuracydf$mae)
  accuracydf$rmse_rank <- rank(accuracydf$rmse)
  accuracydf$sbd_rank <- rank(accuracydf$sbd_median)
  accuracydf$spearman_rank <- rank(accuracydf$spearman)
  accuracydf$avg_rank <- (1/3)*(accuracydf$smape_rank + accuracydf$mae_rank + accuracydf$rmse_rank)
  accuracydf$std_smape <- 1-(accuracydf$smape_no0 - min(accuracydf$smape_no0))/(diff(range(accuracydf$smape_no0)))
  accuracydf$std_mae <- 1-(accuracydf$mae - min(accuracydf$mae))/(diff(range(accuracydf$mae)))
  accuracydf$std_rmse <- 1-(accuracydf$rmse - min(accuracydf$rmse))/(diff(range(accuracydf$rmse)))
  accuracydf$std_comb <- accuracydf$std_mae + accuracydf$std_rmse + accuracydf$std_smape
  accuracydf <- accuracydf[order(accuracydf$mae, decreasing = T),]
  print(accuracydf)
  rm(dat_agg)
  gc()                      
}else{
  ## make the accuracy and print it out
  bigplotoutput$fcst <- pmax(1e-10, bigplotoutput$fcst)
  accuracydf <- ddply(subset(bigplotoutput,type == "fcst"),.(model),summarise,
                      n = length(fcst),
                      smape_no0 = mean( (abs(fcst - truth)/sum(0.5*(abs(fcst) + abs(truth))))[truth>0]   ),#zeros messing this metric up and being unduly influential. Exclude.
                      mae = mean(abs(fcst - truth)),
                      rmse = sqrt(mean(abs(fcst - truth)^2)))
  accuracydf$smape_rank <- rank(accuracydf$smape_no0)
  accuracydf$mae_rank <- rank(accuracydf$mae)
  accuracydf$rmse_rank <- rank(accuracydf$rmse)
  accuracydf$avg_rank <- (1/3)*(accuracydf$smape_rank + accuracydf$mae_rank + accuracydf$rmse_rank)
  accuracydf$std_smape <- 1-(accuracydf$smape_no0 - min(accuracydf$smape_no0))/(diff(range(accuracydf$smape_no0)))
  accuracydf$std_mae <- 1-(accuracydf$mae - min(accuracydf$mae))/(diff(range(accuracydf$mae)))
  accuracydf$std_rmse <- 1-(accuracydf$rmse - min(accuracydf$rmse))/(diff(range(accuracydf$rmse)))
  accuracydf$std_comb <- accuracydf$std_mae + accuracydf$std_rmse + accuracydf$std_smape
  accuracydf <- accuracydf[order(accuracydf$mae, decreasing = T),]
  print(accuracydf)
  gc()      
}

## combine output into a list
output_list <- list(plot_df = bigplotoutput,
                    accuracy_table = accuracydf,
                    features = bigfeatures,
                    weights = bigwts,
                    weights_sd = bigwtssd,
                    components = bigcomponents)

## save the output
saveRDS(object = output_list, file = paste0(savepath, outputname, "_",eval_type,'_',param_type,".RDS"))

