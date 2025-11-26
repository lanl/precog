################################################
################################################
### Get epiFFORMA Predictive Interval Widths ###
################################################
################################################


################################
### Get Global SLURM Options ###
################################

options(warn = -1)

suppressMessages(library("optparse"))

option_list <- list(
  make_option("--eval_key", type="character", default=""),
  make_option("--eval_type", type="character", default="")
)

parser <- OptionParser(usage="%prog [options]", option_list=option_list)
args <- parse_args(parser, positional_arguments = 0)
opt <- args$options

eval_key = as.character(opt$eval_key) 
eval_type = as.character(opt$eval_type) 
print(paste0(eval_type, ' : ', eval_key))


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
my_path = here::here("epifforma", "uq") # this.path::here()

## define paths
# my_path = './uq'
# setwd(my_path)
savepath <- paste0(my_path,"/../evaluate_model/evaluation/")
savepath_new <- paste0(my_path,"/../evaluate_probabilistic/evaluation/")
syntheticpath <- paste0(my_path, '/../raw_data/synthetic/output/') 
savemodelpath <- paste0(my_path, '/features2weights/') 
savetrainpath =  '/../process_data/'

if(!dir.exists(savemodelpath)){dir.create(savemodelpath)}
outputname <- eval_key

## define forecast horizon
h = 4



## get .RDS file names
output_list = readRDS(paste0(savepath, outputname, '_', eval_type, '_multierror.RDS'))
results = output_list['plot_df']$plot_df
results = results[results$type == 'fcst',]

library(dplyr)
results_obs = output_list['plot_df']$plot_df
results_obs = results_obs[results_obs$type == 'obs',]
results_obs = results_obs[!duplicated(paste0(results_obs$x, '_', results_obs$geography, '_', results_obs$disease)),c('x', 'geography', 'disease', 'truth')]

FEATURE_LIST = setdiff(names(output_list['features']$features),c('disease','geography','last_obs_time'))
features_wide = output_list['features']$features
features_wide = merge(features_wide, 
                      data.frame(results[results$model == 'epifforma',c('disease','geography','last_obs_time','h', 'truth')],
                                 fcst_epifforma = unlist(results[results$model == 'epifforma','fcst'])), 
                      by = c('disease','geography','last_obs_time','h'), all.x = T, all.y = F)
features_wide = merge(features_wide, 
                      data.frame(results[results$model == 'equal_wt',c('disease','geography','last_obs_time','h')],
                                 fcst_equalwt = unlist(results[results$model == 'equal_wt','fcst'])), 
                      by = c('disease','geography','last_obs_time','h'), all.x = T, all.y = F)

### Components in Test Data
COMPONENT_LIST = setdiff(names(output_list['components']$components),c('disease','geography','last_obs_time','h'))
components_wide = output_list['components']$components
components_wide$component_sd = apply(output_list['components']$components[,COMPONENT_LIST],1,sd)
weights_wide = output_list['weights']$weights
weights_wide$weight_sd = apply(output_list['weights']$weights[,paste0(COMPONENT_LIST,'_avg') ],1,sd)

features_wide = merge(features_wide, components_wide[,c('component_sd','disease','geography','last_obs_time','h')],
                      by = c('disease','geography','last_obs_time','h'), all.x = T, all.y = F)
features_wide$component_sd = features_wide$component_sd/features_wide$fcst_epifforma #only relevant for epifforma models
features_wide$component_sd[features_wide$component_sd>1 | is.infinite(features_wide$component_sd)]=1
features_wide = merge(features_wide, weights_wide[,c('weight_sd','disease','geography','last_obs_time','h')],
                      by = c('disease','geography','last_obs_time','h'), all.x = T, all.y = F)

rm(list=c('components_wide', 'weights_wide','output_list'))



##########################################
### Get Component Predictive Intervals ###
##########################################

if(!file.exists(paste0(savepath_new,'intervals_',outputname,'_',eval_type,".RDS"))){
ncores = 10
TO_RUN = results_arima[!duplicated(paste0(results_arima$disease,'_',results_arima$geography)),c('disease', 'geography')]
TO_RUN = data.frame(TO_RUN)
cl <- parallel::makeCluster(spec = ncores,type = sockettype)
setDefaultCluster(cl)
registerDoParallel(cl)
print(Sys.time())
myoutput <- foreach(g = 1:length(unlist(TO_RUN[,1])),
                    .verbose = F)%dopar%{
                      source(paste0(my_path,"/../process_data/epi_functions.R"))

                      RESULTS= results_obs[results_obs$disease == TO_RUN$disease[g] &
                                             results_obs$geography == TO_RUN$geography[g],]
                      RESULTS = RESULTS[order(RESULTS$x),]
                      RESULTS = data.frame(RESULTS)
                      fcst_indices <- unique(RESULTS$x)
                      fcst_indices = fcst_indices[fcst_indices>11]
                      library(dplyr)
                      intervaloutput = NULL
                      info_packet = list(ts = NULL, ts_scale = 'counts')
                      for(j in fcst_indices){
                        info_packet$ts=RESULTS$truth[1:j]
                        components = make_component_intervals(info_packet, h = h)
                        components$last_obs_time = j
                        components$disease = TO_RUN$disease[g]
                        components$geography = TO_RUN$geography[g]
                        intervaloutput = rbind(intervaloutput, components)
                      }
                      intervaloutput
                    }
stopCluster(cl)
print(Sys.time())
intervaloutput = do.call('rbind', myoutput)
saveRDS(intervaloutput, file = paste0(savepath_new,'intervals_',outputname,'_',eval_type,".RDS"))
}






##############################
##############################
### Get Forecast Intervals ###
##############################
##############################

TRUNC = 50
input = features_wide[,c(FEATURE_LIST, 'disease', 'geography', 'last_obs_time','fcst_epifforma','fcst_equalwt', 'truth', 'weight_sd', 'component_sd')]

intervaloutput = readRDS(paste0(savepath_new,'intervals_',outputname,'_',eval_type,".RDS"))
intervaloutput$h = rep(1:4, nrow(intervaloutput)/4)
input = merge(input, intervaloutput,
              by = c('disease', 'geography', 'last_obs_time', 'h'), all.x = T, all.y = F)
INTERVAL_LIST = names(input)
INTERVAL_LIST = INTERVAL_LIST[INTERVAL_LIST %in% c('arima', 'theta', 'gam')]
input = data.frame(input)

input$fcst_mult = 4*input$fcst_epifforma
mod_interval = readRDS(paste0(savemodelpath,"interval2.RDS"))
pred_wts_list <- list()
for(j in 1:length(mod_interval)){
  pred_wts_list[[j]] <- predict(mod_interval[[j]], newdata = as.matrix(input[,c(FEATURE_LIST, 'weight_sd', 'component_sd')]))
}
components = data.frame(input)[,c(INTERVAL_LIST, 'fcst_mult')]
pred_wts <- apply(simplify2array(pred_wts_list), 1:2, mean)
pred_wts[pred_wts<0.01] = 0 
pred_wts = t(apply(pred_wts, 1, FUN = function(x){x/sum(x)}))
input$pred_interval_epifforma = rowSums(pred_wts*as.matrix(components))
input$pred_interval_epifforma[input$pred_interval_epifforma>(2*TRUNC*input$fcst_epifforma) & !is.na(input$pred_interval_epifforma)] = 2*TRUNC*input$fcst_epifforma[input$pred_interval_epifforma>(2*TRUNC*input$fcst_epifforma) & !is.na(input$pred_interval_epifforma)]
input$pred_interval_epifforma = pmax(1e-9,input$pred_interval_epifforma) 



input$fcst_mult = 4*input$fcst_equalwt
mod_interval = readRDS(paste0(savemodelpath,"interval2_equal_wt.RDS"))
pred_wts_list <- list()
for(j in 1:length(mod_interval)){
  pred_wts_list[[j]] <- predict(mod_interval[[j]], newdata = as.matrix(input[,c(FEATURE_LIST)]))
}
components = data.frame(input)[,c(INTERVAL_LIST, 'fcst_mult')]
pred_wts <- apply(simplify2array(pred_wts_list), 1:2, mean)
pred_wts[pred_wts<0.01] = 0 
pred_wts = t(apply(pred_wts, 1, FUN = function(x){x/sum(x)}))
input$pred_interval_equalwt = rowSums(pred_wts*as.matrix(components))
input$pred_interval_equalwt[input$pred_interval_equalwt>(2*TRUNC*input$fcst_equalwt) & !is.na(input$pred_interval_equalwt)] = 2*TRUNC*input$fcst_equalwt[input$pred_interval_equalwt>(2*TRUNC*input$fcst_equalwt) & !is.na(input$pred_interval_equalwt)]
input$pred_interval_equalwt = pmax(1e-9,input$pred_interval_equalwt) 




################
### Evaluate ###
################

results = data.frame(input)[,c('disease','geography', 'last_obs_time','h','fcst_epifforma','fcst_equalwt','truth')]
results$x = results$last_obs_time + results$h

results$fcst_epifforma <- pmax(1e-10, results$fcst_epifforma)
results$fcst_equalwt <- pmax(1e-10, results$fcst_equalwt)


results$fcst_interval_epifforma_0.975 = results$fcst_epifforma+0.5*input$pred_interval_epifforma
results$fcst_interval_epifforma_0.025 = pmax(0,results$fcst_epifforma-0.5*input$pred_interval_epifforma)

results$fcst_interval_equalwt_0.975 = results$fcst_equalwt+0.5*input$pred_interval_equalwt
results$fcst_interval_equalwt_0.025 = pmax(0,results$fcst_equalwt-0.5*input$pred_interval_equalwt)



results$covers_interval_epifforma = as.numeric(results$truth <= results$fcst_interval_epifforma_0.975 & results$truth >= results$fcst_interval_epifforma_0.025)
results$width_interval_epifforma = results$fcst_interval_epifforma_0.975 - results$fcst_interval_epifforma_0.025

results$covers_interval_equalwt = as.numeric(results$truth <= results$fcst_interval_equalwt_0.975 & results$truth >= results$fcst_interval_equalwt_0.025)
results$width_interval_equalwt = results$fcst_interval_equalwt_0.975 - results$fcst_interval_equalwt_0.025



results = merge(results, results_rw[,c('geography', 'disease', 'last_obs_time','error_rw')], by = c('geography', 'disease', 'last_obs_time'), all.x = T, all.y = F)

interval_score = function(uppers, lowers, alpha_grid, y){
  (uppers-lowers)+(2/alpha_grid)*(lowers - y)*as.numeric(y<lowers)+
    (2/alpha_grid)*(y - uppers)*as.numeric(y>uppers)
}




results$wis_interval_epifforma = (1/1.5)*(0.5*abs(results$fcst_epifforma-results$truth)+ (0.05/2)*interval_score(uppers = results$fcst_interval_epifforma_0.975, lowers = results$fcst_interval_epifforma_0.025, alpha_grid = 0.05, y = results$truth))
results$wis_interval_equalwt = (1/1.5)*(0.5*abs(results$fcst_equalwt-results$truth)+ (0.05/2)*interval_score(uppers = results$fcst_interval_equalwt_0.975, lowers = results$fcst_interval_equalwt_0.025, alpha_grid = 0.05, y = results$truth))







results_long = rbind(data.frame(type = 'Epifforma (Interval)',             wis = results$wis_interval_epifforma,  width = results$width_interval_epifforma, coverage = results$covers_interval_epifforma, lower = results$fcst_interval_epifforma_0.025, upper = results$fcst_interval_epifforma_0.975, results[,c('geography', 'disease','last_obs_time', 'h','truth')], fcst = results$fcst_epifforma),
                     data.frame(type = 'Equalwt (Interval)',               wis = results$wis_interval_equalwt,   width = results$width_interval_equalwt, coverage = results$covers_interval_equalwt, lower = results$fcst_interval_equalwt_0.025, upper = results$fcst_interval_equalwt_0.975, results[,c('geography', 'disease','last_obs_time', 'h','truth')], fcst = results$fcst_equalwt))

library(plyr)
accuracydf <- ddply(results_long,.(type),summarise,
                    n = length(width),
                    median_width_scaled = median((width/fcst)[fcst>0],na.rm=T),
                    coverage = mean(coverage, na.rm=T), 
                    median_wis = median(wis[fcst>0], na.rm=T))
print(accuracydf)
gc()  
results_big = list(results = results_long, 
                   accuracydf = accuracydf)
saveRDS(results_big, file = paste0(savepath_new,outputname,'_',eval_type,".RDS"))




