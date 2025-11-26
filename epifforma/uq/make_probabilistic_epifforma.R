###################################
###################################
### epiFFORMA UQ Model Training ###
###################################
###################################

library(ggplot2)
library(data.table)
library(GGally)
library(viridis)
library(ggrepel)
library(dplyr)
theme_set(theme_classic())

# library(this.path)
my_path = here::here("epifforma", "uq") #this.path::here()
# setwd(my_path)
savepath <- paste0(my_path,"/../evaluate_model/evaluation/")
savepath_new <- paste0(my_path,"/../evaluate_probabilistic/evaluation/")
savemodelpath <- paste0(my_path, '/features2weights/') 
savevalidpath <- paste0(my_path, '/../uq/')
if(!dir.exists(savemodelpath)){dir.create(savemodelpath)}

## get .RDS file names
FILES = list.files(savepath)
FILES = FILES[grepl('synthetic_test_order_multierror.RDS',FILES)]

output_list = readRDS(paste0(savepath, FILES))

### Epifforma Predictions in Test Data
results = output_list['plot_df']$plot_df#[,-1]
results = results[results$type == 'fcst',]
results = results[results$model == 'epifforma',c('ts_length','ts_id','h','fcst','truth')]
results = results[!duplicated(paste0(results$h,'_', results$ts_length, '_', results$ts_id)),]

### Features in Test Data
FEATURE_LIST = setdiff(names(output_list['features']$features),c('ts_length','ts_id'))
features_wide = output_list['features']$features
features_wide = features_wide[!duplicated(paste0(features_wide$h,'_', features_wide$ts_length, '_', features_wide$ts_id)),]
features_wide = merge(features_wide,results, by = c('ts_length','ts_id','h'), all.x = T, all.y = F)
features_wide$seasonality[is.na(features_wide$seasonality)] = 0


### Components in Test Data
COMPONENT_LIST = setdiff(names(output_list['components']$components),c('ts_length','ts_id','h'))
components_wide = output_list['components']$components
components_wide$component_sd = apply(output_list['components']$components[,COMPONENT_LIST],1,sd)
weights_wide = output_list['weights']$weights
weights_wide$weight_sd = apply(output_list['weights']$weights[,paste0(COMPONENT_LIST,'_avg') ],1,sd)

features_wide = merge(features_wide, components_wide[,c('component_sd','ts_length','ts_id','h')],
                      by = c('ts_length','ts_id','h'), all.x = T, all.y = F)
features_wide$component_sd = features_wide$component_sd/features_wide$fcst
features_wide$component_sd[features_wide$component_sd>1 | is.infinite(features_wide$component_sd)]=1
features_wide = merge(features_wide, weights_wide[,c('weight_sd','ts_length','ts_id','h')],
                      by = c('ts_length','ts_id','h'), all.x = T, all.y = F)

rm(list=c('components_wide', 'weights_wide'))


### Interval Components in Test Data
lf <- list.files(paste0(savevalidpath,"validation_data_packets/"))[grep("interval",list.files(paste0(savevalidpath,"validation_data_packets/")))]
interval_data <- NULL
for(i in 1:length(lf)){
  interval_data <- rbind(interval_data,fread(paste0(paste0(savevalidpath,"validation_data_packets/"),lf[i])))
}
interval_data = interval_data[!duplicated(paste0(interval_data$ts_length, '_', interval_data$ts_id, '_', interval_data$h)),]
if('moa' %in% names(interval_data)){
  interval_data = subset(interval_data, select = -c(moa))
}

INTERVAL_LIST = names(interval_data)[!(names(interval_data) %in% c('ts_length', 'ts_id', 'h','truth'))]
interval_data = data.frame(interval_data)
for(i in 1:length(INTERVAL_LIST)){
  interval_data[interval_data[,INTERVAL_LIST[i]]<0.01,INTERVAL_LIST[i]]=0.01
}


features_wide = merge(features_wide,subset(interval_data, select = -c(truth)), by = c('ts_length','ts_id','h'), all.x = T, all.y = F)

rm(output_list)
rm(results)

input = features_wide[,c(FEATURE_LIST, INTERVAL_LIST,'fcst','truth', 'ts_id','ts_length', 'weight_sd', 'component_sd')]
input$fcst_mult = 4*input$fcst
input$entropy[is.na(input$entropy)] = 0 

set.seed(1234)
TRAIN = sample(0:1, size = nrow(input), replace = T, prob = c(0.25,0.75))
h = 4

NMODELS = 5



##########################
### Interval epiFFORMA ###
##########################

interval_long <- reshape2::melt(subset(interval_data, select = -c(truth)), id.vars = c("h", 'ts_id', 'ts_length'))
names(interval_long)[names(interval_long) == 'variable'] = 'model'
interval_long = merge(interval_long, input[,c(FEATURE_LIST, 'ts_length', 'ts_id','truth', 'fcst')], 
                      by = c( 'ts_length', 'ts_id', 'h'), all.x = T, all.y = F)

SUBSET = interval_long[interval_long$model == 'arima',] 
SUBSET$model = 'fcst_mult'
SUBSET$value = SUBSET$fcst*4
interval_long = rbind(interval_long, data.frame(SUBSET)[,colnames(interval_long)])

interval_long$upper = interval_long$fcst + 0.5*interval_long$value
interval_long$lower = pmax(0,interval_long$fcst - 0.5*interval_long$value)
interval_score = function(uppers, lowers, alpha_grid, y){
  (uppers-lowers)+(2/alpha_grid)*(lowers - y)*as.numeric(y<lowers)+
    (2/alpha_grid)*(y - uppers)*as.numeric(y>uppers)
}
###interval score
interval_long$is = interval_score(uppers = interval_long$upper, lowers = interval_long$lower, alpha_grid = 0.05, y = interval_long$truth)
my_func = function(is){
  if(min(is)==0){
    return(ifelse(is==0, 1/sum(is==0),0))
  }else{
    s1 <- sum(is)/is
    return(s1/sum(s1))
  }
}

interval_long = interval_long %>% dplyr::group_by(ts_length, ts_id, h) %>% dplyr::mutate(class_wt = my_func(is))
interval_long = interval_long[interval_long$class_wt > 0.99/length(unique(interval_long$model)),]
interval_long$class = as.numeric(factor(interval_long$model, levels = c(INTERVAL_LIST,'fcst_mult')))
interval_long = subset(interval_long, select=-c(model, truth, is, upper, lower, fcst, value))

interval_long = merge(interval_long, input[,c('component_sd', 'weight_sd', 'ts_length', 'ts_id', 'h')],
                      by = c('ts_length', 'ts_id', 'h'), all.x = T, all.y = F)

### Training Model
TRAIN_INTERVAL = as.numeric(paste0(interval_long$ts_length, '_', interval_long$ts_id, '_', interval_long$h) %in% paste0(input$ts_length, '_', input$ts_id, '_', input$h)[TRAIN == 1])
source(paste0(my_path, '/../process_data/epi_functions.R'))
mod_interval <- fit_lgbm_wt(input =data.frame(interval_long)[TRAIN_INTERVAL==1,c(FEATURE_LIST, 'weight_sd', 'component_sd', 'class', 'class_wt')], nmodels = NMODELS)
saveRDS(mod_interval, file = paste0(savemodelpath,"interval2.RDS"))

