############################
############################
### Get All UQ Summaries ###
############################
############################




library(data.table)
library(GGally)
lut_path = '~/GitLab/frodo/lookup_tables/'


interval_score = function(uppers, lowers, alpha_grid, y){
  (uppers-lowers)+(2/alpha_grid)*(lowers - y)*as.numeric(y<lowers)+
    (2/alpha_grid)*(y - uppers)*as.numeric(y>uppers)
}



dispersion_score = function(uppers, lowers, alpha_grid, y){
  (uppers-lowers)
}

underprediction_score = function(uppers, lowers, alpha_grid, y){
    (2/alpha_grid)*(y - uppers)*as.numeric(y>uppers)
}
overprediction_score = function(uppers, lowers, alpha_grid, y){
  (2/alpha_grid)*(lowers - y)*as.numeric(y<lowers)
}


#############
### GBM ###
#############

### results
my_path = '~/GitLab/frodo/summaries_uq/results/'
eval_path = '~/GitLab/frodo/uq/frodo/all/'

### modespecific ###
eval_path2 = '~/GitLab/frodo/uq/frodo/modespecific_cold/'

### diseasespecific ###
eval_path3 = '~/GitLab/frodo/uq/frodo/diseasespecific_cold/'

### sourcespecific ###
eval_path4 = '~/GitLab/frodo/uq/frodo/sourcespecific_cold/'


FILES = list.files(paste0(eval_path4))

RESULTS = data.frame(FILES = FILES)
RESULTS$geography = NA
RESULTS$disease = NA
RESULTS$disease_source = NA
RESULTS$N = NA
RESULTS$wis_1 = NA
RESULTS$wis_2 = NA
RESULTS$wis_3 = NA
RESULTS$wis_4 = NA

RESULTS$wis_dispersion_1 = NA
RESULTS$wis_dispersion_2 = NA
RESULTS$wis_dispersion_3 = NA
RESULTS$wis_dispersion_4 = NA

RESULTS$wis_overprediction_1 = NA
RESULTS$wis_overprediction_2 = NA
RESULTS$wis_overprediction_3 = NA
RESULTS$wis_overprediction_4 = NA

RESULTS$wis_underprediction_1 = NA
RESULTS$wis_underprediction_2 = NA
RESULTS$wis_underprediction_3 = NA
RESULTS$wis_underprediction_4 = NA

RESULTS$coverage50_1 = NA
RESULTS$coverage50_2 = NA
RESULTS$coverage50_3 = NA
RESULTS$coverage50_4 = NA
RESULTS$coverage80_1 = NA
RESULTS$coverage80_2 = NA
RESULTS$coverage80_3 = NA
RESULTS$coverage80_4 = NA
RESULTS$coverage95_1 = NA
RESULTS$coverage95_2 = NA
RESULTS$coverage95_3 = NA
RESULTS$coverage95_4 = NA

RESULTS$width50_1 = NA
RESULTS$width50_2 = NA
RESULTS$width50_3 = NA
RESULTS$width50_4 = NA
RESULTS$width80_1 = NA
RESULTS$width80_2 = NA
RESULTS$width80_3 = NA
RESULTS$width80_4 = NA
RESULTS$width95_1 = NA
RESULTS$width95_2 = NA
RESULTS$width95_3 = NA
RESULTS$width95_4 = NA


for(i in 1:length(FILES)){
  if(file.exists(paste0(eval_path4, FILES[i]))){
    output = try(data.frame(data.table::fread(paste0(eval_path, FILES[i]))), silent = T)
    output2 = try(data.frame(data.table::fread(paste0(eval_path2, FILES[i]))), silent = T)
    output3 = try(data.frame(data.table::fread(paste0(eval_path3, FILES[i]))), silent = T)
    output4 = try(data.frame(data.table::fread(paste0(eval_path4, FILES[i]))), silent = T)
    if(grepl('COVID',FILES[i])){output3 = output4} #didn't bother running disease-specific for COVID
    if(grepl('Dengue',FILES[i])){output3 = output4} #didn't bother running disease-specific for Dengue
    if(class(output)[1] != 'try-error' & class(output2)[1] != 'try-error' & class(output3)[1] != 'try-error' & class(output4)[1] != 'try-error'){
      
      
      output$fcst = pmax(0,output$fcst)
      output$fcst[!is.na(output$obs) & output$obs == 0] = 0
      output = output[output$obs > 0,]
      output = output[order(output$quant, output$last_obs_time, output$h),]
      output = output %>% dplyr::group_by(last_obs_time, h) %>% dplyr::mutate(value = sort(value))
      
      
      output2$fcst = pmax(0,output2$fcst)
      output2$fcst[!is.na(output2$obs) & output2$obs == 0] = 0
      output2 = output2[output2$obs > 0,]
      output2 = output2[order(output2$quant, output2$last_obs_time, output2$h),]
      output2 = output2 %>% dplyr::group_by(last_obs_time, h) %>% dplyr::mutate(value = sort(value))
      
      output3$fcst = pmax(0,output3$fcst)
      output3$fcst[!is.na(output3$obs) & output3$obs == 0] = 0
      output3 = output3[output3$obs > 0,] 
      output3 = output3[order(output3$quant, output3$last_obs_time, output3$h),]
      output3 = output3 %>% dplyr::group_by(last_obs_time, h) %>% dplyr::mutate(value = sort(value))
      
      output4$fcst = pmax(0,output4$fcst)
      output4$fcst[!is.na(output4$obs) & output4$obs == 0] = 0
      output4 = output4[output4$obs > 0,] 
      output4 = output4[order(output4$quant, output4$last_obs_time, output4$h),]
      output4 = output4 %>% dplyr::group_by(last_obs_time, h) %>% dplyr::mutate(value = sort(value))

      output2$fcst2 = output2$fcst
      output2$value2 = output2$value
      output3$fcst3 = output3$fcst
      output3$value3 = output3$value
      output4$fcst4 = output4$fcst
      output4$value4 = output4$value

      ####Note this is already satisfied for Frodo
      output$value[output$value > 3*output$obs & output$obs != 0] = 3*output$obs[output$value > 3*output$obs & output$obs != 0]
      output2$value2[output2$value2 > 3*output2$obs & output2$obs != 0] = 3*output2$obs[output2$value2 > 3*output2$obs & output2$obs != 0]
      output3$value3[output3$value3 > 3*output3$obs & output3$obs != 0] = 3*output3$obs[output3$value3 > 3*output3$obs & output3$obs != 0]
      output4$value4[output4$value4 > 3*output4$obs & output4$obs != 0] = 3*output4$obs[output4$value4 > 3*output4$obs & output4$obs != 0]


      #merge by date, last_obs_time got messed up somehow
      output = merge(output, output2[,c('date', 'h', 'fcst2','value2','quant')], by = c('date','h','quant'), all.x = T, all.y = T)
      output = merge(output, output3[,c('date', 'h', 'fcst3','value3','quant')], by = c('date','h','quant'), all.x = T, all.y = T)
      output = merge(output, output4[,c('date', 'h', 'fcst4','value4','quant')], by = c('date','h','quant'), all.x = T, all.y = T)
      output = output[!is.na(output$fcst) & !is.na(output$fcst2) & !is.na(output$fcst3) & !is.na(output$fcst4),]
      output = output[!is.na(output$value) & !is.na(output$value2) & !is.na(output$value3) & !is.na(output$value4),]
   
      output = output[order(output$date, output$h, output$quant),]


      WIS_1 = (abs(output$fcst[output$quant == 0.5] - output$truth[output$quant == 0.5])/2)+
        (0.5/2)*interval_score(uppers = output$value[output$quant == 0.75], lowers = output$value[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*interval_score(uppers = output$value[output$quant == 0.9], lowers = output$value[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*interval_score(uppers = output$value[output$quant == 0.975], lowers = output$value[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_2 = (abs(output$fcst2[output$quant == 0.5] - output$truth[output$quant == 0.5])/2)+
        (0.5/2)*interval_score(uppers = output$value2[output$quant == 0.75], lowers = output$value2[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*interval_score(uppers = output$value2[output$quant == 0.9], lowers = output$value2[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*interval_score(uppers = output$value2[output$quant == 0.975], lowers = output$value2[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_3 = (abs(output$fcst3[output$quant == 0.5] - output$truth[output$quant == 0.5])/2)+
        (0.5/2)*interval_score(uppers = output$value3[output$quant == 0.75], lowers = output$value3[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*interval_score(uppers = output$value3[output$quant == 0.9], lowers = output$value3[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*interval_score(uppers = output$value3[output$quant == 0.975], lowers = output$value3[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_4 = (abs(output$fcst4[output$quant == 0.5] - output$truth[output$quant == 0.5])/2)+
        (0.5/2)*interval_score(uppers = output$value4[output$quant == 0.75], lowers = output$value4[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*interval_score(uppers = output$value4[output$quant == 0.9], lowers = output$value4[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*interval_score(uppers = output$value4[output$quant == 0.975], lowers = output$value4[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      

      
      WIS_DISPERSION_1 = (0.5/2)*dispersion_score(uppers = output$value[output$quant == 0.75], lowers = output$value[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*dispersion_score(uppers = output$value[output$quant == 0.9], lowers = output$value[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*dispersion_score(uppers = output$value[output$quant == 0.975], lowers = output$value[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_DISPERSION_2 = (0.5/2)*dispersion_score(uppers = output$value2[output$quant == 0.75], lowers = output$value2[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*dispersion_score(uppers = output$value2[output$quant == 0.9], lowers = output$value2[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*dispersion_score(uppers = output$value2[output$quant == 0.975], lowers = output$value2[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_DISPERSION_3 = (0.5/2)*dispersion_score(uppers = output$value3[output$quant == 0.75], lowers = output$value3[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*dispersion_score(uppers = output$value3[output$quant == 0.9], lowers = output$value3[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*dispersion_score(uppers = output$value3[output$quant == 0.975], lowers = output$value3[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_DISPERSION_4 = (0.5/2)*dispersion_score(uppers = output$value4[output$quant == 0.75], lowers = output$value4[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*dispersion_score(uppers = output$value4[output$quant == 0.9], lowers = output$value4[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*dispersion_score(uppers = output$value4[output$quant == 0.975], lowers = output$value4[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      

      WIS_OVERPREDICTION_1 = (0.5/2)*overprediction_score(uppers = output$value[output$quant == 0.75], lowers = output$value[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*overprediction_score(uppers = output$value[output$quant == 0.9], lowers = output$value[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*overprediction_score(uppers = output$value[output$quant == 0.975], lowers = output$value[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_OVERPREDICTION_2 = (0.5/2)*overprediction_score(uppers = output$value2[output$quant == 0.75], lowers = output$value2[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*overprediction_score(uppers = output$value2[output$quant == 0.9], lowers = output$value2[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*overprediction_score(uppers = output$value2[output$quant == 0.975], lowers = output$value2[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_OVERPREDICTION_3 = (0.5/2)*overprediction_score(uppers = output$value3[output$quant == 0.75], lowers = output$value3[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*overprediction_score(uppers = output$value3[output$quant == 0.9], lowers = output$value3[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*overprediction_score(uppers = output$value3[output$quant == 0.975], lowers = output$value3[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_OVERPREDICTION_4 = (0.5/2)*overprediction_score(uppers = output$value4[output$quant == 0.75], lowers = output$value4[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*overprediction_score(uppers = output$value4[output$quant == 0.9], lowers = output$value4[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*overprediction_score(uppers = output$value4[output$quant == 0.975], lowers = output$value4[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      

      WIS_UNDERPREDICTION_1 = (0.5/2)*underprediction_score(uppers = output$value[output$quant == 0.75], lowers = output$value[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*underprediction_score(uppers = output$value[output$quant == 0.9], lowers = output$value[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*underprediction_score(uppers = output$value[output$quant == 0.975], lowers = output$value[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_UNDERPREDICTION_2 = (0.5/2)*underprediction_score(uppers = output$value2[output$quant == 0.75], lowers = output$value2[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*underprediction_score(uppers = output$value2[output$quant == 0.9], lowers = output$value2[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*underprediction_score(uppers = output$value2[output$quant == 0.975], lowers = output$value2[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_UNDERPREDICTION_3 = (0.5/2)*underprediction_score(uppers = output$value3[output$quant == 0.75], lowers = output$value3[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*underprediction_score(uppers = output$value3[output$quant == 0.9], lowers = output$value3[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*underprediction_score(uppers = output$value3[output$quant == 0.975], lowers = output$value3[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_UNDERPREDICTION_4 = (0.5/2)*underprediction_score(uppers = output$value4[output$quant == 0.75], lowers = output$value4[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*underprediction_score(uppers = output$value4[output$quant == 0.9], lowers = output$value4[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*underprediction_score(uppers = output$value4[output$quant == 0.975], lowers = output$value4[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])

      RESULTS$N[i] = nrow(output)
      RESULTS$wis_1[i] = (1/(3+0.5))*mean(WIS_1)
      RESULTS$wis_2[i] =  (1/(3+0.5))*mean(WIS_2)
      RESULTS$wis_3[i] =  (1/(3+0.5))*mean(WIS_3)
      RESULTS$wis_4[i] =  (1/(3+0.5))*mean(WIS_4)

      RESULTS$wis_dispersion_1[i] = (1/(3+0.5))*mean(WIS_DISPERSION_1)
      RESULTS$wis_dispersion_2[i] = (1/(3+0.5))*mean(WIS_DISPERSION_2)
      RESULTS$wis_dispersion_3[i] = (1/(3+0.5))*mean(WIS_DISPERSION_3)
      RESULTS$wis_dispersion_4[i] = (1/(3+0.5))*mean(WIS_DISPERSION_4)

      RESULTS$wis_overprediction_1[i] = (1/(3+0.5))*mean(WIS_OVERPREDICTION_1)
      RESULTS$wis_overprediction_2[i] = (1/(3+0.5))*mean(WIS_OVERPREDICTION_2)
      RESULTS$wis_overprediction_3[i] = (1/(3+0.5))*mean(WIS_OVERPREDICTION_3)
      RESULTS$wis_overprediction_4[i] = (1/(3+0.5))*mean(WIS_OVERPREDICTION_4)
      
      RESULTS$wis_underprediction_1[i] = (1/(3+0.5))*mean(WIS_UNDERPREDICTION_1)
      RESULTS$wis_underprediction_2[i] = (1/(3+0.5))*mean(WIS_UNDERPREDICTION_2)
      RESULTS$wis_underprediction_3[i] = (1/(3+0.5))*mean(WIS_UNDERPREDICTION_3)
      RESULTS$wis_underprediction_4[i] = (1/(3+0.5))*mean(WIS_UNDERPREDICTION_4)

      RESULTS$coverage50_1[i] = mean(output$truth[output$quant == 0.5] <= output$value[output$quant == 0.75] & output$truth[output$quant == 0.5] >= output$value[output$quant == 0.25], na.rm=T)
      RESULTS$coverage50_2[i] = mean(output$truth[output$quant == 0.5] <= output$value2[output$quant == 0.75] & output$truth[output$quant == 0.5] >= output$value2[output$quant == 0.25], na.rm=T)
      RESULTS$coverage50_3[i] = mean(output$truth[output$quant == 0.5] <= output$value3[output$quant == 0.75] & output$truth[output$quant == 0.5] >= output$value3[output$quant == 0.25], na.rm=T)
      RESULTS$coverage50_4[i] = mean(output$truth[output$quant == 0.5] <= output$value4[output$quant == 0.75] & output$truth[output$quant == 0.5] >= output$value4[output$quant == 0.25], na.rm=T)

      RESULTS$coverage80_1[i] = mean(output$truth[output$quant == 0.5] <= output$value[output$quant == 0.9] & output$truth[output$quant == 0.5] >= output$value[output$quant == 0.1], na.rm=T)
      RESULTS$coverage80_2[i] = mean(output$truth[output$quant == 0.5] <= output$value2[output$quant == 0.9] & output$truth[output$quant == 0.5] >= output$value2[output$quant == 0.1], na.rm=T)
      RESULTS$coverage80_3[i] = mean(output$truth[output$quant == 0.5] <= output$value3[output$quant == 0.9] & output$truth[output$quant == 0.5] >= output$value3[output$quant == 0.1], na.rm=T)
      RESULTS$coverage80_4[i] = mean(output$truth[output$quant == 0.5] <= output$value4[output$quant == 0.9] & output$truth[output$quant == 0.5] >= output$value4[output$quant == 0.1], na.rm=T)

      RESULTS$coverage95_1[i] = mean(output$truth[output$quant == 0.5] <= output$value[output$quant == 0.975] & output$truth[output$quant == 0.5] >= output$value[output$quant == 0.025], na.rm=T)
      RESULTS$coverage95_2[i] = mean(output$truth[output$quant == 0.5] <= output$value2[output$quant == 0.975] & output$truth[output$quant == 0.5] >= output$value2[output$quant == 0.025], na.rm=T)
      RESULTS$coverage95_3[i] = mean(output$truth[output$quant == 0.5] <= output$value3[output$quant == 0.975] & output$truth[output$quant == 0.5] >= output$value3[output$quant == 0.025], na.rm=T)
      RESULTS$coverage95_4[i] = mean(output$truth[output$quant == 0.5] <= output$value4[output$quant == 0.975] & output$truth[output$quant == 0.5] >= output$value4[output$quant == 0.025], na.rm=T)

      
      RESULTS$width50_1[i] = mean((output$value[output$quant == 0.75]-output$value[output$quant == 0.25])/output$truth, na.rm=T)
      RESULTS$width50_2[i] = mean((output$value2[output$quant == 0.75]-output$value2[output$quant == 0.25])/output$truth, na.rm=T)
      RESULTS$width50_3[i] = mean((output$value3[output$quant == 0.75]-output$value3[output$quant == 0.25])/output$truth, na.rm=T)
      RESULTS$width50_4[i] = mean((output$value4[output$quant == 0.75]-output$value4[output$quant == 0.25])/output$truth, na.rm=T)

      RESULTS$width80_1[i] = mean((output$value[output$quant == 0.9]-output$value[output$quant == 0.1])/output$truth, na.rm=T)
      RESULTS$width80_2[i] = mean((output$value2[output$quant == 0.9]-output$value2[output$quant == 0.1])/output$truth, na.rm=T)
      RESULTS$width80_3[i] = mean((output$value3[output$quant == 0.9]-output$value3[output$quant == 0.1])/output$truth, na.rm=T)
      RESULTS$width80_4[i] = mean((output$value4[output$quant == 0.9]-output$value4[output$quant == 0.1])/output$truth, na.rm=T)

      RESULTS$width95_1[i] = mean((output$value[output$quant == 0.975]-output$value[output$quant == 0.025])/output$truth, na.rm=T)
      RESULTS$width95_2[i] = mean((output$value2[output$quant == 0.975]-output$value2[output$quant == 0.025])/output$truth, na.rm=T)
      RESULTS$width95_3[i] = mean((output$value3[output$quant == 0.975]-output$value3[output$quant == 0.025])/output$truth, na.rm=T)
      RESULTS$width95_4[i] = mean((output$value4[output$quant == 0.975]-output$value4[output$quant == 0.025])/output$truth, na.rm=T)

      
    }else{
      print(paste0('Skipped: ', FILES[i]))
    }
  }
}


RESULTS_ALL = RESULTS[RESULTS$N > 0,]

RESULTS_ALL$disease_source = gsub('real_uq_mat_','',RESULTS_ALL$FILES)
RESULTS_ALL$disease_source = gsub('.csv','',RESULTS_ALL$disease_source)
RESULTS_ALL$disease_source = gsub('(.*)_\\w+', '\\1', RESULTS_ALL$disease_source)
RESULTS_ALL$disease_source = gsub('[0-9]', '', RESULTS_ALL$disease_source)
RESULTS_ALL$disease_source = substr(RESULTS_ALL$disease_source, 2, nchar(RESULTS_ALL$disease_source))
RESULTS_ALL$disease_source = substr(RESULTS_ALL$disease_source, 1, nchar(RESULTS_ALL$disease_source)-1)

RESULTS_ALL$plot_id = RESULTS_ALL$disease_source

write.csv(RESULTS_ALL, file = paste0('~/GitLab/frodo/summaries_uq/tables/',"frodo.csv"), quote = F, row.names = F)


RESULTS_ALL = RESULTS_ALL %>% dplyr::group_by(plot_id) %>% dplyr::mutate(N_agg = sum(N))
RESULTS_AGG = RESULTS_ALL %>% dplyr::group_by(plot_id) %>% dplyr::mutate(wis_1_agg = sum(N*wis_1)/N_agg,
                                                                         wis_2_agg = sum(N*wis_2)/N_agg,
                                                                         wis_3_agg = sum(N*wis_3)/N_agg,
                                                                         wis_4_agg = sum(N*wis_4)/N_agg,
                                                                         wis_dispersion_1_agg = sum(N*wis_dispersion_1)/N_agg,
                                                                         wis_dispersion_2_agg = sum(N*wis_dispersion_2)/N_agg,
                                                                         wis_dispersion_3_agg = sum(N*wis_dispersion_3)/N_agg,
                                                                         wis_dispersion_4_agg = sum(N*wis_dispersion_4)/N_agg,
                                                                         wis_overprediction_1_agg = sum(N*wis_overprediction_1)/N_agg,
                                                                         wis_overprediction_2_agg = sum(N*wis_overprediction_2)/N_agg,
                                                                         wis_overprediction_3_agg = sum(N*wis_overprediction_3)/N_agg,
                                                                         wis_overprediction_4_agg = sum(N*wis_overprediction_4)/N_agg,
                                                                         wis_underprediction_1_agg = sum(N*wis_underprediction_1)/N_agg,
                                                                         wis_underprediction_2_agg = sum(N*wis_underprediction_2)/N_agg,
                                                                         wis_underprediction_3_agg = sum(N*wis_underprediction_3)/N_agg,
                                                                         wis_underprediction_4_agg = sum(N*wis_underprediction_4)/N_agg,
                                                                         coverage50_1_agg = sum(N*coverage50_1)/N_agg,
                                                                         coverage50_2_agg = sum(N*coverage50_2)/N_agg,
                                                                         coverage50_3_agg = sum(N*coverage50_3)/N_agg,
                                                                         coverage50_4_agg = sum(N*coverage50_4)/N_agg,
                                                                         coverage80_1_agg = sum(N*coverage80_1)/N_agg,
                                                                         coverage80_2_agg = sum(N*coverage80_2)/N_agg,
                                                                         coverage80_3_agg = sum(N*coverage80_3)/N_agg,
                                                                         coverage80_4_agg = sum(N*coverage80_4)/N_agg,
                                                                         coverage95_1_agg = sum(N*coverage95_1)/N_agg,
                                                                         coverage95_2_agg = sum(N*coverage95_2)/N_agg,
                                                                         coverage95_3_agg = sum(N*coverage95_3)/N_agg,
                                                                         coverage95_4_agg = sum(N*coverage95_4)/N_agg,
                                                                         width50_1_agg = sum(N*width50_1)/N_agg,
                                                                         width50_2_agg = sum(N*width50_2)/N_agg,
                                                                         width50_3_agg = sum(N*width50_3)/N_agg,
                                                                         width50_4_agg = sum(N*width50_4)/N_agg,
                                                                         width80_1_agg = sum(N*width80_1)/N_agg,
                                                                         width80_2_agg = sum(N*width80_2)/N_agg,
                                                                         width80_3_agg = sum(N*width80_3)/N_agg,
                                                                         width80_4_agg = sum(N*width80_4)/N_agg,
                                                                         width95_1_agg = sum(N*width95_1)/N_agg,
                                                                         width95_2_agg = sum(N*width95_2)/N_agg,
                                                                         width95_3_agg = sum(N*width95_3)/N_agg,
                                                                         width95_4_agg = sum(N*width95_4)/N_agg)

RESULTS_AGG = RESULTS_AGG[!duplicated(RESULTS_AGG$plot_id),]

write.csv(RESULTS_AGG, file = paste0('~/GitLab/frodo/summaries_uq/tables/',"frodo_agg.csv"), quote = F, row.names = F)





############
### LSTM ###
############


### results
my_path = '~/GitLab/frodo/summaries_uq/results/'
eval_path = '~/GitLab/frodo/uq/lstm/all/'

### modespecific ###
eval_path2 = '~/GitLab/frodo/uq/lstm/modespecific_cold/'

### diseasespecific ###
eval_path3 = '~/GitLab/frodo/uq/lstm/diseasespecific_cold/'

### sourcespecific ###
eval_path4 = '~/GitLab/frodo/uq/lstm/sourcespecific_cold/'




FILES = list.files(paste0(eval_path4))

RESULTS = data.frame(FILES = FILES)
RESULTS$geography = NA
RESULTS$disease = NA
RESULTS$disease_source = NA
RESULTS$N = NA
RESULTS$wis_1 = NA
RESULTS$wis_2 = NA
RESULTS$wis_3 = NA
RESULTS$wis_4 = NA

RESULTS$wis_dispersion_1 = NA
RESULTS$wis_dispersion_2 = NA
RESULTS$wis_dispersion_3 = NA
RESULTS$wis_dispersion_4 = NA

RESULTS$wis_overprediction_1 = NA
RESULTS$wis_overprediction_2 = NA
RESULTS$wis_overprediction_3 = NA
RESULTS$wis_overprediction_4 = NA

RESULTS$wis_underprediction_1 = NA
RESULTS$wis_underprediction_2 = NA
RESULTS$wis_underprediction_3 = NA
RESULTS$wis_underprediction_4 = NA

RESULTS$coverage50_1 = NA
RESULTS$coverage50_2 = NA
RESULTS$coverage50_3 = NA
RESULTS$coverage50_4 = NA
RESULTS$coverage80_1 = NA
RESULTS$coverage80_2 = NA
RESULTS$coverage80_3 = NA
RESULTS$coverage80_4 = NA
RESULTS$coverage95_1 = NA
RESULTS$coverage95_2 = NA
RESULTS$coverage95_3 = NA
RESULTS$coverage95_4 = NA

RESULTS$width50_1 = NA
RESULTS$width50_2 = NA
RESULTS$width50_3 = NA
RESULTS$width50_4 = NA
RESULTS$width80_1 = NA
RESULTS$width80_2 = NA
RESULTS$width80_3 = NA
RESULTS$width80_4 = NA
RESULTS$width95_1 = NA
RESULTS$width95_2 = NA
RESULTS$width95_3 = NA
RESULTS$width95_4 = NA

# i=which(grepl('_COVID_jhuowid_935',FILES_ALL))
# i=which(grepl('_COVID_jhuowid_935',FILES))
for(i in 1:length(FILES)){
  if(file.exists(paste0(eval_path4, FILES[i]))){
    output = try(data.frame(data.table::fread(paste0(eval_path, FILES[i]))), silent = T)
    output2 = try(data.frame(data.table::fread(paste0(eval_path2, FILES[i]))), silent = T)
    output3 = try(data.frame(data.table::fread(paste0(eval_path3, FILES[i]))), silent = T)
    output4 = try(data.frame(data.table::fread(paste0(eval_path4, FILES[i]))), silent = T)
    if(grepl('COVID',FILES[i])){output3 = output4} #didn't bother running disease-specific for COVID
    if(grepl('Dengue',FILES[i])){output3 = output4} #didn't bother running disease-specific for Dengue
    if(class(output)[1] != 'try-error' & class(output2)[1] != 'try-error' & class(output3)[1] != 'try-error' & class(output4)[1] != 'try-error'){
      
      
      output$fcst = pmax(0,output$fcst)
      output$fcst[!is.na(output$obs) & output$obs == 0] = 0
      output = output[output$obs > 0,]
      output = output[order(output$quant, output$last_obs_time, output$h),]
      output = output %>% dplyr::group_by(last_obs_time, h) %>% dplyr::mutate(value = sort(value))
      
      
      output2$fcst = pmax(0,output2$fcst)
      output2$fcst[!is.na(output2$obs) & output2$obs == 0] = 0
      output2 = output2[output2$obs > 0,]
      output2 = output2[order(output2$quant, output2$last_obs_time, output2$h),]
      output2 = output2 %>% dplyr::group_by(last_obs_time, h) %>% dplyr::mutate(value = sort(value))
      
      output3$fcst = pmax(0,output3$fcst)
      output3$fcst[!is.na(output3$obs) & output3$obs == 0] = 0
      output3 = output3[output3$obs > 0,] 
      output3 = output3[order(output3$quant, output3$last_obs_time, output3$h),]
      output3 = output3 %>% dplyr::group_by(last_obs_time, h) %>% dplyr::mutate(value = sort(value))
      
      output4$fcst = pmax(0,output4$fcst)
      output4$fcst[!is.na(output4$obs) & output4$obs == 0] = 0
      output4 = output4[output4$obs > 0,] 
      output4 = output4[order(output4$quant, output4$last_obs_time, output4$h),]
      output4 = output4 %>% dplyr::group_by(last_obs_time, h) %>% dplyr::mutate(value = sort(value))
      
      output2$fcst2 = output2$fcst
      output2$value2 = output2$value
      output3$fcst3 = output3$fcst
      output3$value3 = output3$value
      output4$fcst4 = output4$fcst
      output4$value4 = output4$value

      output$value[output$value > 3*output$obs & output$obs != 0] = 3*output$obs[output$value > 3*output$obs & output$obs != 0]
      output2$value2[output2$value2 > 3*output2$obs & output2$obs != 0] = 3*output2$obs[output2$value2 > 3*output2$obs & output2$obs != 0]
      output3$value3[output3$value3 > 3*output3$obs & output3$obs != 0] = 3*output3$obs[output3$value3 > 3*output3$obs & output3$obs != 0]
      output4$value4[output4$value4 > 3*output4$obs & output4$obs != 0] = 3*output4$obs[output4$value4 > 3*output4$obs & output4$obs != 0]
      
 
      #merge by date, last_obs_time got messed up somehow
      output = merge(output, output2[,c('date', 'h', 'fcst2','value2','quant')], by = c('date','h','quant'), all.x = T, all.y = T)
      output = merge(output, output3[,c('date', 'h', 'fcst3','value3','quant')], by = c('date','h','quant'), all.x = T, all.y = T)
      output = merge(output, output4[,c('date', 'h', 'fcst4','value4','quant')], by = c('date','h','quant'), all.x = T, all.y = T)
      output = output[!is.na(output$fcst) & !is.na(output$fcst2) & !is.na(output$fcst3) & !is.na(output$fcst4),]
      output = output[!is.na(output$value) & !is.na(output$value2) & !is.na(output$value3) & !is.na(output$value4),]
      
      output = output[order(output$date, output$h, output$quant),]
      
      
      WIS_1 = (abs(output$fcst[output$quant == 0.5] - output$truth[output$quant == 0.5])/2)+
        (0.5/2)*interval_score(uppers = output$value[output$quant == 0.75], lowers = output$value[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*interval_score(uppers = output$value[output$quant == 0.9], lowers = output$value[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*interval_score(uppers = output$value[output$quant == 0.975], lowers = output$value[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_2 = (abs(output$fcst2[output$quant == 0.5] - output$truth[output$quant == 0.5])/2)+
        (0.5/2)*interval_score(uppers = output$value2[output$quant == 0.75], lowers = output$value2[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*interval_score(uppers = output$value2[output$quant == 0.9], lowers = output$value2[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*interval_score(uppers = output$value2[output$quant == 0.975], lowers = output$value2[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_3 = (abs(output$fcst3[output$quant == 0.5] - output$truth[output$quant == 0.5])/2)+
        (0.5/2)*interval_score(uppers = output$value3[output$quant == 0.75], lowers = output$value3[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*interval_score(uppers = output$value3[output$quant == 0.9], lowers = output$value3[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*interval_score(uppers = output$value3[output$quant == 0.975], lowers = output$value3[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_4 = (abs(output$fcst4[output$quant == 0.5] - output$truth[output$quant == 0.5])/2)+
        (0.5/2)*interval_score(uppers = output$value4[output$quant == 0.75], lowers = output$value4[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*interval_score(uppers = output$value4[output$quant == 0.9], lowers = output$value4[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*interval_score(uppers = output$value4[output$quant == 0.975], lowers = output$value4[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])

      WIS_DISPERSION_1 = (0.5/2)*dispersion_score(uppers = output$value[output$quant == 0.75], lowers = output$value[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*dispersion_score(uppers = output$value[output$quant == 0.9], lowers = output$value[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*dispersion_score(uppers = output$value[output$quant == 0.975], lowers = output$value[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_DISPERSION_2 = (0.5/2)*dispersion_score(uppers = output$value2[output$quant == 0.75], lowers = output$value2[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*dispersion_score(uppers = output$value2[output$quant == 0.9], lowers = output$value2[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*dispersion_score(uppers = output$value2[output$quant == 0.975], lowers = output$value2[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_DISPERSION_3 = (0.5/2)*dispersion_score(uppers = output$value3[output$quant == 0.75], lowers = output$value3[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*dispersion_score(uppers = output$value3[output$quant == 0.9], lowers = output$value3[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*dispersion_score(uppers = output$value3[output$quant == 0.975], lowers = output$value3[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_DISPERSION_4 = (0.5/2)*dispersion_score(uppers = output$value4[output$quant == 0.75], lowers = output$value4[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*dispersion_score(uppers = output$value4[output$quant == 0.9], lowers = output$value4[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*dispersion_score(uppers = output$value4[output$quant == 0.975], lowers = output$value4[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_OVERPREDICTION_1 = (0.5/2)*overprediction_score(uppers = output$value[output$quant == 0.75], lowers = output$value[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*overprediction_score(uppers = output$value[output$quant == 0.9], lowers = output$value[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*overprediction_score(uppers = output$value[output$quant == 0.975], lowers = output$value[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_OVERPREDICTION_2 = (0.5/2)*overprediction_score(uppers = output$value2[output$quant == 0.75], lowers = output$value2[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*overprediction_score(uppers = output$value2[output$quant == 0.9], lowers = output$value2[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*overprediction_score(uppers = output$value2[output$quant == 0.975], lowers = output$value2[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_OVERPREDICTION_3 = (0.5/2)*overprediction_score(uppers = output$value3[output$quant == 0.75], lowers = output$value3[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*overprediction_score(uppers = output$value3[output$quant == 0.9], lowers = output$value3[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*overprediction_score(uppers = output$value3[output$quant == 0.975], lowers = output$value3[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_OVERPREDICTION_4 = (0.5/2)*overprediction_score(uppers = output$value4[output$quant == 0.75], lowers = output$value4[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*overprediction_score(uppers = output$value4[output$quant == 0.9], lowers = output$value4[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*overprediction_score(uppers = output$value4[output$quant == 0.975], lowers = output$value4[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])

      WIS_UNDERPREDICTION_1 = (0.5/2)*underprediction_score(uppers = output$value[output$quant == 0.75], lowers = output$value[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*underprediction_score(uppers = output$value[output$quant == 0.9], lowers = output$value[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*underprediction_score(uppers = output$value[output$quant == 0.975], lowers = output$value[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_UNDERPREDICTION_2 = (0.5/2)*underprediction_score(uppers = output$value2[output$quant == 0.75], lowers = output$value2[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*underprediction_score(uppers = output$value2[output$quant == 0.9], lowers = output$value2[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*underprediction_score(uppers = output$value2[output$quant == 0.975], lowers = output$value2[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_UNDERPREDICTION_3 = (0.5/2)*underprediction_score(uppers = output$value3[output$quant == 0.75], lowers = output$value3[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*underprediction_score(uppers = output$value3[output$quant == 0.9], lowers = output$value3[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*underprediction_score(uppers = output$value3[output$quant == 0.975], lowers = output$value3[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_UNDERPREDICTION_4 = (0.5/2)*underprediction_score(uppers = output$value4[output$quant == 0.75], lowers = output$value4[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*underprediction_score(uppers = output$value4[output$quant == 0.9], lowers = output$value4[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*underprediction_score(uppers = output$value4[output$quant == 0.975], lowers = output$value4[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])

      
      RESULTS$N[i] = nrow(output)
      RESULTS$wis_1[i] = (1/(3+0.5))*mean(WIS_1)
      RESULTS$wis_2[i] =  (1/(3+0.5))*mean(WIS_2)
      RESULTS$wis_3[i] =  (1/(3+0.5))*mean(WIS_3)
      RESULTS$wis_4[i] =  (1/(3+0.5))*mean(WIS_4)
      
      RESULTS$wis_dispersion_1[i] = (1/(3+0.5))*mean(WIS_DISPERSION_1)
      RESULTS$wis_dispersion_2[i] = (1/(3+0.5))*mean(WIS_DISPERSION_2)
      RESULTS$wis_dispersion_3[i] = (1/(3+0.5))*mean(WIS_DISPERSION_3)
      RESULTS$wis_dispersion_4[i] = (1/(3+0.5))*mean(WIS_DISPERSION_4)

      RESULTS$wis_overprediction_1[i] = (1/(3+0.5))*mean(WIS_OVERPREDICTION_1)
      RESULTS$wis_overprediction_2[i] = (1/(3+0.5))*mean(WIS_OVERPREDICTION_2)
      RESULTS$wis_overprediction_3[i] = (1/(3+0.5))*mean(WIS_OVERPREDICTION_3)
      RESULTS$wis_overprediction_4[i] = (1/(3+0.5))*mean(WIS_OVERPREDICTION_4)
      
      RESULTS$wis_underprediction_1[i] = (1/(3+0.5))*mean(WIS_UNDERPREDICTION_1)
      RESULTS$wis_underprediction_2[i] = (1/(3+0.5))*mean(WIS_UNDERPREDICTION_2)
      RESULTS$wis_underprediction_3[i] = (1/(3+0.5))*mean(WIS_UNDERPREDICTION_3)
      RESULTS$wis_underprediction_4[i] = (1/(3+0.5))*mean(WIS_UNDERPREDICTION_4)
      
      RESULTS$coverage50_1[i] = mean(output$truth[output$quant == 0.5] <= output$value[output$quant == 0.75] & output$truth[output$quant == 0.5] >= output$value[output$quant == 0.25], na.rm=T)
      RESULTS$coverage50_2[i] = mean(output$truth[output$quant == 0.5] <= output$value2[output$quant == 0.75] & output$truth[output$quant == 0.5] >= output$value2[output$quant == 0.25], na.rm=T)
      RESULTS$coverage50_3[i] = mean(output$truth[output$quant == 0.5] <= output$value3[output$quant == 0.75] & output$truth[output$quant == 0.5] >= output$value3[output$quant == 0.25], na.rm=T)
      RESULTS$coverage50_4[i] = mean(output$truth[output$quant == 0.5] <= output$value4[output$quant == 0.75] & output$truth[output$quant == 0.5] >= output$value4[output$quant == 0.25], na.rm=T)

      RESULTS$coverage80_1[i] = mean(output$truth[output$quant == 0.5] <= output$value[output$quant == 0.9] & output$truth[output$quant == 0.5] >= output$value[output$quant == 0.1], na.rm=T)
      RESULTS$coverage80_2[i] = mean(output$truth[output$quant == 0.5] <= output$value2[output$quant == 0.9] & output$truth[output$quant == 0.5] >= output$value2[output$quant == 0.1], na.rm=T)
      RESULTS$coverage80_3[i] = mean(output$truth[output$quant == 0.5] <= output$value3[output$quant == 0.9] & output$truth[output$quant == 0.5] >= output$value3[output$quant == 0.1], na.rm=T)
      RESULTS$coverage80_4[i] = mean(output$truth[output$quant == 0.5] <= output$value4[output$quant == 0.9] & output$truth[output$quant == 0.5] >= output$value4[output$quant == 0.1], na.rm=T)

      RESULTS$coverage95_1[i] = mean(output$truth[output$quant == 0.5] <= output$value[output$quant == 0.975] & output$truth[output$quant == 0.5] >= output$value[output$quant == 0.025], na.rm=T)
      RESULTS$coverage95_2[i] = mean(output$truth[output$quant == 0.5] <= output$value2[output$quant == 0.975] & output$truth[output$quant == 0.5] >= output$value2[output$quant == 0.025], na.rm=T)
      RESULTS$coverage95_3[i] = mean(output$truth[output$quant == 0.5] <= output$value3[output$quant == 0.975] & output$truth[output$quant == 0.5] >= output$value3[output$quant == 0.025], na.rm=T)
      RESULTS$coverage95_4[i] = mean(output$truth[output$quant == 0.5] <= output$value4[output$quant == 0.975] & output$truth[output$quant == 0.5] >= output$value4[output$quant == 0.025], na.rm=T)
      
      RESULTS$width50_1[i] = mean((output$value[output$quant == 0.75]-output$value[output$quant == 0.25])/output$truth, na.rm=T)
      RESULTS$width50_2[i] = mean((output$value2[output$quant == 0.75]-output$value2[output$quant == 0.25])/output$truth, na.rm=T)
      RESULTS$width50_3[i] = mean((output$value3[output$quant == 0.75]-output$value3[output$quant == 0.25])/output$truth, na.rm=T)
      RESULTS$width50_4[i] = mean((output$value4[output$quant == 0.75]-output$value4[output$quant == 0.25])/output$truth, na.rm=T)

      RESULTS$width80_1[i] = mean((output$value[output$quant == 0.9]-output$value[output$quant == 0.1])/output$truth, na.rm=T)
      RESULTS$width80_2[i] = mean((output$value2[output$quant == 0.9]-output$value2[output$quant == 0.1])/output$truth, na.rm=T)
      RESULTS$width80_3[i] = mean((output$value3[output$quant == 0.9]-output$value3[output$quant == 0.1])/output$truth, na.rm=T)
      RESULTS$width80_4[i] = mean((output$value4[output$quant == 0.9]-output$value4[output$quant == 0.1])/output$truth, na.rm=T)

      RESULTS$width95_1[i] = mean((output$value[output$quant == 0.975]-output$value[output$quant == 0.025])/output$truth, na.rm=T)
      RESULTS$width95_2[i] = mean((output$value2[output$quant == 0.975]-output$value2[output$quant == 0.025])/output$truth, na.rm=T)
      RESULTS$width95_3[i] = mean((output$value3[output$quant == 0.975]-output$value3[output$quant == 0.025])/output$truth, na.rm=T)
      RESULTS$width95_4[i] = mean((output$value4[output$quant == 0.975]-output$value4[output$quant == 0.025])/output$truth, na.rm=T)

      
    }else{
      print(paste0('Skipped: ', FILES[i]))
    }
  }
}

RESULTS_ALL = RESULTS[RESULTS$N > 0,]

RESULTS_ALL$disease_source = gsub('real_uq_mat_','',RESULTS_ALL$FILES)
RESULTS_ALL$disease_source = gsub('.csv','',RESULTS_ALL$disease_source)
RESULTS_ALL$disease_source = gsub('(.*)_\\w+', '\\1', RESULTS_ALL$disease_source)
RESULTS_ALL$disease_source = gsub('[0-9]', '', RESULTS_ALL$disease_source)
RESULTS_ALL$disease_source = substr(RESULTS_ALL$disease_source, 1,nchar(RESULTS_ALL$disease_source)-2)
RESULTS_ALL$plot_id = RESULTS_ALL$disease_source

write.csv(RESULTS_ALL, file = paste0('~/GitLab/frodo/summaries_uq/tables/',"lstm.csv"), quote = F, row.names = F)

RESULTS_ALL = RESULTS_ALL %>% dplyr::group_by(plot_id) %>% dplyr::mutate(N_agg = sum(N))
RESULTS_AGG = RESULTS_ALL %>% dplyr::group_by(plot_id) %>% dplyr::mutate(wis_1_agg = sum(N*wis_1)/N_agg,
                                                                         wis_2_agg = sum(N*wis_2)/N_agg,
                                                                         wis_3_agg = sum(N*wis_3)/N_agg,
                                                                         wis_4_agg = sum(N*wis_4)/N_agg,
                                                                         wis_dispersion_1_agg = sum(N*wis_dispersion_1)/N_agg,
                                                                         wis_dispersion_2_agg = sum(N*wis_dispersion_2)/N_agg,
                                                                         wis_dispersion_3_agg = sum(N*wis_dispersion_3)/N_agg,
                                                                         wis_dispersion_4_agg = sum(N*wis_dispersion_4)/N_agg,
                                                                         wis_overprediction_1_agg = sum(N*wis_overprediction_1)/N_agg,
                                                                         wis_overprediction_2_agg = sum(N*wis_overprediction_2)/N_agg,
                                                                         wis_overprediction_3_agg = sum(N*wis_overprediction_3)/N_agg,
                                                                         wis_overprediction_4_agg = sum(N*wis_overprediction_4)/N_agg,
                                                                         wis_underprediction_1_agg = sum(N*wis_underprediction_1)/N_agg,
                                                                         wis_underprediction_2_agg = sum(N*wis_underprediction_2)/N_agg,
                                                                         wis_underprediction_3_agg = sum(N*wis_underprediction_3)/N_agg,
                                                                         wis_underprediction_4_agg = sum(N*wis_underprediction_4)/N_agg,
                                                                         coverage50_1_agg = sum(N*coverage50_1)/N_agg,
                                                                         coverage50_2_agg = sum(N*coverage50_2)/N_agg,
                                                                         coverage50_3_agg = sum(N*coverage50_3)/N_agg,
                                                                         coverage50_4_agg = sum(N*coverage50_4)/N_agg,
                                                                         coverage80_1_agg = sum(N*coverage80_1)/N_agg,
                                                                         coverage80_2_agg = sum(N*coverage80_2)/N_agg,
                                                                         coverage80_3_agg = sum(N*coverage80_3)/N_agg,
                                                                         coverage80_4_agg = sum(N*coverage80_4)/N_agg,
                                                                         coverage95_1_agg = sum(N*coverage95_1)/N_agg,
                                                                         coverage95_2_agg = sum(N*coverage95_2)/N_agg,
                                                                         coverage95_3_agg = sum(N*coverage95_3)/N_agg,
                                                                         coverage95_4_agg = sum(N*coverage95_4)/N_agg,
                                                                         width50_1_agg = sum(N*width50_1)/N_agg,
                                                                         width50_2_agg = sum(N*width50_2)/N_agg,
                                                                         width50_3_agg = sum(N*width50_3)/N_agg,
                                                                         width50_4_agg = sum(N*width50_4)/N_agg,
                                                                         width80_1_agg = sum(N*width80_1)/N_agg,
                                                                         width80_2_agg = sum(N*width80_2)/N_agg,
                                                                         width80_3_agg = sum(N*width80_3)/N_agg,
                                                                         width80_4_agg = sum(N*width80_4)/N_agg,
                                                                         width95_1_agg = sum(N*width95_1)/N_agg,
                                                                         width95_2_agg = sum(N*width95_2)/N_agg,
                                                                         width95_3_agg = sum(N*width95_3)/N_agg,
                                                                         width95_4_agg = sum(N*width95_4)/N_agg)

RESULTS_AGG = RESULTS_AGG[!duplicated(RESULTS_AGG$plot_id),]

write.csv(RESULTS_AGG, file = paste0('~/GitLab/frodo/summaries_uq/tables/',"lstm_agg.csv"), quote = F, row.names = F)


############
### SMOA ###
############



### results
my_path = '~/GitLab/frodo/summaries_uq/results/'
eval_path = '~/GitLab/frodo/uq/smoanodiff/all/'

### modespecific ###
eval_path2 = '~/GitLab/frodo/uq/smoanodiff/modespecific/'

### diseasespecific ###
eval_path3 = '~/GitLab/frodo/uq/smoanodiff/diseasespecific/'

### sourcespecific ###
eval_path4 = '~/GitLab/frodo/uq/smoanodiff/sourcespecific/'



FILES = list.files(paste0(eval_path4))

RESULTS = data.frame(FILES = FILES)
RESULTS$geography = NA
RESULTS$disease = NA
RESULTS$disease_source = NA
RESULTS$N = NA
RESULTS$wis_1 = NA
RESULTS$wis_2 = NA
RESULTS$wis_3 = NA
RESULTS$wis_4 = NA

RESULTS$wis_dispersion_1 = NA
RESULTS$wis_dispersion_2 = NA
RESULTS$wis_dispersion_3 = NA
RESULTS$wis_dispersion_4 = NA

RESULTS$wis_overprediction_1 = NA
RESULTS$wis_overprediction_2 = NA
RESULTS$wis_overprediction_3 = NA
RESULTS$wis_overprediction_4 = NA

RESULTS$wis_underprediction_1 = NA
RESULTS$wis_underprediction_2 = NA
RESULTS$wis_underprediction_3 = NA
RESULTS$wis_underprediction_4 = NA

RESULTS$coverage50_1 = NA
RESULTS$coverage50_2 = NA
RESULTS$coverage50_3 = NA
RESULTS$coverage50_4 = NA
RESULTS$coverage80_1 = NA
RESULTS$coverage80_2 = NA
RESULTS$coverage80_3 = NA
RESULTS$coverage80_4 = NA
RESULTS$coverage95_1 = NA
RESULTS$coverage95_2 = NA
RESULTS$coverage95_3 = NA
RESULTS$coverage95_4 = NA

RESULTS$width50_1 = NA
RESULTS$width50_2 = NA
RESULTS$width50_3 = NA
RESULTS$width50_4 = NA
RESULTS$width80_1 = NA
RESULTS$width80_2 = NA
RESULTS$width80_3 = NA
RESULTS$width80_4 = NA
RESULTS$width95_1 = NA
RESULTS$width95_2 = NA
RESULTS$width95_3 = NA
RESULTS$width95_4 = NA

for(i in 1:length(FILES)){
  if(file.exists(paste0(eval_path4, FILES[i]))){
    output = try(data.frame(data.table::fread(paste0(eval_path, FILES[i]))), silent = T)
    output2 = try(data.frame(data.table::fread(paste0(eval_path2, FILES[i]))), silent = T)
    output3 = try(data.frame(data.table::fread(paste0(eval_path3, FILES[i]))), silent = T)
    output4 = try(data.frame(data.table::fread(paste0(eval_path4, FILES[i]))), silent = T)
    if(grepl('COVID',FILES[i])){output3 = output4} #didn't bother running disease-specific for COVID
    if(grepl('Dengue',FILES[i])){output3 = output4} #didn't bother running disease-specific for Dengue
    if(class(output)[1] != 'try-error' & class(output2)[1] != 'try-error' & class(output3)[1] != 'try-error' & class(output4)[1] != 'try-error'){
      
      
      output$fcst = pmax(0,output$fcst)
      output$fcst[!is.na(output$obs) & output$obs == 0] = 0
      output = output[output$obs > 0,]
      output = output[order(output$quant, output$last_obs_time, output$h),]
      output = output %>% dplyr::group_by(last_obs_time, h) %>% dplyr::mutate(value = sort(value))
      
      
      output2$fcst = pmax(0,output2$fcst)
      output2$fcst[!is.na(output2$obs) & output2$obs == 0] = 0
      output2 = output2[output2$obs > 0,]
      output2 = output2[order(output2$quant, output2$last_obs_time, output2$h),]
      output2 = output2 %>% dplyr::group_by(last_obs_time, h) %>% dplyr::mutate(value = sort(value))
      
      output3$fcst = pmax(0,output3$fcst)
      output3$fcst[!is.na(output3$obs) & output3$obs == 0] = 0
      output3 = output3[output3$obs > 0,] 
      output3 = output3[order(output3$quant, output3$last_obs_time, output3$h),]
      output3 = output3 %>% dplyr::group_by(last_obs_time, h) %>% dplyr::mutate(value = sort(value))
      
      output4$fcst = pmax(0,output4$fcst)
      output4$fcst[!is.na(output4$obs) & output4$obs == 0] = 0
      output4 = output4[output4$obs > 0,] 
      output4 = output4[order(output4$quant, output4$last_obs_time, output4$h),]
      output4 = output4 %>% dplyr::group_by(last_obs_time, h) %>% dplyr::mutate(value = sort(value))
      
      output2$fcst2 = output2$fcst
      output2$value2 = output2$value
      output3$fcst3 = output3$fcst
      output3$value3 = output3$value
      output4$fcst4 = output4$fcst
      output4$value4 = output4$value
      
      output$value[output$value > 3*output$obs & output$obs != 0] = 3*output$obs[output$value > 3*output$obs & output$obs != 0]
      output2$value2[output2$value2 > 3*output2$obs & output2$obs != 0] = 3*output2$obs[output2$value2 > 3*output2$obs & output2$obs != 0]
      output3$value3[output3$value3 > 3*output3$obs & output3$obs != 0] = 3*output3$obs[output3$value3 > 3*output3$obs & output3$obs != 0]
      output4$value4[output4$value4 > 3*output4$obs & output4$obs != 0] = 3*output4$obs[output4$value4 > 3*output4$obs & output4$obs != 0]
 
      #merge by date, last_obs_time got messed up somehow
      output = merge(output, output2[,c('date', 'h', 'fcst2','value2','quant')], by = c('date','h','quant'), all.x = T, all.y = T)
      output = merge(output, output3[,c('date', 'h', 'fcst3','value3','quant')], by = c('date','h','quant'), all.x = T, all.y = T)
      output = merge(output, output4[,c('date', 'h', 'fcst4','value4','quant')], by = c('date','h','quant'), all.x = T, all.y = T)
      output = output[!is.na(output$fcst) & !is.na(output$fcst2) & !is.na(output$fcst3) & !is.na(output$fcst4),]
      output = output[!is.na(output$value) & !is.na(output$value2) & !is.na(output$value3) & !is.na(output$value4),]
      
      output = output[order(output$date, output$h, output$quant),]
      
      
      WIS_1 = (abs(output$fcst[output$quant == 0.5] - output$truth[output$quant == 0.5])/2)+
        (0.5/2)*interval_score(uppers = output$value[output$quant == 0.75], lowers = output$value[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*interval_score(uppers = output$value[output$quant == 0.9], lowers = output$value[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*interval_score(uppers = output$value[output$quant == 0.975], lowers = output$value[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_2 = (abs(output$fcst2[output$quant == 0.5] - output$truth[output$quant == 0.5])/2)+
        (0.5/2)*interval_score(uppers = output$value2[output$quant == 0.75], lowers = output$value2[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*interval_score(uppers = output$value2[output$quant == 0.9], lowers = output$value2[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*interval_score(uppers = output$value2[output$quant == 0.975], lowers = output$value2[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_3 = (abs(output$fcst3[output$quant == 0.5] - output$truth[output$quant == 0.5])/2)+
        (0.5/2)*interval_score(uppers = output$value3[output$quant == 0.75], lowers = output$value3[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*interval_score(uppers = output$value3[output$quant == 0.9], lowers = output$value3[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*interval_score(uppers = output$value3[output$quant == 0.975], lowers = output$value3[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_4 = (abs(output$fcst4[output$quant == 0.5] - output$truth[output$quant == 0.5])/2)+
        (0.5/2)*interval_score(uppers = output$value4[output$quant == 0.75], lowers = output$value4[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*interval_score(uppers = output$value4[output$quant == 0.9], lowers = output$value4[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*interval_score(uppers = output$value4[output$quant == 0.975], lowers = output$value4[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_DISPERSION_1 = (0.5/2)*dispersion_score(uppers = output$value[output$quant == 0.75], lowers = output$value[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*dispersion_score(uppers = output$value[output$quant == 0.9], lowers = output$value[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*dispersion_score(uppers = output$value[output$quant == 0.975], lowers = output$value[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_DISPERSION_2 = (0.5/2)*dispersion_score(uppers = output$value2[output$quant == 0.75], lowers = output$value2[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*dispersion_score(uppers = output$value2[output$quant == 0.9], lowers = output$value2[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*dispersion_score(uppers = output$value2[output$quant == 0.975], lowers = output$value2[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_DISPERSION_3 = (0.5/2)*dispersion_score(uppers = output$value3[output$quant == 0.75], lowers = output$value3[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*dispersion_score(uppers = output$value3[output$quant == 0.9], lowers = output$value3[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*dispersion_score(uppers = output$value3[output$quant == 0.975], lowers = output$value3[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_DISPERSION_4 = (0.5/2)*dispersion_score(uppers = output$value4[output$quant == 0.75], lowers = output$value4[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*dispersion_score(uppers = output$value4[output$quant == 0.9], lowers = output$value4[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*dispersion_score(uppers = output$value4[output$quant == 0.975], lowers = output$value4[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_OVERPREDICTION_1 = (0.5/2)*overprediction_score(uppers = output$value[output$quant == 0.75], lowers = output$value[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*overprediction_score(uppers = output$value[output$quant == 0.9], lowers = output$value[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*overprediction_score(uppers = output$value[output$quant == 0.975], lowers = output$value[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_OVERPREDICTION_2 = (0.5/2)*overprediction_score(uppers = output$value2[output$quant == 0.75], lowers = output$value2[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*overprediction_score(uppers = output$value2[output$quant == 0.9], lowers = output$value2[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*overprediction_score(uppers = output$value2[output$quant == 0.975], lowers = output$value2[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_OVERPREDICTION_3 = (0.5/2)*overprediction_score(uppers = output$value3[output$quant == 0.75], lowers = output$value3[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*overprediction_score(uppers = output$value3[output$quant == 0.9], lowers = output$value3[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*overprediction_score(uppers = output$value3[output$quant == 0.975], lowers = output$value3[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_OVERPREDICTION_4 = (0.5/2)*overprediction_score(uppers = output$value4[output$quant == 0.75], lowers = output$value4[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*overprediction_score(uppers = output$value4[output$quant == 0.9], lowers = output$value4[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*overprediction_score(uppers = output$value4[output$quant == 0.975], lowers = output$value4[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])

      WIS_UNDERPREDICTION_1 = (0.5/2)*underprediction_score(uppers = output$value[output$quant == 0.75], lowers = output$value[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*underprediction_score(uppers = output$value[output$quant == 0.9], lowers = output$value[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*underprediction_score(uppers = output$value[output$quant == 0.975], lowers = output$value[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_UNDERPREDICTION_2 = (0.5/2)*underprediction_score(uppers = output$value2[output$quant == 0.75], lowers = output$value2[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*underprediction_score(uppers = output$value2[output$quant == 0.9], lowers = output$value2[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*underprediction_score(uppers = output$value2[output$quant == 0.975], lowers = output$value2[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_UNDERPREDICTION_3 = (0.5/2)*underprediction_score(uppers = output$value3[output$quant == 0.75], lowers = output$value3[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*underprediction_score(uppers = output$value3[output$quant == 0.9], lowers = output$value3[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*underprediction_score(uppers = output$value3[output$quant == 0.975], lowers = output$value3[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      WIS_UNDERPREDICTION_4 = (0.5/2)*underprediction_score(uppers = output$value4[output$quant == 0.75], lowers = output$value4[output$quant == 0.25], alpha_grid = 0.5, y = output$truth[output$quant == 0.5])+
        (0.8/2)*underprediction_score(uppers = output$value4[output$quant == 0.9], lowers = output$value4[output$quant == 0.1], alpha_grid = 0.8, y = output$truth[output$quant == 0.5])+
        (0.95/2)*underprediction_score(uppers = output$value4[output$quant == 0.975], lowers = output$value4[output$quant == 0.025], alpha_grid = 0.95, y = output$truth[output$quant == 0.5])
      
      RESULTS$N[i] = nrow(output)
      RESULTS$wis_1[i] = (1/(3+0.5))*mean(WIS_1)
      RESULTS$wis_2[i] =  (1/(3+0.5))*mean(WIS_2)
      RESULTS$wis_3[i] =  (1/(3+0.5))*mean(WIS_3)
      RESULTS$wis_4[i] =  (1/(3+0.5))*mean(WIS_4)
      
      RESULTS$wis_dispersion_1[i] = (1/(3+0.5))*mean(WIS_DISPERSION_1)
      RESULTS$wis_dispersion_2[i] = (1/(3+0.5))*mean(WIS_DISPERSION_2)
      RESULTS$wis_dispersion_3[i] = (1/(3+0.5))*mean(WIS_DISPERSION_3)
      RESULTS$wis_dispersion_4[i] = (1/(3+0.5))*mean(WIS_DISPERSION_4)

      RESULTS$wis_overprediction_1[i] = (1/(3+0.5))*mean(WIS_OVERPREDICTION_1)
      RESULTS$wis_overprediction_2[i] = (1/(3+0.5))*mean(WIS_OVERPREDICTION_2)
      RESULTS$wis_overprediction_3[i] = (1/(3+0.5))*mean(WIS_OVERPREDICTION_3)
      RESULTS$wis_overprediction_4[i] = (1/(3+0.5))*mean(WIS_OVERPREDICTION_4)
      
      RESULTS$wis_underprediction_1[i] = (1/(3+0.5))*mean(WIS_UNDERPREDICTION_1)
      RESULTS$wis_underprediction_2[i] = (1/(3+0.5))*mean(WIS_UNDERPREDICTION_2)
      RESULTS$wis_underprediction_3[i] = (1/(3+0.5))*mean(WIS_UNDERPREDICTION_3)
      RESULTS$wis_underprediction_4[i] = (1/(3+0.5))*mean(WIS_UNDERPREDICTION_4)
      
      RESULTS$coverage50_1[i] = mean(output$truth[output$quant == 0.5] <= output$value[output$quant == 0.75] & output$truth[output$quant == 0.5] >= output$value[output$quant == 0.25], na.rm=T)
      RESULTS$coverage50_2[i] = mean(output$truth[output$quant == 0.5] <= output$value2[output$quant == 0.75] & output$truth[output$quant == 0.5] >= output$value2[output$quant == 0.25], na.rm=T)
      RESULTS$coverage50_3[i] = mean(output$truth[output$quant == 0.5] <= output$value3[output$quant == 0.75] & output$truth[output$quant == 0.5] >= output$value3[output$quant == 0.25], na.rm=T)
      RESULTS$coverage50_4[i] = mean(output$truth[output$quant == 0.5] <= output$value4[output$quant == 0.75] & output$truth[output$quant == 0.5] >= output$value4[output$quant == 0.25], na.rm=T)

      RESULTS$coverage80_1[i] = mean(output$truth[output$quant == 0.5] <= output$value[output$quant == 0.9] & output$truth[output$quant == 0.5] >= output$value[output$quant == 0.1], na.rm=T)
      RESULTS$coverage80_2[i] = mean(output$truth[output$quant == 0.5] <= output$value2[output$quant == 0.9] & output$truth[output$quant == 0.5] >= output$value2[output$quant == 0.1], na.rm=T)
      RESULTS$coverage80_3[i] = mean(output$truth[output$quant == 0.5] <= output$value3[output$quant == 0.9] & output$truth[output$quant == 0.5] >= output$value3[output$quant == 0.1], na.rm=T)
      RESULTS$coverage80_4[i] = mean(output$truth[output$quant == 0.5] <= output$value4[output$quant == 0.9] & output$truth[output$quant == 0.5] >= output$value4[output$quant == 0.1], na.rm=T)

      RESULTS$coverage95_1[i] = mean(output$truth[output$quant == 0.5] <= output$value[output$quant == 0.975] & output$truth[output$quant == 0.5] >= output$value[output$quant == 0.025], na.rm=T)
      RESULTS$coverage95_2[i] = mean(output$truth[output$quant == 0.5] <= output$value2[output$quant == 0.975] & output$truth[output$quant == 0.5] >= output$value2[output$quant == 0.025], na.rm=T)
      RESULTS$coverage95_3[i] = mean(output$truth[output$quant == 0.5] <= output$value3[output$quant == 0.975] & output$truth[output$quant == 0.5] >= output$value3[output$quant == 0.025], na.rm=T)
      RESULTS$coverage95_4[i] = mean(output$truth[output$quant == 0.5] <= output$value4[output$quant == 0.975] & output$truth[output$quant == 0.5] >= output$value4[output$quant == 0.025], na.rm=T)

      RESULTS$width50_1[i] = mean((output$value[output$quant == 0.75]-output$value[output$quant == 0.25])/output$truth, na.rm=T)
      RESULTS$width50_2[i] = mean((output$value2[output$quant == 0.75]-output$value2[output$quant == 0.25])/output$truth, na.rm=T)
      RESULTS$width50_3[i] = mean((output$value3[output$quant == 0.75]-output$value3[output$quant == 0.25])/output$truth, na.rm=T)
      RESULTS$width50_4[i] = mean((output$value4[output$quant == 0.75]-output$value4[output$quant == 0.25])/output$truth, na.rm=T)

      RESULTS$width80_1[i] = mean((output$value[output$quant == 0.9]-output$value[output$quant == 0.1])/output$truth, na.rm=T)
      RESULTS$width80_2[i] = mean((output$value2[output$quant == 0.9]-output$value2[output$quant == 0.1])/output$truth, na.rm=T)
      RESULTS$width80_3[i] = mean((output$value3[output$quant == 0.9]-output$value3[output$quant == 0.1])/output$truth, na.rm=T)
      RESULTS$width80_4[i] = mean((output$value4[output$quant == 0.9]-output$value4[output$quant == 0.1])/output$truth, na.rm=T)

      RESULTS$width95_1[i] = mean((output$value[output$quant == 0.975]-output$value[output$quant == 0.025])/output$truth, na.rm=T)
      RESULTS$width95_2[i] = mean((output$value2[output$quant == 0.975]-output$value2[output$quant == 0.025])/output$truth, na.rm=T)
      RESULTS$width95_3[i] = mean((output$value3[output$quant == 0.975]-output$value3[output$quant == 0.025])/output$truth, na.rm=T)
      RESULTS$width95_4[i] = mean((output$value4[output$quant == 0.975]-output$value4[output$quant == 0.025])/output$truth, na.rm=T)

      
    }else{
      print(paste0('Skipped: ', FILES[i]))
    }
  }
}


RESULTS_ALL = RESULTS[RESULTS$N > 0,]

RESULTS_ALL$disease_source = gsub('real_uq_mat_','',RESULTS_ALL$FILES)
RESULTS_ALL$disease_source = gsub('.csv','',RESULTS_ALL$disease_source)
RESULTS_ALL$disease_source = gsub('(.*)_\\w+', '\\1', RESULTS_ALL$disease_source)
RESULTS_ALL$disease_source = gsub('[0-9]', '', RESULTS_ALL$disease_source)
RESULTS_ALL$plot_id = RESULTS_ALL$disease_source

write.csv(RESULTS_ALL, file = paste0('~/GitLab/frodo/summaries_uq/tables/',"smoanodiff.csv"), quote = F, row.names = F)


RESULTS_ALL = RESULTS_ALL %>% dplyr::group_by(plot_id) %>% dplyr::mutate(N_agg = sum(N))
RESULTS_AGG = RESULTS_ALL %>% dplyr::group_by(plot_id) %>% dplyr::mutate(wis_1_agg = sum(N*wis_1)/N_agg,
                                                                         wis_2_agg = sum(N*wis_2)/N_agg,
                                                                         wis_3_agg = sum(N*wis_3)/N_agg,
                                                                         wis_4_agg = sum(N*wis_4)/N_agg,
                                                                         wis_dispersion_1_agg = sum(N*wis_dispersion_1)/N_agg,
                                                                         wis_dispersion_2_agg = sum(N*wis_dispersion_2)/N_agg,
                                                                         wis_dispersion_3_agg = sum(N*wis_dispersion_3)/N_agg,
                                                                         wis_dispersion_4_agg = sum(N*wis_dispersion_4)/N_agg,
                                                                         wis_overprediction_1_agg = sum(N*wis_overprediction_1)/N_agg,
                                                                         wis_overprediction_2_agg = sum(N*wis_overprediction_2)/N_agg,
                                                                         wis_overprediction_3_agg = sum(N*wis_overprediction_3)/N_agg,
                                                                         wis_overprediction_4_agg = sum(N*wis_overprediction_4)/N_agg,
                                                                         wis_underprediction_1_agg = sum(N*wis_underprediction_1)/N_agg,
                                                                         wis_underprediction_2_agg = sum(N*wis_underprediction_2)/N_agg,
                                                                         wis_underprediction_3_agg = sum(N*wis_underprediction_3)/N_agg,
                                                                         wis_underprediction_4_agg = sum(N*wis_underprediction_4)/N_agg,
                                                                         coverage50_1_agg = sum(N*coverage50_1)/N_agg,
                                                                         coverage50_2_agg = sum(N*coverage50_2)/N_agg,
                                                                         coverage50_3_agg = sum(N*coverage50_3)/N_agg,
                                                                         coverage50_4_agg = sum(N*coverage50_4)/N_agg,
                                                                         coverage80_1_agg = sum(N*coverage80_1)/N_agg,
                                                                         coverage80_2_agg = sum(N*coverage80_2)/N_agg,
                                                                         coverage80_3_agg = sum(N*coverage80_3)/N_agg,
                                                                         coverage80_4_agg = sum(N*coverage80_4)/N_agg,
                                                                         coverage95_1_agg = sum(N*coverage95_1)/N_agg,
                                                                         coverage95_2_agg = sum(N*coverage95_2)/N_agg,
                                                                         coverage95_3_agg = sum(N*coverage95_3)/N_agg,
                                                                         coverage95_4_agg = sum(N*coverage95_4)/N_agg,
                                                                         width50_1_agg = sum(N*width50_1)/N_agg,
                                                                         width50_2_agg = sum(N*width50_2)/N_agg,
                                                                         width50_3_agg = sum(N*width50_3)/N_agg,
                                                                         width50_4_agg = sum(N*width50_4)/N_agg,
                                                                         width80_1_agg = sum(N*width80_1)/N_agg,
                                                                         width80_2_agg = sum(N*width80_2)/N_agg,
                                                                         width80_3_agg = sum(N*width80_3)/N_agg,
                                                                         width80_4_agg = sum(N*width80_4)/N_agg,
                                                                         width95_1_agg = sum(N*width95_1)/N_agg,
                                                                         width95_2_agg = sum(N*width95_2)/N_agg,
                                                                         width95_3_agg = sum(N*width95_3)/N_agg,
                                                                         width95_4_agg = sum(N*width95_4)/N_agg)

RESULTS_AGG = RESULTS_AGG[!duplicated(RESULTS_AGG$plot_id),]


write.csv(RESULTS_AGG, file = paste0('~/GitLab/frodo/summaries_uq/tables/',"smoanodiff_agg.csv"), quote = F, row.names = F)


