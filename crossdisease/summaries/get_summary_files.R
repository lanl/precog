#########################
#########################
### Get All Summaries ###
#########################
#########################




library(data.table)
library(GGally)
lut_path = '~/GitLab/frodo/lookup_tables/'




#############
### FRODO ###
#############

### results
my_path = '~/GitLab/frodo/summaries/results/'
eval_path = '~/GitLab/frodo/evaluations/frodo/all/'

### results_modespecific ###
eval_path2 = '~/GitLab/frodo/evaluations/frodo/modespecific_cold/'

### results_diseasespecific ###
eval_path3 = '~/GitLab/frodo/evaluations/frodo/diseasespecific_cold/'

### results_sourcespecific ###
eval_path4 = '~/GitLab/frodo/evaluations/frodo/sourcespecific_cold/'


FILES = list.files(paste0(eval_path4))
FILES = FILES[grepl('_0.5',FILES)]


RESULTS = data.frame(FILES = FILES)
RESULTS$geography = NA
RESULTS$disease = NA
RESULTS$disease_source = NA
RESULTS$N = NA
RESULTS$mae = NA
RESULTS$rmse = NA
RESULTS$mae2 = NA
RESULTS$rmse2 = NA
RESULTS$mae3 = NA
RESULTS$rmse3 = NA
RESULTS$mae4 = NA
RESULTS$rmse4 = NA
RESULTS$mae_rw = NA
RESULTS$rmse_rw = NA
RESULTS$truth_avg = NA
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
      output$fcst[!is.na(output$obs_smoothed) & output$obs_smoothed == 0] = 0
      # output = output[output$truth > 0,]
      
      output2$fcst = pmax(0,output2$fcst)
      output2$fcst[!is.na(output2$obs_smoothed) & output2$obs_smoothed == 0] = 0
      # output2 = output2[output2$truth > 0,]
      output2$fcst2 = output2$fcst
      
      output3$fcst = pmax(0,output3$fcst)
      output3$fcst[!is.na(output3$obs_smoothed) & output3$obs_smoothed == 0] = 0
      # output3 = output3[output3$truth > 0,] 
      output3$fcst3 = output3$fcst
      
      output4$fcst = pmax(0,output4$fcst)
      output4$fcst[!is.na(output4$obs_smoothed) & output4$obs_smoothed == 0] = 0
      # output4 = output4[output4$truth > 0,] 
      output4$fcst4 = output4$fcst
      
      output = merge(output, output2[,c('last_obs_time', 'h', 'fcst2')], by = c('last_obs_time','h'), all.x = T, all.y = T)
      output = merge(output, output3[,c('last_obs_time', 'h', 'fcst3')], by = c('last_obs_time','h'), all.x = T, all.y = T)
      output = merge(output, output4[,c('last_obs_time', 'h', 'fcst4')], by = c('last_obs_time','h'), all.x = T, all.y = T)
      output = output[!is.na(output$fcst) & !is.na(output$fcst2) & !is.na(output$fcst3) & !is.na(output$fcst4),]
      
      RESULTS$geography[i] = output$geography[1]
      RESULTS$disease[i] = output$disease[1]
      RESULTS$disease_source[i] = output$disease_source[1]
      RESULTS$N[i] = nrow(output)
      RESULTS$mae[i] = mean(abs(output$truth - output$fcst), na.rm=T)
      RESULTS$rmse[i] = sqrt(mean((output$truth - output$fcst)^2, na.rm=T))
      RESULTS$mae2[i] = mean(abs(output$truth - output$fcst2), na.rm=T)
      RESULTS$rmse2[i] = sqrt(mean((output$truth - output$fcst2)^2, na.rm=T))
      RESULTS$mae3[i] = mean(abs(output$truth - output$fcst3), na.rm=T)
      RESULTS$rmse3[i] = sqrt(mean((output$truth - output$fcst3)^2, na.rm=T))
      RESULTS$mae4[i] = mean(abs(output$truth - output$fcst4), na.rm=T)
      RESULTS$rmse4[i] = sqrt(mean((output$truth - output$fcst4)^2, na.rm=T))
      RESULTS$mae_rw[i] = mean(abs(output$truth - output$obs), na.rm=T)
      RESULTS$rmse_rw[i] = sqrt(mean((output$truth - output$obs)^2, na.rm=T))
      RESULTS$truth_avg[i] = mean(output$truth, na.rm=T)
      
    }else{
      print(paste0('Skipped: ', FILES[i]))
    }
  }
}



RESULTS = RESULTS[!is.na(RESULTS$disease),]

FILENAME = paste0('~/GitLab/frodo/models/results/cat_to_numeric.RDS')
LEVELS_STORAGE <- readRDS(FILENAME)
RESULTS$disease_store = RESULTS$disease
RESULTS$disease = LEVELS_STORAGE$disease[RESULTS$disease]


RESULTS_ALL = RESULTS
RESULTS_ALL$plot_id = RESULTS_ALL$disease_source

write.csv(RESULTS_ALL, file = paste0('~/GitLab/frodo/summaries/tables/',"frodo.csv"), quote = F, row.names = F)

RESULTS_ALL = RESULTS_ALL[!is.na(RESULTS_ALL$mae_rw),] 
RESULTS_ALL = RESULTS_ALL[RESULTS_ALL$N > 0,]


RESULTS_AGG = RESULTS_ALL %>% dplyr::group_by(plot_id) %>% dplyr::mutate(mae_agg = sum(N*mae),
                                                                          rmse_agg = sum(N*(rmse^2)),
                                                                          mae2_agg = sum(N*mae2),
                                                                          rmse2_agg = sum(N*(rmse2^2)),
                                                                          mae3_agg = sum(N*mae3),
                                                                          rmse3_agg = sum(N*(rmse3^2)),
                                                                          mae4_agg = sum(N*mae4),
                                                                          rmse4_agg = sum(N*(rmse4^2)),
                                                                          mae_rw_agg = sum(N*mae_rw),
                                                                          rmse_rw_agg = sum(N*(rmse_rw^2)),
                                                                          N_agg = sum(N))
RESULTS_AGG$mae_agg = RESULTS_AGG$mae_agg/RESULTS_AGG$N_agg
RESULTS_AGG$rmse_agg = sqrt(RESULTS_AGG$rmse_agg/RESULTS_AGG$N_agg)
RESULTS_AGG$mae2_agg = RESULTS_AGG$mae2_agg/RESULTS_AGG$N_agg
RESULTS_AGG$rmse2_agg = sqrt(RESULTS_AGG$rmse2_agg/RESULTS_AGG$N_agg)
RESULTS_AGG$mae3_agg = RESULTS_AGG$mae3_agg/RESULTS_AGG$N_agg
RESULTS_AGG$rmse3_agg = sqrt(RESULTS_AGG$rmse3_agg/RESULTS_AGG$N_agg)
RESULTS_AGG$mae4_agg = RESULTS_AGG$mae4_agg/RESULTS_AGG$N_agg
RESULTS_AGG$rmse4_agg = sqrt(RESULTS_AGG$rmse4_agg/RESULTS_AGG$N_agg)
RESULTS_AGG$mae_rw_agg = RESULTS_AGG$mae_rw_agg/RESULTS_AGG$N_agg
RESULTS_AGG$rmse_rw_agg = sqrt(RESULTS_AGG$rmse_rw_agg/RESULTS_AGG$N_agg)

RESULTS_AGG = RESULTS_AGG[!duplicated(RESULTS_AGG$plot_id),]

write.csv(RESULTS_AGG, file = paste0('~/GitLab/frodo/summaries/tables/',"frodo_agg.csv"), quote = F, row.names = F)









############
### LSTM ###
############


### results
my_path = '~/GitLab/frodo/summaries/lstm/all/'
eval_path = '~/GitLab/frodo/evaluations/lstm/all/'

### results_modespecific ###
eval_path2 = '~/GitLab/frodo/evaluations/lstm/modespecific_cold/'

### results_diseasespecific ###
eval_path3 = '~/GitLab/frodo/evaluations/lstm/diseasespecific_cold/'

### results_sourcespecific ###
eval_path4 = '~/GitLab/frodo/evaluations/lstm/sourcespecific_cold/'

FILES = list.files(paste0(eval_path4))
FILES = FILES[grepl('_0.5',FILES)]

RESULTS = data.frame(FILES = FILES)
RESULTS$N = NA
RESULTS$mae = NA
RESULTS$rmse = NA
RESULTS$mae2 = NA
RESULTS$rmse2 = NA
RESULTS$mae3 = NA
RESULTS$rmse3 = NA
RESULTS$mae4 = NA
RESULTS$rmse4 = NA
RESULTS$mae_rw = NA
RESULTS$rmse_rw = NA
RESULTS$truth_avg = NA
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
      # output = output[output$truth > 0,]
      
      
      output2$fcst = pmax(0,output2$fcst)
      output2$fcst[!is.na(output2$obs) & output2$obs == 0] = 0
      # output2 = output2[output2$truth > 0,]
      output2$fcst2 = output2$fcst
      
      output3$fcst = pmax(0,output3$fcst)
      output3$fcst[!is.na(output3$obs) & output3$obs == 0] = 0
      # output3 = output3[output3$truth > 0,] 
      output3$fcst3 = output3$fcst
      
      output4$fcst = pmax(0,output4$fcst)
      output4$fcst[!is.na(output4$obs) & output4$obs == 0] = 0
      # output4 = output4[output4$truth > 0,] 
      output4$fcst4 = output4$fcst
      
      output = merge(output, output2[,c('date', 'h', 'fcst2')], by = c('date','h'), all.x = T, all.y = T)
      output = merge(output, output3[,c('date', 'h', 'fcst3')], by = c('date','h'), all.x = T, all.y = T)
      output = merge(output, output4[,c('date', 'h', 'fcst4')], by = c('date','h'), all.x = T, all.y = T)
      output = output[!is.na(output$fcst) & !is.na(output$fcst2) & !is.na(output$fcst3) & !is.na(output$fcst4),]

      RESULTS$N[i] = nrow(output)
      RESULTS$mae[i] = mean(abs(output$truth - output$fcst), na.rm=T)
      RESULTS$rmse[i] = sqrt(mean((output$truth - output$fcst)^2, na.rm=T))
      RESULTS$mae2[i] = mean(abs(output$truth - output$fcst2), na.rm=T)
      RESULTS$rmse2[i] = sqrt(mean((output$truth - output$fcst2)^2, na.rm=T))
      RESULTS$mae3[i] = mean(abs(output$truth - output$fcst3), na.rm=T)
      RESULTS$rmse3[i] = sqrt(mean((output$truth - output$fcst3)^2, na.rm=T))
      RESULTS$mae4[i] = mean(abs(output$truth - output$fcst4), na.rm=T)
      RESULTS$rmse4[i] = sqrt(mean((output$truth - output$fcst4)^2, na.rm=T))
      RESULTS$mae_rw[i] = mean(abs(output$truth - output$obs), na.rm=T)
      RESULTS$rmse_rw[i] = sqrt(mean((output$truth - output$obs)^2, na.rm=T))
      RESULTS$truth_avg[i] = mean(output$truth, na.rm=T)
      
    }else{
      print(paste0('Skipped: ', FILES[i]))
    }
  }
}

RESULTS$temp = gsub('.csv','',gsub('real_eval_mat_','',RESULTS$FILES))
RESULTS$temp = gsub('_0.5','',RESULTS$temp)



SPLITS = strsplit(RESULTS$temp, split = '_')
RESULTS$disease_source = NA
RESULTS$disease = NA
for(i in 1:length(SPLITS)){
  if(length(SPLITS[[i]])==3){
    RESULTS$disease_source[i] = paste0(SPLITS[[i]][1],'_',SPLITS[[i]][2])
    RESULTS$disease[i] = SPLITS[[i]][1]
  }else{
    RESULTS$disease_source[i] = paste0(SPLITS[[i]][1],'_',SPLITS[[i]][2],'_',SPLITS[[i]][3])
    RESULTS$disease[i] = paste0(SPLITS[[i]][1],'_',SPLITS[[i]][2])
  }
}


RESULTS_ALL = RESULTS
RESULTS_ALL$plot_id = RESULTS_ALL$disease_source

RESULTS_ALL = RESULTS_ALL[!is.na(RESULTS_ALL$mae_rw),] 
RESULTS_ALL = RESULTS_ALL[RESULTS_ALL$N > 0,]

write.csv(RESULTS_ALL, file = paste0('~/GitLab/frodo/summaries/tables/',"lstm.csv"), quote = F, row.names = F)

RESULTS_META = RESULTS_ALL


RESULTS_AGG = RESULTS_META %>% dplyr::group_by(plot_id) %>% dplyr::mutate(mae_agg = sum(N*mae),
                                                                          rmse_agg = sum(N*(rmse^2)),
                                                                          mae2_agg = sum(N*mae2),
                                                                          rmse2_agg = sum(N*(rmse2^2)),
                                                                          mae3_agg = sum(N*mae3),
                                                                          rmse3_agg = sum(N*(rmse3^2)),
                                                                          mae4_agg = sum(N*mae4),
                                                                          rmse4_agg = sum(N*(rmse4^2)),
                                                                          mae_rw_agg = sum(N*mae_rw),
                                                                          rmse_rw_agg = sum(N*(rmse_rw^2)),
                                                                          N_agg = sum(N))
RESULTS_AGG$mae_agg = RESULTS_AGG$mae_agg/RESULTS_AGG$N_agg
RESULTS_AGG$rmse_agg = sqrt(RESULTS_AGG$rmse_agg/RESULTS_AGG$N_agg)
RESULTS_AGG$mae2_agg = RESULTS_AGG$mae2_agg/RESULTS_AGG$N_agg
RESULTS_AGG$rmse2_agg = sqrt(RESULTS_AGG$rmse2_agg/RESULTS_AGG$N_agg)
RESULTS_AGG$mae3_agg = RESULTS_AGG$mae3_agg/RESULTS_AGG$N_agg
RESULTS_AGG$rmse3_agg = sqrt(RESULTS_AGG$rmse3_agg/RESULTS_AGG$N_agg)
RESULTS_AGG$mae4_agg = RESULTS_AGG$mae4_agg/RESULTS_AGG$N_agg
RESULTS_AGG$rmse4_agg = sqrt(RESULTS_AGG$rmse4_agg/RESULTS_AGG$N_agg)
RESULTS_AGG$mae_rw_agg = RESULTS_AGG$mae_rw_agg/RESULTS_AGG$N_agg
RESULTS_AGG$rmse_rw_agg = sqrt(RESULTS_AGG$rmse_rw_agg/RESULTS_AGG$N_agg)

RESULTS_AGG = RESULTS_AGG[!duplicated(RESULTS_AGG$plot_id),]

write.csv(RESULTS_AGG, file = paste0('~/GitLab/frodo/summaries/tables/',"lstm_agg.csv"), quote = F, row.names = F)







############
### SMOA ###
############


### all
my_path = '~/GitLab/frodo/summaries/smoanodiff/all/'
eval_path = '~/GitLab/frodo/evaluations/smoanodiff/all/'

### modespecific ###
eval_path2 = '~/GitLab/frodo/evaluations/smoanodiff/modespecific/'

### diseasespecific ###
eval_path3 = '~/GitLab/frodo/evaluations/smoanodiff/diseasespecific/'

### sourcespecific ###
eval_path4 = '~/GitLab/frodo/evaluations/smoanodiff/sourcespecific/'

lut_path = '~/GitLab/frodo/lookup_tables/'

FILES = list.files(paste0(eval_path4))

RESULTS = data.frame(FILES = FILES)
RESULTS$N = NA
RESULTS$mae = NA
RESULTS$rmse = NA
RESULTS$mae2 = NA
RESULTS$rmse2 = NA
RESULTS$mae3 = NA
RESULTS$rmse3 = NA
RESULTS$mae4 = NA
RESULTS$rmse4 = NA
RESULTS$mae_rw = NA
RESULTS$rmse_rw = NA
for(i in 1:length(FILES)){
  if(file.exists(paste0(eval_path4, FILES[i]))){
    
    output = try(data.frame(data.table::fread(paste0(eval_path, FILES[i]))), silent = T)
    output2 = try(data.frame(data.table::fread(paste0(eval_path2, FILES[i]))), silent = T)
    output3 = try(data.frame(data.table::fread(paste0(eval_path3, FILES[i]))), silent = T)
    output4 = try(data.frame(data.table::fread(paste0(eval_path4, FILES[i]))), silent = T)
    if(grepl('COVID',FILES[i])){output3 = output4} #didn't bother running disease-specific for COVID
    if(grepl('Dengue',FILES[i])){output3 = output4} #didn't bother running disease-specific for Dengue
    if(class(output)[1] != 'try-error' & class(output2)[1] != 'try-error' & class(output3)[1] != 'try-error' & class(output4)[1] != 'try-error'){
      
      output$fcst = output$fcst
      output$fcst = pmax(0,output$fcst)
      output$fcst[!is.na(output$obs) & output$obs == 0] = 0
    
      # output = output[output$truth > 0,]
      
      output2$fcst = output2$fcst
      output2$fcst = pmax(0,output2$fcst)
      output2$fcst[!is.na(output2$obs) & output2$obs == 0] = 0
      # output2 = output2[output2$truth > 0,] 
      output2$fcst2 = output2$fcst
      
      output3$fcst = output3$fcst
      output3$fcst = pmax(0,output3$fcst)
      output3$fcst[!is.na(output3$obs) & output3$obs == 0] = 0
      # output3 = output3[output3$truth > 0,] 
      output3$fcst3 = output3$fcst
      
      output4$fcst = output4$fcst
      output4$fcst = pmax(0,output4$fcst)
      output4$fcst[!is.na(output4$obs) & output4$obs == 0] = 0
      # output4 = output4[output4$truth > 0,] 
      output4$fcst4 = output4$fcst
      
      output = merge(output, output2[,c('row_num', 'h', 'fcst2')], by = c('row_num','h'), all.x = T, all.y = T)
      output = merge(output, output3[,c('row_num', 'h', 'fcst3')], by = c('row_num','h'), all.x = T, all.y = T)
      output = merge(output, output4[,c('row_num', 'h', 'fcst4')], by = c('row_num','h'), all.x = T, all.y = T)
      output = output[!is.na(output$fcst) & !is.na(output$fcst2) & !is.na(output$fcst3) & !is.na(output$fcst4),]
      
      RESULTS$N[i] = nrow(output)
      RESULTS$mae[i] = mean(abs(output$truth - output$fcst), na.rm=T)
      RESULTS$rmse[i] = sqrt(mean((output$truth - output$fcst)^2, na.rm=T))
      RESULTS$mae2[i] = mean(abs(output$truth - output$fcst2), na.rm=T)
      RESULTS$rmse2[i] = sqrt(mean((output$truth - output$fcst2)^2, na.rm=T))
      RESULTS$mae3[i] = mean(abs(output$truth - output$fcst3), na.rm=T)
      RESULTS$rmse3[i] = sqrt(mean((output$truth - output$fcst3)^2, na.rm=T))
      RESULTS$mae4[i] = mean(abs(output$truth - output$fcst4), na.rm=T)
      RESULTS$rmse4[i] = sqrt(mean((output$truth - output$fcst4)^2, na.rm=T))
      RESULTS$mae_rw[i] = mean(abs(output$truth - output$obs), na.rm=T)
      RESULTS$rmse_rw[i] = sqrt(mean((output$truth - output$obs)^2, na.rm=T))
      RESULTS$truth_avg[i] = mean(output$truth, na.rm=T)
    }else{
      print(paste0('Skipped: ', FILES[i]))
    }
  }
  #print(paste0('Finished: ', i, ' of ', length(FILES)))
}




RESULTS$temp = gsub('.csv','',gsub('real_eval_mat_','',RESULTS$FILES))
SPLITS = strsplit(RESULTS$temp, split = '_')
RESULTS$disease_source = NA
RESULTS$disease = NA
for(i in 1:length(SPLITS)){
  if(length(SPLITS[[i]])==3){
    RESULTS$disease_source[i] = paste0(SPLITS[[i]][1],'_',SPLITS[[i]][2])
    RESULTS$disease[i] = SPLITS[[i]][1]
  }else{
    RESULTS$disease_source[i] = paste0(SPLITS[[i]][1],'_',SPLITS[[i]][2],'_',SPLITS[[i]][3])
    RESULTS$disease[i] = paste0(SPLITS[[i]][1],'_',SPLITS[[i]][2])
  }
}


RESULTS_ALL = RESULTS
RESULTS_ALL$plot_id = RESULTS_ALL$disease_source

RESULTS_ALL = RESULTS_ALL[!is.na(RESULTS_ALL$mae_rw),]
RESULTS_ALL = RESULTS_ALL[RESULTS_ALL$N > 0,]

write.csv(RESULTS_ALL, file = paste0('~/GitLab/frodo/summaries/tables/',"smoanodiff.csv"), quote = F, row.names = F)

RESULTS_META = RESULTS_ALL

table(RESULTS_META$disease, RESULTS_META$ts_time_cadence)
table(RESULTS_META$plot_id)

RESULTS_AGG = RESULTS_META %>% dplyr::group_by(plot_id) %>% dplyr::mutate(mae_agg = sum(N*mae),
                                                                          rmse_agg = sum(N*(rmse^2)),
                                                                          mae2_agg = sum(N*mae2),
                                                                          rmse2_agg = sum(N*(rmse2^2)),
                                                                          mae3_agg = sum(N*mae3),
                                                                          rmse3_agg = sum(N*(rmse3^2)),
                                                                          mae4_agg = sum(N*mae4),
                                                                          rmse4_agg = sum(N*(rmse4^2)),
                                                                          mae_rw_agg = sum(N*mae_rw),
                                                                          rmse_rw_agg = sum(N*(rmse_rw^2)),
                                                                          N_agg = sum(N))
RESULTS_AGG$mae_agg = RESULTS_AGG$mae_agg/RESULTS_AGG$N_agg
RESULTS_AGG$rmse_agg = sqrt(RESULTS_AGG$rmse_agg/RESULTS_AGG$N_agg)
RESULTS_AGG$mae2_agg = RESULTS_AGG$mae2_agg/RESULTS_AGG$N_agg
RESULTS_AGG$rmse2_agg = sqrt(RESULTS_AGG$rmse2_agg/RESULTS_AGG$N_agg)
RESULTS_AGG$mae3_agg = RESULTS_AGG$mae3_agg/RESULTS_AGG$N_agg
RESULTS_AGG$rmse3_agg = sqrt(RESULTS_AGG$rmse3_agg/RESULTS_AGG$N_agg)
RESULTS_AGG$mae4_agg = RESULTS_AGG$mae4_agg/RESULTS_AGG$N_agg
RESULTS_AGG$rmse4_agg = sqrt(RESULTS_AGG$rmse4_agg/RESULTS_AGG$N_agg)
RESULTS_AGG$mae_rw_agg = RESULTS_AGG$mae_rw_agg/RESULTS_AGG$N_agg
RESULTS_AGG$rmse_rw_agg = sqrt(RESULTS_AGG$rmse_rw_agg/RESULTS_AGG$N_agg)

RESULTS_AGG = RESULTS_AGG[!duplicated(RESULTS_AGG$plot_id),]


write.csv(RESULTS_AGG, file = paste0('~/GitLab/frodo/summaries/tables/',"smoanodiff_agg.csv"), quote = F, row.names = F)


