#########################
#########################
### Get All Summaries ###
#########################
#########################




library(ggplot2)
library(data.table)
library(GGally)
library(viridis)
library(ggrepel)
theme_set(theme_classic())


### results
my_path = '~/GitLab/frodo/summaries/'
eval_path = '~/GitLab/frodo/evaluations/smoanodiff/covidfreq/'
lut_path = '~/GitLab/frodo/lookup_tables/'







#######################
#######################
### Aggregate Evals ###
#######################
#######################
FILES = list.files(paste0(eval_path))
RESULTS = replicate(length(FILES), list(NULL))
for(i in 1:length(FILES)){
  if(file.exists(paste0(eval_path, FILES[i]))){
    output = data.frame(data.table::fread(paste0(eval_path, FILES[i])))
    output = output[output$max_dist > 0,]
    if(nrow(output) > 0){
      AGG = aggregate(disease_freq~disease_chosen, FUN = sum, data = output)
      AGG$FILES = FILES[i]
      AGG$N = nrow(output)
      RESULTS[[i]] = AGG[,c('FILES', 'disease_chosen', 'disease_freq')]
    }
  }
  print(paste0('Finished: ', i, ' of ', length(FILES)))
}




RESULTS2 = do.call('rbind', RESULTS)
RESULTS2$disease_chosen2 = RESULTS2$disease_chosen
RESULTS2$disease_chosen2 = gsub('[0-9]','',RESULTS2$disease_chosen)
RESULTS2$disease_chosen2 = gsub('_nndss','',RESULTS2$disease_chosen2)
RESULTS2$disease_chosen2 = gsub('_jhuowid','-19',RESULTS2$disease_chosen2)
RESULTS2$disease_chosen2[RESULTS2$disease_chosen2 == 'HN_whoflunet_'] = 'H1N1'
RESULTS2$disease_chosen2 = gsub('_tycho_', '', RESULTS2$disease_chosen2)
RESULTS2$disease_chosen2 = gsub('_whoflunet_', '', RESULTS2$disease_chosen2)
RESULTS2$disease_chosen2 = gsub('_usflunet_', '', RESULTS2$disease_chosen2)
RESULTS2$disease_chosen2 = gsub('_who_', '', RESULTS2$disease_chosen2)
RESULTS2$disease_chosen2 = gsub('_deSouza_', '', RESULTS2$disease_chosen2)
RESULTS2$disease_chosen2 = gsub('_ushhs_', '', RESULTS2$disease_chosen2)
RESULTS2$disease_chosen2 = gsub('_aido_', '', RESULTS2$disease_chosen2)
RESULTS2$disease_chosen2 = gsub('_jhuowid_', '', RESULTS2$disease_chosen2)
RESULTS2$disease_chosen2 = gsub('_opendengue_', '', RESULTS2$disease_chosen2)
RESULTS2$disease_chosen2 = gsub('InfluenzaA','Influenza',RESULTS2$disease_chosen2)
RESULTS2$disease_chosen2 = gsub('InfluenzaB','Influenza',RESULTS2$disease_chosen2)
RESULTS2$disease_chosen2 = gsub('_',' ',RESULTS2$disease_chosen2)


write.csv(RESULTS2, file = paste0('~/GitLab/frodo/summaries/tables/',"moa_covidfreq.csv"), quote = F, row.names = F)





eval_path2 = '~/GitLab/frodo/evaluations/smoanodiff/marginalfreq/'
FILES = list.files(paste0(eval_path2))
RESULTS_OVERALL = replicate(length(FILES), list(NULL))
for(i in 1:length(FILES)){
  if(file.exists(paste0(eval_path2, FILES[i]))){
    output = data.frame(data.table::fread(paste0(eval_path2, FILES[i])))
    if(nrow(output) > 0){
      AGG = aggregate(disease_freq~disease_chosen, FUN = sum, data = output)
      AGG$FILES = FILES[i]
      AGG$N = nrow(output)
      RESULTS_OVERALL[[i]] = AGG[,c('FILES', 'disease_chosen', 'disease_freq')]
    }
  }
  print(paste0('Finished: ', i, ' of ', length(FILES)))
}



RESULTS_OVERALL2 = do.call('rbind', RESULTS_OVERALL)
RESULTS_OVERALL2$disease_chosen2 = RESULTS_OVERALL2$disease_chosen
RESULTS_OVERALL2$disease_chosen2 = gsub('[0-9]','',RESULTS_OVERALL2$disease_chosen)
RESULTS_OVERALL2$disease_chosen2 = gsub('_nndss','',RESULTS_OVERALL2$disease_chosen2)
RESULTS_OVERALL2$disease_chosen2 = gsub('_jhuowid','-19',RESULTS_OVERALL2$disease_chosen2)
RESULTS_OVERALL2$disease_chosen2[RESULTS_OVERALL2$disease_chosen2 == 'HN_whoflunet_'] = 'H1N1'
RESULTS_OVERALL2$disease_chosen2 = gsub('_tycho_', '', RESULTS_OVERALL2$disease_chosen2)
RESULTS_OVERALL2$disease_chosen2 = gsub('_whoflunet_', '', RESULTS_OVERALL2$disease_chosen2)
RESULTS_OVERALL2$disease_chosen2 = gsub('_usflunet_', '', RESULTS_OVERALL2$disease_chosen2)
RESULTS_OVERALL2$disease_chosen2 = gsub('_who_', '', RESULTS_OVERALL2$disease_chosen2)
RESULTS_OVERALL2$disease_chosen2 = gsub('_deSouza_', '', RESULTS_OVERALL2$disease_chosen2)
RESULTS_OVERALL2$disease_chosen2 = gsub('_ushhs_', '', RESULTS_OVERALL2$disease_chosen2)
RESULTS_OVERALL2$disease_chosen2 = gsub('_aido_', '', RESULTS_OVERALL2$disease_chosen2)
RESULTS_OVERALL2$disease_chosen2 = gsub('_jhuowid_', '', RESULTS_OVERALL2$disease_chosen2)
RESULTS_OVERALL2$disease_chosen2 = gsub('_opendengue_', '', RESULTS_OVERALL2$disease_chosen2)
RESULTS_OVERALL2$disease_chosen2 = gsub('InfluenzaA','Influenza',RESULTS_OVERALL2$disease_chosen2)
RESULTS_OVERALL2$disease_chosen2 = gsub('InfluenzaB','Influenza',RESULTS_OVERALL2$disease_chosen2)
RESULTS_OVERALL2$disease_chosen2 = gsub('_',' ',RESULTS_OVERALL2$disease_chosen2)



write.csv(RESULTS_OVERALL2, file = paste0('~/GitLab/frodo/summaries/tables/',"moa_marginalfreq.csv"), quote = F, row.names = F)



################################
### Second Level Aggregation ###
################################



RESULTS2 = read.csv(file = paste0('~/GitLab/frodo/summaries/tables/',"moa_covidfreq.csv"))
RESULTS_OVERALL2 = read.csv(file = paste0('~/GitLab/frodo/summaries/tables/',"moa_marginalfreq.csv"))





AGG_ALL = aggregate(disease_freq~disease_chosen2, FUN = sum, data = RESULTS2)
AGG_ALL$prop = AGG_ALL$disease_freq / sum(AGG_ALL$disease_freq)
AGG_ALL = AGG_ALL[AGG_ALL$prop > 0.01,]
AGG_ALL = AGG_ALL[order(AGG_ALL$prop),]
AGG_ALL$disease_chosen2 = factor(AGG_ALL$disease_chosen2, levels = AGG_ALL$disease_chosen2)





AGG_ALL_OVERALL = aggregate(disease_freq~disease_chosen2, FUN = sum, data = RESULTS_OVERALL2)
AGG_ALL_OVERALL$prop = AGG_ALL_OVERALL$disease_freq / sum(AGG_ALL_OVERALL$disease_freq)
AGG_ALL_OVERALL = AGG_ALL_OVERALL[AGG_ALL_OVERALL$disease_chosen2 %in% AGG_ALL$disease_chosen2,]

AGG_ALL_OVERALL = AGG_ALL_OVERALL[order(AGG_ALL_OVERALL$prop),]
AGG_ALL_OVERALL$disease_chosen2 = factor(AGG_ALL_OVERALL$disease_chosen2, levels = AGG_ALL_OVERALL$disease_chosen2)





AGG_ALL = merge(AGG_ALL, data.frame(disease_chosen2 = AGG_ALL_OVERALL$disease_chosen2, prop_marginal = AGG_ALL_OVERALL$prop),
                by = c('disease_chosen2'), all.x = T, all.y = F)
AGG_ALL$prop_ratio = AGG_ALL$prop/AGG_ALL$prop_marginal
AGG_ALL = AGG_ALL[order(AGG_ALL$prop_ratio),]
AGG_ALL$disease_chosen2 = factor(AGG_ALL$disease_chosen2, levels = AGG_ALL$disease_chosen2)





write.csv(AGG_ALL, file = paste0('~/GitLab/frodo/summaries/tables/',"moa_covid_AGG_ALL.csv"), quote = F, row.names = F)
write.csv(AGG_ALL_OVERALL, file = paste0('~/GitLab/frodo/summaries/tables/',"moa_covid_AGG_ALL_OVERALL.csv"), quote = F, row.names = F)







###############################
### By Mode of Transmission ###
###############################



disease_lut = readxl::read_xlsx(paste0(lut_path, 'Disease_Mappings.xlsx'))
disease_lut$disease_chosen2 = disease_lut$Disease
disease_lut$transmission = 'other'
disease_lut$transmission[disease_lut$Transmission == 'Respiratory_Secretions'] = 'respiratory'
disease_lut$transmission[disease_lut$Transmission == 'Sexual'] = 'sexual'
disease_lut$transmission[disease_lut$Transmission == 'Vectorborne'] = 'vectorborne'
disease_lut$transmission[disease_lut$Transmission == 'Fecal-Oral' | disease_lut$Transmission == 'Water' | disease_lut$Transmission == 'Fecal-Oral/Bodily Fluids'] = 'fecaloral'

AGG_ALL = aggregate(disease_freq~disease_chosen2, FUN = sum, data = RESULTS2)
AGG_ALL = merge(AGG_ALL, disease_lut[,c('disease_chosen2', 'transmission')], by = 'disease_chosen2')
AGG_ALL2 = aggregate(disease_freq~transmission, FUN = sum, data = AGG_ALL)
AGG_ALL2$prop = AGG_ALL2$disease_freq / sum(AGG_ALL2$disease_freq)




AGG_ALL_OVERALL = aggregate(disease_freq~disease_chosen2, FUN = sum, data = RESULTS_OVERALL2)
AGG_ALL_OVERALL = merge(AGG_ALL_OVERALL, disease_lut[,c('disease_chosen2', 'transmission')], by = 'disease_chosen2')
AGG_ALL_OVERALL2 = aggregate(disease_freq~transmission, FUN = sum, data = AGG_ALL_OVERALL)
AGG_ALL_OVERALL2$prop = AGG_ALL_OVERALL2$disease_freq / sum(AGG_ALL_OVERALL2$disease_freq)

AGG_ALL_OVERALL$prop = AGG_ALL_OVERALL$disease_freq / sum(AGG_ALL_OVERALL$disease_freq)
