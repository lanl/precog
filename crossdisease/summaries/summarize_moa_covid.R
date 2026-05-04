#############################
#############################
### Results Visualization ###
#############################
#############################


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



AGG_ALL = read.csv(file = paste0('~/GitLab/frodo/summaries/tables/',"moa_covid_AGG_ALL.csv"))
AGG_ALL_OVERALL = read.csv(file = paste0('~/GitLab/frodo/summaries/tables/',"moa_covid_AGG_ALL_OVERALL.csv"))




COLORS = as.vector(palette.colors(4, palette = 'Tableau 10'))
AGG_ALL = AGG_ALL[order(AGG_ALL$prop),]
AGG_ALL$disease_chosen2 = factor(AGG_ALL$disease_chosen2, levels = as.character(AGG_ALL$disease_chosen2))
p1 = ggplot(AGG_ALL)+
  geom_bar(aes(x=disease_chosen2, y = prop), stat = 'identity', fill = COLORS[1], color = 'black')+
  theme_classic()+
  xlab('')+
  ylab('Frequency')+
  scale_y_continuous(expand = c(0,0))+
  coord_flip()
ggsave(p1, filename = paste0(my_path, 'freq.pdf'), width = 4, height = 3 )



p1 = ggplot(AGG_ALL_OVERALL)+
  geom_bar(aes(x=disease_chosen2, y = prop), stat = 'identity', fill = COLORS[1], color = 'black')+
  theme_classic()+
  xlab('')+
  ylab('Frequency')+
  scale_y_continuous(expand = c(0,0))+
  coord_flip()
ggsave(p1, filename = paste0(my_path, 'freq_in_library.pdf'), width = 4, height = 3 )



t_shift <- scales::trans_new("shift",
                             transform = function(x) {x-1},
                             inverse = function(x) {x+1})
AGG_ALL = AGG_ALL[order(AGG_ALL$prop_ratio),]
AGG_ALL$disease_chosen2 = factor(AGG_ALL$disease_chosen2, levels = as.character(AGG_ALL$disease_chosen2))
p1 = ggplot(AGG_ALL)+
  geom_bar(aes(x=disease_chosen2, y = prop_ratio), stat = 'identity', fill = COLORS[1], color = 'black')+
  theme_classic()+
  xlab('')+
  ylab('Frequency, Relative to Library')+
  scale_y_continuous(expand = c(0,0), trans = t_shift)+
  coord_flip()
ggsave(p1, filename = paste0(my_path, 'freq_ratio.pdf'), width = 4, height = 3 )











##########################
##########################
### US Overall by Week ###
##########################
##########################

# eval_path = '~/GitLab/frodo/evaluations/frodo/all/'
# FILES = list.files(paste0(eval_path))
# FILES = FILES[grep('COVID',FILES)]
# FILES = FILES[grep('_0.5',FILES)]
# GEOGRAPHIES = rep(NA, length(FILES))
# for(i in 1:length(GEOGRAPHIES)){
#   output = data.frame(data.table::fread(paste0(eval_path, FILES[i])))
#   GEOGRAPHIES[i] = output$geography[1]
# }
# index_of_us_cases = which(GEOGRAPHIES == 'United States')[3]
# FILES[index_of_us_cases]
# evalkey_number_of_us_cases = 825 ### this is the geography number for us cases!!!!!!


eval_path = '~/GitLab/frodo/evaluations/smoanodiff/covidfreq/'
FILES = list.files(paste0(eval_path))
FILES = FILES[!grepl('tycho',FILES)]
FILES = FILES[grepl('_825.',FILES)]
RESULTS = replicate(length(FILES), list(NULL))
for(i in 1:length(FILES)){
  output = data.frame(data.table::fread(paste0(eval_path, FILES[i])))
  RESULTS = output
}

eval_path2 = '~/GitLab/frodo/evaluations/smoanodiff/marginalfreq/'
FILES = list.files(paste0(eval_path2))
FILES = FILES[!grepl('tycho',FILES)]
FILES = FILES[grepl('_825.',FILES)]
RESULTS_OVERALL = replicate(length(FILES), list(NULL))
for(i in 1:length(FILES)){
  output = data.frame(data.table::fread(paste0(eval_path2, FILES[i])))
  RESULTS_OVERALL = output
}



RESULTS2 =RESULTS
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


RESULTS_OVERALL2 = RESULTS_OVERALL
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






AGG_ALL = aggregate(disease_freq~disease_chosen2+row_num+obs, FUN = sum, data = RESULTS2)
AGG_ALL = AGG_ALL %>% dplyr::group_by(row_num) %>% dplyr::mutate(prop = disease_freq/sum(disease_freq))
AGG_ALL_FLU = AGG_ALL[AGG_ALL$disease_chosen2 == 'Influenza',]
AGG_ALL = AGG_ALL[AGG_ALL$disease_chosen2 == 'COVID-19 ',]

AGG_ALL_OVERALL = aggregate(disease_freq~disease_chosen2+row_num+obs, FUN = sum, data = RESULTS_OVERALL2)
AGG_ALL_OVERALL = AGG_ALL_OVERALL %>% dplyr::group_by(row_num) %>% dplyr::mutate(prop_overall = disease_freq/sum(disease_freq))
AGG_ALL_OVERALL_FLU = AGG_ALL_OVERALL[AGG_ALL_OVERALL$disease_chosen2 == 'Influenza',]
AGG_ALL_OVERALL = AGG_ALL_OVERALL[AGG_ALL_OVERALL$disease_chosen2 == 'COVID-19 ',]



AGG_ALL = merge(AGG_ALL,AGG_ALL_OVERALL[,c('row_num', 'prop_overall')], 
                by = c('row_num'), all.x = T, all.y = F)
AGG_ALL$prop_ratio = AGG_ALL$prop/AGG_ALL$prop_overall




AGG_ALL_FLU = merge(AGG_ALL_FLU,AGG_ALL_OVERALL_FLU[,c('row_num', 'prop_overall')], 
                by = c('row_num'), all.x = T, all.y = F)
AGG_ALL_FLU$prop_ratio = AGG_ALL_FLU$prop/AGG_ALL_FLU$prop_overall




TO_PLOT = rbind(data.frame(AGG_ALL[AGG_ALL$row_num < 150,], type = 'COVID-19'),
                data.frame(AGG_ALL_FLU[AGG_ALL_FLU$row_num < 150 & AGG_ALL_FLU$row_num > min(AGG_ALL$row_num),], type = 'Influenza'))

TO_PLOT$type = factor(TO_PLOT$type, levels = rev(c('COVID-19', 'Influenza')))
p1 = ggplot(TO_PLOT)+
  geom_bar(aes(x = row_num, y = prop,  fill = type),linewidth = 1.1,  stat = "identity", position = "stack")+
  geom_line(aes(x=row_num, y = 0.5*obs/max(obs)), color = 'black')+
  geom_point(aes(x=row_num, y = 0.5*obs/max(obs)), color = 'black')+
  scale_fill_manual('Disease', breaks = c('COVID-19', 'Influenza'),values = c('deepskyblue2', 'deepskyblue4'), labels = c('COVID-19', 'Influenza'), guide = 'legend')+
  theme_classic()+
  xlab('Time Index')+
  ylab('Proportion of Neighbors')+
  scale_y_continuous(expand = c(0,0))+
  theme(legend.position = 'top')
ggsave(p1, filename = paste0(my_path, 'moa_freq3_us.pdf'), width = 8, height = 3 )



