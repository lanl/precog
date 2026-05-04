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
library(dplyr)


#########################
### Overall MAE Plots ###
#########################

frodo = read.csv(paste0('~/GitLab/frodo/summaries/tables/',"frodo.csv"), header = T)
lstm = read.csv(paste0('~/GitLab/frodo/summaries/tables/',"lstm.csv"), header = T)
smoa = read.csv(paste0('~/GitLab/frodo/summaries/tables/',"smoanodiff.csv"), header = T)

n = nrow(frodo)
frodo_long = data.frame(disease_source = rep(frodo$disease_source,4),
                        disease = rep(frodo$disease,4),
                        mae = c(frodo$mae, frodo$mae2, frodo$mae3, frodo$mae4),
                        mae_source = rep(frodo$mae4, 4),
                        mae_rw = rep(frodo$mae_rw, 4),
                        method = 'GBM',
                        weight = rep(frodo$N, 4),
                        variation = c(rep('All Available', n),
                                      rep('Mode of Transmission', n),
                                      rep('Disease', n),
                                      rep('Disease x Source', n)))

n = nrow(lstm)
lstm_long = data.frame(disease_source = rep(lstm$disease_source,4),
                       disease = rep(lstm$disease,4),
                       mae = c(lstm$mae, lstm$mae2, lstm$mae3, lstm$mae4),
                       mae_source = rep(lstm$mae4, 4),
                       mae_rw = rep(lstm$mae_rw, 4),
                       method = 'LSTM',
                       weight = rep(lstm$N, 4),
                       variation = c(rep('All Available', n),
                                     rep('Mode of Transmission', n),
                                     rep('Disease', n),
                                     rep('Disease x Source', n)))


n = nrow(smoa)
smoa_long = data.frame(disease_source = rep(smoa$disease_source,4),
                       disease = rep(smoa$disease,4),
                       mae = c(smoa$mae, smoa$mae2, smoa$mae3, smoa$mae4),
                       mae_source = rep(smoa$mae4, 4),
                       mae_rw = rep(smoa$mae_rw, 4),
                       method = 'MOA',
                       weight = rep(smoa$N, 4),
                       variation = c(rep('All Available', n),
                                     rep('Mode of Transmission', n),
                                     rep('Disease', n),
                                     rep('Disease x Source', n)))

dat_long = rbind(frodo_long, lstm_long, smoa_long)
dat_long = dat_long[!grepl('InfluenzaA',dat_long$disease_source),]
dat_long = dat_long[!grepl('InfluenzaB',dat_long$disease_source),]
dat_long = dat_long[dat_long$disease_source != 'Dengue_tycho',]

#####################################
### Summarize Overall Performance ###
#####################################

OTHER = c(as.numeric(frodo$mae4 < apply(cbind(frodo$mae, frodo$mae2, frodo$mae3), 1, min))[frodo$disease != 'COVID' & frodo$disease != 'Dengue'],
                 as.numeric(lstm$mae4 < apply(cbind(lstm$mae, lstm$mae2, lstm$mae3), 1, min))[lstm$disease != 'COVID' & lstm$disease != 'Dengue'],
                 as.numeric(smoa$mae4 < apply(cbind(smoa$mae, smoa$mae2, smoa$mae3), 1, min))[smoa$disease != 'COVID' & smoa$disease != 'Dengue'])
mean(OTHER)

#mae3 = mae4 for covid and dengue 
COVID_DENGUE = c(as.numeric(frodo$mae4 < apply(cbind(frodo$mae, frodo$mae2), 1, min))[frodo$disease == 'COVID' | frodo$disease == 'Dengue'],
                 as.numeric(lstm$mae4 < apply(cbind(lstm$mae, lstm$mae2), 1, min))[lstm$disease == 'COVID' | lstm$disease == 'Dengue'],
                 as.numeric(smoa$mae4 < apply(cbind(smoa$mae, smoa$mae2), 1, min))[smoa$disease == 'COVID' | smoa$disease == 'Dengue'])
mean(COVID_DENGUE)

mean(c(OTHER,COVID_DENGUE ))


SUB = dat_long[dat_long$variation == 'All Available',]
aggregate(as.numeric(mae < mae_source)~method, FUN = mean, data = SUB)
aggregate(as.numeric(mae < mae_source)~1, FUN = mean, data = SUB)

SUB = dat_long[dat_long$variation == 'Mode of Transmission',]
aggregate(as.numeric(mae < mae_source)~method, FUN = mean, data = SUB)
aggregate(as.numeric(mae < mae_source)~1, FUN = mean, data = SUB)

SUB = dat_long[dat_long$variation == 'Disease' & dat_long$disease != 'COVID' & dat_long$disease != 'Dengue',]
aggregate(as.numeric(mae < mae_source)~method, FUN = mean, data = SUB)
aggregate(as.numeric(mae < mae_source)~1, FUN = mean, data = SUB)



####################################
### Plotting Overall Performance ###
####################################


SUBSET = dat_long
SUBSET$variation_pretty = SUBSET$variation
LABELS = rev(c('Disease x Source', 'Disease','Mode of Transmission', 'All Available'))
SUBSET$variation_pretty = factor(SUBSET$variation_pretty, levels = LABELS)
SUBSET2 = SUBSET[(SUBSET$disease == 'COVID' & SUBSET$variation_pretty != 'Disease') | SUBSET$disease != 'COVID',]
SUBSET2 = SUBSET2[SUBSET2$variation_pretty != 'Disease x Source',]
SUBSET2 = SUBSET2[!(SUBSET2$disease == 'Dengue' & SUBSET2$variation_pretty == 'Disease'),]
p1 = ggplot(SUBSET2)+
  geom_hline(yintercept = 1, color = 'gray', linetype = 2)+
  geom_point(aes(y=mae/mae_source, x = variation_pretty, fill = variation_pretty), shape = 21, color = 'black',alpha = 0.3, position = position_jitterdodge())+
  #geom_boxplot(aes(y=mae/mae_source, x = variation_pretty, fill = variation_pretty, weight = weight), alpha = 1, outlier.shape = NA)+
  geom_boxplot(aes(y=mae/mae_source, x = variation_pretty, fill = variation_pretty), alpha = 1, outlier.shape = NA)+
  xlab('')+
  ylab('MAE, Relative to Single Data Stream')+
  labs(title = 'MAE, Relative to Disease x Source')+
  facet_grid(ifelse(disease == 'COVID', 'COVID-19',ifelse(disease == 'Dengue','Dengue','Other') )~method, scales=  'free')+
  #facet_grid(ifelse(disease == 'COVID', 'COVID-19','Other' )~method, scales=  'free')+
  scale_color_manual(values = as.vector(palette.colors(4, palette = 'Tableau 10')),
                     breaks = LABELS)+
  scale_fill_manual(values = as.vector(palette.colors(4, palette = 'Tableau 10')),
                    breaks = LABELS)+
  theme_classic()+
  theme(legend.position = 'top')+
  coord_flip(ylim=c(0.65,1.3))+
  guides(color = 'none', fill = 'none')+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
ggsave(p1,filename = paste0('~/GitLab/frodo/summaries/',"mae_rw_noaggBOTH_source.pdf"), width = 10, height = 4)









### Separating out by mode of transmission
lut_path = '~/GitLab/frodo/lookup_tables/'
disease_lut = readxl::read_xlsx(paste0(lut_path, 'Disease_Mappings.xlsx'))
SUBSET2 = merge(SUBSET, data.frame(disease = disease_lut$Disease, transmission = disease_lut$Transmission ),by = 'disease', all.x = T, all.y = F)
SUBSET2$transmission_cat = NA
SUBSET2$transmission_cat[SUBSET2$transmission == 'Respiratory_Secretions'] = 'respiratory'
SUBSET2$transmission_cat[SUBSET2$transmission == 'Vectorborne'] = 'vectorborne'
SUBSET2$transmission_cat[SUBSET2$transmission == 'Sexual'] = 'sexual'
SUBSET2$transmission_cat[SUBSET2$transmission == 'Fecal-Oral' | SUBSET2$transmission == 'Water' | SUBSET2$transmission == 'Fecal-Oral/Bodily Fluids'] = 'fecaloral'
SUBSET2 = SUBSET2[SUBSET2$variation_pretty != 'Disease x Source',]

SUBSET2$transmission_cat2 = SUBSET2$transmission_cat
SUBSET2$transmission_cat2[SUBSET2$transmission_cat == 'respiratory' & SUBSET2$disease == 'COVID'] = 'Respiratory\n(COVID-19)'
SUBSET2$transmission_cat2[SUBSET2$transmission_cat == 'respiratory' & SUBSET2$disease != 'COVID'] = 'Respiratory\n(other)'
SUBSET2$transmission_cat2[SUBSET2$transmission_cat == 'vectorborne' & SUBSET2$disease == 'Dengue'] = 'Vector-borne\n(Dengue)'
SUBSET2$transmission_cat2[SUBSET2$transmission_cat == 'vectorborne' & SUBSET2$disease != 'Dengue'] = 'Vector-borne\n(Other)'
SUBSET2$transmission_cat2[SUBSET2$transmission_cat == 'sexual'] = 'Sexually-\nTransmitted'
SUBSET2$transmission_cat2[SUBSET2$transmission_cat == 'fecaloral'] = 'Fecal-Oral'

SUBSET2 = SUBSET2[!(SUBSET2$variation_pretty == 'Disease' & SUBSET2$disease == 'COVID'),]
SUBSET2 = SUBSET2[!(SUBSET2$variation_pretty == 'Disease' & SUBSET2$disease == 'Dengue'),]

p1 = ggplot(SUBSET2)+
  geom_hline(yintercept = 1, color = 'gray', linetype = 2)+
  geom_point(aes(y=mae/mae_source, x = variation_pretty, fill = variation_pretty), shape = 21, color = 'black',alpha = 0.3, position = position_jitterdodge())+
  geom_boxplot(aes(y=mae/mae_source, x = variation_pretty, fill = variation_pretty), alpha = 1, outlier.shape = NA)+
  xlab('')+
  ylab('MAE, Relative to Single Data Stream')+
  labs(title = 'MAE, Relative to Disease x Source')+
  facet_grid(transmission_cat2~method, scales=  'free')+
  scale_color_manual(values = as.vector(palette.colors(4, palette = 'Tableau 10')),
                     breaks = LABELS)+
  scale_fill_manual(values = as.vector(palette.colors(4, palette = 'Tableau 10')),
                    breaks = LABELS)+
  theme_classic()+
  theme(legend.position = 'top')+
  coord_flip(ylim=c(0.65,1.3))+
  guides(color = 'none', fill = 'none')+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
ggsave(p1,filename = paste0('~/GitLab/frodo/summaries/',"mae_rw_bymode_source.pdf"), width = 10, height = 6)


SUBSET3 = SUBSET2[SUBSET2$variation_pretty == 'All Available', ]
table(SUBSET3$transmission_cat)
prop.table(table(SUBSET3$transmission_cat))






#########################
### Aggregate Ratios  ###
#########################

frodo = read.csv(paste0('~/GitLab/frodo/summaries/tables/',"frodo.csv"), header = T)
lstm = read.csv(paste0('~/GitLab/frodo/summaries/tables/',"lstm.csv"), header = T)
smoa = read.csv(paste0('~/GitLab/frodo/summaries/tables/',"smoanodiff.csv"), header = T)


sum(smoa$N*4)+sum(frodo$N*4)+sum(lstm$N*4) #total number of forecasts performed

n = nrow(frodo)
frodo_long = data.frame(disease_source = rep(frodo$disease_source,4),
                        disease = rep(frodo$disease,4),
                        mae = c(frodo$mae, frodo$mae2, frodo$mae3, frodo$mae4),
                        mae_source = rep(frodo$mae4, 4),
                        mae_rw = rep(frodo$mae_rw, 4),
                        method = 'GBM',
                        weight = rep(frodo$N, 4),
                        variation = c(rep('All Available', n),
                                      rep('Mode of Transmission', n),
                                      rep('Disease', n),
                                      rep('Disease x Source', n)))

n = nrow(lstm)
lstm_long = data.frame(disease_source = rep(lstm$disease_source,4),
                       disease = rep(lstm$disease,4),
                       mae = c(lstm$mae, lstm$mae2, lstm$mae3, lstm$mae4),
                       mae_source = rep(lstm$mae4, 4),
                       mae_rw = rep(lstm$mae_rw, 4),
                       method = 'LSTM',
                       weight = rep(lstm$N, 4),
                       variation = c(rep('All Available', n),
                                     rep('Mode of Transmission', n),
                                     rep('Disease', n),
                                     rep('Disease x Source', n)))


n = nrow(smoa)
smoa_long = data.frame(disease_source = rep(smoa$disease_source,4),
                       disease = rep(smoa$disease,4),
                       mae = c(smoa$mae, smoa$mae2, smoa$mae3, smoa$mae4),
                       mae_source = rep(smoa$mae4, 4),
                       mae_rw = rep(smoa$mae_rw, 4),
                       method = 'MOA',
                       weight = rep(smoa$N, 4),
                       variation = c(rep('All Available', n),
                                     rep('Mode of Transmission', n),
                                     rep('Disease', n),
                                     rep('Disease x Source', n)))

dat_long = rbind(frodo_long, lstm_long, smoa_long)
dat_long = dat_long[!grepl('InfluenzaA',dat_long$disease_source),]
dat_long = dat_long[!grepl('InfluenzaB',dat_long$disease_source),]
dat_long = dat_long[dat_long$disease_source != 'Dengue_tycho',]
summary(dat_long$mae/dat_long$mae_source)


dat_long[dat_long$mae_source == 0,]


SUBSET = dat_long %>% dplyr::group_by(disease_source, method, variation) %>% dplyr::mutate(mae_ratio = mean(mae/mae_source, na.rm=T))
SUBSET = SUBSET[!duplicated(paste0(SUBSET$disease_source, '_' ,SUBSET$method, '_', SUBSET$variation)),]
SUBSET$variation_pretty = SUBSET$variation
LABELS = rev(c('Disease x Source', 'Disease','Mode of Transmission', 'All Available'))
SUBSET$variation_pretty = factor(SUBSET$variation_pretty, levels = LABELS)
SUBSET$disease_source_pretty = SUBSET$disease_source
SUBSET$disease_source_pretty = gsub('_nndss','',SUBSET$disease_source_pretty)
SUBSET$disease_source_pretty = gsub('_usflunet',' (US FluNet)',SUBSET$disease_source_pretty)
SUBSET$disease_source_pretty = gsub('_whoflunet',' (WHO FluNet)',SUBSET$disease_source_pretty)
SUBSET$disease_source_pretty = gsub('_ushhs',' (US HHS)',SUBSET$disease_source_pretty)
SUBSET$disease_source_pretty = gsub('_jhuowid','-19',SUBSET$disease_source_pretty)
SUBSET$disease_source_pretty = gsub('_opendengue','',SUBSET$disease_source_pretty)
SUBSET$disease_source_pretty = gsub('_',' ',SUBSET$disease_source_pretty)
SUBSET$disease_source_pretty = gsub('InfluenzaA','Influenza A',SUBSET$disease_source_pretty)
SUBSET$disease_source_pretty = gsub('InfluenzaB','Influenza B',SUBSET$disease_source_pretty)
SUBSET$disease_source_pretty = factor(SUBSET$disease_source_pretty, levels = rev(sort(unique(as.character(SUBSET$disease_source_pretty)))))

SUBSET$disease_source_pretty = as.character(SUBSET$disease_source_pretty)

lut_path = '~/GitLab/frodo/lookup_tables/'
disease_lut = readxl::read_xlsx(paste0(lut_path, 'Disease_Mappings.xlsx'))
disease_lut$disease = disease_lut$Disease
SUBSET = merge(SUBSET, disease_lut[,c('disease', 'Transmission')], by = 'disease', all.x = T, all.y = F)
SUBSET$Transmission[SUBSET$Transmission == 'Sexual'] = 'Sexually-\nTransmitted'
SUBSET$Transmission[SUBSET$Transmission == 'Respiratory_Secretions'] = 'Respiratory'
SUBSET$Transmission[SUBSET$Transmission == 'Vectorborne'] = 'Vector-borne'
SUBSET$Transmission[SUBSET$Transmission == 'Fecal-Oral'|SUBSET$Transmission == 'Water'] = 'Fecal-Oral'



embedding_path = '~/GitLab/frodo/features/embeddings_lstm/'
dat = data.frame(data.table::fread(paste0(embedding_path,"/real_train_mat.csv")))
dat$disease_source = sub("_[^_]*$", "", dat$FILES)
TAB = table(dat$disease_source)


SUBSET = merge(SUBSET, data.frame(disease_source = names(TAB), value = as.numeric(TAB)), by = 'disease_source', all.x = T, all.y = F)

library(ggh4x)
A= SUBSET[!duplicated(paste0(SUBSET$disease_source_pretty,'_', SUBSET$method)),]
A = A[A$method == 'MOA',]

TO_PLOT = SUBSET[SUBSET$variation_pretty != 'Disease x Source',]
TO_PLOT = TO_PLOT[!(TO_PLOT$disease == 'COVID' & TO_PLOT$variation_pretty == 'Disease'),]
TO_PLOT = TO_PLOT[!(TO_PLOT$disease == 'Dengue' & TO_PLOT$variation_pretty == 'Disease'),]

p1 = ggplot(TO_PLOT)+#[SUBSET$variation %in% c('Disease x Source, cold', 'All'),])+
  theme_classic()+
  theme(legend.position = 'top')+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  geom_hline(yintercept = 1, color = 'gray', linetype = 1, linewidth = 1)+
  geom_line(aes(x=disease_source_pretty, y = mae_ratio, color = variation_pretty, group = variation_pretty), linewidth = 1.5)+
  #geom_point(aes(x=disease_source, y = mae/mae_rw, group = variation_pretty, shape = variation_pretty), size = 2.5, color = 'black')+
  geom_point(aes(x=disease_source_pretty, y = mae_ratio, fill = variation_pretty, shape = variation_pretty), size = 2)+
  xlab('')+
  ylab('MAE, Relative to Single Data Stream')+
  labs(title = 'MAE / Disease x Source MAE')+
  scale_color_manual('',values = rev(as.vector(palette.colors(4, palette = 'Tableau 10'))),
                     breaks = rev(LABELS))+
  scale_fill_manual('',values = rev(as.vector(palette.colors(4, palette = 'Tableau 10'))),
                    breaks = rev(LABELS))+
  scale_shape_manual('',values = c(21, 22, 23, 24),
                     breaks = rev(LABELS))+
  # scale_color_brewer('',palette = 'Spectral')+
  facet_grid(Transmission~method, scales = 'free_y')+
  # geom_text(aes(y=1.2, x = disease_source_pretty, label = value), size = 3, hjust = 1,
  #           data = A)+
  coord_flip(ylim=c(0.8,1.2))
ggsave(p1,filename = paste0('~/GitLab/frodo/summaries/',"mae_ratio_new.pdf"), width = 8, height = 6)



library(ggh4x)


TO_PLOT = SUBSET[SUBSET$variation_pretty != 'Disease x Source',]
TO_PLOT = TO_PLOT[!(TO_PLOT$disease == 'COVID' & TO_PLOT$variation_pretty == 'Disease'),]
TO_PLOT = TO_PLOT[!(TO_PLOT$disease == 'Dengue' & TO_PLOT$variation_pretty == 'Disease'),]
INCLUDE = c('COVID_jhuowid', 'Influenza_usflunet', 'Influenza_ushhs', 'Influenza_whoflunet', 'Dengue_opendengue')
TO_PLOT = TO_PLOT[TO_PLOT$disease_source %in% INCLUDE,]

p1 = ggplot(TO_PLOT)+#[SUBSET$variation %in% c('Disease x Source, cold', 'All'),])+
  theme_classic()+
  theme(legend.position = 'top')+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  geom_hline(yintercept = 1, color = 'gray', linetype = 1, linewidth = 1)+
  geom_line(aes(x=disease_source_pretty, y = mae_ratio, color = variation_pretty, group = variation_pretty), linewidth = 1.5)+
  #geom_point(aes(x=disease_source, y = mae/mae_rw, group = variation_pretty, shape = variation_pretty), size = 2.5, color = 'black')+
  geom_point(aes(x=disease_source_pretty, y = mae_ratio, fill = variation_pretty, shape = variation_pretty), size = 2)+
  xlab('')+
  ylab('MAE, Relative to Single Data Stream')+
  labs(title = 'MAE / Disease x Source MAE')+
  scale_color_manual('',values = rev(as.vector(palette.colors(4, palette = 'Tableau 10'))),
                     breaks = rev(LABELS))+
  scale_fill_manual('',values = rev(as.vector(palette.colors(4, palette = 'Tableau 10'))),
                    breaks = rev(LABELS))+
  scale_shape_manual('',values = c(21, 22, 23, 24),
                     breaks = rev(LABELS))+
  # scale_color_brewer('',palette = 'Spectral')+
  facet_grid(.~method, scales = 'free_y')+
  # geom_text(aes(y=1.2, x = disease_source_pretty, label = value), size = 3, hjust = 1,
  #           data = A)+
  coord_flip(ylim=c(0.8,1.2))
ggsave(p1,filename = paste0('~/GitLab/frodo/summaries/',"mae_ratio_new_subset.pdf"), width = 8, height = 2.5)






TO_PLOT = SUBSET[SUBSET$value < 20000 & SUBSET$variation_pretty != 'Disease x Source',]
TO_PLOT$mae_ratio[TO_PLOT$mae_ratio>1.05] = 1.05
p1=ggplot(TO_PLOT)+
  geom_hline(yintercept = 1)+
  geom_point(aes(x=value,y=mae_ratio,fill = variation_pretty,  shape = variation_pretty), size = 2)+
  #geom_point(aes(x=log10(value),y=mae_ratio,fill = variation_pretty,  shape = variation_pretty), size = 2)+
  scale_fill_manual('',values = rev(as.vector(palette.colors(4, palette = 'Tableau 10'))),
                     breaks = rev(LABELS))+
  scale_shape_manual('',values = c(21, 22, 23, 24),
                     breaks = rev(LABELS))+
  coord_cartesian(ylim=c(0.75,1.05))+
  theme_classic()+
  theme(legend.position = 'top')+
  ylab('MAE, Relative to Single Data Stream')+
  xlab('Number of Single Data Stream Embeddings')
ggsave(p1,filename = paste0('~/GitLab/frodo/summaries/',"mae_ratio_scatter.pdf"), width = 6, height = 4)









###################################
### Forecastability Exploration ###
###################################

frodo = read.csv(paste0('~/GitLab/frodo/summaries/tables/',"frodo.csv"), header = T)
lstm = read.csv(paste0('~/GitLab/frodo/summaries/tables/',"lstm.csv"), header = T)
smoa = read.csv(paste0('~/GitLab/frodo/summaries/tables/',"smoanodiff.csv"), header = T)

frodo = frodo[!grepl('InfluenzaA',frodo$disease_source),]
frodo = frodo[!grepl('InfluenzaB',frodo$disease_source),]
frodo = frodo[!grepl('tycho',frodo$disease_source),]
frodo$eval_key_temp = gsub('_0.5.csv','',gsub('real_eval_mat_','',frodo$FILES))

frodo$eval_key = NA
for(i in 1:length(frodo[,1])){
  SPLIT = strsplit(frodo$eval_key_temp[i], split = '_')
  frodo$eval_key[i] = paste0(paste0(SPLIT[[1]][-1], collapse = '_'), '_',SPLIT[[1]][1])
}

lstm = lstm[!grepl('InfluenzaA',lstm$disease_source),]
lstm = lstm[!grepl('InfluenzaB',lstm$disease_source),]
lstm = lstm[!grepl('tycho',lstm$disease_source),]
lstm$eval_key = gsub('_0.5.csv','',gsub('real_eval_mat_','',lstm$FILES))
lstm$eval_key = gsub('.csv','',lstm$eval_key)

smoa = smoa[!grepl('InfluenzaA',smoa$disease_source),]
smoa = smoa[!grepl('InfluenzaB',smoa$disease_source),]
smoa = smoa[!grepl('tycho',smoa$disease_source),]
smoa$eval_key = gsub('_0.5.csv','',gsub('real_eval_mat_','',smoa$FILES))
smoa$eval_key = gsub('.csv','',smoa$eval_key)


FILES_TO_RUN = gsub('_0.5.csv','',gsub('real_eval_mat_','',frodo$FILES))

dat = NULL
for(i in 1:length(FILES_TO_RUN)){
  dat_temp = read.csv(paste0('~/GitLab/frodo/features/forecastability/real_mat_',FILES_TO_RUN[i],'.csv' ))
  dat = rbind(dat, data.frame(dat_temp, eval_key = frodo$eval_key[i]))
}




n = nrow(frodo)
frodo_long = data.frame(disease_source = rep(frodo$disease_source,4),
                        disease = rep(frodo$disease,4),
                        eval_key = rep(frodo$eval_key,4),
                        mae = c(frodo$mae, frodo$mae2, frodo$mae3, frodo$mae4),
                        mae_rw = rep(frodo$mae_rw, 4),
                        mae_source = rep(frodo$mae4, 4),
                        method = 'GBM',
                        variation = c(rep('All Available', n),
                                      rep('Mode of Transmission', n),
                                      rep('Disease', n),
                                      rep('Disease x Source', n)))

n = nrow(lstm)
lstm_long = data.frame(disease_source = rep(lstm$disease_source,4),
                       disease = rep(lstm$disease,4),
                       eval_key = rep(lstm$eval_key,4),
                       mae = c(lstm$mae, lstm$mae2, lstm$mae3, lstm$mae4),
                       mae_rw = rep(lstm$mae_rw, 4),
                       mae_source = rep(lstm$mae4, 4),
                       method = 'LSTM',
                       variation = c(rep('All Available', n),
                                     rep('Mode of Transmission', n),
                                     rep('Disease', n),
                                     rep('Disease x Source', n)))


n = nrow(smoa)
smoa_long = data.frame(disease_source = rep(smoa$disease_source,4),
                       disease = rep(smoa$disease,4),
                       eval_key = rep(smoa$eval_key,4),
                       mae = c(smoa$mae, smoa$mae2, smoa$mae3, smoa$mae4),
                       mae_rw = rep(smoa$mae_rw, 4),
                       mae_source = rep(smoa$mae4, 4),
                       method = 'MOA',
                       variation = c(rep('All Available', n),
                                     rep('Mode of Transmission', n),
                                     rep('Disease', n),
                                     rep('Disease x Source', n)))

dat_long = rbind(frodo_long, lstm_long, smoa_long)
dat_long = merge(dat_long, dat[,c('eval_key', 'hurst', 'hurst_empirical','sample_entropy','approx_entropy')],
                 by = c('eval_key'), all.x = T, all.y = F)

dat_long$hurst = pmin(1,pmax(0,dat_long$hurst))
dat_long$hurst_empirical = pmin(1,pmax(0,dat_long$hurst_empirical))
dat_long$sample_entropy = pmin(1,pmax(0,dat_long$sample_entropy))
dat_long$approx_entropy = pmin(1,pmax(0,dat_long$approx_entropy))
dat_long = dat_long[!is.na(dat_long$sample_entropy),]
dat_long$entropy_cat = cut(dat_long$sample_entropy, breaks = c(0, 0.33, 0.67,1))
dat_long$entropy_cat2 = cut(dat_long$sample_entropy, breaks = seq(0,1,0.1))




SUBSET = dat_long
SUBSET$variation_pretty = SUBSET$variation
SUBSET$variation_pretty = factor(SUBSET$variation_pretty, levels = LABELS)
SUBSET = SUBSET[SUBSET$variation_pretty != 'Disease x Source',]
SUBSET = SUBSET[SUBSET$disease != 'COVID' | (SUBSET$disease == 'COVID' & SUBSET$variation_pretty != 'Disease'),]
SUBSET = SUBSET[SUBSET$disease != 'Dengue' | (SUBSET$disease == 'Dengue' & SUBSET$variation_pretty != 'Disease'),]
SUBSET$mae_ratio= SUBSET$mae/SUBSET$mae_source
MEDIANS = aggregate(mae_ratio~entropy_cat2+variation_pretty, FUN = median, data = SUBSET)


p1 = ggplot(SUBSET)+
  geom_hline(yintercept = 1, color = 'gray', linetype = 2)+
  geom_boxplot(aes(y=mae/mae_source, x = factor(entropy_cat2), fill = variation_pretty, group = factor(entropy_cat2)))+
  geom_point(aes(y=mae_ratio, x=factor(entropy_cat2), group = factor(entropy_cat2)), size = 2, data = MEDIANS)+
  xlab('Sample Entropy (Higher = Harder to Forecast)')+
  ylab('MAE, Relative to Single Data Stream')+
  labs(title = 'MAE, Relative to Disease x Source')+
  #facet_grid(variation_pretty~method, scales=  'free_x')+
  facet_grid(.~factor(variation_pretty, levels = c('Disease', 'Mode of Transmission', 'All Available')), scales=  'free_y')+
  scale_color_manual(values = as.vector(palette.colors(4, palette = 'Tableau 10')),
                     breaks = LABELS)+
  scale_fill_manual(values = as.vector(palette.colors(4, palette = 'Tableau 10')),
                    breaks = LABELS)+
  theme_classic()+
  theme(legend.position = 'top')+
  coord_cartesian(ylim=c(0.65,1.3))+
  guides(color = 'none', fill = 'none')+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  theme(axis.text.x = element_text(angle = -90, vjust = 0, hjust = 0.5))
ggsave(p1,filename = paste0('~/GitLab/frodo/summaries/',"mae_forecastability3.pdf"), width = 8, height = 4)





