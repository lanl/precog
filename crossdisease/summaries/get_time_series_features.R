#########################
#########################
### Get All Summaries ###
#########################
#########################




library(data.table)
library(GGally)

embedding_path_moa = '~/GitLab/frodo/features/embeddings/'
embedding_path_lstm = '~/GitLab/frodo/features/embeddings_lstm/'
feature_path_gbm = '~/GitLab/frodo/features/results/'


lut_path = '~/GitLab/frodo/lookup_tables/'
DISEASES = readxl::read_xlsx(paste0(lut_path, 'Disease_ML.xlsx'))
DISEASES = data.frame(DISEASES)
DISEASES = DISEASES$FILES

DISEASE_SOURCES = readxl::read_xlsx(paste0(lut_path, 'Evaluations_ML.xlsx'))
DISEASE_SOURCES = data.frame(DISEASE_SOURCES)
DISEASE_SOURCES = DISEASE_SOURCES$FILES
DISEASE_SOURCES = gsub('.RDS','',DISEASE_SOURCES)





#######################
### Forecastability ###
#######################

FILES_TO_RUN = list.files(paste0('~/GitLab/frodo/features/forecastability/' ))
dat = NULL
for(i in 1:length(FILES_TO_RUN)){
  dat_temp = read.csv(paste0('~/GitLab/frodo/features/forecastability/',FILES_TO_RUN[i]))
  dat = rbind(dat, data.frame(dat_temp, FILES = FILES_TO_RUN[i]))
}
dat$sample_entropy = pmin(1,pmax(0,dat$sample_entropy))

disease_lut = readxl::read_xlsx(paste0(lut_path, 'Disease_Mappings.xlsx'))
disease_lut$disease = disease_lut$Disease
disease_lut$transmission = 'other'
disease_lut$transmission[disease_lut$Transmission == 'Respiratory_Secretions'] = 'respiratory'
disease_lut$transmission[disease_lut$Transmission == 'Sexual'] = 'sexual'
disease_lut$transmission[disease_lut$Transmission == 'Vectorborne'] = 'vectorborne'
disease_lut$transmission[disease_lut$Transmission == 'Fecal-Oral' | disease_lut$Transmission == 'Water' | disease_lut$Transmission == 'Fecal-Oral/Bodily Fluids'] = 'fecaloral'

dat = merge(dat, disease_lut[,c('disease', 'transmission')], by = 'disease')


SUMMARIZE = data.frame(disease_source = DISEASE_SOURCES, 
                       disease = sub("_[^_]*$", "", DISEASE_SOURCES),
                       min = NA, mean = NA, max = NA,
                       min2 = NA, mean2 = NA, max2 = NA,
                       min3 = NA, mean3 = NA, max3 = NA,
                       min4 = NA, mean4 = NA, max4 = NA)
SUMMARIZE = merge(SUMMARIZE, disease_lut[,c('disease', 'transmission')], by = 'disease')

for(i in 1:nrow(SUMMARIZE)){
  
  ### All 
  dat_sub = dat
  SUMMARIZE$min[i] = min(dat_sub$sample_entropy, na.rm=T)
  SUMMARIZE$mean[i] = mean(dat_sub$sample_entropy, na.rm=T)
  SUMMARIZE$max[i] = max(dat_sub$sample_entropy, na.rm=T)
  
  ### Mode 
  dat_sub = dat[dat$transmission == SUMMARIZE$transmission[i], ]
  SUMMARIZE$min2[i] = min(dat_sub$sample_entropy, na.rm=T)
  SUMMARIZE$mean2[i] = mean(dat_sub$sample_entropy, na.rm=T)
  SUMMARIZE$max2[i] = max(dat_sub$sample_entropy, na.rm=T)
  
  ### Disease 
  dat_sub = dat[dat$disease == SUMMARIZE$disease[i], ]
  SUMMARIZE$min3[i] = min(dat_sub$sample_entropy, na.rm=T)
  SUMMARIZE$mean3[i] = mean(dat_sub$sample_entropy, na.rm=T)
  SUMMARIZE$max3[i] = max(dat_sub$sample_entropy, na.rm=T)
  
  ### Disease x source 
  dat_sub = dat[dat$disease_source == SUMMARIZE$disease_source[i], ]
  SUMMARIZE$min4[i] = min(dat_sub$sample_entropy, na.rm=T)
  SUMMARIZE$mean4[i] = mean(dat_sub$sample_entropy, na.rm=T)
  SUMMARIZE$max4[i] = max(dat_sub$sample_entropy, na.rm=T)
  
  
}
n = nrow(SUMMARIZE)
SUMMARIZE_LONG = data.frame(disease_source = rep(SUMMARIZE$disease_source,4),
                        disease = rep(SUMMARIZE$disease,4),
                        transmission = rep(SUMMARIZE$transmission,4),
                        min = c(SUMMARIZE$min, SUMMARIZE$min2, SUMMARIZE$min3, SUMMARIZE$min4),
                        mean = c(SUMMARIZE$mean, SUMMARIZE$mean2, SUMMARIZE$mean3, SUMMARIZE$mean4),
                        max = c(SUMMARIZE$max, SUMMARIZE$max2, SUMMARIZE$max3, SUMMARIZE$max4),
                        variation = c(rep('All Available', n),
                                      rep('Mode of Transmission', n),
                                      rep('Disease', n),
                                      rep('Disease x Source', n)))



SUMMARIZE_LONG$variation_pretty = SUMMARIZE_LONG$variation
LABELS = rev(c('Disease x Source', 'Disease','Mode of Transmission', 'All Available'))
SUMMARIZE_LONG$variation_pretty = factor(SUMMARIZE_LONG$variation_pretty, levels = LABELS)
SUMMARIZE_LONG$disease_source_pretty = SUMMARIZE_LONG$disease_source
SUMMARIZE_LONG$disease_source_pretty = gsub('_nndss','',SUMMARIZE_LONG$disease_source_pretty)
SUMMARIZE_LONG$disease_source_pretty = gsub('_usflunet',' (US FluNet)',SUMMARIZE_LONG$disease_source_pretty)
SUMMARIZE_LONG$disease_source_pretty = gsub('_whoflunet',' (WHO FluNet)',SUMMARIZE_LONG$disease_source_pretty)
SUMMARIZE_LONG$disease_source_pretty = gsub('_ushhs',' (US HHS)',SUMMARIZE_LONG$disease_source_pretty)
SUMMARIZE_LONG$disease_source_pretty = gsub('_jhuowid','-19',SUMMARIZE_LONG$disease_source_pretty)
SUMMARIZE_LONG$disease_source_pretty = gsub('_opendengue','',SUMMARIZE_LONG$disease_source_pretty)
SUMMARIZE_LONG$disease_source_pretty = gsub('_',' ',SUMMARIZE_LONG$disease_source_pretty)
SUMMARIZE_LONG$disease_source_pretty = gsub('InfluenzaA','Influenza A',SUMMARIZE_LONG$disease_source_pretty)
SUMMARIZE_LONG$disease_source_pretty = gsub('InfluenzaB','Influenza B',SUMMARIZE_LONG$disease_source_pretty)
SUMMARIZE_LONG$disease_source_pretty = factor(SUMMARIZE_LONG$disease_source_pretty, levels = rev(sort(unique(as.character(SUMMARIZE_LONG$disease_source_pretty)))))
SUMMARIZE_LONG$disease_source_pretty = as.character(SUMMARIZE_LONG$disease_source_pretty)



SUMMARIZE_LONG$transmission[SUMMARIZE_LONG$transmission == 'sexual'] = 'Sexually-\nTransmitted'
SUMMARIZE_LONG$transmission[SUMMARIZE_LONG$transmission == 'respiratory'] = 'Respiratory'
SUMMARIZE_LONG$transmission[SUMMARIZE_LONG$transmission == 'vectorborne'] = 'Vector-borne'
SUMMARIZE_LONG$transmission[SUMMARIZE_LONG$transmission == 'fecaloral'] = 'Fecal-Oral'

SUMMARIZE_LONG = SUMMARIZE_LONG[!grepl('Cocci',SUMMARIZE_LONG$disease_source),]
SUMMARIZE_LONG = SUMMARIZE_LONG[!grepl('Malaria',SUMMARIZE_LONG$disease_source),]



SAVE_FORECASTABILITY_LONG = SUMMARIZE_LONG







#######################
### GBM Predictors ###
#######################

dat = data.table::fread(file =paste0( "~/GitLab/frodo/models","/../features/results/real_train_mat.csv"))

disease_lut = readxl::read_xlsx(paste0(lut_path, 'Disease_Mappings.xlsx'))
disease_lut$disease = disease_lut$Disease
disease_lut$transmission = 'other'
disease_lut$transmission[disease_lut$Transmission == 'Respiratory_Secretions'] = 'respiratory'
disease_lut$transmission[disease_lut$Transmission == 'Sexual'] = 'sexual'
disease_lut$transmission[disease_lut$Transmission == 'Vectorborne'] = 'vectorborne'
disease_lut$transmission[disease_lut$Transmission == 'Fecal-Oral' | disease_lut$Transmission == 'Water' | disease_lut$Transmission == 'Fecal-Oral/Bodily Fluids'] = 'fecaloral'

dat = subset(dat, select = -c(transmission))
dat = merge(dat, disease_lut[,c('disease', 'transmission')], by = 'disease')


prop.table(table(dat$transmission)) #training data mode breakdown



dat = dat[dat$h == 1,]
dat = data.frame(dat)
dat$disease_source = paste0(dat$disease,'_', dat$source)

dat$abs_gradient =abs(atanh(dat$gr12_div_23)+1)

INCLUDE = c('propzerototal', 'obs', 'coefvar', 'seasonality_spearman', 'avg_recent_div_avg_global','cur_relative_max', 'abs_gradient')


SUMMARIZE = expand.grid(disease_source = DISEASE_SOURCES, variation_pretty = c('Disease x Source', 'Disease','Mode of Transmission', 'All Available'))
SUMMARIZE[,INCLUDE] = NA
SUMMARIZE = merge(SUMMARIZE, dat[!duplicated(dat$disease_source),c('disease_source', 'disease', 'transmission')], by = 'disease_source', all.x = T, all.y = F)
SUMMARIZE = SUMMARIZE[!grepl('Cocci', SUMMARIZE$disease_source),]
SUMMARIZE = SUMMARIZE[!grepl('Malaria', SUMMARIZE$disease_source),]
SUMMARIZE = SUMMARIZE[!(SUMMARIZE$disease == 'COVID' & SUMMARIZE$variation_pretty == 'Disease'),]
SUMMARIZE = SUMMARIZE[!(SUMMARIZE$disease == 'Dengue' & SUMMARIZE$variation_pretty == 'Disease'),]
SUMMARIZE$N = NA
for(i in 1:nrow(SUMMARIZE)){
  if(SUMMARIZE$variation_pretty[i] == 'Disease x Source'){
    SUMMARIZE[i,INCLUDE] = apply(dat[dat$disease_source == SUMMARIZE$disease_source[i],c(INCLUDE)],2,mean, na.rm=T)
    SUMMARIZE$N[i] = sum(dat$disease_source == SUMMARIZE$disease_source[i])
  }
  if(SUMMARIZE$variation_pretty[i] == 'Disease'){
    SUMMARIZE[i,INCLUDE] = apply(dat[dat$disease ==SUMMARIZE$disease[i],c(INCLUDE)],2,mean, na.rm=T)
    SUMMARIZE$N[i] = sum(dat$disease == SUMMARIZE$disease[i])
  }
  if(SUMMARIZE$variation_pretty[i] == 'Mode of Transmission'){
    SUMMARIZE[i,INCLUDE] = apply(dat[dat$transmission == SUMMARIZE$transmission[i] ,c(INCLUDE)],2,mean, na.rm=T)
    SUMMARIZE$N[i] = sum(dat$transmission == SUMMARIZE$transmission[i] )
  }
  if(SUMMARIZE$variation_pretty[i] == 'All Available'){
    SUMMARIZE[i,INCLUDE] = apply(data.frame(dat)[,c(INCLUDE)],2,mean, na.rm=T)
    SUMMARIZE$N[i] = nrow(dat)
  }
  print(i)
}


SUMMARIZE_LONG = SUMMARIZE
LABELS = rev(c('Disease x Source', 'Disease','Mode of Transmission', 'All Available'))
SUMMARIZE_LONG$variation_pretty = factor(SUMMARIZE_LONG$variation_pretty, levels = LABELS)
SUMMARIZE_LONG$disease_source_pretty = SUMMARIZE_LONG$disease_source
SUMMARIZE_LONG$disease_source_pretty = gsub('_nndss','',SUMMARIZE_LONG$disease_source_pretty)
SUMMARIZE_LONG$disease_source_pretty = gsub('_usflunet',' (US FluNet)',SUMMARIZE_LONG$disease_source_pretty)
SUMMARIZE_LONG$disease_source_pretty = gsub('_whoflunet',' (WHO FluNet)',SUMMARIZE_LONG$disease_source_pretty)
SUMMARIZE_LONG$disease_source_pretty = gsub('_ushhs',' (US HHS)',SUMMARIZE_LONG$disease_source_pretty)
SUMMARIZE_LONG$disease_source_pretty = gsub('_jhuowid','-19',SUMMARIZE_LONG$disease_source_pretty)
SUMMARIZE_LONG$disease_source_pretty = gsub('_opendengue','',SUMMARIZE_LONG$disease_source_pretty)
SUMMARIZE_LONG$disease_source_pretty = gsub('_',' ',SUMMARIZE_LONG$disease_source_pretty)
SUMMARIZE_LONG$disease_source_pretty = gsub('InfluenzaA','Influenza A',SUMMARIZE_LONG$disease_source_pretty)
SUMMARIZE_LONG$disease_source_pretty = gsub('InfluenzaB','Influenza B',SUMMARIZE_LONG$disease_source_pretty)
SUMMARIZE_LONG$disease_source_pretty = factor(SUMMARIZE_LONG$disease_source_pretty, levels = rev(sort(unique(as.character(SUMMARIZE_LONG$disease_source_pretty)))))
SUMMARIZE_LONG$disease_source_pretty = as.character(SUMMARIZE_LONG$disease_source_pretty)


SUMMARIZE_LONG$transmission[SUMMARIZE_LONG$transmission == 'sexual'] = 'Sexually-\nTransmitted'
SUMMARIZE_LONG$transmission[SUMMARIZE_LONG$transmission == 'respiratory'] = 'Respiratory'
SUMMARIZE_LONG$transmission[SUMMARIZE_LONG$transmission == 'vectorborne'] = 'Vector-borne'
SUMMARIZE_LONG$transmission[SUMMARIZE_LONG$transmission == 'fecaloral'] = 'Fecal-Oral'
SUMMARIZE_LONG$coefvar2 = (atanh(SUMMARIZE_LONG$coefvar)/0.1)+1

SUMMARIZE_LONG = merge(SUMMARIZE_LONG, data.frame(SAVE_FORECASTABILITY_LONG[,c('disease_source', 'variation_pretty','mean')]),
                       by = c('disease_source', 'variation_pretty'), all.x = T, all.y = F)




TO_PLOT = data.frame(SUMMARIZE_LONG[,c('disease_source_pretty', 'variation_pretty', 'transmission')], value = log10(SUMMARIZE_LONG$N), label = 'log10(N)')
TO_PLOT = rbind(TO_PLOT,data.frame(SUMMARIZE_LONG[,c('disease_source_pretty', 'variation_pretty', 'transmission')], value = SUMMARIZE_LONG$mean, label = 'Sample Entropy'))
TO_PLOT = rbind(TO_PLOT,data.frame(SUMMARIZE_LONG[,c('disease_source_pretty', 'variation_pretty', 'transmission')], value = SUMMARIZE_LONG$coefvar2, label = 'Coefficient of Variation'))

TO_PLOT$variation_pretty = as.character(TO_PLOT$variation_pretty)
TO_PLOT$variation_pretty[TO_PLOT$variation_pretty == 'Disease x Source'] = 'Single Data Stream'
LABELS = rev(c('Single Data Stream', 'Disease','Mode of Transmission', 'All Available'))
TO_PLOT$variation_pretty = factor(TO_PLOT$variation_pretty, levels = LABELS)
p1 = ggplot(TO_PLOT)+#[SUBSET$variation %in% c('Disease x Source, cold', 'All'),])+
  theme_classic()+
  theme(legend.position = 'top')+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  geom_line(aes(x=disease_source_pretty, y = value, color = variation_pretty, group = variation_pretty), linewidth = 1.5)+
  geom_point(aes(x=disease_source_pretty, y = value, fill = variation_pretty, shape = variation_pretty), size = 2)+
  xlab('')+
  ylab('Value')+
  labs(title = 'Summaries per Training Dataset (Cumulative)')+
  scale_color_manual('',values = rev(as.vector(palette.colors(4, palette = 'Tableau 10'))),
                     breaks = rev(LABELS))+
  scale_fill_manual('',values = rev(as.vector(palette.colors(4, palette = 'Tableau 10'))),
                    breaks = rev(LABELS))+
  scale_shape_manual('',values = c(21, 22, 23, 24),
                     breaks = rev(LABELS))+
  facet_grid(transmission~label, scales = 'free')+
  coord_flip()
ggsave(p1,filename = paste0('~/GitLab/frodo/summaries/',"summaries_both2.pdf"), width = 8, height = 6)




TO_PLOT = data.frame(SUMMARIZE_LONG[,c('disease_source_pretty', 'variation_pretty', 'transmission')], value = log10(SUMMARIZE_LONG$N), label = 'log10(N)')
TO_PLOT = rbind(TO_PLOT,data.frame(SUMMARIZE_LONG[,c('disease_source_pretty', 'variation_pretty', 'transmission')], value = SUMMARIZE_LONG$mean, label = 'Sample Entropy'))
TO_PLOT = rbind(TO_PLOT,data.frame(SUMMARIZE_LONG[,c('disease_source_pretty', 'variation_pretty', 'transmission')], value = SUMMARIZE_LONG$coefvar2, label = 'Coefficient of Variation'))
INCLUDE = c('COVID-19', 'Influenza (US FluNet)', 'Influenza (US HHS)', 'Influenza (WHO FluNet)', 'Dengue')
TO_PLOT = TO_PLOT[TO_PLOT$disease_source %in% INCLUDE,]

TO_PLOT$variation_pretty = as.character(TO_PLOT$variation_pretty)
TO_PLOT$variation_pretty[TO_PLOT$variation_pretty == 'Disease x Source'] = 'Single Data Stream'
LABELS = rev(c('Single Data Stream', 'Disease','Mode of Transmission', 'All Available'))
TO_PLOT$variation_pretty = factor(TO_PLOT$variation_pretty, levels = LABELS)
p1 = ggplot(TO_PLOT)+#[SUBSET$variation %in% c('Disease x Source, cold', 'All'),])+
  theme_classic()+
  theme(legend.position = 'top')+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))+
  geom_line(aes(x=disease_source_pretty, y = value, color = variation_pretty, group = variation_pretty), linewidth = 1.5)+
  geom_point(aes(x=disease_source_pretty, y = value, fill = variation_pretty, shape = variation_pretty), size = 2)+
  xlab('')+
  ylab('Value')+
  labs(title = 'Summaries per Training Dataset (Cumulative)')+
  scale_color_manual('',values = rev(as.vector(palette.colors(4, palette = 'Tableau 10'))),
                     breaks = rev(LABELS))+
  scale_fill_manual('',values = rev(as.vector(palette.colors(4, palette = 'Tableau 10'))),
                    breaks = rev(LABELS))+
  scale_shape_manual('',values = c(21, 22, 23, 24),
                     breaks = rev(LABELS))+
  facet_grid(.~label, scales = 'free')+
  coord_flip()
ggsave(p1,filename = paste0('~/GitLab/frodo/summaries/',"summaries_both2_sub.pdf"), width = 8, height = 2.5)






#########################
### Aggregate Ratios  ###
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
disease_lut$transmission = disease_lut$Transmission
SUBSET = merge(SUBSET, disease_lut[,c('disease', 'transmission')], by = 'disease', all.x = T, all.y = F)
SUBSET$transmission[SUBSET$transmission == 'Sexual'] = 'Sexually-\nTransmitted'
SUBSET$transmission[SUBSET$transmission == 'Respiratory_Secretions'] = 'Respiratory'
SUBSET$transmission[SUBSET$transmission == 'Vectorborne'] = 'Vector-borne'
SUBSET$transmission[SUBSET$transmission == 'Fecal-Oral'|SUBSET$transmission == 'Water'] = 'Fecal-Oral'



embedding_path = '~/GitLab/frodo/features/embeddings_lstm/'
dat = data.frame(data.table::fread(paste0(embedding_path,"/real_train_mat.csv")))
dat$disease_source = sub("_[^_]*$", "", dat$FILES)
TAB = table(dat$disease_source)


SUBSET = merge(SUBSET, data.frame(disease_source = names(TAB), value = as.numeric(TAB)), by = 'disease_source', all.x = T, all.y = F)





SUBSET = merge(SUBSET, SUMMARIZE_LONG[,c('disease_source', 'N', 'coefvar2', 'variation_pretty')], by = c('disease_source','variation_pretty'))




unique(SUBSET$coefvar2[SUBSET$variation_pretty == 'Disease x Source' & SUBSET$disease_source == 'Anaplasmosis_nndss'])
SUBSET = SUBSET %>% dplyr::group_by(disease_source) %>% dplyr::mutate(coefvar2_source = unique(coefvar2[variation_pretty == 'Disease x Source']), 
                                                                      N_source = unique(N[variation_pretty == 'Disease x Source']))


p1=ggplot(SUBSET[SUBSET$variation_pretty != 'Disease x Source',])+
  geom_hline(yintercept = 1)+
  #geom_point(aes(x=coefvar2/coefvar2_source,y=mae_ratio,fill = variation_pretty,  shape = variation_pretty), size = 2)+
  geom_point(aes(x=log10(coefvar2/coefvar2_source),y=mae_ratio,fill = variation_pretty,  shape = variation_pretty, size = log10(N)))+
  #geom_point(aes(x=log10(value),y=mae_ratio,fill = variation_pretty,  shape = variation_pretty), size = 2)+
  scale_fill_manual('',values = rev(as.vector(palette.colors(4, palette = 'Tableau 10'))),
                    breaks = rev(LABELS))+
  scale_shape_manual('',values = c(21, 22, 23, 24),
                     breaks = rev(LABELS))+
  #coord_cartesian(ylim=c(0.75,1.05))+
  theme_classic()+
  theme(legend.position = 'top')+
  ylab('Coefficient of Variation, Relative to Disease x Source')+
  xlab('Relative Coefficient of Variation')+
  facet_grid(.~method)+
  scale_size_continuous('',range = c(0.1,3))+
  guides(size = 'none')
ggsave(p1,filename = paste0('~/GitLab/frodo/summaries/',"mae_coefvar.pdf"), width = 6, height = 4)




SUBSET$coefvar2_ratio = SUBSET$coefvar2/SUBSET$coefvar2_source
SUBSET$log10N = log10(SUBSET$N_source)
SUBSET$mae_ratio2 = SUBSET$mae_ratio
SUBSET$mae_ratio2[SUBSET$mae_ratio2>1.05]=1.05

TEMP = SUBSET[SUBSET$variation_pretty != 'Disease x Source',]
TEMP[TEMP$mae_ratio2 == 1.05,]
library(viridis)
p1=ggplot(SUBSET[SUBSET$variation_pretty != 'Disease x Source' & SUBSET$mae_ratio2 > 1,])+# & SUBSET$method == 'MOA',])+
  geom_point(aes(x=coefvar2_ratio,y = log10N,  fill= mae_ratio2,shape = variation_pretty, group = variation_pretty), size = 2)+
  scale_shape_manual('',values = c(21, 22, 23, 24),
                     breaks = rev(LABELS))+
  theme_classic()+
  theme(legend.position = 'top')+
  xlab('Coefficient of Variation, Relative to Disease x Source')+
  ylab('log10(N)')+
  scale_fill_viridis('',  limits = c(0.8,1.05), oob = scales::oob_squish)+
  guides(shape = 'none')
ggsave(p1,filename = paste0('~/GitLab/frodo/summaries/',"mae_coefvar2.pdf"), width = 4, height = 4)





### Of time series, what percent benefitted from training using other data sources? 
TEMP = dat_long[dat_long$variation == 'All Available',]
mean(TEMP$mae <= TEMP$mae_source) # 0.7443196

TEMP = dat_long[dat_long$variation == 'Mode of Transmission',]
mean(TEMP$mae <= TEMP$mae_source) #0.805974

TEMP = dat_long[dat_long$variation == 'Disease' & dat_long$disease != 'COVID' & dat_long$disease != 'Dengue',]
mean(TEMP$mae <= TEMP$mae_source) #0.6962843













