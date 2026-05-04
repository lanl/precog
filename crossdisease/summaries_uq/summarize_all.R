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



#######################
### Read in Results ###
#######################

frodo = read.csv(paste0('~/GitLab/frodo/summaries_uq/tables/',"frodo.csv"), header = T)
lstm = read.csv(paste0('~/GitLab/frodo/summaries_uq/tables/',"lstm.csv"), header = T)
smoa = read.csv(paste0('~/GitLab/frodo/summaries_uq/tables/',"smoanodiff.csv"), header = T)

frodo = frodo[!is.na(frodo$N),]
lstm = lstm[!is.na(lstm$N),]
smoa = smoa[!is.na(smoa$N),]

n = nrow(frodo)
frodo_long = data.frame(disease_source = rep(frodo$disease_source,4),
                        #disease = rep(frodo$disease,4),
                        wis = c(frodo$wis_1, frodo$wis_2, frodo$wis_3, frodo$wis_4),
                        wis_source = rep(frodo$wis_4,4),
                        
                        wis_dispersion = c(frodo$wis_dispersion_1, frodo$wis_dispersion_2, frodo$wis_dispersion_3, frodo$wis_dispersion_4),
                        wis_dispersion_source = rep(frodo$wis_dispersion_4,4),
                        wis_overprediction = c(frodo$wis_overprediction_1, frodo$wis_overprediction_2, frodo$wis_overprediction_3, frodo$wis_overprediction_4),
                        wis_overprediction_source = rep(frodo$wis_overprediction_4,4),
                        wis_underprediction = c(frodo$wis_underprediction_1, frodo$wis_underprediction_2, frodo$wis_underprediction_3, frodo$wis_underprediction_4),
                        wis_underprediction_source = rep(frodo$wis_underprediction_4,4),
                        
                        coverage50_source = rep(frodo$coverage50_4,4),
                        coverage50 = c(frodo$coverage50_1, frodo$coverage50_2, frodo$coverage50_3, frodo$coverage50_4),
                        coverage80_source = rep(frodo$coverage80_4,4),
                        coverage80 = c(frodo$coverage80_1, frodo$coverage80_2, frodo$coverage80_3, frodo$coverage80_4),
                        coverage95_source = rep(frodo$coverage95_4,4),
                        coverage95 = c(frodo$coverage95_1, frodo$coverage95_2, frodo$coverage95_3, frodo$coverage95_4),
                        width50_source = rep(frodo$width50_4,4),
                        width50 = c(frodo$width50_1, frodo$width50_2, frodo$width50_3, frodo$width50_4),
                        width80_source = rep(frodo$width80_4,4),
                        width80 = c(frodo$width80_1, frodo$width80_2, frodo$width80_3, frodo$width80_4),
                        width95_source = rep(frodo$width95_4,4),
                        width95 = c(frodo$width95_1, frodo$width95_2, frodo$width95_3, frodo$width95_4),
                        method = 'GBM',
                        variation = c(rep('All Available', n),
                                      rep('Mode of Transmission', n),
                                      rep('Disease', n),
                                      rep('Disease x Source', n)))

n = nrow(lstm)
lstm_long = data.frame(disease_source = rep(lstm$disease_source,4),
                       #disease = rep(lstm$disease,4),
                       wis = c(lstm$wis_1, lstm$wis_2, lstm$wis_3, lstm$wis_4),
                       wis_source = rep(lstm$wis_4,4),

                       wis_dispersion = c(lstm$wis_dispersion_1, lstm$wis_dispersion_2, lstm$wis_dispersion_3, lstm$wis_dispersion_4),
                       wis_dispersion_source = rep(lstm$wis_dispersion_4,4),
                       wis_overprediction = c(lstm$wis_overprediction_1, lstm$wis_overprediction_2, lstm$wis_overprediction_3, lstm$wis_overprediction_4),
                       wis_overprediction_source = rep(lstm$wis_overprediction_4,4),
                       wis_underprediction = c(lstm$wis_underprediction_1, lstm$wis_underprediction_2, lstm$wis_underprediction_3, lstm$wis_underprediction_4),
                       wis_underprediction_source = rep(lstm$wis_underprediction_4,4),

                       coverage50_source = rep(lstm$coverage50_4,4),
                       coverage50 = c(lstm$coverage50_1, lstm$coverage50_2, lstm$coverage50_3, lstm$coverage50_4),
                       coverage80_source = rep(lstm$coverage80_4,4),
                       coverage80 = c(lstm$coverage80_1, lstm$coverage80_2, lstm$coverage80_3, lstm$coverage80_4),
                       coverage95_source = rep(lstm$coverage95_4,4),
                       coverage95 = c(lstm$coverage95_1, lstm$coverage95_2, lstm$coverage95_3, lstm$coverage95_4),
                       width50_source = rep(lstm$width50_4,4),
                       width50 = c(lstm$width50_1, lstm$width50_2, lstm$width50_3, lstm$width50_4),
                       width80_source = rep(lstm$width80_4,4),
                       width80 = c(lstm$width80_1, lstm$width80_2, lstm$width80_3, lstm$width80_4),
                       width95_source = rep(lstm$width95_4,4),
                       width95 = c(lstm$width95_1, lstm$width95_2, lstm$width95_3, lstm$width95_4),
                       method = 'LSTM',
                       variation = c(rep('All Available', n),
                                     rep('Mode of Transmission', n),
                                     rep('Disease', n),
                                     rep('Disease x Source', n)))

n = nrow(smoa)
smoa_long = data.frame(disease_source = rep(smoa$disease_source,4),
                       #disease = rep(smoa$disease,4),
                       wis = c(smoa$wis_1, smoa$wis_2, smoa$wis_3, smoa$wis_4),
                       wis_source = rep(smoa$wis_4,4),
                       
                       wis_dispersion = c(smoa$wis_dispersion_1, smoa$wis_dispersion_2, smoa$wis_dispersion_3, smoa$wis_dispersion_4),
                       wis_dispersion_source = rep(smoa$wis_dispersion_4,4),
                       wis_overprediction = c(smoa$wis_overprediction_1, smoa$wis_overprediction_2, smoa$wis_overprediction_3, smoa$wis_overprediction_4),
                       wis_overprediction_source = rep(smoa$wis_overprediction_4,4),
                       wis_underprediction = c(smoa$wis_underprediction_1, smoa$wis_underprediction_2, smoa$wis_underprediction_3, smoa$wis_underprediction_4),
                       wis_underprediction_source = rep(smoa$wis_underprediction_4,4),
                       
                       coverage50_source = rep(smoa$coverage50_4,4),
                       coverage50 = c(smoa$coverage50_1, smoa$coverage50_2, smoa$coverage50_3, smoa$coverage50_4),
                       coverage80_source = rep(smoa$coverage80_4,4),
                       coverage80 = c(smoa$coverage80_1, smoa$coverage80_2, smoa$coverage80_3, smoa$coverage80_4),
                       coverage95_source = rep(smoa$coverage95_4,4),
                       coverage95 = c(smoa$coverage95_1, smoa$coverage95_2, smoa$coverage95_3, smoa$coverage95_4),
                       width50_source = rep(smoa$width50_4,4),
                       width50 = c(smoa$width50_1, smoa$width50_2, smoa$width50_3, smoa$width50_4),
                       width80_source = rep(smoa$width80_4,4),
                       width80 = c(smoa$width80_1, smoa$width80_2, smoa$width80_3, smoa$width80_4),
                       width95_source = rep(smoa$width95_4,4),
                       width95 = c(smoa$width95_1, smoa$width95_2, smoa$width95_3, smoa$width95_4),
                       method = 'MOA',
                       variation = c(rep('All Available', n),
                                     rep('Mode of Transmission', n),
                                     rep('Disease', n),
                                     rep('Disease x Source', n)))

dat_long = rbind(frodo_long, lstm_long, smoa_long)

dat_long = dat_long[!is.na(dat_long$wis),]
dat_long = dat_long[!grepl('InfluenzaA',dat_long$disease_source),]
dat_long = dat_long[!grepl('InfluenzaB',dat_long$disease_source),]
dat_long = dat_long[dat_long$disease_source != 'Dengue_tycho',]
SPLITS = strsplit(dat_long$disease)
dat_long$disease = unlist(lapply(dat_long$disease_source,  FUN = function(x){paste(strsplit(x, split = '_')[[1]][-length(strsplit(x, split = '_')[[1]])], collapse = '_')}))



library(dplyr)



# Average Coverage Error, measuring prediction interval calibration
dat_long$ace = (1/3)*( abs(dat_long$coverage50-0.5) + abs(dat_long$coverage80-0.8) + abs(dat_long$coverage95-0.95))
dat_long$ace_source = (1/3)*( abs(dat_long$coverage50_source-0.5) + abs(dat_long$coverage80_source-0.8) + abs(dat_long$coverage95_source-0.95))




SUBSET = dat_long
SUBSET$variation_pretty = SUBSET$variation
LABELS = rev(c('Disease x Source', 'Disease','Mode of Transmission', 'All Available'))
SUBSET$variation_pretty = factor(SUBSET$variation_pretty, levels = LABELS)
SUBSET2 = SUBSET[(SUBSET$disease == 'COVID' & SUBSET$variation_pretty != 'Disease') | SUBSET$disease != 'COVID',]
SUBSET2 = SUBSET2[SUBSET2$variation_pretty != 'Disease x Source',]
SUBSET2 = SUBSET2[!(SUBSET2$disease == 'Dengue' & SUBSET2$variation_pretty == 'Disease'),]

p1 = ggplot(SUBSET2)+
  geom_hline(yintercept = 1, color = 'gray', linetype = 2)+
  geom_point(aes(y=wis/wis_source, x = variation_pretty, fill = variation_pretty), shape = 21, color = 'black',alpha = 0.3, position = position_jitterdodge())+
  geom_boxplot(aes(y=wis/wis_source, x = variation_pretty, fill = variation_pretty), alpha = 1, outlier.shape = NA)+
  xlab('')+
  ylab('WIS, Relative to Single Data Stream')+
  labs(title = 'WIS, Relative to Disease x Source')+
  facet_grid(ifelse(disease == 'COVID', 'COVID-19',ifelse(disease == 'Dengue','Dengue','Other') )~method, scales=  'free')+
  scale_color_manual(values = as.vector(palette.colors(4, palette = 'Tableau 10')),
                     breaks = LABELS)+
  scale_fill_manual(values = as.vector(palette.colors(4, palette = 'Tableau 10')),
                    breaks = LABELS)+
  theme_classic()+
  theme(legend.position = 'top')+
  coord_flip(ylim=c(0.65,1.3))+
  guides(color = 'none', fill = 'none')+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
ggsave(p1,filename = paste0('~/GitLab/frodo/summaries_uq/',"wis_source_noaggBOTH.pdf"), width = 10, height = 4)





SUBSET = dat_long
SUBSET$variation_pretty = SUBSET$variation
LABELS = rev(c('Disease x Source', 'Disease','Mode of Transmission', 'All Available'))
SUBSET$variation_pretty = factor(SUBSET$variation_pretty, levels = LABELS)
SUBSET2 = SUBSET[(SUBSET$disease == 'COVID' & SUBSET$variation_pretty != 'Disease') | SUBSET$disease != 'COVID',]
SUBSET2 = SUBSET2[SUBSET2$variation_pretty != 'Disease x Source',]
SUBSET2 = SUBSET2[!(SUBSET2$disease == 'Dengue' & SUBSET2$variation_pretty == 'Disease'),]

SUBSET2$variation_pretty = as.character(SUBSET2$variation_pretty)
SUBSET2$variation_pretty[SUBSET2$variation_pretty == 'Disease x Source'] = 'Single Data Stream'
LABELS = rev(c('Single Data Stream', 'Disease','Mode of Transmission', 'All Available'))
SUBSET2$variation_pretty = factor(SUBSET2$variation_pretty, levels = LABELS)

p1 = ggplot(SUBSET2)+
  geom_hline(yintercept = 0.95, color = 'gray', linetype = 2)+
  geom_point(aes(y=coverage95, x = variation_pretty, fill = variation_pretty), shape = 21, color = 'black',alpha = 0.3, position = position_jitterdodge())+
  geom_boxplot(aes(y=coverage95, x = variation_pretty, fill = variation_pretty), alpha = 1, outlier.shape = NA)+
  geom_point(aes(y=coverage95_source, x = 'Single Data Stream', fill = 'Single Data Stream'), shape = 21, color = 'black',alpha = 0.3, position = position_jitterdodge())+
  geom_boxplot(aes(y=coverage95_source, x = 'Single Data Stream', fill = 'Single Data Stream'), alpha = 1, outlier.shape = NA)+
  xlab('')+
  ylab('95% Prediction Interval Coverage')+
  labs(title = '95% Prediction Interval Coverage')+
  facet_grid(ifelse(disease == 'COVID', 'COVID-19',ifelse(disease == 'Dengue','Dengue','Other') )~method, scales=  'free')+
  scale_color_manual(values = as.vector(palette.colors(4, palette = 'Tableau 10')),
                     breaks = LABELS)+
  scale_fill_manual(values = as.vector(palette.colors(4, palette = 'Tableau 10')),
                    breaks = LABELS)+
  theme_classic()+
  theme(legend.position = 'top')+
  coord_flip(ylim=c(0.6,1))+
  guides(color = 'none', fill = 'none')+
  theme(panel.border = element_rect(color = "black", fill = NA, linewidth = 1))
ggsave(p1,filename = paste0('~/GitLab/frodo/summaries_uq/',"coverage_source_noagg.pdf"), width = 10, height = 5)


