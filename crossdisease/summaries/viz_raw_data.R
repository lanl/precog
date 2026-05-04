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

####################
### Read in Data ###
####################

FILES = list.files('~/GitLab/frodo/raw_data/')
FILES = FILES[!grepl('foodnet', FILES)] #monthly only
FILES = FILES[!grepl('Polio', FILES)] #duplicated 

LENS = rep(NA, length(FILES))
for(i in 1:length(FILES)){
  dat = readRDS(paste0('~/GitLab/frodo/raw_data/', FILES[i]))
  CADENCE = unlist(lapply(dat, FUN = function(x){x$ts_time_cadence}))
  LENS[i] = sum(CADENCE == 'weekly')
}

FILES = FILES[LENS>0]
LENS = LENS[LENS>0]

SPLITS = strsplit(FILES, split = '_')

SOURCES = apply(cbind(1:length(SPLITS)),1,FUN = function(x,SPLITS){SPLITS[[x]][length(SPLITS[[x]])]} , SPLITS)
SOURCES = gsub('.RDS', '', SOURCES)
SOURCES[SOURCES == 'aido'] = 'AIDO'
SOURCES[SOURCES == 'jhuowid'] = 'JHU CSSEGIS/OWID'
SOURCES[SOURCES == 'opendengue'] = 'OpenDengue/NOAA'
SOURCES[SOURCES == 'deSouza'] = 'de Souza et al. (2023)'
SOURCES[SOURCES == 'ushhs'] = 'US HHS'
SOURCES[SOURCES == 'usflunet'] = 'US FluNet/WHO NREVSS'
SOURCES[SOURCES == 'whoflunet'] = 'WHO FluNet'
SOURCES[SOURCES == 'nndss'] = 'US NNDSS'
SOURCES[SOURCES == 'tycho'] = 'Project Tycho'
SOURCES[SOURCES == 'who'] = 'WHO'

DISEASES_ALL = apply(cbind(1:length(SPLITS)),1,FUN = function(x,SPLITS){paste0(SPLITS[[x]][1:(length(SPLITS[[x]])-1)], collapse = ' ')} , SPLITS)
DISEASES_ALL[DISEASES_ALL == 'COVID'] = 'COVID-19'
DISEASES_ALL[DISEASES_ALL == 'Influenza1918'] = 'Influenza'

length(DISEASES_ALL) #101, Number of Data Streams
length(unique(DISEASES_ALL)) #76, number of unique disease variations (influenza, influenza A, influenza B, etc.)


### Ignoring Subtypes ### 
DISEASES = DISEASES_ALL
DISEASES[DISEASES == 'InfluenzaA'] = 'Influenza'
DISEASES[DISEASES == 'InfluenzaB'] = 'Influenza'
DISEASES[DISEASES == 'DengueSero1'] = 'Dengue'
DISEASES[DISEASES == 'DengueSero2'] = 'Dengue'
DISEASES[DISEASES == 'DengueSero3'] = 'Dengue'
DISEASES[DISEASES == 'DengueSero4'] = 'Dengue'
DISEASES[DISEASES == 'Ecoli Infection Shiga'] = 'Ecoli Infection'
DISEASES[DISEASES == 'Hepatitis A'] = 'Hepatitis'
DISEASES[DISEASES == 'Hepatitis B'] = 'Hepatitis'
DISEASES[DISEASES == 'Hepatitis NonAB'] = 'Hepatitis'

length(unique(DISEASES)) #66, number of "diseases"


###################
### Make Sankey ###
###################

library(networkD3)
library("curl")


DATA = data.frame(Source = SOURCES, Disease = DISEASES, Lengths = LENS)

nodes = data.frame("name" = c(unique(SOURCES),unique(DISEASES)))
nodes$short_name = c(unique(SOURCES),
                     rep('stage3',length(unique(DISEASES))))
### Plotting parameter, defines how far away the nodes are from one another vertically
BASELINE_SANKEY_HEIGHT = 500
### Define links to time scales
links = data.frame(source = NULL, target = NULL, value = NULL, link_group = NULL)
### Define links to variables
for(i in c(1:length(DATA[,1]))){
  links = rbind(links, data.frame(source =which(nodes$name == DATA$Source[i])-1, target = which(nodes$name == DATA$Disease[i])-1, value = DATA[i,'Lengths'], link_group = DATA$Source[i]))
}

### Note: this figure is toggle-able
p1 = sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "source", Target = "target",
                   Value = "value", NodeID = "name",
                   nodeWidth = 20, NodeGroup = 'short_name',
                   LinkGroup = 'link_group',
                   fontSize= 12,  margin = list(left = 260))
p2 = htmlwidgets::onRender(p1, '
  function(el) {
    var cols_x = this.sankey.nodes().map(d => d.x).filter((v, i, a) => a.indexOf(v) === i).sort(function(a, b){return a - b});
    var labels = [ " ", " "];
    cols_x.forEach((d, i) => {
      d3.select(el).select("svg")
        .append("text")
        .attr("x", d)
        .attr("y", 12)
        .attr("font-family", "Arial")
        .attr("font-size", 30)
        .text(labels[i]);
      d3.select(el).selectAll(".link")
        .style("stroke-opacity", 0.3); 
    })
  }
')
p2[["sizingPolicy"]][["defaultWidth"]] <- "90%" # Example using percentage
p2[["sizingPolicy"]][["defaultHeight"]] <- "400px" # Example using pixels

require(htmlwidgets)
saveWidget(p2, file=paste0('~/GitLab/frodo/summaries/',"sankey.html"), knitrOptions = list(width = 600, height = 700))


#######################
### Other Summaries ###
#######################



dat = data.frame(FILES = FILES, SOURCES = SOURCES, DISEASES = DISEASES)
dat$length_of_list = NA
dat$total_obs_in_list = NA
dat$total_obs_in_list_nonzero = NA
dat$median_obs_per_ts = NA

for(i in 1:length(FILES)){
  data_list = readRDS(paste0('~/GitLab/frodo/raw_data/', FILES[i]))
  CADENCE = unlist(lapply(data_list, FUN = function(x){x$ts_time_cadence}))
  data_list = data_list[which(CADENCE == 'weekly')]
  dat$length_of_list[i] = length(data_list) 
  nums = unlist(lapply(data_list, FUN = function(x){length(x$ts)}))
  dat$total_obs_in_list[i] = sum(nums)
  dat$median_obs_per_ts[i] = median(nums)
  nums = unlist(lapply(data_list, FUN = function(x){length(x$ts[x$ts>0])}))
  dat$total_obs_in_list_nonzero[i] = sum(nums)
}

sum(dat$length_of_list) #6875, number of time series
sum(dat$total_obs_in_list) #5605772, total observations in list
median(dat$median_obs_per_ts) #458, median number of observations per time series

sum(dat$total_obs_in_list_nonzero) #2269070, number of nonzero observations



AGG = aggregate(total_obs_in_list_nonzero~DISEASES, data = dat, FUN = sum)
AGG = AGG[order(AGG$total_obs_in_list_nonzero),]
AGG$DISEASES = factor(AGG$DISEASES, levels = rev(AGG$DISEASES))
AGG$fill_cat = round(log10(AGG$total_obs_in_list_nonzero),0)
p1 = ggplot(AGG)+
  geom_bar(aes(y=DISEASES, x= log10(total_obs_in_list_nonzero), fill = factor(fill_cat)), linewidth = 0.1, color = 'black',stat = 'identity')+
  annotate('label', label = rev(c('~100', '~1K', '~10K', '~100K')), x = rep(1,4), y = c(11,36,57,64))+
  theme_classic()+
  scale_x_continuous(expand = c(0,0))+
  xlab('Nonzero Values (log10)')+
  ylab('')+
  scale_fill_brewer('Order', palette = 'Spectral', breaks = c(2:6), labels = c('100','1K','10K','100K','1M'))+
  guides(fill = 'none', color = 'none')
ggsave(p1, filename = paste0('~/GitLab/frodo/summaries/',"histogram.png"), width = 8, height = 10)





#####################################
### Visualizing Disease Timelines ###
#####################################


dat = data.frame(FILES = FILES, SOURCES = SOURCES, DISEASES = DISEASES, DISEASES_ALL = DISEASES_ALL)

DATES = replicate(length(FILES), list(NULL))
for(i in 1:length(FILES)){
  data_list = readRDS(paste0('~/GitLab/frodo/raw_data/', FILES[i]))
  CADENCE = unlist(lapply(data_list, FUN = function(x){x$ts_time_cadence}))
  data_list = data_list[which(CADENCE == 'weekly')]
  DATES_SINGLE = sort(unique(unlist(lapply(data_list, FUN = function(x){as.character(x$ts_dates)}))))
  
  TIME_GAPS = as.numeric(difftime(as.Date(DATES_SINGLE[-1], format = '%Y-%m-%d'),as.Date(DATES_SINGLE[-length(DATES_SINGLE)], format = '%Y-%m-%d'), units = 'days'))
  TIME_GAPS[TIME_GAPS<7] = 7
  if(length(unique(TIME_GAPS[TIME_GAPS])) == 1){
    DATES[[i]] = data.frame(dat[i,c('FILES', 'SOURCES', 'DISEASES', 'DISEASES_ALL')],dates_lower = min(DATES_SINGLE), dates_upper = max(DATES_SINGLE))
  }else{
    OFF = which(TIME_GAPS > 7)
    temp = NULL
    for(u in 1:length(OFF)){
      if(u == 1){
        temp = rbind(temp, data.frame(dat[i,c('FILES', 'SOURCES', 'DISEASES', 'DISEASES_ALL')],dates_lower = min(DATES_SINGLE[1:OFF[u]]), dates_upper = DATES_SINGLE[OFF[u]]))
      }
      if(u > 1 & u != length(OFF)){
        temp = rbind(temp, data.frame(dat[i,c('FILES', 'SOURCES', 'DISEASES', 'DISEASES_ALL')],dates_lower = DATES_SINGLE[OFF[u]+1], dates_upper = DATES_SINGLE[OFF[u+1]]))
      }
      if(u == length(OFF)){
        temp = rbind(temp, data.frame(dat[i,c('FILES', 'SOURCES', 'DISEASES', 'DISEASES_ALL')],dates_lower = DATES_SINGLE[OFF[u]+1], dates_upper = max(DATES_SINGLE)))
      }
    }
    DATES[[i]] = temp
    print(FILES[i])
  }
}
DATES_LONG = do.call('rbind', DATES)

DATES_LONG = DATES_LONG[DATES_LONG$DISEASES == DATES_LONG$DISEASES_ALL,]

nrow(DATES_LONG) #91, diseases x sources (so not counting Influenza A and Influenza B differently)

DATES_LONG$dates_lower = as.Date(DATES_LONG$dates_lower, format = '%Y-%m-%d')
DATES_LONG$dates_upper = as.Date(DATES_LONG$dates_upper, format = '%Y-%m-%d')
DATES_LONG = DATES_LONG[order(DATES_LONG$dates_lower),]
DATES_LONG$DISEASES = factor(DATES_LONG$DISEASES, levels = unique(DATES_LONG$DISEASES))
DATES_LONG$KEY = paste0(DATES_LONG$DISEASES, '_', DATES_LONG$SOURCES, '_', DATES_LONG$dates_lower)
DATES_LONG$KEY = factor(DATES_LONG$KEY)


DATES_LONG$in_eval = 0
DATES_LONG$in_eval[DATES_LONG$SOURCES == 'US NNDSS' & DATES_LONG$DISEASES == 'Anaplasmosis'] = 1
DATES_LONG$in_eval[DATES_LONG$SOURCES == 'US NNDSS' & DATES_LONG$DISEASES == 'Babesiosis'] = 1
DATES_LONG$in_eval[DATES_LONG$SOURCES == 'US NNDSS' & DATES_LONG$DISEASES == 'Campylobacteriosis'] = 1
DATES_LONG$in_eval[DATES_LONG$SOURCES == 'US NNDSS' & DATES_LONG$DISEASES == 'Chlamydia'] = 1
DATES_LONG$in_eval[DATES_LONG$SOURCES == 'JHU CSSEGIS/OWID' & DATES_LONG$DISEASES == 'COVID-19'] = 1
DATES_LONG$in_eval[DATES_LONG$SOURCES == 'US NNDSS' & DATES_LONG$DISEASES == 'Cryptosporidiosis'] = 1
DATES_LONG$in_eval[DATES_LONG$SOURCES == 'OpenDengue/NOAA' & DATES_LONG$DISEASES == 'Dengue'] = 1
DATES_LONG$in_eval[DATES_LONG$SOURCES == 'US NNDSS' & DATES_LONG$DISEASES == 'Ecoli Infection'] = 1
DATES_LONG$in_eval[DATES_LONG$SOURCES == 'US NNDSS' & DATES_LONG$DISEASES == 'Ehrlichiosis'] = 1
DATES_LONG$in_eval[DATES_LONG$SOURCES == 'US NNDSS' & DATES_LONG$DISEASES == 'Giardiasis'] = 1
DATES_LONG$in_eval[DATES_LONG$SOURCES == 'US NNDSS' & DATES_LONG$DISEASES == 'Gonorrhea'] = 1
DATES_LONG$in_eval[DATES_LONG$SOURCES == 'US NNDSS' & DATES_LONG$DISEASES == 'Haemophilus Influenzae'] = 1
DATES_LONG$in_eval[DATES_LONG$SOURCES == 'US FluNet/WHO NREVSS' & DATES_LONG$DISEASES == 'Influenza'] = 1
DATES_LONG$in_eval[DATES_LONG$SOURCES == 'US HHS' & DATES_LONG$DISEASES == 'Influenza'] = 1
DATES_LONG$in_eval[DATES_LONG$SOURCES == 'WHO FluNet' & DATES_LONG$DISEASES == 'Influenza'] = 1
DATES_LONG$in_eval[DATES_LONG$SOURCES == 'US NNDSS' & DATES_LONG$DISEASES == 'Legionellosis'] = 1
DATES_LONG$in_eval[DATES_LONG$SOURCES == 'US NNDSS' & DATES_LONG$DISEASES == 'Pertussis'] = 1
DATES_LONG$in_eval[DATES_LONG$SOURCES == 'US NNDSS' & DATES_LONG$DISEASES == 'Salmonellosis'] = 1
DATES_LONG$in_eval[DATES_LONG$SOURCES == 'US NNDSS' & DATES_LONG$DISEASES == 'Shigellosis'] = 1
DATES_LONG$in_eval[DATES_LONG$SOURCES == 'US NNDSS' & DATES_LONG$DISEASES == 'Tuberculosis'] = 1

p1 = ggplot(DATES_LONG)+
  geom_hline(yintercept = as.Date(c('2010-01-01'), format = '%Y-%m-%d'), linetype = 2, color = 'gray')+
  geom_vline(xintercept = seq(0.5,74.5,1), linetype = 1, color = 'lightgray', linewidth = 0.5)+
  #geom_rect(aes(ymin=dates_lower-1, ymax = dates_upper+1, xmin = as.numeric(DISEASES)-0.5, xmax = as.numeric(DISEASES)+0.5), color = 'black',fill = NA)+
  geom_rect(aes(ymin=dates_lower, ymax = dates_upper, xmin = as.numeric(DISEASES)-0.4, xmax = as.numeric(DISEASES)+0.4, fill = SOURCES), color = 'black',alpha = 0.7)+
  geom_segment(aes(y=dates_lower, yend = dates_upper, x = DISEASES,  group = KEY),alpha = 0)+
  scale_y_date()+
  xlab('')+
  ylab('Dates')+
  theme(legend.position = 'top')+
  coord_flip()+
  #scale_color_manual('In Evaluation', breaks = factor(c(0,1)), labels = c('No','Yes'), values = c('gray','black'))+
  scale_fill_brewer('Source', palette = 'Spectral')+
  scale_y_date(expand = c(0,0))+
  guides(fill = guide_legend(nrow = 3))
ggsave(filename = paste0('~/GitLab/frodo/summaries/',"date_intervals.pdf"), width = 10, height = 10)



DISEASES_SUB = c('Shigellosis', 'Salmonellosis', 'Pertussis',  'Legionellosis', 
                 'Influenza', 'Dengue', 'COVID-19', 'Haemophilus Influenzae',
                 'Gonorrhea', 'Giardiasis', 'Ehrlichiosis', 'Cryptosporidiosis',
                 'Chlamydia', 'Campylobacteriosis',
                 'Babesiosis', 'Anaplasmosis', 'Tuberculosis', 'Ecoli Infection')

SUBSET = DATES_LONG[DATES_LONG$DISEASES %in% DISEASES_SUB,]
SUBSET = SUBSET[order(SUBSET$dates_lower),]
SUBSET$DISEASES = factor(SUBSET$DISEASES, levels = unique(SUBSET$DISEASES))
SUBSET$KEY = paste0(SUBSET$DISEASES, '_', SUBSET$SOURCES, '_', SUBSET$dates_lower)
SUBSET$KEY = factor(SUBSET$KEY)
p1 = ggplot(SUBSET)+
  geom_hline(yintercept = as.Date(c('2010-01-01'), format = '%Y-%m-%d'), linetype = 2, color = 'gray')+
  geom_vline(xintercept = seq(0.5,22.5,1), linetype = 1, color = 'lightgray', linewidth = 0.5)+
  geom_rect(aes(ymin=dates_lower, ymax = dates_upper, xmin = as.numeric(DISEASES)-0.4, xmax = as.numeric(DISEASES)+0.4, fill = SOURCES), color = NA, alpha = 0.7)+
  geom_rect(aes(ymin=dates_lower, ymax = dates_upper, xmin = as.numeric(DISEASES)-0.4, xmax = as.numeric(DISEASES)+0.4), fill = NA, color = 'black', alpha = 0.7, data = SUBSET[SUBSET$in_eval == 1,])+
  geom_segment(aes(y=dates_lower, yend = dates_upper, x = DISEASES,  group = KEY),alpha = 0)+
  scale_y_date()+
  xlab('')+
  ylab('Dates')+
  theme(legend.position = 'top')+
  coord_flip()+
  #scale_color_manual('In Evaluation', breaks = factor(c(0,1)), labels = c('No','Yes'), values = c('gray','black'))+
  scale_fill_brewer('Source', palette = 'Spectral')+
  scale_y_date(expand = c(0,0))+
  guides(fill = guide_legend(nrow = 3), alpha = 'none')
ggsave(p1,filename = paste0('~/GitLab/frodo/summaries/',"date_intervals_subset.pdf"), width = 7, height = 5)


