###########################################
###########################################
### Visualize and Summarize Time Series ###
###########################################
###########################################


library(ggplot2)
library(data.table)
library(GGally)
library(viridis)
library(ggrepel)
theme_set(theme_classic())

savepath <- "./_Organized_Lists/"

FILES = list.files(savepath)
SOURCES = gsub('.RDS','',unlist(lapply(strsplit(FILES, split = '_'), FUN = function(x){x[length(x)]})))
DISEASES = unlist(lapply(strsplit(FILES, split = '_'), FUN = function(x){paste(x[-length(x)], collapse = '_')}))


LENS = rep(NA, length(FILES))
for(i in 1:length(FILES)){
  dat = readRDS(paste0(savepath, FILES[i]))
  LENS[i] = length(dat) 
}


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
SOURCES[SOURCES == 'foodnet'] = 'US FoodNet'

DISEASES = apply(cbind(1:length(SPLITS)),1,FUN = function(x,SPLITS){paste0(SPLITS[[x]][1:(length(SPLITS[[x]])-1)], collapse = ' ')} , SPLITS)
DISEASES[DISEASES == 'COVID'] = 'COVID-19'




########################
### Plot Time Series ###
########################

INDS = 1:length(FILES)
INDS = INDS[!grepl('COVID',FILES)]


for(i in INDS){
  data_list = readRDS(paste0(savepath, FILES[i]))

  #lags = NULL
  for(j in 1:length(data_list)){
    A = table(diff(as.Date(data_list[[j]]$ts_dates_actual, format = '%Y-%m-%d')))
    if(length(A)==0){
      A = table(diff(as.Date(data_list[[j]]$ts_dates, format = '%Y-%m-%d')))
    }
    #lags = rbind(lags,data.frame(diff = names(A), freq = as.numeric(A)))
    if(data_list[[j]]$ts_time_cadence %in% c('daily', 'weekly')){
      DATES = as.Date(data_list[[j]]$ts_dates, format = '%Y-%m-%d')
    }else{
      DATES = as.Date(paste0(gsub('_','-',data_list[[j]]$ts_dates),'-01'), format = '%Y-%m-%d')
    }
    if(j == 1){
      p1 = ggplot(data.frame(start_date = DATES, ts = data_list[[j]]$ts, loc_id = j))+
        geom_line(aes(x=start_date, y = ts, group = loc_id, color = loc_id))+
        geom_point(aes(x=start_date, y = ts, group = loc_id, color = loc_id))
    }else{
      p1 = p1 + geom_line(aes(x=start_date, y = ts, group = loc_id, color = loc_id), data = data.frame(start_date = DATES, ts = data_list[[j]]$ts, loc_id = j))+
        geom_point(aes(x=start_date, y = ts, group = loc_id, color = loc_id), data = data.frame(start_date = DATES, ts = data_list[[j]]$ts, loc_id = j))
    }
    p2 = p1 + guides(color = 'none')
  }
  #AGG = aggregate(freq~diff, FUN = sum, data = lags)
  ggsave(p2, filename = paste0(savepath,'../_Organized_Figs/Time_Series_Plots/',gsub('.RDS','',FILES[i]),'.png'), width = 10, height = 6)
  gc()
}



for(i in INDS){
  data_list = readRDS(paste0(savepath, FILES[i]))
  
  for(j in 1:length(data_list)){
    if(data_list[[j]]$ts_time_cadence %in% c('daily', 'weekly')){
      DATES = as.Date(data_list[[j]]$ts_dates, format = '%Y-%m-%d')
    }else{
      DATES = as.Date(paste0(gsub('_','-',data_list[[j]]$ts_dates),'-01'), format = '%Y-%m-%d')
    }
    if(j == 1){
      p1 = ggplot(data.frame(start_date = DATES, ts = data_list[[j]]$ts, loc_id = j))+
        geom_line(aes(x=start_date, y = ts/max(ts), group = loc_id, color = loc_id))+
        geom_point(aes(x=start_date, y = ts/max(ts), group = loc_id, color = loc_id))
    }else{
      p1 = p1 + geom_line(aes(x=start_date, y = ts/max(ts), group = loc_id, color = loc_id), data = data.frame(start_date = DATES, ts = data_list[[j]]$ts, loc_id = j))+
        geom_point(aes(x=start_date, y = ts/max(ts), group = loc_id, color = loc_id), data = data.frame(start_date = DATES, ts = data_list[[j]]$ts, loc_id = j))
    }
    p2 = p1 + guides(color = 'none')
  }
  ggsave(p2, filename = paste0(savepath,'../_Organized_Figs/Time_Series_Plots/',gsub('.RDS','',FILES[i]),'_scaled.png'), width = 10, height = 6)
  gc()
}


####################################
### Obtain Time Series Summaries ###
####################################

dat = data.frame(FILES = FILES, SOURCES = SOURCES, DISEASES = DISEASES)
dat$length_of_list = NA
dat$total_obs_in_list = NA
dat$total_obs_in_list_nonzero = NA
dat$median_obs_per_ts = NA

OBS_VEC = c()
for(i in 1:length(FILES)){
  data_list = readRDS(paste0(savepath, FILES[i]))
  CADENCE = unlist(lapply(data_list, FUN = function(x){x$ts_time_cadence}))
  if(length(unique(CADENCE))>1){
    data_list = data_list[which(CADENCE != 'daily')] #every daily time series has a weekly analog
  }
  dat$length_of_list[i] = length(data_list) 
  nums = unlist(lapply(data_list, FUN = function(x){length(x$ts)}))
  OBS_VEC = c(OBS_VEC, nums)
  dat$total_obs_in_list[i] = sum(nums)
  dat$median_obs_per_ts[i] = median(nums)
  nums = unlist(lapply(data_list, FUN = function(x){length(x$ts[x$ts>0])}))
  dat$total_obs_in_list_nonzero[i] = sum(nums)
}


sum(dat$length_of_list)
sum(dat$total_obs_in_list)
median(dat$median_obs_per_ts)

sum(dat$total_obs_in_list_nonzero)
length(unique(DISEASES))

AGG = aggregate(total_obs_in_list_nonzero~DISEASES, data = dat, FUN = sum)
AGG = AGG[order(AGG$total_obs_in_list_nonzero),]
AGG$DISEASES = factor(AGG$DISEASES, levels = rev(AGG$DISEASES))
AGG$fill_cat = round(log10(AGG$total_obs_in_list_nonzero),0)
p1 = ggplot(AGG)+
  geom_bar(aes(y=DISEASES, x= log10(total_obs_in_list_nonzero), fill = factor(fill_cat)), linewidth = 0.1, color = 'black',stat = 'identity')+
  annotate('label', label = rev(c('~100', '~1K', '~10K', '~100K')), x = rep(1,4), y = c(13,42,67,78))+
  theme_classic()+
  scale_x_continuous(expand = c(0,0))+
  xlab('Nonzero Values (log10)')+
  ylab('')+
  scale_fill_brewer('Order', palette = 'Spectral', breaks = c(2:6), labels = c('100','1K','10K','100K','1M'))+
  guides(fill = 'none', color = 'none')
ggsave(p1, filename = paste0('./_Organized_Figs/',"histogram_diseases.png"), width = 8, height = 10)




install.packages("wordcloud")
library(wordcloud)

AGG = aggregate(total_obs_in_list_nonzero~DISEASES, data = dat, FUN = sum)
AGG = AGG[order(AGG$total_obs_in_list_nonzero),]
AGG$DISEASES = factor(AGG$DISEASES, levels = rev(AGG$DISEASES))
AGG$fill_cat = round(log10(AGG$total_obs_in_list_nonzero),0)
AGG = AGG[order(AGG$total_obs_in_list_nonzero),]
library(RColorBrewer)

COLORS = rev(viridis(10))[-1]

pdf(paste0('./_Organized_Figs/',"wordcloud.pdf"), width = 10, height = 10)
wordcloud(
  words = AGG$DISEASES,
  freq = AGG$total_obs_in_list_nonzero/100000,
  scale = c(2,0.5),      # largest and smallest size
  colors = COLORS
)
dev.off()


p1 = ggplot(data.frame(value = OBS_VEC))+
  geom_histogram(aes(x=value), bins = 100, fill = 'lightblue', color = 'black', lwd = 0.1)+
  #geom_histogram(aes(x=log10(value)), bins = 100, fill = 'lightblue', color = 'black', lwd = 0.1)+
  #geom_histogram(aes(x=log10(value)), bins = 100)+
  xlab('Time Series Lengths')+
  ylab('Number of Time Series')#+
  #scale_y_continuous(expand = c(0,0))+
  #scale_x_continuous(expand = c(0,0))
ggsave(p1, filename = paste0('./_Organized_Figs/',"ts_lengths.png"), width = 4, height = 3)

###################
### Make Sankey ###
###################

library(networkD3)
library("curl")

DISEASES_EXCLUDE = c('InfluenzaA', 'InfluenzaB', 'DengueSero1',
                     'DengueSero2', 'DengueSero3', 'DengueSero4',
                     'Hepatitis A','Hepatitis B','Hepatitis NonAB')


DATA = data.frame(Source = SOURCES, Disease = DISEASES, Lengths = LENS)
DATA$Disease[DATA$Disease == 'Influenza1918'] = 'Influenza'
DATA = DATA[!(DATA$Disease %in% DISEASES_EXCLUDE),]

nodes = data.frame("name" = c(unique(DATA$Source),unique(DATA$Disease)))
nodes$short_name = c(unique(DATA$Source),
                     rep('stage3',length(unique(DATA$Disease))))
### Plotting parameter, defines how far away the nodes are from one another vertically
BASELINE_SANKEY_HEIGHT = 500
### Define links to time scales
links = data.frame(source = NULL, target = NULL, value = NULL, link_group = NULL)
### Define links to variables
for(i in c(1:length(DATA[,1]))){
  links = rbind(links, data.frame(source =which(nodes$name == DATA$Source[i])-1, target = which(nodes$name == DATA$Disease[i])-1, value = DATA[i,'Lengths'], link_group = DATA$Source[i]))
}


### Note: this figure is toggle-able
p1 = networkD3::sankeyNetwork(Links = links, Nodes = nodes,
                   Source = "source", Target = "target",
                   Value = "value", NodeID = "name",
                   nodeWidth = 20, NodeGroup = 'short_name',
                   LinkGroup = 'link_group',
                   fontSize= 12,  margin = list(left = 260))#,

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
saveWidget(p2, file=paste0('./_Organized_Figs/',"sankey.html"), knitrOptions = list(width = 600, height = 700))


##########################################
### Visualizing Disease Time Coverages ###
##########################################


dat = data.frame(FILES = FILES, SOURCES = SOURCES, DISEASES = DISEASES)

DATES = replicate(length(FILES), list(NULL))
for(i in 1:length(FILES)){
  data_list = readRDS(paste0(savepath, FILES[i]))
  CADENCE = unlist(lapply(data_list, FUN = function(x){x$ts_time_cadence}))
  data_list = data_list[which(CADENCE == 'weekly')]
  if(length(data_list)>0){
    DATES_SINGLE = sort(unique(unlist(lapply(data_list, FUN = function(x){as.character(x$ts_dates)}))))
    
    TIME_GAPS = as.numeric(difftime(as.Date(DATES_SINGLE[-1], format = '%Y-%m-%d'),as.Date(DATES_SINGLE[-length(DATES_SINGLE)], format = '%Y-%m-%d'), units = 'days'))
    TIME_GAPS[TIME_GAPS<7] = 7
    if(length(unique(TIME_GAPS[TIME_GAPS])) == 1){
      DATES[[i]] = data.frame(dat[i,c('FILES', 'SOURCES', 'DISEASES')],dates_lower = min(DATES_SINGLE), dates_upper = max(DATES_SINGLE))
    }else{
      OFF = which(TIME_GAPS > 7)
      temp = NULL
      for(u in 1:length(OFF)){
        if(u == 1){
          temp = rbind(temp, data.frame(dat[i,c('FILES', 'SOURCES', 'DISEASES')],dates_lower = min(DATES_SINGLE[1:OFF[u]]), dates_upper = DATES_SINGLE[OFF[u]]))
        }
        if(u > 1 & u != length(OFF)){
          temp = rbind(temp, data.frame(dat[i,c('FILES', 'SOURCES', 'DISEASES')],dates_lower = DATES_SINGLE[OFF[u]+1], dates_upper = DATES_SINGLE[OFF[u+1]]))
        }
        if(u == length(OFF)){
          temp = rbind(temp, data.frame(dat[i,c('FILES', 'SOURCES', 'DISEASES')],dates_lower = DATES_SINGLE[OFF[u]+1], dates_upper = max(DATES_SINGLE)))
        }
      }
      DATES[[i]] = temp
      print(FILES[i]) #gaps may occur because of different locations having different time intervals
    }
  }
}
DATES_LONG = do.call('rbind', DATES)



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
DATES_LONG$in_eval[DATES_LONG$SOURCES == 'US NNDSS' & DATES_LONG$DISEASES == 'Ehrlichiosis'] = 1
DATES_LONG$in_eval[DATES_LONG$SOURCES == 'US NNDSS' & DATES_LONG$DISEASES == 'Giardiasis'] = 1
DATES_LONG$in_eval[DATES_LONG$SOURCES == 'US NNDSS' & DATES_LONG$DISEASES == 'Gonorrhea'] = 1
DATES_LONG$in_eval[DATES_LONG$SOURCES == 'US NNDSS' & DATES_LONG$DISEASES == 'Haemophilus Influenzae'] = 1
DATES_LONG$in_eval[DATES_LONG$SOURCES == 'US-FluNet/WHO NREVSS' & DATES_LONG$DISEASES == 'Influenza'] = 1
DATES_LONG$in_eval[DATES_LONG$SOURCES == 'USHHS' & DATES_LONG$DISEASES == 'Influenza'] = 1
DATES_LONG$in_eval[DATES_LONG$SOURCES == 'WHO-FluNet' & DATES_LONG$DISEASES == 'Influenza'] = 1
DATES_LONG$in_eval[DATES_LONG$SOURCES == 'US NNDSS' & DATES_LONG$DISEASES == 'Legionellosis'] = 1
DATES_LONG$in_eval[DATES_LONG$SOURCES == 'US NNDSS' & DATES_LONG$DISEASES == 'Malaria'] = 1
DATES_LONG$in_eval[DATES_LONG$SOURCES == 'US NNDSS' & DATES_LONG$DISEASES == 'Pertussis'] = 1
DATES_LONG$in_eval[DATES_LONG$SOURCES == 'US NNDSS' & DATES_LONG$DISEASES == 'Salmonellosis'] = 1
DATES_LONG$in_eval[DATES_LONG$SOURCES == 'US NNDSS' & DATES_LONG$DISEASES == 'Shigellosis'] = 1

DATES_TO_PLOT = c(as.Date(c('1900-01-01'), format = '%Y-%m-%d'),
                  as.Date(c('1910-01-01'), format = '%Y-%m-%d'),
                  as.Date(c('1920-01-01'), format = '%Y-%m-%d'),
                  as.Date(c('1930-01-01'), format = '%Y-%m-%d'),
                  as.Date(c('1940-01-01'), format = '%Y-%m-%d'),
                  as.Date(c('1950-01-01'), format = '%Y-%m-%d'),
                  as.Date(c('1960-01-01'), format = '%Y-%m-%d'),
                  as.Date(c('1970-01-01'), format = '%Y-%m-%d'),
                  as.Date(c('1980-01-01'), format = '%Y-%m-%d'),
                  as.Date(c('1990-01-01'), format = '%Y-%m-%d'),
                  as.Date(c('2000-01-01'), format = '%Y-%m-%d'),
                  as.Date(c('2010-01-01'), format = '%Y-%m-%d'),
                  as.Date(c('2020-01-01'), format = '%Y-%m-%d'))
p1 = ggplot(DATES_LONG)+
  geom_hline(yintercept = DATES_TO_PLOT, linetype = 2, color = 'gray')+
  geom_vline(xintercept = seq(0.5,78.5,1), linetype = 1, color = 'lightgray', linewidth = 0.5)+
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
  scale_y_date(expand = c(0,0), breaks = DATES_TO_PLOT, labels = seq(1900,2020,10))
ggsave(filename = paste0('./_Organized_Figs/',"date_intervals.pdf"), width = 11, height = 10)




# A = DATES_LONG
# A = A %>% dplyr::group_by(DISEASES) %>% dplyr::mutate(num_sources = length(SOURCES))
# A = A[A$dates_upper > as.Date('01-01-2010', format = '%m-%d-%Y'),]
# A = A[A$SOURCES != 'Project Tycho',]
# A = A[(A$num_sources > 1)|(A$DISEASES == 'COVID-19'),]
# A = A[!grepl('InfluenzaA',A$DISEASES) & !grepl('InfluenzaB',A$DISEASES),]
# library(dplyr)
# A = A %>% dplyr::group_by(DISEASES) %>% dplyr::mutate(num_sources = length(SOURCES))
# 
# data.frame(A[,c('FILES', 'dates_lower', 'dates_upper')])
# 
# table(A$SOURCES, A$DISEASES)

