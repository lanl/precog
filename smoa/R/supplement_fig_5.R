## Dave Osthus and AC Murph
## 12-17-24
## view embeddings
library(ggplot2)
library(data.table)
library(umap)
library(dplyr)
library(ggExtra)
library(gridExtra)
library(this.path)
library(grid)
setwd(paste0(this.path::here(),"/../"))
theme_set(theme_bw())

savepath <- "data/embed_synthetic_w_data"
filepath_covid <- "data/forecasting_23-26/embed_synthetic_w_data/"
filepath_other <- "data/forecasting_23-26/embed_synthetic_w_data/embeddings_for_murph/"
realpath <- "data/forecasting_23-26/embed_synthetic_w_data/embeddings_for_murph/"
load("data/forecasting_23-26/embed_synthetic_w_data/synthetic_X.RData")

## pick 100k rows of synthetic data
set.seed(75600)
nsample <- 10000
set.seed(534300)
X <- X[sample(1:nrow(X),nsample,replace=F),]
synthetic_summary <- data.frame(type = "synthetic",
                                id = 1:nrow(X),
                                avg = rowMeans(X))
smallX <- data.frame(t(apply(X,1,diff)))
names(smallX) <- paste0("X",1:ncol(smallX))
rm(X)


#########################################################################
## load chik data
lf <- list.files(realpath)
chik_lf <- lf[grepl("Chikungunya", lf) & grepl("_X", lf)]
chikX <- NULL
for(i in 1:length(chik_lf)){
  tempdf <- fread(paste0(realpath,chik_lf[i]))
  chikX <- rbind(chikX, tempdf)
}

## remove negatives
chikX[chikX < 0] <- 0

## compute chik summary
chik_summary <- data.frame(type = "Chikungunya",
                           id = 1:nrow(chikX),
                           avg = rowMeans(chikX))

## get differences
chikXdiff <- data.frame(t(apply(chikX,1,diff)))
names(chikXdiff) <- paste0("X",1:ncol(chikXdiff))

## subset chikXdiff 
if(nrow(chikXdiff) > nsample){
  set.seed(9878)
  chikXdiff <- chikXdiff[sample(1:nrow(chikXdiff), nsample, replace = F),]
}

## combine with synthetic
chikdf <- as.matrix(rbind(smallX, chikXdiff))

## fit UMAP (2-dimensions)
set.seed(1128)
chik_umap <- umap(chikdf, n_components = 2)
chik_umap_df <- data.frame(chik_umap$layout)
names(chik_umap_df) <- paste0("X",1:ncol(chik_umap_df))
chik_umap_df$type <- "Chikungunya"
chik_umap_df[1:nrow(smallX),]$type <- "Synthetic"


#########################################################################
# load covid data
covid_lf <- list.files('data/embeddings_by_location')
covidX <- NULL
for(i in 1:length(covid_lf)){
  load(paste0('data/embeddings_by_location/',covid_lf[i]))
  covidX <- rbind(covidX, state_X)
}

## remove negatives
covidX[covidX < 0] <- 0

## compute covid summary
covid_summary <- data.frame(type = "COVID-19",
                            id = 1:nrow(covidX),
                            avg = rowMeans(covidX))

## get differences
covidXdiff <- data.frame(t(apply(covidX,1,diff)))
names(covidXdiff) <- paste0("X",1:ncol(covidXdiff))

## subset covidXdiff 
if(nrow(covidXdiff) > nsample){
  set.seed(9878)
  covidXdiff <- covidXdiff[sample(1:nrow(covidXdiff), nsample, replace = F),]
}

## combine with synthetic
coviddf <- as.matrix(rbind(smallX, covidXdiff))

## fit UMAP (2-dimensions)
set.seed(1128)
covid_umap <- umap(coviddf, n_components = 2)

## reformat
covid_umap_df <- data.frame(covid_umap$layout)
names(covid_umap_df) <- paste0("X",1:ncol(covid_umap_df))
covid_umap_df$type <- "COVID"
covid_umap_df[1:nrow(smallX),]$type <- "Synthetic"


#########################################################################
## load influenza ushhs data
ushhs_lf <- lf[grepl("ushhs", lf) & grepl("_X", lf)]
ushhsX <- NULL
for(i in 1:length(ushhs_lf)){
  tempdf <- fread(paste0(realpath,ushhs_lf[i]))
  ushhsX <- rbind(ushhsX, tempdf)
}

## remove negatives
ushhsX[ushhsX < 0] <- 0

## compute ushhs summary
ushhs_summary <- data.frame(type = "Influenza Hospitalizations",
                            id = 1:nrow(ushhsX),
                            avg = rowMeans(ushhsX))


## get differences
ushhsXdiff <- data.frame(t(apply(ushhsX,1,diff)))
names(ushhsXdiff) <- paste0("X",1:ncol(ushhsXdiff))

## subset ushhsXdiff 
if(nrow(ushhsXdiff) > nsample){
  set.seed(9878)
  ushhsXdiff <- ushhsXdiff[sample(1:nrow(ushhsXdiff), nsample, replace = F),]
}

## combine with synthetic
ushhsdf <- as.matrix(rbind(smallX, ushhsXdiff))

## fit UMAP (2-dimensions)
set.seed(1128)
ushhs_umap <- umap(ushhsdf, n_components = 2)

## reformat
ushhs_umap_df <- data.frame(ushhs_umap$layout)
names(ushhs_umap_df) <- paste0("X",1:ncol(ushhs_umap_df))
ushhs_umap_df$type <- "Influenza Hospitalizations"
ushhs_umap_df[1:nrow(smallX),]$type <- "Synthetic"


#########################################################################
## load influenza usflunet data
usflunet_lf <- lf[grepl("usflunet", lf) & grepl("_X", lf)]
usflunetX <- NULL
for(i in 1:length(usflunet_lf)){
  tempdf <- fread(paste0(realpath,usflunet_lf[i]))
  usflunetX <- rbind(usflunetX, tempdf)
}

## remove negatives
usflunetX[usflunetX < 0] <- 0

## compute usflunet summary
usflunet_summary <- data.frame(type = "ILI Incidence",
                               id = 1:nrow(usflunetX),
                               avg = rowMeans(usflunetX))


## get differences
usflunetXdiff <- data.frame(t(apply(usflunetX,1,diff)))
names(usflunetXdiff) <- paste0("X",1:ncol(usflunetXdiff))

## subset usflunetXdiff 
if(nrow(usflunetXdiff) > nsample){
  set.seed(9878)
  usflunetXdiff <- usflunetXdiff[sample(1:nrow(usflunetXdiff), nsample, replace = F),]
}

## combine with synthetic
usflunetdf <- as.matrix(rbind(smallX, usflunetXdiff))

## fit UMAP (2-dimensions)
set.seed(1128)
usflunet_umap <- umap(usflunetdf, n_components = 2)

## reformat
usflunet_umap_df <- data.frame(usflunet_umap$layout)
names(usflunet_umap_df) <- paste0("X",1:ncol(usflunet_umap_df))
usflunet_umap_df$type <- "ILI Incidence"
usflunet_umap_df[1:nrow(smallX),]$type <- "Synthetic"


#########################################################################
## load dengue data
dengue_lf <- lf[grepl("Dengue", lf) & grepl("_X", lf)]
dengueX <- NULL
for(i in 1:length(dengue_lf)){
  tempdf <- fread(paste0(realpath,dengue_lf[i]))
  dengueX <- rbind(dengueX, tempdf)
}

## remove negatives
dengueX[dengueX < 0] <- 0

## compute dengue summary
dengue_summary <- data.frame(type = "Dengue",
                             id = 1:nrow(dengueX),
                             avg = rowMeans(dengueX))


## get differences
dengueXdiff <- data.frame(t(apply(dengueX,1,diff)))
names(dengueXdiff) <- paste0("X",1:ncol(dengueXdiff))

## subset dengueXdiff 
if(nrow(dengueXdiff) > nsample){
  set.seed(9878)
  dengueXdiff <- dengueXdiff[sample(1:nrow(dengueXdiff), nsample, replace = F),]
}

## combine with synthetic
denguedf <- as.matrix(rbind(smallX, dengueXdiff))

## fit UMAP (2-dimensions)
set.seed(1128)
dengue_umap <- umap(denguedf, n_components = 2)

## reformat
dengue_umap_df <- data.frame(dengue_umap$layout)
names(dengue_umap_df) <- paste0("X",1:ncol(dengue_umap_df))
dengue_umap_df$type <- "Dengue"
dengue_umap_df[1:nrow(smallX),]$type <- "Synthetic"


###########################################################################
#### make FIGURE 5
###########################################################################

## plot synthetic on top of real
patall <- grid.arrange(
  ggplot()+
    geom_point(aes(x=X1, y=X2), data=subset(covid_umap_df, type != "Synthetic"), color = I('#d95f02'))+
    geom_point(aes(x=X1, y=X2), data=subset(covid_umap_df, type == "Synthetic"), color = I("black"))+
    ggtitle("COVID-19")+
    xlim(range(covid_umap_df$X1))+
    ylim(range(covid_umap_df$X2)) +
    theme(plot.title = element_text(hjust = 0.5)),
  ggplot()+
    geom_point(aes(x=X1, y=X2), data=subset(dengue_umap_df, type != "Synthetic"), color = I('#66a61e'))+
    geom_point(aes(x=X1, y=X2), data=subset(dengue_umap_df, type == "Synthetic"), color = I("black"))+
    ggtitle("Dengue")+
    xlim(range(dengue_umap_df$X1))+
    ylim(range(dengue_umap_df$X2)) +
    theme(plot.title = element_text(hjust = 0.5)),
  ggplot()+
    geom_point(aes(x=X1, y=X2), data=subset(usflunet_umap_df, type != "Synthetic"), color = I('#7570b3'))+
    geom_point(aes(x=X1, y=X2), data=subset(usflunet_umap_df, type == "Synthetic"), color = I("black"))+
    ggtitle("ILI Incidence")+
    xlim(range(usflunet_umap_df$X1))+
    ylim(range(usflunet_umap_df$X2)) +
    theme(plot.title = element_text(hjust = 0.5)),
  ggplot()+
    geom_point(aes(x=X1, y=X2), data=subset(ushhs_umap_df, type != "Synthetic"), color = I('#e6ab02'))+
    geom_point(aes(x=X1, y=X2), data=subset(ushhs_umap_df, type == "Synthetic"), color = I("black"))+
    ggtitle("Influenza Hospitalizations")+
    xlim(range(ushhs_umap_df$X1))+
    ylim(range(ushhs_umap_df$X2)) +
    theme(plot.title = element_text(hjust = 0.5)),
  ggplot()+
    geom_point(aes(x=X1, y=X2), data=subset(chik_umap_df, type != "Synthetic"), color = I('#e7298a'))+
    geom_point(aes(x=X1, y=X2), data=subset(chik_umap_df, type == "Synthetic"), color = I("black"))+
    ggtitle("Chikungunya")+
    xlim(range(chik_umap_df$X1))+
    ylim(range(chik_umap_df$X2)) +
    theme(plot.title = element_text(hjust = 0.5)),ncol=1, top=textGrob("UMAP: Synthetic over Real"))


## plot real on top of synthetic
pbtall <- grid.arrange(
  ggplot()+
    geom_point(aes(x=X1, y=X2), data=subset(covid_umap_df, type == "Synthetic"), color = I("black"))+
    geom_point(aes(x=X1, y=X2), data=subset(covid_umap_df, type != "Synthetic"), color = I('#d95f02'))+
    ggtitle("COVID-19")+
    xlim(range(covid_umap_df$X1))+
    ylim(range(covid_umap_df$X2)) +
    theme(plot.title = element_text(hjust = 0.5)),
  ggplot()+
    geom_point(aes(x=X1, y=X2), data=subset(dengue_umap_df, type == "Synthetic"), color = I("black"))+
    geom_point(aes(x=X1, y=X2), data=subset(dengue_umap_df, type != "Synthetic"), color = I('#66a61e'))+
    ggtitle("Dengue")+
    xlim(range(dengue_umap_df$X1))+
    ylim(range(dengue_umap_df$X2)) +
    theme(plot.title = element_text(hjust = 0.5)),
  ggplot()+
    geom_point(aes(x=X1, y=X2), data=subset(usflunet_umap_df, type == "Synthetic"), color = I("black"))+
    geom_point(aes(x=X1, y=X2), data=subset(usflunet_umap_df, type != "Synthetic"), color = I('#7570b3'))+
    ggtitle("ILI Incidence")+
    xlim(range(usflunet_umap_df$X1))+
    ylim(range(usflunet_umap_df$X2)) +
    theme(plot.title = element_text(hjust = 0.5)),
  ggplot()+
    geom_point(aes(x=X1, y=X2), data=subset(ushhs_umap_df, type == "Synthetic"), color = I("black"))+
    geom_point(aes(x=X1, y=X2), data=subset(ushhs_umap_df, type != "Synthetic"), color = I('#e6ab02'))+
    ggtitle("Influenza Hospitalizations")+
    xlim(range(ushhs_umap_df$X1))+
    ylim(range(ushhs_umap_df$X2)) +
    theme(plot.title = element_text(hjust = 0.5)),
  ggplot()+
    geom_point(aes(x=X1, y=X2), data=subset(chik_umap_df, type == "Synthetic"), color = I("black"))+
    geom_point(aes(x=X1, y=X2), data=subset(chik_umap_df, type != "Synthetic"), color = I('#e7298a'))+
    ggtitle("Chikungunya")+
    xlim(range(chik_umap_df$X1))+
    ylim(range(chik_umap_df$X2)) +
    theme(plot.title = element_text(hjust = 0.5)),ncol=1, top=textGrob("UMAP: Real over Synthetic"))


## average of segment
pctall <- grid.arrange(
  ggplot() +
    geom_histogram(data = synthetic_summary, aes(x = avg, y = ..density..), 
                   fill = "black",  bins = 150, position = "identity") +
    geom_histogram(data = covid_summary, aes(x = avg, y = ..density..), 
                   fill = "#d95f02", alpha = 0.75, bins = 150, position = "identity") +
    scale_x_log10()+
    labs(x = "Average of Segment", y = "Density") +
    ggtitle("COVID-19")+
    theme(plot.title = element_text(hjust = 0.5)),
  ##
  ggplot() +
    geom_histogram(data = synthetic_summary, aes(x = avg, y = ..density..), 
                   fill = "black",  bins = 150, position = "identity") +
    geom_histogram(data = dengue_summary, aes(x = avg, y = ..density..), 
                   fill = "#66a61e", alpha = 0.75, bins = 150, position = "identity") +
    scale_x_log10()+
    labs(x = "Average of Segment", y = "Density") +
    ggtitle("Dengue")+
    theme(plot.title = element_text(hjust = 0.5)),
  ##
  ggplot() +
    geom_histogram(data = synthetic_summary, aes(x = avg, y = ..density..), 
                   fill = "black",  bins = 150, position = "identity") +
    geom_histogram(data = usflunet_summary, aes(x = avg, y = ..density..), 
                   fill = "#7570b3", alpha = 0.75, bins = 150, position = "identity") +
    scale_x_log10()+
    labs(x = "Average of Segment", y = "Density") +
    ggtitle("ILI Incidence")+
    theme(plot.title = element_text(hjust = 0.5)),
  ##
  ggplot() +
    geom_histogram(data = synthetic_summary, aes(x = avg, y = ..density..), 
                   fill = "black",  bins = 150, position = "identity") +
    geom_histogram(data = ushhs_summary, aes(x = avg, y = ..density..), 
                   fill = "#e6ab02", alpha = 0.75, bins = 150, position = "identity") +
    scale_x_log10()+
    labs(x = "Average of Segment", y = "Density") +
    ggtitle("Influenza Hospitalizations")+
    theme(plot.title = element_text(hjust = 0.5)),
  ##
  ggplot() +
    geom_histogram(data = synthetic_summary, aes(x = avg, y = ..density..), 
                   fill = "black",  bins = 150, position = "identity") +
    geom_histogram(data = chik_summary, aes(x = avg, y = ..density..), 
                   fill = "#e7298a", alpha = 0.75, bins = 150, position = "identity") +
    scale_x_log10()+
    labs(x = "Average of Segment", y = "Density") +
    ggtitle("Chikungunya")+
    theme(plot.title = element_text(hjust = 0.5)),ncol=1, top=textGrob("Average of Segment"))


## combine the plots
grid.arrange(patall, pbtall, pctall, nrow=1)