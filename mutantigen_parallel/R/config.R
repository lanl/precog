## Dave Osthus
## 6-20-25
## 4 quadrants paper: define quantities common to multiple scripts for the analysis

################################################################################
# set the theme for ggplot2
theme_set(theme_bw())


################################################################################
## define quantities need for make_training_data.R

## define max forecast horizon (max steps ahead)
max_horizon <- 4

## define max lag for making features
max_lag <- 9

## define the number of synthetic rows you want made
nrow_synthetic = 500000 

## define the names of the epi columns
epi_cols <- "cases"


################################################################################
## define quantities need for make_persistence_model_forecasts.R
n_samples <- 100000 # should be 100k


################################################################################
## define quantities need for run_forecasts.R

## define the models
models <- paste0("M",1:6)

## synthetic only models
synthetic_only_models <- c("M2","M5")

## node fraction for parallelization
core_frac <- 1/3

## define sequence over dates for analysis
date_vec <- ymd((seq(ymd("2020-05-31"), ymd("2022-12-31"), by = "week")))
# date_vec <- date_vec[length(date_vec)-57]

## define quantiles for analysis
alphavec <- c(.025,.1,.25,.5,.75,.9,.975)

## how many replicates should LightGBM be fit?
nfits <- 10 ## should be 5? 10? 20?
n_lhs <- 15 ## number of lhs design points for hyperparameter selection

## LightGBM parameters
early_stopping_rounds <- 50 
learning_rate <- .05 # should be [0.01, 0.1]
num_leaves <- 62
nrounds = 50000 # should be big
feature_fraction <- 1
bagging_fraction <- .8
bagging_freq <- 1

################################################################################
## define model colors for plotting
m0_color <- "#bebada"
m1_color <- "#8dd3c7"
m2_color <- "#fdb462"
m3_color <- "#e39dac"
m4_color <- "#fb8072"
m5_color <- "#b3de69"
m6_color <- "#80b1d3"
color_vec <- c(m0_color, m1_color, m2_color, m3_color, m4_color, m5_color, m6_color)
color_df <- data.frame(model = paste0("M",0:(length(color_vec)-1)),
                       color = color_vec)
# real_us_epi <- "#d95f02"
# synthetic_epi <- "#1b9e77"


real_us_epi <- "#2b7d2b"
synthetic_epi <- "#762a83"
umap_low <- "#1b2a49"
umap_high <- "#946c00"



