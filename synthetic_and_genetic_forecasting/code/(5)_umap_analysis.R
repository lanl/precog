
## Perform UMAP on the train/test data to see if the synthetic training data is "closer" to the test data than the real training data.

library(data.table)
library(matrixStats)
library(ggplot2)
library(patchwork)
library(uwot)
library(lightgbm)
library(lubridate)


# ----------------------------------------------------------
# 0) Knobs
# ----------------------------------------------------------
datapath <- paste0(here::here("synthetic_and_genetic_forecasting", "data_clean"), "/") 
configpath <- paste0(here::here("synthetic_and_genetic_forecasting", "trained_models"), "/") 
savepath <- paste0(here::here("synthetic_and_genetic_forecasting", "figs"), "/") 


# ----------------------------------------------------------
# read in a config file
cfg <- readRDS(file = paste0(configpath,"cfg_all.rds"))

# get context length
context_len <- cfg$context_length_in


# ----------------------------------------------------------
# 0.5) Load your lists
# ----------------------------------------------------------

# window length
M <- context_len
n_train <- 150000
n_valid <- 25000
n_test  <- 8000

# helper: sliding windows as rows, in ascending order within each row
make_windows <- function(x, M) {
  if (length(x) < M) return(matrix(numeric(0), nrow = 0, ncol = M))
  w <- stats::embed(x, M)       # rows = length(x) - M + 1; cols = M (reversed)
  w[, M:1, drop = FALSE]        # reverse columns to get ascending order within each row
}

## row standardize W
row_standardize_fast <- function(W, na.rm = FALSE) {
  stopifnot(is.matrix(W))
  storage.mode(W) <- "double"              # be explicit
  
  mu <- rowMeans2(W, na.rm = na.rm)
  # compute row sds relative to the SAME means to avoid any mismatch
  s  <- rowSds(W, center = mu, na.rm = na.rm)
  
  nz <- s > 0 & is.finite(s)
  
  if (any(nz)) {
    # direct row-wise broadcasting; base R recycles a length-rows vector down columns
    W[nz, ] <- (W[nz, , drop = FALSE] - mu[nz]) / s[nz]
  }
  W
}



# ----------------------------------------------------------
# Non-COVID-19, real respiratory data
# ----------------------------------------------------------
dfresp <- readRDS(paste0(datapath,"train_real.RDS"))
length(dfresp)
noncovidresp <- which(unlist(lapply(dfresp, function(x){return(x$ts_disease[1])})) != "covid")
dfresp <- dfresp[noncovidresp]
length(dfresp)
mats_resp <- lapply(dfresp, function(dt) {
  setorder(dt, time)       # in-place, ascending
  make_windows(dt$cases, M)
})
rm(dfresp)

# rbind all series together
W_resp_full <- do.call(rbind, mats_resp)
rm(mats_resp)
invisible(gc()) 

## break out ids
resp_id       <- sort(sample(1:nrow(W_resp_full), n_train, replace = F))
resp_id_valid <- sort(sample(setdiff(1:nrow(W_resp_full),resp_id), n_valid, replace = F))
resp_id_test  <- sort(sample(setdiff(1:nrow(W_resp_full),c(resp_id_valid,resp_id)), n_test, replace = F))

Xtrain_resp <- row_standardize_fast(W_resp_full[resp_id,])
Xvalid_resp <- row_standardize_fast(W_resp_full[resp_id_valid,])
Xtest_resp  <- row_standardize_fast(W_resp_full[resp_id_test,])




# ----------------------------------------------------------
# Total cases (tc) Synthetic
# ----------------------------------------------------------
dfsyn_tc <- readRDS(paste0(datapath,"train_syn_tc.RDS"))
mats_syn_tc <- lapply(dfsyn_tc, function(dt) {
  setorder(dt, time)       # in-place, ascending
  make_windows(dt$cases, M)
})
rm(dfsyn_tc)
invisible(gc())

# rbind all series together
W_syn_tc_full <- do.call(rbind, mats_syn_tc)
rm(mats_syn_tc)
invisible(gc())

## break out ids
syn_tc_id <- sort(sample(1:nrow(W_syn_tc_full), n_train, replace = F))
syn_tc_id_valid <- sort(sample(setdiff(1:nrow(W_syn_tc_full),syn_tc_id), n_valid, replace = F))
syn_tc_id_test <- sort(sample(setdiff(1:nrow(W_syn_tc_full),c(syn_tc_id_valid,syn_tc_id)), n_test, replace = F))

Xtrain_syn_tc <- row_standardize_fast(W_syn_tc_full[syn_tc_id,])
Xvalid_syn_tc <- row_standardize_fast(W_syn_tc_full[syn_tc_id_valid,])
Xtest_syn_tc  <- row_standardize_fast(W_syn_tc_full[syn_tc_id_test,])



# ----------------------------------------------------------
# US state COVID-19 
# ----------------------------------------------------------
df    <- readRDS(paste0(datapath,"test_covid.RDS"))

## remove all dates in in the analysis range
min_date <- ymd("2020-06-01") - 7*20
max_date <- ymd("2022-12-31")

df <- lapply(df, function(dt) {
  dt <- copy(dt)                 # optional: avoid modifying originals
  dt[, date := as.IDate(date)]   # ensure proper type
  dt[date %between% c(min_date, max_date)]
})


mats_covid <- lapply(df, function(dt) {
  setorder(dt, time)       # in-place, ascending
  make_windows(dt$cases, M)
})

## I need to make a state/end date look up table for the Xtrain_covid inputs
mats_covid_dates <- lapply(df, function(dt) {
  setorder(dt, time)       # in-place, ascending
  data.table(series_id = dt$series_id[1],
             time = make_windows(dt$time, M)[,M])
})

# make the lookup tble for matrices
lookup_covid1 <- do.call(rbind, mats_covid_dates)
lookup_covid1$rowid <- 1:nrow(lookup_covid1)

## date lookup
mats_covid_realdates <- lapply(df, function(dt){
  setorder(dt, time)       # in-place, ascending
  data.table(series_id = dt$series_id[1],
            date = dt$date,
            avg_cases = mean(dt$cases),
            max_cases = max(dt$cases),
            cases = dt$cases,
            time = dt$time)
})

# make the lookup tble for matrices
lookup_covid2 <- do.call(rbind, mats_covid_realdates)

## merge
lookup_covid <- merge(lookup_covid1, lookup_covid2, by=c("series_id","time"), all.y=T)


rm(df)
invisible(gc())





# rbind all series together
W_covid <- do.call(rbind, mats_covid)
rm(mats_covid)
invisible(gc())

Xtrain_covid <- row_standardize_fast(W_covid)

# ----------------------------------------------------------
# Put the real, syn_tc, syn_vac, and test_covid together
# ----------------------------------------------------------
Wz <- rbind(Xtrain_resp, Xtrain_syn_tc, Xtrain_covid)
W_type <- c(rep(1,nrow(Xtrain_resp)),
            rep(2,nrow(Xtrain_syn_tc)),
            rep(3,nrow(Xtrain_covid)))



# ----------------------------------------------------------
# Perform UMAP on Wz
# ----------------------------------------------------------

# Wz: your row-standardized matrix (rows = windows, cols = M)
stopifnot(is.matrix(Wz))
storage.mode(Wz) <- "double"


# choose a fast ANN method
nn_method <- if (requireNamespace("RcppHNSW", quietly = TRUE)) "hnsw" else "annoy"

# sensible defaults; tweak as needed
n_neighbors <- 15
min_dist    <- 0.1
metric      <- "euclidean"  # try "cosine" if you prefer angle similarity
myseed <- 560

n_threads <- max(1L, parallel::detectCores() - 1L)

Wz <- as.matrix(Wz)  # your input matrix

## 1) Make a stable key per row (round if you have tiny numeric noise)
row_key <- apply(round(Wz, 5), 1, paste, collapse = "\r")

keep1 <- !duplicated(row_key)
keep2 <- (apply(Wz,1,sd) > 0)
keep  <- (keep1 == TRUE & keep2 == TRUE)

Wz_unique <- Wz[keep, , drop = FALSE]
W_type_unique <- W_type[keep]

set.seed(42)  # R-side RNG
um_model <- umap(
  Wz_unique,
  n_components = 2,
  n_neighbors  = n_neighbors,
  min_dist     = min_dist,
  metric       = metric,
  nn_method    = nn_method,
  n_threads    = n_threads,
  ret_model    = TRUE,
  verbose      = TRUE,
  seed         = myseed
)

# 2) allocate full embedding and fill training part
U <- um_model$embedding
colnames(U) <- c("UMAP1", "UMAP2")
U <- data.table(U)
U$type <- as.numeric(W_type[keep])


# ------------------------------------------------------------------------------

traincolors <- data.frame(variable = c("real","syn_tc","covid"),
                          type = c(1,2,3),
                         color = c("#ef7271","#7baacd","#86c368"))


plot_bins <- 150

p_real <- ggplot(U, aes(UMAP1, UMAP2)) +
  stat_bin2d(bins = plot_bins, alpha = I(1), fill = I("white"), color = I("white")) +
  stat_bin2d(aes(alpha = after_stat(count)), bins = plot_bins, fill = traincolors[traincolors$type == 1,]$color, data = U[U$type == 1, ])+
  scale_alpha(range = c(0.5, 2), guide = "none") +
  guides(fill = "none") +
  theme_minimal(base_size = 12) +
  theme(
    plot.background  = element_rect(fill = "black", color = NA),
    panel.background = element_rect(fill = "black", color = NA),
    panel.grid       = element_line(color = "grey20"),
    axis.text        = element_text(color = "grey85"),
    axis.title       = element_text(color = "grey85"),
    strip.text       = element_text(color = "grey85"),
    legend.position  = "none"
  )+
  ggtitle(traincolors[traincolors$type == 1,]$variable) 

p_syn_tc <- ggplot(U, aes(UMAP1, UMAP2)) +
  stat_bin2d(bins = plot_bins, alpha = I(1), fill = I("white"), color = I("white")) +
  stat_bin2d(aes(alpha = after_stat(count)), bins = plot_bins, fill = traincolors[traincolors$type == 2,]$color, data = U[U$type == 2, ])+
  scale_alpha(range = c(0.5, 2), guide = "none") +
  guides(fill = "none") +
  theme_minimal(base_size = 12) +
  theme(
    plot.background  = element_rect(fill = "black", color = NA),
    panel.background = element_rect(fill = "black", color = NA),
    panel.grid       = element_line(color = "grey20"),
    axis.text        = element_text(color = "grey85"),
    axis.title       = element_text(color = "grey85"),
    strip.text       = element_text(color = "grey85"),
    legend.position  = "none"
  )+
  ggtitle(traincolors[traincolors$type == 2,]$variable) 

p_test_covid <- ggplot(U, aes(UMAP1, UMAP2)) +
  stat_bin2d(bins = plot_bins, alpha = I(1), fill = I("white"), color = I("white")) +
  stat_bin2d(aes(alpha = after_stat(count)), bins = plot_bins, fill = traincolors[traincolors$type == 3,]$color, data = U[U$type == 3, ])+
  scale_alpha(range = c(0.5, 2), guide = "none") +
  guides(fill = "none") +
  theme_minimal(base_size = 12) +
  theme(
    plot.background  = element_rect(fill = "black", color = NA),
    panel.background = element_rect(fill = "black", color = NA),
    panel.grid       = element_line(color = "grey20"),
    axis.text        = element_text(color = "grey85"),
    axis.title       = element_text(color = "grey85"),
    strip.text       = element_text(color = "grey85"),
    legend.position  = "none"
  )+
  ggtitle(traincolors[traincolors$type == 4,]$variable)






# -----------------------------------------------------------
# Version 2
# -----------------------------------------------------------

p_real2 <- ggplot(U, aes(UMAP1, UMAP2)) +
  stat_bin2d(bins = plot_bins, alpha = I(1), fill = I("white"), color = I("white")) +
  stat_bin2d(aes(alpha = after_stat(count)), bins = plot_bins, fill = traincolors[traincolors$type == 1,]$color, data = U[U$type == 1, ])+
  stat_bin2d(aes(alpha = after_stat(count)), bins = plot_bins, fill = traincolors[traincolors$type == 4,]$color, data = U[U$type == 4, ])+
  scale_alpha(range = c(0.5, 2), guide = "none") +
  guides(fill = "none") +
  theme_minimal(base_size = 12) +
  theme(
    plot.background  = element_rect(fill = "black", color = NA),
    panel.background = element_rect(fill = "black", color = NA),
    panel.grid       = element_line(color = "grey20"),
    axis.text        = element_text(color = "grey85"),
    axis.title       = element_text(color = "grey85"),
    strip.text       = element_text(color = "grey85"),
    legend.position  = "none"
  )+
  ggtitle(traincolors[traincolors$type == 1,]$variable) 

p_syn_tc2 <- ggplot(U, aes(UMAP1, UMAP2)) +
  stat_bin2d(bins = plot_bins, alpha = I(1), fill = I("white"), color = I("white")) +
  stat_bin2d(aes(alpha = after_stat(count)), bins = plot_bins, fill = traincolors[traincolors$type == 2,]$color, data = U[U$type == 2, ])+
  stat_bin2d(aes(alpha = after_stat(count)), bins = plot_bins, fill = traincolors[traincolors$type == 4,]$color, data = U[U$type == 4, ])+
  scale_alpha(range = c(0.5, 2), guide = "none") +
  guides(fill = "none") +
  theme_minimal(base_size = 12) +
  theme(
    plot.background  = element_rect(fill = "black", color = NA),
    panel.background = element_rect(fill = "black", color = NA),
    panel.grid       = element_line(color = "grey20"),
    axis.text        = element_text(color = "grey85"),
    axis.title       = element_text(color = "grey85"),
    strip.text       = element_text(color = "grey85"),
    legend.position  = "none"
  )+
  ggtitle(traincolors[traincolors$type == 2,]$variable) 

p_test_covid2 <- ggplot(U, aes(UMAP1, UMAP2)) +
  stat_bin2d(bins = plot_bins, alpha = I(1), fill = I("white"), color = I("white")) +
  stat_bin2d(aes(alpha = after_stat(count)), bins = plot_bins, fill = traincolors[traincolors$type == 3,]$color, data = U[U$type == 3, ])+
  scale_alpha(range = c(0.5, 2), guide = "none") +
  guides(fill = "none") +
  theme_minimal(base_size = 12) +
  theme(
    plot.background  = element_rect(fill = "black", color = NA),
    panel.background = element_rect(fill = "black", color = NA),
    panel.grid       = element_line(color = "grey20"),
    axis.text        = element_text(color = "grey85"),
    axis.title       = element_text(color = "grey85"),
    strip.text       = element_text(color = "grey85"),
    legend.position  = "none"
  )+
  ggtitle(traincolors[traincolors$type == 3,]$variable)







# -------------------------------------------------------------
## TRAIN A CLASSIFICATION MODEL TO SEE IF COVID EXAMPLES "LOOK MORE LIKE" SYNTHETIC EXAMPLES OR REAL EXAMPLES
# install.packages("lightgbm")   # if needed (use the precompiled binary for your OS)
# -------------------------------------------------------------


## make training data for LGBM
Xtrain <- rbind(Xtrain_syn_tc, Xtrain_resp)
Xvalid <- rbind(Xvalid_syn_tc, Xvalid_resp)
Ytrain <- c(rep(1, nrow(Xtrain_syn_tc)),
            rep(0, nrow(Xtrain_resp)))
Yvalid <- c(rep(1, nrow(Xvalid_syn_tc)),
            rep(0, nrow(Xvalid_resp)))
Xtest_resp  <- Xtest_resp
Xtest_syn_tc   <- Xtest_syn_tc
Xtest_covid <- Xtrain_covid



## --- prepare data ---
Xtrain_m <- as.matrix(Xtrain)
Xvalid_m <- as.matrix(Xvalid)

y <- as.numeric(Ytrain)
yvalid <- as.numeric(Yvalid)

dtrain <- lgb.Dataset(data = Xtrain_m, label = y)
dvalid <- lgb.Dataset(data = Xvalid_m, label = yvalid)

## --- model params (sensible defaults) ---
params <- list(
  objective = "binary",
  metric = c("binary_logloss"),
  learning_rate = 0.05,
  num_leaves = 64,
  feature_fraction = 0.9,
  bagging_fraction = 0.8,
  bagging_freq = 1,
  min_data_in_leaf = 20,
  max_depth = -1,
  verbose = 1
)

set.seed(123)
## Train for a fixed number of rounds (simple & robust)
model <- lgb.train(
  params  = params,
  data    = dtrain,
  valids = list(val = dvalid),
  early_stopping_rounds = 500,
  nrounds = 30000,
  verbose = 1
)


# -----------------------------------------------------------
# --- predict probabilities on non-COVID-19, real data
# -----------------------------------------------------------
umaptest_resp <- data.frame(umap_transform(
  Xtest_resp,
  um_model,
  verbose   = TRUE,
  seed = myseed
))
names(umaptest_resp) <- c("UMAP1","UMAP2")
dftest_resp <- data.frame(Xtest_resp)
dftest_resp <- cbind(dftest_resp,
                     umaptest_resp)
dfpred_resp <- data.table(predict(model, Xtest_resp))
names(dfpred_resp) <- c("prob_syn_tc")
dftest_resp<- data.table(dftest_resp, dfpred_resp)
dftest_resp$type <- "Non-COVID-19 Real Respiratory"


# -----------------------------------------------------------
# --- predict probabilities on synthetic total cases
# -----------------------------------------------------------
umaptest_syn_tc <- data.frame(umap_transform(
  Xtest_syn_tc,
  um_model,
  verbose   = TRUE,
  seed = myseed
))
names(umaptest_syn_tc) <- c("UMAP1","UMAP2")
dftest_syn_tc <- data.frame(Xtest_syn_tc)
dftest_syn_tc <- cbind(dftest_syn_tc,
                       umaptest_syn_tc)
dfpred_syn_tc <- data.table(predict(model, Xtest_syn_tc))
names(dfpred_syn_tc) <- c("prob_syn_tc")
dftest_syn_tc<- data.table(dftest_syn_tc, dfpred_syn_tc)
dftest_syn_tc$type <- "Synthetic Total Cases"


# -----------------------------------------------------------
# --- predict probabilities on COVID test set
# -----------------------------------------------------------
umaptest_covid <- data.frame(umap_transform(
  Xtest_covid,
  um_model,
  verbose   = TRUE,
  seed = myseed
))
names(umaptest_covid) <- c("UMAP1","UMAP2")
dftest_covid <- data.frame(Xtest_covid)
dftest_covid <- cbind(dftest_covid,
                      umaptest_covid)
dfpred_covid <- data.table(predict(model, Xtest_covid))
names(dfpred_covid) <- c("prob_syn_tc")
dftest_covid<- data.table(dftest_covid, dfpred_covid)
dftest_covid$type <- "COVID-19"


# -----------------------------------------------------------
# --- combine predictions
# -----------------------------------------------------------

dftest_resp$rowid <- 1:nrow(dftest_resp)
dftest_syn_tc$rowid <- 1:nrow(dftest_syn_tc)
dftest_covid$rowid <- 1:nrow(dftest_covid)

dftest <- rbind(dftest_resp,
                dftest_syn_tc,
                dftest_covid)

dftest$type <- factor(as.factor(dftest$type), levels = c("Non-COVID-19 Real Respiratory","COVID-19","Synthetic Total Cases"))

## some quick summaries of classifier performance:
dftest[,list(n_syn = sum(prob_syn_tc > .5), prob_syn = mean(prob_syn_tc > .5), n = length(UMAP1)), by="type"]

## plot the results
dftest_melt <- data.table(reshape2::melt(subset(dftest, select = grep("X",names(dftest), value = T, invert = T)),
                              id.vars = c("type","UMAP1","UMAP2","rowid")))


# plot results
setEPS()
postscript(paste0(savepath,"prob_syn_tc_by_data_types.eps"), height = 5, width = 10)
(ggplot(dftest_melt, aes(x = value)) +
    geom_histogram(
      aes(y = after_stat(density), fill = after_stat(x)),
      binwidth = 0.05,
      boundary = 0,
      color = "black"
    ) +
    scale_fill_gradient2(
      low = "#e31a1c",
      mid = "white",
      high = "#1f78b4",
      midpoint = 0.5,
      name="Probability\nSynthetic\nTotal\nCases",
      limits = c(0, 1)
    ) +
    scale_x_continuous(limits = c(0, 1)) +
    ylab("Density") +
    xlab("Probability Input is Synthetic TC") +
    facet_wrap(~type, nrow = 1))/
(qplot(UMAP1, UMAP2,  data = subset(dftest_melt), color = value)+
   scale_color_gradient2(
    low = "#e31a1c",      # min color
    mid = "white",     # midpoint color
    high = "#1f78b4",# max color
    midpoint = 0.5,
    limits = c(0,1),
    name="Probability\nSynthetic\nTotal\nCases")+
   facet_wrap( ~ type, nrow = 1))
dev.off()



# ---------------- title ----------------
state_from_series_id <- function(x) {
  if (identical(x, "all")) return("US")
  core <- sub("^unitedstates_", "", x)
  
  map <- c(
    newyork = "New York",
    newhampshire = "New Hampshire",
    newjersey = "New Jersey",
    newmexico = "New Mexico",
    northcarolina = "North Carolina",
    northdakota = "North Dakota",
    southcarolina = "South Carolina",
    southdakota = "South Dakota",
    westvirginia = "West Virginia",
    rhodeisland = "Rhode Island"
  )
  
  if (core %in% names(map)) return(unname(map[[core]]))
  core <- gsub("_", " ", core)
  tools::toTitleCase(tolower(core))
}



# lookup_covid$rowid <- 1:nrow(lookup_covid)
dftest_covid2 <- merge(dftest_covid,
                       lookup_covid,
                       by="rowid",
                       all.y=T)

state_order <- dftest_covid2[,list(avg_p = min(prob_syn_tc)),by="series_id"]
state_order <- state_order[order(state_order$avg_p),]
dftest_covid2$series_id <- factor(as.factor(dftest_covid2$series_id), levels = state_order$series_id)


dftest_covid2$state <- NA
for(i in 1:nrow(dftest_covid2)){
  dftest_covid2$state[i] <- state_from_series_id(dftest_covid2$series_id[i])
}


# assuming df$date is Date
start_jan <- as.Date(sprintf("%d-01-01", as.integer(format(min(dftest_covid2$date), "%Y"))))
end_date  <- max(dftest_covid2$date)



## plot it
setEPS()
postscript(paste0(savepath,"prob_syn_tc_allstates.eps"), height = 12, width = 10)
(qplot(date, cases/avg_cases, color = prob_syn_tc, group = state, geom=c("line","point"), data = dftest_covid2)+
  ylab("Cases / Average Cases")+
    xlab("")+
    scale_color_gradient2(
      low = "#e31a1c",      # min color
      mid = "white",     # midpoint color
      high = "#1f78b4",# max color
      midpoint = 0.5,
      na.value = "grey20",
      limits = c(0,1),
      name="Probability Synthetic Total Cases")+
    theme(legend.position="none",
          plot.title = element_text(hjust = 0.5))+
    ggtitle("COVID-19 Cases")+
  scale_x_date(
    breaks = seq(start_jan, end_date, by = "3 months"),
    date_labels = "%b %Y",
    expand = c(0, 0)))/
      
(qplot(date, state, fill=prob_syn_tc, data = dftest_covid2, geom="tile", color=I("black"))+
  scale_fill_gradient2(
    low = "#e31a1c",      # min color
    mid = "white",     # midpoint color
    high = "#1f78b4",# max color
    midpoint = 0.5,
    na.value = "grey20",
    limits = c(0,1),
    name="Probability Synthetic Total Cases")+
  ylab("")+
  xlab("")+
  scale_y_discrete(expand = c(0, 0),
                   limits = rev) +
  scale_x_date(
    breaks = seq(start_jan, end_date, by = "3 months"),
    date_labels = "%b %Y",
    expand = c(0, 0)
  )+
  theme(
    panel.grid = element_blank(),
    legend.position = "bottom",
    legend.box.margin = margin(t = -20, r = 0, b = 0, l = 0),  # pull legend up
    legend.margin = margin(t = 0, r = 0, b = 0, l = 0),
    axis.title.x = element_text(margin = margin(t = 4))        # reduce space above x title
  )) + plot_layout(heights = c(1, 3))
dev.off()
