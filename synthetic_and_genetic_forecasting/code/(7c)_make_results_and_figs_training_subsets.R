
# ==================================================================================================
# - Assumes `preds` is already in memory (data.table or data.frame)
# - Assumes paths (e.g., output_root) already defined (only used for fig output)
# - Uses your column names + your model mapping (train_mod + model_input -> M(0)..M(a,v) with a/b splits)
# - Computes: MAE, WIS, interval coverage, relative metrics vs M(0), bootstrap CIs, and key plots
# ==================================================================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(lubridate)
  library(scales)
  library(patchwork)
  library(colorspace)
})
theme_set(theme_bw())

## define paths
fcstpath <- paste0(here::here("synthetic_and_genetic_forecasting", "output"), "/") #"/Users/dosthus/Documents/DR24_PreCog/projects/four_quadrants/quadrants_complete_2-11-26/output/"
figpath <- paste0(here::here("synthetic_and_genetic_forecasting", "figs"), "/") #"/Users/dosthus/Documents/DR24_PreCog/projects/four_quadrants/quadrants_complete_2-11-26/figs/"
datapath <- paste0(here::here("synthetic_and_genetic_forecasting", "data_clean"), "/") #"/Users/dosthus/Documents/DR24_PreCog/projects/four_quadrants/quadrants_complete_2-11-26/data_clean/"
modelpath <- paste0(here::here("synthetic_and_genetic_forecasting", "trained_models_subsets"), "/") #"/Users/dosthus/Documents/DR24_PreCog/projects/four_quadrants/quadrants_complete_2-11-26/trained_models_subsets/"

## Read in test data
test_covid <- readRDS(paste0(datapath,"test_covid.RDS"))

## read in forecasts
lf <- list.files(fcstpath, pattern="_subsets_")
preds <- NULL
for(i in 1:length(lf)){
  preds <- rbind(preds,
                 fread(paste0(fcstpath,lf[i])))
}

# Compute the max cases per series_id
max_cases_df <- data.table(rbindlist(lapply(test_covid, function(x){data.frame(series_id = x$series_id[1], max_cases = max(x$cases))})))

## append the max cases to each state in preds
preds <- merge(preds, 
               max_cases_df,
               by = "series_id",
               all.x=T)

## get training data size
preds$train_size <- as.numeric(sub(".*_(\\d+)$", "\\1", preds$train_mod))
preds[is.na(preds$train_size),]$train_size <- 0
preds$model <- paste0(sub("_[0-9]+$", "", preds$train_mod),"_",preds$model_input)



# --------------------------------------------------------
# Back to preds
# --------------------------------------------------------

## dates that align with the PLOS paper
min_plos_date <- ymd("2020-07-28")
max_plos_date <- ymd("2021-12-21")

## define model colors
MODEL_COLORS <- c(
  # "M(0)" = "darkgrey",
  "M(r,t)" = "#fb9a99",  
  "M(st,t)" = "#a6cee3", 
  # "M(sv,t)" = "#b2df8a", 
  "M(a,t)" = "#cab2d6",  
  "M(r,v)" = "#e31a1c",  
  # "M(st,v)" = "#1f78b4",  
  "M(sv,v)" = "#33a02c",  
  "M(a,v)" = "#6a3d9a"
)

# ------------------------------------------------------------------------------
# Harmonize preds (NEW schema -> old schema expected by the rest of this script)
# ------------------------------------------------------------------------------
preds <- as.data.table(preds)

# sanity: required columns in new schema
req <- c("series_id","time","cases","date","horizon","last_obs","ref_date",
         "train_mod","model_input",
         "q0.0250","q0.1000","q0.2500","q0.5000","q0.7500","q0.9000","q0.9750")
miss <- setdiff(req, names(preds))
if (length(miss)) stop("preds is missing required columns: ", paste(miss, collapse = ", "))

# core renames / aliases used downstream
preds[, `:=`(
  y_true       = as.numeric(cases),
  step_ahead   = as.integer(horizon),
  last_obs_date= as.Date(ref_date),
  persist      = as.numeric(last_obs),
  
  # `end_t` existed in old scripts; safest compatibility alias:
  end_t        = suppressWarnings(as.integer(time))
)]

# rename models
preds[preds$model == "mod_baseline_tc"]$model <- "M(0)"
preds[preds$model == "mod_real_tc"]$model <- "M(r,t)"
preds[preds$model == "mod_real_vac"]$model <- "M(r,v)"
preds[preds$model == "mod_syn_tc_tc"]$model <- "M(st,t)"
preds[preds$model == "mod_syn_tc_vac"]$model <- "M(st,v)"
preds[preds$model == "mod_syn_vac_tc"]$model <- "M(sv,t)"
preds[preds$model == "mod_syn_vac_vac"]$model <- "M(sv,v)"
preds[preds$model == "mod_all_tc"]$model <- "M(a,t)"
preds[preds$model == "mod_all_vac"]$model <- "M(a,v)"

model_order <- c("M(0)","M(r,t)","M(st,t)","M(sv,t)","M(a,t)","M(r,v)","M(st,v)","M(sv,v)","M(a,v)")
preds$model <- factor(as.factor(preds$model), levels = model_order)

# your new quantile columns are q0.0250 etc (dots)
# old code uses q0_025 etc (underscores)
preds[, `:=`(
  q0_025 = as.numeric(q0.0250),
  q0_100 = as.numeric(q0.1000),
  q0_250 = as.numeric(q0.2500),
  q0_500 = as.numeric(q0.5000),
  q0_750 = as.numeric(q0.7500),
  q0_900 = as.numeric(q0.9000),
  q0_975 = as.numeric(q0.9750)
)]

# plos_date flag
preds[, plos_date := as.integer(last_obs_date >= min_plos_date & last_obs_date <= max_plos_date)]

# ------------------------------------------------------------------------------
# Everything below here can remain (almost) unchanged,
# because y_true/step_ahead/last_obs_date/model/q0_* now exist.
# ------------------------------------------------------------------------------

## add in wis
preds$K <- 3
preds$mult <- 1/(preds$K + .5)
preds$w0 <- 0.5
preds$ae <- abs(preds$y_true - preds$q0_500)

preds$w.5 <- .5/2
preds$w.8 <- .2/2
preds$w.95 <- .05/2

preds$IS.5 <- (preds$q0_750 - preds$q0_250) + 
  (2/.5)*ifelse(preds$y_true < preds$q0_250, preds$q0_250 - preds$y_true, 0) + 
  (2/.5)*ifelse(preds$y_true > preds$q0_750, preds$y_true - preds$q0_750, 0)

preds$IS.8 <- (preds$q0_900 - preds$q0_100) + 
  (2/.2)*ifelse(preds$y_true < preds$q0_100, preds$q0_100 - preds$y_true, 0) + 
  (2/.2)*ifelse(preds$y_true > preds$q0_900, preds$y_true - preds$q0_900, 0)

preds$IS.95 <- (preds$q0_975 - preds$q0_025) + 
  (2/.05)*ifelse(preds$y_true < preds$q0_025, preds$q0_025 - preds$y_true, 0) + 
  (2/.05)*ifelse(preds$y_true > preds$q0_975, preds$y_true - preds$q0_975, 0)

preds$wis <- preds$mult*(preds$w0*preds$ae + preds$w.5*preds$IS.5 + preds$w.8*preds$IS.8 + preds$w.95*preds$IS.95)


## add in wis on standardized data
# preds$K <- 3
# preds$mult <- 1/(preds$K + .5)
# preds$w0 <- 0.5
preds$ae_std <- abs(preds$y_true - preds$q0_500)/preds$max_cases

# preds$w.5 <- .5/2
# preds$w.8 <- .2/2
# preds$w.95 <- .05/2

preds$IS.5_std <- (preds$q0_750 - preds$q0_250)/preds$max_cases + 
  (2/.5)*ifelse(preds$y_true < preds$q0_250, (preds$q0_250 - preds$y_true)/preds$max_cases, 0) + 
  (2/.5)*ifelse(preds$y_true > preds$q0_750, (preds$y_true - preds$q0_750)/preds$max_cases, 0)

preds$IS.8_std  <- (preds$q0_900 - preds$q0_100)/preds$max_cases + 
  (2/.2)*ifelse(preds$y_true < preds$q0_100, (preds$q0_100 - preds$y_true)/preds$max_cases, 0) + 
  (2/.2)*ifelse(preds$y_true > preds$q0_900, (preds$y_true - preds$q0_900)/preds$max_cases, 0)

preds$IS.95_std  <- (preds$q0_975 - preds$q0_025)/preds$max_cases + 
  (2/.05)*ifelse(preds$y_true < preds$q0_025, (preds$q0_025 - preds$y_true)/preds$max_cases, 0) + 
  (2/.05)*ifelse(preds$y_true > preds$q0_975, (preds$y_true - preds$q0_975)/preds$max_cases, 0)

preds$wis_std  <- preds$mult*(preds$w0*preds$ae_std  + preds$w.5*preds$IS.5_std  + preds$w.8*preds$IS.8_std  + preds$w.95*preds$IS.95_std )



## make coverage
preds$cov95 <- 0
preds[preds$q0_025 <= preds$y_true & preds$q0_975 >= preds$y_true,]$cov95 <- 1
preds$cov80 <- 0
preds[preds$q0_100 <= preds$y_true & preds$q0_900 >= preds$y_true,]$cov80 <- 1
preds$cov50 <- 0
preds[preds$q0_250 <= preds$y_true & preds$q0_750 >= preds$y_true,]$cov50 <- 1






## subset to the randomly partitioned data
preds_full <- preds
preds <- subset(preds_full, model %in% grep("_rp",unique(preds$model),invert=T,value=T))

## unique states and dates
unq_state_date <- preds_full[,list(n=1),by=c("series_id","last_obs_date","step_ahead")]





## --------------------------------------------------------------------------------------
## metrics

# overall
metrics_overall      <- preds[,list(mae = mean(abs(q0_500 - y_true)),
                                    wis = mean(wis)), by = c("model","train_size")]
metrics_overall$mae_m0 <- subset(metrics_overall, model == "M(0)")$mae
metrics_overall$rmae <- metrics_overall$mae / metrics_overall$mae_m0

metrics_overall$wis_m0 <- subset(metrics_overall, model == "M(0)")$wis
metrics_overall$rwis <- metrics_overall$wis / metrics_overall$wis_m0

((qplot(train_size, rmae, data = subset(metrics_overall, model %in% c("M(r,t)","M(st,t)","M(a,t)")), color = model, group = model, size = I(5)) + 
  scale_x_log10() + 
    scale_color_manual(values = MODEL_COLORS, drop = FALSE) +
  geom_hline(aes(yintercept = 1)))|
(qplot(train_size, rwis, data = subset(metrics_overall, model %in% c("M(r,t)","M(st,t)","M(a,t)")), color = model, group = model, size = I(5)) + 
  scale_x_log10() + 
  geom_hline(aes(yintercept = 1)))) /
((qplot(train_size, rmae, data = subset(metrics_overall, model %in% c("M(r,v)","M(sv,v)","M(a,v)")), color = model, group = model, size = I(5)) + 
    scale_x_log10() + 
    geom_hline(aes(yintercept = 1)))|
  (qplot(train_size, rwis, data = subset(metrics_overall, model %in% c("M(r,v)","M(sv,v)","M(a,v)")), color = model, group = model, size = I(5)) + 
     scale_x_log10() + 
     geom_hline(aes(yintercept = 1))))




# # --------------------------------------------------------------------------------------
## bootstrapping
bootstrap_rmse_like_summary_fast <- function(preds, Nb = 1L, Nwindow = 15L,
                                             baseline_model = "mod_baseline",
                                             progress = TRUE) {
  stopifnot(is.data.table(preds))
  stopifnot(all(c("series_id", "ref_date", "train_mod", "model_input", "q0_500", "y_true",
                  "max_cases", "wis", "wis_std", "cov95", "cov80", "cov50") %in% names(preds)))
  
  # Ensure fast joins
  if (!haskey(preds) || !identical(key(preds), c("series_id", "ref_date"))) {
    setkey(preds, series_id, ref_date)
  }
  
  unqstates <- sort(unique(preds$series_id))
  unqdates  <- sort(unique(preds$ref_date))
  
  nd <- length(unqdates)
  if (nd < Nwindow) stop("Nwindow is larger than number of unique ref_date values.")
  
  # Same as your code: sample start indices from 1:(nd - Nwindow)
  nstarts <- nd - Nwindow
  if (nstarts < 1L) stop("Need at least one possible window start: length(unqdates) must be > Nwindow.")
  
  # Same as your code: length(unqdates) / Nwindow (not necessarily integer)
  # sample() will coerce size to integer; we mimic that explicitly
  n_start_draws <- as.integer(nd / Nwindow)
  if (n_start_draws < 1L) n_start_draws <- 1L
  
  offsets <- 0:(Nwindow - 1L)
  
  out_list <- vector("list", Nb)
  
  for (b in seq_len(Nb)) {
    if (progress) message("bootstrap ", b, "/", Nb)
    
    # Resample states with replacement (same as your code)
    states_b <- sort(sample(unqstates, length(unqstates), replace = TRUE))
    
    # For EACH resampled state, resample window start indices with replacement
    # and then expand each start into Nwindow consecutive dates.
    #
    # Vectorized equivalent of:
    # for each j in states_b:
    #   temp_start_dates <- sample(1:nstarts, n_start_draws, replace=TRUE)
    #   for each k:
    #     my_ref_dates <- unqdates[start:(start+Nwindow-1)]
    #     rbind subset(preds, series_id==state & ref_date %in% my_ref_dates)
    #
    # IMPORTANT: We preserve multiplicity exactly by creating repeated key rows.
    state_rep <- rep(states_b, each = n_start_draws)
    starts    <- sample.int(nstarts, size = length(states_b) * n_start_draws, replace = TRUE)
    
    # Expand starts to date indices (preserves duplicates like the original rbind loop)
    date_idx <- rep(starts, each = Nwindow) + rep(offsets, times = length(starts))
    keys_dt <- data.table(
      series_id = rep(state_rep, each = Nwindow),
      ref_date  = unqdates[date_idx]
    )
    
    # Single join to pull all rows needed; allow.cartesian keeps duplicates
    temp_df <- preds[keys_dt, on = .(series_id, ref_date), nomatch = 0L, allow.cartesian = TRUE]
    
    # Summarize exactly as you did
    temp_sum_overall <- temp_df[, .(
      mae     = mean(abs(q0_500 - y_true)),
      mae_std = mean(abs(q0_500 - y_true) / max_cases),
      wis     = mean(wis),
      wis_std = mean(wis_std),
      cov95   = mean(cov95),
      cov80   = mean(cov80),
      cov50   = mean(cov50)
    ), by = c("train_mod", "model_input")]
    
    # Baseline model values (M(0))
    base <- temp_sum_overall[train_mod == baseline_model,
                             .(mae_M0 = mae, mae_std_M0 = mae_std, wis_M0 = wis, wis_std_M0 = wis_std)]
    
    # If baseline missing for some reason, fill with NA to avoid errors
    if (nrow(base) == 0L) {
      base <- data.table(mae_M0 = NA_real_, mae_std_M0 = NA_real_, wis_M0 = NA_real_, wis_std_M0 = NA_real_)
    }
    
    temp_sum_overall[, `:=`(
      mae_M0     = base$mae_M0[1],
      mae_std_M0 = base$mae_std_M0[1],
      wis_M0     = base$wis_M0[1],
      wis_std_M0 = base$wis_std_M0[1]
    )]
    
    temp_sum_overall[, `:=`(
      rmae     = mae / mae_M0,
      rmae_std = mae_std / mae_std_M0,
      rwis     = wis / wis_M0,
      rwis_std = wis_std / wis_std_M0,
      bs_id    = b
    )]
    
    out_list[[b]] <- temp_sum_overall
  }
  
  rbindlist(out_list, use.names = TRUE, fill = TRUE)
}



## get the bootstraps
Nb <- 5000
preds_bs <- bootstrap_rmse_like_summary_fast(preds = preds, 
                                             Nb = Nb, 
                                             Nwindow = 15)


## summarize
bs_uq_preds <- preds_bs[,list(rmae_q0_025 = quantile(rmae, probs=.025),
                              rmae_q0_975 = quantile(rmae, probs=.975),
                              rmae_std_q0_025 = quantile(rmae_std, probs=.025),
                              rmae_std_q0_975 = quantile(rmae_std, probs=.975),
                              rwis_q0_025 = quantile(rwis, probs = .025),
                              rwis_q0_975 = quantile(rwis, probs = .975),
                              rwis_std_q0_025 = quantile(rwis_std, probs = .025),
                              rwis_std_q0_975 = quantile(rwis_std, probs = .975),
                              cov95_q0_025 = quantile(cov95, probs = .025),
                              cov95_q0_975 = quantile(cov95, probs = .975),
                              cov80_q0_025 = quantile(cov80, probs = .025),
                              cov80_q0_975 = quantile(cov80, probs = .975),
                              cov50_q0_025 = quantile(cov50, probs = .025),
                              cov50_q0_975 = quantile(cov50, probs = .975)),by=c("train_mod","model_input")]







## ===============================================================================================================
## 4) MAE overall, by horizon, by state, by plos dates

# overall
metrics_overall      <- preds[,list(mae = mean(abs(q0_500 - y_true)),
                                wis = mean(wis)), by = c("train_mod","model","train_size","model_input")]
metrics_overall$mae_m0 <- subset(metrics_overall, train_mod == "mod_baseline")$mae
metrics_overall$rmae <- metrics_overall$mae / metrics_overall$mae_m0

metrics_overall$wis_m0 <- subset(metrics_overall, train_mod == "mod_baseline")$wis
metrics_overall$rwis <- metrics_overall$wis / metrics_overall$wis_m0

metrics_overall <- merge(metrics_overall,
                     subset(bs_uq_preds, select=c("train_mod","model_input","rmae_q0_025","rmae_q0_975","rwis_q0_025","rwis_q0_975")),
                     by=c("train_mod","model_input"))

pd <- position_dodge(width = 0.6)


## which 
card_sizes <- sort(unique(metrics_overall$train_size)[which(unique(metrics_overall$train_size) %% 100 == 0)])
all_sizes <- setdiff(unique(metrics_overall$train_size), card_sizes)

metrics_overall$train_size_fact <- as.character(metrics_overall$train_size)
metrics_overall[metrics_overall$train_size %in% all_sizes,]$train_size_fact <- "All"
metrics_overall$train_size_fact <- factor(as.factor(metrics_overall$train_size_fact), levels = c(as.character(card_sizes),"All"))

# metrics_overall <- subset(metrics_overall, train_size > 3000)
metrics_overall <- subset(metrics_overall, train_size <= 25000000)

# unq_train_size <- sort(unique(metrics_overall$train_size))
# unq_train_size <- unq_train_size[which(unq_train_size %% 100 == 0)]
# unq_train_size <- sort(unq_train_size[which(unq_train_size > 0)])
# metrics_overall <- subset(metrics_overall, train_size %in% unq_train_size)
# metrics_overall$train_size_fact <- factor(as.factor(metrics_overall$train_size), levels = sort(unq_train_size))

brks <- c(2.5e3, 2.5e4, 2.5e5, 1e6, 2.5e6, 5e6, 1e7)

labs <- scales::label_number(scale_cut = scales::cut_short_scale(), accuracy=.1)(brks)
names(labs) <- as.character(brks)
labs


setEPS()
postscript(paste0(figpath,"training_subsets_metrics.eps"), width = 10, height = 4)
print((ggplot(subset(metrics_overall, model != "M(0)" & model != "M(st,v)" & model != "M(sv,t)"), aes(x=train_size_fact, y = rmae, color = model)) +
   geom_errorbar(aes(ymin = rmae_q0_025, ymax = rmae_q0_975),
                position = pd, width = 0.5) +
   geom_point(aes(fill = model), position = pd, size = 2, shape=I(21), color = I('black'))+
   scale_fill_manual(values = MODEL_COLORS, drop = FALSE, name="") +
   scale_color_manual(values = MODEL_COLORS, drop = FALSE, name="") +
   geom_hline(aes(yintercept = 1), linetype=I(2))+
    theme(plot.title = element_text(hjust = 0.5),
          legend.position="bottom",
          axis.title.x = element_text(margin = margin(t=10)),
          axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 10))+ 
   ylab("")+
    ggtitle("rMAE")+
    scale_y_continuous(breaks = seq(0,100,.2))+
   scale_x_discrete(labels = labs)+
    guides(colour = guide_legend(nrow = 1),
           fill   = guide_legend(nrow = 1))+
    xlab("Number of Training Data Examples")) |
(ggplot(subset(metrics_overall, model != "M(0)" & model != "M(st,v)" & model != "M(sv,t)"), aes(x=train_size_fact, y = rwis, color = model)) +
   geom_errorbar(aes(ymin = rwis_q0_025, ymax = rwis_q0_975),
                position = pd, width = 0.5) +
   geom_point(aes(fill = model), position = pd, size = 2, shape=I(21), color = I('black'))+
   scale_fill_manual(values = MODEL_COLORS, drop = FALSE, name="") +
   scale_color_manual(values = MODEL_COLORS, drop = FALSE, name="") +
   geom_hline(aes(yintercept = 1), linetype=I(2))+
   ggtitle("rWIS")+
   scale_y_continuous(breaks = seq(0,100,.2))+
   theme(plot.title = element_text(hjust = 0.5),
         legend.position="bottom",
         axis.title.x = element_text(margin = margin(t=10)),
         axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 10))+ 
   ylab("")+
   guides(colour = guide_legend(nrow = 1),
          fill   = guide_legend(nrow = 1))+
   scale_x_discrete(labels = labs)+
   xlab("Number of Training Data Examples")))
dev.off()






## --------------------------------------------------------------------------------
## --------------------------------------------------------------------------------

# keep model columns exactly as-is
preds_bs <- subset(preds_bs, train_mod != "mod_baseline")

preds_bs$train_size <- as.numeric(sub(".*_(\\d+)$", "\\1", preds_bs$train_mod))
preds_bs$model <- paste0(sub("_[0-9]+$", "", preds_bs$train_mod),"_",preds_bs$model_input)

## which 
card_sizes <- sort(unique(preds_bs$train_size)[which(unique(preds_bs$train_size) %% 100 == 0)])
all_sizes <- setdiff(unique(preds_bs$train_size), card_sizes)

preds_bs$train_size_fact <- as.character(preds_bs$train_size)
preds_bs[preds_bs$train_size %in% all_sizes,]$train_size_fact <- "All"
preds_bs$train_size_fact <- factor(as.factor(preds_bs$train_size_fact), levels = c(as.character(card_sizes),"All"))
# preds_bs <- subset(preds_bs, train_size %in% unq_train_size)


# rename models
preds_bs[preds_bs$model == "mod_real_tc"]$model <- "M(r,t)"
preds_bs[preds_bs$model == "mod_real_vac"]$model <- "M(r,v)"
preds_bs[preds_bs$model == "mod_syn_tc_tc"]$model <- "M(st,t)"
preds_bs[preds_bs$model == "mod_syn_tc_vac"]$model <- "M(st,v)"
preds_bs[preds_bs$model == "mod_syn_vac_tc"]$model <- "M(sv,t)"
preds_bs[preds_bs$model == "mod_syn_vac_vac"]$model <- "M(sv,v)"
preds_bs[preds_bs$model == "mod_all_tc"]$model <- "M(a,t)"
preds_bs[preds_bs$model == "mod_all_vac"]$model <- "M(a,v)"

model_order <- c("M(r,t)","M(st,t)","M(sv,t)","M(a,t)","M(r,v)","M(st,v)","M(sv,v)","M(a,v)")
preds_bs$model <- factor(as.factor(preds_bs$model), levels = model_order)


tab_rmae <- data.frame(reshape2::dcast(preds_bs, bs_id + train_size_fact ~ model, value.var="rmae"), check.names = FALSE)
tab_rwis <- data.frame(reshape2::dcast(preds_bs, bs_id + train_size_fact ~ model, value.var="rwis"), check.names = FALSE)



# =========================================================
## Qtrain
mae_qtrain <- data.table(
  type = "rMAE",
  train_size_fact = tab_rmae$train_size_fact,
  `M(r,t) vs M(st,t)` = tab_rmae$`M(r,t)` - tab_rmae$`M(st,t)`,
  `M(r,t) vs M(a,t)` = tab_rmae$`M(r,t)` - tab_rmae$`M(a,t)`,
  `M(st,t) vs M(a,t)` = tab_rmae$`M(st,t)` - tab_rmae$`M(a,t)`,
  ##
  `M(r,v) vs M(sv,v)` = tab_rmae$`M(r,v)` - tab_rmae$`M(sv,v)`,
  `M(r,v) vs M(a,v)` = tab_rmae$`M(r,v)` - tab_rmae$`M(a,v)`,
  `M(sv,v) vs M(a,v)` = tab_rmae$`M(sv,v)` - tab_rmae$`M(a,v)`
)

wis_qtrain <- data.table(
  type = "rWIS",
  train_size_fact = tab_rwis$train_size_fact,
  `M(r,t) vs M(st,t)` = tab_rwis$`M(r,t)` - tab_rwis$`M(st,t)`,
  `M(r,t) vs M(a,t)` = tab_rwis$`M(r,t)` - tab_rwis$`M(a,t)`,
  `M(st,t) vs M(a,t)` = tab_rwis$`M(st,t)` - tab_rwis$`M(a,t)`,
  ##
  `M(r,v) vs M(sv,v)` = tab_rwis$`M(r,v)` - tab_rwis$`M(sv,v)`,
  `M(r,v) vs M(a,v)` = tab_rwis$`M(r,v)` - tab_rwis$`M(a,v)`,
  `M(sv,v) vs M(a,v)` = tab_rwis$`M(sv,v)` - tab_rwis$`M(a,v)`
)

q1df <- rbind(mae_qtrain, wis_qtrain)
q1dfmelt <- reshape2::melt(q1df, id.vars=c("type","train_size_fact"))

q1dfmelt <- subset(q1dfmelt, train_size_fact != "2500")
q1dfmelt <- subset(q1dfmelt, !is.na(value))

q1dfmelt_summary <- data.table(q1dfmelt)[,list(
  q.025 = quantile(value, probs = .025, na.rm=T),
  avg = mean(value, na.rm=T),
  q.975 = quantile(value, probs = .975, na.rm=T)
), by=c("type","train_size_fact","variable")]


q1dfcolors <- data.frame(
  variable = c(
    "M(r,t) vs M(st,t)",
    "M(r,t) vs M(a,t)",
    "M(st,t) vs M(a,t)",
    "M(r,v) vs M(sv,v)",
    "M(r,v) vs M(a,v)",
    "M(sv,v) vs M(a,v)"
  ),
  color = c("#D1B4BE","#E3A6B8", "#B8C0DD","#8B5D24","#A72C5B","#4E6E63")
)

q1dfmelt <- merge(q1dfmelt, q1dfcolors, by="variable", all.x=TRUE)

q1dfmelt$train_size_fact <- factor(as.factor(q1dfmelt$train_size_fact), levels = c(as.character(card_sizes),"All"))

brks <- c(2.5e3, 2.5e4, 2.5e5, 1e6, 2.5e6, 5e6, 1e7)

fmt_km <- function(x) {
  lab <- scales::label_number(
    scale_cut = scales::cut_short_scale(),
    accuracy  = 1
  )(x)
  
  lab <- sub("\\.0(?=[A-Za-z]|$)", "", lab, perl = TRUE)  # drop trailing .0
  lab <- sub("K$", "k", lab)                              # K -> k (keep M)
  lab
}

labs <- fmt_km(brks)
names(labs) <- as.character(brks)
labs


setEPS()
postscript(paste0(figpath,"training_subsets_paired_differences_boxplot.eps"), width = 16, height = 7)
(qplot(train_size_fact, value, data = subset(q1dfmelt, type %in% c("rMAE","rWIS")),
       geom="boxplot", fill = variable, show.legend = FALSE) +
    facet_grid(type~variable, scales="free") +
    scale_fill_manual(values = q1dfcolors$color, drop = FALSE, name="") +
    geom_hline(aes(yintercept = 0), linetype=I(2)) +
    xlab("") + ylab("Model A - Model B") + 
    ggtitle("Paired Differences") +
    xlab("Number of Training Data Examples")+
    scale_x_discrete(labels = labs)+
    theme(plot.title = element_text(hjust = 0.5),
          axis.title.x = element_text(margin = margin(t=10)),
          axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 10))) 
dev.off()


setEPS()
postscript(paste0(figpath,"training_subsets_paired_differences.eps"), width = 12, height = 4)
ggplot(aes(x=train_size_fact, y=avg,color = variable),data = q1dfmelt_summary, show.legend = FALSE)+
  geom_hline(aes(yintercept = 0), linetype=I(2))+
  geom_errorbar(aes(ymin = q.025, ymax = q.975),size=I(1.15),
                position = pd, width = .5, show.legend = FALSE) +
  facet_grid(type~variable, scales="free") +
  geom_point(aes(fill = variable), position = pd, size = 4, shape=I(21), color = I('black'), show.legend = FALSE)+
  scale_fill_manual(values = q1dfcolors$color, drop = FALSE, name="") +
  scale_color_manual(values = q1dfcolors$color, drop = FALSE, name="") +
  xlab("") + ylab("Model A - Model B") + 
  ggtitle("Paired Differences") +
  xlab("Number of Training Data Examples")+
  scale_x_discrete(labels = labs)+
  theme(plot.title = element_text(hjust = 0.5),
        axis.title.x = element_text(margin = margin(t=10), size = 12),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 9)) 
dev.off()






