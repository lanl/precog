
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
fcstpath <- paste0(here::here("synthetic_and_genetic_forecasting", "output"), "/")  #"/Users/dosthus/Documents/DR24_PreCog/projects/four_quadrants/quadrants_complete_2-11-26/output/"
figpath <- paste0(here::here("synthetic_and_genetic_forecasting", "figs"), "/")  #"/Users/dosthus/Documents/DR24_PreCog/projects/four_quadrants/quadrants_complete_2-11-26/figs/"
datapath <- paste0(here::here("synthetic_and_genetic_forecasting", "data_clean"), "/") # "/Users/dosthus/Documents/DR24_PreCog/projects/four_quadrants/quadrants_complete_2-11-26/data_clean/"
modelpath <- paste0(here::here("synthetic_and_genetic_forecasting", "trained_models"), "/")  #"/Users/dosthus/Documents/DR24_PreCog/projects/four_quadrants/quadrants_complete_2-11-26/trained_models/"

## Read in test data
test_covid <- readRDS(paste0(datapath,"test_covid.RDS"))

## read in forecasts
preds <- fread(file = paste0(fcstpath,"forecasts.csv"))

# Compute the max cases per series_id
max_cases_df <- data.table(rbindlist(lapply(test_covid, function(x){data.frame(series_id = x$series_id[1], max_cases = max(x$cases))})))

## append the max cases to each state in preds
preds <- merge(preds, 
               max_cases_df,
               by = "series_id",
               all.x=T)


# --------------------------------------------------------
# Read in config files
# --------------------------------------------------------
lf_rds <- list.files(modelpath, pattern = ".rds")
df_cfg <- NULL
for(i in 1:length(lf_rds)){
  temp_cfg <- readRDS(paste0(modelpath,lf_rds[i]))
  df_cfg <- rbind(df_cfg,
                  data.frame(name = lf_rds[i],
                             run_time = temp_cfg$run_time,
                             total_possible_flashcards = temp_cfg$total_possible_flashcards,
                             cards_seen = temp_cfg$cards_seen,
                             cards_seen_no_noise_injection = temp_cfg$cards_seen/2,
                             weight_updates = temp_cfg$weight_updates))
  
}
df_cfg <- df_cfg[order(df_cfg$total_possible_flashcards),]
df_cfg$n_unique_cards_seen <- df_cfg$cards_seen_no_noise_injection / df_cfg$total_possible_flashcards
df_cfg

# --------------------------------------------------------
# PLOT AN EXAMPLE OF THE DATA
# --------------------------------------------------------

df_covid <- rbind(test_covid[[1]],
                  test_covid[[5]])
df_covid <- subset(df_covid, date <= max(preds$ref_date))
df_covid_melt <- data.table(reshape2::melt(df_covid, id.vars = c("series_id","time","date","cases")))
df_covid_melt$state <- tools::toTitleCase(gsub("unitedstates_","",df_covid_melt$series_id))
df_covid_melt$variant <- tools::toTitleCase(gsub("var_attr_","",df_covid_melt$variable))
df_covid_melt[df_covid_melt$variant == "Ba1",]$variant <- "Omicron (BA.1)"
df_covid_melt[df_covid_melt$variant == "Ba2",]$variant <- "BA.2"
df_covid_melt[df_covid_melt$variant == "Ba45",]$variant <- "BA.4/BA.5"
df_covid_melt[df_covid_melt$variant == "Xbb",]$variant <- "XBB"

var_order <- df_covid_melt[,list(date = date[min(which(value > 0))]), by="variant"]
var_order <- var_order[order(var_order$date),]
df_covid_melt$variant <- factor(as.factor(df_covid_melt$variant), levels = var_order$variant)

p1 <- qplot(date, cases, data = df_covid_melt, geom="line")+
  facet_wrap(~state, scale = "free", ncol = 1)+
  scale_y_continuous(
    labels = function(x) {
      ifelse(abs(x) >= 1000,
             paste0(number(x / 1000, accuracy = 1), "k"),
             number(x))
    }
  )+
  ylab("")+
  xlab("")+
  ggtitle("Total Cases")+
  theme(plot.title = element_text(hjust = 0.5))

p2 <- ggplot(data = df_covid_melt)+
       geom_line(aes(x = date, y = value/cases, group = variant, color = variant))+
       facet_wrap(~state, scale = "free", ncol = 1)+
       theme(legend.position = "bottom")+
       scale_color_manual(values = qualitative_hcl(12, palette = "Dark 3"))+
       guides(color = guide_legend(title = NULL, override.aes = list(linewidth = 2), nrow = 1, byrow = TRUE))+
       ylab("")+
       xlab("")+
       ggtitle("Variant Proportions")+
       theme(plot.title = element_text(hjust = 0.5))

p3 <- ggplot(data = df_covid_melt)+
        geom_line(aes(x = date, y = value, group = variant, color = variant), show.legend = F)+
        facet_wrap(~state, scale = "free", ncol = 1)+
        scale_y_sqrt(
          labels = function(x) {
            ifelse(abs(x) >= 1000,
                   paste0(number(x / 1000, accuracy = 1), "k"),
                   number(x))
          }
        )+
        theme(legend.position = "right")+
        scale_color_manual(values = qualitative_hcl(12, palette = "Dark 3"))+
        guides(color = guide_legend(override.aes = list(linewidth = 2), ncol = 1, byrow = TRUE))+
        ggtitle("Variant Attributable Cases")+
        ylab("")+
        xlab("")+
        theme(plot.title = element_text(hjust = 0.5))

setEPS() 
postscript(paste0(figpath,"example_data.eps"),width = 11, height = 5)
  (p1 | p2 | p3) +
  plot_annotation(tag_levels = "a", tag_prefix = "(", tag_suffix = ")")
dev.off()




# --------------------------------------------------------
# Back to preds
# --------------------------------------------------------

## dates that align with the PLOS paper
min_plos_date <- ymd("2020-07-28")
max_plos_date <- ymd("2021-12-21")

## define model colors
MODEL_COLORS <- c(
  "M(0)" = "darkgrey",
  "M(r,t)" = "#fb9a99",  
  "M(st,t)" = "#a6cee3", 
  "M(sv,t)" = "#b2df8a", 
  "M(a,t)" = "#cab2d6",  
  "M(r,v)" = "#e31a1c",  
  "M(st,v)" = "#1f78b4",  
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
  
  # old scripts sometimes used `model_name` (and sometimes `model_code`)
  model   = paste0(as.character(train_mod), "_", as.character(model_input)),
  
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





# # --------------------------------------------------------------------------------------
## bootstrapping
bootstrap_rmse_like_summary_fast <- function(preds, Nb = 1L, Nwindow = 15L,
                                             baseline_model = "M(0)",
                                             progress = TRUE) {
  stopifnot(is.data.table(preds))
  stopifnot(all(c("series_id", "ref_date", "model", "q0_500", "y_true",
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
    ), by = model]
    
    # Baseline model values (M(0))
    base <- temp_sum_overall[model == baseline_model,
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
preds_bs <- bootstrap_rmse_like_summary_fast(preds, Nb = Nb, Nwindow = 15)


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
                              cov50_q0_975 = quantile(cov50, probs = .975)),by=c("model")]



# --------------------------------------------------------------------
## get bootstraps assuming iid (the wrong way to do this)
## this is for a comparison in the Supplemental Materials

df_iid_mae <- data.frame(
  reshape2::dcast(preds, ref_date + series_id + horizon ~ model, value.var = "ae"),
  check.names = FALSE   # <-- critical: keep names like "M(0)" intact
)
df_iid_mae <- df_iid_mae[order(df_iid_mae$ref_date, df_iid_mae$series_id, df_iid_mae$horizon),]

df_iid_wis <- data.frame(
  reshape2::dcast(preds, ref_date + series_id + horizon ~ model, value.var = "wis"),
  check.names = FALSE   # <-- critical
)
df_iid_wis <- df_iid_wis[order(df_iid_wis$ref_date, df_iid_wis$series_id, df_iid_wis$horizon),]

preds_bs_iid <- NULL

# model columns: anything that starts with M(
model_cols <- grep("^M\\(", names(df_iid_mae), value = TRUE)

for(i in 1:Nb){
  print(i)
  
  ## get iid indices
  temp_idx <- sort(sample(seq_len(nrow(df_iid_mae)), nrow(df_iid_mae), replace = TRUE))
  
  ## subset data
  temp_mae <- df_iid_mae[temp_idx, , drop = FALSE]
  temp_wis <- df_iid_wis[temp_idx, , drop = FALSE]
  
  ## summarise (avoid subset(select=...) to handle non-syntactic names safely)
  temp_summary <- data.table(
    model = model_cols,
    mae   = as.numeric(colMeans(temp_mae[, model_cols, drop = FALSE], na.rm = TRUE)),
    wis   = as.numeric(colMeans(temp_wis[, model_cols, drop = FALSE], na.rm = TRUE)),
    bs_id = i
  )
  
  # baseline model: M(0)
  base_mae <- temp_summary[model == "M(0)", mae]
  base_wis <- temp_summary[model == "M(0)", wis]
  
  temp_summary[, rmae := mae / base_mae]
  temp_summary[, rwis := wis / base_wis]
  
  ## concatenate
  preds_bs_iid <- rbind(preds_bs_iid, temp_summary)
}








# --------------------------------------------------------------------
# plot the blocked bootstrap way vs. the iid boot strapped way
bs_comp <- rbind(data.table(type = "block", subset(preds_bs, select = c("model","rmae","rwis"))),
                 data.table(type = "iid", subset(preds_bs_iid, select = c("model","rmae","rwis"))))
bs_comp$type <- factor(as.factor(bs_comp$type), levels = c("iid","block"))
bs_comp$group <- paste0(bs_comp$model, "_", bs_comp$type)
bs_comp <- subset(bs_comp,model != "M(0)")


## get the overall mean
df_overall      <- preds[,list(mae = mean(abs(q0_500 - y_true)),
                               wis = mean(wis)), by = "model"]
df_overall$mae_m0 <- subset(df_overall, model == "M(0)")$mae
df_overall$rmae <- df_overall$mae / df_overall$mae_m0
df_overall$wis_m0 <- subset(df_overall, model == "M(0)")$wis
df_overall$rwis <- df_overall$wis / df_overall$wis_m0
df_overall <- subset(df_overall, model != "M(0)")

## plot it
setEPS() 
postscript(paste0(figpath,"bootstrap_comparison.eps"),width = 11, height = 5)
  (ggplot(bs_comp, aes(x = type, y = rmae, fill = model)) +
    geom_boxplot(width = 0.8, show.legend=F, color = I("black")) +
    facet_grid(. ~ model, switch = "x", space = "free_x") +
    labs(x = NULL, fill = NULL) +
    theme(
      strip.placement = "outside",
      strip.background = element_blank()
    )+
    ylab("rMAE")+
    geom_hline(aes(yintercept = 1), linetype = I(2))+
    geom_hline(aes(yintercept = rmae), data = df_overall, linetype = I(3), size=I(.5), color =I("white"))+
    scale_fill_manual(values = MODEL_COLORS, drop = FALSE)) /
    
  (ggplot(bs_comp, aes(x = type, y = rwis, fill = model)) +
    geom_boxplot(width = 0.8, show.legend=F) +
    facet_grid(. ~ model, switch = "x", space = "free_x") +
    labs(x = NULL, fill = NULL) +
    theme(
      strip.placement = "outside",
      strip.background = element_blank()
    )+
    ylab("rWIS")+
    geom_hline(aes(yintercept = 1), linetype = I(2))+
    geom_hline(aes(yintercept = rwis), data = df_overall, linetype = I(3), size=I(.5), color =I("white"))+
    scale_fill_manual(values = MODEL_COLORS, drop = FALSE))
dev.off()




## ===============================================================================================================
## 4) MAE overall, by horizon, by state, by plos dates

# overall
mae_overall      <- preds[,list(mae = mean(abs(q0_500 - y_true)),
                                mae_std = mean(abs(q0_500 - y_true)/max_cases)), by = "model"]
mae_overall$mae_m0 <- subset(mae_overall, model == "M(0)")$mae
mae_overall$mae_std_m0 <- subset(mae_overall, model == "M(0)")$mae_std
mae_overall$rmae <- mae_overall$mae / mae_overall$mae_m0
mae_overall$rmae_std <- mae_overall$mae_std / mae_overall$mae_std_m0
mae_overall <- merge(mae_overall,
                     subset(bs_uq_preds, select=c("model","rmae_q0_025","rmae_q0_975","rmae_std_q0_025","rmae_std_q0_975")),
                     by="model")

# --------------------------------------------
# plot
setEPS() 
postscript(paste0(figpath,"mae_overall.eps"),width = 7, height = 4)
ggplot(mae_overall, aes(x = model, y = rmae, fill = model)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.6, alpha = 0.8) +
  geom_col(width = 0.65, color = "black", alpha = 0.95) +
  scale_fill_manual(values = MODEL_COLORS, drop = FALSE) +
  scale_y_continuous("Relative MAE", labels = scales::number_format(accuracy = 0.01),
                     breaks = seq(0, 10, .1)) +
  labs(x = "", title = "Relative MAE by Model") +
  geom_text(aes(y = 0.8 * rmae, label = round(rmae, 3))) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    axis.text.x = element_text(angle = 0, hjust = .5, size = 12)
  )
dev.off()

setEPS() 
postscript(paste0(figpath,"mae_std_overall.eps"),width = 7, height = 4)
ggplot(mae_overall, aes(x = model, y = rmae_std, fill = model)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.6, alpha = 0.8) +
  geom_col(width = 0.65, color = "black", alpha = 0.95) +
  scale_fill_manual(values = MODEL_COLORS, drop = FALSE) +
  scale_y_continuous("Relative MAE", labels = scales::number_format(accuracy = 0.01),
                     breaks = seq(0, 10, .1)) +
  labs(x = "", title = "Relative MAE for Standardized Cases by Model") +
  geom_text(aes(y = 0.8 * rmae_std, label = round(rmae_std, 3))) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    axis.text.x = element_text(angle = 0, hjust = .5, size = 12)
  )
dev.off()

# --------------------------------------------
# plot
setEPS() 
postscript(paste0(figpath,"mae_overall_errorbar.eps"),width = 7, height = 4)
ggplot(mae_overall, aes(x = model, y = rmae, fill = model)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.6) +
  geom_col(width = 0.65, color = "black") +
  geom_errorbar(aes(ymin = rmae_q0_025, ymax = rmae_q0_975),
                width = 0.15, linewidth = 0.6, color = "black") +
  scale_fill_manual(values = MODEL_COLORS, drop = FALSE) +
  scale_y_continuous("Relative MAE", labels = scales::number_format(accuracy = 0.01),
                     breaks = seq(0, 10, .1)) +
  labs(x = "", title = "Relative MAE by Model") +
  geom_text(aes(y = 0.8 * rmae, label = round(rmae, 3))) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    axis.text.x = element_text(angle = 0, hjust = .5, size = 12)
  )
dev.off()


setEPS() 
postscript(paste0(figpath,"mae_std_overall_errorbar.eps"),width = 7, height = 4)
ggplot(mae_overall, aes(x = model, y = rmae_std, fill = model)) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.6, alpha = 0.8) +
  geom_col(width = 0.65, color = "black", alpha = 0.95) +
  geom_errorbar(aes(ymin = rmae_std_q0_025, ymax = rmae_std_q0_975),
                width = 0.15, linewidth = 0.6, color = "black") +
  scale_fill_manual(values = MODEL_COLORS, drop = FALSE) +
  scale_y_continuous("Relative MAE", labels = scales::number_format(accuracy = 0.01),
                     breaks = seq(0, 10, .1)) +
  labs(x = "", title = "Relative MAE for Standardized Cases by Model") +
  geom_text(aes(y = 0.8 * rmae_std, label = round(rmae_std, 3))) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    axis.text.x = element_text(angle = 0, hjust = .5, size = 12)
  )
dev.off()


## subset to comparison dates
mae_overall_comp      <- preds[plos_date == 1,list(mae = mean(abs(q0_500 - y_true)),
                                                   mae_std = mean(abs(q0_500 - y_true)/max_cases)), by = "model"]
mae_overall_comp$mae_m0 <- subset(mae_overall_comp , model == "M(0)")$mae
mae_overall_comp$mae_std_m0 <- subset(mae_overall_comp , model == "M(0)")$mae_std
mae_overall_comp$rmae <- mae_overall_comp$mae / mae_overall_comp$mae_m0
mae_overall_comp$rmae_std <- mae_overall_comp$mae_std / mae_overall_comp$mae_std_m0
mae_overall_comp


# --------------------------------------------
# plot
setEPS() 
postscript(paste0(figpath,"mae_overall_comp.eps"),width = 7, height = 4)
ggplot(mae_overall_comp, aes(x = model, y = rmae, fill = model)) +
  geom_col(width = 0.65, color = "black", alpha = 0.95) +
  scale_fill_manual(values = MODEL_COLORS, drop = FALSE) +
  scale_y_continuous("Relative MAE", labels = number_format(accuracy = 0.01), breaks = seq(0,10,.1)) +
  labs(x = "", title = paste0("Relative MAE by Model: ", min_plos_date," - ",max_plos_date)) +
  theme_bw(base_size = 12) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.6, alpha = 0.8) +
  geom_text(aes(x = model, y = .8*rmae, label = round(rmae,3)))+
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    axis.text.x = element_text(angle = 0, hjust = .5, size = 12))
dev.off()


setEPS() 
postscript(paste0(figpath,"mae_std_overall_comp.eps"),width = 7, height = 4)
ggplot(mae_overall_comp, aes(x = model, y = rmae_std, fill = model)) +
  geom_col(width = 0.65, color = "black", alpha = 0.95) +
  scale_fill_manual(values = MODEL_COLORS, drop = FALSE) +
  scale_y_continuous("Relative MAE", labels = number_format(accuracy = 0.01), breaks = seq(0,10,.1)) +
  labs(x = "", title = paste0("Relative Standardized MAE by Model: ", min_plos_date," - ",max_plos_date)) +
  theme_bw(base_size = 12) +
  geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.6, alpha = 0.8) +
  geom_text(aes(x = model, y = .8*rmae_std, label = round(rmae_std,3)))+
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    axis.text.x = element_text(angle = 0, hjust = .5, size = 12))
dev.off()




# by horizon
mae_horizon      <- preds[,list(mae = mean(abs(q0_500 - y_true)),
                                mae_std = mean(abs(q0_500 - y_true)/max_cases)), by = c("model","step_ahead")]
mae_horizon_m0 <- subset(mae_horizon, model == "M(0)")
names(mae_horizon_m0) <- c("model","step_ahead","mae_m0","mae_std_m0")
mae_horizon <- merge(mae_horizon, subset(mae_horizon_m0, select = c("step_ahead","mae_m0","mae_std_m0")), by=c("step_ahead"), all.x=T)
mae_horizon$rmae <- mae_horizon$mae / mae_horizon$mae_m0
mae_horizon$rmae_std <- mae_horizon$mae_std / mae_horizon$mae_std_m0
mae_horizon

# --------------------------------------------
# plot
setEPS() 
postscript(paste0(figpath,"mae_horizon.eps"),width = 7, height = 4)
ggplot(mae_horizon, aes(x = step_ahead, y = rmae, color = model, group = model)) +
  geom_line(size=I(2)) +
  scale_color_manual(values = MODEL_COLORS, drop = FALSE, name="") +
  scale_fill_manual(values = MODEL_COLORS, drop = FALSE, name="") +
  scale_y_continuous("Relative MAE", labels = number_format(accuracy = 0.01), breaks = seq(0,10,.05)) +
  labs(x = "Horizon", title = "Relative MAE by Model") +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    axis.text.x = element_text(angle = 0, hjust = .5, size = 12))
dev.off()


setEPS() 
postscript(paste0(figpath,"mae_std_horizon.eps"),width = 7, height = 4)
ggplot(mae_horizon, aes(x = step_ahead, y = rmae_std, color = model, group = model)) +
  geom_line(size=I(2)) +
  scale_color_manual(values = MODEL_COLORS, drop = FALSE, name="") +
  scale_fill_manual(values = MODEL_COLORS, drop = FALSE, name="") +
  scale_y_continuous("Relative MAE", labels = number_format(accuracy = 0.01), breaks = seq(0,10,.05)) +
  labs(x = "Horizon", title = "Relative Standardized MAE by Model") +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    axis.text.x = element_text(angle = 0, hjust = .5, size = 12))
dev.off()



mae_horizon_comp      <- preds[plos_date == 1,list(mae = mean(abs(q0_500 - y_true)),
                                                   mae_std = mean(abs(q0_500 - y_true)/max_cases)), by = c("model","step_ahead")]
mae_horizon_m0_comp <- subset(mae_horizon_comp, model == "M(0)")
names(mae_horizon_m0_comp) <- c("model","step_ahead","mae_m0","mae_std_m0")
mae_horizon_comp <- merge(mae_horizon_comp, subset(mae_horizon_m0_comp, select = c("step_ahead","mae_m0","mae_std_m0")), by=c("step_ahead"), all.x=T)
mae_horizon_comp$rmae <- mae_horizon_comp$mae / mae_horizon_comp$mae_m0
mae_horizon_comp$rmae_std <- mae_horizon_comp$mae_std / mae_horizon_comp$mae_std_m0
mae_horizon_comp

# --------------------------------------------
# plot
setEPS() 
postscript(paste0(figpath,"mae_horizon_comp.eps"),width = 7, height = 4)
ggplot(mae_horizon_comp, aes(x = step_ahead, y = rmae, color = model, group = model)) +
  geom_line(size=I(2)) +
  scale_color_manual(values = MODEL_COLORS, drop = FALSE, name="") +
  scale_fill_manual(values = MODEL_COLORS, drop = FALSE, name="") +  scale_y_continuous("Relative MAE", labels = number_format(accuracy = 0.01), breaks = seq(0,10,.05)) +
  labs(x = "Horizon", title = paste0("Relative MAE by Model: ", min_plos_date," - ",max_plos_date)) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    axis.text.x = element_text(angle = 0, hjust = .5, size = 12))
dev.off()


setEPS() 
postscript(paste0(figpath,"mae_std_horizon_comp.eps"),width = 7, height = 4)
ggplot(mae_horizon_comp, aes(x = step_ahead, y = rmae_std, color = model, group = model)) +
  geom_line(size=I(2)) +
  scale_color_manual(values = MODEL_COLORS, drop = FALSE, name="") +
  scale_fill_manual(values = MODEL_COLORS, drop = FALSE, name="") +  scale_y_continuous("Relative MAE", labels = number_format(accuracy = 0.01), breaks = seq(0,10,.05)) +
  labs(x = "Horizon", title = paste0("Relative MAE by Model: ", min_plos_date," - ",max_plos_date)) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    axis.text.x = element_text(angle = 0, hjust = .5, size = 12))
dev.off()


# by state
mae_state      <- preds[,list(mae = mean(abs(q0_500 - y_true))), by = c("model","series_id")]
mae_state_m0 <- subset(mae_state, model == "M(0)")
names(mae_state_m0) <- c("model","series_id","mae_m0")
mae_state <- merge(mae_state, subset(mae_state_m0, select = c("series_id","mae_m0")), by=c("series_id"), all.x=T)
mae_state$rmae <- mae_state$mae / mae_state$mae_m0
mae_state$state <- gsub("unitedstates_","",mae_state$series_id)
mae_state_order <- mae_state[,list(avg_mae = mean(rmae)),by="state"]
mae_state_order <- mae_state_order[order(mae_state_order$avg_mae),]
mae_state$state <- factor(as.factor(mae_state$state), levels = mae_state_order$state)

mae_state_rank <- data.table::copy(mae_state)
mae_state_rank[, model := factor(model, levels = names(MODEL_COLORS))]
setorder(mae_state_rank, state, rmae, model)
mae_state_rank[, rank := seq_len(.N), by = state]   

# --------------------------------------------
# plot
setEPS() 
postscript(paste0(figpath,"mae_state.eps"),width = 10, height = 5)
ggplot(mae_state, aes(x = state, y = rmae, color = model, group = model)) +
  geom_line(size=I(1)) +
  scale_color_manual(values = MODEL_COLORS, drop = FALSE, name="") +
  scale_y_continuous("Relative MAE", labels = number_format(accuracy = 0.01), breaks = seq(0,10,.05)) +
  labs(x = "", title = "Relative MAE by State") +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    axis.text.x = element_text(angle = 90, vjust = .5, hjust=1))
dev.off()



# --------------------------------------------
# plot
setEPS() 
postscript(paste0(figpath,"mae_rank_state.eps"),width = 10, height = 5)
ggplot(mae_state_rank, aes(x = state, y = rank, fill = model)) +
  geom_tile(color = "white", linewidth = 0.3, height = 0.9) +
  scale_fill_manual(values = MODEL_COLORS, drop = FALSE, name = "") +
  scale_y_continuous(limits = c(0, length(MODEL_COLORS)+1), breaks = 1:length(MODEL_COLORS), minor_breaks = NULL, expand = c(0,0)) +
  labs(
    x = "",
    y = "Rank (1 = best, 9 = worst)",
    title = "Relative MAE Rank by State"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))
dev.off()



mae_state_comp      <- preds[plos_date == 1,list(mae = mean(abs(q0_500 - y_true))), by = c("model","series_id")]
mae_state_m0_comp <- subset(mae_state_comp, model == "M(0)")
names(mae_state_m0_comp) <- c("model","series_id","mae_m0")
mae_state_comp <- merge(mae_state_comp, subset(mae_state_m0_comp, select = c("series_id","mae_m0")), by=c("series_id"), all.x=T)
mae_state_comp$rmae <- mae_state_comp$mae / mae_state_comp$mae_m0
mae_state_comp$state <- gsub("unitedstates_","",mae_state_comp$series_id)
mae_state_order_comp <- mae_state_comp[,list(avg_mae = mean(rmae)),by="state"]
mae_state_order_comp <- mae_state_order_comp[order(mae_state_order_comp$avg_mae),]
mae_state_comp$state <- factor(as.factor(mae_state_comp$state), levels = mae_state_order_comp$state)

mae_state_rank_comp <- data.table::copy(mae_state_comp)
mae_state_rank_comp[, model := factor(model, levels = names(MODEL_COLORS))]
setorder(mae_state_rank_comp, state, rmae, model)
mae_state_rank_comp[, rank := seq_len(.N), by = state]   # 1..7 per state

# --------------------------------------------
# plot
setEPS() 
postscript(paste0(figpath,"mae_state_comp.eps"),width = 10, height = 5)
ggplot(mae_state_comp, aes(x = state, y = rmae, color = model, group = model)) +
  geom_line(size=I(1)) +
  scale_color_manual(values = MODEL_COLORS, drop = FALSE, name="") +
  scale_fill_manual(values = MODEL_COLORS, drop = FALSE, name="") +  scale_y_continuous("Relative MAE", labels = number_format(accuracy = 0.01), breaks = seq(0,10,.05)) +
  labs(x = "", title = paste0("Relative MAE by State: ", min_plos_date," - ",max_plos_date)) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))
dev.off()


# --------------------------------------------
# plot
setEPS() 
postscript(paste0(figpath,"mae_rank_state_comp.eps"),width = 10, height = 5)
ggplot(mae_state_rank_comp, aes(x = state, y = rank, fill = model)) +
  geom_tile(color = "white", linewidth = 0.3, height = 0.9) +
  scale_fill_manual(values = MODEL_COLORS, drop = FALSE, name = "") +
  scale_y_continuous(limits = c(0, length(MODEL_COLORS)+1), breaks = 1:length(MODEL_COLORS), minor_breaks = NULL, expand = c(0,0)) +
  labs(
    x = "",
    y = "Rank (1 = best, 9 = worst)",
    title = paste0("Relative MAE Rank by State: ", min_plos_date," - ",max_plos_date)
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))
dev.off()




## Running relative error vs persistence (M(0)) — MAE or WIS
plot_running_rel_mae_base_only <- function(
    preds,
    metric             = c("mae","wis","mae_std","wis_std"),
    point_col          = "q0_500",
    wis_col            = "wis",
    persist_model      = "M(0)",
    include_conformal  = FALSE,
    exclude_models     = c("M(0)"),
    model_colors       = c(
      "M(0)"    = "darkgrey",
      "M(r,t)"  = "#fb9a99",
      "M(st,t)" = "#a6cee3",
      "M(sv,t)" = "#b2df8a",
      "M(a,t)"  = "#cab2d6",
      "M(r,v)"  = "#e31a1c",
      "M(st,v)" = "#1f78b4",
      "M(sv,v)" = "#33a02c",
      "M(a,v)"  = "#6a3d9a"
    ),
    line_size          = 0.7,
    date_labels        = "%b\n%Y",
    title              = NULL,
    subtitle           = NULL,
    eval_range         = c("all","comp"),
    min_plos_date      = ymd("2020-07-28"),
    max_plos_date      = ymd("2021-12-21")
) {
  metric <- match.arg(metric)
  DT <- data.table::as.data.table(preds)
  
  if (!inherits(DT$last_obs_date, "Date"))
    stop("`last_obs_date` must be a Date column.")
  
  req_base <- c("series_id","last_obs_date","step_ahead","end_t","model","y_true")
  miss <- setdiff(req_base, names(DT))
  if (length(miss)) stop("Missing required columns: ", paste(miss, collapse = ", "))
  
  if (metric == "mae" && !point_col %in% names(DT))
    stop("For metric='mae', missing column: ", point_col)
  if (metric == "wis" && !wis_col %in% names(DT))
    stop("For metric='wis', missing column: ", wis_col)
  
  key_cols <- c("series_id","last_obs_date","step_ahead","end_t")
  
  ## --- Baseline (persistence M(0)) ---
  if (metric == "mae") {
    base_dt <- unique(DT[model == persist_model, c(key_cols, "y_true", point_col), with = FALSE])
    if (!nrow(base_dt)) stop("No rows found for baseline: model == '", persist_model, "'")
    data.table::setnames(base_dt, point_col, "persist_pred")
  }
  if (metric == "mae_std") {
    base_dt <- unique(DT[model == persist_model, c(key_cols, "y_true", point_col,"max_cases"), with = FALSE])
    if (!nrow(base_dt)) stop("No rows found for baseline: model == '", persist_model, "'")
    data.table::setnames(base_dt, point_col, "persist_pred")
  }
  if (metric %in% c("wis","wis_std")) {
    base_dt <- unique(DT[model == persist_model, c(key_cols, "y_true", wis_col), with = FALSE])
    if (!nrow(base_dt)) stop("No rows found for baseline: model == '", persist_model, "'")
    data.table::setnames(base_dt, wis_col, "persist_wis")
  }
  
  data.table::setkeyv(base_dt, key_cols)
  
  ## --- Candidate models ---
  if (metric == "mae")      cand_cols <- c(key_cols, "model", "y_true", point_col)
  if (metric == "mae_std")  cand_cols <- c(key_cols, "model", "y_true", point_col,"max_cases")
  if (metric %in% c("wis","wis_std")) cand_cols <- c(key_cols, "model", "y_true", wis_col)
  
  cand <- DT[
    !model %in% exclude_models &
      (include_conformal | !grepl("_conform$", model)),
    cand_cols, with = FALSE
  ]
  
  need_col <- if (metric %in% c("wis","wis_std")) wis_col else point_col
  cand <- cand[!is.na(get(need_col)) & !is.na(y_true)]
  if (!nrow(cand)) stop("No candidate rows after filtering (check include_conformal/exclude_models).")
  data.table::setkeyv(cand, key_cols)
  
  ## --- Join to baseline ---
  joined <- cand[base_dt, nomatch = 0L]
  if (!nrow(joined)) stop("No overlap between candidates and baseline keys.")
  
  if (metric == "mae") {
    joined[, ae_model   := abs(get(point_col) - y_true)]
    joined[, ae_persist := abs(persist_pred    - y_true)]
    err_model_col <- "ae_model"; err_persist_col <- "ae_persist"
    ylab <- "Relative MAE"
  }
  if (metric == "mae_std") {
    joined[, ae_model   := abs(get(point_col) - y_true)/max_cases]
    joined[, ae_persist := abs(persist_pred    - y_true)/max_cases]
    err_model_col <- "ae_model"; err_persist_col <- "ae_persist"
    ylab <- "Relative Standardized MAE"
  }
  if (metric == "wis") {
    joined[, wis_model   := get(wis_col)]
    joined[, wis_persist := persist_wis]
    err_model_col <- "wis_model"; err_persist_col <- "wis_persist"
    ylab <- "Relative WIS"
  }
  if (metric == "wis_std") {
    joined[, wis_model   := get(wis_col)]
    joined[, wis_persist := persist_wis]
    err_model_col <- "wis_model"; err_persist_col <- "wis_persist"
    ylab <- "Relative Standardized WIS"
  }
  
  ## --- Daily aggregation ---
  daily <- joined[, .(
    ERR_model   = sum(get(err_model_col),   na.rm = TRUE),
    ERR_persist = sum(get(err_persist_col), na.rm = TRUE)
  ), by = .(last_obs_date, model)]
  daily$model <- factor(
    as.factor(daily$model),
    levels = c("M(0)","M(r,t)","M(st,t)","M(sv,t)","M(a,t)","M(r,v)","M(st,v)","M(sv,v)","M(a,v)")
  )
  
  data.table::setorder(daily, model, last_obs_date)
  daily[, cum_ERR_model   := cumsum(ERR_model),   by = model]
  daily[, cum_ERR_persist := cumsum(ERR_persist), by = model]
  daily[, rel_running := data.table::fifelse(cum_ERR_persist > 0, cum_ERR_model / cum_ERR_persist, NA_real_)]
  daily <- daily[!is.na(rel_running)]
  
  ## Palette
  missing_cols <- setdiff(levels(daily$model), names(model_colors))
  if (length(missing_cols)) {
    warning("No color specified for: ", paste(missing_cols, collapse = ", "), ". Using fallback gray.")
    model_colors <- c(model_colors, stats::setNames(rep("#7f7f7f", length(missing_cols)), missing_cols))
  }
  pal <- model_colors[levels(daily$model)]
  
  ## Titles
  if (is.null(title) & metric == "mae")     title <- "Running Relative MAE"
  if (is.null(title) & metric == "mae_std") title <- "Running Relative Standardized MAE"
  if (is.null(title) & metric == "wis")     title <- "Running Relative WIS"
  if (is.null(title) & metric == "wis_std") title <- "Running Relative Standardized WIS"
  if (eval_range == "comp") title <- paste0(title, ": ", min_plos_date, " - ", max_plos_date)
  
  ## --- Quarter breaks anchored to Jan/Apr/Jul/Oct ---
  q_breaks <- function(lims) {
    d0 <- as.Date(lims[1]); d1 <- as.Date(lims[2])
    start_q <- as.Date(lubridate::floor_date(d0, unit = "quarter"))
    end_q   <- as.Date(lubridate::ceiling_date(d1, unit = "quarter"))
    seq(start_q, end_q, by = "3 months")
  }
  
  ## Plot
  ggplot2::ggplot(daily, ggplot2::aes(x = last_obs_date, y = rel_running, color = model, group = model)) +
    ggplot2::geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.4) +
    ggplot2::geom_line(linewidth = line_size) +
    ggplot2::scale_y_continuous(ylab, labels = scales::number_format(accuracy = 0.01)) +
    ggplot2::scale_x_date(
      NULL,
      breaks = q_breaks,            # <-- Jan/Apr/Jul/Oct
      date_labels = date_labels,
      expand = ggplot2::expansion(mult = c(0.01, 0.03))
    ) +
    ggplot2::scale_color_manual(values = pal, breaks = names(pal),
                                guide = ggplot2::guide_legend(title = "", override.aes = list(linewidth = 1))) +
    ggplot2::labs(title = title) +
    ggplot2::theme_bw() +
    ggplot2::theme(
      plot.title.position = "plot",
      legend.position = "right",
      legend.box = "vertical",
      plot.title = element_text(hjust = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      axis.title.x = ggplot2::element_text(margin = ggplot2::margin(t = 6)),
      axis.title.y = ggplot2::element_text(margin = ggplot2::margin(r = 6))
    )
}



# --------------------------------------------
# plot
setEPS() 
postscript(paste0(figpath,"mae_running.eps"),width = 10, height = 3.5)
  plot_running_rel_mae_base_only(preds, line_size = 1, eval_range = "all", metric = "mae")
dev.off()

setEPS() 
postscript(paste0(figpath,"mae_std_running.eps"),width = 10, height = 3.5)
  plot_running_rel_mae_base_only(preds, line_size = 1, eval_range = "all", metric = "mae_std")
dev.off()


# --------------------------------------------
# plot
setEPS() 
postscript(paste0(figpath,"mae_running_comp.eps"),width = 10, height = 3.5)
  plot_running_rel_mae_base_only(subset(preds,plos_date == 1), line_size = 1,  eval_range = "comp", metric = "mae")
dev.off()

setEPS() 
postscript(paste0(figpath,"mae_std_running_comp.eps"),width = 10, height = 3.5)
  plot_running_rel_mae_base_only(subset(preds,plos_date == 1), line_size = 1, eval_range = "comp", metric = "mae_std")
dev.off()






# --------------------------------------------
# overall
wis_overall      <- preds[,list(wis = mean(wis),
                                wis_std = mean(wis_std)), by = "model"]
wis_overall$wis_m0 <- subset(wis_overall, model == "M(0)")$wis
wis_overall$wis_std_m0 <- subset(wis_overall, model == "M(0)")$wis_std
wis_overall$rwis <- wis_overall$wis / wis_overall$wis_m0
wis_overall$rwis_std <- wis_overall$wis_std / wis_overall$wis_std_m0
wis_overall <- merge(wis_overall,
                     subset(bs_uq_preds, select = c("model","rwis_q0_025","rwis_q0_975","rwis_std_q0_025","rwis_std_q0_975")),
                     by="model")


# --------------------------------------------
# plot
setEPS() 
postscript(paste0(figpath,"wis_overall.eps"),width = 7, height = 4)
  ggplot(wis_overall, aes(x = model, y = rwis, fill = model)) +
    geom_col(width = 0.65, color = "black", alpha = 0.95) +
    scale_fill_manual(values = MODEL_COLORS, drop = FALSE) +
    scale_y_continuous("Relative WIS", labels = number_format(accuracy = 0.01), breaks = seq(0,10,.1)) +
    labs(x = "", title = "Relative WIS by Model") +
    geom_text(aes(x = model, y = .8*rwis, label = round(rwis,3)))+
    theme_bw(base_size = 12) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.6, alpha = 0.8) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      plot.title.position = "plot",
      axis.text.x = element_text(angle = 0, vjust = 0, hjust = .5, size = 12))
dev.off()

setEPS() 
postscript(paste0(figpath,"wis_std_overall.eps"),width = 7, height = 4)
  ggplot(wis_overall, aes(x = model, y = rwis_std, fill = model)) +
    geom_col(width = 0.65, color = "black", alpha = 0.95) +
    scale_fill_manual(values = MODEL_COLORS, drop = FALSE) +
    scale_y_continuous("Relative WIS", labels = number_format(accuracy = 0.01), breaks = seq(0,10,.1)) +
    labs(x = "", title = "Relative Standardized WIS by Model") +
    geom_text(aes(x = model, y = .8*rwis_std, label = round(rwis_std,3)))+
    theme_bw(base_size = 12) +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.6, alpha = 0.8) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      plot.title.position = "plot",
      axis.text.x = element_text(angle = 0, vjust = 0, hjust = .5, size = 12))
dev.off()


# --------------------------------------------
# plot
setEPS() 
postscript(paste0(figpath,"wis_overall_errorbar.eps"),width = 7, height = 4)
  ggplot(wis_overall, aes(x = model, y = rwis, fill = model)) +
    geom_col(width = 0.65, color = "black") +
    scale_fill_manual(values = MODEL_COLORS, drop = FALSE) +
    scale_y_continuous("Relative WIS", labels = number_format(accuracy = 0.01), breaks = seq(0,10,.1)) +
    labs(x = "", title = "Relative WIS by Model") +
    geom_text(aes(x = model, y = .8*rwis, label = round(rwis,3)))+
    theme_bw(base_size = 12) +
    geom_errorbar(aes(ymin = rwis_q0_025, ymax = rwis_q0_975),
                  width = 0.15, linewidth = 0.6, color = "black") +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.6) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      plot.title.position = "plot",
      axis.text.x = element_text(angle = 0, vjust = 0, hjust = .5, size = 12))
dev.off()


setEPS() 
postscript(paste0(figpath,"wis_std_overall_errorbar.eps"),width = 7, height = 4)
  ggplot(wis_overall, aes(x = model, y = rwis_std, fill = model)) +
    geom_col(width = 0.65, color = "black", alpha = 0.95) +
    scale_fill_manual(values = MODEL_COLORS, drop = FALSE) +
    scale_y_continuous("Relative WIS", labels = number_format(accuracy = 0.01), breaks = seq(0,10,.1)) +
    labs(x = "", title = "Relative Standardized WIS by Model") +
    geom_text(aes(x = model, y = .8*rwis_std, label = round(rwis_std,3)))+
    theme_bw(base_size = 12) +
    geom_errorbar(aes(ymin = rwis_std_q0_025, ymax = rwis_std_q0_975),
                  width = 0.15, linewidth = 0.6, color = "black") +
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.6, alpha = 0.8) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      plot.title.position = "plot",
      axis.text.x = element_text(angle = 0, vjust = 0, hjust = .5, size = 12))
dev.off()

## comparison dates
wis_overall_comp      <- preds[plos_date == 1,list(wis = mean(wis),
                                                   wis_std = mean(wis_std)), by = "model"]
wis_overall_comp$wis_m0 <- subset(wis_overall_comp, model == "M(0)")$wis
wis_overall_comp$wis_std_m0 <- subset(wis_overall_comp, model == "M(0)")$wis_std
wis_overall_comp$rwis <- wis_overall_comp$wis / wis_overall_comp$wis_m0
wis_overall_comp$rwis_std <- wis_overall_comp$wis_std / wis_overall_comp$wis_std_m0
wis_overall_comp

# --------------------------------------------
# plot
setEPS() 
postscript(paste0(figpath,"wis_overall_comp.eps"),width = 7, height = 4)
  ggplot(subset(wis_overall_comp, model != "M(0)"), aes(x = model, y = rwis, fill = model)) +
    geom_col(width = 0.65, color = "black", alpha = 0.95) +
    scale_fill_manual(values = MODEL_COLORS, drop = FALSE) +
    scale_y_continuous("Relative WIS", labels = number_format(accuracy = 0.01), breaks = seq(0,10,.1)) +
    labs(x = "", title = paste0("Relative WIS by Model: ", min_plos_date," - ",max_plos_date)) +
    theme_bw(base_size = 12) +
    geom_text(aes(x = model, y = .8*rwis, label = round(rwis,3)))+
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.6, alpha = 0.8) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      plot.title.position = "plot",
      axis.text.x = element_text(angle = 0, vjust = 0, hjust = .5, size = 12))
dev.off()

setEPS() 
postscript(paste0(figpath,"wis_overall_comp_with_comp_line.eps"),width = 7, height = 4)
ggplot(subset(wis_overall_comp, model != "M(0)"), aes(x = model, y = rwis, fill = model)) +
  geom_col(width = 0.65, color = "black", alpha = 0.95) +
  scale_fill_manual(values = MODEL_COLORS, drop = FALSE) +
  scale_y_continuous("Relative WIS", labels = number_format(accuracy = 0.01), breaks = seq(0,10,.1),limits=c(0,max(1,max(wis_overall_comp$rwis)))) +
  labs(x = "", title = paste0("Relative WIS by Model: ", min_plos_date," - ",max_plos_date)) +
  theme_bw(base_size = 12) +
  geom_text(aes(x = model, y = .8*rwis, label = round(rwis,3)))+
  geom_hline(yintercept = 0.81, linetype = "dashed", linewidth = 0.6, alpha = 0.8) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    axis.text.x = element_text(angle = 0, vjust = 0, hjust = .5, size = 12))
dev.off()


setEPS() 
postscript(paste0(figpath,"wis_std_overall_comp.eps"),width = 7, height = 4)
  ggplot(subset(wis_overall_comp, model != "M(0)"), aes(x = model, y = rwis_std, fill = model)) +
    geom_col(width = 0.65, color = "black", alpha = 0.95) +
    scale_fill_manual(values = MODEL_COLORS, drop = FALSE) +
    scale_y_continuous("Relative WIS", labels = number_format(accuracy = 0.01), breaks = seq(0,10,.1)) +
    labs(x = "", title = paste0("Relative Standardized WIS by Model: ", min_plos_date," - ",max_plos_date)) +
    theme_bw(base_size = 12) +
    geom_text(aes(x = model, y = .8*rwis_std, label = round(rwis_std,3)))+
    geom_hline(yintercept = 1, linetype = "dashed", linewidth = 0.6, alpha = 0.8) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      plot.title.position = "plot",
      axis.text.x = element_text(angle = 0, vjust = 0, hjust = .5, size = 12))
dev.off()



# by horizon
wis_horizon      <- preds[,list(wis = mean(wis),
                                wis_std = mean(wis_std)), by = c("model","step_ahead")]
wis_horizon_m0 <- subset(wis_horizon, model == "M(0)")
names(wis_horizon_m0) <- c("model","step_ahead","wis_m0","wis_std_m0")
wis_horizon <- merge(wis_horizon, subset(wis_horizon_m0, select = c("step_ahead","wis_m0","wis_std_m0")), by=c("step_ahead"), all.x=T)
wis_horizon$rwis <- wis_horizon$wis / wis_horizon$wis_m0
wis_horizon$rwis_std <- wis_horizon$wis_std / wis_horizon$wis_std_m0
wis_horizon


# --------------------------------------------
# plot
setEPS() 
postscript(paste0(figpath,"wis_horizon.eps"),width = 7, height = 4)
  ggplot(wis_horizon, aes(x = step_ahead, y = rwis, color = model, group = model)) +
    geom_line(size=I(2)) +
    # geom_point(size=I(4), color=I("black"), shape=I(21), aes(fill = model)) +
    scale_color_manual(values = MODEL_COLORS, drop = FALSE, name="") +
    scale_fill_manual(values = MODEL_COLORS, drop = FALSE, name="") +
    scale_y_continuous("Relative WIS", labels = number_format(accuracy = 0.01), breaks = seq(0,10,.05)) +
    labs(x = "Horizon", title = "Relative WIS by Model") +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5),
      panel.grid.minor = element_blank(),
      plot.title.position = "plot",
      axis.text.x = element_text(angle = 0, hjust = .5, size = 12))
dev.off()

setEPS() 
postscript(paste0(figpath,"wis_std_horizon.eps"),width = 7, height = 4)
ggplot(wis_horizon, aes(x = step_ahead, y = rwis_std, color = model, group = model)) +
  geom_line(size=I(2)) +
  # geom_point(size=I(4), color=I("black"), shape=I(21), aes(fill = model)) +
  scale_color_manual(values = MODEL_COLORS, drop = FALSE, name="") +
  scale_fill_manual(values = MODEL_COLORS, drop = FALSE, name="") +
  scale_y_continuous("Relative WIS", labels = number_format(accuracy = 0.01), breaks = seq(0,10,.05)) +
  labs(x = "Horizon", title = "Relative WIS by Model") +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    axis.text.x = element_text(angle = 0, hjust = .5, size = 12))
dev.off()



wis_horizon_comp      <- preds[plos_date == 1,list(wis = mean(wis),
                                                   wis_std = mean(wis_std)), by = c("model","step_ahead")]
wis_horizon_m0_comp <- subset(wis_horizon_comp, model == "M(0)")
names(wis_horizon_m0_comp) <- c("model","step_ahead","wis_m0","wis_std_m0")
wis_horizon_comp <- merge(wis_horizon_comp, subset(wis_horizon_m0_comp, select = c("step_ahead","wis_m0","wis_std_m0")), by=c("step_ahead"), all.x=T)
wis_horizon_comp$rwis <- wis_horizon_comp$wis / wis_horizon_comp$wis_m0
wis_horizon_comp$rwis_std <- wis_horizon_comp$wis_std / wis_horizon_comp$wis_std_m0
wis_horizon_comp

# --------------------------------------------
# plot
setEPS() 
postscript(paste0(figpath,"wis_horizon_comp.eps"),width = 7, height = 4)
ggplot(wis_horizon_comp, aes(x = step_ahead, y = rwis, color = model, group = model)) +
  geom_line(size=I(2)) +
  scale_color_manual(values = MODEL_COLORS, drop = FALSE, name="") +
  scale_fill_manual(values = MODEL_COLORS, drop = FALSE, name="") +  scale_y_continuous("Relative WIS", labels = number_format(accuracy = 0.01), breaks = seq(0,10,.05)) +
  labs(x = "Horizon", title = paste0("Relative WIS by Model: ", min_plos_date," - ",max_plos_date)) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    axis.text.x = element_text(angle = 0, hjust = .5, size = 12))
dev.off()


setEPS() 
postscript(paste0(figpath,"wis_std_horizon_comp.eps"),width = 7, height = 4)
ggplot(wis_horizon_comp, aes(x = step_ahead, y = rwis_std, color = model, group = model)) +
  geom_line(size=I(2)) +
  scale_color_manual(values = MODEL_COLORS, drop = FALSE, name="") +
  scale_fill_manual(values = MODEL_COLORS, drop = FALSE, name="") +  scale_y_continuous("Relative WIS", labels = number_format(accuracy = 0.01), breaks = seq(0,10,.05)) +
  labs(x = "Horizon", title = paste0("Relative WIS by Model: ", min_plos_date," - ",max_plos_date)) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    axis.text.x = element_text(angle = 0, hjust = .5, size = 12))
dev.off()



# by state
wis_state      <- preds[,list(wis = mean(wis)), by = c("model","series_id")]
wis_state_m0 <- subset(wis_state, model == "M(0)")
names(wis_state_m0) <- c("model","series_id","wis_m0")
wis_state <- merge(wis_state, subset(wis_state_m0, select = c("series_id","wis_m0")), by=c("series_id"), all.x=T)
wis_state$rwis <- wis_state$wis / wis_state$wis_m0
wis_state$state <- gsub("unitedstates_","",wis_state$series_id)
wis_state_order <- wis_state[,list(avg_wis = mean(rwis)),by="state"]
wis_state_order <- wis_state_order[order(wis_state_order$avg_wis),]
wis_state$state <- factor(as.factor(wis_state$state), levels = wis_state_order$state)

wis_state_rank <- data.table::copy(wis_state)
wis_state_rank[, model := factor(model, levels = names(MODEL_COLORS))]
setorder(wis_state_rank, state, rwis, model)
wis_state_rank[, rank := seq_len(.N), by = state]   # 1..7 per state

# --------------------------------------------
# plot
setEPS() 
postscript(paste0(figpath,"wis_state.eps"),width = 10, height = 5)
ggplot(wis_state, aes(x = state, y = rwis, color = model, group = model)) +
  geom_line(size=I(1)) +
  scale_color_manual(values = MODEL_COLORS, drop = FALSE, name="") +
  scale_fill_manual(values = MODEL_COLORS, drop = FALSE, name="") +  scale_y_continuous("Relative WIS", labels = number_format(accuracy = 0.01), breaks = seq(0,10,.05)) +
  labs(x = "", title = "Relative WIS by State") +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    axis.text.x = element_text(angle = 90, vjust = .5, hjust=1))
dev.off()

# --------------------------------------------
# plot
setEPS() 
postscript(paste0(figpath,"wis_rank_state.eps"),width = 10, height = 5)
ggplot(wis_state_rank, aes(x = state, y = rank, fill = model)) +
  geom_tile(color = "white", linewidth = 0.3, height = 0.9) +
  scale_fill_manual(values = MODEL_COLORS, drop = FALSE, name = "") +
  scale_y_continuous(limits = c(0, length(MODEL_COLORS)+1), breaks = 1:length(MODEL_COLORS), minor_breaks = NULL, expand = c(0,0)) +
  labs(
    x = "",
    y = "Rank (1 = best, 9 = worst)",
    title = "Relative WIS Rank by State"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))
dev.off()




## comparison dates by state
wis_state_comp      <- preds[plos_date == 1,list(wis = mean(wis)), by = c("model","series_id")]
wis_state_m0_comp <- subset(wis_state_comp, model == "M(0)")
names(wis_state_m0_comp) <- c("model","series_id","wis_m0")
wis_state_comp <- merge(wis_state_comp, subset(wis_state_m0_comp, select = c("series_id","wis_m0")), by=c("series_id"), all.x=T)
wis_state_comp$rwis <- wis_state_comp$wis / wis_state_comp$wis_m0
wis_state_comp$state <- gsub("unitedstates_","",wis_state_comp$series_id)
wis_state_order_comp <- wis_state_comp[,list(avg_wis = mean(rwis)),by="state"]
wis_state_order_comp <- wis_state_order_comp[order(wis_state_order_comp$avg_wis),]
wis_state_comp$state <- factor(as.factor(wis_state_comp$state), levels = wis_state_order_comp$state)

wis_state_rank_comp <- data.table::copy(wis_state_comp)
wis_state_rank_comp[, model := factor(model, levels = names(MODEL_COLORS))]
setorder(wis_state_rank_comp, state, rwis, model)
wis_state_rank_comp[, rank := seq_len(.N), by = state]   # 1..7 per state

# --------------------------------------------
# plot
setEPS() 
postscript(paste0(figpath,"wis_state_comp.eps"),width = 10, height = 5)
ggplot(wis_state_comp, aes(x = state, y = rwis, color = model, group = model)) +
  geom_line(size=I(1)) +
  scale_color_manual(values = MODEL_COLORS, drop = FALSE, name="") +
  scale_fill_manual(values = MODEL_COLORS, drop = FALSE, name="") +  scale_y_continuous("Relative WIS", labels = number_format(accuracy = 0.01), breaks = seq(0,10,.05)) +
  labs(x = "", title = paste0("Relative WIS by State: ", min_plos_date," - ",max_plos_date)) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))
dev.off()

# --------------------------------------------
# plot
setEPS() 
postscript(paste0(figpath,"wis_rank_state_comp.eps"),width = 10, height = 5)
ggplot(wis_state_rank_comp, aes(x = state, y = rank, fill = model)) +
  geom_tile(color = "white", linewidth = 0.3, height = 0.9) +
  scale_fill_manual(values = MODEL_COLORS, drop = FALSE, name = "") +
  scale_y_continuous(limits = c(0, length(MODEL_COLORS)+1), breaks = 1:length(MODEL_COLORS), minor_breaks = NULL, expand = c(0,0)) +
  labs(
    x = "",
    y = "Rank (1 = best, 9 = worst)",
    title = paste0("Relative WIS Rank by State: ", min_plos_date," - ",max_plos_date)
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    panel.grid = element_blank(),
    axis.text.x = element_text(angle = 90, vjust = .5, hjust = 1))
dev.off()


# running results

# --------------------------------------------
# plot
setEPS() 
postscript(paste0(figpath,"wis_running.eps"),width = 10, height = 3.5)
  plot_running_rel_mae_base_only(preds, line_size = 1, eval_range = "all", metric = "wis", wis_col = "wis")
dev.off()

setEPS() 
postscript(paste0(figpath,"wis_std_running.eps"),width = 10, height = 3.5)
  plot_running_rel_mae_base_only(preds, line_size = 1, eval_range = "all", metric = "wis_std", wis_col = "wis_std")
dev.off()

# --------------------------------------------
# plot
setEPS() 
postscript(paste0(figpath,"wis_running_comp.eps"),width = 10, height = 3.5)
  plot_running_rel_mae_base_only(subset(preds,plos_date == 1), line_size = 1, eval_range = "comp", metric = "wis")
dev.off()

setEPS() 
postscript(paste0(figpath,"wis_std_running_comp.eps"),width = 10, height = 3.5)
  plot_running_rel_mae_base_only(subset(preds,plos_date == 1), line_size = 1, eval_range = "comp", metric = "wis_std", wis_col = "wis_std")
dev.off()



## =========================================================================================================
## 6) Coverage
cov_overall <- preds[,list(cov95 = mean(cov95),
                           cov80 = mean(cov80),
                           cov50 = mean(cov50)), by = "model"]
cov_overall_melt <- reshape2::melt(cov_overall, id.vars  = "model")
cov_overall_melt$nom_cov <- as.numeric(paste0(".",gsub("cov","",cov_overall_melt$variable)))

bs_uq_preds_melt_upper <- melt(subset(bs_uq_preds, select = c("model","cov95_q0_975","cov80_q0_975","cov50_q0_975")),
                               id.vars="model")
bs_uq_preds_melt_upper$variable <- gsub("_q0_975","",bs_uq_preds_melt_upper$variable)
names(bs_uq_preds_melt_upper)[which(names(bs_uq_preds_melt_upper) == "value")] <- "upper"

bs_uq_preds_melt_lower <- melt(subset(bs_uq_preds, select = c("model","cov95_q0_025","cov80_q0_025","cov50_q0_025")),
                               id.vars="model")
bs_uq_preds_melt_lower$variable <- gsub("_q0_025","",bs_uq_preds_melt_lower$variable)
names(bs_uq_preds_melt_lower)[which(names(bs_uq_preds_melt_lower) == "value")] <- "lower"

cov_overall_melt <- merge(cov_overall_melt,
                          bs_uq_preds_melt_lower,
                          by = c("model","variable"))
cov_overall_melt <- merge(cov_overall_melt,
                          bs_uq_preds_melt_upper,
                          by = c("model","variable"))

# --------------------------------------------
# plot
setEPS() 
postscript(paste0(figpath,"coverage_overall.eps"),width = 14, height = 4.25)
ggplot(cov_overall_melt, aes(x = model, y = value, fill = model)) +
  facet_wrap(~nom_cov)+
  geom_hline(aes(yintercept = nom_cov), linetype = "dashed", linewidth = 0.6, alpha = 0.8) +
  geom_col(width = 0.65, color = "black", alpha = 0.95) +
  scale_fill_manual(values = MODEL_COLORS, drop = FALSE) +
  scale_y_continuous("Coverage", labels = number_format(accuracy = 0.01), breaks = seq(0,10,.1), limits = c(0,1)) +
  labs(x = "", title = "Coverage by Model") +
  geom_text(aes(x = model, y = .8*value, label = round(value,2)), size=I(3))+
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    axis.text.x = element_text(angle = 0, hjust = .5))
dev.off()

# --------------------------------------------
# plot
pcomperror <- ggplot(cov_overall_melt, aes(x = model, y = value, fill = model)) +
  facet_wrap(~nom_cov)+
  geom_hline(aes(yintercept = nom_cov), linetype = "dashed", linewidth = 0.6) +
  geom_col(width = 0.65, color = "black") +
  scale_fill_manual(values = MODEL_COLORS, drop = FALSE) +
  scale_y_continuous("Coverage", labels = number_format(accuracy = 0.01), breaks = seq(0,10,.1), limits = c(0,1)) +
  labs(x = "", title = "Coverage by Model") +
  geom_text(aes(x = model, y = .8*value, label = round(value,2)), size=I(2.0))+
  geom_errorbar(aes(ymin = lower, ymax = upper),
                width = 0.15, linewidth = 0.6, color = "black") +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    axis.text.x = element_text(angle = 60, hjust = 1))

setEPS() 
postscript(paste0(figpath,"coverage_overall_errorbar.eps"),width = 9.5, height = 4)
  pcomperror
dev.off()


# comparison dates
cov_overall_comp <- preds[plos_date == 1,list(cov95 = mean(cov95),
                                              cov80 = mean(cov80),
                                              cov50 = mean(cov50)), by = "model"]
cov_overall_melt_comp <- reshape2::melt(cov_overall_comp, id.vars  = "model")
cov_overall_melt_comp$nom_cov <- as.numeric(paste0(".",gsub("cov","",cov_overall_melt_comp$variable)))

## this 0.8 coverage is from Lopez et. al. (2024) "Challenges of COVID-19 Case Forecasting in the US, 2020–2021"
cov_overall_melt_comp$covidhub <- 0.8
cov_overall_melt_comp[cov_overall_melt_comp$variable %in% c("cov50","cov80"),]$covidhub <- NA

# --------------------------------------------
# plot
pcompcov <- ggplot(subset(cov_overall_melt_comp, model != "M(0)" & nom_cov == 0.95), aes(x = model, y = value, fill = model)) +
  # facet_wrap(~nom_cov,nrow=1)+
  geom_hline(aes(yintercept = nom_cov), linetype = I(1), linewidth = 0.6) +
  geom_col(width = 0.65, color = "black") +
  scale_fill_manual(values = MODEL_COLORS, drop = FALSE) +
  scale_y_continuous("95% Coverage", labels = number_format(accuracy = 0.01), breaks = seq(0,10,.25), limits = c(0,1)) +
  labs(x = "", title = paste0("Coverage by Model: ", min_plos_date," - ",max_plos_date)) +
  theme_bw(base_size = 12) +
  geom_hline(aes(yintercept = covidhub), linetype = "dashed", linewidth = 0.6) +
  geom_text(aes(x = model, y = .8*value, label = round(value,2)), size=I(4))+
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "none",
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    axis.text.x = element_text(angle = 0, hjust = .5, size = 12))

setEPS() 
postscript(paste0(figpath,"coverage_overall_comp.eps"),width = 6, height = 4.25)
  pcompcov
dev.off()





# plot coverage comparison
setEPS() 
postscript(paste0(figpath,"wis_and_coverage_overall_comp.eps"),width = 12, height = 4.25)
  (ggplot(subset(wis_overall_comp, model != "M(0)"), aes(x = model, y = rwis, fill = model)) +
    geom_col(width = 0.65, color = "black") +
    scale_fill_manual(values = MODEL_COLORS, drop = FALSE) +
    scale_y_continuous("Relative WIS", labels = number_format(accuracy = 0.01), breaks = seq(0,10,.1),limits=c(0,max(1,max(wis_overall_comp$rwis)))) +
    labs(x = "", title = paste0("Relative WIS by Model: ", min_plos_date," - ",max_plos_date)) +
    theme_bw(base_size = 12) +
    geom_text(aes(x = model, y = .8*rwis, label = round(rwis,3)))+
    geom_hline(yintercept = 0.81, linetype = "dashed", linewidth = 0.6) +
    theme(
      plot.title = element_text(hjust = 0.5),
      legend.position = "none",
      panel.grid.minor = element_blank(),
      plot.title.position = "plot",
      axis.text.x = element_text(angle = 0, vjust = 0, hjust = .5, size = 12))) |
  pcompcov
dev.off()





## by horizon
cov_horizon <- preds[,list(cov95 = mean(cov95),
                           cov80 = mean(cov80),
                           cov50 = mean(cov50)), by = c("model","step_ahead")]

cov_horizon_melt <- reshape2::melt(cov_horizon, id.vars  = c("model","step_ahead"))
cov_horizon_melt$nom_cov <- as.numeric(paste0(".",gsub("cov","",cov_horizon_melt$variable)))

# --------------------------------------------
# plot
setEPS() 
postscript(paste0(figpath,"coverage_horizon.eps"),width = 9.5, height = 4)
ggplot(cov_horizon_melt, aes(x = step_ahead, y = value, color = model, group=model)) +
  facet_wrap(~nom_cov)+
  geom_line(size=I(1.5)) +
  geom_hline(aes(yintercept = nom_cov), linetype = "dashed", linewidth = 0.6, alpha = 0.8) +
  scale_color_manual(values = MODEL_COLORS, drop = FALSE, name="") +
  scale_fill_manual(values = MODEL_COLORS, drop = FALSE, name="") +
  scale_y_continuous("Coverage", labels = number_format(accuracy = 0.01), breaks = seq(0,10,.1), limits = c(0,1)) +
  labs(x = "Horizon", title = "Coverage by Model") +
  theme_bw(base_size = 12) +
  guides(color = guide_legend(nrow = 1))+
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "bottom",
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    axis.text.x = element_text(angle = 0, hjust = .5))
dev.off()



cov_horizon_comp <- preds[plos_date == 1,list(cov95 = mean(cov95),
                                              cov80 = mean(cov80),
                                              cov50 = mean(cov50)), by = c("model","step_ahead")]

cov_horizon_melt_comp <- reshape2::melt(cov_horizon_comp, id.vars  = c("model","step_ahead"))
cov_horizon_melt_comp$nom_cov <- as.numeric(paste0(".",gsub("cov","",cov_horizon_melt_comp$variable)))

# --------------------------------------------
# plot
setEPS() 
postscript(paste0(figpath,"coverage_horizon_comp.eps"),width = 12, height = 4.5)
ggplot(cov_horizon_melt_comp, aes(x = step_ahead, y = value, color = model, group=model)) +
  facet_wrap(~nom_cov)+
  geom_line(size=I(1.5)) +
  geom_hline(aes(yintercept = nom_cov), linetype = "dashed", linewidth = 0.6, alpha = 0.8) +
  scale_color_manual(values = MODEL_COLORS, drop = FALSE, name="") +
  scale_fill_manual(values = MODEL_COLORS, drop = FALSE, name="") +
  scale_y_continuous("Coverage", labels = number_format(accuracy = 0.01), breaks = seq(0,10,.1), limits=c(0,1)) +
  labs(x = "Horizon", title = paste0("Coverage by Model: ", min_plos_date," - ",max_plos_date)) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    axis.text.x = element_text(angle = 0, hjust = .5))
dev.off()




## by state
cov_state <- preds[,list(cov95 = mean(cov95),
                         cov80 = mean(cov80),
                         cov50 = mean(cov50)), by = c("model","series_id")]

cov_state_melt <- reshape2::melt(cov_state, id.vars  = c("model","series_id"))
cov_state_melt$nom_cov <- as.numeric(paste0(".",gsub("cov","",cov_state_melt$variable)))

cov_state_melt$state <- gsub("unitedstates_","",cov_state_melt$series_id)
cov_state_order <- data.table(cov_state_melt)[,list(avg_cov = mean(value)),by="state"]
cov_state_order <- cov_state_order[order(cov_state_order$avg_cov),]
cov_state_melt$state <- factor(as.factor(cov_state_melt$state), levels = cov_state_order$state)

# --------------------------------------------
# plot
setEPS() 
postscript(paste0(figpath,"coverage_state.eps"),width = 16, height = 8)
ggplot(cov_state_melt, aes(x = state, y = value, color = model, group=model)) +
  facet_wrap(~nom_cov, scales="free")+
  geom_line(size=I(1)) +
  geom_hline(aes(yintercept = nom_cov), linetype = "dashed", linewidth = 0.6, alpha = 0.8) +
  scale_color_manual(values = MODEL_COLORS, drop = FALSE, name="") +
  scale_fill_manual(values = MODEL_COLORS, drop = FALSE, name="") +
  scale_y_continuous("Coverage", labels = number_format(accuracy = 0.01), breaks = seq(0,10,.05)) +
  labs(x = "", title = "Coverage by Model") +
  theme_bw(base_size = 12) +
  coord_flip()+
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    axis.text.x = element_text(angle = 0, hjust = .5))
dev.off()



cov_state_comp <- preds[plos_date == 1,list(cov95 = mean(cov95),
                                            cov80 = mean(cov80),
                                            cov50 = mean(cov50)), by = c("model","series_id")]

cov_state_melt_comp <- reshape2::melt(cov_state_comp, id.vars  = c("model","series_id"))
cov_state_melt_comp$nom_cov <- as.numeric(paste0(".",gsub("cov","",cov_state_melt_comp$variable)))

cov_state_melt_comp$state <- gsub("unitedstates_","",cov_state_melt_comp$series_id)
cov_state_order_comp <- data.table(cov_state_melt_comp)[,list(avg_cov = mean(value)),by="state"]
cov_state_order_comp <- cov_state_order_comp[order(cov_state_order_comp$avg_cov),]
cov_state_melt_comp$state <- factor(as.factor(cov_state_melt_comp$state), levels = cov_state_order_comp$state)

# --------------------------------------------
# plot
setEPS() 
postscript(paste0(figpath,"coverage_state_comp.eps"),width = 16, height = 8)
ggplot(cov_state_melt_comp, aes(x = state, y = value, color = model, group=model)) +
  facet_wrap(~nom_cov, scales="free")+
  geom_line(size=I(1)) +
  geom_hline(aes(yintercept = nom_cov), linetype = "dashed", linewidth = 0.6, alpha = 0.8) +
  scale_color_manual(values = MODEL_COLORS, drop = FALSE, name="") +
  scale_fill_manual(values = MODEL_COLORS, drop = FALSE, name="") +
  scale_y_continuous("Coverage", labels = number_format(accuracy = 0.01), breaks = seq(0,10,.05)) +
  labs(x = "", title = paste0("Coverage by Model: ", min_plos_date," - ",max_plos_date)) +
  theme_bw(base_size = 12) +
  coord_flip()+
  theme(
    plot.title = element_text(hjust = 0.5),
    legend.position = "right",
    panel.grid.minor = element_blank(),
    plot.title.position = "plot",
    axis.text.x = element_text(angle = 0, hjust = .5))
dev.off()





##########################################################################
## Tables answering questions (updated for model names like M(r,t), etc.)
## - keeps non-syntactic model names from being mangled (check.names = FALSE)
## - uses three-line x-axis labels: "M(...)\nvs\nM(...)" (horizontal)
## - keeps axis labels horizontal everywhere

# keep model columns exactly as-is
tab_rmae <- data.frame(reshape2::dcast(preds_bs, bs_id ~ model, value.var="rmae"), check.names = FALSE)
tab_rmae_std <- data.frame(reshape2::dcast(preds_bs, bs_id ~ model, value.var="rmae_std"), check.names = FALSE)
tab_rwis <- data.frame(reshape2::dcast(preds_bs, bs_id ~ model, value.var="rwis"), check.names = FALSE)
tab_rwis_std <- data.frame(reshape2::dcast(preds_bs, bs_id ~ model, value.var="rwis_std"), check.names = FALSE)

# helper: 3-line "A \nvs\n B" labels
vs3 <- function(a, b) paste0(a, "\nvs\n", b)

# =========================================================
## Q1
mae_q1 <- data.table(
  type = "rMAE",
  `M(r,t) vs M(st,t)` = tab_rmae$`M(r,t)` - tab_rmae$`M(st,t)`,
  # `M(r,t) vs M(sv,t)` = tab_rmae$`M(r,t)` - tab_rmae$`M(sv,t)`,
  # `M(r,v) vs M(st,v)` = tab_rmae$`M(r,v)` - tab_rmae$`M(st,v)`,
  `M(r,v) vs M(sv,v)` = tab_rmae$`M(r,v)` - tab_rmae$`M(sv,v)`
)

wis_q1 <- data.table(
  type = "rWIS",
  `M(r,t) vs M(st,t)` = tab_rwis$`M(r,t)` - tab_rwis$`M(st,t)`,
  # `M(r,t) vs M(sv,t)` = tab_rwis$`M(r,t)` - tab_rwis$`M(sv,t)`,
  # `M(r,v) vs M(st,v)` = tab_rwis$`M(r,v)` - tab_rwis$`M(st,v)`,
  `M(r,v) vs M(sv,v)` = tab_rwis$`M(r,v)` - tab_rwis$`M(sv,v)`
)

q1df <- rbind(mae_q1, wis_q1)
q1dfmelt <- reshape2::melt(q1df, id.vars="type")

q1dfmelt_summary <- data.table(q1dfmelt)[,list(
  q.025 = quantile(value, probs = .025),
  q.975 = quantile(value, probs = .975)
), by=c("type","variable")]

q1dfcolors <- data.frame(
  variable = c(
    "M(r,t) vs M(st,t)",
    # "M(r,t) vs M(sv,t)",
    # "M(r,t) vs M(a,t)",
    # "M(r,v) vs M(st,v)",
    "M(r,v) vs M(sv,v)"
    # "M(r,v) vs M(a,v)"
  ),
  color = c("#D0B4BE","#8B5D24") #"#e4a7bb","#a85985","#aa7625","#b42f72")
)

q1dfmelt <- merge(q1dfmelt, q1dfcolors, by="variable", all.x=TRUE)

setEPS() 
postscript(paste0(figpath,"q1_boxplot.eps"), width = 9.25, height = 3)
qplot(variable, value, data = subset(q1dfmelt, type %in% c("rMAE","rWIS")),
      geom="boxplot", fill = variable, show.legend = FALSE) +
  facet_wrap(~type, nrow = 1) +
  scale_fill_manual(values = q1dfcolors$color, drop = FALSE, name="") +
  geom_point(aes(x=variable, y=q.025, fill = variable), shape = I(25), color=I("black"), size=I(2),
             data = subset(q1dfmelt_summary, type %in% c("rMAE","rWIS")), show.legend=FALSE) +
  geom_point(aes(x=variable, y=q.975, fill = variable), shape = I(24), color=I("black"), size=I(2),
             data = subset(q1dfmelt_summary, type %in% c("rMAE","rWIS")), show.legend=FALSE) +
  geom_hline(aes(yintercept = 0), linetype=I(2)) +
  xlab("") + ylab("Model A - Model B") + ggtitle("Paired Differences") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5, size = 12))
dev.off()

# =========================================================
## Q2
mae_q2 <- data.table(
  type = "rMAE",
  `M(r,t) vs M(a,t)`  = tab_rmae$`M(r,t)`  - tab_rmae$`M(a,t)`,
  `M(st,t) vs M(a,t)` = tab_rmae$`M(st,t)` - tab_rmae$`M(a,t)`,
  `M(r,v) vs M(a,v)`  = tab_rmae$`M(r,v)`  - tab_rmae$`M(a,v)`,
  `M(sv,v) vs M(a,v)` = tab_rmae$`M(sv,v)` - tab_rmae$`M(a,v)`
)

wis_q2 <- data.table(
  type = "rWIS",
  `M(r,t) vs M(a,t)`  = tab_rwis$`M(r,t)`  - tab_rwis$`M(a,t)`,
  `M(st,t) vs M(a,t)` = tab_rwis$`M(st,t)` - tab_rwis$`M(a,t)`,
  `M(r,v) vs M(a,v)`  = tab_rwis$`M(r,v)`  - tab_rwis$`M(a,v)`,
  `M(sv,v) vs M(a,v)` = tab_rwis$`M(sv,v)` - tab_rwis$`M(a,v)`
)

q2df <- rbind(mae_q2, wis_q2)
q2dfmelt <- reshape2::melt(q2df, id.vars="type")

q2dfmelt_summary <- data.table(q2dfmelt)[,list(
  q.025 = quantile(value, probs = .025),
  q.975 = quantile(value, probs = .975)
), by=c("type","variable")]

q2dfcolors <- data.frame(
  variable = c(
    "M(r,t) vs M(a,t)",
    "M(st,t) vs M(a,t)",
    "M(r,v) vs M(a,v)",
    "M(sv,v) vs M(a,v)"
  ),
  color = c("#e4a7bb","#b9c1dd","#b42f72","#547b74")
)

setEPS() 
postscript(paste0(figpath,"q2_boxplot.eps"), width = 9.25, height = 3)
qplot(variable, value, data = subset(q2dfmelt, type %in% c("rMAE","rWIS")),
      geom="boxplot", fill = variable, show.legend = FALSE) +
  facet_wrap(~type, nrow = 1) +
  geom_hline(aes(yintercept = 0), linetype=I(2)) +
  scale_fill_manual(values = q2dfcolors$color, drop = FALSE, name="") +
  scale_color_manual(values = q2dfcolors$color, drop = FALSE, name="") +
  xlab("") +
  geom_point(aes(x=variable, y=q.025, fill = variable), shape = I(25), color=I("black"), size=I(2),
             data = subset(q2dfmelt_summary, type %in% c("rMAE","rWIS")), show.legend=FALSE) +
  geom_point(aes(x=variable, y=q.975, fill = variable), shape = I(24), color=I("black"), size=I(2),
             data = subset(q2dfmelt_summary, type %in% c("rMAE","rWIS")), show.legend=FALSE) +
  ylab("Model A - Model B") + ggtitle("Paired Differences") +
  theme(plot.title = element_text(hjust = 0.5))
dev.off()

# =========================================================
## Q3
mae_q3 <- data.table(
  type = "rMAE",
  `M(st,t) vs M(sv,v)` = tab_rmae$`M(st,t)` - tab_rmae$`M(sv,v)`,
  `M(sv,t) vs M(st,t)` = tab_rmae$`M(sv,t)` - tab_rmae$`M(st,t)`,
  `M(st,v) vs M(sv,v)` = tab_rmae$`M(st,v)` - tab_rmae$`M(sv,v)`
  
)

wis_q3 <- data.table(
  type = "rWIS",
  `M(st,t) vs M(sv,v)` = tab_rwis$`M(st,t)` - tab_rwis$`M(sv,v)`,
  `M(sv,t) vs M(st,t)` = tab_rwis$`M(sv,t)` - tab_rwis$`M(st,t)`,
  `M(st,v) vs M(sv,v)` = tab_rwis$`M(st,v)` - tab_rwis$`M(sv,v)`
)

q3df <- rbind(mae_q3, wis_q3)
q3dfmelt <- reshape2::melt(q3df, id.vars="type")

q3dfmelt_summary <- data.table(q3dfmelt)[,list(
  q.025 = quantile(value, probs = .025),
  q.975 = quantile(value, probs = .975)
), by=c("type","variable")]

q3dfcolors <- data.frame(
  variable = c("M(st,t) vs M(sv,v)","M(sv,t) vs M(st,t)","M(st,v) vs M(sv,v)"),
  color = c("#6DB788","#ACD7B7","#298C70")
)


setEPS() 
postscript(paste0(figpath,"q3_boxplot.eps"), width = 9.25, height = 3)
qplot(variable, value, data = subset(q3dfmelt, type %in% c("rMAE","rWIS")),
      geom="boxplot", fill = variable, show.legend = FALSE) +
  facet_wrap(~type, nrow = 1) +
  geom_hline(aes(yintercept = 0), linetype=I(2)) +
  xlab("") +
  geom_point(aes(x=variable, y=q.025, fill = variable), shape = I(25), color=I("black"), size=I(2),
             data = subset(q3dfmelt_summary, type %in% c("rMAE","rWIS")), show.legend=FALSE) +
  geom_point(aes(x=variable, y=q.975, fill = variable), shape = I(24), color=I("black"), size=I(2),
             data = subset(q3dfmelt_summary, type %in% c("rMAE","rWIS")), show.legend=FALSE) +
  scale_fill_manual(values = q3dfcolors$color, drop = FALSE, name="") +
  ylab("Model A - Model B") + ggtitle("Paired Differences") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))
dev.off()




# =========================================================
## Q4
mae_q4 <- data.table(
  type = "rMAE",
  `M(r,t) vs M(r,v)`   = tab_rmae$`M(r,t)`  - tab_rmae$`M(r,v)`,
  `M(st,t) vs M(st,v)` = tab_rmae$`M(st,t)` - tab_rmae$`M(st,v)`,
  `M(sv,t) vs M(sv,v)` = tab_rmae$`M(sv,t)` - tab_rmae$`M(sv,v)`,
  `M(a,t) vs M(a,v)`   = tab_rmae$`M(a,t)`  - tab_rmae$`M(a,v)`
)


wis_q4 <- data.table(
  type = "rWIS",
  `M(r,t) vs M(r,v)`   = tab_rwis$`M(r,t)`  - tab_rwis$`M(r,v)`,
  `M(st,t) vs M(st,v)` = tab_rwis$`M(st,t)` - tab_rwis$`M(st,v)`,
  `M(sv,t) vs M(sv,v)` = tab_rwis$`M(sv,t)` - tab_rwis$`M(sv,v)`,
  `M(a,t) vs M(a,v)`   = tab_rwis$`M(a,t)`  - tab_rwis$`M(a,v)`
)

q4df <- rbind(mae_q4, wis_q4)
q4dfmelt <- reshape2::melt(q4df, id.vars="type")

q4dfmelt_summary <- data.table(q4dfmelt)[,list(
  q.025 = quantile(value, probs = .025),
  q.975 = quantile(value, probs = .975)
), by=c("type","variable")]

q4dfcolors <- data.frame(
  variable = c(
    "M(r,t) vs M(r,v)",
    "M(st,t) vs M(st,v)",
    "M(sv,t) vs M(sv,v)",
    "M(a,t) vs M(a,v)"
  ),
  color = c("#ef7271","#7baacd","#86c368","#a388bb")
)


setEPS() 
postscript(paste0(figpath,"q4_boxplot.eps"), width = 9.25, height = 3)
qplot(variable, value, data = subset(q4dfmelt, type %in% c("rMAE","rWIS")),
      geom="boxplot", fill = variable, show.legend = FALSE) +
  facet_wrap(~type, nrow = 1) +
  geom_hline(aes(yintercept = 0), linetype=I(2)) +
  geom_point(aes(x=variable, y=q.025, fill = variable), shape = I(25), color=I("black"), size=I(2),
             data = subset(q4dfmelt_summary, type %in% c("rMAE","rWIS")), show.legend=FALSE) +
  geom_point(aes(x=variable, y=q.975, fill = variable), shape = I(24), color=I("black"), size=I(2),
             data = subset(q4dfmelt_summary, type %in% c("rMAE","rWIS")), show.legend=FALSE) +
  xlab("") +
  scale_fill_manual(values = q4dfcolors$color, drop = FALSE, name="") +
  ylab("Model A - Model B") + ggtitle("Paired Differences") +
  theme(plot.title = element_text(hjust = 0.5),
        axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5))
dev.off()




# ------------------------------------------------------------------------------
# Make individual state/date forecasts
# ------------------------------------------------------------------------------

panel_forecast_plot <- function(dat, preds, last_date, state, history_window = 0L, show_truth = T) {
  suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(patchwork)
  })
  
  MODEL_COLORS <- c(
    "M(0)"  = "darkgrey",
    "M(r,t)"  = "#fb9a99",
    "M(st,t)" = "#a6cee3",
    "M(sv,t)" = "#b2df8a",
    "M(a,t)"  = "#cab2d6",
    "M(r,v)"  = "#e31a1c",
    "M(st,v)" = "#1f78b4",
    "M(sv,v)" = "#33a02c",
    "M(a,v)"  = "#6a3d9a"
  )
  
  # ---- checks ----
  if (!is.list(dat) || length(dat) == 0) stop("`dat` must be a non-empty list of data.tables.")
  preds <- as.data.table(preds)
  
  last_date <- as.IDate(last_date)
  if (is.na(last_date)) stop("`last_date` must be coercible to IDate/Date.")
  
  history_window <- as.integer(history_window)
  if (!is.finite(history_window) || history_window < 0L) stop("`history_window` must be an integer >= 0.")
  
  # find state's historical dt
  dat_dt <- NULL
  for (i in seq_along(dat)) {
    dt_i <- as.data.table(dat[[i]])
    if ("series_id" %in% names(dt_i) && any(dt_i$series_id == state)) {
      dat_dt <- dt_i[series_id == state]
      break
    }
  }
  if (is.null(dat_dt) || nrow(dat_dt) == 0) stop("Could not find `state` in `dat` list.")
  
  req_dat_cols <- c("series_id", "date", "cases")
  miss_dat <- setdiff(req_dat_cols, names(dat_dt))
  if (length(miss_dat)) stop("`dat` is missing columns: ", paste(miss_dat, collapse = ", "))
  
  req_pred_cols <- c(
    "series_id", "ref_date", "date", "y_true",
    "q0.0250","q0.1000","q0.2500","q0.5000","q0.7500","q0.9000","q0.9750",
    "model"
  )
  miss_pred <- setdiff(req_pred_cols, names(preds))
  if (length(miss_pred)) stop("`preds` is missing columns: ", paste(miss_pred, collapse = ", "))
  
  # coerce dates
  dat_dt[, date := as.IDate(date)]
  preds[,  date := as.IDate(date)]
  preds[, ref_date := as.IDate(ref_date)]
  
  # ---- historical window ----
  hist_start <- last_date - history_window * 7L
  hist_dt <- dat_dt[date >= hist_start & date <= last_date]
  
  # ---- choose forecast ref_date: nearest not later than requested ----
  p_state <- preds[series_id == state & !is.na(ref_date)]
  if (nrow(p_state) == 0) stop("No rows in `preds` for this state.")
  
  eligible_dates <- p_state[ref_date <= last_date, unique(ref_date)]
  if (length(eligible_dates) == 0) {
    stop("No forecast ref_date in `preds` is <= requested last_date for this state.")
  }
  ref_date_fcst <- max(eligible_dates)
  
  # subset forecasts to that chosen forecast ref date
  fc <- p_state[ref_date == ref_date_fcst]
  
  # dedupe (keep 1 row per model x forecast date)
  fc <- unique(
    fc[, .(model, series_id, ref_date, date, y_true,
           q0.0250, q0.1000, q0.2500, q0.5000, q0.7500, q0.9000, q0.9750)],
    by = c("model", "series_id", "ref_date", "date")
  )
  
  # enforce model ordering
  model_order <- c("M(0)","M(r,t)","M(st,t)","M(sv,t)","M(a,t)","M(r,v)","M(st,v)","M(sv,v)","M(a,v)")
  fc[, model := factor(as.character(model), levels = model_order)]
  
  # ---- MAE per model (for titles) ----
  mae_by_model <- fc[
    is.finite(q0.5000) & is.finite(y_true),
    .(mae = mean(abs(q0.5000 - y_true))),
    by = model
  ]
  fc <- mae_by_model[fc, on = "model"]   # adds column `mae` to fc
  
  
  # ---- per-panel plot ----
  if(show_truth){myalpha <- 1}else{myalpha <- 0}
  make_one <- function(model_code, myalpha = 1) {
    fcm <- fc[model == model_code][order(date)]
    if (nrow(fcm) == 0) {
      return(
        ggplot() +
          theme_void() +
          ggtitle(as.character(model_code)) +
          theme(plot.title = element_text(hjust = 0.5))
      )
    }
    
    col <- MODEL_COLORS[[as.character(model_code)]]
    
    mae_val <- fcm$mae[1]
    mae_txt <- if (is.finite(mae_val)) {
      paste0("MAE=", scales::comma(round(mae_val, 2), accuracy = 0.01))
    } else {
      "MAE=NA"
    }
    
    if(myalpha == 0){mae_txt = ""}
      
    ggplot() +
      geom_line(
        data = hist_dt,
        aes(x = date, y = cases),
        linewidth = 0.7,
        color = "black"
      ) +
      geom_ribbon(
        data = fcm,
        aes(x = date, ymin = q0.0250, ymax = q0.9750),
        fill = col, alpha = 0.4
      ) +
      geom_ribbon(
        data = fcm,
        aes(x = date, ymin = q0.1000, ymax = q0.9000),
        fill = col, alpha = 0.4
      ) +
      geom_ribbon(
        data = fcm,
        aes(x = date, ymin = q0.2500, ymax = q0.7500),
        fill = col, alpha = 0.4
      ) +
      geom_line(
        data = fcm,
        aes(x = date, y = q0.5000),
        linewidth = 0.9,
        color = col
      ) +
      geom_point(
        data = fcm[is.finite(y_true)],
        aes(x = date, y = y_true),
        color = "black",
        size = 1.6,
        alpha = myalpha,
      ) +
      geom_vline(xintercept = as.Date(last_date), linetype = "dashed", linewidth = 0.5, alpha = 0.8) +
      ggtitle(paste0(as.character(model_code), " | ", mae_txt)) +
      theme_bw(base_size = 11) +
      theme(
        plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank()
      ) +
      labs(x = NULL, y = NULL)
  }
  
  
  blank_panel <- ggplot() + theme_void()
  
  p_top <- make_one("M(0)", myalpha) | make_one("M(r,t)", myalpha) | make_one("M(st,t)", myalpha) | make_one("M(sv,t)", myalpha) | make_one("M(a,t)", myalpha)
  p_bot <- blank_panel   | make_one("M(r,v)", myalpha) | make_one("M(st,v)", myalpha) | make_one("M(sv,v)", myalpha) | make_one("M(a,v)", myalpha)
  
  (p_top / p_bot) +
    plot_annotation(
      title = paste0(
        state,
        " | requested last_date = ", as.character(last_date),
        " | forecast ref_date used = ", as.character(ref_date_fcst),
        " | history_window = ", history_window
      ),
      theme = theme(plot.title = element_text(hjust = 0.5))
    )
}



## spot check some forecasts
unqdates <- sort(unique(preds$ref_date))
unqstates <- unique(preds$series_id)


## (formerly) problematic forecasts
panel_forecast_plot(test_covid, preds, last_date = "2022-02-14",
                    state = "unitedstates_maine", history_window = 20)
panel_forecast_plot(test_covid, preds, last_date = "2022-11-07",
                    state = "unitedstates_florida", history_window = 20)
panel_forecast_plot(test_covid, preds, last_date = "2022-11-14",
                    state = "unitedstates_florida", history_window = 20)
panel_forecast_plot(test_covid, preds, last_date = "2021-12-27",
                    state = "unitedstates_newyork", history_window = 20)
panel_forecast_plot(test_covid, preds, last_date = "2022-01-03",
                    state = "unitedstates_newyork", history_window = 20)




## the problem of Omicron
panel_forecast_plot(test_covid, preds, last_date = "2021-12-27",
                    state = "unitedstates_newyork", history_window = 20)
panel_forecast_plot(test_covid, preds, last_date = "2022-01-03",
                    state = "unitedstates_california", history_window = 20)
panel_forecast_plot(test_covid, preds, last_date = "2022-01-03",
                    state = "unitedstates_texas", history_window = 20)
panel_forecast_plot(test_covid, preds, last_date = "2021-12-27",
                    state = "unitedstates_florida", history_window = 20)



## interesting TC vs VAC difference
panel_forecast_plot(test_covid, preds, last_date = "2022-01-03",
                    state = "unitedstates_iowa", history_window = 20)
panel_forecast_plot(test_covid, preds, last_date = "2022-01-10",
                    state = "unitedstates_iowa", history_window = 20)
panel_forecast_plot(test_covid, preds, last_date = "2022-01-15",
                    state = "unitedstates_minnesota", history_window = 20)




# --------------------------------------------------------------------------
# show example forecasts for all times
# --------------------------------------------------------------------------

panel_forecast_plot_5wk_grid <- function(
    dat,
    preds,
    last_date,
    state,
    sqrt_y = TRUE,
    layout_2col = TRUE,
    history_window = 0L,
    show_truth = TRUE,
    step_weeks = 5L
) {
  suppressPackageStartupMessages({
    library(data.table)
    library(ggplot2)
    library(patchwork)
    library(lubridate)
    library(scales)
  })
  
  blend_hex <- function(fg, alpha, bg = "#FFFFFF") {
    if (!is.numeric(alpha) || length(alpha) != 1 || is.na(alpha) || alpha < 0 || alpha > 1) {
      stop("`alpha` must be a single number between 0 and 1.")
    }
    
    fg_rgb <- grDevices::col2rgb(fg)[, 1]
    bg_rgb <- grDevices::col2rgb(bg)[, 1]
    
    blended <- round(alpha * fg_rgb + (1 - alpha) * bg_rgb)
    
    grDevices::rgb(
      red = blended[1],
      green = blended[2],
      blue = blended[3],
      maxColorValue = 255
    )
  }
  
  MODEL_COLORS <- c(
    "M(r,t)"  = "#fb9a99",
    "M(st,t)" = "#a6cee3",
    "M(sv,t)" = "#b2df8a",
    "M(a,t)"  = "#cab2d6",
    "M(r,v)"  = "#e31a1c",
    "M(st,v)" = "#1f78b4",
    "M(sv,v)" = "#33a02c",
    "M(a,v)"  = "#6a3d9a"
  )
  
  if (!is.list(dat) || length(dat) == 0) stop("`dat` must be a non-empty list of data.tables.")
  preds <- as.data.table(preds)
  
  last_date <- as.IDate(last_date)
  if (is.na(last_date)) stop("`last_date` must be coercible to IDate/Date.")
  
  history_window <- as.integer(history_window)
  if (!is.finite(history_window) || history_window < 0L) {
    stop("`history_window` must be an integer >= 0.")
  }
  
  step_weeks <- as.integer(step_weeks)
  if (!is.finite(step_weeks) || step_weeks <= 0L) {
    stop("`step_weeks` must be an integer >= 1.")
  }
  
  req_pred_cols <- c(
    "series_id", "ref_date", "date", "y_true",
    "q0.0250", "q0.1000", "q0.2500", "q0.5000", "q0.7500", "q0.9000", "q0.9750",
    "model"
  )
  miss_pred <- setdiff(req_pred_cols, names(preds))
  if (length(miss_pred)) stop("`preds` is missing columns: ", paste(miss_pred, collapse = ", "))
  
  preds[, date := as.IDate(date)]
  preds[, ref_date := as.IDate(ref_date)]
  
  req_dat_cols <- c("series_id", "date", "cases")
  
  if (identical(state, "all")) {
    dat_all <- rbindlist(lapply(dat, as.data.table), fill = TRUE)
    miss_dat <- setdiff(req_dat_cols, names(dat_all))
    if (length(miss_dat)) stop("`dat` is missing columns: ", paste(miss_dat, collapse = ", "))
    
    dat_all[, date := as.IDate(date)]
    dat_dt <- dat_all[
      is.finite(cases),
      .(cases = sum(cases, na.rm = TRUE)),
      by = .(date)
    ][order(date)]
    dat_dt[, series_id := "all"]
  } else {
    dat_dt <- NULL
    for (i in seq_along(dat)) {
      dt_i <- as.data.table(dat[[i]])
      if ("series_id" %in% names(dt_i) && any(dt_i$series_id == state)) {
        dat_dt <- dt_i[series_id == state]
        break
      }
    }
    if (is.null(dat_dt) || nrow(dat_dt) == 0) stop("Could not find `state` in `dat` list.")
    
    miss_dat <- setdiff(req_dat_cols, names(dat_dt))
    if (length(miss_dat)) stop("`dat` is missing columns: ", paste(miss_dat, collapse = ", "))
    
    dat_dt[, date := as.IDate(date)]
  }
  
  hist_dt <- dat_dt
  
  if (identical(state, "all")) {
    p_state <- preds[!is.na(ref_date)]
    if (nrow(p_state) == 0) stop("No rows in `preds` (after removing NA ref_date).")
    
    fc0 <- p_state[
      ,
      .(
        y_true  = sum(y_true, na.rm = TRUE),
        q0.0250 = sum(q0.0250, na.rm = TRUE),
        q0.1000 = sum(q0.1000, na.rm = TRUE),
        q0.2500 = sum(q0.2500, na.rm = TRUE),
        q0.5000 = sum(q0.5000, na.rm = TRUE),
        q0.7500 = sum(q0.7500, na.rm = TRUE),
        q0.9000 = sum(q0.9000, na.rm = TRUE),
        q0.9750 = sum(q0.9750, na.rm = TRUE)
      ),
      by = .(model, ref_date, date)
    ]
    fc0[, series_id := "all"]
  } else {
    p_state <- preds[series_id == state & !is.na(ref_date)]
    if (nrow(p_state) == 0) stop("No rows in `preds` for this state.")
    fc0 <- p_state
  }
  
  available_ref <- sort(unique(fc0$ref_date))
  if (!length(available_ref)) stop("No available `ref_date` values in `preds` for this selection.")
  
  available_ref <- available_ref[available_ref <= last_date]
  if (!length(available_ref)) stop("No available `ref_date` values on/before `last_date` for this selection.")
  
  step_days <- step_weeks * 7L
  min_ref <- min(available_ref)
  k_back <- as.integer(floor(as.numeric((last_date - min_ref) / step_days)))
  
  cand_ref <- as.IDate(sort(unique(last_date - (0:k_back) * step_days)))
  
  ref_grid <- cand_ref[cand_ref %in% available_ref]
  if (!length(ref_grid)) {
    stop(
      "No ref_dates on the 5-week grid were found in `preds` for this selection. ",
      "Try a different `last_date` or `step_weeks`."
    )
  }
  
  fc <- fc0[ref_date %in% ref_grid]
  
  fc <- unique(
    fc[, .(
      model, series_id, ref_date, date, y_true,
      q0.0250, q0.1000, q0.2500, q0.5000, q0.7500, q0.9000, q0.9750
    )],
    by = c("model", "series_id", "ref_date", "date")
  )
  
  model_order <- c("M(r,t)", "M(st,t)", "M(sv,t)", "M(a,t)", "M(r,v)", "M(st,v)", "M(sv,v)", "M(a,v)")
  fc[, model := factor(as.character(model), levels = model_order)]
  
  make_one <- function(model_code) {
    fcm <- fc[model == model_code][order(ref_date, date)]
    if (nrow(fcm) == 0) {
      return(
        ggplot() +
          theme_void() +
          ggtitle(as.character(model_code)) +
          theme(plot.title = element_text(hjust = 0.5))
      )
    }
    
    col <- MODEL_COLORS[[as.character(model_code)]]
    
    ribbon_col_95 <- blend_hex(col, 0.20, bg = "#FFFFFF")
    ribbon_col_80 <- blend_hex(col, 0.35, bg = "#FFFFFF")
    ribbon_col_50 <- blend_hex(col, 0.50, bg = "#FFFFFF")
    
    sqrt_breaks <- function(n = 6) {
      function(limits) {
        r <- range(pmax(limits, 0), finite = TRUE)
        if (!all(is.finite(r))) return(NULL)
        b <- seq(sqrt(r[1]), sqrt(r[2]), length.out = n)
        unique(round(pmax(0, b^2), -3))
      }
    }
    
    mxdate <- max(fcm$date)
    
    p <- ggplot() +
      geom_ribbon(
        data = fcm,
        aes(x = date, ymin = q0.0250, ymax = q0.9750, group = ref_date),
        fill = ribbon_col_95,
        color = NA
      ) +
      geom_ribbon(
        data = fcm,
        aes(x = date, ymin = q0.1000, ymax = q0.9000, group = ref_date),
        fill = ribbon_col_80,
        color = NA
      ) +
      geom_ribbon(
        data = fcm,
        aes(x = date, ymin = q0.2500, ymax = q0.7500, group = ref_date),
        fill = ribbon_col_50,
        color = NA
      ) +
      geom_line(
        data = data.frame(hist_dt)[data.frame(hist_dt)$date >= (ymd("2020-06-01") - 140) &
                                     data.frame(hist_dt)$date <= mxdate, ],
        aes(x = date, y = cases),
        linewidth = 0.7,
        color = "black",
        inherit.aes = FALSE,
        show.legend = FALSE
      ) +
      geom_point(
        data = fcm,
        aes(x = date, y = q0.5000, group = ref_date),
        fill = col,
        shape = 21,
        color = "black"
      ) +
      ggtitle(as.character(model_code)) +
      theme_bw(base_size = 11) +
      theme(
        plot.title = element_text(hjust = 0.5),
        panel.grid.minor = element_blank()
      ) +
      labs(
        x = NULL,
        y = "Total Cases"
      )
    
    if (sqrt_y) {
      p <- p +
        scale_y_sqrt(
          breaks = sqrt_breaks(6),
          labels = label_comma()
        )
    }
    
    p
  }
  
  if (layout_2col) {
    p <- (make_one("M(r,t)")  | make_one("M(r,v)")) /
      (make_one("M(st,t)") | make_one("M(st,v)")) /
      (make_one("M(sv,t)") | make_one("M(sv,v)")) /
      (make_one("M(a,t)")  | make_one("M(a,v)"))
  } else {
    p <- (make_one("M(r,t)")  | make_one("M(st,t)") | make_one("M(sv,t)") | make_one("M(a,t)")) /
      (make_one("M(r,v)")  | make_one("M(st,v)") | make_one("M(sv,v)") | make_one("M(a,v)"))
  }
  
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
  
  p +
    plot_annotation(
      title = state_from_series_id(state),
      theme = theme(plot.title = element_text(hjust = 0.5, size = 20))
    )
}




# --------------------------------------
# make the plot
# --------------------------------------

setEPS()
postscript(paste0(figpath,"visualize_fcsts_nm.eps"),width = 10, height = 9)
  panel_forecast_plot_5wk_grid(dat = test_covid, 
                               preds = preds, 
                               last_date = max(unique(preds$ref_date)[unique(preds$ref_date) < ymd("2022-12-24")]),
                               state = "unitedstates_newmexico",
                               sqrt_y = T,
                               step_weeks = 6,
                               layout_2col = T)
dev.off()


setEPS() 
postscript(paste0(figpath,"visualize_fcsts_usa.eps"),width = 10, height = 9)
panel_forecast_plot_5wk_grid(dat = test_covid, 
                             preds = preds, 
                             last_date = max(unique(preds$ref_date)[unique(preds$ref_date) < ymd("2021-12-01")]),
                             state = "all",
                             sqrt_y = F,
                             step_weeks = 6,
                             layout_2col = T)
dev.off()




