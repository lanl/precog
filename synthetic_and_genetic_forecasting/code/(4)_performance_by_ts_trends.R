## Where do models outperform? Inflection points?
## I need to compute an MAE, rMAE, WIS, and rWIS by state and forecast date.
## I then need to assign the class label 

## -----------------------------------------------------------------------------
## load libraries
library(ggplot2)
library(data.table)
library(patchwork)
library(lubridate)
theme_set(theme_bw())

## -----------------------------------------------------------------------------
## set paths
datatrendspath <- paste0(here::here("synthetic_and_genetic_forecasting", "data_raw", "ts_phases"), "/") 
figpath <- paste0(here::here("synthetic_and_genetic_forecasting", "figs"), "/")  
fcstpath <- paste0(here::here("synthetic_and_genetic_forecasting", "output"), "/") 
datapath <- paste0(here::here("synthetic_and_genetic_forecasting", "data_clean"), "/")

 
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


cols <- c(
  "Impending Rise" = "#80CDC1",  # blue
  "Rise" = "#01665E",  # light blue
  "Impending Fall" = "#F1B6DA",  # red
  "Fall" = "#C51B7D",  # light red
  "Start" = "#ffffbf"   # light grey
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
preds$ae_std <- abs(preds$y_true - preds$q0_500)/preds$max_cases

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







## make preds for m0
preds_m0 <- subset(preds, model == "M(0)", select = c("series_id","end_t","step_ahead","last_obs_date","ae","wis"))
names(preds_m0) <- c("series_id","end_t","step_ahead","last_obs_date","ae_m0","wis_m0")

## merge with preds
preds <- merge(preds, preds_m0, by=c("series_id","end_t","step_ahead","last_obs_date"), all.x = T)


## -----------------------------------------------------------------------------
## start the 

lf_raw <- sort(list.files(datatrendspath))
lf_raw <- grep(".csv",lf_raw,value=T)


## loop over states
results_list <- list()
for(i in 1:length(lf_raw)){
  
  ## get the state
  low_state <- gsub(" ","",tolower(gsub(".csv","",lf_raw[i])))
  
  ## read in raw fake data
  df_raw <- fread(paste0(datatrendspath,lf_raw[i]))
  df_raw$Time <- df_raw$Time + days(1)
  
  
  v <- df_raw$FAKE_03
  
  label <- rep("Start", length(v))
  label[v ==  1] <- "Impending Rise"
  label[v == -1] <- "Impending Fall"
  
  # index of most recent non-zero at or before each position (0 if none)
  last_nz <- cummax(ifelse(v != 0, seq_along(v), 0L))
  
  # safely look up the most recent non-zero value (0 if none)
  prev_nz_val <- integer(length(v))
  idx <- last_nz > 0
  prev_nz_val[idx] <- v[last_nz[idx]]
  
  # assign labels for zeros based on previous non-zero
  is0 <- v == 0
  label[is0 & prev_nz_val ==  1] <- "Rise"
  label[is0 & prev_nz_val == -1] <- "Fall"
  
  df_raw$Phase <- label
  
  
  ## subset preds
  temp_preds <- subset(preds, series_id == paste0("unitedstates_",low_state))
  
  ## combine:
  ## df_raw $Time matches temp_preds$last_obs_date
  temp_preds2 <- merge(temp_preds,
                       subset(df_raw, select = c("Time","FAKE_03","Phase")),
                       by.x = "last_obs_date",
                       by.y = "Time",
                       all.x=T)
  
  results_list[[i]] <- temp_preds2
  

}
results <- rbindlist(results_list)





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
    ), by = c("model", "Phase")]
    
    # Baseline model values (M0)
    base <- temp_sum_overall[model == baseline_model,
                             .(Phase = Phase, mae_M0 = mae, mae_std_M0 = mae_std, wis_M0 = wis, wis_std_M0 = wis_std)]
    
    # If baseline missing for some reason, fill with NA to avoid errors
    if (nrow(base) == 0L) {
      base <- data.table(mae_M0 = NA_real_, mae_std_M0 = NA_real_, wis_M0 = NA_real_, wis_std_M0 = NA_real_)
    }
    
    temp_sum_overall <- merge(temp_sum_overall,
                              base,
                              by="Phase")
    
    
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



## get results by model and phase
results$Phase <- results$Phase
results_summary <- results[,list(mae = mean(ae),
                                 wis = mean(wis)), by=c("model","Phase")]


## get the bootstraps
Nb <- 5000
preds_bs <- bootstrap_rmse_like_summary_fast(results, Nb = Nb, Nwindow = 15)
preds_bs$Phase <- preds_bs$Phase


## reshape
preds_bs_mae <- reshape2::dcast(preds_bs, Phase + bs_id ~ model, value.var = "rmae")
preds_bs_mae$type <- "rMAE"
preds_bs_wis <- reshape2::dcast(preds_bs, Phase + bs_id ~ model, value.var = "rwis")
preds_bs_wis$type <- "rWIS"
preds_bs2 <- rbind(preds_bs_mae, preds_bs_wis)

preds_bs2$`M(r,t) vs M(r,v)` <- (preds_bs2$`M(r,t)` - preds_bs2$`M(r,v)`)
preds_bs2$`M(st,t) vs M(st,v)` <- (preds_bs2$`M(st,t)` - preds_bs2$`M(st,v)`)
preds_bs2$`M(sv,t) vs M(sv,v)` <- (preds_bs2$`M(sv,t)` - preds_bs2$`M(sv,v)`)
preds_bs2$`M(a,t) vs M(a,v)` <- (preds_bs2$`M(a,t)` - preds_bs2$`M(a,v)`)


preds_bs_melt <- data.table(reshape2::melt(subset(preds_bs2, select = c("Phase","bs_id","type","M(r,t) vs M(r,v)",
                                                                        "M(st,t) vs M(st,v)","M(sv,t) vs M(sv,v)","M(a,t) vs M(a,v)")),
                                id.vars=c("Phase","bs_id","type")))

preds_bs_melt2 <- preds_bs_melt[,list(q.025 = quantile(value,probs=.025),
                                      med = median(value),
                                      avg = mean(value),
                                      q.975 = quantile(value, probs=.975)), by=c("Phase","type","variable")]

preds_bs_melt2$Phase <- factor(as.factor(preds_bs_melt2$Phase), levels = c("Start","Impending Rise","Rise","Impending Fall","Fall"))


q4dfcolors <- data.frame(variable = c("M(r,t) vs M(r,v)","M(st,t) vs M(st,v)","M(sv,t) vs M(sv,v)","M(a,t) vs M(a,v)"),
                         color = c("#ef7271","#7baacd","#86c368","#a388bb"))



## plot it
pd <- position_dodge(width = 0.6)

boxplot_plot <- ggplot(subset(preds_bs_melt2, Phase != "Start")) +
  geom_errorbar(
    aes(x = Phase, ymin = q.025, ymax = q.975, color = variable, group = variable),
    width = .25,
    size = I(.8),
    position = pd,
    show.legend=F
  ) +
  geom_point(
    aes(x = Phase, y = avg, fill = variable, group = variable),
    color=I("black"),
    shape=I(21),
    size = I(3),
    position = pd
  ) +
  facet_wrap(~type, nrow = 1)+
  ylab("Total Cases -\nVariant Attributable Cases")+
  xlab("")+
  geom_hline(aes(yintercept = 0), linetype = I(2))+
  scale_color_manual(values = q4dfcolors$color, drop = FALSE, name="") +
  scale_fill_manual(values = q4dfcolors$color, drop = FALSE, name="") +
  theme(legend.position = "bottom")





## make example figure
results$Phase <- factor(as.factor(results$Phase), levels = c("Start","Impending Rise","Rise","Impending Fall","Fall"))
phase_ex <- qplot(last_obs_date, persist, data = subset(results, Phase != "Start" & series_id == "unitedstates_alabama" & model =="M(0)"), fill = Phase, shape = I(21), color = I("black"), size = I(2),show.legend = T)+
  scale_fill_manual(values = cols, name = "")+
  xlab("")+
  ylab("Total Cases")+
  ggtitle("Alabama")+
  ylim(c(-100, 1.1*max(subset(results, Phase != "Start" & series_id == "unitedstates_alabama" & model =="M(0)")$persist)))+
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")

# ## save the plot
setEPS()
postscript(paste0(figpath,"phase_of_outbreak.eps"),width = 8, height = 5.5)
   (phase_ex / boxplot_plot) + plot_layout(heights = c(1.5, 2)) 
dev.off()

 
