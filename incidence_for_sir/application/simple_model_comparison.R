library(deSolve)
library(ggplot2)
library(dplyr)
library(reshape2)
library(data.table)
library(plyr)
library(gridExtra)
library(tmvtnorm)
library(latex2exp)
library(this.path)
library(EpiModel)
library(tidyr)
library(purrr)
library(parallel)
setwd(paste0(this.path::here(), "/code"))
theme_set(theme_bw())

flu_seasons = c(2010, 2010) #, 2014, 2014
lengths = c(13, 22, 12, 22)
graphic_idx = 1
# -------------------------
# 1. Load or simulate data
# -------------------------
flu_season = flu_seasons[graphic_idx]
# Set forecast horizon
length_of_time_series = lengths[graphic_idx]

## define paths
## Note: to recreate this experiment, get the data from Osthus et al 2017 and
## place this data in a directory located at the path below.
datapath <- paste0(this.path::here(), "/data/")
jagspath <- paste0(this.path::here(), "/code/")

## load ILI and surveillance data
dfili <- read.csv(paste0(datapath,"EW08-2020_nat.csv"))
dfili <- subset(dfili, region == "nat", select=c("epi_year","epi_week","wili","epi_time","epi_season"))
dfsurveillance <- read.csv(paste0(datapath,"EW08-2020_who_combined_prior_to_2015_16.csv"))
dfsurveillance <- subset(dfsurveillance, select=c("year","week","percent_positive"), region == "National")

## merge ILI data with surveillance data
df <- merge(dfili, dfsurveillance, by.x=c("epi_year","epi_week"), by.y=c("year","week"))
df$iliplus <- (df$wili/100)*(df$percent_positive/100)

## focus on 2010, national wILI (Fig 6 of Osthus et. al. AoAS)
final_epi_time = 35
df_forecast_year <- subset(df, epi_season == flu_season & epi_time <= final_epi_time)
df_forecast_year <- df_forecast_year[order(df_forecast_year$epi_time),]

# --- Data (expects your df_forecast_year already in the environment) ---
# Keep one influenza season and the first 35 weeks of that season
season_year <- unique(df_forecast_year$epi_season)[1]
dat_season <- df_forecast_year %>%
  filter(epi_season == season_year) %>%
  arrange(epi_time) %>%
  filter(epi_time <= 35)

# Use iliplus as weekly incidence proportion (new infections per capita)
inc_obs <- dat_season$iliplus
time_obs <- dat_season$epi_time
stopifnot(length(inc_obs) >= length_of_time_series)

# Data slice (weeks 1..35 within season)
season_year <- unique(df_forecast_year$epi_season)[1]
dat_season <- df_forecast_year %>%
  filter(epi_season == season_year) %>%
  arrange(epi_time) %>%
  filter(epi_time <= 35)

inc_obs <- dat_season$iliplus
stopifnot(length(inc_obs) >= length_of_time_series)
inc_fit <- inc_obs[1:length_of_time_series]

# Deterministic incidence from dcm (per-capita, using N=1)
sir_dcm_incidence <- function(beta, gamma, i0, nsteps = 35) {
  i0 <- min(max(i0, 1e-10), 1 - 1e-10)
  par <- param.dcm(inf.prob = 1, act.rate = beta, rec.rate = gamma)
  ini <- init.dcm(s.num = 1 - i0, i.num = i0, r.num = 0)
  ctl <- control.dcm(type = "SIR", nsteps = nsteps, dt = 1)
  df  <- as.data.frame(dcm(par, ini, ctl))
  df$si.flow
}

# Reparameterize to enforce R0 > 1: beta = gamma * exp(eta), eta > 0
# Also estimate reporting fraction rho in (0,1]
log_eps <- 1e-12 # to avoid log(0)
obj <- function(theta) {
  # theta = (log_gamma, log_eta, logit_rho_raw)
  gamma <- exp(theta[1])
  eta   <- exp(theta[2])                # positive
  beta  <- gamma * eta                  # => R0 = beta/gamma = eta > 1
  rho   <- 1 / (1 + exp(-theta[3]))     # (0,1)
  
  # tie i0 to first observation so week-1 is on-scale
  i0 <- max(1e-10, min(1e-2, inc_fit[1] / (rho * beta)))
  
  inc_hat <- sir_dcm_incidence(beta, gamma, i0, nsteps = 35)[1:length_of_time_series]
  
  # fit on log-scale with a small floor
  ll <- (log(pmax(inc_fit, log_eps)) - log(pmax(rho * inc_hat, log_eps)))^2
  sum(ll)
}

theta0 <- c(log(1/3), log(1.5), qlogis(0.1)) # gamma ~ 1/3wk, R0 ≈ 1.5, rho ≈ 0.1
fit <- optim(theta0, obj, method = "Nelder-Mead",
             control = list(maxit = 4000, reltol = 1e-12))

gamma_hat <- exp(fit$par[1])
eta_hat   <- exp(fit$par[2])
beta_hat  <- gamma_hat * eta_hat
rho_hat   <- 1 / (1 + exp(-fit$par[3]))
i0_hat    <- max(1e-10, min(1e-2, inc_fit[1] / (rho_hat * beta_hat)))

cat(sprintf("Fitted (log-scale, R0>1): R0=%.3f, beta=%.5f, gamma=%.5f, rho=%.3f, i0=%.6g\n",
            beta_hat/gamma_hat, beta_hat, gamma_hat, rho_hat, i0_hat))

# Deterministic forecast (per-capita)
inc_det <- sir_dcm_incidence(beta_hat, gamma_hat, i0_hat, nsteps = 35)

# Compare the fit vs observed in first 10 weeks
compare_fit <- tibble(
  epi_time = 1:length_of_time_series,
  obs      = inc_fit,
  model    = rho_hat * inc_det[1:length_of_time_series]
)
print(compare_fit)
## Plot
ggplot(compare_fit, aes(x = epi_time)) +
  geom_line(aes(y = model), linewidth = 1) +
  geom_point(aes(y = obs), size = 2) +
  theme_bw()

# 1) deterministic run to get week-10 state
df_dcm <- {
  par <- param.dcm(inf.prob = 1, act.rate = beta_hat, rec.rate = gamma_hat)
  ini <- init.dcm(s.num = 1 - i0_hat, i.num = i0_hat, r.num = 0)
  ctl <- control.dcm(type = "SIR", nsteps = 35, dt = 1)
  as.data.frame(dcm(par, ini, ctl))
}
s10 <- df_dcm$s.num[length_of_time_series]; i10 <- df_dcm$i.num[length_of_time_series]; r10 <- df_dcm$r.num[length_of_time_series]

# 2) map to icm contact process: choose plausible contacts and compute per-contact prob
contacts_week <- 20
inf_prob <- min(0.5, beta_hat / contacts_week)

# 3) initialize counts at week 10 and simulate only weeks 11..35
Npop <- 1e5
s0 <- max(0L, round(s10 * Npop))
i0 <- max(1L, round(i10 * Npop))
r0 <- max(0L, Npop - s0 - i0)

par_icm <- param.icm(inf.prob = inf_prob, act.rate = contacts_week, rec.rate = gamma_hat)
ini_icm <- init.icm(s.num = s0, i.num = i0, r.num = r0)

steps_after <- 35 - length_of_time_series
ncores_use <- 10
ctl_icm <- control.icm(type = "SIR", nsteps = steps_after, nsims = 300,
                       ncores = ncores_use, verbose = FALSE)

set.seed(123)
mod_icm <- icm(par_icm, ini_icm, ctl_icm)

# Set your conditioning week and horizon explicitly
week_cond <- length_of_time_series
horizon   <- 35

# Pull out data frame and make the relative time index robust
df_icm <- as.data.frame(mod_icm)   # columns like: sim, time, s.num, i.num, num, r.num, si.flow, ir.flow

# Population used for per-capita scaling
Npop_used <- if ("num" %in% names(df_icm)) max(df_icm$num) else if (exists("Npop")) Npop else 1

# Build per-sim incidence, aligning weeks to 11..35 regardless of 'time' starting at 0 or 1
df_inc <- df_icm %>%
  dplyr::arrange(sim, time) %>%
  dplyr::group_by(sim) %>%
  dplyr::mutate(
    t_rel    = seq_along(time),        # 1..(horizon - week_cond) per sim
    epi_time = week_cond + t_rel  - 1     # => 11..35
  ) %>%
  dplyr::ungroup() %>%
  dplyr::transmute(sim, epi_time, incidence = si.flow / Npop_used) %>%
  subset(epi_time > length_of_time_series)

# Summarize predictive intervals by week (reported scale if you estimated rho_hat)
scale_to_reported <- function(x) {
  if (exists("rho_hat")) rho_hat * x else x
}

pred_summary <- df_inc %>%
  group_by(epi_time) %>%
  dplyr::summarize(
    pred_median = mean(scale_to_reported(incidence)),
    pred_lo80   = quantile(scale_to_reported(incidence), 0.10),
    pred_hi80   = quantile(scale_to_reported(incidence), 0.90),
    pred_lo95   = quantile(scale_to_reported(incidence), 0.025),
    pred_hi95   = quantile(scale_to_reported(incidence), 0.975),
    .groups = "drop"
  )

# Deterministic line on the plotting scale:
# prefer 'det_rep' if you already computed it; else use rho_hat*inc_det if rho_hat exists; else inc_det
y_det <- if (exists("det_rep")) {
  det_rep
} else if (exists("inc_det") && exists("rho_hat")) {
  rho_hat * inc_det
} else if (exists("inc_det")) {
  inc_det
} else {
  rep(NA_real_, horizon)
}

# Base frame for plotting
pred_df <- tibble(epi_time = 1:horizon, y_det = y_det) %>%
  dplyr::left_join(pred_summary, by = "epi_time")

# Reattach observed ILI+ for the same season (assumes df_forecast_year & season_year exist)
obs_df <- df_forecast_year %>%
  filter(epi_season == unique(df_forecast_year$epi_season)[1],
         epi_time <= horizon) %>%
  select(epi_time, iliplus)

pred_df <- pred_df %>%
  left_join(obs_df, by = "epi_time") %>%
  # Force weeks 1..week_cond to hug the deterministic/observed path
  mutate(
    pred_median = ifelse(epi_time <= week_cond, y_det, pred_median),
    pred_lo80   = ifelse(epi_time <= week_cond, y_det, pred_lo80),
    pred_hi80   = ifelse(epi_time <= week_cond, y_det, pred_hi80),
    pred_lo95   = ifelse(epi_time <= week_cond, y_det, pred_lo95),
    pred_hi95   = ifelse(epi_time <= week_cond, y_det, pred_hi95)
  )

# --- Sanity checks (optional) ---
# range(unique(df_inc$epi_time))     # should be (week_cond+1)..horizon
# dplyr::count(df_inc, epi_time)     # ~ nsims rows per week

# --- Plot ---
ggplot(pred_df, aes(x = epi_time)) +
  geom_ribbon(
    data = subset(pred_df, epi_time > week_cond),
    aes(ymin = pred_lo95, ymax = pred_hi95, fill = "95% PI"),
    alpha = 0.15
  ) +
  geom_ribbon(
    data = subset(pred_df, epi_time > week_cond),
    aes(ymin = pred_lo80, ymax = pred_hi80, fill = "80% PI"),
    alpha = 0.25
  ) +
  geom_line(
    data = subset(pred_df, epi_time > week_cond),
    aes(y = pred_median, linetype = "Median forecast"),
    linewidth = 0.9
  ) +
  geom_line(aes(y = y_det, linetype = "Deterministic fit"), linewidth = 1) +
  geom_point(aes(y = iliplus, shape = "Observed ILI+"), size = 2) +
  geom_vline(xintercept = week_cond + 0.5, linetype = "dashed") +
  annotate("text",
           x = week_cond + 0.5,
           y = max(pred_df$y_det, na.rm = TRUE),
           label = paste0("conditioning\n(after week", length_of_time_series,")"),
           vjust = -0.5, hjust = 0, size = 3) +
  scale_fill_manual(NULL, values = c("95% PI" = "grey70", "80% PI" = "grey50")) +
  scale_linetype_manual(NULL, values = c("Deterministic fit" = "solid",
                                         "Median forecast" = "dashed")) +
  scale_shape_manual(NULL, values = c("Observed ILI+" = 16)) +
  labs(x = "Epidemic week",
       y = "Weekly incidence (ILI+ per capita)",
       title = "Conditioned stochastic forecasts: exact match on weeks 1–13, uncertainty on weeks 14–35",
       subtitle = "Line = deterministic fit; dots = observed ILI+; ribbons = 80% and 95% predictive intervals") +
  theme_bw() +
  theme(legend.position = "top",
        legend.key.width = unit(1.6, "lines"))

