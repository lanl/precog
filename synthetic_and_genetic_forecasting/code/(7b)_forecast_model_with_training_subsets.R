# ============================================================
# Forecast script compatible with NEW training script:
# - normalization: divide by max(y_in) (cases >= 0)
# - model outputs residual quantiles that are:
#     * monotone across quantiles
#     * satisfy residual lower bound r >= -z_last
#   => implied z_out_hat >= 0 automatically
# - inverse transform:
#     y_hat = z_out_hat * scale
#     where scale = max(y_in) if >0 else 0 (and z_out_hat == y_out when scale==0 case)
# ============================================================

suppressPackageStartupMessages({
  library(torch)
  library(data.table)
  library(lubridate)
  library(ggplot2)
  library(foreach)
  library(doParallel)
})

# ---------------------------
# define paths
# ---------------------------
datapath <- paste0(here::here("synthetic_and_genetic_forecasting", "data_clean"), "/") 
modelpath <- paste0(here::here("synthetic_and_genetic_forecasting", "trained_models_subsets"), "/") 
savepath <- paste0(here::here("synthetic_and_genetic_forecasting", "output"), "/") 

# ---------------------------
# Read in test data
# ---------------------------
test_covid <- readRDS(paste0(datapath,"test_covid.RDS"))

# ---------------------------
# Forecast evaluation dates
# ---------------------------
last_obs_vec <- seq(ymd("2020-06-01"), ymd("2022-12-31"), by = "7 days")

# ---------------------------
# Hyperparameters
# ---------------------------
Nmcmc <- 100000L  # for CovidHub baseline + VAC bottom-up Monte Carlo

# ---------------------------
# Small utilities
# ---------------------------
`%||%` <- function(a, b) if (!is.null(a)) a else b

.assert_cols <- function(dt, cols, where = "data") {
  miss <- setdiff(cols, names(dt))
  if (length(miss)) stop(sprintf("Missing columns in %s: %s", where, paste(miss, collapse = ", ")))
}

# Build input exactly like UPDATED training:
# - take last min(n, C) observations as y_in
# - scale = max(y_in) (fallback 0)
# - z_in = y_in/scale (or identity if scale==0)
# - left-pad to length C using first available z (or 0)
.make_context_x <- function(y, C) {
  C <- as.integer(C)
  if (!is.finite(C) || C < 1L) stop("C must be a positive integer.")
  
  y <- as.numeric(y)
  y <- y[is.finite(y)]
  n <- length(y)
  if (n < 1L) stop("Empty 'cases' after filtering non-finite values.")
  
  y_in <- if (n >= C) y[(n - C + 1L):n] else y
  
  scale <- max(y_in, na.rm = TRUE)
  if (!is.finite(scale)) scale <- 0
  
  z <- if (scale > 0) (y_in / scale) else y_in
  
  if (length(z) < C) {
    pad_len <- C - length(z)
    pad_val <- if (length(z) >= 1L && is.finite(z[1])) z[1] else 0
    z <- c(rep(pad_val, pad_len), z)
  }
  
  list(x = z, z_last = z[C], scale = scale)
}

# Future time stamps from a time vector (numeric/integer/Date/POSIXct).
.make_future_time <- function(time_vec, H) {
  H <- as.integer(H)
  if (H < 1L) return(rep(NA, 0))
  
  t <- time_vec
  t <- t[!is.na(t)]
  if (!length(t)) return(rep(NA, H))
  t <- t[order(t)]
  last_t <- t[length(t)]
  
  if (length(t) >= 2) {
    if (inherits(t, "Date")) {
      diffs <- as.numeric(diff(t))
      step <- stats::median(diffs[is.finite(diffs) & diffs > 0], na.rm = TRUE)
      if (!is.finite(step) || step <= 0) step <- 7
      return(last_t + seq_len(H) * step)
    }
    if (inherits(t, "POSIXct") || inherits(t, "POSIXt")) {
      diffs <- as.numeric(diff(t))
      step <- stats::median(diffs[is.finite(diffs) & diffs > 0], na.rm = TRUE)
      if (!is.finite(step) || step <= 0) step <- 86400
      return(last_t + seq_len(H) * step)
    }
    if (is.numeric(t) || is.integer(t)) {
      diffs <- as.numeric(diff(t))
      step <- stats::median(diffs[is.finite(diffs) & diffs > 0], na.rm = TRUE)
      if (!is.finite(step) || step <= 0) step <- 1
      return(last_t + seq_len(H) * step)
    }
  }
  
  rep(NA, H)
}

# Split each series into:
#   le: rows with date <= ref_date
#   gt: next H rows after that (by time index), with date > ref_date, and time within (max_time+H)
split_list_by_date_with_horizon <- function(x, ref_date, H) {
  if (!is.list(x)) stop("x must be a list of data.tables/data.frames.")
  H <- as.integer(H)
  if (!is.finite(H) || H < 1L) stop("H must be a positive integer.")
  ref_date <- as.Date(ref_date)
  if (is.na(ref_date)) stop("ref_date could not be converted to Date.")
  
  le <- vector("list", length(x))
  gt <- vector("list", length(x))
  
  for (i in seq_along(x)) {
    dt <- as.data.table(x[[i]])
    .assert_cols(dt, c("date", "time"), where = sprintf("x[[%d]]", i))
    
    d <- as.Date(dt[["date"]])
    t <- dt[["time"]]
    
    ok <- !is.na(d) & is.finite(t)
    dt_ok <- dt[ok]
    
    le_i <- dt_ok[d[ok] <= ref_date]
    le[[i]] <- le_i
    
    if (nrow(le_i) == 0L) {
      gt[[i]] <- dt_ok[0]
      next
    }
    
    max_time <- max(le_i$time, na.rm = TRUE)
    if (!is.finite(max_time)) {
      gt[[i]] <- dt_ok[0]
      next
    }
    
    gt[[i]] <- dt_ok[d[ok] > ref_date & dt_ok$time > max_time & dt_ok$time <= (max_time + H)]
  }
  
  list(le = le, gt = gt)
}

# Attach future truth (and any columns present) to a forecast table by "time".
# Also adds series_id, last_obs, ref_date, train_mod, model_input.
.attach_meta <- function(fcst_dt, past_dt, future_dt, ref_date, train_mod, model_input) {
  fcst_dt <- as.data.table(fcst_dt)
  past_dt <- as.data.table(past_dt)
  future_dt <- as.data.table(future_dt)
  
  # ensure time exists
  .assert_cols(fcst_dt, c("time", "horizon"), where = "forecast")
  .assert_cols(past_dt, c("time", "cases", "series_id"), where = "past_dt")
  
  # last observed value
  tmax <- max(past_dt$time, na.rm = TRUE)
  last_val <- past_dt[time == tmax, cases][1]
  if (!is.finite(last_val)) last_val <- NA_real_
  
  # merge in truth/future columns when available
  if (nrow(future_dt) > 0L && "time" %in% names(future_dt)) {
    setkey(future_dt, time)
    setkey(fcst_dt, time)
    fcst_dt <- future_dt[fcst_dt]  # keep fcst rows, add future cols
  }
  
  fcst_dt[, series_id := past_dt$series_id[1]]
  fcst_dt[, last_obs := last_val]
  fcst_dt[, ref_date := as.IDate(ref_date)]
  fcst_dt[, train_mod := train_mod]
  fcst_dt[, model_input := model_input]
  
  # standardize truth column name for downstream: keep "cases" from future if present
  # (your downstream code expects "cases" to exist as y_true)
  fcst_dt[]
}

# Build a list of {time, cases} series for each var_attr_* column in a dt
make_var_attr_list <- function(dt, time_col = "time", pattern = "^var_attr") {
  dt <- as.data.table(dt)
  .assert_cols(dt, time_col, where = "dt")
  var_cols <- grep(pattern, names(dt), value = TRUE)
  if (!length(var_cols)) stop(sprintf("No columns matched pattern '%s'.", pattern))
  
  out <- lapply(var_cols, function(v) {
    x <- dt[, .(time = get(time_col), cases = get(v))]
    x
  })
  names(out) <- var_cols
  out
}

# ---------------------------
# Torch model definition (must match training)
# ---------------------------
positional_embedding <- nn_module(
  initialize = function(max_len, d_model) {
    self$max_len <- max_len
    self$pos_emb <- nn_embedding(num_embeddings = max_len, embedding_dim = d_model)
  },
  forward = function(x) {
    T <- x$size(2)
    pos <- torch_tensor(seq_len(T), device = x$device, dtype = torch_long())
    if (T > self$max_len) pos <- pos$clamp(min = 1, max = self$max_len)
    pe <- self$pos_emb(pos)$unsqueeze(1)  # [1,T,D]
    x + pe
  }
)

TransformerForecaster <- nn_module(
  initialize = function(
    context_length,
    horizon,
    n_quantiles,
    d_model = 256,
    n_head = 4,
    num_layers = 2,
    dim_feedforward = 512,
    dropout = 0
  ) {
    self$C <- context_length
    self$H <- horizon
    self$Q <- n_quantiles
    
    self$in_proj <- nn_linear(1, d_model)
    self$pos <- positional_embedding(max_len = context_length, d_model = d_model)
    
    enc_layer <- nn_transformer_encoder_layer(
      d_model = d_model,
      nhead = n_head,
      dim_feedforward = dim_feedforward,
      dropout = dropout,
      batch_first = TRUE
    )
    self$encoder <- nn_transformer_encoder(enc_layer, num_layers = num_layers)
    
    self$out_proj <- nn_linear(d_model, horizon * n_quantiles)
  },
  
  forward = function(x) {
    B <- x$size(1)
    
    z_last <- x[, self$C]
    base <- (-z_last)$view(c(B, 1, 1))
    
    x3 <- x$unsqueeze(3)
    h <- self$in_proj(x3)
    h <- self$pos(h)
    h <- self$encoder(h)
    
    last <- h[, self$C, ]
    raw <- self$out_proj(last)$view(c(B, self$H, self$Q))
    
    if (self$Q == 1) return(base + nnf_softplus(raw))
    
    first_raw <- raw[, , 1]$unsqueeze(3)
    first <- base + nnf_softplus(first_raw)
    
    deltas <- raw[, , 2:self$Q]
    deltas_pos <- nnf_softplus(deltas)
    cum <- torch_cumsum(deltas_pos, dim = 3)
    rest <- first + cum
    
    torch_cat(list(first, rest), dim = 3)
  }
)

# Infer transformer context length from checkpoint (safest way to avoid mismatches)
.infer_context_length_from_state <- function(state) {
  key <- "pos.pos_emb.weight"
  if (!key %in% names(state)) stop("Checkpoint missing 'pos.pos_emb.weight'; can't infer context_length.")
  as.integer(state[[key]]$size(1))
}

# ---------------------------
# Forecast with torch model (vectorized over series list)
# ---------------------------
forecast_mod <- function(model, cfg, test, batch_size = 2048L, enforce_monotone = TRUE) {
  if (!is.list(test) || length(test) == 0L) stop("test must be a non-empty list.")
  test <- lapply(test, as.data.table)
  
  H <- as.integer(cfg$horizon %||% cfg$H %||% model$H)
  if (!is.finite(H) || H < 1L) stop("Invalid horizon in cfg/model.")
  
  quantile_vec <- cfg$quantile_vec
  quantile_vec <- sort(unique(as.numeric(quantile_vec)))
  Q <- length(quantile_vec)
  if (Q < 1L) stop("cfg$quantile_vec missing or empty.")
  
  C <- as.integer(model$C)
  
  # device from model params
  params <- model$parameters
  device <- params[[1]]$device
  
  n_series <- length(test)
  x_mat <- matrix(0, nrow = n_series, ncol = C)
  z_last_vec <- numeric(n_series)
  scale_vec  <- numeric(n_series)
  time_next_list <- vector("list", n_series)
  
  for (i in seq_len(n_series)) {
    dt <- test[[i]]
    .assert_cols(dt, c("time", "cases"), where = sprintf("test[[%d]]", i))
    
    ord <- order(dt$time)
    y <- dt$cases[ord]
    t <- dt$time[ord]
    
    ok <- is.finite(y)
    y <- y[ok]
    t <- t[ok]
    
    if (!length(y)) stop(sprintf("Series %d has no finite cases.", i))
    
    ctx <- .make_context_x(y, C = C)
    x_mat[i, ] <- ctx$x
    z_last_vec[i] <- ctx$z_last
    scale_vec[i]  <- ctx$scale
    time_next_list[[i]] <- .make_future_time(t, H)
  }
  
  preds_array <- array(NA_real_, dim = c(n_series, H, Q))
  
  model$eval()
  with_no_grad({
    idx <- 1L
    while (idx <= n_series) {
      end_idx <- min(n_series, idx + as.integer(batch_size) - 1L)
      b <- end_idx - idx + 1L
      
      x <- torch_tensor(
        x_mat[idx:end_idx, , drop = FALSE],
        dtype = torch_float(),
        device = device
      )  # [b, C]
      
      r_hat <- model(x)                   # [b, H, Q]
      r_hat_np <- as.array(r_hat$to(device = "cpu"))
      
      for (j in seq_len(b)) {
        ii <- idx + j - 1L
        z_last <- z_last_vec[ii]
        scale  <- scale_vec[ii]
        
        z_out_hat <- z_last + r_hat_np[j, , ]  # [H, Q]
        z_out_hat <- pmax(z_out_hat, 0)
        
        y_hat <- if (scale > 0) z_out_hat * scale else matrix(0, nrow = H, ncol = Q)
        
        if (isTRUE(enforce_monotone) && Q > 1L) {
          for (h in seq_len(H)) y_hat[h, ] <- cummax(y_hat[h, ])
        }
        
        preds_array[ii, , ] <- y_hat
      }
      
      idx <- end_idx + 1L
    }
  })
  
  q_names <- paste0("q", formatC(quantile_vec, format = "f", digits = 4))
  out <- vector("list", n_series)
  for (i in seq_len(n_series)) {
    dt_out <- data.table(horizon = seq_len(H), time = time_next_list[[i]])
    for (q in seq_len(Q)) dt_out[[q_names[q]]] <- preds_array[i, , q]
    out[[i]] <- dt_out
  }
  
  out
}

# ---------------------------
# CovidHub-style baseline (quantile forecast of random-walk diffs)
# ---------------------------
covidhub_baseline <- function(y, quantiles, h, N = 100000L, seed = NULL, force_median = TRUE) {
  y <- as.numeric(y)
  y <- y[is.finite(y)]
  if (length(y) < 1L) stop("y has no finite values.")
  
  h <- as.integer(h)
  if (!is.finite(h) || h < 1L) stop("h must be positive integer.")
  
  N <- as.integer(N)
  if (!is.finite(N) || N < 1L) stop("N must be positive integer.")
  
  quantiles <- sort(unique(as.numeric(quantiles)))
  if (any(!is.finite(quantiles)) || any(quantiles <= 0) || any(quantiles >= 1)) {
    stop("quantiles must be in (0,1).")
  }
  
  if (!is.null(seed)) set.seed(seed)
  
  y_last <- y[length(y)]
  d <- if (length(y) >= 2L) diff(y) else numeric(0)
  d <- d[is.finite(d)]
  pool <- c(d, -d)
  if (!length(pool)) pool <- 0
  
  v <- sort(pool)
  n <- length(v)
  p <- if (n == 1L) 0.5 else (0:(n - 1L)) / (n - 1L)
  
  u <- stats::runif(N * h)
  draws <- stats::approx(x = p, y = v, xout = u, rule = 2, ties = "ordered")$y
  diff_mat <- matrix(draws, nrow = N, ncol = h)
  
  samp_mat <- diff_mat
  if (h >= 2L) {
    for (k in 2:h) samp_mat[, k] <- samp_mat[, k - 1L] + diff_mat[, k]
  }
  samp_mat <- y_last + samp_mat
  samp_mat <- pmax(samp_mat, 0)
  
  qmat <- apply(samp_mat, 2, stats::quantile, probs = quantiles, names = FALSE, type = 7)
  qmat <- t(qmat) # h x Q
  
  if (isTRUE(force_median)) {
    mid <- which(abs(quantiles - 0.5) < 1e-12)
    if (length(mid) == 1L) qmat[, mid] <- y_last
  }
  
  out <- data.table(horizon = seq_len(h))
  q_names <- paste0("q", formatC(quantiles, format = "f", digits = 4))
  for (j in seq_along(q_names)) out[[q_names[j]]] <- qmat[, j]
  out
}







# ---------------------------
# Bottom-up aggregation from many quantile tables (VACs)
# ---------------------------
mc_sum_from_quantile_tables <- function(fcst_list, N = 5000L, seed = NULL) {
  if (!is.list(fcst_list) || length(fcst_list) < 1L) stop("fcst_list must be a non-empty list.")
  N <- as.integer(N)
  if (!is.finite(N) || N < 1L) stop("N must be positive integer.")
  if (!is.null(seed)) set.seed(seed)
  
  fcst_list <- lapply(fcst_list, as.data.table)
  
  dt0 <- copy(fcst_list[[1]])
  base_cols <- intersect(c("horizon", "time"), names(dt0))
  q_cols <- setdiff(names(dt0), base_cols)
  if (length(q_cols) < 2L) stop("Need >=2 q-columns (e.g., q0.0050 ... q0.9950).")
  
  probs <- suppressWarnings(as.numeric(sub("^q", "", q_cols)))
  if (any(!is.finite(probs)) || any(probs <= 0) || any(probs >= 1)) {
    stop("Couldn't parse probs from q-cols. Expected names like 'q0.0250'.")
  }
  
  ordp <- order(probs)
  probs <- probs[ordp]
  q_cols_ord <- q_cols[ordp]
  
  R <- nrow(dt0)
  M <- length(fcst_list)
  
  comp_mat_list <- lapply(fcst_list, function(dt) as.matrix(dt[, ..q_cols_ord]))
  
  sample_from_quantile_fn <- function(qvals, u) {
    qv <- as.numeric(qvals)
    qv[!is.finite(qv)] <- NA_real_
    if (all(is.na(qv))) return(rep(NA_real_, length(u)))
    
    if (anyNA(qv)) {
      ok <- which(!is.na(qv))
      if (length(ok) == 1L) {
        qv <- rep(qv[ok], length(qv))
      } else {
        qv <- approx(x = probs[ok], y = qv[ok], xout = probs, rule = 2, ties = "ordered")$y
      }
    }
    
    qv <- cummax(qv)
    approx(x = probs, y = qv, xout = u, rule = 2, ties = "ordered")$y
  }
  
  out_mat <- matrix(NA_real_, nrow = R, ncol = length(q_cols_ord))
  
  for (r in seq_len(R)) {
    sum_draws <- numeric(N)
    any_na_component <- FALSE
    
    for (m in seq_len(M)) {
      u <- stats::runif(N)
      draws_m <- sample_from_quantile_fn(comp_mat_list[[m]][r, ], u)
      if (all(is.na(draws_m))) any_na_component <- TRUE
      draws_m[is.na(draws_m)] <- 0
      sum_draws <- sum_draws + draws_m
    }
    
    out_mat[r, ] <- if (any_na_component) NA_real_ else
      as.numeric(stats::quantile(sum_draws, probs = probs, type = 7, names = FALSE))
  }
  
  out <- dt0[, ..base_cols]
  out_q <- as.data.table(out_mat)
  setnames(out_q, q_cols_ord)
  out_q <- out_q[, ..q_cols]  # restore original q col ordering
  cbind(out, out_q)
}




# fcst_list: list of data.tables like your example
# Nmcmc: number of draws to sample per row
# Returns: list of data.tables, each with Nmcmc rows and nrow(fcst_dt) columns
sample_fcst_list_draws <- function(fcst_list, Nmcmc, seed = NULL) {
  stopifnot(is.list(fcst_list), length(fcst_list) > 0)
  stopifnot(length(Nmcmc) == 1, is.numeric(Nmcmc), Nmcmc >= 1)
  Nmcmc <- as.integer(Nmcmc)
  
  if (!is.null(seed)) set.seed(seed)
  
  lapply(fcst_list, function(fcst_dt) {
    stopifnot(is.data.table(fcst_dt))
    
    # quantile columns are those starting with "q"
    qcols <- grep("^q", names(fcst_dt), value = TRUE)
    if (length(qcols) < 2) stop("Need at least 2 quantile columns (names like q0.0500).")
    
    # probabilities extracted from colnames, e.g. "q0.0050" -> 0.005
    probs <- as.numeric(sub("^q", "", qcols))
    if (any(!is.finite(probs))) stop("Could not parse probabilities from quantile column names.")
    o <- order(probs)
    probs <- probs[o]
    qcols <- qcols[o]
    
    nr <- nrow(fcst_dt)
    
    # nice column names for each row (one column per horizon/time row)
    if (all(c("horizon", "time") %in% names(fcst_dt))) {
      out_names <- paste0("h", fcst_dt[["horizon"]], "_t", fcst_dt[["time"]])
    } else {
      out_names <- paste0("row", seq_len(nr))
    }
    
    # sample uniform values once; reuse across rows (each row still gets its own interpolated values)
    u <- runif(Nmcmc)
    
    out_mat <- matrix(NA_real_, nrow = Nmcmc, ncol = nr)
    
    # for each row, treat provided quantiles as a piecewise-linear quantile function Q(p)
    for (r in seq_len(nr)) {
      yq <- as.numeric(fcst_dt[r, ..qcols])
      # enforce nondecreasing to avoid quantile crossing issues
      yq <- cummax(yq)
      
      out_mat[, r] <- stats::approx(
        x = probs, y = yq, xout = u,
        method = "linear",
        rule = 2,            # flat extrapolation beyond ends
        ties = "ordered"
      )$y
    }
    
    out <- as.data.table(out_mat)
    setnames(out, out_names)
    out
  })
}






# ======================================================================================
# 1) Forecast baseline (non-torch) on total cases
# ======================================================================================
message("Forecasting baseline...")

# For baseline we just need quantiles + horizon; reuse a real cfg for those values.
cfg_baseline <- readRDS(file.path(modelpath, "cfg_real_2500.rds"))

fcst_output_baseline <- NULL

for (k in seq_along(last_obs_vec)) {
  ref_date <- last_obs_vec[k]
  split_test <- split_list_by_date_with_horizon(test_covid, ref_date = ref_date, H = cfg_baseline$horizon)
  past_list <- split_test$le
  fut_list  <- split_test$gt
  
  for (i in seq_along(past_list)) {
    past_dt <- past_list[[i]]
    fut_dt  <- fut_list[[i]]
    
    if (nrow(past_dt) == 0L) next
    
    fcst <- covidhub_baseline(
      y = past_dt$cases,
      quantiles = cfg_baseline$quantile_vec,
      h = cfg_baseline$horizon,
      N = Nmcmc,
      seed = 1,
      force_median = TRUE
    )
    
    # attach time: prefer actual future times; otherwise infer
    if (nrow(fut_dt) > 0L) {
      fcst[, time := fut_dt$time[match(horizon, seq_len(nrow(fut_dt)))]]
    } else {
      fcst[, time := .make_future_time(past_dt$time, cfg_baseline$horizon)]
    }
    
    fcst_dt <- .attach_meta(
      fcst_dt = fcst,
      past_dt = past_dt,
      future_dt = fut_dt,
      ref_date = ref_date,
      train_mod = "mod_baseline",
      model_input = "tc"
    )
    
    fcst_output_baseline <- rbind(fcst_output_baseline, fcst_dt, fill = TRUE)
  }
  
  if (k %% 25 == 0) message(sprintf("  baseline: %d / %d", k, length(last_obs_vec)))
}


# ======================================================================================
# 2) Forecast torch models with total cases input (tc): PARALLELIZED like #3
# ======================================================================================
message("Forecasting torch models (tc, parallel)...")

n_workers <- max(1L, parallel::detectCores() - 1L)
cl <- parallel::makeCluster(n_workers)
doParallel::registerDoParallel(cl)

# Put required definitions on workers (nn_modules must exist in worker env)
parallel::clusterEvalQ(cl, {
  library(torch)
  library(data.table)
  
  `%||%` <- function(a, b) if (!is.null(a)) a else b
  
  positional_embedding <- nn_module(
    initialize = function(max_len, d_model) {
      self$max_len <- max_len
      self$pos_emb <- nn_embedding(num_embeddings = max_len, embedding_dim = d_model)
    },
    forward = function(x) {
      T <- x$size(2)
      pos <- torch_tensor(seq_len(T), device = x$device, dtype = torch_long())
      if (T > self$max_len) pos <- pos$clamp(min = 1, max = self$max_len)
      pe <- self$pos_emb(pos)$unsqueeze(1)
      x + pe
    }
  )
  
  TransformerForecaster <- nn_module(
    initialize = function(context_length, horizon, n_quantiles,
                          d_model = 256, n_head = 4, num_layers = 2,
                          dim_feedforward = 512, dropout = 0) {
      self$C <- context_length
      self$H <- horizon
      self$Q <- n_quantiles
      
      self$in_proj <- nn_linear(1, d_model)
      self$pos <- positional_embedding(max_len = context_length, d_model = d_model)
      
      enc_layer <- nn_transformer_encoder_layer(
        d_model = d_model,
        nhead = n_head,
        dim_feedforward = dim_feedforward,
        dropout = dropout,
        batch_first = TRUE
      )
      self$encoder <- nn_transformer_encoder(enc_layer, num_layers = num_layers)
      self$out_proj <- nn_linear(d_model, horizon * n_quantiles)
    },
    forward = function(x) {
      B <- x$size(1)
      z_last <- x[, self$C]
      base <- (-z_last)$view(c(B, 1, 1))
      
      x3 <- x$unsqueeze(3)
      h <- self$in_proj(x3)
      h <- self$pos(h)
      h <- self$encoder(h)
      
      last <- h[, self$C, ]
      raw <- self$out_proj(last)$view(c(B, self$H, self$Q))
      
      if (self$Q == 1) return(base + nnf_softplus(raw))
      
      first_raw <- raw[, , 1]$unsqueeze(3)
      first <- base + nnf_softplus(first_raw)
      
      deltas <- raw[, , 2:self$Q]
      deltas_pos <- nnf_softplus(deltas)
      cum <- torch_cumsum(deltas_pos, dim = 3)
      rest <- first + cum
      
      torch_cat(list(first, rest), dim = 3)
    }
  )
  
  NULL
})


# .tc_models <- data.table(
#   cfg_name = c("cfg_real.rds","cfg_syn_tc.rds","cfg_syn_vac.rds","cfg_all.rds"),   # cfg files inside modelpath
#   ckpt     = c("mod_real.pt","mod_syn_tc.pt","mod_syn_vac.pt","mod_all.pt")     # checkpoint files inside modelpath
# )

.tc_models <- data.table(
  cfg_name = sort(list.files(modelpath, pattern=".rds")),   # cfg files inside modelpath
  ckpt     = sort(list.files(modelpath, pattern=".pt"))     # checkpoint files inside modelpath
)

parallel::clusterExport(
  cl,
  varlist = c(
    "modelpath", "test_covid", "last_obs_vec",
    ".tc_models",
    "split_list_by_date_with_horizon",
    "forecast_mod", ".attach_meta",
    ".infer_context_length_from_state"
  ),
  envir = environment()
)

id_setup <- data.table(
  j = rep(seq_len(nrow(.tc_models)), each = length(last_obs_vec)),
  k = rep(seq_along(last_obs_vec), times = nrow(.tc_models))
)

fcst_output_tc <- foreach(
  idx = seq_len(nrow(id_setup)),
  .combine = "rbind",
  .multicombine = TRUE,
  .inorder = FALSE,
  .packages = c("torch", "data.table", "lubridate")
) %dopar% {
  
  j <- id_setup$j[idx]
  k <- id_setup$k[idx]
  
  cfg_name   <- .tc_models$cfg_name[j]
  model_ckpt <- .tc_models$ckpt[j]
  
  cfg <- readRDS(file.path(modelpath, cfg_name))
  state <- torch::torch_load(file.path(modelpath, model_ckpt))
  
  # infer context length from checkpoint to guarantee compatibility
  C_state <- .infer_context_length_from_state(state)
  
  # NOTE: If 1 GPU only, consider forcing CPU here:
  # device <- "cpu"
  device <- if (torch::cuda_is_available()) "cuda" else "cpu"
  
  model <- TransformerForecaster(
    context_length  = C_state,
    horizon         = cfg$horizon,
    n_quantiles     = length(cfg$quantile_vec),
    d_model         = cfg$d_model,
    n_head          = cfg$n_head,
    num_layers      = cfg$num_layers,
    dim_feedforward = cfg$dim_feedforward,
    dropout         = cfg$dropout
  )$to(device = torch_device(device))
  
  model$load_state_dict(state)
  model$eval()
  
  ref_date <- last_obs_vec[k]
  split_test <- split_list_by_date_with_horizon(test_covid, ref_date = ref_date, H = cfg$horizon)
  past_list <- split_test$le
  fut_list  <- split_test$gt
  
  # Forecast all series at once (tc input)
  fcst_list <- forecast_mod(model = model, cfg = cfg, test = past_list, batch_size = 2048L)
  
  # Attach truth/meta per series
  out_worker <- NULL
  for (i in seq_along(fcst_list)) {
    if (nrow(past_list[[i]]) == 0L) next
    fcst_i <- fcst_list[[i]]
    
    # attach actual future time stamps if present
    if (nrow(fut_list[[i]]) > 0L) {
      fcst_i[, time := fut_list[[i]]$time[match(horizon, seq_len(nrow(fut_list[[i]])))]]
    }
    
    fcst_i <- .attach_meta(
      fcst_dt = fcst_i,
      past_dt = past_list[[i]],
      future_dt = fut_list[[i]],
      ref_date = ref_date,
      train_mod = gsub("\\.pt$", "", model_ckpt),
      model_input = "tc"
    )
    
    out_worker <- rbind(out_worker, fcst_i, fill = TRUE)
  }
  out_worker
}

parallel::stopCluster(cl)

message("Done forecasting torch models (tc, parallel).")



print(head(fcst_output_tc))




# ======================================================================================
# 3) Forecast torch models on VACs (bottom-up): parallelized
# ======================================================================================
message("Forecasting torch models (vac bottom-up, parallel)...")

n_workers <- max(1L, parallel::detectCores() - 1L)
cl <- parallel::makeCluster(n_workers)
doParallel::registerDoParallel(cl)

# Put required definitions on workers (nn_modules must exist in worker env)
parallel::clusterEvalQ(cl, {
  library(torch)
  library(data.table)

  `%||%` <- function(a, b) if (!is.null(a)) a else b

  positional_embedding <- nn_module(
    initialize = function(max_len, d_model) {
      self$max_len <- max_len
      self$pos_emb <- nn_embedding(num_embeddings = max_len, embedding_dim = d_model)
    },
    forward = function(x) {
      T <- x$size(2)
      pos <- torch_tensor(seq_len(T), device = x$device, dtype = torch_long())
      if (T > self$max_len) pos <- pos$clamp(min = 1, max = self$max_len)
      pe <- self$pos_emb(pos)$unsqueeze(1)
      x + pe
    }
  )

  TransformerForecaster <- nn_module(
    initialize = function(context_length, horizon, n_quantiles,
                          d_model = 256, n_head = 4, num_layers = 2,
                          dim_feedforward = 512, dropout = 0) {
      self$C <- context_length
      self$H <- horizon
      self$Q <- n_quantiles

      self$in_proj <- nn_linear(1, d_model)
      self$pos <- positional_embedding(max_len = context_length, d_model = d_model)

      enc_layer <- nn_transformer_encoder_layer(
        d_model = d_model,
        nhead = n_head,
        dim_feedforward = dim_feedforward,
        dropout = dropout,
        batch_first = TRUE
      )
      self$encoder <- nn_transformer_encoder(enc_layer, num_layers = num_layers)
      self$out_proj <- nn_linear(d_model, horizon * n_quantiles)
    },
    forward = function(x) {
      B <- x$size(1)
      z_last <- x[, self$C]
      base <- (-z_last)$view(c(B, 1, 1))

      x3 <- x$unsqueeze(3)
      h <- self$in_proj(x3)
      h <- self$pos(h)
      h <- self$encoder(h)

      last <- h[, self$C, ]
      raw <- self$out_proj(last)$view(c(B, self$H, self$Q))

      if (self$Q == 1) return(base + nnf_softplus(raw))

      first_raw <- raw[, , 1]$unsqueeze(3)
      first <- base + nnf_softplus(first_raw)

      deltas <- raw[, , 2:self$Q]
      deltas_pos <- nnf_softplus(deltas)
      cum <- torch_cumsum(deltas_pos, dim = 3)
      rest <- first + cum

      torch_cat(list(first, rest), dim = 3)
    }
  )

  NULL
})

parallel::clusterExport(
  cl,
  varlist = c(
    "modelpath", "test_covid", "last_obs_vec", "Nmcmc",
    "split_list_by_date_with_horizon", "make_var_attr_list",
    "forecast_mod", "mc_sum_from_quantile_tables", ".attach_meta",
    ".infer_context_length_from_state"
  ),
  envir = environment()
)

id_setup <- data.table(
  j = rep(seq_len(nrow(.tc_models)), each = length(last_obs_vec)),
  k = rep(seq_along(last_obs_vec), times = nrow(.tc_models))
)

fcst_output_vac <- foreach(
  idx = seq_len(nrow(id_setup)),
  .combine = "rbind",
  .multicombine = TRUE,
  .inorder = FALSE,
  .packages = c("torch", "data.table", "lubridate")
) %dopar% {

  j <- id_setup$j[idx]
  k <- id_setup$k[idx]

  cfg_name   <- .tc_models$cfg_name[j]
  model_ckpt <- .tc_models$ckpt[j]

  cfg <- readRDS(file.path(modelpath, cfg_name))
  state <- torch::torch_load(file.path(modelpath, model_ckpt))
  C_state <- .infer_context_length_from_state(state)

  device <- if (torch::cuda_is_available()) "cuda" else "cpu"
  model <- TransformerForecaster(
    context_length  = C_state,
    horizon         = cfg$horizon,
    n_quantiles     = length(cfg$quantile_vec),
    d_model         = cfg$d_model,
    n_head          = cfg$n_head,
    num_layers      = cfg$num_layers,
    dim_feedforward = cfg$dim_feedforward,
    dropout         = cfg$dropout
  )$to(device = torch_device(device))

  model$load_state_dict(state)
  model$eval()

  ref_date <- last_obs_vec[k]
  split_test <- split_list_by_date_with_horizon(test_covid, ref_date = ref_date, H = cfg$horizon)
  past_list <- split_test$le
  fut_list  <- split_test$gt

  out_worker <- NULL

  for (i in seq_along(past_list)) {
    past_dt <- past_list[[i]]
    fut_dt  <- fut_list[[i]]
    if (nrow(past_dt) == 0L) next

    # Make per-variant series list (each is time/cases where cases = var_attr_*)
    vac_list <- make_var_attr_list(past_dt, time_col = "time", pattern = "^var_attr")

    # Forecast each variant series (one pass; returns list of q-tables)
    fcst_vacs <- forecast_mod(model = model, cfg = cfg, test = vac_list, batch_size = 2048L)

    # Monte Carlo sum across variants, then re-quantile
    fcst_mc <- mc_sum_from_quantile_tables(fcst_vacs, N = Nmcmc, seed = 1)

    # attach time from true future if present
    if (nrow(fut_dt) > 0L) {
      fcst_mc[, time := fut_dt$time[match(horizon, seq_len(nrow(fut_dt)))]]
    }

    fcst_dt <- .attach_meta(
      fcst_dt = fcst_mc,
      past_dt = past_dt,
      future_dt = fut_dt,
      ref_date = ref_date,
      train_mod = gsub("\\.pt$", "", model_ckpt),
      model_input = "vac"
    )

    out_worker <- rbind(out_worker, fcst_dt, fill = TRUE)
  }

  out_worker
}

parallel::stopCluster(cl)



# ======================================================================================
# 4) Concatenate and save (keep only 7 WIS quantiles)
# ======================================================================================
message("Combining + saving...")

fcst_output <- rbindlist(
  list(fcst_output_baseline, fcst_output_tc,  fcst_output_vac),
  use.names = TRUE,
  fill = TRUE
)

# Keep only the 7 quantiles used downstream for WIS
keep_q <- c(0.025, 0.10, 0.25, 0.50, 0.75, 0.90, 0.975)
keep_q_cols <- paste0("q", formatC(keep_q, format = "f", digits = 4))

# Ensure standard columns exist
needed_meta <- c("series_id", "time", "horizon", "ref_date", "last_obs", "train_mod", "model_input")
missing_meta <- setdiff(needed_meta, names(fcst_output))
if (length(missing_meta)) {
  stop("Final fcst_output is missing required columns: ", paste(missing_meta, collapse = ", "))
}

# Keep all non-quantile columns + the 7 quantiles
non_q_cols <- grep("^q0\\.", names(fcst_output), value = TRUE, invert = TRUE)
keep_cols <- unique(c(non_q_cols, keep_q_cols))
keep_cols <- intersect(keep_cols, names(fcst_output))

fcst_output_q <- fcst_output[, ..keep_cols]

## get model name
fcst_output_q$just_model <- sub("_[0-9]+$", "", fcst_output_q$train_mod)

## save these by just_model
dir.create(savepath, recursive = TRUE, showWarnings = FALSE)

unq_just_model <- unique(fcst_output_q$just_model)
for(i in 1:length(unq_just_model)){
  fwrite(
    subset(fcst_output_q, just_model == unq_just_model[i], select = setdiff(names(fcst_output_q),c("just_model",grep("var_",names(fcst_output_q), value=T)))),
    file = paste0(savepath,"forecasts_model_subsets_",unq_just_model[i],".csv"),
    row.names = FALSE)
}


