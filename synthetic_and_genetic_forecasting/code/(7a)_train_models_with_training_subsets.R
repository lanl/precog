# ============================================================
# Transformer time-series quantile forecaster (torch)
# - NORMALIZATION: divide by max(y_in)
#     z_in  = y_in / max(y_in) if max(y_in) > 0 else z_in = y_in
#     z_out = y_out / max(y_in) if max(y_in) > 0 else z_out = y_out
# - TARGET: residuals on normalized scale
#     r_out = z_out - z_last, where z_last = z_in[C]
# - CONSTRAINT (NEW): since y >= 0 => z_out >= 0, therefore
#     z_last + r_out >= 0  ==>  r_out >= -z_last
#   We enforce this LOWER BOUND for *every* predicted quantile by
#   parameterizing the model output as:
#     r_q1 = -z_last + softplus(raw_q1)
#     r_qk = r_q(k-1) + softplus(raw_qk)   for k>=2
#   This simultaneously guarantees:
#     (1) non-crossing quantiles (monotone in q)
#     (2) residual lower bound r >= -z_last  (=> z_out_hat >= 0)
# ============================================================


# ---------------------------
# define paths
# ---------------------------
datapath <- paste0(here::here("synthetic_and_genetic_forecasting", "data_clean"), "/") #"/Users/dosthus/Documents/DR24_PreCog/projects/four_quadrants/quadrants_complete_2-11-26/data_clean/"
modelpath <- paste0(here::here("synthetic_and_genetic_forecasting", "trained_models_subsets"), "/") #"/Users/dosthus/Documents/DR24_PreCog/projects/four_quadrants/quadrants_complete_2-11-26/trained_models_subsets/"


# ---------------------------
# load packages
# ---------------------------
suppressPackageStartupMessages({
  library(torch)
  library(data.table)
  library(reshape2)
})

# ---------------------------
# Read in training data
# ---------------------------
train_real <- readRDS(paste0(datapath,"train_real.RDS"))
train_syn_tc <- readRDS(paste0(datapath,"train_syn_tc.RDS"))
train_syn_vac <- readRDS(paste0(datapath,"train_syn_vac.RDS"))


# ---------------------------
# some training data summaries
# ---------------------------

sum_train <- NULL
sum_train <- rbind(sum_train,
                   data.frame(train_type = "real",
                              n_ts = length(train_real),
                              total_obs = sum(unlist(lapply(train_real, function(x){nrow(x)}))),
                              avg_length = mean(unlist(lapply(train_real, function(x){nrow(x)}))),
                              median_length = median(unlist(lapply(train_real, function(x){nrow(x)})))))

sum_train <- rbind(sum_train,
                   data.frame(train_type = "syn_tc",
                              n_ts = length(train_syn_tc),
                              total_obs = sum(unlist(lapply(train_syn_tc, function(x){nrow(x)}))),
                              avg_length = mean(unlist(lapply(train_syn_tc, function(x){nrow(x)}))),
                              median_length = median(unlist(lapply(train_syn_tc, function(x){nrow(x)})))))

sum_train <- rbind(sum_train,
                   data.frame(train_type = "syn_vac",
                              n_ts = length(train_syn_vac),
                              total_obs = sum(unlist(lapply(train_syn_vac, function(x){nrow(x)}))),
                              avg_length = mean(unlist(lapply(train_syn_vac, function(x){nrow(x)}))),
                              median_length = median(unlist(lapply(train_syn_vac, function(x){nrow(x)})))))
sum_train

# ---------------------------
# Define hyperparameters
# ---------------------------
context_length <- 20
horizon <- 4

# This dense set of quantiles is used for the bottom-up forecasting of VACs later
# The quantiles for WIS computation is only c(0.025, 0.10, 0.25, 0.5, 0.75, 0.90, 0.975)
quantile_vec <- c(0.0005, 0.005, 0.01, 0.025,
                  seq(0.05, 0.95, by = 0.05),
                  0.975, 0.99, 0.995, 0.9995)

# Cosine LR schedule hyperparameters (start -> end over weight_updates)
lr_start <- 5e-4
lr_end   <- 5e-5

# Mini-batching and weight updates
batch_size <- 256
weight_updates <- 25000 #this is smaller due to the number of models needing to be fit for this exercise (the number in the paper is closer to 100k)

n_head <- 4
num_layers <- 2
d_model <- 256
dropout <- 0

# EMA hyperparameter
ema_alpha <- 0.98

print_every <- 250
val_n <- 25000
val_seed <- 1122233

# ---------------------------
# input-noise multiplier hyperparameter
mult_alpha <- 0.15  # y_in' = pmax(0, y_in * runif(., 1-mult_alpha, 1+mult_alpha))

# ---------------------------
# Utilities: normalization + loss
# ---------------------------

.normalize_by_input_max <- function(y_in, y_out) {
  m <- max(y_in, na.rm = TRUE)
  if (!is.finite(m)) m <- 0
  
  if (m > 0) {
    z_in  <- y_in  / m
    z_out <- y_out / m
    scale <- m
  } else {
    z_in  <- y_in
    z_out <- y_out
    scale <- 0
  }
  
  list(scale = scale, z_in = z_in, z_out = z_out)
}

.pinball_loss <- function(pred, target, quantiles_tensor) {
  q <- quantiles_tensor$view(c(1, 1, -1))     # [1,1,Q]
  e <- target - pred                          # [B,H,Q]
  torch_mean(torch_maximum(q * e, (q - 1) * e))
}

.clone_state_dict <- function(state) {
  lapply(state, function(t) t$clone())
}

.ema_update_state <- function(ema_state, cur_state, alpha) {
  one_minus <- 1 - alpha
  for (nm in names(ema_state)) {
    ema_state[[nm]] <- ema_state[[nm]] * alpha + cur_state[[nm]] * one_minus
  }
  ema_state
}

.cosine_lr <- function(step, T, lr_start, lr_end) {
  t <- (step - 1) / max(1, (T - 1))  # 0..1
  lr_end + 0.5 * (lr_start - lr_end) * (1 + cos(pi * t))
}

.set_optimizer_lr <- function(optimizer, lr_value) {
  for (pg in optimizer$param_groups) {
    pg$lr <- lr_value
  }
}

# ---------------------------
# NEW: helper to perturb y_in
# ---------------------------
.perturb_y_in <- function(y_in, alpha = 0.1) {
  if (!is.finite(alpha) || alpha < 0 || alpha > 1) stop("alpha must be in [0,1].")
  mult <- stats::runif(length(y_in), min = 1 - alpha, max = 1 + alpha)
  pmax(0, y_in * mult)
}

.get_y <- function(df) {
  ord <- order(df[["time"]])
  y <- as.numeric(df[["cases"]][ord])
  y[is.finite(y)]
}

# ---------------------------
# Data preparation: indexing
# ---------------------------

.prepare_series_index <- function(train, C, H) {
  n_series <- length(train)
  if (n_series == 0) stop("train is empty.")
  
  n_cards <- integer(n_series)
  for (i in seq_len(n_series)) {
    df <- train[[i]]
    if (!("time" %in% names(df)) || !("cases" %in% names(df))) {
      stop(sprintf("Series %d missing required columns 'time' and/or 'cases'.", i))
    }
    
    y <- .get_y(df)
    
    n <- length(y)
    n_cards[i] <- max(0L, n - C - H + 1L)
  }
  
  total_cards <- sum(n_cards)
  if (total_cards <= 0) stop("No flash cards can be formed. Need length(cases) >= C + H.")
  
  probs <- n_cards / total_cards
  list(n_cards = n_cards, probs = probs, total_cards = total_cards)
}


## new for exercise
.make_flashcard_batch <- function(train, C, H, batch_size, series_probs,
                                  mult_alpha = 0.1,
                                  device = "cpu",
                                  flashcard_pool = NULL) {
  
  B2 <- 2L * as.integer(batch_size)
  
  x_mat <- matrix(0, nrow = B2, ncol = C)
  y_out_mat <- matrix(0, nrow = B2, ncol = H)
  scale_vec <- numeric(B2)
  z_last_vec <- numeric(B2)
  
  if (!is.null(flashcard_pool)) {
    pool_n <- as.integer(flashcard_pool$n)
  } else {
    n_series <- length(train)
  }
  
  for (b in seq_len(batch_size)) {
    
    # keep trying until we get a valid card
    repeat {
      if (!is.null(flashcard_pool)) {
        pick <- sample.int(pool_n, size = 1L)
        s_idx <- flashcard_pool$series_idx[pick]
        start <- flashcard_pool$start_idx[pick]
      } else {
        s_idx <- sample.int(length(train), size = 1L, prob = series_probs$probs)
        start <- NA_integer_
      }
      
      df <- train[[s_idx]]
      y <- .get_y(df)
      
      n <- length(y)
      max_start <- n - C - H + 1L
      if (max_start < 1L) next
      
      if (is.na(start)) start <- sample.int(max_start, size = 1L)
      if (start < 1L || start > max_start) next
      
      end_ctx <- start + C - 1L
      y_in  <- y[start:end_ctx]
      y_out <- y[(end_ctx + 1L):(end_ctx + H)]
      
      if (length(y_in) != C || length(y_out) != H) next
      
      break
    }
    
    # ---- original card ----
    st <- .normalize_by_input_max(y_in, y_out)
    z_in <- st$z_in
    m <- st$scale
    
    x_mat[b, ] <- z_in
    z_last_vec[b] <- z_in[C]
    y_out_mat[b, ] <- y_out
    scale_vec[b] <- m
    
    # ---- noisy twin card ----
    y_in_p <- .perturb_y_in(y_in, alpha = mult_alpha)
    st2 <- .normalize_by_input_max(y_in_p, y_out)
    
    b2 <- b + batch_size
    x_mat[b2, ] <- st2$z_in
    z_last_vec[b2] <- st2$z_in[C]
    y_out_mat[b2, ] <- y_out
    scale_vec[b2] <- st2$scale
  }
  
  list(
    x = torch_tensor(x_mat, dtype = torch_float(), device = device),
    y_out = torch_tensor(y_out_mat, dtype = torch_float(), device = device),
    scale = torch_tensor(scale_vec, dtype = torch_float(), device = device),
    z_last = torch_tensor(z_last_vec, dtype = torch_float(), device = device)
  )
}

# ---------------------------
# NEW: Flashcard pool utilities
# ---------------------------

# Given series_index$n_cards and a vector of global card ids in 1..total_cards,
# map each id to (series_idx, start_idx_within_series)
.card_ids_to_pairs <- function(n_cards, card_ids) {
  n_cards <- as.integer(n_cards)
  cs <- cumsum(n_cards)
  total <- cs[length(cs)]
  
  card_ids <- as.integer(card_ids)
  if (any(card_ids < 1L) || any(card_ids > total)) {
    stop("card_ids out of range.")
  }
  
  # Safe mapping (handles card_ids == total correctly)
  series_idx <- findInterval(card_ids, cs, rightmost.closed = TRUE) + 1L
  
  # Clamp (extra safety)
  series_idx <- pmin(series_idx, length(n_cards))
  
  prev_cs <- c(0L, cs)[series_idx]
  start_idx <- card_ids - prev_cs
  
  # start_idx should be 1..n_cards[series_idx]; clamp as safety
  start_idx <- pmax(1L, pmin(start_idx, n_cards[series_idx]))
  
  list(
    series_idx = as.integer(series_idx),
    start_idx  = as.integer(start_idx),
    n = length(card_ids),
    total_cards = as.integer(total)
  )
}

# Build a nested pool: one random permutation of card ids up to K_max
.make_nested_card_id_perm <- function(total_cards, K_max, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  total_cards <- as.integer(total_cards)
  K_max <- as.integer(K_max)
  if (K_max < 1L) stop("K_max must be >= 1")
  if (K_max > total_cards) K_max <- total_cards
  sample.int(total_cards, size = K_max, replace = FALSE)
}

# ---------------------------
# Model components
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
    
    # x has length self$C (exactly context_length observations)
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


# ---------------------------
# Training: main entry point
# ---------------------------

train_model <- function(
    train,
    context_length = 20,
    horizon = 4,
    quantile_vec = seq(0.01, 0.99, 0.01),
    batch_size = 256,
    lr_start = 5e-3,
    lr_end = 1e-4,
    weight_updates = 100000,
    ema_alpha = 0.98,
    mult_alpha = 0.1,
    d_model = 256,
    n_head = 4,
    num_layers = 2,
    dim_feedforward = 512,
    dropout = 0,
    val_n = NULL,
    val_frac = 0.01,
    val_seed = 123,
    print_every = 10,
    device = NULL,
    # NEW:
    flashcard_pool = NULL,
    total_possible_flashcards_override = NULL
) {
  if (!is.list(train)) stop("train must be a list of data.frames/data.tables/tibbles.")
  C_obs <- as.integer(context_length)
  H <- as.integer(horizon)
  
  if (C_obs < 1 || H < 1) stop("context_length and horizon must be >= 1.")
  if (batch_size < 1) stop("batch_size must be >= 1.")
  if (weight_updates < 1) stop("weight_updates must be >= 1.")
  if (!is.finite(ema_alpha) || ema_alpha < 0 || ema_alpha >= 1) stop("ema_alpha must be in [0, 1).")
  if (!is.finite(lr_start) || !is.finite(lr_end) || lr_start <= 0 || lr_end <= 0) stop("lr_start and lr_end must be > 0.")
  if (!is.finite(mult_alpha) || mult_alpha < 0 || mult_alpha > 1) {
    stop("mult_alpha must be in [0,1].")
  }
  
  if (any(!is.finite(quantile_vec)) || any(quantile_vec <= 0) || any(quantile_vec >= 1)) {
    stop("quantile_vec must be in (0,1).")
  }
  quantile_vec <- sort(unique(as.numeric(quantile_vec)))
  Q <- length(quantile_vec)
  if (Q < 1) stop("quantile_vec must contain at least one quantile.")
  
  if (is.null(device)) device <- if (cuda_is_available()) "cuda" else "cpu"
  device <- torch_device(device)
  
  series_index <- .prepare_series_index(train, C_obs, H)
  
  # If using a pool, redefine what "total cards" means for this training run
  effective_total_cards <- if (!is.null(flashcard_pool)) as.integer(flashcard_pool$n) else as.integer(series_index$total_cards)
  true_total_cards <- as.integer(series_index$total_cards)
  
  if (is.null(val_n)) {
    val_n <- max(1L, as.integer(round(effective_total_cards * val_frac)))
  } else {
    val_n <- as.integer(val_n)
  }
  if (val_n < 1L) val_n <- 1L
  
  set.seed(val_seed)
  if (!is.null(flashcard_pool)) {
    # Validation samples from the SAME pool of allowed cards
    val_pick <- sample.int(flashcard_pool$n, size = val_n, replace = TRUE)
    val_series_idx <- flashcard_pool$series_idx[val_pick]
    val_start_idx  <- flashcard_pool$start_idx[val_pick]
  } else {
    # Original behavior
    val_series_idx <- sample.int(length(train), size = val_n, replace = TRUE, prob = series_index$probs)
    val_start_idx  <- integer(val_n)
    
    for (i in seq_len(val_n)) {
      df <- train[[val_series_idx[i]]]
      y <- .get_y(df)
      n <- length(y)
      max_start <- n - C_obs - H + 1L
      if (max_start < 1L) {
        while (max_start < 1L) {
          j <- sample.int(length(train), size = 1L, prob = series_index$probs)
          df <- train[[j]]
          y <- .get_y(df)
          n <- length(y)
          max_start <- n - C_obs - H + 1L
          val_series_idx[i] <- j
        }
      }
      val_start_idx[i] <- sample.int(max_start, size = 1L)
    }
  }
  
  eval_validation <- function(model, chunk_size = 2048L) {
    model$eval()
    on.exit(model$train(), add = TRUE)
    
    quantiles_t <- torch_tensor(quantile_vec, dtype = torch_float(), device = device)
    total_loss <- 0
    total_n <- 0L
    
    draw_valid_val_card <- function() {
      repeat {
        if (!is.null(flashcard_pool)) {
          pick <- sample.int(flashcard_pool$n, size = 1L)
          s_idx <- flashcard_pool$series_idx[pick]
          s0    <- flashcard_pool$start_idx[pick]
        } else {
          s_idx <- sample.int(length(train), size = 1L, prob = series_index$probs)
          s0 <- NA_integer_
        }
        
        # <-- add this guard right here
        if (s_idx < 1L || s_idx > length(train)) next
        
        df <- train[[s_idx]]
        y <- .get_y(df)
        
        n <- length(y)
        max_start <- n - C_obs - H + 1L
        if (max_start < 1L) next
        
        if (is.na(s0)) s0 <- sample.int(max_start, size = 1L)
        if (s0 < 1L || s0 > max_start) next
        
        end_ctx <- s0 + C_obs - 1L
        y_in  <- y[s0:end_ctx]
        y_out <- y[(end_ctx + 1L):(end_ctx + H)]
        if (length(y_in) != C_obs || length(y_out) != H) next
        
        return(list(y_in = y_in, y_out = y_out))
      }
    }
    
    with_no_grad({
      idx <- 1L
      while (idx <= val_n) {
        end_idx <- min(val_n, idx + chunk_size - 1L)
        m <- end_idx - idx + 1L
        
        x_mat <- matrix(0, nrow = m, ncol = C_obs)
        y_out_mat <- matrix(0, nrow = m, ncol = H)
        scale_vec <- numeric(m)
        z_last_vec <- numeric(m)
        
        for (k in seq_len(m)) {
          card <- draw_valid_val_card()
          
          st <- .normalize_by_input_max(card$y_in, card$y_out)
          z_in <- st$z_in
          mscale <- st$scale
          
          x_mat[k, ] <- z_in
          z_last_vec[k] <- z_in[C_obs]
          y_out_mat[k, ] <- card$y_out
          scale_vec[k] <- mscale
        }
        
        x <- torch_tensor(x_mat, dtype = torch_float(), device = device)              # [m, C]
        y_true <- torch_tensor(y_out_mat, dtype = torch_float(), device = device)     # [m, H]
        y_true <- y_true$unsqueeze(3)$expand(c(m, H, Q))                               # [m, H, Q]
        
        scale_t <- torch_tensor(scale_vec, dtype = torch_float(), device = device)     # [m]
        zlast_t  <- torch_tensor(z_last_vec, dtype = torch_float(), device = device)  # [m]
        
        r_hat <- model(x)
        z_hat <- r_hat + zlast_t$view(c(m, 1, 1))
        y_hat <- z_hat * scale_t$view(c(m, 1, 1))
        
        loss <- .pinball_loss(y_hat, y_true, quantiles_t)
        
        total_loss <- total_loss + as.numeric(loss$item()) * m
        total_n <- total_n + m
        idx <- end_idx + 1L
      }
    })
    
    total_loss / max(1L, total_n)
  }
  
  model <- TransformerForecaster(
    context_length = C_obs,
    horizon = H,
    n_quantiles = Q,
    d_model = d_model,
    n_head = n_head,
    num_layers = num_layers,
    dim_feedforward = dim_feedforward,
    dropout = dropout
  )$to(device = device)
  
  optimizer <- optim_adam(model$parameters, lr = lr_start)
  quantiles_t <- torch_tensor(quantile_vec, dtype = torch_float(), device = device)
  
  model$train()
  
  # RAW best
  best_raw_val <- Inf
  best_raw_step <- 0L
  best_raw_state <- NULL
  
  # EMA (delayed start) + best EMA
  ema_active <- FALSE
  ema_started_step <- NA_integer_
  ema_state <- NULL
  
  best_ema_val <- Inf
  best_ema_step <- NA_integer_
  best_ema_state <- NULL
  
  # Init validation (EMA is OFF)
  raw_val <- eval_validation(model)
  cur_state <- model$state_dict()
  
  best_raw_val <- raw_val
  best_raw_step <- 0L
  best_raw_state <- .clone_state_dict(cur_state)
  
  cat(sprintf("Init: val_raw=%.6f * | val_ema=NA (ema_off)\n", raw_val))
  
  for (step in seq_len(weight_updates)) {
    lr_t <- .cosine_lr(step, weight_updates, lr_start, lr_end)
    .set_optimizer_lr(optimizer, lr_t)
    
    batch <- .make_flashcard_batch(
      train = train,
      C = C_obs,
      H = H,
      batch_size = batch_size,
      series_probs = series_index,
      mult_alpha = mult_alpha,
      device = device,
      flashcard_pool = flashcard_pool
    )
    
    # NEW: effective batch size is doubled
    B2 <- 2L * as.integer(batch_size)
    
    x <- batch$x
    y_true <- batch$y_out$unsqueeze(3)$expand(c(B2, H, Q))
    scale_t <- batch$scale
    zlast_t <- batch$z_last
    
    optimizer$zero_grad()
    
    r_hat <- model(x)
    z_hat <- r_hat + zlast_t$view(c(B2, 1, 1))
    y_hat <- z_hat * scale_t$view(c(B2, 1, 1))
    
    loss <- .pinball_loss(y_hat, y_true, quantiles_t)
    loss$backward()
    optimizer$step()
    
    cur_state <- model$state_dict()
    
    if (ema_active) {
      ema_state <- .ema_update_state(ema_state, cur_state, ema_alpha)
    }
    
    if (step %% print_every == 0L || step == 1L || step == weight_updates) {
      raw_val <- eval_validation(model)
      
      is_new_best_raw <- is.finite(raw_val) && (raw_val < best_raw_val)
      if (is_new_best_raw) {
        best_raw_val <- raw_val
        best_raw_step <- step
        best_raw_state <- .clone_state_dict(cur_state)
      }
      
      if (!ema_active && !is_new_best_raw) {
        ema_active <- TRUE
        ema_started_step <- step
        ema_state <- .clone_state_dict(cur_state)
        
        model$load_state_dict(ema_state)
        ema_val_init <- eval_validation(model)
        model$load_state_dict(cur_state)
        
        best_ema_val <- ema_val_init
        best_ema_step <- step
        best_ema_state <- .clone_state_dict(ema_state)
      }
      
      if (ema_active) {
        model$load_state_dict(ema_state)
        ema_val <- eval_validation(model)
        model$load_state_dict(cur_state)
        
        is_new_best_ema <- is.finite(ema_val) && (ema_val < best_ema_val)
        if (is_new_best_ema) {
          best_ema_val <- ema_val
          best_ema_step <- step
          best_ema_state <- .clone_state_dict(ema_state)
        }
        
        cat(sprintf(
          "step=%d/%d | lr=%.3e | train=%.6f | val_raw=%.6f%s | val_ema=%.6f%s\n",
          step, weight_updates, lr_t, as.numeric(loss$item()),
          raw_val, if (is_new_best_raw) " *" else "",
          ema_val, if (is_new_best_ema) " *" else ""
        ))
      } else {
        cat(sprintf(
          "step=%d/%d | lr=%.3e | train=%.6f | val_raw=%.6f%s | val_ema=NA (ema_off)\n",
          step, weight_updates, lr_t, as.numeric(loss$item()),
          raw_val, if (is_new_best_raw) " *" else ""
        ))
      }
    }
  }
  
  # Ensure best EMA exists (if EMA never started)
  if (!ema_active) {
    ema_active <- TRUE
    ema_started_step <- weight_updates
    ema_state <- .clone_state_dict(model$state_dict())
    
    model$load_state_dict(ema_state)
    ema_val_final <- eval_validation(model)
    model$load_state_dict(model$state_dict())
    
    best_ema_val <- ema_val_final
    best_ema_step <- weight_updates
    best_ema_state <- .clone_state_dict(ema_state)
  }
  
  list(
    best_raw_state = best_raw_state,
    best_ema_state = best_ema_state,
    config = list(
      normalization = "divide_by_max(y_in) with fallback max=0 => identity",
      extra_feature = NULL,
      context_length_obs = C_obs,
      context_length_in = C_obs,
      horizon = H,
      quantile_vec = quantile_vec,
      batch_size = batch_size,
      effective_batch_size = 2L * as.integer(batch_size),
      lr_start = lr_start,
      lr_end = lr_end,
      weight_updates = weight_updates,
      ema_alpha = ema_alpha,
      ema_started_step = ema_started_step,
      d_model = d_model,
      n_head = n_head,
      num_layers = num_layers,
      dim_feedforward = dim_feedforward,
      dropout = dropout,
      mult_alpha = mult_alpha,
      noise_sd = NULL,
      noise_transform = "y_in' = pmax(0, y_in * runif(., 1-alpha, 1+alpha))",
      val_n = val_n,
      val_frac = val_frac,
      val_seed = val_seed,
      device = as.character(device),
      cards_seen = as.integer(weight_updates) * as.integer(batch_size) * 2L,
      total_possible_flashcards_true = true_total_cards,
      total_possible_flashcards_used = if (is.null(total_possible_flashcards_override)) effective_total_cards else as.integer(total_possible_flashcards_override),
      series_n_flashcards = series_index$n_cards,
      best_raw_val_pinball_loss = best_raw_val,
      best_raw_step = best_raw_step,
      best_ema_val_pinball_loss = best_ema_val,
      best_ema_step = best_ema_step
    )
  )
}


# ---------------------------
# Fit and save the models for nested flashcard subsets
# ---------------------------

card_sizes <- c(2500L, 25000L, 250000L, 1000000L, 2500000L, 5000000L, 10000000L, 25000000L)

# Lookup by ms = 1:4
train_list_by_ms <- list(
  train_real,
  train_syn_tc,
  train_syn_vac,
  c(train_real, train_syn_tc, train_syn_vac)
)

name_by_ms <- c("real", "syn_tc", "syn_vac", "all")
seed_by_ms <- c(111, 222, 333, 444)

for (ms in 4) {
  
  temp_train <- train_list_by_ms[[ms]]
  
  # Optional: keep your original series-order randomization
  set.seed(seed_by_ms[ms])
  perm <- sample.int(length(temp_train))
  temp_train <- temp_train[perm]
  
  base_model_name <- paste0("mod_", name_by_ms[ms])
  
  # Full series index for this training set
  series_index_full <- .prepare_series_index(temp_train, context_length, horizon)
  total_cards_full <- as.integer(series_index_full$total_cards)
  
  # Which K to run (<= total_cards_full), plus "ALL" (total_cards_full)
  K_list <- card_sizes[card_sizes <= total_cards_full]
  if (!any(K_list == total_cards_full)) K_list <- c(K_list, total_cards_full)
  K_list <- unique(as.integer(K_list))
  K_list <- K_list[order(K_list)]
  
  # Build one nested permutation up to max(K_list)
  K_max <- max(K_list)
  ids_perm <- .make_nested_card_id_perm(total_cards_full, K_max, seed = seed_by_ms[ms] + 999)
  
  cat(sprintf("\n== %s: total_cards=%d | runs: %s ==\n",
              base_model_name, total_cards_full, paste(K_list, collapse = ", ")))
  
  for (K in K_list) {
    
    # Nested subset of card ids
    ids_K <- ids_perm[seq_len(K)]
    
    # Map ids -> (series_idx, start_idx)
    pool <- .card_ids_to_pairs(series_index_full$n_cards, ids_K)
    
    # Train
    start_time <- Sys.time()
    
    fit <- train_model(
      train = temp_train,
      context_length = context_length,
      horizon = horizon,
      quantile_vec = quantile_vec,
      batch_size = batch_size,
      lr_start = lr_start,
      lr_end = lr_end,
      weight_updates = weight_updates,
      ema_alpha = ema_alpha,
      mult_alpha = mult_alpha,
      n_head = n_head,
      num_layers = num_layers,
      d_model = d_model,
      print_every = print_every,
      val_n = val_n,
      val_seed = val_seed,
      
      flashcard_pool = pool,
      total_possible_flashcards_override = K
    )
    
    end_time <- Sys.time()
    run_time <- difftime(end_time, start_time, units = "mins")
    
    # Config
    cfg <- fit$config
    cfg$run_time <- run_time
    cfg$model_family <- name_by_ms[ms]
    cfg$subset_cards_used <- as.integer(K)
    
    # Save best EMA only
    device_save <- if (torch::cuda_is_available()) "cuda" else "cpu"
    model_obj <- TransformerForecaster(
      context_length = context_length,
      horizon = horizon,
      n_quantiles = length(sort(unique(as.numeric(quantile_vec)))),
      d_model = d_model,
      n_head = n_head,
      num_layers = num_layers,
      dim_feedforward = 512,
      dropout = dropout
    )$to(device = torch::torch_device(device_save))
    
    model_obj$load_state_dict(fit$best_ema_state)
    
    model_file <- file.path(modelpath, paste0(base_model_name, "_", K, ".pt"))
    cfg_file   <- file.path(modelpath, paste0("cfg_", name_by_ms[ms], "_", K, ".rds"))
    
    torch::torch_save(model_obj$state_dict(), model_file)
    saveRDS(cfg, cfg_file)
    
    cat(sprintf("Saved: %s and %s\n", basename(model_file), basename(cfg_file)))
  }
}
