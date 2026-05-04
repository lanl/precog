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