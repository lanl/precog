# depends: estimate_hi_from_outputs, predict_feas_glm, osfd_ei_acquisition
#' @export
propose_batch <- function(
  X_succ, Y_succ,
  remaining_idx, CAND,
  feas_model,
  B = 8,
  cand_batch = 20000,
  beta = 2.0,
  p_floor = 0.05,
  tau_hard = NULL,
  repel_radius = 0.05   # in [0,1]^p units; tune
) {
  # --- compute OSFD quantities from successful outputs ---
  sc <- minmax_scale(Y_succ)
  hi <- estimate_hi_from_outputs(sc$Y_sc)

  # --- sample a big evaluation pool from remaining candidates ---
  k <- min(cand_batch, length(remaining_idx))
  pool_idx <- sample(remaining_idx, k, replace = FALSE)
  X_pool <- CAND[pool_idx, , drop = FALSE]

  # --- feasibility probs ---
  p_feas <- predict_feas_glm(feas_model, X_pool, fallback_p = 1.0)
  p_feas <- pmin(1, pmax(p_floor, p_feas))

  if (!is.null(tau_hard)) {
    ok <- p_feas >= tau_hard
    if (any(ok)) {
      X_pool <- X_pool[ok, , drop = FALSE]
      pool_idx <- pool_idx[ok]
      p_feas <- p_feas[ok]
    }
  }

  # --- base EI score ---
  EI <- osfd_ei_acquisition(X_obs = X_succ, hi = hi, X_cand = X_pool)
  base_score <- EI * (p_feas ^ beta)

  # --- greedy within-batch selection with repulsion ---
  chosen <- integer(0)
  chosen_X <- matrix(numeric(0), ncol = ncol(CAND))

  score <- base_score

  for (b in seq_len(B)) {
    if (length(score) == 0 || all(score <= 0)) break

    j <- which.max(score)
    chosen <- c(chosen, pool_idx[j])
    chosen_X <- rbind(chosen_X, X_pool[j, , drop = FALSE])

    # repulsion: downweight points near newly chosen point in input space
    dx <- sqrt(rowSums((X_pool - matrix(X_pool[j,], nrow(X_pool), ncol(X_pool), byrow=TRUE))^2))
    score <- score * pmin(1, dx / repel_radius)  # ~0 near, ~1 far

    # avoid reselecting the same point
    score[j] <- -Inf
  }

  list(idx = chosen, X = chosen_X)
}
