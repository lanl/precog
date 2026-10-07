# depends: 
#' @export
scale_rows_X_then_apply_to_y <- function(X, y, eps = 1e-8, zero_sd = c("eps", "skip")) {
  zero_sd <- match.arg(zero_sd)

  # row mean + row sd of X (fast)
  mu <- matrixStats::rowMeans2(X)
  s  <- matrixStats::rowSds(X)

  if (zero_sd == "eps") {
    # avoid division by ~0 by adding eps
    s_adj <- s + eps
    Xs <- sweep(X, 1, mu, FUN = "-")
    Xs <- sweep(Xs, 1, s_adj, FUN = "/")

    ys <- sweep(y, 1, mu, FUN = "-")
    ys <- sweep(ys, 1, s_adj, FUN = "/")

  } else {
    # "skip": center always; only divide rows with nonzero sd
    good <- s > 0

    Xs <- sweep(X, 1, mu, FUN = "-")
    ys <- sweep(y, 1, mu, FUN = "-")

    Xs[good, ] <- Xs[good, , drop = FALSE] / s[good]
    ys[good, ] <- ys[good, , drop = FALSE] / s[good]
    # rows with sd==0 are just mean-centered (no scaling)
  }

  list(X = Xs, y = ys, mu = mu, sd = s)
}