# depends: 
#' @export
nn1_excluding_self <- function(X) {
  # for each row i, nearest neighbor among rows != i
  if (nrow(X) < 2) stop("Need at least 2 points for nearest-neighbor LOOCV.")
  if (requireNamespace("RANN", quietly = TRUE)) {
    ans <- RANN::nn2(data = X, query = X, k = 2)  # first is self
    list(idx = ans$nn.idx[,2], dist = ans$nn.dists[,2])
  } else {
    idx <- integer(nrow(X))
    dist <- numeric(nrow(X))
    D2 <- as.matrix(dist(X))^2
    diag(D2) <- Inf
    for (i in seq_len(nrow(X))) {
      j <- which.min(D2[i,])
      idx[i] <- j
      dist[i] <- sqrt(D2[i,j])
    }
    list(idx = idx, dist = dist)
  }
}