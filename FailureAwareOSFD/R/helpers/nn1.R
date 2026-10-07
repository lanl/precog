# depends: osfd_helpers
#' @export
nn1 <- function(data, query) {
  # We'll assume it returns a list with fields like idx / dist (or similar).
  osfd_knnx_my   <- get_osfd_fn("knnx_my")
  ans <- osfd_knnx_my(data, query, k = 1)

  # Make it robust to small naming differences
  idx <- NULL; dist <- NULL

  idx <- ans[["nn_index"]]
  dist <- ans[["nn_dist"]]

  if (is.null(idx) || is.null(dist)) {
    stop("Couldn't parse output of OSFD::knnx_my. Run str(OSFD::knnx_my(data, query, k=1)).")
  }

  # Ensure vectors
  idx <- as.integer(idx) + 1
  dist <- as.numeric(dist)

  list(idx = idx, dist = dist)
}
