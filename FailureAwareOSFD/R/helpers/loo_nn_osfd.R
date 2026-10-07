# depends: osfd_helpers
#' @export
loo_nn_osfd <- function(X) {
  # knn_my likely computes neighbor lists within the same set (includes self as nearest),
  # so we request k=2 and take the second neighbor.
  osfd_knn_my    <- get_osfd_fn("knn_my")
  ans <- osfd_knn_my(X, k = 2)

  # Parse
  idx_mat <- NULL; dist_mat <- NULL

  idx_mat <- ans[["nn_index"]]
  dist_mat <- ans[["nn_dist"]]

  list(
    idx  = as.integer(idx_mat[,2]) + 1,
    dist = as.numeric(dist_mat[,2])
  )
}
