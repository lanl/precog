# depends: osfd_helpers, RcppExports
#' @name 
#' spanfill_lite
#' @title
#' Generate points to approximate the space spanned by the existing points.  Murph experiments to speed it up.
#
#' @description
#' \code{spanfill} generates points to approximate a space based on existing points.
#' These approximate points can be used to find local fill distance in the space or be used as candidate points in active learning.
#' 
#' @details 
#' \code{spanfill} generates points to approximate the space spanned by the existing points. Details can be found in  Wang et al. (2024).
#' 
#' @param X a matrix specifying the existing points
#' @param bound a binary variable indicating whether to bound the generated points to 0 to 1 in each dimension. If bound=TRUE, all the generated points will be projected to the unit hypercube. The default value is FALSE.
#' 
#' @return a matrix of the generated points to approximate the space.
#' 
#' @references 
#' Wang, Shangkun, Adam P. Generale, Surya R. Kalidindi, and V. Roshan Joseph. (2024). "Sequential designs for filling output spaces", Technometrics, 66, 65–76.
#'
#' @export
#' 
#' @examples
#' 
#' X = matrix(runif(20), ncol=2)
#' spanfill_points = spanfill(X)
#' plot(spanfill_points, type='p')
#' 
spanfill_lite = function (X, bound=FALSE){
  print("Doing more progressing printing for the spanfill function that I've been modifiying.")
  q = ncol(X)
  X.u = unique(X)
  no.u = dim(X.u)[1]

  # Check if we have enough unique points for meaningful KNN
  # Need at least q+2 points: one query point + at least q+1 neighbors for proper subspace construction
  min_required <- q + 2
  if(no.u < min_required) {
    warning(sprintf("spanfill_lite: Only %d unique points, need at least %d for q=%d. Returning input.",
                    no.u, min_required, q))
    return(X.u)
  }

  no.nb = min(2*q, no.u-1)

  # Additional safety check: no.nb must be positive
  if(no.nb < 1) {
    warning(sprintf("spanfill_lite: no.nb=%d is too small (no.u=%d, q=%d). Returning input.",
                    no.nb, no.u, q))
    return(X.u)
  }

  osfd_knn_my                <- get_osfd_fn("knn_my")
  osfd_runif_in_sphere_cpp   <- get_osfd_fn("runif_in_sphere_cpp")

  # ---- timing helpers (inline, no function) ----
  t0 <- proc.time()[["elapsed"]]

  # print("Performing the KNN on the data...")
  t_knn0 <- proc.time()[["elapsed"]]
  knn.result <- osfd_knn_my(X.u, no.nb)
  t_knn1 <- proc.time()[["elapsed"]]
  # cat(sprintf("KNN step took %.3f sec (%.2f min; %.2f hr)\n",
  #             t_knn1 - t_knn0, (t_knn1 - t_knn0)/60, (t_knn1 - t_knn0)/3600))

  # points in balls
  n_ball <- 2*q + 2*(q + 1) + 1  # number of points in each ball

  # print("Sampling uniform points from the ball...")
  t_ball0 <- proc.time()[["elapsed"]]
  twinsample <- osfd_runif_in_sphere_cpp(100*n_ball, q)
  t_ball1 <- proc.time()[["elapsed"]]
  # cat(sprintf("Ball sampling step took %.3f sec (%.2f min; %.2f hr)\n",
  #             t_ball1 - t_ball0, (t_ball1 - t_ball0)/60, (t_ball1 - t_ball0)/3600))

  # print("Twinning (sub-sampling) the data...")
  t_twin0 <- proc.time()[["elapsed"]]
  twinsample <- twinsample[twinning::twin(twinsample, r = 100), , drop = FALSE]
  t_twin1 <- proc.time()[["elapsed"]]
  # cat(sprintf("Twinning step took %.3f sec (%.2f min; %.2f hr)\n",
  #             t_twin1 - t_twin0, (t_twin1 - t_twin0)/60, (t_twin1 - t_twin0)/3600))

  print("Calculating the approximating points (this now is in parallel on the C++ level).")
  cat(sprintf("DEBUG: X.u has %d rows, %d cols; no.nb=%d; q=%d\n",
              nrow(X.u), ncol(X.u), no.nb, q))
  t_apx0 <- proc.time()[["elapsed"]]

  approx_points <- approx_gen(X.u, knn.result, q, FALSE, twinsample)
  cat(sprintf("DEBUG: approx_gen returned %d points\n", nrow(approx_points)))
  if (bound) {
    approx_points <- pmax(pmin(approx_points, 1), 0)
  }
  t_apx1 <- proc.time()[["elapsed"]]
  cat(sprintf("Determining candidate output space points took %.3f sec (%.2f min; %.2f hr)\n",
              t_apx1 - t_apx0, (t_apx1 - t_apx0)/60, (t_apx1 - t_apx0)/3600))

  t1 <- proc.time()[["elapsed"]]
  # cat(sprintf("TOTAL spanfill sections took %.3f sec (%.2f min; %.2f hr)\n",
  #             t1 - t0, (t1 - t0)/60, (t1 - t0)/3600))

  return (approx_points)
}