# depends: loo_nn_osfd, nn1
#' Compute Expected Improvement Acquisition Function for OSFD
#'
#' This function computes the Expected Improvement (EI) acquisition function for
#' selecting new input points in Output-Space Filling Design. The EI measures
#' the expected reduction in the maximum hole size if we were to evaluate the
#' function at each candidate input point.
#'
#' The acquisition function balances:
#'   - Exploitation: prefer candidate points predicted to have high hole index
#'   - Exploration: prefer candidate points far from observed inputs (high uncertainty)
#'
#' @param X_obs Matrix (m x p): observed input points where f has been successfully evaluated
#' @param hi Numeric vector of length m: hole index for each observed output (from estimate_hi_from_outputs)
#' @param X_cand Matrix (n x p): candidate input points at which to compute EI
#'
#' @return Numeric vector of length n: Expected Improvement for each candidate point.
#'   Higher values indicate more promising candidates for filling output space holes.
#'
#' @details
#' Algorithm (based on Gaussian Process Expected Improvement):
#'
#'   1. Estimate the prediction variance using a distance-based model:
#'      - For each observed point X_obs[i, ], find its leave-one-out (LOO)
#'        nearest neighbor in X_obs
#'      - Compute the squared prediction error: (hi[i] - hi[LOO neighbor])^2
#'      - Normalize by the distance to the LOO neighbor
#'      - Average across all observed points to get sigma2_hat (noise variance)
#'
#'   2. For each candidate point X_cand[j, ]:
#'      - Find its nearest neighbor in X_obs
#'      - Predict the hole index: mu[j] = hi[nearest neighbor]
#'      - Predict the standard deviation: s[j] = sqrt(sigma2_hat * distance)
#'
#'   3. Compute Expected Improvement:
#'      - hmax = max(hi) is the current maximum hole size
#'      - u = (mu - hmax) / s is the standardized improvement
#'      - EI = s * (u * Phi(u) + phi(u)), where Phi is CDF and phi is PDF of standard normal
#'
#' Interpretation:
#'   - EI is high when mu is high (predicted to be in a hole) AND s is high (uncertain)
#'   - EI is zero when s = 0 (candidate point equals an observed input point)
#'   - EI balances mean prediction and uncertainty, preferring high-value uncertain regions
#'
#' Technical notes:
#'   - We need at least 2 observed points to estimate variance (otherwise error)
#'   - The LOO nearest neighbor distance is clamped to >= 1e-12 to avoid division by zero
#'   - Candidate point distances are clamped to >= 0 to avoid sqrt of negative values
#'   - When s = 0, we set EI = 0 (no improvement expected from exact re-evaluation)
#'
#' @export
osfd_ei_acquisition <- function(X_obs, hi, X_cand) {

  # Number of observed points
  m <- nrow(X_obs)

  # Need at least 2 points to estimate prediction variance
  if (m < 2) {
    stop("Need at least 2 successful points before EI is meaningful.")
  }

  # Current maximum hole index (we want to fill the largest holes)
  hmax <- max(hi)

  # ==========================================================================
  # Step 1: Estimate prediction variance using LOO cross-validation
  # ==========================================================================
  # For each observed point, find its nearest neighbor among the other observed points
  # This gives us a sense of how much the hole index varies between nearby input points

  loo <- loo_nn_osfd(X_obs)

  # idx_r[i] = index of the LOO nearest neighbor for X_obs[i, ]
  idx_r <- loo$idx

  # Hole index of each point's LOO nearest neighbor
  hi_loo <- hi[idx_r]

  # Estimate variance of hole index predictions
  # For each observed point, compute: (actual hi - predicted hi)^2 / distance^2
  # This gives variance per unit distance squared
  # We clamp distance to >= 1e-12 to avoid division by zero for duplicate points
  sigma2_hat <- mean((hi - hi_loo)^2 / pmax(loo$dist, 1e-12))

  # ==========================================================================
  # Step 2: Predict hole index for each candidate point
  # ==========================================================================
  # For each candidate point, use its nearest neighbor in X_obs to predict hole index

  nnc <- nn1(X_obs, X_cand)

  # Predicted mean hole index for each candidate
  # mu[j] = hole index of the nearest observed point to X_cand[j, ]
  mu <- hi[nnc$idx]

  # Predicted standard deviation for each candidate
  # s[j] = sqrt(sigma2_hat * distance to nearest observed point)
  # The farther a candidate is from observed points, the more uncertain we are
  # We clamp distance to >= 0 to handle any numerical issues
  s <- sqrt(sigma2_hat * pmax(nnc$dist, 0))

  # ==========================================================================
  # Step 3: Compute Expected Improvement
  # ==========================================================================
  # EI formula from Jones et al. (1998) for Gaussian processes:
  #   If s > 0: EI = s * (u * Phi(u) + phi(u))
  #   If s = 0: EI = 0
  # where u = (mu - best) / s, Phi is standard normal CDF, phi is standard normal PDF

  # Initialize EI vector
  EI <- numeric(nrow(X_cand))

  # Only compute EI for points with s > 0 (points not coinciding with observed points)
  ok <- s > 0

  # Standardized improvement: how many standard deviations above current max?
  u <- (mu[ok] - hmax) / s[ok]

  # Expected Improvement formula
  # pnorm(u) = Phi(u) = P(Z <= u) for Z ~ N(0,1)
  # dnorm(u) = phi(u) = probability density at u
  EI[ok] <- s[ok] * (u * pnorm(u) + dnorm(u))

  # Points with s = 0 have EI = 0 (already evaluated or coincide with observed point)
  EI[!ok] <- 0

  return(EI)
}
