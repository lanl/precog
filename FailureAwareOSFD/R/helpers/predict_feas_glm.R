#' Predict Feasibility Probabilities for New Input Points
#'
#' This function uses a fitted feasibility model (from fit_feas_model_glm) to
#' predict the probability that each new input point will result in a successful
#' function evaluation. These probabilities are used to guide the sequential
#' design toward feasible regions.
#'
#' @param model A fitted GLM object from fit_feas_model_glm, or NULL
#' @param X_new Matrix (n x p): new input points at which to predict feasibility
#' @param fallback_p Numeric in [0,1]: probability to return when model is NULL
#'   or when predictions are non-finite. Default is 1.0 (assume all points feasible).
#'
#' @return Numeric vector of length n: predicted probability of success for each
#'   row of X_new. All probabilities are in [0, 1], with non-finite values
#'   replaced by fallback_p.
#'
#' @details
#' If model is NULL (because all past attempts had the same outcome), this
#' function returns fallback_p for all points. The fallback should typically
#' be set to the empirical success rate: mean(success_flag).
#'
#' Non-finite predictions can occur due to:
#' - Numerical overflow/underflow in the GLM predictions
#' - Extrapolation far outside the range of the training data
#' These are replaced with fallback_p as a conservative estimate.
#'
#' Note: The returned probabilities are NOT clipped to [p_floor, 1] - that
#' clipping happens in the main OSFD function.
#'
#' @export
predict_feas_glm <- function(model, X_new, fallback_p = 1.0) {

  # If no model was fitted (all past attempts had same outcome),
  # return the fallback probability for all points
  if (is.null(model)) {
    return(rep(fallback_p, nrow(X_new)))
  }

  # Prepare data frame for prediction (must match structure of training data)
  df <- as.data.frame(X_new)

  # Get predicted probabilities using logistic regression
  # type = "response" returns probabilities (not log-odds)
  # We suppress warnings because extrapolation warnings are common and expected
  p <- suppressWarnings(
    predict(model, newdata = df, type = "response")
  )

  # Ensure predictions are numeric (handles edge cases)
  p <- as.numeric(p)

  # Replace any non-finite predictions (NA, NaN, Inf) with fallback
  p[!is.finite(p)] <- fallback_p

  return(p)
}
