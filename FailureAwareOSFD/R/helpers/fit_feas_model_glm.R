#' Fit a Logistic Regression Model to Predict Feasibility
#'
#' This function fits a generalized linear model (GLM) with binomial family to
#' predict the probability that an input point will result in a successful
#' evaluation of the black-box function. The model learns the shape of the
#' feasible region in input space based on past successes and failures.
#'
#' @param X_all Matrix (n x p): all input points that have been attempted so far
#' @param success Logical vector of length n: TRUE if evaluation succeeded, FALSE if failed
#'
#' @return A fitted GLM object (from stats::glm), or NULL if fitting is not possible.
#'   Returns NULL when:
#'   - All attempts succeeded (no failures to learn from)
#'   - All attempts failed (no successes to learn from)
#'
#' @details
#' The model is: logit(P(success | x)) = beta_0 + beta_1 * x_1 + ... + beta_p * x_p
#'
#' This is a simple linear model in the input space. In practice, the feasible
#' region may have a complex shape, but the linear approximation provides
#' useful guidance for avoiding infeasible regions, especially when combined
#' with the p_floor parameter in the main OSFD function.
#'
#' Warning messages from glm (e.g., about perfect separation) are suppressed
#' because:
#' - They are common when the feasible region has sharp boundaries
#' - The predictions are clipped to [p_floor, 1] in the main OSFD function
#' - The model is used as a heuristic guide, not for precise inference
#'
#' @export
fit_feas_model_glm <- function(X_all, success) {

  # Check if we have both successes and failures
  # If all attempts have the same outcome, we can't fit a discriminative model
  if (length(unique(success)) < 2) {
    return(NULL)
  }

  # Prepare data frame for glm
  df <- as.data.frame(X_all)
  df$success <- as.integer(success)  # Convert logical to 0/1

  # Fit logistic regression: success ~ all input dimensions
  # family = binomial() uses logit link by default
  # We suppress warnings about separation because:
  #   - Perfect/quasi separation is common with sharp feasibility boundaries
  #   - We handle edge cases by clipping predictions to [p_floor, 1]
  suppressWarnings(
    glm(success ~ ., data = df, family = binomial())
  )
}
