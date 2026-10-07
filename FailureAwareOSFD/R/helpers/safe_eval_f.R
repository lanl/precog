#' Safely Evaluate a Black-Box Function with Error Handling
#'
#' This function wraps the evaluation of f(x) with error handling to catch any
#' runtime errors, and validates that the output has the expected dimension and
#' contains only finite values. This is critical for OSFD because some input
#' points may lie in infeasible regions where the function fails or returns
#' invalid outputs.
#'
#' @param f Function to evaluate: should accept a length-p numeric vector and
#'   return a length-q numeric vector
#' @param x Numeric vector of length p: the input point at which to evaluate f
#' @param q Integer: the expected dimension of the output (used for validation)
#'
#' @return List with three components:
#'   - success: Logical, TRUE if evaluation succeeded and returned valid output
#'   - y: Numeric vector of length q, the output values (or NA if failed)
#'   - msg: Character string, error message if failed (empty string if succeeded)
#'
#' @details
#' The function is considered to have succeeded if ALL of the following hold:
#'   1. f(x) executes without throwing an error
#'   2. The output can be coerced to numeric
#'   3. The output has exactly length q
#'   4. All output values are finite (not NA, NaN, Inf, or -Inf)
#'
#' If any of these conditions fail, success = FALSE and the error message
#' is captured in msg.
#'
#' @export
safe_eval_f <- function(f, x, q) {

  # Wrap evaluation in tryCatch to handle any runtime errors
  out <- tryCatch(
    {
      # Attempt to evaluate the function at x
      y <- f(x)

      # Coerce output to numeric (handles cases where f returns integer, etc.)
      y <- as.numeric(y)

      # Validate output: must have correct length and all finite values
      ok <- (length(y) == q) && all(is.finite(y))

      if (!ok) {
        # Output failed validation
        list(
          success = FALSE,
          y = rep(NA_real_, q),  # Return vector of NAs with correct length
          msg = "Returned wrong length or non-finite values."
        )
      } else {
        # Output passed validation
        list(
          success = TRUE,
          y = y,
          msg = ""  # Empty message indicates success
        )
      }
    },
    error = function(e) {
      # Function threw an error during evaluation
      list(
        success = FALSE,
        y = rep(NA_real_, q),  # Return vector of NAs with correct length
        msg = conditionMessage(e)  # Capture the error message
      )
    }
  )

  return(out)
}
