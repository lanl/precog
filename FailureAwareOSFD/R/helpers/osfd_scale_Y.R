# depends: osfd_helpers
#' Scale Output Matrix to Unit Hypercube [0,1]^q
#'
#' This function scales each column (output dimension) of the output matrix Y
#' to the range [0, 1] using min-max scaling. This is essential for OSFD
#' because:
#'   1. Different output dimensions may have vastly different scales
#'   2. Distance calculations in output space must treat all dimensions equally
#'   3. The hole-filling metric (hi) should not be dominated by large-scale outputs
#'
#' @param Y Matrix (n x q): observed output points, one row per observation
#'
#' @return Matrix (n x q): scaled output points with each column in [0, 1]
#'
#' @details
#' For each column j:
#'   Y_scaled[, j] = (Y[, j] - min(Y[, j])) / (max(Y[, j]) - min(Y[, j]))
#'
#' This is implemented via a C++ function (OSFD::scale_cpp) for efficiency.
#' The C++ function may return either:
#'   - A matrix directly, or
#'   - A list containing the scaled matrix under various possible names
#'
#' This wrapper handles both cases and extracts the scaled matrix.
#'
#' If a column has zero range (all values identical), the scaled values will
#' be 0/0 = NaN. The C++ function should handle this gracefully, but users
#' should be aware that constant output dimensions are not useful for OSFD.
#'
#' @export
osfd_scale_Y <- function(Y) {

  # Get the C++ scaling function from OSFD package
  osfd_scale_cpp <- get_osfd_fn("scale_cpp")

  # Call the C++ function to perform min-max scaling
  out <- osfd_scale_cpp(Y)

  # Handle return value based on what the C++ function returns:

  # Case 1: Returns a matrix directly (expected behavior)
  if (is.matrix(out)) {
    return(out)
  }

  # Case 2: Returns a list containing the scaled matrix
  # The matrix might be stored under various names depending on the OSFD version
  if (is.list(out)) {
    # Check common field names used in different OSFD versions
    for (nm in c("Y", "X", "Ys", "Xs", "Y_sc", "X_sc")) {
      if (!is.null(out[[nm]]) && is.matrix(out[[nm]])) {
        return(out[[nm]])
      }
    }
  }

  # If we get here, the return value is neither a matrix nor a list with a matrix
  # This should not happen with the current OSFD package, but we error informatively
  stop(
    "OSFD::scale_cpp returned an unexpected type; ",
    "inspect it with str(scale_cpp(Y))."
  )
}
