#' Check the axis limits of a plot method
#'
#' @param lim \code{NULL}, or the limits given by the user.
#' @param name Name of the argument, for the error message.
#' @param positive Logical. The axis is logarithmic, so that the limits must
#'   be positive.
#'
#' @return \code{lim}, invisibly, when it is \code{NULL} or two finite
#'   increasing numbers (positive when \code{positive} is \code{TRUE});
#'   otherwise an error.
#'
#' @noRd
check_limits <- function(lim, name, positive = FALSE) {
  if (is.null(lim)) return(invisible(lim))
  ok <- is.numeric(lim) && length(lim) == 2L && all(is.finite(lim)) &&
    lim[1L] < lim[2L] && (!positive || lim[1L] > 0)
  if (!ok) {
    stop("'", name, "' must be two increasing ",
         if (positive) "positive " else "", "numbers",
         if (positive) " (the axis is logarithmic)" else "")
  }
  invisible(lim)
}
