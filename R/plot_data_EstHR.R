#' Values plotted by plot.EstHR()
#'
#' @param x An object of class \code{"EstHR"} returned by \code{EstHR()}.
#' @param labels \code{NULL} for the default labels of \code{plot_labels()},
#'   or one label per method.
#' @param value_digits Number of significant digits of the values written on
#'   the right axis, a positive integer. Digits before the decimal point are
#'   all written, and a value without decimals has no trailing point.
#'
#' @return A data frame with one row per method, in the order of the methods
#'   of \code{x}, and the columns \code{method}, \code{label}, \code{HR},
#'   \code{lower}, \code{upper}, \code{text} (the hazard ratio and the
#'   confidence limits as written on the right axis, or \code{"NA"}), and
#'   \code{y} (vertical position, from the number of methods at the top down
#'   to 1 at the bottom). The attribute \code{cox_hr} is the Cox estimate of
#'   the hazard ratio (from \code{"Cox"}, or \code{"Cox_NR"} without
#'   \code{"Cox"}), and is absent when there is no Cox estimate.
#'
#' @noRd
plot_data_EstHR <- function(x, labels = NULL, value_digits = 3L) {
  method <- attr(x, "method")
  if (is.null(method) || length(x) != 5L * length(method)) {
    stop("'x' must be the output of EstHR()")
  }
  cols <- if (is.matrix(x)) colnames(x) else names(x)
  est <- matrix(as.vector(x), nrow = length(method),
                dimnames = list(method, cols))
  k <- length(method)
  # unname() also for one method, where est[, j] is a single number
  hr <- unname(est[, 2L])
  lower <- unname(est[, 4L])
  upper <- unname(est[, 5L])
  ok <- is.finite(hr)
  if (!any(ok)) stop("there is no estimate to plot")
  if (!is.numeric(value_digits) || length(value_digits) != 1L ||
      !is.finite(value_digits) || value_digits < 1 ||
      value_digits != round(value_digits)) {
    stop("'value_digits' must be a positive integer")
  }

  # formatC() pads its result with blanks, which trimws() removes, and with
  # flag "#" it ends a value without decimals (such as "1234.") with a point,
  # which sub() removes
  fmt <- function(v) {
    s <- trimws(formatC(v, digits = value_digits, format = "fg", flag = "#"))
    sub("\\.$", "", s)
  }
  text <- ifelse(ok, paste0(fmt(hr), " (", fmt(lower), ", ", fmt(upper),
                            ")"), "NA")
  out <- data.frame(method = method, label = plot_labels(method, labels),
                    HR = hr, lower = lower, upper = upper, text = text,
                    y = rev(seq_len(k)), stringsAsFactors = FALSE)
  cox <- intersect(c("Cox", "Cox_NR"), method[ok])
  if (length(cox) > 0L) attr(out, "cox_hr") <- hr[match(cox[1L], method)]
  out
}
