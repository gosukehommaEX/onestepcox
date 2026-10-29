#' Values plotted by plot.LimitHR()
#'
#' @param x An object of class \code{"LimitHR"} returned by
#'   \code{LimitHR()}.
#' @param methods Methods to plot, in the order of the legend.
#' @param scale \code{"abs"} for the absolute bias, or \code{"signed"} for
#'   the bias.
#' @param labels \code{NULL} for the default labels of \code{plot_labels()},
#'   or one label per method of \code{methods}.
#' @param floor Smallest absolute bias kept with \code{scale = "abs"}, a
#'   positive number; smaller values (including 0) are set to \code{floor}.
#'
#' @return A data frame with one row per method and hazard ratio (hazard
#'   ratios in increasing order, varying fastest) and the columns \code{HR},
#'   \code{method}, \code{bias} (signed, before \code{floor} is applied),
#'   \code{value} (the plotted value), and \code{label}. The attributes are
#'   \code{methods} and \code{labels} (one per method).
#'
#' @noRd
plot_data_LimitHR <- function(x, methods, scale, labels = NULL,
                              floor = 1e-16) {
  if (!all(c("HR", "method", "bias") %in% names(x))) {
    stop("'x' must have the columns HR, method, and bias of LimitHR()")
  }
  if (length(methods) == 0L || !all(methods %in% x$method)) {
    stop("'methods' must be methods of 'x': ",
         paste(unique(x$method), collapse = ", "))
  }
  # The absolute bias is drawn on a logarithmic axis
  if (!is.numeric(floor) || length(floor) != 1L || !is.finite(floor) ||
      floor <= 0) {
    stop("'floor' must be a single positive number")
  }
  hr <- sort(unique(x$HR))
  lab <- plot_labels(methods, labels)
  out <- do.call(rbind, lapply(methods, function(m) {
    sel <- x$method == m
    data.frame(HR = hr, method = m,
               bias = x$bias[sel][match(hr, x$HR[sel])],
               stringsAsFactors = FALSE)
  }))
  out$value <- if (scale == "abs") pmax(abs(out$bias), floor) else out$bias
  if (!any(is.finite(out$value))) stop("there is no value to plot")
  out$label <- lab[match(out$method, methods)]
  rownames(out) <- NULL
  attr(out, "methods") <- methods
  attr(out, "labels") <- lab
  out
}
