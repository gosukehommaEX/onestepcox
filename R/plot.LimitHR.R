#' Plot the large-sample bias of hazard ratio estimators
#'
#' Plots the large-sample bias of the log hazard ratio estimators computed by
#' \code{\link{LimitHR}} against the true hazard ratio, on a logarithmic
#' horizontal axis, with one line per method.
#'
#' @param x An object of class \code{"LimitHR"} returned by
#'   \code{LimitHR()}.
#' @param methods Methods to plot, in the order of the legend. Default: the
#'   methods whose limit differs from the true value, \code{"Pike"},
#'   \code{"Peto"}, \code{"Pike_Newton"}, \code{"Pike_Chebyshev"}, and
#'   \code{"Pike_Halley"}.
#' @param scale \code{"abs"} (default) for the absolute bias on a
#'   logarithmic vertical axis, or \code{"signed"} for the bias on a linear
#'   axis with a dashed line at 0.
#' @param labels Labels of the methods in the legend, one per method of
#'   \code{methods}. Default \code{NULL}: the names of the methods with
#'   underscores replaced by hyphens, such as \code{"Pike-Halley"}.
#' @param col Colours of the methods, recycled. Default \code{NULL} (the
#'   Okabe-Ito palette).
#' @param lty,lwd,pch,cex Line type, line width, plotting symbol, and symbol
#'   size of the methods, recycled.
#' @param base_size Character expansion of the axis labels, axis titles, and
#'   legend (the title is \code{1.2 * base_size}). Default 1.
#' @param legend_pos Position of the legend: a keyword of
#'   \code{graphics::legend()} such as \code{"top"} (default) or
#'   \code{"bottomright"}, numeric coordinates \code{c(x, y)} of its top left
#'   corner, or \code{"none"}.
#' @param legend_inset Inset of the legend from the margins as a fraction of
#'   the plot region. Default 0.02.
#' @param legend_ncol Number of columns of the legend. Default 1.
#' @param xlab,ylab,main Axis titles and main title. \code{ylab = NULL} gives
#'   a title that describes \code{scale}.
#' @param xlim,ylim Limits of the axes. Default \code{NULL} (the range of the
#'   values).
#' @param floor Smallest absolute bias shown with \code{scale = "abs"};
#'   smaller values (including 0 at a hazard ratio of 1) are drawn at
#'   \code{floor}. Default \code{1e-16}.
#' @param ... Further graphical parameters passed to \code{graphics::plot()},
#'   such as \code{las}.
#'
#' @return A matrix of the bias (signed, before \code{floor} is applied) with
#'   one row per hazard ratio, in increasing order, and one column per
#'   method, invisibly. The graphical parameters are restored on exit.
#'
#' @seealso \code{\link{LimitHR}}, \code{\link{print.LimitHR}}.
#'
#' @examples
#' lim <- LimitHR(HR = exp(seq(log(0.1), log(10), length.out = 21)),
#'                hazard_C = log(2) / 18, n_T = 100, n_C = 100, a = 12,
#'                E = 160)
#' plot(lim)
#' plot(lim, scale = "signed", methods = c("Pike", "Peto", "MH_Halley"),
#'      legend_pos = "topright", pch = c(1, 2, 4), lwd = 3, base_size = 1.2)
#'
#' @export
plot.LimitHR <- function(x, methods = c("Pike", "Peto", "Pike_Newton",
                                        "Pike_Chebyshev", "Pike_Halley"),
                         scale = c("abs", "signed"), labels = NULL,
                         col = NULL, lty = 1, lwd = 2, pch = 19, cex = 0.8,
                         base_size = 1, legend_pos = "top",
                         legend_inset = 0.02, legend_ncol = 1L,
                         xlab = "Hazard ratio",
                         ylab = NULL, main = NULL, xlim = NULL, ylim = NULL,
                         floor = 1e-16, ...) {
  scale <- match.arg(scale)
  if (!all(c("HR", "method", "bias") %in% names(x))) {
    stop("'x' must have the columns HR, method, and bias of LimitHR()")
  }
  if (length(methods) == 0L || !all(methods %in% x$method)) {
    stop("'methods' must be methods of 'x': ",
         paste(unique(x$method), collapse = ", "))
  }
  hr <- sort(unique(x$HR))
  k <- length(methods)
  mat <- matrix(NA_real_, nrow = length(hr), ncol = k,
                dimnames = list(as.character(signif(hr, 6)), methods))
  for (j in seq_len(k)) {
    sel <- x$method == methods[j]
    mat[, j] <- x$bias[sel][match(hr, x$HR[sel])]
  }
  if (!any(is.finite(mat))) stop("there is no value to plot")
  labels <- plot_labels(methods, labels)
  col <- plot_colors(k, col)
  lty <- rep_len(lty, k)
  lwd <- rep_len(lwd, k)
  pch <- rep_len(pch, k)

  y <- mat
  if (scale == "abs") {
    y <- abs(y)
    y[which(y < floor)] <- floor
  }
  if (is.null(ylim)) {
    ylim <- range(y, finite = TRUE)
    if (scale == "signed") ylim <- range(c(ylim, 0))
    if (ylim[1L] == ylim[2L]) {
      ylim <- if (scale == "abs") ylim * c(0.1, 10) else ylim + c(-0.5, 0.5)
    }
  }
  if (is.null(xlim)) xlim <- range(hr)
  if (is.null(ylab)) {
    ylab <- if (scale == "abs") {
      "|Bias| of the log hazard ratio"
    } else {
      "Bias of the log hazard ratio"
    }
  }

  op <- graphics::par(cex.axis = base_size, cex.lab = base_size,
                      cex.main = 1.2 * base_size)
  on.exit(graphics::par(op))
  graphics::plot(NA, xlim = xlim, ylim = ylim,
                 log = if (scale == "abs") "xy" else "x", xlab = xlab,
                 ylab = ylab, main = main, ...)
  if (scale == "signed") graphics::abline(h = 0, col = "grey50", lty = 2)
  for (j in seq_len(k)) {
    graphics::lines(hr, y[, j], type = "b", col = col[j], lty = lty[j],
                    lwd = lwd[j], pch = pch[j], cex = cex)
  }
  draw_legend(legend_pos, labels, inset = legend_inset, ncol = legend_ncol,
              cex = base_size, col = col, lty = lty, lwd = lwd, pch = pch)
  invisible(mat)
}
