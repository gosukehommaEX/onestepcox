#' Plot the large-sample bias of hazard ratio estimators
#'
#' Plots the large-sample bias of the log hazard ratio estimators computed by
#' \code{\link{LimitHR}} against the true hazard ratio with \code{ggplot2},
#' on a logarithmic horizontal axis, with one line per method.
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
#' @param col Colors of the methods, recycled. Default \code{NULL} (the
#'   Okabe-Ito palette).
#' @param shape Point shapes of the methods, with the numbers of \code{pch}
#'   of base graphics, recycled. Default \code{c(16, 17, 15, 18, 1, 2, 0, 5)}
#'   (filled circle, triangle, square, and diamond, then the open ones), so
#'   that methods with nearly equal values remain distinguishable.
#' @param size Size of the points. Default 1.5.
#' @param linetype Line types of the methods, recycled. Default
#'   \code{c("solid", "dashed", "dotted", "dotdash", "longdash", "twodash")}.
#' @param linewidth Width of the lines. Default 0.8.
#' @param base_size Base font size in points, as in
#'   \code{ggplot2::theme_bw()}. Default 11.
#' @param legend_pos Position of the legend: \code{"bottom"} (default),
#'   \code{"top"}, \code{"right"}, or \code{"left"} (outside the plot
#'   panel), \code{"none"}, or a numeric vector \code{c(x, y)} with elements
#'   between 0 and 1 for a position inside the panel (\code{c(1, 1)} is the
#'   top right corner).
#' @param legend_ncol Number of columns of the legend, filled row by row in
#'   the order of the methods. Default \code{NULL}: at most 3 for a legend at
#'   the bottom or top, and 1 otherwise.
#' @param xlab,ylab,main Axis titles and main title. \code{ylab = NULL} gives
#'   a title that describes \code{scale}.
#' @param xlim,ylim Limits of the axes, two increasing numbers (positive on
#'   a logarithmic axis: always for \code{xlim}, and for \code{ylim} with
#'   \code{scale = "abs"}). Default \code{NULL} (the range of the values).
#' @param floor Smallest absolute bias shown with \code{scale = "abs"}, a
#'   positive number; smaller values (including 0 at a hazard ratio of 1) are
#'   drawn at \code{floor}. Default \code{1e-16}.
#' @param ... Not used.
#'
#' @return A \code{ggplot} object. It is drawn when it is printed, which
#'   happens automatically at the top level; inside functions and loops, call
#'   \code{print()} or save it with \code{ggplot2::ggsave()}. It can be
#'   modified further with \code{ggplot2}.
#'
#' @seealso \code{\link{LimitHR}}, \code{\link{print.LimitHR}}.
#'
#' @examples
#' lim <- LimitHR(HR = exp(seq(log(0.1), log(10), length.out = 21)),
#'                hazard_C = log(2) / 18, n_T = 100, n_C = 100, a = 12,
#'                E = 160)
#' plot(lim)
#' plot(lim, scale = "signed", methods = c("Pike", "Peto", "MH_Halley"),
#'      legend_pos = "right", shape = c(1, 2, 4), linewidth = 1.2,
#'      base_size = 13)
#'
#' @export
plot.LimitHR <- function(x, methods = c("Pike", "Peto", "Pike_Newton",
                                        "Pike_Chebyshev", "Pike_Halley"),
                         scale = c("abs", "signed"), labels = NULL,
                         col = NULL, shape = c(16, 17, 15, 18, 1, 2, 0, 5),
                         size = 1.5,
                         linetype = c("solid", "dashed", "dotted", "dotdash",
                                      "longdash", "twodash"),
                         linewidth = 0.8, base_size = 11,
                         legend_pos = "bottom", legend_ncol = NULL,
                         xlab = "Hazard ratio", ylab = NULL, main = NULL,
                         xlim = NULL, ylim = NULL, floor = 1e-16, ...) {
  scale <- match.arg(scale)
  check_limits(xlim, "xlim", positive = TRUE)
  check_limits(ylim, "ylim", positive = scale == "abs")
  d <- plot_data_LimitHR(x, methods, scale, labels, floor)
  methods <- attr(d, "methods")
  lab <- attr(d, "labels")
  k <- length(methods)
  if (is.null(ylab)) {
    ylab <- if (scale == "abs") {
      "|Bias| of the log hazard ratio"
    } else {
      "Bias of the log hazard ratio"
    }
  }
  d$method <- factor(d$method, levels = methods)
  dd <- d[is.finite(d$value), , drop = FALSE]

  p <- ggplot2::ggplot(dd, ggplot2::aes(x = .data$HR, y = .data$value,
                                        colour = .data$method,
                                        shape = .data$method,
                                        linetype = .data$method,
                                        group = .data$method))
  if (scale == "signed") {
    p <- p + ggplot2::geom_hline(yintercept = 0, colour = "grey50",
                                 linetype = "dashed")
  }
  # Lines only when some method has two or more points to join
  if (anyDuplicated(as.character(dd$method)) > 0L) {
    p <- p + ggplot2::geom_line(linewidth = linewidth)
  }
  p <- p + ggplot2::geom_point(size = size) + ggplot2::scale_x_log10()
  if (scale == "abs") p <- p + ggplot2::scale_y_log10()
  p +
    plot_scales(methods, lab, col = col, shape = shape, linetype = linetype) +
    ggplot2::coord_cartesian(xlim = xlim, ylim = ylim) +
    ggplot2::labs(x = xlab, y = ylab, title = main) +
    plot_theme(base_size) +
    plot_legend(legend_pos, legend_ncol, k)
}
