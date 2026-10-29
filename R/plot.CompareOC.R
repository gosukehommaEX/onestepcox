#' Plot operating characteristics of hazard ratio estimators
#'
#' Plots a result of \code{\link{CompareOC}} with \code{ggplot2}. With
#' \code{what = "deviation"}, box plots of the absolute deviation of each
#' method from the reference estimator in the replicates, grouped by
#' scenario, on a logarithmic axis. With the other values of \code{what}, a
#' summary of each method (bias, RMSE, coverage, median ratio of the standard
#' errors to that of the reference, or computing time per replicate) across
#' the scenarios, as points joined by lines, with a dashed reference line at
#' the target value (0 for the bias, \code{conf.level} for the coverage, and
#' 1 for the ratio of the standard errors).
#'
#' @param x An object of class \code{"CompareOC"} returned by
#'   \code{CompareOC()}.
#' @param what Quantity to plot: \code{"deviation"} (default), \code{"bias"},
#'   \code{"rmse"}, \code{"coverage"}, \code{"se_ratio"}, or \code{"time"}.
#' @param methods Methods to plot, in the order of the legend. Default
#'   \code{NULL}: all methods, except the reference for \code{"deviation"}
#'   and \code{"se_ratio"}.
#' @param labels Labels of the methods in the legend, one per method of
#'   \code{methods}. Default \code{NULL}: the names of the methods with
#'   underscores replaced by hyphens, such as \code{"MH-Halley"}.
#' @param col Colours of the methods, recycled. Default \code{NULL} (the
#'   Okabe-Ito palette).
#' @param shape Point shapes of the methods, with the numbers of \code{pch}
#'   of base graphics, recycled; for \code{"deviation"}, the shapes of the
#'   outlying points of the box plots. Default 19 (filled circle).
#' @param size Size of the points. Default 2.
#' @param linetype Line types of the methods, recycled (not used for
#'   \code{"deviation"}). Default 1 (solid).
#' @param linewidth Width of the lines, and of the outlines of the box plots.
#'   Default 0.8.
#' @param base_size Base font size in points, as in
#'   \code{ggplot2::theme_bw()}. Default 11.
#' @param legend_pos Position of the legend: \code{"bottom"} (default),
#'   \code{"top"}, \code{"right"}, or \code{"left"} (outside the plot
#'   panel), \code{"none"}, or a numeric vector \code{c(x, y)} with elements
#'   between 0 and 1 for a position inside the panel (\code{c(1, 1)} is the
#'   top right corner).
#' @param legend_ncol Number of columns of the legend. Default \code{NULL}:
#'   at most 3 for a legend at the bottom or top, and 1 otherwise.
#' @param xlab,ylab,main Axis titles and main title. \code{ylab = NULL} gives
#'   a title that describes \code{what}.
#' @param ylim Limits of the vertical axis, two increasing numbers (positive
#'   on a logarithmic axis). Default \code{NULL} (the range of the values and
#'   of the reference line). Values beyond the limits are kept in the box
#'   plots and cut at the edge of the panel.
#' @param floor Smallest deviation shown with \code{"deviation"}, a positive
#'   number; smaller deviations (including 0) are drawn at \code{floor}.
#'   Default \code{1e-16}.
#' @param ... Not used.
#'
#' @return A \code{ggplot} object. It is drawn when it is printed, which
#'   happens automatically at the top level; inside functions and loops, call
#'   \code{print()} or save it with \code{ggplot2::ggsave()}. It can be
#'   modified further with \code{ggplot2}. The vertical axis is logarithmic
#'   for \code{"deviation"}, and for \code{"time"} when all times are
#'   positive.
#'
#' @seealso \code{\link{CompareOC}}, \code{\link{print.CompareOC}}.
#'
#' @examples
#' scenarios <- data.frame(label = c("HR 0.5", "HR 2"), n_T = 100,
#'                         n_C = 100, hazard_C = log(2) / 18,
#'                         HR = c(0.5, 2), a = 12, E = 110)
#' res <- CompareOC(scenarios, c("Cox", "Peto", "Pike", "MH_Halley"),
#'                  nsim = 30, seed = 1)
#' plot(res)
#' plot(res, what = "bias", legend_pos = "right", base_size = 13)
#' plot(res, what = "coverage", methods = c("Cox", "MH_Halley"),
#'      col = c("black", "red"), linetype = c(1, 2), ylim = c(0.8, 1),
#'      legend_pos = c(0.02, 0.02))
#'
#' @export
plot.CompareOC <- function(x, what = c("deviation", "bias", "rmse",
                                       "coverage", "se_ratio", "time"),
                           methods = NULL, labels = NULL, col = NULL,
                           shape = 19, size = 2, linetype = 1,
                           linewidth = 0.8, base_size = 11,
                           legend_pos = "bottom", legend_ncol = NULL,
                           xlab = "Scenario", ylab = NULL, main = NULL,
                           ylim = NULL, floor = 1e-16, ...) {
  what <- match.arg(what)
  d <- plot_data_CompareOC(x, what, methods, labels, floor)
  ref <- attr(d, "reference")
  methods <- attr(d, "methods")
  lab <- attr(d, "labels")
  k <- length(methods)
  d$scenario <- factor(d$scenario, levels = attr(d, "scenarios"))
  d$method <- factor(d$method, levels = methods)

  if (what == "deviation") {
    if (is.null(ylab)) ylab <- paste0("|coef - coef of ", ref, "|")
    aes_used <- c("colour", "shape")
    log_y <- TRUE
    p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$scenario, y = .data$value,
                                         colour = .data$method,
                                         shape = .data$method)) +
      # outlier.shape = NULL takes the shape of the method (the default of
      # ggplot2 3.5 is a fixed shape)
      ggplot2::geom_boxplot(fill = "white", linewidth = linewidth,
                            outlier.shape = NULL, outlier.size = size) +
      ggplot2::scale_y_log10()
  } else {
    target <- switch(what, bias = 0, coverage = attr(x, "conf.level"),
                     se_ratio = 1, NULL)
    if (is.null(ylab)) {
      ylab <- switch(what, bias = "Bias of the log hazard ratio",
                     rmse = "RMSE of the log hazard ratio",
                     coverage = "Coverage of the confidence interval",
                     se_ratio = paste0("Median ratio of se to that of ", ref),
                     time = "Seconds per replicate")
    }
    aes_used <- c("colour", "shape", "linetype")
    dd <- d[is.finite(d$value), , drop = FALSE]
    p <- ggplot2::ggplot(dd, ggplot2::aes(x = .data$scenario,
                                          y = .data$value,
                                          colour = .data$method,
                                          shape = .data$method,
                                          linetype = .data$method,
                                          group = .data$method))
    if (!is.null(target)) {
      p <- p + ggplot2::geom_hline(yintercept = target, colour = "grey50",
                                   linetype = "dashed")
    }
    # Lines only when some method has two or more points to join
    if (anyDuplicated(as.character(dd$method)) > 0L) {
      p <- p + ggplot2::geom_line(linewidth = linewidth)
    }
    p <- p + ggplot2::geom_point(size = size)
    log_y <- what == "time" && all(dd$value > 0)
    if (log_y) p <- p + ggplot2::scale_y_log10()
  }
  check_limits(ylim, "ylim", positive = log_y)
  p +
    ggplot2::scale_x_discrete(drop = FALSE) +
    plot_scales(methods, lab, col = col, shape = shape, linetype = linetype,
                aesthetics = aes_used) +
    ggplot2::coord_cartesian(ylim = ylim) +
    ggplot2::labs(x = xlab, y = ylab, title = main) +
    plot_theme(base_size) +
    plot_legend(legend_pos, legend_ncol, k, aesthetics = aes_used)
}
