#' Plot hazard ratio estimates
#'
#' Draws a forest plot of the estimates of \code{\link{EstHR}} with
#' \code{ggplot2}: for each method, the hazard ratio (point) and its
#' confidence interval (segment) on a logarithmic axis, with the methods from
#' top to bottom in the order of \code{method}. Vertical reference lines mark
#' the hazard ratio 1 and, when a Cox method is among the methods, the Cox
#' estimate. The estimate and interval of each method are written on the
#' right axis.
#'
#' @param x An object of class \code{"EstHR"} returned by \code{EstHR()}.
#' @param labels Labels of the methods on the vertical axis and in the
#'   legend, one per method. Default \code{NULL}: the names of the methods
#'   with underscores replaced by hyphens, such as \code{"MH-Halley"}.
#' @param col Colours of the methods, recycled. Default \code{NULL} (the
#'   Okabe-Ito palette).
#' @param shape Point shapes of the methods, with the numbers of \code{pch}
#'   of base graphics, recycled. Default 15 (filled square).
#' @param size Size of the points. Default 3.
#' @param linetype Line types of the interval segments, recycled. Default 1
#'   (solid).
#' @param linewidth Width of the interval segments. Default 0.8.
#' @param base_size Base font size in points, as in
#'   \code{ggplot2::theme_bw()}. Default 11.
#' @param show_values Logical. Write the hazard ratio and the confidence
#'   limits on the right axis. Default \code{TRUE}.
#' @param value_digits Number of significant digits of the values on the
#'   right axis. Default 3.
#' @param ref_lines Logical. Draw the reference lines at the hazard ratio 1
#'   (solid) and at the Cox estimate (dashed; from \code{"Cox"}, or
#'   \code{"Cox_NR"} without \code{"Cox"}). Default \code{TRUE}.
#' @param ref_col Colour of the reference lines. Default \code{"grey50"}.
#' @param legend_pos Position of a legend of the methods: \code{"none"}
#'   (default, as the methods label the vertical axis), \code{"bottom"},
#'   \code{"top"}, \code{"right"}, or \code{"left"} (outside the plot panel),
#'   or a numeric vector \code{c(x, y)} with elements between 0 and 1 for a
#'   position inside the panel (\code{c(1, 1)} is the top right corner).
#' @param legend_ncol Number of columns of the legend. Default \code{NULL}:
#'   at most 3 for a legend at the bottom or top, and 1 otherwise.
#' @param xlab,main Axis title and main title.
#' @param xlim Limits of the horizontal axis (hazard ratio), two increasing
#'   positive numbers. Default \code{NULL}: the range of the estimates, the
#'   intervals, and the reference line at 1, widened by 4 per cent on the
#'   logarithmic scale. Intervals beyond the limits are cut at the edge of
#'   the panel.
#' @param ... Not used.
#'
#' @return A \code{ggplot} object. It is drawn when it is printed, which
#'   happens automatically at the top level; inside functions and loops, call
#'   \code{print()} or save it with \code{ggplot2::ggsave()}. It can be
#'   modified further with \code{ggplot2}, for example by adding
#'   \code{ggplot2::theme()}.
#'
#' @seealso \code{\link{EstHR}}, \code{\link{print.EstHR}}.
#'
#' @examples
#' est <- EstHR(Surv(rfstime, status) ~ hormon + strata(grade),
#'              data = survival::gbsg,
#'              method = c("MH_Halley", "MH", "Pike", "Peto", "Cox"))
#' plot(est)
#' plot(est, col = "black", shape = 19, linewidth = 1.2, base_size = 13,
#'      main = "GBSG: hormonal therapy", show_values = FALSE)
#' plot(est, legend_pos = "bottom", xlim = c(0.5, 1))
#' # Further changes with ggplot2
#' plot(est) + ggplot2::theme(axis.text = ggplot2::element_text(size = 12))
#'
#' @export
plot.EstHR <- function(x, labels = NULL, col = NULL, shape = 15, size = 3,
                       linetype = 1, linewidth = 0.8, base_size = 11,
                       show_values = TRUE, value_digits = 3L,
                       ref_lines = TRUE, ref_col = "grey50",
                       legend_pos = "none", legend_ncol = NULL,
                       xlab = "Hazard ratio", main = NULL, xlim = NULL, ...) {
  d <- plot_data_EstHR(x, labels, value_digits)
  cox_hr <- attr(d, "cox_hr")
  k <- nrow(d)
  methods <- d$method

  check_limits(xlim, "xlim", positive = TRUE)
  if (is.null(xlim)) {
    v <- c(d$HR, d$lower, d$upper, if (ref_lines) 1)
    v <- v[is.finite(v) & v > 0]
    xlim <- range(v)
    if (xlim[1L] == xlim[2L]) xlim <- xlim * c(0.8, 1.25)
    pad <- 0.04 * diff(log(xlim))
    xlim <- exp(log(xlim) + c(-pad, pad))
  }

  # Intervals are cut at the limits, so that no value is infinite on the
  # logarithmic axis
  d$lo <- pmin(pmax(d$lower, xlim[1L]), xlim[2L])
  d$up <- pmin(pmax(d$upper, xlim[1L]), xlim[2L])
  d$method <- factor(d$method, levels = methods)
  pts <- d[is.finite(d$HR) & d$HR > 0, , drop = FALSE]
  seg <- pts[is.finite(pts$lo) & is.finite(pts$up), , drop = FALSE]

  p <- ggplot2::ggplot()
  if (ref_lines) {
    p <- p + ggplot2::geom_vline(xintercept = 1, colour = ref_col,
                                 linetype = "solid")
    if (!is.null(cox_hr)) {
      p <- p + ggplot2::geom_vline(xintercept = cox_hr, colour = ref_col,
                                   linetype = "dashed")
    }
  }
  if (nrow(seg) > 0L) {
    p <- p + ggplot2::geom_segment(
      data = seg,
      mapping = ggplot2::aes(x = .data$lo, xend = .data$up, y = .data$y,
                             yend = .data$y, colour = .data$method,
                             linetype = .data$method),
      linewidth = linewidth
    )
  }
  p <- p +
    ggplot2::geom_point(
      data = pts,
      mapping = ggplot2::aes(x = .data$HR, y = .data$y,
                             colour = .data$method, shape = .data$method),
      size = size
    ) +
    ggplot2::scale_x_log10() +
    ggplot2::scale_y_continuous(
      breaks = d$y, labels = d$label,
      sec.axis = if (show_values) {
        ggplot2::sec_axis(function(v) v, breaks = d$y, labels = d$text)
      } else {
        ggplot2::waiver()
      }
    ) +
    plot_scales(methods, d$label, col = col, shape = shape,
                linetype = linetype) +
    ggplot2::coord_cartesian(xlim = xlim, ylim = c(0.5, k + 0.5),
                             expand = FALSE) +
    ggplot2::labs(x = xlab, y = NULL, title = main) +
    plot_theme(base_size) +
    plot_legend(legend_pos, legend_ncol, k) +
    ggplot2::theme(axis.ticks.y = ggplot2::element_blank(),
                   panel.grid.major.y = ggplot2::element_blank())
  p
}
