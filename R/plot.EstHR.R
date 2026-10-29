#' Plot hazard ratio estimates
#'
#' Draws a forest plot of the estimates of \code{\link{EstHR}}: for each
#' method, the hazard ratio (point) and its confidence interval (segment) on a
#' logarithmic axis, with the methods from top to bottom in the order of
#' \code{method}. Vertical reference lines mark the hazard ratio 1 and, when
#' a Cox method is among the methods, the Cox estimate. The estimate and
#' interval of each method are written in the right margin.
#'
#' @param x An object of class \code{"EstHR"} returned by \code{EstHR()}.
#' @param labels Labels of the methods on the vertical axis and in the
#'   legend, one per method. Default \code{NULL}: the names of the methods
#'   with underscores replaced by hyphens, such as \code{"MH-Halley"}.
#' @param col Colours of the methods, recycled. Default \code{NULL} (the
#'   Okabe-Ito palette).
#' @param pch,cex Plotting symbol and its size for the hazard ratios.
#' @param lwd,lty Width and type of the interval segments.
#' @param base_size Character expansion of the axis labels, axis titles,
#'   values, and legend (the title is \code{1.2 * base_size}). Default 1.
#' @param show_values Logical. Write the hazard ratio and the confidence
#'   limits in the right margin. Default \code{TRUE}.
#' @param value_digits Number of significant digits of the values in the
#'   right margin. Default 3.
#' @param ref_lines Logical. Draw the reference lines at the hazard ratio 1
#'   (solid) and at the Cox estimate (dashed; from \code{"Cox"}, or
#'   \code{"Cox_NR"} without \code{"Cox"}). Default \code{TRUE}.
#' @param ref_col Colour of the reference lines. Default \code{"grey50"}.
#' @param legend_pos Position of a legend of the methods: \code{"none"}
#'   (default, as the methods label the vertical axis), a keyword of
#'   \code{graphics::legend()} such as \code{"topright"} or
#'   \code{"bottomleft"}, or numeric coordinates \code{c(x, y)} of its top
#'   left corner.
#' @param legend_inset Inset of the legend from the margins as a fraction of
#'   the plot region. Default 0.02.
#' @param legend_ncol Number of columns of the legend. Default 1.
#' @param xlab,main Axis title and main title.
#' @param xlim Limits of the horizontal axis (hazard ratio). Default
#'   \code{NULL} (the range of the intervals and 1).
#' @param mar Margins in lines, as in \code{graphics::par()}. Default
#'   \code{NULL} (wide enough for the method names and the values).
#' @param ... Further graphical parameters passed to \code{graphics::plot()},
#'   such as \code{las} or \code{family}.
#'
#' @return A data frame of the plotted values, with the columns
#'   \code{method}, \code{label}, \code{HR}, \code{lower}, \code{upper},
#'   \code{text} (the values written in the right margin), and \code{y} (the
#'   vertical position), invisibly. The graphical parameters are restored on
#'   exit.
#'
#' @seealso \code{\link{EstHR}}, \code{\link{print.EstHR}}.
#'
#' @examples
#' est <- EstHR(Surv(rfstime, status) ~ hormon + strata(grade),
#'              data = survival::gbsg,
#'              method = c("MH_Halley", "MH", "Pike", "Peto", "Cox"))
#' plot(est)
#' plot(est, col = "black", pch = 19, lwd = 3, base_size = 1.2,
#'      main = "GBSG: hormonal therapy", show_values = FALSE)
#' plot(est, legend_pos = "topright", xlim = c(0.5, 1))
#'
#' @export
plot.EstHR <- function(x, labels = NULL, col = NULL, pch = 15, cex = 1.2,
                       lwd = 2, lty = 1, base_size = 1, show_values = TRUE,
                       value_digits = 3L, ref_lines = TRUE,
                       ref_col = "grey50", legend_pos = "none",
                       legend_inset = 0.02, legend_ncol = 1L,
                       xlab = "Hazard ratio", main = NULL, xlim = NULL,
                       mar = NULL, ...) {
  method <- attr(x, "method")
  if (is.null(method) || length(x) != 5L * length(method)) {
    stop("'x' must be the output of EstHR()")
  }
  cols <- if (is.matrix(x)) colnames(x) else names(x)
  est <- matrix(as.vector(x), nrow = length(method),
                dimnames = list(method, cols))
  k <- length(method)
  # Named by method also for one method (a 1 x 1 result loses its names)
  hr <- stats::setNames(est[, 2L], method)
  lower <- stats::setNames(est[, 4L], method)
  upper <- stats::setNames(est[, 5L], method)
  ok <- is.finite(hr)
  if (!any(ok)) stop("there is no estimate to plot")
  y <- rev(seq_len(k))
  labels <- plot_labels(method, labels)
  col <- plot_colors(k, col)
  pch <- rep_len(pch, k)
  lwd <- rep_len(lwd, k)
  lty <- rep_len(lty, k)

  # Reference value of the Cox estimate
  cox <- intersect(c("Cox", "Cox_NR"), method[ok])
  cox_hr <- if (length(cox) > 0L) hr[[cox[1L]]] else NULL

  if (is.null(xlim)) {
    xlim <- range(c(lower, upper, hr, 1), finite = TRUE)
    pad <- 0.04 * diff(log(xlim))
    xlim <- exp(log(xlim) + c(-pad, pad))
  }
  # formatC() pads its result with blanks, which trimws() removes
  fmt <- function(v) {
    trimws(formatC(v, digits = value_digits, format = "fg", flag = "#"))
  }
  values <- ifelse(ok, paste0(fmt(hr), " (", fmt(lower), ", ", fmt(upper),
                              ")"), "NA")
  if (is.null(mar)) {
    # About 0.45 lines per character of the labels at cex = 1
    left <- 0.45 * base_size * max(nchar(labels)) + 1.5
    right <- if (show_values) {
      0.45 * base_size * max(nchar(values)) + 1.5
    } else {
      1
    }
    mar <- c(4.1, left, if (is.null(main)) 1.5 else 3.5, right)
  }
  op <- graphics::par(mar = mar)
  on.exit(graphics::par(op))

  graphics::plot(NA, xlim = xlim, ylim = c(0.5, k + 0.5), log = "x",
                 yaxt = "n", xlab = xlab, ylab = "", main = main,
                 cex.axis = base_size, cex.lab = base_size,
                 cex.main = 1.2 * base_size, ...)
  if (ref_lines) {
    graphics::abline(v = 1, col = ref_col, lty = 1)
    if (!is.null(cox_hr)) graphics::abline(v = cox_hr, col = ref_col, lty = 2)
  }
  # Intervals clipped to the axis, then the points
  lo <- pmin(pmax(lower, xlim[1L]), xlim[2L])
  up <- pmin(pmax(upper, xlim[1L]), xlim[2L])
  ci <- ok & is.finite(lo) & is.finite(up)
  if (any(ci)) {
    graphics::segments(lo[ci], y[ci], up[ci], y[ci], col = col[ci],
                       lwd = lwd[ci], lty = lty[ci])
  }
  graphics::points(hr[ok], y[ok], pch = pch[ok], col = col[ok], cex = cex)
  graphics::axis(2, at = y, labels = labels, las = 1, tick = FALSE,
                 cex.axis = base_size)
  if (show_values) {
    graphics::axis(4, at = y, labels = values, las = 1, tick = FALSE,
                   cex.axis = base_size)
  }
  draw_legend(legend_pos, labels, inset = legend_inset, ncol = legend_ncol,
              cex = base_size, col = col, pch = pch, lty = lty, lwd = lwd)

  invisible(data.frame(method = method, label = labels, HR = unname(hr),
                       lower = unname(lower), upper = unname(upper),
                       text = unname(values), y = y,
                       stringsAsFactors = FALSE))
}
