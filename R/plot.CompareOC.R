#' Plot operating characteristics of hazard ratio estimators
#'
#' Plots a result of \code{\link{CompareOC}}. With
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
#' @param lty,lwd,pch,cex Line type, line width, plotting symbol, and symbol
#'   size of the methods (recycled), for the lines and points; \code{pch} and
#'   \code{cex} are also used for the outlying points of the box plots.
#' @param base_size Character expansion of the axis labels, axis titles, and
#'   legend (the title is \code{1.2 * base_size}). Default 1.
#' @param legend_pos Position of the legend: a keyword of
#'   \code{graphics::legend()} such as \code{"topright"} (default) or
#'   \code{"bottomleft"}, numeric coordinates \code{c(x, y)} of its top left
#'   corner, or \code{"none"}.
#' @param legend_inset Inset of the legend from the margins as a fraction of
#'   the plot region. Default 0.02.
#' @param legend_ncol Number of columns of the legend. Default 1.
#' @param xlab,ylab,main Axis titles and main title. \code{ylab = NULL} gives
#'   a title that describes \code{what}.
#' @param ylim Limits of the vertical axis. Default \code{NULL} (the range
#'   of the values and of the reference line).
#' @param floor Smallest deviation shown with \code{"deviation"}; smaller
#'   deviations (including 0) are drawn at \code{floor}. Default
#'   \code{1e-16}.
#' @param ... Further arguments passed to \code{graphics::boxplot()}
#'   (\code{"deviation"}), such as \code{boxwex}, or graphical parameters
#'   passed to \code{graphics::plot()} (the others), such as \code{las}.
#'
#' @return The plotted values, invisibly: for \code{"deviation"}, a list of
#'   the deviations by scenario and method, named
#'   \code{"<scenario>:<method>"}; otherwise a matrix with one row per
#'   scenario and one column per method. The vertical axis is logarithmic for
#'   \code{"deviation"}, and for \code{"time"} when all times are positive.
#'   The graphical parameters are restored on exit.
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
#' plot(res, what = "bias", legend_pos = "bottomleft", base_size = 1.2)
#' plot(res, what = "coverage", methods = c("Cox", "MH_Halley"),
#'      col = c("black", "red"), lty = c(1, 2), ylim = c(0.8, 1))
#'
#' @export
plot.CompareOC <- function(x, what = c("deviation", "bias", "rmse",
                                       "coverage", "se_ratio", "time"),
                           methods = NULL, labels = NULL, col = NULL,
                           lty = 1, lwd = 2,
                           pch = 19, cex = 1, base_size = 1,
                           legend_pos = "topright", legend_inset = 0.02,
                           legend_ncol = 1L, xlab = "Scenario", ylab = NULL,
                           main = NULL, ylim = NULL, floor = 1e-16, ...) {
  what <- match.arg(what)
  s <- x$summary
  est <- x$estimates
  if (!is.data.frame(s) || !is.data.frame(est)) {
    stop("'x' must be the output of CompareOC()")
  }
  ref <- attr(x, "reference")
  if (is.null(ref)) ref <- s$method[1L]
  scen <- unique(as.character(s$label))
  all_methods <- unique(s$method)
  if (is.null(methods)) {
    methods <- if (what %in% c("deviation", "se_ratio")) {
      setdiff(all_methods, ref)
    } else {
      all_methods
    }
  }
  if (length(methods) == 0L) stop("there is no method to plot")
  if (!all(methods %in% all_methods)) {
    stop("'methods' must be methods of 'x': ",
         paste(all_methods, collapse = ", "))
  }
  k <- length(methods)
  n_lab <- length(scen)
  legend_labels <- plot_labels(methods, labels)
  col <- plot_colors(k, col)
  lty <- rep_len(lty, k)
  lwd <- rep_len(lwd, k)
  pch <- rep_len(pch, k)

  op <- graphics::par(cex.axis = base_size, cex.lab = base_size,
                      cex.main = 1.2 * base_size)
  on.exit(graphics::par(op))

  if (what == "deviation") {
    # Deviation of each replicate from the reference, matched by simID
    vals <- list()
    for (lab in scen) {
      e <- est[as.character(est$label) == lab, , drop = FALSE]
      r <- e[e$method == ref, , drop = FALSE]
      for (m in methods) {
        em <- e[e$method == m, , drop = FALSE]
        d <- abs(em$coef - r$coef[match(em$simID, r$simID)])
        d <- d[!is.na(d)]
        d[d < floor] <- floor
        vals[[paste(lab, m, sep = ":")]] <- d
      }
    }
    v <- unlist(vals)
    if (length(v) == 0L) stop("there is no deviation to plot")
    if (is.null(ylim)) {
      ylim <- range(v)
      if (ylim[1L] == ylim[2L]) ylim <- ylim * c(0.1, 10)
    }
    if (is.null(ylab)) ylab <- paste0("|coef - coef of ", ref, "|")
    at <- rep((seq_len(n_lab) - 1L) * (k + 1L), each = k) +
      rep(seq_len(k), n_lab)
    fill <- grDevices::adjustcolor(col, alpha.f = 0.35)
    graphics::boxplot(vals, at = at, log = "y", col = rep(fill, n_lab),
                      border = rep(col, n_lab), outpch = rep(pch, n_lab),
                      outcex = cex, xaxt = "n", xlim = c(0.5, max(at) + 0.5),
                      ylim = ylim, xlab = xlab, ylab = ylab, main = main,
                      names = rep("", length(at)), ...)
    graphics::axis(1, at = (seq_len(n_lab) - 1L) * (k + 1L) + (k + 1) / 2,
                   labels = scen, tick = FALSE)
    draw_legend(legend_pos, legend_labels, inset = legend_inset,
                ncol = legend_ncol, cex = base_size, fill = fill,
                border = col)
    return(invisible(vals))
  }

  column <- c(bias = "bias", rmse = "rmse", coverage = "coverage",
              se_ratio = "se_ratio_median", time = "time_per_sim")[[what]]
  mat <- matrix(NA_real_, nrow = n_lab, ncol = k,
                dimnames = list(scen, methods))
  for (j in seq_len(k)) {
    sm <- s[s$method == methods[j], , drop = FALSE]
    mat[, j] <- sm[[column]][match(scen, as.character(sm$label))]
  }
  if (!any(is.finite(mat))) stop("there is no value to plot")
  target <- switch(what, bias = 0, coverage = attr(x, "conf.level"),
                   se_ratio = 1, NULL)
  log_y <- what == "time" && all(mat[is.finite(mat)] > 0)
  if (is.null(ylim)) {
    ylim <- range(c(mat, target), finite = TRUE)
    if (ylim[1L] == ylim[2L]) {
      ylim <- if (log_y) ylim * c(0.5, 2) else ylim + c(-0.5, 0.5)
    }
  }
  if (is.null(ylab)) {
    ylab <- switch(what, bias = "Bias of the log hazard ratio",
                   rmse = "RMSE of the log hazard ratio",
                   coverage = "Coverage of the confidence interval",
                   se_ratio = paste0("Median ratio of se to that of ", ref),
                   time = "Seconds per replicate")
  }
  graphics::plot(NA, xlim = c(0.75, n_lab + 0.25), ylim = ylim,
                 log = if (log_y) "y" else "", xaxt = "n", xlab = xlab,
                 ylab = ylab, main = main, ...)
  graphics::axis(1, at = seq_len(n_lab), labels = scen)
  for (j in seq_len(k)) {
    graphics::lines(seq_len(n_lab), mat[, j], type = "b", col = col[j],
                    lty = lty[j], lwd = lwd[j], pch = pch[j], cex = cex)
  }
  if (!is.null(target)) {
    graphics::abline(h = target, col = "grey50", lty = 2)
  }
  draw_legend(legend_pos, legend_labels, inset = legend_inset,
              ncol = legend_ncol,
              cex = base_size, col = col, lty = lty, lwd = lwd, pch = pch)
  invisible(mat)
}
