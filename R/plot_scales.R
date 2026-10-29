#' Manual scales of the methods in the plot methods
#'
#' @param methods Character vector of methods, in the order of the legend.
#' @param labels Character vector of labels, one per method.
#' @param col \code{NULL} for the default palette of \code{plot_colors()}, or
#'   a vector of colours, recycled.
#' @param shape Point shapes, recycled.
#' @param linetype Line types, recycled.
#' @param aesthetics Aesthetics for which a scale is returned, among
#'   \code{"colour"}, \code{"shape"}, and \code{"linetype"}.
#'
#' @return A list of \code{ggplot2} manual scales, one per element of
#'   \code{aesthetics}, which map each method to its value and share the
#'   limits, breaks, labels, and (empty) title, so that their legends are
#'   combined even when a layer lacks some of the methods.
#'
#' @noRd
plot_scales <- function(methods, labels, col = NULL, shape = 19,
                        linetype = 1,
                        aesthetics = c("colour", "shape", "linetype")) {
  k <- length(methods)
  if (length(shape) == 0L) stop("'shape' must have at least one value")
  if (length(linetype) == 0L) {
    stop("'linetype' must have at least one value")
  }
  values <- list(colour = plot_colors(k, col),
                 shape = rep_len(shape, k),
                 linetype = rep_len(linetype, k))
  makers <- list(colour = ggplot2::scale_colour_manual,
                 shape = ggplot2::scale_shape_manual,
                 linetype = ggplot2::scale_linetype_manual)
  lapply(aesthetics, function(a) {
    makers[[a]](values = stats::setNames(values[[a]], methods),
                limits = methods, breaks = methods, labels = labels,
                name = NULL)
  })
}
