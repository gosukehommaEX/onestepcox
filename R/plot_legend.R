#' Position and layout of the legend of the plot methods
#'
#' @param legend_pos Position of the legend: \code{"bottom"}, \code{"top"},
#'   \code{"right"}, or \code{"left"} (outside the plot panel),
#'   \code{"none"} (no legend), or a numeric vector \code{c(x, y)} with
#'   elements between 0 and 1 for a position inside the panel (\code{c(0, 0)}
#'   is the bottom left corner and \code{c(1, 1)} the top right corner).
#' @param legend_ncol Number of columns of the legend, or \code{NULL}: at
#'   most 3 columns (and at most \code{k}) at the bottom or top, and 1 column
#'   otherwise.
#' @param k Number of entries of the legend.
#' @param aesthetics Aesthetics whose legends are combined into one.
#'
#' @return A list of two elements, a \code{ggplot2} theme that sets the
#'   position of the legend and the guides of \code{aesthetics} with
#'   \code{legend_ncol} columns, to be added to a plot after
#'   \code{plot_theme()}.
#'
#' @noRd
plot_legend <- function(legend_pos, legend_ncol = NULL, k = 1L,
                        aesthetics = c("colour", "shape", "linetype")) {
  sides <- c("bottom", "top", "right", "left", "none")
  if (is.character(legend_pos) && length(legend_pos) == 1L &&
      !is.na(legend_pos) && legend_pos %in% sides) {
    pos <- ggplot2::theme(legend.position = legend_pos)
    if (is.null(legend_ncol)) {
      legend_ncol <- if (legend_pos %in% c("bottom", "top")) min(k, 3L) else 1L
    }
  } else if (is.numeric(legend_pos) && length(legend_pos) == 2L &&
             all(is.finite(legend_pos)) && all(legend_pos >= 0) &&
             all(legend_pos <= 1)) {
    # The legend is placed with its own corner at the same relative position,
    # so that c(1, 1) puts it in the top right corner of the panel
    pos <- ggplot2::theme(
      legend.position = "inside",
      legend.position.inside = legend_pos,
      legend.justification = legend_pos,
      legend.background = ggplot2::element_rect(colour = "grey50")
    )
    if (is.null(legend_ncol)) legend_ncol <- 1L
  } else {
    stop("'legend_pos' must be one of ",
         paste0("\"", sides, "\"", collapse = ", "),
         ", or a numeric vector c(x, y) with elements between 0 and 1")
  }
  if (!is.numeric(legend_ncol) || length(legend_ncol) != 1L ||
      !is.finite(legend_ncol) || legend_ncol < 1) {
    stop("'legend_ncol' must be a positive integer")
  }
  g <- ggplot2::guide_legend(ncol = as.integer(legend_ncol))
  guides <- do.call(ggplot2::guides,
                    stats::setNames(rep(list(g), length(aesthetics)),
                                    aesthetics))
  list(pos, guides)
}
