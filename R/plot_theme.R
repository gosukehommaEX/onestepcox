#' Theme of the plot methods
#'
#' @param base_size Base font size in points, passed to
#'   \code{ggplot2::theme_bw()}.
#'
#' @return A complete \code{ggplot2} theme: \code{theme_bw()} without the
#'   minor grid lines and the legend title, with legend keys wide enough to
#'   show the line types. It is added to a plot before \code{plot_legend()},
#'   whose position of the legend would otherwise be reset.
#'
#' @noRd
plot_theme <- function(base_size = 11) {
  if (!is.numeric(base_size) || length(base_size) != 1L ||
      !is.finite(base_size) || base_size <= 0) {
    stop("'base_size' must be a single positive number")
  }
  ggplot2::theme_bw(base_size = base_size) +
    ggplot2::theme(panel.grid.minor = ggplot2::element_blank(),
                   legend.title = ggplot2::element_blank(),
                   legend.key.width = ggplot2::unit(2, "lines"))
}
