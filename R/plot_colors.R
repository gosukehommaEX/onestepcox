#' Colors of the lines and points of the plot methods
#'
#' @param n Number of colors.
#' @param col \code{NULL} for the default palette, or a vector of colors
#'   recycled to length \code{n}.
#'
#' @return A character vector of \code{n} colors. The default palette is
#'   that of Okabe and Ito, which is distinguishable under common forms of
#'   color blindness, starting with blue and vermilion, and recycled when
#'   \code{n} exceeds 8.
#'
#' @noRd
plot_colors <- function(n, col = NULL) {
  if (is.null(col)) {
    col <- c("#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00",
             "#56B4E9", "#000000", "#F0E442")
  }
  if (length(col) == 0L) stop("'col' must have at least one color")
  rep_len(col, n)
}
