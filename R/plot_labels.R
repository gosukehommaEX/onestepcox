#' Labels of the methods in the plot methods
#'
#' @param methods Character vector of methods of \code{EstHR()}.
#' @param labels \code{NULL} for the default labels, or a character vector
#'   with one label per method.
#'
#' @return A character vector of labels. The default label of a method is its
#'   name with underscores replaced by hyphens, such as \code{"MH-Halley"}
#'   for \code{"MH_Halley"}.
#'
#' @noRd
plot_labels <- function(methods, labels = NULL) {
  if (is.null(labels)) return(gsub("_", "-", methods, fixed = TRUE))
  if (length(labels) != length(methods)) {
    stop("'labels' must have one label per method (", length(methods), ")")
  }
  as.character(labels)
}
