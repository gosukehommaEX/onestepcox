#' Draw the legend of a plot method
#'
#' @param pos Position of the legend: \code{"none"} (no legend), one of the
#'   keywords of \code{graphics::legend()} (\code{"topright"}, \code{"top"},
#'   \code{"topleft"}, \code{"left"}, \code{"center"}, \code{"right"},
#'   \code{"bottomright"}, \code{"bottom"}, \code{"bottomleft"}), or a
#'   numeric vector \code{c(x, y)} of the coordinates of the top left corner
#'   in the units of the axes.
#' @param legend Character vector of labels.
#' @param inset Inset from the margins as a fraction of the plot region, used
#'   with a keyword.
#' @param ncol Number of columns of the legend.
#' @param cex Character expansion of the legend text.
#' @param ... Further arguments of \code{graphics::legend()}, such as
#'   \code{col}, \code{lty}, \code{lwd}, \code{pch}, \code{fill}, and
#'   \code{border}.
#'
#' @return The value of \code{graphics::legend()}, invisibly, or \code{NULL}
#'   for \code{pos = "none"}.
#'
#' @noRd
draw_legend <- function(pos, legend, inset = 0.02, ncol = 1L, cex = 1, ...) {
  keywords <- c("bottomright", "bottom", "bottomleft", "left", "topleft",
                "top", "topright", "right", "center")
  if (identical(pos, "none")) return(invisible(NULL))
  if (is.character(pos) && length(pos) == 1L && pos %in% keywords) {
    out <- graphics::legend(pos, legend = legend, inset = inset, ncol = ncol,
                            cex = cex, bg = "white", ...)
  } else if (is.numeric(pos) && length(pos) == 2L && all(is.finite(pos))) {
    out <- graphics::legend(pos[1L], pos[2L], legend = legend, ncol = ncol,
                            cex = cex, bg = "white", ...)
  } else {
    stop("'legend_pos' must be \"none\", one of ",
         paste0("\"", keywords, "\"", collapse = ", "),
         ", or numeric coordinates c(x, y)")
  }
  invisible(out)
}
