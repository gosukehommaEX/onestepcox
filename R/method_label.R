#' Descriptions of the methods of EstHR()
#'
#' @param method Character vector of methods of \code{EstHR()}.
#'
#' @return A character vector of one-line descriptions, named by the
#'   methods; \code{NA} for an unknown method.
#'
#' @noRd
method_label <- function(method) {
  labels <- c(
    MH_Halley      = "Halley step from the Mantel-Haenszel estimate (proposed)",
    MH_Newton      = "Newton step from the Mantel-Haenszel estimate",
    Pike_Halley    = "Halley step from the Pike estimate",
    Pike_Newton    = "Newton step from the Pike estimate",
    MH_Chebyshev   = "Chebyshev step from the Mantel-Haenszel estimate",
    Pike_Chebyshev = "Chebyshev step from the Pike estimate",
    MH             = "Mantel-Haenszel estimate",
    Pike           = "Pike estimate",
    Peto           = "Peto estimate",
    Cox            = "Cox estimate by survival::coxph()",
    Cox_NR         = "Cox estimate by Newton-Raphson iteration"
  )
  stats::setNames(unname(labels[method]), method)
}
