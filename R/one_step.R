#' One correction step on the Cox score equation
#'
#' From an anchor \code{eta0} with score \code{U}, information \code{I}, and
#' higher derivatives \code{J} and \code{L} evaluated there (see
#' \code{score_derivs()}), takes one step towards the root of the score:
#' \describe{
#'   \item{Newton}{\code{delta = U / I}.}
#'   \item{Halley}{\code{delta = 2 U I / (2 I^2 + U J)}. \code{NA} when the
#'     denominator is not positive.}
#'   \item{Chebyshev}{\code{delta = U / I - J U^2 / (2 I^3)}, the
#'     second-order truncation of the Halley step in \code{U / I}.}
#' }
#' The standard error is updated in the same pass from the second-order
#' Taylor expansion of the information about the anchor,
#' \code{(I + J delta + L delta^2 / 2)^(-1/2)}, and is \code{NA} when this
#' expansion is not positive.
#'
#' @param eta0 Anchor of the log hazard ratio.
#' @param d0 Named vector \code{c(U, I, J, L)} at \code{eta0}.
#' @param type One of \code{"Newton"}, \code{"Halley"}, \code{"Chebyshev"}.
#'
#' @return A named numeric vector \code{c(eta, se)}.
#'
#' @noRd
one_step <- function(eta0, d0, type) {
  U <- d0[["U"]]
  I <- d0[["I"]]
  J <- d0[["J"]]
  L <- d0[["L"]]
  if (type == "Newton") {
    delta <- U / I
  } else if (type == "Halley") {
    den <- 2 * I * I + U * J
    delta <- if (isTRUE(den > 0)) 2 * U * I / den else NA_real_
  } else {
    delta <- U / I - J * U * U / (2 * I * I * I)
  }
  info <- I + J * delta + L * delta * delta / 2
  c(eta = eta0 + delta,
    se = if (isTRUE(info > 0)) 1 / sqrt(info) else NA_real_)
}
