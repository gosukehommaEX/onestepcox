#' Score and its derivatives for the two-group Cox model
#'
#' Evaluates, at the log hazard ratio \code{eta}, the score
#' \code{U = O_T - sum_m w_m p_m}, the observed information
#' \code{I = -U' = sum_m w_m p_m q_m}, and the higher derivatives
#' \code{J = -U'' = sum_m w_m p_m q_m (q_m - p_m)} and
#' \code{L = -U''' = sum_m w_m p_m q_m (1 - 6 p_m q_m)}, where
#' \code{p_m = n_T_m exp(eta) / (n_C_m + n_T_m exp(eta))} and
#' \code{q_m = 1 - p_m}. All four are obtained in one pass over the terms.
#'
#' @param eta Log hazard ratio (scalar).
#' @param dt List returned by \code{deriv_table()}.
#'
#' @return A named numeric vector \code{c(U, I, J, L)}.
#'
#' @noRd
score_derivs <- function(eta, dt) {
  a <- dt$n_T * exp(eta)
  den <- dt$n_C + a
  p <- a / den
  q <- dt$n_C / den
  wpq <- dt$w * p * q
  c(U = dt$O_T - sum(dt$w * p),
    I = sum(wpq),
    J = sum(wpq * (q - p)),
    L = sum(wpq * (1 - 6 * p * q)))
}
