#' Closed-form anchors: Mantel-Haenszel, Pike, and Peto estimators
#'
#' Computes one of the closed-form estimators of the log hazard ratio and its
#' standard error from the risk-set summaries. With
#' \code{n_k = n_T_k + n_C_k} and \code{O_C_k = O_k - O_T_k},
#' \describe{
#'   \item{MH}{\code{log(R / S)} with \code{R = sum O_T_k n_C_k / n_k} and
#'     \code{S = sum O_C_k n_T_k / n_k}. Variance
#'     \code{sum(O_k n_T_k n_C_k / n_k^2) / (R S)}.}
#'   \item{Pike}{\code{log((O_T E_C) / (O_C E_T))} with
#'     \code{E_T = sum n_T_k O_k / n_k} and \code{E_C = sum n_C_k O_k / n_k}.
#'     Variance \code{1 / E_T + 1 / E_C}.}
#'   \item{Peto}{\code{(O_T - E_T) / V} with the hypergeometric variance
#'     \code{V = sum n_T_k n_C_k O_k (n_k - O_k) / (n_k^2 (n_k - 1))} over
#'     the times with \code{n_k > 1}. Variance \code{1 / V}.}
#' }
#'
#' @param rs List returned by \code{risk_sets()}.
#' @param anchor One of \code{"MH"}, \code{"Pike"}, \code{"Peto"}.
#'
#' @return A named numeric vector \code{c(eta, se)}. The values are not
#'   finite when the estimator is undefined (for example, when no event
#'   occurs in one of the groups).
#'
#' @noRd
anchor_estimate <- function(rs, anchor) {
  n_k <- rs$n_T + rs$n_C
  if (anchor == "MH") {
    R <- sum(rs$O_T * rs$n_C / n_k)
    S <- sum((rs$O - rs$O_T) * rs$n_T / n_k)
    v <- sum(rs$O * rs$n_T * rs$n_C / (n_k * n_k)) / (R * S)
    return(c(eta = log(R / S), se = sqrt(v)))
  }
  O_T <- sum(rs$O_T)
  E_T <- sum(rs$n_T * rs$O / n_k)
  if (anchor == "Pike") {
    O_C <- sum(rs$O) - O_T
    E_C <- sum(rs$n_C * rs$O / n_k)
    return(c(eta = log((O_T * E_C) / (O_C * E_T)),
             se = sqrt(1 / E_T + 1 / E_C)))
  }
  ok <- n_k > 1
  V <- sum(rs$n_T[ok] * rs$n_C[ok] * rs$O[ok] * (n_k[ok] - rs$O[ok]) /
             (n_k[ok] * n_k[ok] * (n_k[ok] - 1)))
  c(eta = (O_T - E_T) / V, se = 1 / sqrt(V))
}
