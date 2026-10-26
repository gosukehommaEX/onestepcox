#' Weighted risk sets for the derivatives of the log partial likelihood
#'
#' Converts the risk-set summaries of \code{risk_sets()} into the form used by
#' \code{score_derivs()}, in which the log partial likelihood of the two-group
#' Cox model is
#' \code{O_T eta - sum_m w_m log(n_C_m + n_T_m exp(eta))}.
#' Under the Breslow handling of ties, there is one term per distinct event
#' time with \code{w_m = O_k}. Under the Efron handling, the tie block at event
#' time k with \code{O_k} events is split into \code{O_k} terms with
#' \code{w_m = 1}; in the r-th term (\code{r = 0, ..., O_k - 1}), a fraction
#' \code{r / O_k} of the treated and control events at that time is removed
#' from the numbers at risk. Without ties, both give the same terms.
#'
#' @param rs List returned by \code{risk_sets()}.
#' @param ties Either \code{"efron"} or \code{"breslow"}.
#'
#' @return A list with elements \code{n_T}, \code{n_C}, \code{w} (vectors of
#'   equal length), and \code{O_T} (total number of treated events).
#'
#' @noRd
deriv_table <- function(rs, ties) {
  O_T_tot <- sum(rs$O_T)
  if (ties == "breslow" || all(rs$O == 1)) {
    return(list(n_T = rs$n_T, n_C = rs$n_C, w = rs$O, O_T = O_T_tot))
  }
  idx <- rep.int(seq_along(rs$O), rs$O)
  frac <- (sequence(rs$O) - 1) / rs$O[idx]
  list(n_T = rs$n_T[idx] - frac * rs$O_T[idx],
       n_C = rs$n_C[idx] - frac * (rs$O[idx] - rs$O_T[idx]),
       w = rep.int(1, length(idx)),
       O_T = O_T_tot)
}
