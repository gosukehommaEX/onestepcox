#' Newton-Raphson iteration for the Cox estimate on risk-set summaries
#'
#' Solves the score equation of the two-group Cox model by Newton-Raphson
#' iteration from \code{eta0}, using the same terms as the closed-form
#' estimators. The score is decreasing in \code{eta}, so a step that does not
#' reduce the absolute score is halved (at most 30 times). Iteration stops
#' when the Newton step is smaller than \code{tol} in absolute value; the
#' information is then evaluated at the final value.
#'
#' @param eta0 Starting value of the log hazard ratio.
#' @param dt List returned by \code{deriv_table()}.
#' @param tol Positive convergence tolerance on the log hazard ratio scale.
#' @param maxit Maximum number of iterations.
#'
#' @return A list with elements \code{eta}, \code{info} (information at
#'   \code{eta}), \code{iter} (number of iterations), and \code{converged}.
#'
#' @noRd
newton_cox <- function(eta0, dt, tol, maxit = 100L) {
  eta <- eta0
  d <- score_derivs(eta, dt)
  for (it in seq_len(maxit)) {
    step <- d[["U"]] / d[["I"]]
    if (!is.finite(step)) break
    small <- abs(step) < tol
    d_new <- score_derivs(eta + step, dt)
    if (!small) {
      h <- 0L
      while (!isTRUE(abs(d_new[["U"]]) <= abs(d[["U"]])) && h < 30L) {
        step <- step / 2
        d_new <- score_derivs(eta + step, dt)
        h <- h + 1L
      }
    }
    eta <- eta + step
    d <- d_new
    if (small) {
      return(list(eta = eta, info = d[["I"]], iter = it, converged = TRUE))
    }
  }
  list(eta = eta, info = d[["I"]], iter = it, converged = FALSE)
}
