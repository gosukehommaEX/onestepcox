#' Analysis time of the large-sample limit of GenSimData()
#'
#' Finds the calendar time \code{tau} at which the expected fraction of the
#' planned subjects with an observed event reaches \code{e_frac}, for
#' uniform accrual on \code{(0, a)}, exponential event times with hazard
#' \code{lambda_j} in group j (allocation fraction \code{pi_j}), and
#' independent exponential dropout with hazard \code{lambda_D}. With
#' \code{r_j = lambda_j + lambda_D} and \code{w = min(tau, a)}, the fraction is
#' \code{sum_j pi_j (lambda_j / r_j)
#'   (w - (exp(-r_j (tau - w)) - exp(-r_j tau)) / r_j) / a},
#' which increases in \code{tau}.
#'
#' @param hazard_T,hazard_C Event hazards of the two groups.
#' @param pi_T Allocation fraction of the treatment group.
#' @param a Length of the accrual period.
#' @param e_frac Target event fraction, \code{E / (n_T + n_C)}.
#' @param dropout_hazard Dropout hazard.
#'
#' @return The analysis time \code{tau}. An error is raised when
#'   \code{e_frac} is not below the limit of the fraction as \code{tau} grows,
#'   \code{sum_j pi_j lambda_j / r_j}.
#'
#' @noRd
limit_analysis_time <- function(hazard_T, hazard_C, pi_T, a, e_frac,
                                dropout_hazard) {
  lam <- c(hazard_T, hazard_C)
  pi_j <- c(pi_T, 1 - pi_T)
  r <- lam + dropout_hazard
  frac_max <- sum(pi_j * lam / r)
  if (!(e_frac < frac_max)) {
    stop("the event fraction E / (n_T + n_C) cannot be reached; ",
         "it must be below ", signif(frac_max, 6))
  }
  frac <- function(tau) {
    w <- min(tau, a)
    sum(pi_j * lam / r * (w - (exp(-r * (tau - w)) - exp(-r * tau)) / r)) / a
  }
  upper <- a
  while (frac(upper) <= e_frac) upper <- 2 * upper
  stats::uniroot(function(tau) frac(tau) - e_frac, c(0, upper),
                 tol = 1e-10, maxiter = 1000L)$root
}
