#' Risk-set summaries at the distinct event times
#'
#' Tabulates, for every distinct event time (within strata when
#' \code{strata} is given), the numbers of treated and control subjects at
#' risk just before that time and the numbers of treated and total events at
#' that time. All estimators of \code{\link{EstHR}} except \code{"Cox"} are
#' computed from these summaries, and this function returns them from the
#' same internal code.
#'
#' Subjects censored at an event time are counted as at risk at that time, as
#' in \code{survival::coxph()}. Times are tied only when they are exactly
#' equal.
#'
#' @inheritParams EstHR
#'
#' @return A data frame with one row per distinct event time within stratum,
#'   ordered by stratum and time, and the columns
#'   \describe{
#'     \item{stratum}{Stratum label (only when \code{strata} is given).}
#'     \item{time}{Event time.}
#'     \item{n_T, n_C}{Numbers at risk in the treatment and control groups
#'       just before \code{time}.}
#'     \item{O_T}{Number of events in the treatment group at \code{time}.}
#'     \item{O}{Total number of events at \code{time}.}
#'   }
#'   The counts are stored as doubles.
#'
#' @seealso \code{\link{EstHR}}.
#'
#' @examples
#' time  <- c(1, 1, 2, 2, 3, 3, 4, 5, 5, 6)
#' event <- c(1, 1, 1, 0, 1, 1, 1, 0, 1, 1)
#' group <- c("C", "T", "C", "T", "C", "T", "T", "C", "C", "T")
#' RiskTable(time, event, group, control = "C")
#' RiskTable(time, event, group, control = "C",
#'           strata = rep(c("a", "b"), each = 5))
#'
#' @export
RiskTable <- function(time, event, group, control, strata = NULL,
                      is_sorted = FALSE) {
  trt <- check_inputs(time, event, group, control, strata)
  rs <- risk_sets(time, event, trt, strata, is_sorted)
  out <- data.frame(time = rs$time, n_T = rs$n_T, n_C = rs$n_C,
                    O_T = rs$O_T, O = rs$O)
  if (!is.null(strata)) {
    out <- cbind(data.frame(stratum = rs$stratum, stringsAsFactors = FALSE),
                 out)
  }
  out
}
