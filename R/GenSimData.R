#' Generate simulated time-to-event data for two-group trials
#'
#' Generates Monte Carlo replicates of a two-group parallel-group trial with
#' exponential survival times, uniform accrual, optional independent
#' dropout, and analysis triggered when the E-th overall event is observed.
#' Unequal allocation between treatment and control groups is supported via
#' separate \code{n_T} and \code{n_C} arguments. Optionally, subjects are
#' assigned to strata with their own baseline hazards, and the recorded
#' times are rounded up to a time unit, which produces tied times.
#'
#' @param nsim Integer. Number of Monte Carlo replicates.
#' @param n_T Integer. Number of subjects on the treatment group
#'   (group label \code{"T"}).
#' @param n_C Integer. Number of subjects on the control group
#'   (group label \code{"C"}).
#' @param hazard_C Positive numeric. Hazard rate for the control group.
#' @param HR Positive numeric. True hazard ratio,
#'   \code{hazard_T = HR * hazard_C}.
#' @param a Positive numeric. Length of the accrual period; subject
#'   accrual times are drawn from \code{Uniform(0, a)}.
#' @param E Positive integer. Required total number of observed events
#'   (pooled across both groups) at which the analysis is triggered.
#' @param dropout_hazard Non-negative numeric. Hazard rate of independent
#'   dropout. Set to 0 (default) to disable dropout.
#' @param strata_prop Optional numeric vector of stratum probabilities
#'   (positive, summing to 1). Each subject is assigned to stratum \code{s}
#'   with probability \code{strata_prop[s]}, independently of the group.
#'   Default \code{NULL} (no strata).
#' @param strata_hazard Optional positive numeric vector of the same length
#'   as \code{strata_prop}: the hazard of a subject in stratum \code{s} is
#'   multiplied by \code{strata_hazard[s]} in both groups, so that the hazard
#'   ratio is common to all strata. Default \code{NULL}, which means 1 for
#'   every stratum.
#' @param time_unit Non-negative numeric. If positive, the recorded times
#'   are rounded up to multiples of \code{time_unit}, which produces tied
#'   times; whether a subject has an event and the analysis time are
#'   determined from the exact times. Default 0 (no rounding).
#' @param seed Optional integer seed for reproducibility.
#'
#' @return A \code{data.table} with one row per subject per simulation,
#'   carrying the columns \code{simID}, \code{group} (factor with levels
#'   \code{"T"} and \code{"C"}), \code{time}, and \code{event}, and, when
#'   \code{strata_prop} is given, \code{stratum} (integer stratum number).
#'   Subjects accrued after the E-th event are dropped, and subjects whose
#'   total calendar time exceeds the analysis time are administratively
#'   censored. The rows are sorted by \code{simID} and \code{time}, or by
#'   \code{simID}, \code{stratum}, and \code{time} when there are strata, so
#'   that the rows of each simulation can be passed to \code{\link{EstHR}}
#'   with \code{is_sorted = TRUE}.
#'
#' @details
#' The data-generating process is:
#' \itemize{
#'   \item Accrual times \code{u_{i,j} ~ Uniform(0, a)};
#'   \item Survival times \code{t_{i,j} ~ Exp(lambda_j)} with
#'     \code{lambda_T = HR * lambda_C};
#'   \item Optional independent dropout
#'     \code{d_{i,j} ~ Exp(lambda_D)};
#'   \item Analysis triggered when the E-th overall event is observed.
#' }
#' With strata, \code{lambda_C} is multiplied by \code{strata_hazard[s]}
#' for the subjects of stratum \code{s}, and \code{lambda_T} by the same
#' factor. The \code{time} and \code{event} columns (and \code{stratum}) are
#' the variables to be passed to \code{\link{EstHR}}. Without
#' \code{strata_prop} and with \code{time_unit = 0}, the random numbers are
#' drawn in the same order as in earlier versions, so that a given
#' \code{seed} gives the same data.
#'
#' @examples
#' d <- GenSimData(nsim = 100, n_T = 100, n_C = 100,
#'                 hazard_C = log(2) / 18, HR = 0.7,
#'                 a = 12, E = 110, seed = 1)
#' head(d)
#'
#' # Two strata with different baseline hazards, times in whole months
#' d2 <- GenSimData(nsim = 10, n_T = 100, n_C = 100,
#'                  hazard_C = log(2) / 18, HR = 0.7, a = 12, E = 110,
#'                  strata_prop = c(0.4, 0.6), strata_hazard = c(1, 2),
#'                  time_unit = 1, seed = 1)
#' head(d2)
#'
#' @import data.table
#' @importFrom stats runif rexp
#' @export
GenSimData <- function(nsim = 1e3, n_T, n_C, hazard_C, HR, a, E,
                       dropout_hazard = 0, strata_prop = NULL,
                       strata_hazard = NULL, time_unit = 0, seed = NULL) {

  # Input validation
  if (!is.null(seed)) set.seed(seed)
  stopifnot(
    nsim >= 1, n_T >= 1, n_C >= 1, hazard_C > 0, HR > 0, a > 0, E >= 1,
    dropout_hazard >= 0
  )
  N_per_sim <- as.integer(n_T) + as.integer(n_C)
  if (E > N_per_sim) {
    stop("'E' is larger than the total sample size per simulation ",
         "(n_T + n_C)")
  }
  if (!is.numeric(time_unit) || length(time_unit) != 1L ||
      !(time_unit >= 0)) {
    stop("'time_unit' must be a single non-negative number")
  }
  if (is.null(strata_prop)) {
    if (!is.null(strata_hazard)) {
      stop("'strata_hazard' requires 'strata_prop'")
    }
  } else {
    if (!is.numeric(strata_prop) || anyNA(strata_prop) ||
        any(strata_prop <= 0) || abs(sum(strata_prop) - 1) > 1e-8) {
      stop("'strata_prop' must be positive probabilities summing to 1")
    }
    if (is.null(strata_hazard)) {
      strata_hazard <- rep(1, length(strata_prop))
    }
    if (!is.numeric(strata_hazard) ||
        length(strata_hazard) != length(strata_prop) ||
        anyNA(strata_hazard) || any(strata_hazard <= 0)) {
      stop("'strata_hazard' must be positive numbers, one per stratum")
    }
  }

  hazard_T <- HR * hazard_C
  total_n  <- N_per_sim * nsim

  # Vectorized generation across all simulations and subjects.
  # Within each simulation, the first n_T subjects are treatment and the
  # next n_C subjects are control.
  group_one_sim <- rep.int(c("T", "C"), times = c(n_T, n_C))
  group_all     <- rep.int(group_one_sim, times = nsim)
  haz_all       <- ifelse(group_all == "T", hazard_T, hazard_C)

  # Optional strata: drawn before the other random numbers, so that the
  # draws without strata are unchanged
  if (!is.null(strata_prop)) {
    stratum_all <- sample.int(length(strata_prop), total_n, replace = TRUE,
                              prob = strata_prop)
    haz_all <- haz_all * strata_hazard[stratum_all]
  }

  dt <- data.table::data.table(
    simID        = rep.int(seq_len(nsim), rep.int(N_per_sim, nsim)),
    group        = group_all,
    accrual_time = stats::runif(total_n, min = 0, max = a),
    surv_time    = stats::rexp(total_n, rate = haz_all)
  )
  if (!is.null(strata_prop)) {
    dt[, stratum := stratum_all]
  }

  # Optional independent dropout
  if (dropout_hazard > 0) {
    dt[, dropout_time := stats::rexp(total_n, rate = dropout_hazard)]
    dt[, tte        := pmin(surv_time, dropout_time)]
    dt[, surv_event := as.integer(surv_time <= dropout_time)]
  } else {
    dt[, tte        := surv_time]
    dt[, surv_event := 1L]
  }

  dt[, total := accrual_time + tte]

  # Determine analysis_time per simID = calendar time of the E-th event
  data.table::setkey(dt, simID, total)
  cutoff <- dt[surv_event == 1L,
               .(analysis_time = total[E]),
               by = simID]
  if (anyNA(cutoff$analysis_time)) {
    stop("Some simulations have fewer than 'E' events; reduce E or ",
         "increase n_T / n_C / a")
  }
  dt <- cutoff[dt, on = "simID"]

  # Apply the cutoff: drop subjects accrued after analysis_time
  # (their on-trial time would be negative), and censor subjects whose
  # total calendar time exceeds analysis_time.
  dt <- dt[accrual_time <= analysis_time]
  dt[, event := as.integer(surv_event == 1L & total <= analysis_time)]
  dt[, time  := pmin(total, analysis_time) - accrual_time]

  # Optional rounding of the recorded times up to multiples of time_unit
  if (time_unit > 0) {
    dt[, time := ceiling(time / time_unit) * time_unit]
  }

  # Output: minimal columns
  if (is.null(strata_prop)) {
    out <- dt[, list(simID, group = factor(group, levels = c("T", "C")),
                     time, event)]
    data.table::setkey(out, simID, time)
  } else {
    out <- dt[, list(simID, stratum,
                     group = factor(group, levels = c("T", "C")),
                     time, event)]
    data.table::setkey(out, simID, stratum, time)
  }
  out[]
}
