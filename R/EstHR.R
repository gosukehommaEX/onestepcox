#' Estimate the hazard ratio of a two-group comparison
#'
#' Estimates the hazard ratio (HR) of the treatment group relative to the
#' control group, with its standard error on the log scale and a Wald
#' confidence interval, by one of eleven methods: the proposed closed-form
#' estimator (\code{"MH_Halley"}), five other one-step corrections, the
#' Mantel-Haenszel, Pike, and Peto estimators, and the Cox partial likelihood
#' estimator computed either by \code{survival::coxph()} or by Newton-Raphson
#' iteration on the same risk-set summaries as the closed-form estimators.
#' Tied event times are handled by the Efron or the Breslow method, and the
#' baseline hazard may differ between strata.
#'
#' @param time Numeric vector of follow-up times.
#' @param event Vector of event indicators, 1 (event) or 0 (censored), or
#'   logical, aligned with \code{time}.
#' @param group Vector of group labels aligned with \code{time}. Subjects with
#'   \code{group != control} form the treatment group.
#' @param control Single value of \code{group} that identifies the control
#'   group.
#' @param method Character string, one of \code{"MH_Halley"} (proposed),
#'   \code{"MH_Newton"}, \code{"Pike_Halley"}, \code{"Pike_Newton"},
#'   \code{"MH_Chebyshev"}, \code{"Pike_Chebyshev"}, \code{"MH"},
#'   \code{"Pike"}, \code{"Peto"}, \code{"Cox"}, and \code{"Cox_NR"}. See
#'   Details.
#' @param ties Character string, \code{"efron"} (default, as in
#'   \code{survival::coxph()}) or \code{"breslow"}: the handling of tied
#'   event times in the partial likelihood. It affects the one-step methods,
#'   \code{"Cox"}, and \code{"Cox_NR"}, and not \code{"MH"}, \code{"Pike"},
#'   or \code{"Peto"}.
#' @param strata Optional vector of stratum labels aligned with \code{time}.
#'   Each stratum has its own baseline hazard and the HR is common to all
#'   strata. Default \code{NULL} (no stratification).
#' @param conf.level Confidence level of the interval. Default 0.95.
#' @param is_sorted Logical. If \code{TRUE}, the rows are taken to be sorted
#'   in ascending order of \code{time} within strata, with the rows of each
#'   stratum contiguous, and the internal sorting is skipped. An error is
#'   raised if they are not. Default \code{FALSE}.
#' @param tol Positive convergence tolerance of \code{"Cox"} (passed to
#'   \code{survival::coxph.control()} as \code{eps}) and of \code{"Cox_NR"}
#'   (on the Newton step of the log HR). Default \code{1e-12}.
#'
#' @return A named numeric vector with elements \code{coef} (log HR),
#'   \code{exp(coef)} (HR), \code{se(coef)} (standard error of the log HR),
#'   and the lower and upper confidence limits of the HR, named as in
#'   \code{summary()} of a \code{coxph} fit (\code{"lower .95"} and
#'   \code{"upper .95"} for \code{conf.level = 0.95}). All elements are
#'   \code{NA} when the estimate cannot be computed, for example when there is
#'   no event or no event in one of the groups.
#'
#' @details
#' Let \code{t_k} (\code{k = 1, ..., K}) be the distinct event times, within
#' strata when \code{strata} is given. At \code{t_k}, let \code{n_T_k} and
#' \code{n_C_k} be the numbers at risk in the treatment and control groups
#' just before \code{t_k}, \code{O_T_k} and \code{O_C_k} the numbers of
#' events, \code{n_k = n_T_k + n_C_k}, and \code{O_k = O_T_k + O_C_k}. These
#' summaries are returned by \code{\link{RiskTable}} and are shared by all
#' methods except \code{"Cox"}. Let \code{O_T = sum O_T_k},
#' \code{O_C = sum O_C_k}, \code{E_T = sum n_T_k O_k / n_k}, and
#' \code{E_C = sum n_C_k O_k / n_k}.
#'
#' The closed-form estimators of the log HR and their standard errors are
#' \describe{
#'   \item{MH}{Mantel-Haenszel estimator \code{log(R / S)}, where
#'     \code{R = sum O_T_k n_C_k / n_k} and \code{S = sum O_C_k n_T_k / n_k}.
#'     The variance is \code{sum(O_k n_T_k n_C_k / n_k^2) / (R S)}.}
#'   \item{Pike}{\code{log((O_T E_C) / (O_C E_T))}, with variance
#'     \code{1 / E_T + 1 / E_C}.}
#'   \item{Peto}{\code{(O_T - E_T) / V}, with variance \code{1 / V}, where
#'     \code{V = sum n_T_k n_C_k O_k (n_k - O_k) / (n_k^2 (n_k - 1))} is the
#'     hypergeometric variance of the log-rank statistic.}
#' }
#'
#' The one-step methods \code{"<anchor>_<step>"} start from the
#' Mantel-Haenszel (\code{"MH"}) or the Pike (\code{"Pike"}) estimate
#' \code{eta_0} and take one step on the score equation \code{U(eta) = 0} of
#' the Cox partial likelihood. With \code{U_0 = U(eta_0)},
#' \code{I_0 = -U'(eta_0)}, \code{J_0 = -U''(eta_0)}, and
#' \code{L_0 = -U'''(eta_0)}, all computed in one pass over the risk sets,
#' the step \code{delta} is
#' \describe{
#'   \item{Newton}{\code{U_0 / I_0};}
#'   \item{Halley}{\code{2 U_0 I_0 / (2 I_0^2 + U_0 J_0)} (the estimate is
#'     \code{NA} if the denominator is not positive);}
#'   \item{Chebyshev}{\code{U_0 / I_0 - J_0 U_0^2 / (2 I_0^3)}, the
#'     second-order truncation of the Halley step.}
#' }
#' The estimate of the log HR is \code{eta_0 + delta}, and its standard error
#' is \code{(I_0 + J_0 delta + L_0 delta^2 / 2)^(-1/2)}, the information
#' expanded about \code{eta_0} to the second order and evaluated at the
#' estimate. The proposed estimator is \code{"MH_Halley"}.
#'
#' Under the Breslow method, the score is
#' \code{U(eta) = sum (O_T_k - O_k p_k)} with
#' \code{p_k = n_T_k exp(eta) / (n_C_k + n_T_k exp(eta))}. Under the Efron
#' method, each tie block contributes \code{O_k} terms,
#' \code{r = 0, ..., O_k - 1}, in which \code{n_T_k} and \code{n_C_k} are
#' reduced by \code{(r / O_k) O_T_k} and \code{(r / O_k) O_C_k}.
#'
#' The Cox estimate is computed in two ways.
#' \describe{
#'   \item{Cox}{\code{survival::coxph()} with the convergence tolerance
#'     \code{eps = tol}, at most 100 iterations, and \code{timefix = FALSE},
#'     so that times are tied only when they are exactly equal, as in
#'     \code{\link{RiskTable}}. The standard error is that of
#'     \code{coxph()}.}
#'   \item{Cox_NR}{Newton-Raphson iteration on the risk-set summaries,
#'     started from the Mantel-Haenszel estimate. A step that does not reduce
#'     the absolute score is halved. Iteration stops when the Newton step is
#'     smaller than \code{tol}, and the standard error is
#'     \code{I(eta)^(-1/2)} at the final value. If the iteration does not
#'     converge in 100 steps, a warning is given and \code{NA} is returned.}
#' }
#'
#' Swapping the roles of the groups changes the sign of the log HR of every
#' method and leaves the standard error unchanged.
#'
#' @seealso \code{\link{RiskTable}} for the risk-set summaries.
#'
#' @examples
#' set.seed(1)
#' n <- 200
#' time  <- c(rexp(n / 2, rate = 0.10), rexp(n / 2, rate = 0.15))
#' event <- rbinom(n, 1, 0.8)
#' group <- rep(c("C", "T"), each = n / 2)
#' EstHR(time, event, group, control = "C", method = "MH_Halley")
#' EstHR(time, event, group, control = "C", method = "Cox")
#'
#' # All methods, Breslow handling of ties
#' methods <- c("MH_Halley", "MH_Newton", "Pike_Halley", "Pike_Newton",
#'              "MH_Chebyshev", "Pike_Chebyshev", "MH", "Pike", "Peto",
#'              "Cox", "Cox_NR")
#' t(sapply(methods, function(m) {
#'   EstHR(ceiling(time), event, group, "C", m, ties = "breslow")
#' }))
#'
#' @importFrom survival Surv strata
#' @export
EstHR <- function(time, event, group, control, method,
                  ties = "efron", strata = NULL, conf.level = 0.95,
                  is_sorted = FALSE, tol = 1e-12) {

  valid_methods <- c("MH_Halley", "MH_Newton", "Pike_Halley", "Pike_Newton",
                     "MH_Chebyshev", "Pike_Chebyshev", "MH", "Pike", "Peto",
                     "Cox", "Cox_NR")
  if (length(method) != 1L || !method %in% valid_methods) {
    stop(sprintf("'method' must be one of: %s",
                 paste(valid_methods, collapse = ", ")))
  }
  if (length(ties) != 1L || !ties %in% c("efron", "breslow")) {
    stop("'ties' must be \"efron\" or \"breslow\"")
  }
  if (!is.numeric(conf.level) || length(conf.level) != 1L ||
      !(conf.level > 0 && conf.level < 1)) {
    stop("'conf.level' must be a single number between 0 and 1")
  }
  if (!is.numeric(tol) || length(tol) != 1L || !(tol > 0)) {
    stop("'tol' must be a single positive number")
  }
  trt <- check_inputs(time, event, group, control, strata)

  # Output in the form of summary.coxph(): coef, exp(coef), se(coef), limits
  z <- stats::qnorm((1 + conf.level) / 2)
  out_names <- c("coef", "exp(coef)", "se(coef)",
                 paste0("lower .", round(100 * conf.level, 2)),
                 paste0("upper .", round(100 * conf.level, 2)))
  na_out <- stats::setNames(rep(NA_real_, 5L), out_names)
  make_out <- function(b, se) {
    if (!is.finite(b)) return(na_out)
    if (!is.finite(se)) se <- NA_real_
    stats::setNames(c(b, exp(b), se, exp(b - z * se), exp(b + z * se)),
                    out_names)
  }

  if (length(time) == 0L || sum(event) == 0) return(na_out)

  # Cox by survival::coxph()
  if (method == "Cox") {
    x <- as.numeric(trt)
    ctrl <- survival::coxph.control(
      eps = tol,
      toler.chol = min(.Machine$double.eps^0.875, tol / 10),
      iter.max = 100L,
      timefix = FALSE
    )
    if (is.null(strata)) {
      fit <- survival::coxph(Surv(time, event) ~ x, ties = ties,
                             control = ctrl)
    } else {
      s <- strata
      fit <- survival::coxph(Surv(time, event) ~ x + strata(s), ties = ties,
                             control = ctrl)
    }
    return(make_out(unname(fit$coefficients[1L]),
                    sqrt(unname(fit$var[1L, 1L]))))
  }

  # All other methods share the risk-set summaries
  rs <- risk_sets(time, event, trt, strata, is_sorted)
  O_T <- sum(rs$O_T)
  if (O_T == 0 || O_T == sum(rs$O)) return(na_out)

  if (method %in% c("MH", "Pike", "Peto")) {
    a <- anchor_estimate(rs, method)
    return(make_out(a[["eta"]], a[["se"]]))
  }

  anchor <- if (startsWith(method, "Pike")) "Pike" else "MH"
  a <- anchor_estimate(rs, anchor)
  if (!is.finite(a[["eta"]])) return(na_out)
  dt <- deriv_table(rs, ties)

  if (method == "Cox_NR") {
    fit <- newton_cox(a[["eta"]], dt, tol)
    if (!fit$converged) {
      warning("Newton-Raphson iteration of 'Cox_NR' did not converge; ",
              "NA is returned")
      return(na_out)
    }
    return(make_out(fit$eta, 1 / sqrt(fit$info)))
  }

  d0 <- score_derivs(a[["eta"]], dt)
  if (!isTRUE(d0[["I"]] > 0)) return(na_out)
  st <- one_step(a[["eta"]], d0, sub("^.*_", "", method))
  make_out(st[["eta"]], st[["se"]])
}
