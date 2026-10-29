#' Estimate the hazard ratio of a two-group comparison
#'
#' Estimates the hazard ratio (HR) of the treatment group relative to the
#' control group, with its standard error on the log scale and a Wald
#' confidence interval, by one or more of eleven methods: the proposed
#' closed-form estimator (\code{"MH_Halley"}), five other one-step
#' corrections, the
#' Mantel-Haenszel, Pike, and Peto estimators, and the Cox partial likelihood
#' estimator computed either by \code{survival::coxph()} or by Newton-Raphson
#' iteration on the same risk-set summaries as the closed-form estimators.
#' Tied event times are handled by the Efron or the Breslow method, and the
#' baseline hazard may differ between strata. The data are given either as
#' vectors or by a formula \code{Surv(time, event) ~ group + strata(s)} with
#' a data frame, as in \code{survival::coxph()}.
#'
#' @param time Numeric vector of follow-up times, or a formula
#'   \code{Surv(time, event) ~ group}, optionally with \code{+ strata(s)}
#'   (see Formula interface in Details).
#' @param event Vector of event indicators, 1 (event) or 0 (censored), or
#'   logical, aligned with \code{time}. With a formula, the data frame may be
#'   given here instead of in \code{data}, as in \code{coxph(formula, data)}.
#' @param group Vector of group labels aligned with \code{time}. Subjects with
#'   \code{group != control} form the treatment group. Not used with a
#'   formula.
#' @param control Single value of \code{group} that identifies the control
#'   group. With a formula, the default is the first level of the group
#'   variable if it is a factor and its smallest value otherwise, which is
#'   the reference level of \code{coxph()}.
#' @param method Character vector of one or more of \code{"MH_Halley"}
#'   (proposed), \code{"MH_Newton"}, \code{"Pike_Halley"},
#'   \code{"Pike_Newton"}, \code{"MH_Chebyshev"}, \code{"Pike_Chebyshev"},
#'   \code{"MH"}, \code{"Pike"}, \code{"Peto"}, \code{"Cox"}, and
#'   \code{"Cox_NR"}, without duplicates. See Details.
#' @param ties Character string, \code{"efron"} (default, as in
#'   \code{survival::coxph()}) or \code{"breslow"}: the handling of tied
#'   event times in the partial likelihood. It affects the one-step methods,
#'   \code{"Cox"}, and \code{"Cox_NR"}, and not \code{"MH"}, \code{"Pike"},
#'   or \code{"Peto"}.
#' @param strata Optional vector of stratum labels aligned with \code{time}.
#'   Each stratum has its own baseline hazard and the HR is common to all
#'   strata. Default \code{NULL} (no stratification). With a formula, strata
#'   are given by \code{strata()} in the formula instead.
#' @param conf.level Confidence level of the interval. Default 0.95.
#' @param is_sorted Logical. If \code{TRUE}, the rows are taken to be sorted
#'   in ascending order of \code{time} within strata, with the rows of each
#'   stratum contiguous, and the internal sorting is skipped. An error is
#'   raised if they are not. Default \code{FALSE}.
#' @param tol Positive convergence tolerance of \code{"Cox"} (passed to
#'   \code{survival::coxph.control()} as \code{eps}) and of \code{"Cox_NR"}
#'   (on the Newton step of the log HR). Default \code{1e-12}.
#' @param data A data frame (or list) in which the variables of a formula are
#'   evaluated. Used only with a formula. Default \code{NULL} (the
#'   environment of the formula).
#'
#' @return An object of class \code{"EstHR"}. For one method, a named numeric
#'   vector with elements \code{coef} (log HR), \code{exp(coef)} (HR),
#'   \code{se(coef)} (standard error of the log HR), and the lower and upper
#'   confidence limits of the HR, named as in \code{summary()} of a
#'   \code{coxph} fit (\code{"lower .95"} and \code{"upper .95"} for
#'   \code{conf.level = 0.95}). For several methods, a matrix with one row
#'   per method, in the order of \code{method}, and these five columns. The
#'   elements are \code{NA} when the estimate cannot be computed, for example
#'   when there is no event or no event in one of the groups.
#'
#'   The object has the attributes \code{method}, \code{ties}, and
#'   \code{conf.level}; \code{groups} (the labels of the treatment and
#'   control groups; the treatment label is the group of the first subject
#'   of the treatment group); \code{n} and \code{events} (numbers of subjects
#'   and events by group, named \code{treatment} and \code{control});
#'   \code{n_strata} (1 without strata); and \code{n_removed} (rows removed
#'   by the formula interface for missing values). Subsetting with \code{[}
#'   or \code{[[} returns plain numbers. See \code{\link{print.EstHR}} and
#'   \code{\link{plot.EstHR}} for the methods of \code{print()} and
#'   \code{plot()}.
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
#' Formula interface. With a formula \code{Surv(time, event) ~ group} as the
#' first argument, the variables are taken from \code{data} (or from the
#' environment of the formula), as in \code{survival::coxph()}. The response
#' must be \code{Surv(time, event)} for right-censored data; an event
#' indicator coded 1 and 2 (2 for an event) is accepted, as by
#' \code{Surv()}. The right-hand side has one group term, a variable or an
#' expression such as \code{factor(x)} with exactly two distinct values, and
#' optionally \code{strata()} terms; several variables in \code{strata()} (or
#' several \code{strata()} terms) define strata by their combinations. Other
#' terms, interactions, and \code{.} are not allowed. Rows with a missing
#' value in any of the variables are removed, and their number is stored in
#' the attribute \code{n_removed}. The arguments \code{group} and
#' \code{strata} are not used with a formula, and \code{method} and the
#' other arguments are given by name.
#'
#' @seealso \code{\link{RiskTable}} for the risk-set summaries,
#'   \code{\link{print.EstHR}} and \code{\link{plot.EstHR}}.
#'
#' @examples
#' # Formula interface, as in survival::coxph()
#' EstHR(Surv(rfstime, status) ~ hormon, data = survival::gbsg,
#'       method = "MH_Halley")
#'
#' # Several methods, stratified by tumour grade
#' est <- EstHR(Surv(rfstime, status) ~ hormon + strata(grade),
#'              data = survival::gbsg,
#'              method = c("MH_Halley", "Pike", "Peto", "Cox"))
#' est
#' plot(est)
#'
#' # Vectors
#' set.seed(1)
#' n <- 200
#' time  <- c(rexp(n / 2, rate = 0.10), rexp(n / 2, rate = 0.15))
#' event <- rbinom(n, 1, 0.8)
#' group <- rep(c("C", "T"), each = n / 2)
#' EstHR(time, event, group, control = "C", method = "MH_Halley")
#'
#' # All methods, Breslow handling of ties; subsetting gives a plain matrix
#' methods <- c("MH_Halley", "MH_Newton", "Pike_Halley", "Pike_Newton",
#'              "MH_Chebyshev", "Pike_Chebyshev", "MH", "Pike", "Peto",
#'              "Cox", "Cox_NR")
#' est_all <- EstHR(ceiling(time), event, group, "C", methods,
#'                  ties = "breslow")
#' est_all[, c("coef", "se(coef)")]
#'
#' @importFrom survival Surv strata
#' @export
EstHR <- function(time, event, group, control, method,
                  ties = "efron", strata = NULL, conf.level = 0.95,
                  is_sorted = FALSE, tol = 1e-12, data = NULL) {

  # Formula interface: Surv(time, event) ~ group + strata(s)
  n_removed <- 0L
  if (inherits(time, "formula")) {
    f <- parse_formula(time, data,
                       event = if (missing(event)) NULL else event,
                       group = if (missing(group)) NULL else group,
                       strata = strata)
    time <- f$time
    event <- f$event
    group <- f$group
    strata <- f$strata
    if (missing(control)) control <- f$control
    n_removed <- f$n_removed
  } else if (!is.null(data)) {
    stop("'data' is used only with a formula")
  }

  valid_methods <- c("MH_Halley", "MH_Newton", "Pike_Halley", "Pike_Newton",
                     "MH_Chebyshev", "Pike_Chebyshev", "MH", "Pike", "Peto",
                     "Cox", "Cox_NR")
  if (!is.character(method) || length(method) == 0L || anyNA(method) ||
      !all(method %in% valid_methods)) {
    stop(sprintf("'method' must be one or more of: %s",
                 paste(valid_methods, collapse = ", ")))
  }
  if (anyDuplicated(method) > 0L) {
    stop("'method' must not contain duplicates")
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
  est <- estimate_hr(time, event, trt, strata, method, ties, is_sorted, tol)
  b <- est[, "coef"]
  b[!is.finite(b)] <- NA_real_
  se <- est[, "se"]
  se[is.na(b) | !is.finite(se)] <- NA_real_
  out <- cbind(b, exp(b), se, exp(b - z * se), exp(b + z * se))
  dimnames(out) <- list(method, out_names)

  if (length(method) == 1L) {
    out <- out[1L, ]
    class(out) <- c("EstHR", "numeric")
  } else {
    class(out) <- c("EstHR", "matrix", "array")
  }
  n_T <- sum(trt)
  e_T <- sum(event[trt])
  attr(out, "method") <- method
  attr(out, "ties") <- ties
  attr(out, "conf.level") <- conf.level
  attr(out, "groups") <- c(
    treatment = if (n_T > 0) as.character(group[match(TRUE, trt)]) else NA,
    control = as.character(control)
  )
  attr(out, "n") <- c(treatment = n_T, control = length(trt) - n_T)
  attr(out, "events") <- c(treatment = e_T, control = sum(event) - e_T)
  attr(out, "n_strata") <- if (is.null(strata)) 1L else length(unique(strata))
  attr(out, "n_removed") <- n_removed
  out
}
