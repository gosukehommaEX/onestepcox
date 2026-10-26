#' onestepcox: One-Step Closed-Form Hazard Ratio Estimators
#'
#' Closed-form estimators of the hazard ratio of a two-group comparison that
#' approximate the Cox partial likelihood estimator without iteration. The
#' proposed estimator takes one Halley step on the Cox score equation from the
#' Mantel-Haenszel estimator, and its standard error is updated in the same
#' pass over the risk sets. The Efron and Breslow handling of tied event times
#' and stratification are supported.
#'
#' The main functions are
#' \describe{
#'   \item{\code{\link{EstHR}}}{Hazard ratio, standard error, and confidence
#'     interval by the proposed estimator, other one-step corrections, the
#'     Mantel-Haenszel, Pike, and Peto estimators, and the Cox estimator.}
#'   \item{\code{\link{RiskTable}}}{Risk-set summaries shared by the
#'     estimators.}
#'   \item{\code{\link{GenSimData}}}{Simulated two-group trial data.}
#'   \item{\code{\link{CompareOC}}}{Operating characteristics of the
#'     estimators across simulation scenarios.}
#' }
#'
#' @importFrom utils globalVariables
#' @keywords internal
"_PACKAGE"

# Suppress R CMD check NOTEs of the form
# "no visible binding for global variable ..." emitted by data.table's
# non-standard evaluation inside GenSimData() and CompareOC().
utils::globalVariables(c(
  # data.table's .() syntactic sugar for list()
  ".",
  # GenSimData.R
  "accrual_time", "surv_time", "dropout_time", "tte", "surv_event",
  "total", "analysis_time", "event", "time", "simID", "group",
  # CompareOC.R
  "HR"
))
