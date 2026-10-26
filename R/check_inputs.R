#' Check the data arguments shared by EstHR() and RiskTable()
#'
#' @param time,event,group,control,strata As in \code{EstHR()}.
#'
#' @return The logical treatment indicator \code{group != control}. An error
#'   is raised for unequal lengths, missing values, event codes other than 0
#'   and 1, and a \code{control} that matches no value of \code{group}.
#'
#' @noRd
check_inputs <- function(time, event, group, control, strata) {
  n <- length(time)
  if (length(event) != n || length(group) != n) {
    stop("'time', 'event', and 'group' must have the same length")
  }
  if (!is.null(strata) && length(strata) != n) {
    stop("'strata' must have the same length as 'time'")
  }
  if (length(control) != 1L) {
    stop("'control' must be a single value")
  }
  if (!is.numeric(time)) {
    stop("'time' must be numeric")
  }
  if (anyNA(time) || anyNA(event) || anyNA(group) ||
      (!is.null(strata) && anyNA(strata))) {
    stop("missing values in 'time', 'event', 'group', or 'strata' ",
         "are not allowed")
  }
  if (any(event != 0 & event != 1)) {
    stop("'event' must be coded 1 (event) and 0 (censored)")
  }
  trt <- group != control
  if (n > 0L && all(trt)) {
    stop("'control' does not match any value of 'group'")
  }
  trt
}
