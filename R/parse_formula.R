#' Evaluate a formula of the form Surv(time, event) ~ group + strata(s)
#'
#' Internal helper of the formula interface of \code{EstHR()} and
#' \code{RiskTable()}. The response is read from the call to \code{Surv()}
#' without calling it, and an event indicator coded 1 and 2 (2 for an event)
#' is recoded to 0 and 1, as by \code{Surv()}. The right-hand side must have
#' one group term (a variable or an expression such as \code{factor(x)}) and
#' may add one or more \code{strata()} terms; several variables in
#' \code{strata()} (or several \code{strata()} terms) define strata by their
#' combinations. Rows with a missing value in any of the variables are
#' removed, as in \code{survival::coxph()}.
#'
#' @param formula A formula \code{Surv(time, event) ~ group}, optionally with
#'   \code{+ strata(s)}.
#' @param data A data frame or list in which the variables are evaluated, or
#'   \code{NULL} for the environment of the formula.
#' @param event The value given as the second argument of the calling
#'   function, or \code{NULL}; a data frame or list given there is used as
#'   \code{data}.
#' @param group,strata The values given as \code{group} and \code{strata},
#'   which must be \code{NULL} with a formula.
#'
#' @return A list with elements \code{time}, \code{event}, \code{group},
#'   \code{strata} (\code{NULL} without \code{strata()}), \code{control} (the
#'   first level of the group variable, or its smallest value), and
#'   \code{n_removed} (number of rows removed for missing values).
#'
#' @noRd
parse_formula <- function(formula, data = NULL, event = NULL, group = NULL,
                          strata = NULL) {
  # A data frame given as the second argument is the data, as in
  # coxph(formula, data)
  if (!is.null(event)) {
    if (is.null(data) && is.list(event)) {
      data <- event
    } else {
      stop("with a formula, the second argument must be the data frame ",
           "or be omitted")
    }
  }
  if (!is.null(group)) {
    stop("with a formula, the group variable is given in the formula, ",
         "and the other arguments are given by name")
  }
  if (!is.null(strata)) {
    stop("with a formula, strata are given by strata() in the formula")
  }
  if (!is.null(data) && !is.list(data)) {
    stop("'data' must be a data frame or a list")
  }
  if (length(formula) != 3L) {
    stop("the formula must have the form Surv(time, event) ~ group")
  }

  env <- environment(formula)
  if (is.null(env)) env <- parent.frame(2L)
  # Name of the function called by an expression, also for pkg::fun
  head_name <- function(e) {
    h <- e[[1L]]
    if (is.name(h)) return(as.character(h))
    if (is.call(h) && length(h) == 3L && is.name(h[[1L]]) &&
        as.character(h[[1L]]) %in% c("::", ":::")) {
      return(as.character(h[[3L]]))
    }
    ""
  }

  # Response: Surv(time, event), with positional or named arguments
  lhs <- formula[[2L]]
  if (!is.call(lhs) || head_name(lhs) != "Surv") {
    stop("the response of the formula must be Surv(time, event)")
  }
  s_args <- as.list(lhs)[-1L]
  s_nm <- names(s_args)
  if (is.null(s_nm)) s_nm <- rep("", length(s_args))
  if (length(s_args) != 2L || !all(s_nm %in% c("", "time", "event")) ||
      anyDuplicated(s_nm[s_nm != ""]) > 0L) {
    stop("the response must be Surv(time, event) for right-censored data")
  }
  i_time <- match("time", s_nm)
  i_event <- match("event", s_nm)
  free <- which(s_nm == "")
  if (is.na(i_time)) {
    i_time <- free[1L]
    free <- free[-1L]
  }
  if (is.na(i_event)) i_event <- free[1L]

  # Right-hand side: one group term and optional strata() terms. Formula
  # operators other than + and the special terms of coxph() are not allowed.
  split_plus <- function(e) {
    if (is.call(e) && identical(e[[1L]], as.name("+")) && length(e) == 3L) {
      c(split_plus(e[[2L]]), split_plus(e[[3L]]))
    } else {
      list(e)
    }
  }
  not_allowed <- c("+", "*", ":", "-", "^", "/", "%in%", "|", "cluster",
                   "tt", "frailty", "offset", "pspline", "ridge")
  bad_rhs <- paste("the right-hand side of the formula must be one group",
                   "variable, optionally with strata()")
  g_expr <- NULL
  s_exprs <- list()
  for (e in split_plus(formula[[3L]])) {
    if (is.call(e) && head_name(e) == "strata") {
      a <- as.list(e)[-1L]
      if (length(a) == 0L || (!is.null(names(a)) && any(names(a) != ""))) {
        stop("strata() in the formula must list variables only")
      }
      s_exprs <- c(s_exprs, a)
    } else if ((is.name(e) && !identical(as.character(e), ".")) ||
               (is.call(e) && !head_name(e) %in% not_allowed)) {
      if (!is.null(g_expr)) stop(bad_rhs)
      g_expr <- e
    } else {
      stop(bad_rhs)
    }
  }
  if (is.null(g_expr)) {
    stop("the right-hand side of the formula must name the group variable")
  }

  # Values of the variables; the class "AsIs" of I() is removed
  eval_var <- function(e) {
    v <- if (is.null(data)) eval(e, env) else eval(e, data, env)
    if (inherits(v, "AsIs")) oldClass(v) <- setdiff(oldClass(v), "AsIs")
    v
  }
  time <- eval_var(s_args[[i_time]])
  event <- eval_var(s_args[[i_event]])
  # Event indicator coded 1 and 2, as accepted by Surv()
  if (is.numeric(event) && any(!is.na(event)) &&
      max(event, na.rm = TRUE) == 2) {
    event <- event - 1
  }
  group <- eval_var(g_expr)
  s_vals <- lapply(s_exprs, eval_var)
  n <- length(time)
  if (length(event) != n || length(group) != n ||
      any(vapply(s_vals, length, integer(1)) != n)) {
    stop("the variables of the formula must have the same length")
  }

  # Remove rows with missing values
  keep <- !(is.na(time) | is.na(event) | is.na(group))
  for (v in s_vals) keep <- keep & !is.na(v)
  time <- time[keep]
  event <- event[keep]
  group <- group[keep]
  s_vals <- lapply(s_vals, function(v) v[keep])
  strata <- if (length(s_vals) == 0L) {
    NULL
  } else if (length(s_vals) == 1L) {
    s_vals[[1L]]
  } else {
    do.call(paste, c(s_vals, sep = ":"))
  }

  lev <- if (is.factor(group)) {
    levels(droplevels(group))
  } else {
    sort(unique(group))
  }
  if (length(lev) != 2L) {
    stop("the group variable must have exactly two values")
  }
  list(time = time, event = event, group = group, strata = strata,
       control = lev[1L], n_removed = sum(!keep))
}
