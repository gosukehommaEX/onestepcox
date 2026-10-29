#' Log hazard ratio estimates and standard errors of several methods
#'
#' Internal workhorse of \code{EstHR()}. The risk-set summaries, the anchor
#' estimates, and the table of tie blocks are computed once and shared by all
#' requested methods.
#'
#' @param time,event,strata,ties,is_sorted,tol As in \code{EstHR()}, after
#'   the checks of \code{check_inputs()}.
#' @param trt Logical treatment indicator returned by \code{check_inputs()}.
#' @param method Character vector of distinct valid methods.
#'
#' @return A matrix with one row per method (named by the method) and the
#'   columns \code{coef} and \code{se}, \code{NA} where the estimate cannot
#'   be computed.
#'
#' @noRd
estimate_hr <- function(time, event, trt, strata, method, ties, is_sorted,
                        tol) {
  out <- matrix(NA_real_, nrow = length(method), ncol = 2L,
                dimnames = list(method, c("coef", "se")))
  if (length(time) == 0L || sum(event) == 0) return(out)

  # Cox by survival::coxph()
  if ("Cox" %in% method) {
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
    out["Cox", ] <- c(unname(fit$coefficients[1L]),
                      sqrt(unname(fit$var[1L, 1L])))
  }
  others <- setdiff(method, "Cox")
  if (length(others) == 0L) return(out)

  # All other methods share the risk-set summaries
  rs <- risk_sets(time, event, trt, strata, is_sorted)
  O_T <- sum(rs$O_T)
  if (O_T == 0 || O_T == sum(rs$O)) return(out)

  closed <- c("MH", "Pike", "Peto")
  anchor_of <- function(m) {
    if (m %in% closed) m else if (startsWith(m, "Pike")) "Pike" else "MH"
  }
  needed <- unique(vapply(others, anchor_of, character(1)))
  anchors <- lapply(stats::setNames(nm = needed),
                    function(a) anchor_estimate(rs, a))
  dt <- NULL

  for (m in others) {
    a <- anchors[[anchor_of(m)]]
    if (m %in% closed) {
      out[m, ] <- c(a[["eta"]], a[["se"]])
      next
    }
    if (!is.finite(a[["eta"]])) next
    if (is.null(dt)) dt <- deriv_table(rs, ties)
    if (m == "Cox_NR") {
      fit <- newton_cox(a[["eta"]], dt, tol)
      if (!fit$converged) {
        warning("Newton-Raphson iteration of 'Cox_NR' did not converge; ",
                "NA is returned")
        next
      }
      out[m, ] <- c(fit$eta, 1 / sqrt(fit$info))
      next
    }
    d0 <- score_derivs(a[["eta"]], dt)
    if (!isTRUE(d0[["I"]] > 0)) next
    st <- one_step(a[["eta"]], d0, sub("^.*_", "", m))
    out[m, ] <- c(st[["eta"]], st[["se"]])
  }
  out
}
