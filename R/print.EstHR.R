#' Print hazard ratio estimates
#'
#' Prints the estimates of \code{\link{EstHR}} in the layout of
#' \code{summary()} of a \code{coxph} fit: the numbers of subjects and events
#' by group, a table of \code{coef}, \code{exp(coef)}, \code{se(coef)}, the
#' Wald statistic \code{z = coef / se(coef)} and its two-sided p-value, a
#' table of the hazard ratio and its confidence limits, and a description of
#' each method.
#'
#' @param x An object of class \code{"EstHR"} returned by \code{EstHR()}.
#' @param digits Number of significant digits. Default
#'   \code{max(3, getOption("digits") - 3)}.
#' @param ... Not used.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{EstHR}}, \code{\link{plot.EstHR}}.
#'
#' @examples
#' est <- EstHR(Surv(rfstime, status) ~ hormon, data = survival::gbsg,
#'              method = c("MH_Halley", "Cox"))
#' print(est)
#' print(est, digits = 6)
#'
#' @export
print.EstHR <- function(x, digits = max(3L, getOption("digits") - 3L), ...) {
  method <- attr(x, "method")
  groups <- attr(x, "groups")
  n <- attr(x, "n")
  events <- attr(x, "events")
  cols <- if (is.matrix(x)) colnames(x) else names(x)
  if (is.null(method) || is.null(groups) || is.null(n) || is.null(events) ||
      length(cols) != 5L || length(x) != 5L * length(method)) {
    # Not the complete output of EstHR(): print the numbers only
    keep <- intersect(names(attributes(x)), c("names", "dim", "dimnames"))
    y <- x
    attributes(y) <- attributes(x)[keep]
    print(y, digits = digits)
    return(invisible(x))
  }
  est <- matrix(as.vector(x), nrow = length(method),
                dimnames = list(method, cols))

  cat("Hazard ratio of the treatment group relative to the control group\n\n")
  lab <- format(c(paste0("Treatment (", groups[["treatment"]], "):"),
                  paste0("Control (", groups[["control"]], "):")))
  cat("  ", lab[1L], "  n = ", n[["treatment"]], ", events = ",
      events[["treatment"]], "\n", sep = "")
  cat("  ", lab[2L], "  n = ", n[["control"]], ", events = ",
      events[["control"]], "\n", sep = "")
  n_strata <- attr(x, "n_strata")
  n_removed <- attr(x, "n_removed")
  info <- paste0("Ties: ", attr(x, "ties"))
  if (!is.null(n_strata) && n_strata > 1L) {
    info <- c(info, paste0("strata: ", n_strata))
  }
  if (!is.null(n_removed) && n_removed > 0L) {
    info <- c(info, paste0(n_removed, " observation",
                           if (n_removed > 1L) "s", " deleted due to ",
                           "missingness"))
  }
  cat("  ", paste(info, collapse = "; "), "\n\n", sep = "")

  if (all(is.na(est[, 1L]))) {
    cat("Estimate not available (no event, or no event in one group)\n")
  } else {
    z <- est[, "coef"] / est[, "se(coef)"]
    tab <- cbind(est[, 1:3, drop = FALSE], z = z,
                 "Pr(>|z|)" = 2 * stats::pnorm(-abs(z)))
    stats::printCoefmat(tab, digits = digits, signif.stars = FALSE,
                        P.values = TRUE, has.Pvalue = TRUE, na.print = "NA")
    cat("\n")
    print(est[, c(2L, 4L, 5L), drop = FALSE], digits = digits)
  }

  cat("\nMethods:\n")
  cat(paste0("  ", format(method), "  ", method_label(method), "\n"),
      sep = "")
  invisible(x)
}
