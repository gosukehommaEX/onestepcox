#' Print the large-sample limits of hazard ratio estimators
#'
#' Prints the design and, for each true hazard ratio, the calendar time of
#' the analysis in the limit and the large-sample bias of each estimator of
#' the log hazard ratio (\code{limit - log(HR)}), as a table with one row
#' per method and one column per hazard ratio.
#'
#' @param x An object of class \code{"LimitHR"} returned by
#'   \code{LimitHR()}.
#' @param digits Number of significant digits. Default 3.
#' @param ... Not used.
#'
#' @return \code{x}, invisibly. \code{as.data.frame(x)} gives the data in
#'   long format.
#'
#' @seealso \code{\link{LimitHR}}, \code{\link{plot.LimitHR}}.
#'
#' @examples
#' lim <- LimitHR(HR = c(0.5, 2, 10), hazard_C = log(2) / 18,
#'                n_T = 100, n_C = 100, a = 12, E = 160)
#' print(lim)
#' print(lim, digits = 6)
#'
#' @export
print.LimitHR <- function(x, digits = 3L, ...) {
  if (!all(c("HR", "analysis_time", "method", "bias") %in% names(x))) {
    # Not the complete output of LimitHR(): print as a data frame
    y <- x
    class(y) <- "data.frame"
    print(y, digits = digits)
    return(invisible(x))
  }
  hr <- sort(unique(x$HR))
  methods <- unique(x$method)
  hr_lab <- paste("HR =", formatC(hr, digits = digits, format = "fg"))
  bias <- matrix(NA_real_, nrow = length(methods), ncol = length(hr),
                 dimnames = list(methods, hr_lab))
  for (i in seq_along(methods)) {
    sel <- x$method == methods[i]
    bias[i, ] <- x$bias[sel][match(hr, x$HR[sel])]
  }
  tau <- x$analysis_time[match(hr, x$HR)]

  cat("Large-sample limits of the log hazard ratio estimators\n")
  design <- attr(x, "design")
  if (!is.null(design)) {
    f <- function(v) formatC(v, digits = digits, format = "fg")
    cat("Design: hazard_C = ", f(design[["hazard_C"]]), ", n_T : n_C = ",
        f(design[["n_T"]]), " : ", f(design[["n_C"]]), ", a = ",
        f(design[["a"]]), ", E = ", f(design[["E"]]), ", dropout hazard = ",
        f(design[["dropout_hazard"]]), "\n", sep = "")
  }
  cat("\nAnalysis time:\n")
  print(noquote(stats::setNames(
    formatC(tau, digits = digits, format = "g"), hr_lab)), right = TRUE)
  cat("\nBias (limit - log HR):\n")
  print(noquote(formatC(bias, digits = digits, format = "g")), right = TRUE)
  invisible(x)
}
