#' Print operating characteristics of hazard ratio estimators
#'
#' Prints, for each scenario of \code{\link{CompareOC}}, a table with one row
#' per method of the number of replicates without an estimate, the median
#' and maximum deviation from the reference estimator, the median ratio of
#' the standard errors, the bias, RMSE, and coverage, and the computing time
#' per replicate. All columns are in the component \code{summary}.
#'
#' @param x An object of class \code{"CompareOC"} returned by
#'   \code{CompareOC()}.
#' @param digits Number of significant digits. Default 3.
#' @param ... Not used.
#'
#' @return \code{x}, invisibly.
#'
#' @seealso \code{\link{CompareOC}}, \code{\link{plot.CompareOC}}.
#'
#' @examples
#' scenarios <- data.frame(label = "demo", n_T = 100, n_C = 100,
#'                         hazard_C = log(2) / 18, HR = 0.7, a = 12,
#'                         E = 110)
#' res <- CompareOC(scenarios, c("Cox", "Pike", "MH_Halley"), nsim = 20,
#'                  seed = 1)
#' print(res, digits = 2)
#'
#' @export
print.CompareOC <- function(x, digits = 3L, ...) {
  s <- x$summary
  cols <- c("n_NA", "dev_median", "dev_max", "se_ratio_median", "bias",
            "rmse", "coverage", "time_per_sim")
  if (!is.data.frame(s) ||
      !all(c("label", "HR_true", "method", cols) %in% names(s))) {
    # Not the complete output of CompareOC(): print the list itself
    print(x[seq_along(x)], digits = digits)
    return(invisible(x))
  }
  labels <- unique(as.character(s$label))
  cat("Operating characteristics of hazard ratio estimators\n")
  info <- c(paste0(length(labels), " scenario",
                   if (length(labels) > 1L) "s"),
            if (!is.null(attr(x, "nsim"))) {
              paste0(attr(x, "nsim"), " replicates each")
            },
            if (!is.null(attr(x, "reference"))) {
              paste0("reference: ", attr(x, "reference"))
            },
            if (!is.null(attr(x, "ties"))) paste0("ties: ", attr(x, "ties")),
            if (!is.null(attr(x, "conf.level"))) {
              paste0("conf.level: ", attr(x, "conf.level"))
            })
  cat(paste(info, collapse = "; "), "\n", sep = "")

  for (lab in labels) {
    r <- s[as.character(s$label) == lab, , drop = FALSE]
    cat("\nScenario ", lab, " (true HR = ", format(r$HR_true[1L]), ")\n",
        sep = "")
    tab <- r[, cols]
    names(tab) <- c("n_NA", "dev_med", "dev_max", "se_ratio", "bias",
                    "rmse", "coverage", "time")
    rownames(tab) <- r$method
    print(tab, digits = digits)
  }
  cat("\ndev: |coef - coef of the reference| (median and maximum);",
      "se_ratio: median of\nse / se of the reference; time: seconds per",
      "replicate. All summaries are in x$summary.\n")
  invisible(x)
}
