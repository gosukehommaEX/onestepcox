#' Values plotted by plot.CompareOC()
#'
#' @param x An object of class \code{"CompareOC"} returned by
#'   \code{CompareOC()}.
#' @param what Quantity to plot, one of \code{"deviation"}, \code{"bias"},
#'   \code{"rmse"}, \code{"coverage"}, \code{"se_ratio"}, and \code{"time"}.
#' @param methods Methods to plot, or \code{NULL} for all methods, except the
#'   reference for \code{"deviation"} and \code{"se_ratio"}.
#' @param labels \code{NULL} for the default labels of \code{plot_labels()},
#'   or one label per method of \code{methods}.
#' @param floor Smallest deviation kept with \code{"deviation"}, a positive
#'   number; smaller deviations (including 0) are set to \code{floor}.
#'
#' @return A data frame with the columns \code{scenario}, \code{method},
#'   \code{value}, and \code{label}. For \code{"deviation"}, one row per
#'   scenario, method, and replicate with a deviation, where \code{value} is
#'   the absolute difference between \code{coef} and that of the reference in
#'   the same replicate; otherwise one row per scenario and method, where
#'   \code{value} is the summary (\code{NA} when it is missing), with the
#'   scenarios varying fastest. The attributes are \code{reference},
#'   \code{methods}, \code{labels} (one per method), and \code{scenarios}
#'   (all scenarios in the order of \code{x$summary}).
#'
#' @noRd
plot_data_CompareOC <- function(x, what, methods = NULL, labels = NULL,
                                floor = 1e-16) {
  s <- x$summary
  est <- x$estimates
  if (!is.data.frame(s) || !is.data.frame(est)) {
    stop("'x' must be the output of CompareOC()")
  }
  # The deviations are drawn on a logarithmic axis
  if (!is.numeric(floor) || length(floor) != 1L || !is.finite(floor) ||
      floor <= 0) {
    stop("'floor' must be a single positive number")
  }
  ref <- attr(x, "reference")
  if (is.null(ref)) ref <- s$method[1L]
  scen <- unique(as.character(s$label))
  all_methods <- unique(s$method)
  if (is.null(methods)) {
    methods <- if (what %in% c("deviation", "se_ratio")) {
      setdiff(all_methods, ref)
    } else {
      all_methods
    }
  }
  if (length(methods) == 0L) stop("there is no method to plot")
  if (!all(methods %in% all_methods)) {
    stop("'methods' must be methods of 'x': ",
         paste(all_methods, collapse = ", "))
  }
  lab <- plot_labels(methods, labels)

  if (what == "deviation") {
    # Deviation of each replicate from the reference, matched by simID
    rows <- list()
    for (sc in scen) {
      e <- est[as.character(est$label) == sc, , drop = FALSE]
      r <- e[e$method == ref, , drop = FALSE]
      for (m in methods) {
        em <- e[e$method == m, , drop = FALSE]
        dv <- abs(em$coef - r$coef[match(em$simID, r$simID)])
        dv <- dv[!is.na(dv)]
        dv[dv < floor] <- floor
        if (length(dv) > 0L) {
          rows[[length(rows) + 1L]] <- data.frame(
            scenario = sc, method = m, value = dv, stringsAsFactors = FALSE
          )
        }
      }
    }
    if (length(rows) == 0L) stop("there is no deviation to plot")
    out <- do.call(rbind, rows)
  } else {
    column <- c(bias = "bias", rmse = "rmse", coverage = "coverage",
                se_ratio = "se_ratio_median", time = "time_per_sim")[[what]]
    out <- do.call(rbind, lapply(methods, function(m) {
      sm <- s[s$method == m, , drop = FALSE]
      data.frame(scenario = scen, method = m,
                 value = sm[[column]][match(scen, as.character(sm$label))],
                 stringsAsFactors = FALSE)
    }))
    if (!any(is.finite(out$value))) stop("there is no value to plot")
  }
  out$label <- lab[match(out$method, methods)]
  rownames(out) <- NULL
  attr(out, "reference") <- ref
  attr(out, "methods") <- methods
  attr(out, "labels") <- lab
  attr(out, "scenarios") <- scen
  out
}
