#' Compare operating characteristics of HR estimators across scenarios
#'
#' Runs Monte Carlo simulations across one or more two-group trial scenarios
#' and compares the estimators of \code{\link{EstHR}}: their deviation from a
#' reference estimator (by default the Cox estimate) on the log scale, the
#' ratio of their standard errors to that of the reference, their bias, root
#' mean squared error (RMSE), and the coverage of their confidence intervals
#' relative to the true hazard ratio, and the computing time. Each scenario
#' specifies its own design (sample sizes, hazard, true HR, accrual length,
#' target number of events, and optionally dropout, strata, and the rounding
#' of times). The function uses \code{\link{GenSimData}} to generate the
#' replicate datasets and \code{\link{EstHR}} to compute the estimates.
#'
#' @param scenarios A data frame with one row per scenario. Required
#'   columns: \code{label} (character or factor), \code{n_T}
#'   (treatment-group size), \code{n_C} (control-group size),
#'   \code{hazard_C}, \code{HR}, \code{a}, \code{E}, as in
#'   \code{GenSimData()}. Optional columns: \code{dropout_hazard} and
#'   \code{time_unit} (both default to 0), and the list columns
#'   \code{strata_prop} and \code{strata_hazard}, whose elements are passed
#'   to \code{GenSimData()} (an element \code{NULL} means no strata). When a
#'   scenario has strata, the estimates are stratified.
#' @param methods Character vector of methods to compare, each a value of
#'   \code{method} in \code{EstHR()}.
#' @param nsim Integer. Number of Monte Carlo replicates per scenario.
#' @param seed Optional integer. Base seed; scenario \code{i} uses
#'   \code{seed + 100 * (i - 1)}.
#' @param reference Character string naming the reference method, which must
#'   be one of \code{methods}. Default \code{"Cox"}.
#' @param ties Handling of tied event times passed to \code{EstHR()},
#'   \code{"efron"} (default) or \code{"breslow"}.
#' @param conf.level Confidence level of the intervals. Default 0.95.
#'
#' @return A list with two components.
#'   \describe{
#'     \item{estimates}{A data frame in long format with one row per
#'       scenario, replicate, and method, and the columns \code{label},
#'       \code{simID}, \code{method}, \code{coef} (log HR), \code{se}
#'       (standard error of the log HR), and \code{lower} and \code{upper}
#'       (confidence limits of the HR).}
#'     \item{summary}{A data frame with one row per scenario and method,
#'       and the columns \code{label}, \code{HR_true}, \code{method},
#'       \code{n_NA} (replicates without an estimate); \code{dev_median},
#'       \code{dev_q90}, \code{dev_q99}, and \code{dev_max} (median, 90th
#'       and 99th percentiles, and maximum of the absolute difference between
#'       \code{coef} and that of the reference); \code{se_ratio_median},
#'       \code{se_ratio_min}, and \code{se_ratio_max} (of \code{se} divided
#'       by that of the reference); \code{bias} and \code{rmse} (of
#'       \code{coef} relative to \code{log(HR_true)}); \code{coverage}
#'       (proportion of intervals containing \code{HR_true}); and
#'       \code{time_per_sim} (seconds per replicate).}
#'   }
#'
#' @details
#' All replicates are summarized; no estimate is trimmed. Replicates in
#' which a method gives no estimate are counted in \code{n_NA} and left out
#' of that method's summaries, and those in which the reference gives no
#' estimate are left out of the deviations and standard-error ratios. The
#' percentiles are those of \code{stats::quantile()} with its default
#' \code{type = 7}. The computing time is the elapsed time of all replicates
#' of a method divided by \code{nsim}; it includes the grouping of the data
#' by replicate and is meant as a rough guide only.
#'
#' @examples
#' scenarios <- data.frame(
#'   label    = "demo",
#'   n_T      = 100,
#'   n_C      = 100,
#'   hazard_C = log(2) / 18,
#'   HR       = 0.7,
#'   a        = 12,
#'   E        = 110
#' )
#' res <- CompareOC(
#'   scenarios = scenarios,
#'   methods   = c("Cox", "Peto", "Pike", "MH_Halley"),
#'   nsim      = 50,
#'   seed      = 1
#' )
#' res$summary
#'
#' # A stratified scenario with times in whole months
#' scenarios$strata_prop <- list(c(0.5, 0.5))
#' scenarios$strata_hazard <- list(c(1, 3))
#' scenarios$time_unit <- 1
#' CompareOC(scenarios, c("Cox", "MH_Halley"), nsim = 20, seed = 1)$summary
#'
#' @import data.table
#' @export
CompareOC <- function(scenarios, methods, nsim, seed = NULL,
                      reference = "Cox", ties = "efron", conf.level = 0.95) {

  required_cols <- c("label", "n_T", "n_C", "hazard_C", "HR", "a", "E")
  missing_cols  <- setdiff(required_cols, names(scenarios))
  if (length(missing_cols) > 0L) {
    stop("scenarios is missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }
  if (!reference %in% methods) {
    stop("'reference' (", reference, ") must be one of 'methods'")
  }

  # Optional scenario settings: a column, or a list column, or the default
  get_opt <- function(i, col, default) {
    if (!col %in% names(scenarios)) return(default)
    v <- scenarios[[col]]
    if (is.list(v)) v[[i]] else v[i]
  }

  # Move the reference method to the front so it is computed first
  methods_ord <- c(reference, setdiff(methods, reference))

  base_seed <- if (is.null(seed)) NULL else as.integer(seed)

  est_rows  <- vector("list", nrow(scenarios))
  summ_rows <- vector("list", nrow(scenarios))

  for (i in seq_len(nrow(scenarios))) {
    sc <- scenarios[i, ]
    seed_i <- if (is.null(base_seed)) NULL else base_seed + 100L * (i - 1L)

    dat <- GenSimData(
      nsim           = nsim,
      n_T            = as.integer(sc$n_T),
      n_C            = as.integer(sc$n_C),
      hazard_C       = sc$hazard_C,
      HR             = sc$HR,
      a              = sc$a,
      E              = as.integer(sc$E),
      dropout_hazard = get_opt(i, "dropout_hazard", 0),
      strata_prop    = get_opt(i, "strata_prop", NULL),
      strata_hazard  = get_opt(i, "strata_hazard", NULL),
      time_unit      = get_opt(i, "time_unit", 0),
      seed           = seed_i
    )
    has_strata <- "stratum" %in% names(dat)

    # Run each method on every replicate; the reference goes first
    by_method <- vector("list", length(methods_ord))
    names(by_method) <- methods_ord
    time_by_method <- numeric(length(methods_ord))
    names(time_by_method) <- methods_ord
    for (m in methods_ord) {
      tm <- system.time({
        ests <- dat[, {
          e <- EstHR(time = time, event = event, group = group,
                     control = "C", method = m, ties = ties,
                     strata = if (has_strata) stratum else NULL,
                     conf.level = conf.level, is_sorted = TRUE)
          list(coef = e[[1L]], se = e[[3L]], lower = e[[4L]],
               upper = e[[5L]])
        }, by = simID]
      })[["elapsed"]]
      by_method[[m]] <- ests
      time_by_method[m] <- tm / nsim
    }

    est_rows[[i]] <- do.call(rbind, lapply(methods_ord, function(m) {
      data.frame(label = as.character(sc$label), simID = by_method[[m]]$simID,
                 method = m, coef = by_method[[m]]$coef,
                 se = by_method[[m]]$se, lower = by_method[[m]]$lower,
                 upper = by_method[[m]]$upper, stringsAsFactors = FALSE)
    }))

    # Summary table (one row per method)
    ref <- by_method[[reference]]
    eta_true <- log(sc$HR)
    summ_rows[[i]] <- do.call(rbind, lapply(methods_ord, function(m) {
      x <- by_method[[m]]
      ok <- !is.na(x$coef)
      dev <- abs(x$coef - ref$coef)
      dev <- dev[!is.na(dev)]
      ratio <- x$se / ref$se
      ratio <- ratio[!is.na(ratio)]
      covered <- x$lower <= sc$HR & sc$HR <= x$upper
      q <- if (length(dev) > 0L) {
        stats::quantile(dev, c(0.5, 0.9, 0.99), names = FALSE)
      } else {
        rep(NA_real_, 3L)
      }
      data.frame(
        label           = as.character(sc$label),
        HR_true         = sc$HR,
        method          = m,
        n_NA            = sum(!ok),
        dev_median      = q[1L],
        dev_q90         = q[2L],
        dev_q99         = q[3L],
        dev_max         = if (length(dev) > 0L) max(dev) else NA_real_,
        se_ratio_median = if (length(ratio) > 0L) stats::median(ratio) else NA_real_,
        se_ratio_min    = if (length(ratio) > 0L) min(ratio) else NA_real_,
        se_ratio_max    = if (length(ratio) > 0L) max(ratio) else NA_real_,
        bias            = if (any(ok)) mean(x$coef[ok] - eta_true) else NA_real_,
        rmse            = if (any(ok)) sqrt(mean((x$coef[ok] - eta_true)^2)) else NA_real_,
        coverage        = if (any(!is.na(covered))) mean(covered, na.rm = TRUE) else NA_real_,
        time_per_sim    = unname(time_by_method[m]),
        stringsAsFactors = FALSE
      )
    }))
  }

  out_est <- do.call(rbind, est_rows)
  out_summ <- do.call(rbind, summ_rows)
  rownames(out_est) <- NULL
  rownames(out_summ) <- NULL
  list(estimates = out_est, summary = out_summ)
}
