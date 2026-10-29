#' Summary of one method in one scenario of CompareOC()
#'
#' Internal workhorse of \code{CompareOC()}: summarizes the estimates of one
#' method in the replicates of one scenario, against the reference method and
#' the true hazard ratio.
#'
#' @param x Data frame (or data.table) of the estimates of the method, one
#'   row per replicate, with the columns \code{coef}, \code{se},
#'   \code{lower}, and \code{upper}.
#' @param ref Data frame of the estimates of the reference method, with the
#'   rows in the same order of the replicates as \code{x}.
#' @param HR_true True hazard ratio.
#' @param time_per_sim Computing time per replicate in seconds.
#'
#' @return A data frame with one row and the columns \code{n_NA},
#'   \code{dev_median}, \code{dev_q90}, \code{dev_q99}, \code{dev_max},
#'   \code{devse_median}, \code{devse_q99}, \code{devse_max},
#'   \code{mismatch_2dp}, \code{se_ratio_median}, \code{se_ratio_min},
#'   \code{se_ratio_max}, \code{bias}, \code{sd}, \code{rmse},
#'   \code{se_mean}, \code{coverage}, and \code{time_per_sim}, as described
#'   in \code{CompareOC()}. A summary without any value is \code{NA}.
#'
#' @noRd
summarize_oc <- function(x, ref, HR_true, time_per_sim) {
  eta_true <- log(HR_true)
  ok <- !is.na(x$coef)
  # Replicates in which both the method and the reference give an estimate
  dev <- abs(x$coef - ref$coef)
  devse <- dev / ref$se
  both <- ok & !is.na(ref$coef)
  dev <- dev[!is.na(dev)]
  devse <- devse[!is.na(devse)]
  ratio <- x$se / ref$se
  ratio <- ratio[!is.na(ratio)]
  # Hazard ratios written with two decimals; strings are compared, so that
  # no tolerance is involved
  mismatch <- sprintf("%.2f", exp(x$coef[both])) !=
    sprintf("%.2f", exp(ref$coef[both]))
  covered <- x$lower <= HR_true & HR_true <= x$upper
  se_ok <- x$se[ok & !is.na(x$se)]
  quant <- function(v, p) {
    if (length(v) > 0L) {
      stats::quantile(v, p, names = FALSE)
    } else {
      rep(NA_real_, length(p))
    }
  }
  q <- quant(dev, c(0.5, 0.9, 0.99))
  qs <- quant(devse, c(0.5, 0.99))
  data.frame(
    n_NA            = sum(!ok),
    dev_median      = q[1L],
    dev_q90         = q[2L],
    dev_q99         = q[3L],
    dev_max         = if (length(dev) > 0L) max(dev) else NA_real_,
    devse_median    = qs[1L],
    devse_q99       = qs[2L],
    devse_max       = if (length(devse) > 0L) max(devse) else NA_real_,
    mismatch_2dp    = if (any(both)) mean(mismatch) else NA_real_,
    se_ratio_median = if (length(ratio) > 0L) stats::median(ratio) else NA_real_,
    se_ratio_min    = if (length(ratio) > 0L) min(ratio) else NA_real_,
    se_ratio_max    = if (length(ratio) > 0L) max(ratio) else NA_real_,
    bias            = if (any(ok)) mean(x$coef[ok] - eta_true) else NA_real_,
    sd              = if (sum(ok) >= 2L) stats::sd(x$coef[ok]) else NA_real_,
    rmse            = if (any(ok)) sqrt(mean((x$coef[ok] - eta_true)^2)) else NA_real_,
    se_mean         = if (length(se_ok) > 0L) mean(se_ok) else NA_real_,
    coverage        = if (any(!is.na(covered))) mean(covered, na.rm = TRUE) else NA_real_,
    time_per_sim    = time_per_sim
  )
}
