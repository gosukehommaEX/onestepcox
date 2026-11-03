# Large-sample limits and simulations with event times grouped into whole
# months, in the design of grid B2 (1:1 allocation, no dropout, no strata,
# recorded times rounded up to months): Section S9 and Table S8 of the
# Supporting Information. With grouped times, ties persist as n grows, so
# that the Mantel-Haenszel estimator and the Cox estimator with the Efron or
# Breslow approximation converge to different limits, and the difference
# between the MH-Halley and Cox estimates does not vanish.
#
# Run from the root of the package (open onestepcox.Rproj) after
# devtools::load_all(), or with the package installed:
#   source("inst/sup_info/data_generate_grouped.R")
# For a quick test run, first run options(onestepcox.quick = TRUE); the
# results are then written to inst/sup_info/results/quick/. This script only
# computes and saves the results; table_and_figure_supplement.R makes
# Table S8 from them.

source(file.path("inst", "scripts", "settings.R"))
out <- results_dir("sup_info")
run_log <- list()
unit <- gridB_settings$B2$time_unit
pi_T <- 0.5

# 1. Large-sample limits ------------------------------------------------------
# Subjects enter uniformly over (0, a), and the analysis takes place at the
# calendar time tau at which the expected fraction of subjects with an event
# is event_frac (as in LimitHR()). Divided by n, the number at risk of group j
# in the k-th interval of time since entry, ((k - 1) u, k u] with u = unit,
# tends to
#   Y_j(k) = pi_j G((k - 1) u) exp(-lambda_j (k - 1) u),
#   G(t) = min(1, (tau - t) / a) for t <= tau,
# and the number of events of group j with recorded time k u to
#   O_j(k) = pi_j int_{(k - 1) u}^{min(k u, tau)} lambda_j exp(-lambda_j t) G(t) dt,
# because the recorded times of the events and of the censoring are rounded
# up to multiples of u (GenSimData()). The sums over the event times of
# Section 2 of the article become sums over the intervals. In the Efron
# approximation, the r-th of the O(k) terms of an interval reduces the
# numbers at risk by r / O(k) times the events of each group; as n grows,
# r / O(k) fills (0, 1), and the sum over r becomes an integral over this
# fraction s.
K_H <- function(x) (2 + 3 * x * exp(x)) / (6 * (3 - exp(x)))

grouped_tables <- function(theta) {
  lam <- c(theta * hazard_C, hazard_C)
  pi_j <- c(pi_T, 1 - pi_T)
  tau <- onestepcox:::limit_analysis_time(lam[1L], lam[2L], pi_T, accrual,
                                          event_frac, 0)
  G <- function(t) pmax(0, pmin(1, (tau - t) / accrual))
  # Integral over (lo, hi), split where G changes its form (t = tau - a)
  integ <- function(f, lo, hi) {
    br <- sort(unique(c(lo, hi, tau - accrual)))
    br <- br[br >= lo & br <= hi]
    sum(vapply(seq_len(length(br) - 1L), function(k) {
      stats::integrate(f, br[k], br[k + 1L], rel.tol = 1e-10,
                       abs.tol = 0)$value
    }, numeric(1)))
  }
  K <- ceiling(tau / unit)
  start <- (seq_len(K) - 1) * unit
  Y <- O <- matrix(0, nrow = K, ncol = 2L, dimnames = list(NULL, c("T", "C")))
  for (j in 1:2) {
    Y[, j] <- pi_j[j] * G(start) * exp(-lam[j] * start)
    O[, j] <- vapply(seq_len(K), function(k) {
      pi_j[j] * integ(function(t) lam[j] * exp(-lam[j] * t) * G(t),
                      start[k], min(start[k] + unit, tau))
    }, numeric(1))
  }
  list(theta = theta, tau = tau, Y = Y, O = O)
}

# Score U and its derivatives I = -U' and J = -U'' of the limit, divided by n
grouped_derivs <- function(eta, tb, ties) {
  Y <- tb$Y
  O <- tb$O
  o <- O[, "T"] + O[, "C"]
  if (ties == "breslow") {
    p <- Y[, "T"] * exp(eta) / (Y[, "C"] + Y[, "T"] * exp(eta))
    return(c(U = sum(O[, "T"]) - sum(o * p), I = sum(o * p * (1 - p)),
             J = sum(o * p * (1 - p) * (1 - 2 * p))))
  }
  # Efron: average over the fraction s of the events removed from the risk set
  avg <- function(k, f) {
    stats::integrate(function(s) {
      yt <- Y[k, "T"] - s * O[k, "T"]
      yc <- Y[k, "C"] - s * O[k, "C"]
      p <- yt * exp(eta) / (yc + yt * exp(eta))
      f(p)
    }, 0, 1, rel.tol = 1e-10, abs.tol = 0)$value
  }
  ks <- which(o > 0)
  ep <- vapply(ks, function(k) avg(k, function(p) p), numeric(1))
  epq <- vapply(ks, function(k) avg(k, function(p) p * (1 - p)), numeric(1))
  epqd <- vapply(ks, function(k) {
    avg(k, function(p) p * (1 - p) * (1 - 2 * p))
  }, numeric(1))
  c(U = sum(O[, "T"]) - sum(o[ks] * ep), I = sum(o[ks] * epq),
    J = sum(o[ks] * epqd))
}

# Root of the limiting score by Newton-Raphson iteration from eta0
grouped_root <- function(eta0, tb, ties) {
  eta <- eta0
  for (it in 1:100) {
    d <- grouped_derivs(eta, tb, ties)
    step <- d[["U"]] / d[["I"]]
    eta <- eta + step
    if (abs(step) < 1e-15) break
  }
  if (abs(step) >= 1e-12) stop("the limit of the Cox estimate did not converge")
  eta
}

grouped_limit_row <- function(theta) {
  tb <- grouped_tables(theta)
  Y <- tb$Y
  O <- tb$O
  y <- Y[, "T"] + Y[, "C"]
  mh <- log(sum(O[, "T"] * Y[, "C"] / y) / sum(O[, "C"] * Y[, "T"] / y))
  res <- list(theta = theta, tau = tb$tau, intervals = nrow(Y),
              events = sum(O), mh = mh)
  for (ties in c("efron", "breslow")) {
    cox <- grouped_root(mh, tb, ties)
    d0 <- grouped_derivs(mh, tb, ties)
    halley <- mh + 2 * d0[["U"]] * d0[["I"]] /
      (2 * d0[["I"]]^2 + d0[["U"]] * d0[["J"]])
    b <- mh - cox
    res[[paste0("cox_", ties)]] <- cox
    res[[paste0("halley_", ties)]] <- halley
    res[[paste0("b_", ties)]] <- b
    res[[paste0("diff_", ties)]] <- halley - cox
    res[[paste0("bound_", ties)]] <- K_H(abs(b)) * abs(b)^3
  }
  as.data.frame(res)
}

t0 <- proc.time()[["elapsed"]]
message("[grouped_limits] large-sample limits with grouped times")
grouped_theta <- c(0.2, 0.5, 2, 5)
limits <- do.call(rbind, lapply(grouped_theta, grouped_limit_row))
run_log[[length(run_log) + 1L]] <-
  log_entry("grouped_limits", length(grouped_theta), 0L, t0)

# 2. Simulations with growing n -----------------------------------------------
# The design of grid B2 with the hazard ratio 5 and n_T = n_C, with the
# estimates of the Mantel-Haenszel anchor and the MH-Halley estimator and the
# reference Cox_NR (the root of the Cox score by Newton-Raphson iteration on
# the same risk-set summaries), for the Efron and the Breslow approximation on
# the same datasets. Scenario i uses the seed base + 100 (i - 1), and chunk k
# of a scenario run in chunks of 10 uses base + 100 (i - 1) + (k - 1), as in
# run_growth() of settings.R.
grouped_sim_theta <- 5
grouped_n <- c(200, 2000, 20000, 200000)
grouped_nsim <- c(2000L, 500L, 200L, 50L)
grouped_methods <- c("Cox_NR", "MH", "MH_Halley")
sim_sc <- data.frame(
  label = sprintf("grouped_HR%s_n%.0f", grouped_sim_theta, grouped_n),
  n = grouped_n, n_T = grouped_n / 2, n_C = grouped_n / 2,
  hazard_C = hazard_C, HR = grouped_sim_theta, a = accrual,
  E = E_of(grouped_n / 2, grouped_n / 2), dropout_hazard = 0,
  time_unit = unit, nsim = grouped_nsim, stringsAsFactors = FALSE
)

t0 <- proc.time()[["elapsed"]]
est_rows <- list()
sum_rows <- list()
for (i in seq_len(nrow(sim_sc))) {
  nsim <- n_rep(sim_sc$nsim[i])
  chunk <- if (sim_sc$n[i] >= 2e4) min(10L, nsim) else nsim
  starts <- seq(1L, nsim, by = chunk)
  for (ties in c("efron", "breslow")) {
    message(sprintf("[grouped_sim] %s, %s ties (%d replicates)",
                    sim_sc$label[i], ties, nsim))
    pieces <- vector("list", length(starts))
    for (k in seq_along(starts)) {
      size <- min(chunk, nsim - starts[k] + 1L)
      r <- onestepcox::CompareOC(
        sim_sc[i, c("label", "n_T", "n_C", "hazard_C", "HR", "a", "E",
                    "dropout_hazard", "time_unit")],
        grouped_methods, nsim = size,
        seed = seeds[["grouped"]] + 100L * (i - 1L) + (k - 1L),
        reference = "Cox_NR", ties = ties
      )
      e <- r$estimates[, c("label", "simID", "method", "coef")]
      e$rep <- starts[k] - 1L + e$simID
      e$simID <- NULL
      pieces[[k]] <- e
    }
    e <- do.call(rbind, pieces)
    e$ties <- ties
    est_rows[[length(est_rows) + 1L]] <- e
    w <- stats::reshape(e[, c("rep", "method", "coef")], idvar = "rep",
                        timevar = "method", direction = "wide")
    d_mh <- w$coef.MH - w$coef.Cox_NR
    d_h <- w$coef.MH_Halley - w$coef.Cox_NR
    # Replicates without a finite Cox estimate are left out; a Halley step that
    # is not defined (denominator not positive) is counted in n_halley_na
    ok <- !is.na(w$coef.Cox_NR)
    ok_h <- ok & !is.na(d_h)
    sum_rows[[length(sum_rows) + 1L]] <- data.frame(
      label = sim_sc$label[i], n = sim_sc$n[i], HR = sim_sc$HR[i],
      ties = ties, nsim = nsim, n_undefined = sum(!ok),
      n_halley_na = sum(ok & is.na(d_h)),
      mh_mean = mean(d_mh[ok]), mh_median_abs = stats::median(abs(d_mh[ok])),
      halley_mean = mean(d_h[ok_h]),
      halley_median_abs = stats::median(abs(d_h[ok_h])),
      halley_max_abs = max(abs(d_h[ok_h])), stringsAsFactors = FALSE
    )
  }
}
simulation <- do.call(rbind, sum_rows)
estimates <- do.call(rbind, est_rows)
rownames(simulation) <- NULL
rownames(estimates) <- NULL
# Each replicate is analyzed with both approximations
run_log[[length(run_log) + 1L]] <-
  log_entry("grouped_sim", nrow(sim_sc),
            2L * sum(vapply(sim_sc$nsim, n_rep, integer(1))), t0)

saveRDS(list(limits = limits, scenarios = sim_sc, simulation = simulation,
             estimates = estimates, time_unit = unit),
        file.path(out, "data", "grouped.rds"))
write_run_log(run_log, out, "run_log_grouped")
message("Finished. Results in ", file.path(out, "data", "grouped.rds"))
