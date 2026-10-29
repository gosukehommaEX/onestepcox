# Data of the figures and tables of the Supporting Information.
#
# Run from the root of the package (open onestepcox.Rproj) after
# devtools::load_all(), or with the package installed:
#   source("inst/sup_info/data_generate_supplement.R")
# For a quick test run, first run options(onestepcox.quick = TRUE); the
# results are then written to inst/sup_info/results/quick/. The settings are
# shared with the main text (inst/scripts/settings.R) and described in
# inst/scripts/README.md. This script only computes and saves the results;
# table_and_figure_supplement.R makes the tables and figures from them.

source(file.path("inst", "scripts", "settings.R"))
# Checked first, so that a missing package does not stop the run at the end
if (!requireNamespace("microbenchmark", quietly = TRUE)) {
  stop("The package 'microbenchmark' is needed for the computing times")
}
out <- results_dir("sup_info")
run_log <- list()

# 1. Large-sample limits with 2:1 allocation and with dropout of 10 per
#    cent per year (Figure S1)
t0 <- proc.time()[["elapsed"]]
message("[limit_sup] large-sample limits, 2:1 allocation and dropout")
hr_grid <- exp(seq(log(0.1), log(10), length.out = if (quick) 9L else 81L))
limit_sup <- list(
  r2 = onestepcox::LimitHR(HR = hr_grid, hazard_C = hazard_C, n_T = 2,
                           n_C = 1, a = accrual, E = 3 * event_frac),
  dropout = onestepcox::LimitHR(HR = hr_grid, hazard_C = hazard_C, n_T = 1,
                                n_C = 1, a = accrual, E = 2 * event_frac,
                                dropout_hazard = dropout_10)
)
saveRDS(limit_sup, file.path(out, "data", "limit_sup.rds"))
run_log[[length(run_log) + 1L]] <-
  log_entry("limit_sup", 2L * length(hr_grid), 0L, t0)

# 2. Grid B with ties (B2, B5) and the Breslow handling of ties for all
#    methods, including the reference Cox estimate (Table S3)
t0 <- proc.time()[["elapsed"]]
sc <- gridB_scenarios(c("B2", "B5"))
gridB_breslow <- list(scenarios = sc,
                      summary = run_grid(sc, seeds[["gridB_breslow"]],
                                         "gridB_breslow", ties = "breslow"))
saveRDS(gridB_breslow, file.path(out, "data", "gridB_breslow.rds"))
run_log[[length(run_log) + 1L]] <-
  log_entry("gridB_breslow", nrow(sc), nrow(sc) * n_rep(10000L), t0)

# 3. Grid C: small samples and extreme hazard ratios (Table S4). In some
#    replicates the Cox estimate is not finite (for example, when a group
#    has no event); all methods are then left out of the summaries
#    (n_undefined in the summary; see summarize_methods() in settings.R)
t0 <- proc.time()[["elapsed"]]
sc <- gridC_scenarios()
gridC <- list(scenarios = sc,
              summary = run_grid(sc, seeds[["gridC"]], "gridC"))
saveRDS(gridC, file.path(out, "data", "gridC.rds"))
run_log[[length(run_log) + 1L]] <-
  log_entry("gridC", nrow(sc), nrow(sc) * n_rep(10000L), t0)

# 4. Computing time (Table S5): one dataset for each n, HR 0.5, 1:1. The
#    risk-set summaries (aggregation) are timed apart from the closed-form
#    step from the Mantel-Haenszel anchor and from the Newton-Raphson
#    iteration from the same anchor on the same summaries, and the total time
#    of EstHR() and of survival::coxph() is timed as well. Internal functions
#    of the package are called with ':::'.
t0 <- proc.time()[["elapsed"]]
n_time <- c(200, 1000, 5000, 20000)
times <- if (quick) 5L else 200L
timing_rows <- list()
iter_rows <- list()
for (j in seq_along(n_time)) {
  n <- n_time[j]
  message(sprintf("[timing] n = %.0f", n))
  d <- onestepcox::GenSimData(nsim = 1, n_T = n / 2, n_C = n / 2,
                              hazard_C = hazard_C, HR = 0.5, a = accrual,
                              E = E_of(n / 2, n / 2),
                              seed = seeds[["timing"]] + (j - 1L))
  time <- d$time
  event <- d$event
  group <- d$group
  trt <- group == "T"
  x <- as.numeric(trt)
  rs <- onestepcox:::risk_sets(time, event, trt, NULL, is_sorted = TRUE)
  dt <- onestepcox:::deriv_table(rs, "efron")
  eta0 <- onestepcox:::anchor_estimate(rs, "MH")[["eta"]]
  mb <- microbenchmark::microbenchmark(
    aggregation = onestepcox:::deriv_table(
      onestepcox:::risk_sets(time, event, trt, NULL, is_sorted = TRUE),
      "efron"
    ),
    mh_halley_step = {
      e0 <- onestepcox:::anchor_estimate(rs, "MH")[["eta"]]
      onestepcox:::one_step(e0, onestepcox:::score_derivs(e0, dt), "Halley")
    },
    newton_1e8 = onestepcox:::newton_cox(eta0, dt, tol = 1e-8),
    newton_1e12 = onestepcox:::newton_cox(eta0, dt, tol = 1e-12),
    EstHR_MH_Halley = onestepcox::EstHR(time, event, group, "C",
                                        "MH_Halley", is_sorted = TRUE),
    EstHR_Cox_NR = onestepcox::EstHR(time, event, group, "C", "Cox_NR",
                                     is_sorted = TRUE),
    EstHR_Cox = onestepcox::EstHR(time, event, group, "C", "Cox",
                                  is_sorted = TRUE),
    coxph_default = survival::coxph(survival::Surv(time, event) ~ x),
    times = times
  )
  # include_cld = FALSE: no comparison letters (which would need 'multcomp');
  # lq and uq are the hinges of fivenum()
  s <- summary(mb, unit = "us", include_cld = FALSE)
  timing_rows[[j]] <- data.frame(n = n, events = sum(event),
                                 event_times = length(rs$time),
                                 expr = as.character(s$expr),
                                 median_us = s$median, lq_us = s$lq,
                                 uq_us = s$uq, neval = s$neval,
                                 stringsAsFactors = FALSE)
  iter_rows[[j]] <- data.frame(
    n = n,
    iter_1e8 = onestepcox:::newton_cox(eta0, dt, tol = 1e-8)$iter,
    iter_1e12 = onestepcox:::newton_cox(eta0, dt, tol = 1e-12)$iter
  )
}
timing <- list(times = do.call(rbind, timing_rows),
               iterations = do.call(rbind, iter_rows),
               machine = Sys.info()[c("sysname", "release", "machine")])
saveRDS(timing, file.path(out, "data", "timing.rds"))
run_log[[length(run_log) + 1L]] <-
  log_entry("timing", length(n_time), 0L, t0)

write_run_log(run_log, out, "run_log_supplement")
message("Finished. Results in ", file.path(out, "data"))
