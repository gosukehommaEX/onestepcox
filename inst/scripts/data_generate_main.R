# Data of the figures and tables of the main text of the article.
#
# Run from the root of the package (open onestepcox.Rproj) after
# devtools::load_all(), or with the package installed:
#   source("inst/scripts/data_generate_main.R")
# For a quick test run with 1/100 of the replicates, written to
# inst/scripts/results/quick/, first run options(onestepcox.quick = TRUE).
# The sections, seeds, and expected computing times are described in
# inst/scripts/README.md. This script only computes and saves the results;
# table_and_figure_manuscript.R makes the tables and figures from them.

source(file.path("inst", "scripts", "settings.R"))
# Checked first, so that a missing package does not stop the run at the end
if (!requireNamespace("asaur", quietly = TRUE)) {
  stop("The package 'asaur' is needed for the pharmacoSmoking data")
}
out <- results_dir("scripts")
run_log <- list()

# 1. Large-sample limits under the design of GenSimData() with 1:1
#    allocation (Figure 1), and at the hazard ratios of the growth scenarios
#    (the limits of the Pike-anchored estimators in Figures 2 and 3)
t0 <- proc.time()[["elapsed"]]
message("[limit] large-sample limits")
hr_grid <- exp(seq(log(0.1), log(10), length.out = if (quick) 9L else 81L))
limit <- list(
  grid = onestepcox::LimitHR(HR = hr_grid, hazard_C = hazard_C, n_T = 1,
                             n_C = 1, a = accrual, E = 2 * event_frac),
  fixed = onestepcox::LimitHR(HR = theta_fixed, hazard_C = hazard_C,
                              n_T = 1, n_C = 1, a = accrual,
                              E = 2 * event_frac)
)
saveRDS(limit, file.path(out, "data", "limit.rds"))
run_log[[length(run_log) + 1L]] <-
  log_entry("limit", length(hr_grid) + length(theta_fixed), 0L, t0)

# 2. Fixed alternatives with n from 100 to 100,000 (Figures 2 to 4)
t0 <- proc.time()[["elapsed"]]
sc <- growth_scenarios("fixed")
fixed <- run_growth(sc, seeds[["fixed"]], "fixed")
saveRDS(fixed, file.path(out, "data", "fixed.rds"))
run_log[[length(run_log) + 1L]] <-
  log_entry("fixed", nrow(sc), sum(vapply(sc$nsim, n_rep, integer(1))), t0)

# 3. Local alternatives exp(c / sqrt(n)) with the same n (Figures 2 to 4)
t0 <- proc.time()[["elapsed"]]
sc <- growth_scenarios("local")
local_alt <- run_growth(sc, seeds[["local"]], "local")
saveRDS(local_alt, file.path(out, "data", "local.rds"))
run_log[[length(run_log) + 1L]] <-
  log_entry("local", nrow(sc), sum(vapply(sc$nsim, n_rep, integer(1))), t0)

# 4. Grid A: operating characteristics with continuous times (Table 1)
t0 <- proc.time()[["elapsed"]]
sc <- gridA_scenarios()
gridA <- list(scenarios = sc,
              summary = run_grid(sc, seeds[["gridA"]], "gridA"))
saveRDS(gridA, file.path(out, "data", "gridA.rds"))
run_log[[length(run_log) + 1L]] <-
  log_entry("gridA", nrow(sc), nrow(sc) * n_rep(10000L), t0)

# 5. Grid B: dropout, ties, and strata (Table 1)
t0 <- proc.time()[["elapsed"]]
sc <- gridB_scenarios()
gridB <- list(scenarios = sc,
              summary = run_grid(sc, seeds[["gridB"]], "gridB"))
saveRDS(gridB, file.path(out, "data", "gridB.rds"))
run_log[[length(run_log) + 1L]] <-
  log_entry("gridB", nrow(sc), nrow(sc) * n_rep(10000L), t0)

# 6. Applications (Table 2): gbsg of 'survival' (hormonal therapy, control
#    hormon = 0) and pharmacoSmoking of 'asaur' (control "patchOnly", as in
#    the earlier version of the article), without and with strata, with the
#    Efron and Breslow handling of ties
t0 <- proc.time()[["elapsed"]]
message("[application] gbsg and pharmacoSmoking")
analyses <- list(
  gbsg = list(formula = Surv(rfstime, status) ~ hormon,
              data = survival::gbsg, control = 0),
  gbsg_grade = list(formula = Surv(rfstime, status) ~ hormon + strata(grade),
                    data = survival::gbsg, control = 0),
  pharmacoSmoking = list(formula = Surv(ttr, relapse) ~ grp,
                         data = asaur::pharmacoSmoking,
                         control = "patchOnly"),
  pharmacoSmoking_ageGroup2 = list(
    formula = Surv(ttr, relapse) ~ grp + strata(ageGroup2),
    data = asaur::pharmacoSmoking, control = "patchOnly"
  )
)
app_rows <- list()
for (a in names(analyses)) {
  an <- analyses[[a]]
  rt <- onestepcox::RiskTable(an$formula, data = an$data,
                              control = an$control)
  for (ties in c("efron", "breslow")) {
    est <- onestepcox::EstHR(an$formula, data = an$data,
                             control = an$control, method = methods,
                             ties = ties)
    app_rows[[length(app_rows) + 1L]] <- data.frame(
      analysis = a, ties = ties, method = methods,
      coef = unname(est[, 1L]), HR = unname(est[, 2L]),
      se = unname(est[, 3L]), lower = unname(est[, 4L]),
      upper = unname(est[, 5L]),
      treatment = attr(est, "groups")[["treatment"]],
      control = attr(est, "groups")[["control"]],
      n_T = attr(est, "n")[["treatment"]],
      n_C = attr(est, "n")[["control"]],
      events_T = attr(est, "events")[["treatment"]],
      events_C = attr(est, "events")[["control"]],
      n_strata = attr(est, "n_strata"), n_removed = attr(est, "n_removed"),
      event_times = nrow(rt), tied_events = sum(rt$O[rt$O > 1]),
      stringsAsFactors = FALSE
    )
  }
}
application <- do.call(rbind, app_rows)
rownames(application) <- NULL
saveRDS(application, file.path(out, "data", "application.rds"))
# The same estimates with 17 significant digits, for the check with an
# independent implementation outside R
app_csv <- application
for (v in c("coef", "HR", "se", "lower", "upper")) {
  app_csv[[v]] <- sprintf("%.17g", app_csv[[v]])
}
utils::write.csv(app_csv, file.path(out, "data", "application.csv"),
                 row.names = FALSE)
run_log[[length(run_log) + 1L]] <-
  log_entry("application", length(analyses), 0L, t0)

write_run_log(run_log, out, "run_log_main")
message("Finished. Results in ", file.path(out, "data"))
