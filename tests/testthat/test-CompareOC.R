# Tests for CompareOC()

scenarios_demo <- function(extra = NULL) {
  base <- data.frame(
    label    = c("S1", "S2"),
    n_T      = c(80, 100),
    n_C      = c(80, 100),
    hazard_C = c(log(2) / 18, log(2) / 12),
    HR       = c(0.7, 0.8),
    a        = c(12, 12),
    E        = c(60, 90),
    stringsAsFactors = FALSE
  )
  if (!is.null(extra)) base <- cbind(base, extra)
  base
}

test_that("Output has the documented structure", {
  res <- CompareOC(
    scenarios = scenarios_demo(),
    methods   = c("Cox", "Peto", "Pike", "MH_Halley"),
    nsim      = 30,
    seed      = 1
  )
  expect_s3_class(res, "CompareOC")
  expect_named(res, c("estimates", "summary"))
  expect_identical(attr(res, "reference"), "Cox")
  expect_identical(attr(res, "ties"), "efron")
  expect_identical(attr(res, "conf.level"), 0.95)
  expect_identical(attr(res, "nsim"), 30)
  expect_named(res$estimates,
               c("label", "simID", "method", "coef", "se", "lower", "upper"))
  expect_named(res$summary,
               c("label", "HR_true", "method", "n_NA", "dev_median",
                 "dev_q90", "dev_q99", "dev_max", "se_ratio_median",
                 "se_ratio_min", "se_ratio_max", "bias", "rmse", "coverage",
                 "time_per_sim"))
  # 2 scenarios * 4 methods (* 30 replicates)
  expect_equal(nrow(res$summary), 2 * 4)
  expect_equal(nrow(res$estimates), 2 * 4 * 30)
  expect_true(all(res$summary$n_NA == 0))
})

test_that("The reference has zero deviation and unit SE ratio", {
  res <- CompareOC(scenarios_demo(), c("Peto", "Cox", "MH_Halley"),
                   nsim = 30, seed = 1, reference = "Cox")
  ref_rows <- res$summary[res$summary$method == "Cox", ]
  expect_true(all(ref_rows$dev_max == 0))
  expect_true(all(ref_rows$se_ratio_min == 1 & ref_rows$se_ratio_max == 1))
  # The reference is computed first
  expect_identical(res$summary$method[1:3], c("Cox", "Peto", "MH_Halley"))
})

test_that("Summaries agree with the estimates", {
  res <- CompareOC(scenarios_demo(), c("Cox", "MH_Halley", "Pike"),
                   nsim = 30, seed = 2)
  for (lab in c("S1", "S2")) {
    r <- res$estimates[res$estimates$label == lab &
                         res$estimates$method == "Cox", ]
    for (m in c("MH_Halley", "Pike")) {
      e <- res$estimates[res$estimates$label == lab &
                           res$estimates$method == m, ]
      s <- res$summary[res$summary$label == lab & res$summary$method == m, ]
      eta <- log(s$HR_true)
      dev <- abs(e$coef - r$coef)
      info <- paste(lab, m)
      expect_equal(s$bias, mean(e$coef - eta), info = info)
      expect_equal(s$rmse, sqrt(mean((e$coef - eta)^2)), info = info)
      expect_equal(s$coverage,
                   mean(e$lower <= s$HR_true & s$HR_true <= e$upper),
                   info = info)
      expect_equal(s$dev_median, stats::median(dev), info = info)
      expect_equal(s$dev_max, max(dev), info = info)
      expect_equal(s$dev_q99, unname(stats::quantile(dev, 0.99)), info = info)
      expect_equal(s$se_ratio_median, stats::median(e$se / r$se), info = info)
    }
  }
})

test_that("Estimates equal EstHR() on the generated data with strata and ties", {
  sc <- scenarios_demo()[1, ]
  sc$strata_prop <- list(c(0.4, 0.6))
  sc$strata_hazard <- list(c(1, 2))
  sc$time_unit <- 1
  res <- CompareOC(sc, c("Cox", "MH_Halley"), nsim = 5, seed = 3,
                   ties = "breslow", conf.level = 0.9)
  dat <- GenSimData(nsim = 5, n_T = 80, n_C = 80, hazard_C = log(2) / 18,
                    HR = 0.7, a = 12, E = 60, strata_prop = c(0.4, 0.6),
                    strata_hazard = c(1, 2), time_unit = 1, seed = 3)
  d2 <- dat[simID == 2L]
  for (m in c("Cox", "MH_Halley")) {
    e <- EstHR(d2$time, d2$event, d2$group, "C", m, ties = "breslow",
               strata = d2$stratum, conf.level = 0.9)
    got <- res$estimates[res$estimates$simID == 2L &
                           res$estimates$method == m, ]
    expect_equal(c(got$coef, got$se, got$lower, got$upper), unname(e[-2]),
                 tolerance = 1e-12, info = m)
  }
})

test_that("MH_Halley deviates much less from Cox than Pike", {
  # In 400 runs of these scenarios in the Python implementation, the 90th
  # percentile for MH_Halley was at most 1e-4 times the median for Pike
  res <- CompareOC(scenarios_demo(), c("Cox", "Pike", "MH_Halley"),
                   nsim = 50, seed = 11)
  for (lab in c("S1", "S2")) {
    s <- res$summary[res$summary$label == lab, ]
    expect_lt(s$dev_q90[s$method == "MH_Halley"],
              s$dev_median[s$method == "Pike"] / 100)
  }
})

test_that("Coverage, bias, and RMSE are in range", {
  res <- CompareOC(scenarios_demo(), c("Cox", "Peto", "MH_Halley"),
                   nsim = 30, seed = 4)
  expect_true(all(res$summary$coverage >= 0 & res$summary$coverage <= 1))
  expect_true(all(res$summary$rmse >= abs(res$summary$bias)))
})

test_that("Optional columns default to no dropout and no rounding", {
  res_no <- CompareOC(scenarios_demo(), c("Cox", "MH_Halley"),
                      nsim = 20, seed = 7)
  res_with <- CompareOC(scenarios_demo(extra = data.frame(
    dropout_hazard = c(0, 0), time_unit = c(0, 0))),
    c("Cox", "MH_Halley"), nsim = 20, seed = 7)
  expect_equal(res_no$estimates, res_with$estimates, tolerance = 1e-12)
})

test_that("Reference must be one of the supplied methods", {
  expect_error(
    CompareOC(scenarios_demo(), methods = c("Peto", "Pike", "MH_Halley"),
              nsim = 10, seed = 1, reference = "Cox"),
    "reference"
  )
})

test_that("Missing required scenario columns triggers an error", {
  bad <- scenarios_demo()
  bad$E <- NULL
  expect_error(
    CompareOC(scenarios = bad, methods = c("Cox", "MH_Halley"),
              nsim = 10, seed = 1),
    "missing required columns"
  )
})
