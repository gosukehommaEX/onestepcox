# Tests for EstHR()

# Compare all methods except "Cox" with the values of the independent
# Python implementation (helper-expected.R)
check_against_python <- function(name) {
  d <- test_data(name)
  ref <- expected_estimates[expected_estimates$data == name, ]
  # 10 methods (all but "Cox") times 2 handlings of ties
  expect_equal(nrow(ref), 20L, info = name)
  for (i in seq_len(nrow(ref))) {
    est <- EstHR(d$time, d$event, d$group, d$control, ref$method[i],
                 ties = ref$ties[i], strata = d$strata)
    lab <- paste(name, ref$ties[i], ref$method[i])
    expect_equal(est[["coef"]], ref$coef[i], tolerance = 1e-10,
                 info = paste(lab, "coef"))
    expect_equal(est[["se(coef)"]], ref$se[i], tolerance = 1e-10,
                 info = paste(lab, "se"))
  }
}

# ---- Output -------------------------------------------------------------

test_that("Output has the names of summary.coxph() and a Wald interval", {
  d <- test_data("ovarian")
  est <- EstHR(d$time, d$event, d$group, d$control, "MH_Halley")
  expect_named(est, c("coef", "exp(coef)", "se(coef)",
                      "lower .95", "upper .95"))
  z <- stats::qnorm(0.975)
  b <- est[["coef"]]
  se <- est[["se(coef)"]]
  expect_equal(est[["exp(coef)"]], exp(b), tolerance = 1e-14)
  expect_equal(est[["lower .95"]], exp(b - z * se), tolerance = 1e-14)
  expect_equal(est[["upper .95"]], exp(b + z * se), tolerance = 1e-14)

  est90 <- EstHR(d$time, d$event, d$group, d$control, "MH_Halley",
                 conf.level = 0.9)
  expect_named(est90, c("coef", "exp(coef)", "se(coef)",
                        "lower .90", "upper .90"))
  df <- data.frame(time = d$time, event = d$event, x = d$group)
  fit <- survival::coxph(Surv(time, event) ~ x, data = df)
  expect_identical(colnames(summary(fit, conf.int = 0.9)$conf.int)[3:4],
                   names(est90)[4:5])
})

# ---- Agreement with the independent Python implementation ---------------

test_that("Estimates and SEs match the Python values (survival data)", {
  for (nm in c("ovarian", "gbsg", "gbsg_grade")) {
    check_against_python(nm)
  }
})

test_that("Estimates and SEs match the Python values (pharmacoSmoking)", {
  skip_if_not_installed("asaur")
  for (nm in c("pharmacoSmoking", "pharmacoSmoking_ageGroup2")) {
    check_against_python(nm)
  }
})

# ---- Cox ---------------------------------------------------------------

test_that("Cox matches survival::coxph() with its default settings", {
  for (nm in c("ovarian", "gbsg", "gbsg_grade")) {
    d <- test_data(nm)
    df <- data.frame(time = d$time, event = d$event,
                     x = as.numeric(d$group != d$control))
    for (ti in c("efron", "breslow")) {
      est <- EstHR(d$time, d$event, d$group, d$control, "Cox",
                   ties = ti, strata = d$strata)
      if (is.null(d$strata)) {
        fit <- survival::coxph(Surv(time, event) ~ x, data = df, ties = ti)
      } else {
        df$s <- d$strata
        fit <- survival::coxph(Surv(time, event) ~ x + strata(s),
                               data = df, ties = ti)
      }
      # coxph() stops at eps = 1e-9 by default, so the tolerance is loose;
      # the next test checks the accuracy of "Cox"
      expect_equal(est[["coef"]], unname(stats::coef(fit)),
                   tolerance = 1e-7, info = paste(nm, ti))
      expect_equal(est[["se(coef)"]], unname(sqrt(fit$var[1, 1])),
                   tolerance = 1e-7, info = paste(nm, ti))
    }
  }
})

test_that("Cox and Cox_NR agree with the Python root of the score", {
  for (nm in c("ovarian", "gbsg", "gbsg_grade")) {
    d <- test_data(nm)
    for (ti in c("efron", "breslow")) {
      ref <- expected_estimates[expected_estimates$data == nm &
                                  expected_estimates$ties == ti &
                                  expected_estimates$method == "Cox_NR", ]
      for (m in c("Cox", "Cox_NR")) {
        est <- EstHR(d$time, d$event, d$group, d$control, m,
                     ties = ti, strata = d$strata)
        expect_equal(est[["coef"]], ref$coef, tolerance = 1e-9,
                     info = paste(nm, ti, m))
        expect_equal(est[["se(coef)"]], ref$se, tolerance = 1e-9,
                     info = paste(nm, ti, m))
      }
    }
  }
})

# ---- Hand calculations ---------------------------------------------------

test_that("MH, Pike, and Peto reproduce hand calculations", {
  # No ties. Sorted by time, the risk sets are
  # k:    n_T  n_C  O_T  O_C
  # 1:     3    3    0    1
  # 2:     3    2    1    0
  # 3:     2    2    0    1
  # 4:     2    1    1    0
  # 5:     1    1    0    1
  # 6:     1    0    1    0
  time  <- c(1, 2, 3, 4, 5, 6)
  event <- c(1, 1, 1, 1, 1, 1)
  group <- c("C", "T", "C", "T", "C", "T")
  n_T <- c(3, 3, 2, 2, 1, 1)
  n_C <- c(3, 2, 2, 1, 1, 0)
  O_T <- c(0, 1, 0, 1, 0, 1)
  O_C <- 1 - O_T
  n <- n_T + n_C
  E_T <- sum(n_T / n)
  E_C <- sum(n_C / n)

  R <- sum(O_T * n_C / n)
  S <- sum(O_C * n_T / n)
  mh <- EstHR(time, event, group, "C", "MH")
  expect_equal(mh[["coef"]], log(R / S), tolerance = 1e-12)
  expect_equal(mh[["se(coef)"]], sqrt(sum(n_T * n_C / n^2) / (R * S)),
               tolerance = 1e-12)

  pike <- EstHR(time, event, group, "C", "Pike")
  expect_equal(pike[["coef"]], log(3 * E_C / (3 * E_T)), tolerance = 1e-12)
  expect_equal(pike[["se(coef)"]], sqrt(1 / E_T + 1 / E_C), tolerance = 1e-12)

  # Without ties, V = sum n_T n_C / n^2 over the times with n > 1
  ok <- n > 1
  V <- sum(n_T[ok] * n_C[ok] / n[ok]^2)
  peto <- EstHR(time, event, group, "C", "Peto")
  expect_equal(peto[["coef"]], (3 - E_T) / V, tolerance = 1e-12)
  expect_equal(peto[["se(coef)"]], 1 / sqrt(V), tolerance = 1e-12)
})

# ---- Invariances ---------------------------------------------------------

test_that("Swapping the groups negates coef and keeps se", {
  d <- test_data("gbsg_grade")
  for (ti in c("efron", "breslow")) {
    for (m in all_methods) {
      a <- EstHR(d$time, d$event, d$group, 0, m, ties = ti, strata = d$strata)
      b <- EstHR(d$time, d$event, d$group, 1, m, ties = ti, strata = d$strata)
      tol <- if (m == "Cox") 1e-8 else 1e-10
      expect_equal(b[["coef"]], -a[["coef"]], tolerance = tol,
                   info = paste(ti, m))
      expect_equal(b[["se(coef)"]], a[["se(coef)"]], tolerance = tol,
                   info = paste(ti, m))
    }
  }
})

test_that("Row order does not matter and is_sorted = TRUE skips sorting", {
  d <- test_data("gbsg_grade")
  set.seed(3)
  p <- sample(length(d$time))
  o <- order(d$strata, d$time)
  for (m in setdiff(all_methods, "Cox")) {
    a <- EstHR(d$time, d$event, d$group, d$control, m, strata = d$strata)
    b <- EstHR(d$time[p], d$event[p], d$group[p], d$control, m,
               strata = d$strata[p])
    s <- EstHR(d$time[o], d$event[o], d$group[o], d$control, m,
               strata = d$strata[o], is_sorted = TRUE)
    expect_equal(b, a, tolerance = 1e-12, info = m)
    expect_equal(s, a, tolerance = 1e-12, info = m)
  }
  expect_error(EstHR(d$time[p], d$event[p], d$group[p], d$control,
                     "MH_Halley", strata = d$strata[p], is_sorted = TRUE),
               "is_sorted")
  expect_error(EstHR(d$time[p], d$event[p], d$group[p], d$control,
                     "MH_Halley", is_sorted = TRUE),
               "is_sorted")
})

test_that("Relabelling the strata does not change the estimates", {
  d <- test_data("gbsg_grade")
  relabelled <- c("z", "x", "y")[d$strata]
  for (m in all_methods) {
    a <- EstHR(d$time, d$event, d$group, d$control, m, strata = d$strata)
    b <- EstHR(d$time, d$event, d$group, d$control, m, strata = relabelled)
    tol <- if (m == "Cox") 1e-8 else 1e-12
    expect_equal(b, a, tolerance = tol, info = m)
  }
})

test_that("A single stratum gives the unstratified estimates", {
  d <- test_data("gbsg")
  one <- rep(1, length(d$time))
  for (m in all_methods) {
    a <- EstHR(d$time, d$event, d$group, d$control, m)
    b <- EstHR(d$time, d$event, d$group, d$control, m, strata = one)
    tol <- if (m == "Cox") 1e-8 else 1e-10
    expect_equal(b, a, tolerance = tol, info = m)
  }
})

test_that("Efron and Breslow differ only for methods using the score", {
  d <- test_data("gbsg")   # tied event times
  for (m in all_methods) {
    e <- EstHR(d$time, d$event, d$group, d$control, m, ties = "efron")
    b <- EstHR(d$time, d$event, d$group, d$control, m, ties = "breslow")
    if (m %in% c("MH", "Pike", "Peto")) {
      expect_identical(e, b, info = m)
    } else {
      expect_gt(abs(e[["coef"]] - b[["coef"]]), 1e-6)
    }
  }
  d <- test_data("ovarian")   # no tied event times
  for (m in all_methods) {
    e <- EstHR(d$time, d$event, d$group, d$control, m, ties = "efron")
    b <- EstHR(d$time, d$event, d$group, d$control, m, ties = "breslow")
    tol <- if (m == "Cox") 1e-8 else 1e-10
    expect_equal(e, b, tolerance = tol, info = m)
  }
})

# ---- Accuracy -----------------------------------------------------------

test_that("MH_Halley is much closer to the Cox estimate than Pike_Chebyshev", {
  # HR = 4, where the Pike anchor is biased; the largest ratio of the mean
  # deviations over 200 sets of 20 replicates in the Python implementation
  # was 0.04
  dat <- GenSimData(nsim = 20, n_T = 100, n_C = 100,
                    hazard_C = log(2) / 18, HR = 4, a = 12, E = 160,
                    seed = 1)
  dev <- dat[, list(
    mhh = abs(EstHR(time, event, group, "C", "MH_Halley")[["coef"]] -
                EstHR(time, event, group, "C", "Cox_NR")[["coef"]]),
    pkc = abs(EstHR(time, event, group, "C", "Pike_Chebyshev")[["coef"]] -
                EstHR(time, event, group, "C", "Cox_NR")[["coef"]])
  ), by = simID]
  expect_lt(mean(dev$mhh), mean(dev$pkc) / 5)
})

# ---- Cases without an estimate -----------------------------------------

test_that("NA is returned when the estimate cannot be computed", {
  is_na5 <- function(x) length(x) == 5L && all(is.na(x))
  for (m in all_methods) {
    expect_true(is_na5(EstHR(numeric(0), integer(0), character(0),
                             "C", m)), info = m)
    expect_true(is_na5(EstHR(c(1, 2, 3, 4), c(0, 0, 0, 0),
                             c("C", "T", "C", "T"), "C", m)), info = m)
  }
  # Events in the control group only: the Cox estimate is infinite
  for (m in setdiff(all_methods, "Cox")) {
    expect_true(is_na5(EstHR(c(1, 2, 3, 4), c(1, 0, 1, 0),
                             c("C", "T", "C", "T"), "C", m)), info = m)
  }
})

# ---- Input checks ---------------------------------------------------------

test_that("Invalid arguments give informative errors", {
  d <- test_data("ovarian")
  expect_error(EstHR(d$time, d$event, d$group, d$control, "PiHE"), "method")
  expect_error(EstHR(d$time, d$event, d$group, d$control, "MH_Halley",
                     ties = "exact"), "ties")
  expect_error(EstHR(d$time, d$event, d$group, d$control, "MH_Halley",
                     conf.level = 95), "conf.level")
  expect_error(EstHR(d$time, d$event, d$group, d$control, "MH_Halley",
                     tol = 0), "tol")
  expect_error(EstHR(c(1, 2, 3), c(1, 0), c("C", "T", "C"), "C", "MH"),
               "same length")
  expect_error(EstHR(d$time, d$event, d$group, 3, "MH_Halley"), "control")
})
