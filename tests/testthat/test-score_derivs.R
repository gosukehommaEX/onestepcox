# Tests for the internal function score_derivs()

# Log partial likelihood of the two-group model, written per event time from
# RiskTable(), independently of deriv_table()
loglik <- function(eta, rt, ties) {
  ll <- 0
  for (k in seq_len(nrow(rt))) {
    O_C <- rt$O[k] - rt$O_T[k]
    ll <- ll + rt$O_T[k] * eta
    if (ties == "breslow") {
      ll <- ll - rt$O[k] * log(rt$n_C[k] + rt$n_T[k] * exp(eta))
    } else {
      r <- seq_len(rt$O[k]) - 1
      ll <- ll - sum(log(rt$n_C[k] - r / rt$O[k] * O_C +
                           (rt$n_T[k] - r / rt$O[k] * rt$O_T[k]) * exp(eta)))
    }
  }
  ll
}

derivs_at <- function(eta, d, ties) {
  rs <- risk_sets(d$time, d$event, d$group != d$control, d$strata)
  score_derivs(eta, deriv_table(rs, ties))
}

# Compare U, I, J, and L with the Python values for both handlings of ties
check_derivs <- function(name) {
  d <- test_data(name)
  for (ti in c("efron", "breslow")) {
    ref <- expected_derivs[expected_derivs$data == name &
                             expected_derivs$ties == ti, ]
    expect_equal(nrow(ref), 1L, info = paste(name, ti))
    got <- derivs_at(ref$eta, d, ti)
    expect_equal(unname(got), c(ref$U, ref$I, ref$J, ref$L),
                 tolerance = 1e-10, info = paste(name, ti))
    expect_named(got, c("U", "I", "J", "L"))
  }
}

test_that("Derivatives match the Python values (survival data)", {
  for (nm in c("ovarian", "gbsg", "gbsg_grade")) {
    check_derivs(nm)
  }
})

test_that("Derivatives match the Python values (pharmacoSmoking)", {
  skip_if_not_installed("asaur")
  for (nm in c("pharmacoSmoking", "pharmacoSmoking_ageGroup2")) {
    check_derivs(nm)
  }
})

test_that("Derivatives match numerical derivatives of the log likelihood", {
  d <- test_data("gbsg_grade")
  h <- 1e-4
  for (ti in c("efron", "breslow")) {
    rt <- RiskTable(d$time, d$event, d$group, d$control, strata = d$strata)
    for (eta in c(-0.8, 0.4)) {
      dd <- derivs_at(eta, d, ti)
      up <- derivs_at(eta + h, d, ti)
      dn <- derivs_at(eta - h, d, ti)
      U_num <- (loglik(eta + h, rt, ti) - loglik(eta - h, rt, ti)) / (2 * h)
      expect_equal(dd[["U"]], U_num, tolerance = 1e-6, info = ti)
      expect_equal(dd[["I"]], -(up[["U"]] - dn[["U"]]) / (2 * h),
                   tolerance = 1e-6, info = ti)
      expect_equal(dd[["J"]], (up[["I"]] - dn[["I"]]) / (2 * h),
                   tolerance = 1e-6, info = ti)
      expect_equal(dd[["L"]], (up[["J"]] - dn[["J"]]) / (2 * h),
                   tolerance = 1e-5, info = ti)
    }
  }
})
