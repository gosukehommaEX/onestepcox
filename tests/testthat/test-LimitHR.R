# Tests for LimitHR()

limit_args <- function(design) {
  d <- expected_limit_designs[expected_limit_designs$design == design, ]
  list(hazard_C = d$hazard_C, n_T = d$n_T, n_C = d$n_C, a = d$a, E = d$E,
       dropout_hazard = d$dropout_hazard)
}

test_that("Limits and analysis times match the Python values", {
  # Design A has its analysis after the end of accrual, design B before it
  for (des in c("A", "B")) {
    ref <- expected_limits[expected_limits$design == des, ]
    got <- do.call(LimitHR, c(list(HR = unique(ref$HR)), limit_args(des)))
    expect_equal(nrow(got), nrow(ref), info = des)
    pos <- match(paste(ref$HR, ref$method), paste(got$HR, got$method))
    expect_false(anyNA(pos), info = des)
    expect_lt(max(abs(got$limit[pos] - ref$limit)), 1e-9)
    expect_lt(max(abs(got$analysis_time[pos] / ref$analysis_time - 1)), 1e-9)
    expect_equal(got$bias, got$limit - log(got$HR), tolerance = 1e-14)
  }
})

test_that("The MH-anchored and Cox limits are the true log hazard ratio", {
  lim <- LimitHR(HR = c(0.3, 1, 4), hazard_C = 0.05, n_T = 150, n_C = 100,
                 a = 10, E = 120, dropout_hazard = 0.01)
  for (m in c("MH", "MH_Newton", "MH_Halley", "MH_Chebyshev", "Cox")) {
    expect_lt(max(abs(lim$bias[lim$method == m])), 1e-9)
  }
  # The Pike anchor is biased away from HR = 1
  expect_gt(min(abs(lim$bias[lim$method == "Pike" & lim$HR != 1])), 1e-3)
})

test_that("All limits are zero when HR = 1", {
  lim <- LimitHR(1, hazard_C = 0.05, n_T = 150, n_C = 100, a = 10, E = 120)
  expect_lt(max(abs(lim$limit)), 1e-10)
})

test_that("The one-step limits from the Pike anchor are ordered", {
  # Newton removes most of the Pike bias, and Halley and Chebyshev more
  lim <- LimitHR(c(0.2, 5), hazard_C = log(2) / 18, n_T = 100, n_C = 100,
                 a = 12, E = 160)
  for (hr in c(0.2, 5)) {
    b <- abs(stats::setNames(lim$bias, lim$method)[lim$HR == hr])
    expect_lt(b[["Pike_Newton"]], b[["Pike"]])
    expect_lt(b[["Pike_Chebyshev"]], b[["Pike_Newton"]])
    expect_lt(b[["Pike_Halley"]], b[["Pike_Newton"]])
  }
})

test_that("Output has one row per hazard ratio and method", {
  lim <- LimitHR(c(0.5, 2, 3), hazard_C = 0.05, n_T = 100, n_C = 100,
                 a = 10, E = 120)
  expect_identical(class(lim), c("LimitHR", "data.frame"))
  expect_named(lim, c("HR", "analysis_time", "method", "limit", "bias"))
  expect_identical(attr(lim, "design"),
                   c(hazard_C = 0.05, n_T = 100, n_C = 100, a = 10, E = 120,
                     dropout_hazard = 0))
  expect_identical(class(as.data.frame(lim)), "data.frame")
  expect_equal(nrow(lim), 3 * 10)
  expect_identical(lim$method[1:10],
                   c("MH_Halley", "MH_Newton", "Pike_Halley", "Pike_Newton",
                     "MH_Chebyshev", "Pike_Chebyshev", "MH", "Pike", "Peto",
                     "Cox"))
  expect_equal(lim$HR, rep(c(0.5, 2, 3), each = 10))
})

test_that("Invalid arguments give errors", {
  expect_error(LimitHR(c(1, -1), 0.05, 100, 100, 10, 120), "HR")
  # Without dropout every subject has an event eventually, so E must be
  # below n_T + n_C
  expect_error(LimitHR(2, 0.05, 100, 100, 10, 200), "cannot be reached")
  # With dropout only part of the subjects has an event
  expect_error(LimitHR(2, 0.05, 100, 100, 10, 150, dropout_hazard = 0.05),
               "cannot be reached")
})
