# Tests for the internal function summarize_oc()

# Estimates written by hand. The method has no estimate in the third
# replicate and the reference none in the fifth. The expected values were
# computed beforehand by an independent implementation in numpy,
# tools/test_expected/summarize_oc_expected.py
oc_inputs <- function() {
  x <- data.frame(coef = log(c(0.50, 0.70, NA, 1.25, 0.92, 2.10)),
                  se = c(0.20, 0.25, NA, 0.30, 0.22, 0.35),
                  lower = c(0.35, 0.45, NA, 0.70, 0.60, 1.20),
                  upper = c(0.72, 1.10, NA, 2.20, 1.40, 3.60))
  ref <- data.frame(coef = c(log(0.503), log(0.70) + 1e-4, log(0.80),
                             log(1.25) - 0.02, NA, log(2.104)),
                    se = c(0.21, 0.24, 0.26, 0.31, NA, 0.33),
                    lower = c(0.33, 0.44, 0.48, 0.67, NA, 1.11),
                    upper = c(0.77, 1.12, 1.33, 2.25, NA, 4.01))
  list(x = x, ref = ref)
}

test_that("The summaries agree with an independent implementation", {
  d <- oc_inputs()
  s <- summarize_oc(d$x, d$ref, HR_true = 0.8, time_per_sim = 0.0012)
  expected <- c(
    n_NA = 1,
    dev_median = 0.0039425109118167989,
    dev_q90 = 0.015794621503264245,
    dev_q99 = 0.019579462150326413,
    dev_max = 0.01999999999999999,
    devse_median = 0.01712628560078467,
    devse_q99 = 0.06343522682951136,
    devse_max = 0.064516129032258035,
    mismatch_2dp = 0.25,
    se_ratio_median = 1.0047043010752688,
    se_ratio_min = 0.95238095238095244,
    se_ratio_max = 1.0606060606060606,
    bias = 0.18951898383538138,
    sd = 0.549825775261507,
    rmse = 0.52703334974888394,
    se_mean = 0.26399999999999996,
    coverage = 0.59999999999999998,
    time_per_sim = 0.0011999999999999999
  )
  expect_identical(names(s), names(expected))
  expect_equal(unlist(s), expected, tolerance = 1e-12)
  # The RMSE, bias, and SD (divisor n - 1) are related as expected
  n <- 5
  expect_equal(s$rmse^2, s$bias^2 + s$sd^2 * (n - 1) / n, tolerance = 1e-12)
})

test_that("The reference against itself has no deviation or mismatch", {
  d <- oc_inputs()
  s <- summarize_oc(d$ref, d$ref, HR_true = 0.8, time_per_sim = 0)
  expect_equal(c(s$dev_max, s$devse_max, s$mismatch_2dp), c(0, 0, 0))
  expect_equal(c(s$se_ratio_min, s$se_ratio_max), c(1, 1))
  expect_equal(s$n_NA, 1)
})

test_that("A method without any estimate gives NA summaries", {
  d <- oc_inputs()
  none <- d$x
  none[] <- NA_real_
  s <- summarize_oc(none, d$ref, HR_true = 0.8, time_per_sim = 0.5)
  expect_equal(s$n_NA, 6)
  expect_equal(s$time_per_sim, 0.5)
  expect_true(all(is.na(unlist(s[setdiff(names(s), c("n_NA",
                                                     "time_per_sim"))]))))
  # One estimate: no standard deviation, but the other summaries
  one <- none
  one[2L, ] <- d$x[2L, ]
  s <- summarize_oc(one, d$ref, HR_true = 0.8, time_per_sim = 0.5)
  expect_true(is.na(s$sd))
  expect_equal(s$bias, log(0.70) - log(0.8))
  expect_equal(s$mismatch_2dp, 0)
})
