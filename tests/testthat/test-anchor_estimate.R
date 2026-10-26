# Tests for the internal function anchor_estimate()

# Two event times with ties: n_T, n_C, O_T, O
rs_hand <- list(n_T = c(10, 6), n_C = c(10, 8), O_T = c(2, 1), O = c(3, 3))

test_that("MH reproduces a hand calculation", {
  n <- c(20, 14)
  R <- 2 * 10 / 20 + 1 * 8 / 14
  S <- 1 * 10 / 20 + 2 * 6 / 14
  v <- (3 * 10 * 10 / 20^2 + 3 * 6 * 8 / 14^2) / (R * S)
  a <- anchor_estimate(rs_hand, "MH")
  expect_equal(a, c(eta = log(R / S), se = sqrt(v)), tolerance = 1e-14)
  expect_equal(n, rs_hand$n_T + rs_hand$n_C)
})

test_that("Pike reproduces a hand calculation", {
  E_T <- 10 * 3 / 20 + 6 * 3 / 14
  E_C <- 10 * 3 / 20 + 8 * 3 / 14
  a <- anchor_estimate(rs_hand, "Pike")
  expect_equal(a, c(eta = log(3 * E_C / (3 * E_T)),
                    se = sqrt(1 / E_T + 1 / E_C)), tolerance = 1e-14)
})

test_that("Peto uses the hypergeometric variance with ties", {
  E_T <- 10 * 3 / 20 + 6 * 3 / 14
  V <- 10 * 10 * 3 * 17 / (20^2 * 19) + 6 * 8 * 3 * 11 / (14^2 * 13)
  a <- anchor_estimate(rs_hand, "Peto")
  expect_equal(a, c(eta = (3 - E_T) / V, se = 1 / sqrt(V)), tolerance = 1e-14)
})

test_that("A risk set of one subject adds nothing to the Peto variance", {
  rs <- list(n_T = c(10, 1), n_C = c(10, 0), O_T = c(1, 1), O = c(2, 1))
  V <- 10 * 10 * 2 * 18 / (20^2 * 19)
  expect_equal(anchor_estimate(rs, "Peto")[["se"]], 1 / sqrt(V),
               tolerance = 1e-14)
})

test_that("Undefined estimators are not finite", {
  rs <- list(n_T = c(5, 4), n_C = c(5, 4), O_T = c(0, 0), O = c(1, 1))
  expect_false(is.finite(anchor_estimate(rs, "MH")[["eta"]]))
  expect_false(is.finite(anchor_estimate(rs, "Pike")[["eta"]]))
})
