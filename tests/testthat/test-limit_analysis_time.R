# Tests for the internal function limit_analysis_time()

# Expected event fraction at calendar time tau, by numerical integration over
# the accrual time (independent of the closed form in the function)
event_fraction <- function(tau, hazard_T, hazard_C, pi_T, a, dropout) {
  lam <- c(hazard_T, hazard_C)
  pi_j <- c(pi_T, 1 - pi_T)
  total <- 0
  for (j in 1:2) {
    r <- lam[j] + dropout
    total <- total + pi_j[j] * stats::integrate(
      function(u) lam[j] / r * (1 - exp(-r * (tau - u))) / a,
      0, min(tau, a), rel.tol = 1e-10
    )$value
  }
  total
}

test_that("The event fraction at the analysis time is the target", {
  cases <- list(
    list(hT = 0.08, hC = 0.04, piT = 0.5, a = 12, e = 0.6, d = 0),
    list(hT = 0.2, hC = 0.1, piT = 2 / 3, a = 24, e = 0.3, d = 0.02)
  )
  taus <- numeric(2)
  for (k in 1:2) {
    cs <- cases[[k]]
    taus[k] <- limit_analysis_time(cs$hT, cs$hC, cs$piT, cs$a, cs$e, cs$d)
    expect_equal(event_fraction(taus[k], cs$hT, cs$hC, cs$piT, cs$a, cs$d),
                 cs$e, tolerance = 1e-9)
  }
  # The first case is analysed after the end of accrual, the second before
  expect_gt(taus[1], 12)
  expect_lt(taus[2], 24)
})

test_that("An unreachable event fraction gives an error", {
  expect_error(limit_analysis_time(0.1, 0.05, 0.5, 10, 0.75, 0.05),
               "cannot be reached")
  expect_error(limit_analysis_time(0.1, 0.05, 0.5, 10, 1, 0),
               "cannot be reached")
})
