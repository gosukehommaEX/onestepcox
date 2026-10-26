# Tests for the internal function newton_cox()

gbsg_table <- function(ties) {
  d <- test_data("gbsg")
  deriv_table(risk_sets(d$time, d$event, d$group != d$control), ties)
}

test_that("The iteration converges to the Python root from several starts", {
  for (ti in c("efron", "breslow")) {
    ref <- expected_estimates[expected_estimates$data == "gbsg" &
                                expected_estimates$ties == ti &
                                expected_estimates$method == "Cox_NR", ]
    dt <- gbsg_table(ti)
    for (eta0 in c(-0.36, 0, 3, -5)) {
      fit <- newton_cox(eta0, dt, tol = 1e-12)
      expect_true(fit$converged, info = paste(ti, eta0))
      expect_equal(fit$eta, ref$coef, tolerance = 1e-10,
                   info = paste(ti, eta0))
      expect_equal(1 / sqrt(fit$info), ref$se, tolerance = 1e-10,
                   info = paste(ti, eta0))
    }
  }
})

test_that("Few iterations are needed from a nearby start", {
  fit <- newton_cox(-0.36, gbsg_table("efron"), tol = 1e-12)
  expect_lte(fit$iter, 4L)
})

test_that("Non-convergence is reported", {
  fit <- newton_cox(3, gbsg_table("efron"), tol = 1e-12, maxit = 2L)
  expect_false(fit$converged)
  expect_equal(fit$iter, 2L)
})
