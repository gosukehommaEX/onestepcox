# Tests for the internal function estimate_hr()

test_that("The values of all methods match the Python values", {
  d <- test_data("gbsg_grade")
  trt <- d$group != d$control
  est <- estimate_hr(d$time, d$event, trt, d$strata, all_methods, "efron",
                     FALSE, 1e-12)
  expect_identical(dimnames(est), list(all_methods, c("coef", "se")))
  ref <- expected_estimates[expected_estimates$data == "gbsg_grade" &
                              expected_estimates$ties == "efron", ]
  pos <- match(ref$method, all_methods)
  expect_false(anyNA(pos))
  expect_equal(unname(est[pos, "coef"]), ref$coef, tolerance = 1e-10)
  expect_equal(unname(est[pos, "se"]), ref$se, tolerance = 1e-10)
})

test_that("The order of the methods does not change the values", {
  d <- test_data("gbsg")
  trt <- d$group != d$control
  a <- estimate_hr(d$time, d$event, trt, NULL,
                   c("Pike_Newton", "MH", "Cox_NR", "Peto"), "breslow",
                   FALSE, 1e-12)
  b <- estimate_hr(d$time, d$event, trt, NULL,
                   c("Cox_NR", "Peto", "Pike_Newton", "MH"), "breslow",
                   FALSE, 1e-12)
  expect_identical(a[rownames(b), ], b)
})

test_that("NA is returned for every method without an estimate", {
  # Events in the control group only
  est <- estimate_hr(c(1, 2, 3, 4), c(1, 0, 1, 0),
                     c(FALSE, TRUE, FALSE, TRUE), NULL,
                     c("MH", "Pike_Halley", "Cox_NR"), "efron", FALSE, 1e-12)
  expect_true(all(is.na(est)))
  # No subject: coxph() is not called
  est0 <- estimate_hr(numeric(0), numeric(0), logical(0), NULL,
                      c("MH", "Cox"), "efron", FALSE, 1e-12)
  expect_true(all(is.na(est0)))
  expect_identical(dim(est0), c(2L, 2L))
})
