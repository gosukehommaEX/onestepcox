# Tests for print.EstHR()

test_that("print() shows the groups, the estimates, and the methods", {
  est <- EstHR(Surv(rfstime, status) ~ hormon + strata(grade),
               data = survival::gbsg, method = c("MH_Halley", "Cox"))
  out <- utils::capture.output(v <- withVisible(print(est)))
  expect_false(v$visible)
  expect_identical(v$value, est)
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "Treatment (1):  n = 246, events = 94", fixed = TRUE)
  expect_match(txt, "n = 440, events = 205", fixed = TRUE)
  expect_match(txt, "Ties: efron; strata: 3", fixed = TRUE)
  expect_match(txt, "Pr(>|z|)", fixed = TRUE)
  expect_match(txt, "lower .95", fixed = TRUE)
  expect_match(txt, "Halley step from the Mantel-Haenszel estimate",
               fixed = TRUE)
  expect_match(txt, "Cox estimate by survival::coxph()", fixed = TRUE)
})

test_that("print() reports removed rows and a missing estimate", {
  est <- EstHR(Surv(time, status) ~ sex + strata(ph.ecog),
               data = survival::lung, method = "MH_Halley")
  expect_output(print(est), "1 observation deleted due to missingness",
                fixed = TRUE)
  na_est <- EstHR(c(1, 2, 3, 4), c(0, 0, 0, 0), c("C", "T", "C", "T"), "C",
                  "MH")
  expect_output(print(na_est), "Estimate not available")
})

test_that("print() uses digits and prints the numbers without attributes", {
  d <- test_data("ovarian")
  est <- EstHR(d$time, d$event, d$group, d$control, "MH_Halley")
  a <- utils::capture.output(print(est, digits = 3))
  b <- utils::capture.output(print(est, digits = 8))
  expect_false(identical(a, b))
  bare <- est
  attr(bare, "method") <- NULL
  out <- utils::capture.output(print(bare))
  expect_match(out[1], "coef")
  expect_false(any(grepl("Treatment", out)))
})
