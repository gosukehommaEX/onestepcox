# Tests for the internal function plot_data_EstHR()

# A result of EstHR() written by hand, with the column names of EstHR()
hand_est <- function(values, method) {
  x <- matrix(values, nrow = length(method),
              dimnames = list(method, c("coef", "exp(coef)", "se(coef)",
                                        "lower .95", "upper .95")))
  attr(x, "method") <- method
  x
}

test_that("The plotted values follow the estimates of EstHR()", {
  est <- small_est()
  d <- plot_data_EstHR(est)
  expect_identical(d$method, attr(est, "method"))
  expect_identical(d$label, c("MH-Halley", "MH", "Pike", "Peto", "Cox"))
  expect_equal(d$HR, unname(est[, "exp(coef)"]))
  expect_equal(d$lower, unname(est[, "lower .95"]))
  expect_equal(d$upper, unname(est[, "upper .95"]))
  expect_equal(d$y, 5:1)
  expect_true(all(grepl("^[0-9.]+ \\([0-9.]+, [0-9.]+\\)$", d$text)))
  expect_equal(attr(d, "cox_hr"), est["Cox", "exp(coef)"])
  d <- plot_data_EstHR(est, labels = c("A", "B", "C", "D", "E"))
  expect_identical(d$label, c("A", "B", "C", "D", "E"))
})

test_that("The values are written with the given significant digits", {
  x <- hand_est(c(log(2), 2, 0.1, 1.5, 2.5), "MH")
  expect_identical(plot_data_EstHR(x)$text, "2.00 (1.50, 2.50)")
  expect_identical(plot_data_EstHR(x, value_digits = 2)$text,
                   "2.0 (1.5, 2.5)")
  expect_equal(plot_data_EstHR(x)$y, 1)
  expect_null(attr(plot_data_EstHR(x), "cox_hr"))
  # Values without decimals have no trailing point
  x <- hand_est(c(log(1234), 1234, 0.1, 1000, 1500), "MH")
  expect_identical(plot_data_EstHR(x)$text, "1234 (1000, 1500)")
})

test_that("A missing estimate is kept as NA, and Cox_NR is the reference", {
  x <- hand_est(c(log(2), NA, 2, NA, 0.1, NA, 1.5, NA, 2.5, NA),
                c("Cox_NR", "Pike"))
  d <- plot_data_EstHR(x)
  expect_identical(d$text[2L], "NA")
  expect_true(is.na(d$HR[2L]))
  expect_equal(d$y, 2:1)
  expect_equal(attr(d, "cox_hr"), 2)
  # "Cox" is preferred to "Cox_NR"
  x <- hand_est(c(log(2), log(3), 2, 3, 0.1, 0.1, 1.5, 2, 2.5, 4),
                c("Cox_NR", "Cox"))
  expect_equal(attr(plot_data_EstHR(x), "cox_hr"), 3)
})

test_that("Invalid input gives informative errors", {
  na_est <- EstHR(c(1, 2, 3, 4), c(0, 0, 0, 0), c("C", "T", "C", "T"), "C",
                  "MH")
  expect_error(plot_data_EstHR(na_est), "no estimate")
  expect_error(plot_data_EstHR(1:5), "output of EstHR")
  expect_error(plot_data_EstHR(small_est(), labels = c("A", "B")), "labels")
  expect_error(plot_data_EstHR(small_est(), value_digits = 0),
               "value_digits")
  expect_error(plot_data_EstHR(small_est(), value_digits = 2.5),
               "value_digits")
})
