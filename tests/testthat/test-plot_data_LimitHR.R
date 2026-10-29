# Tests for the internal function plot_data_LimitHR()

test_that("The bias is arranged by method and hazard ratio", {
  lim <- small_limit()
  d <- plot_data_LimitHR(lim, c("Pike", "MH_Halley"), "signed")
  expect_identical(attr(d, "methods"), c("Pike", "MH_Halley"))
  expect_identical(attr(d, "labels"), c("Pike", "MH-Halley"))
  expect_equal(d$HR, rep(c(0.25, 0.5, 1, 2, 4), 2))
  expect_identical(d$method, rep(c("Pike", "MH_Halley"), each = 5))
  expect_equal(d$bias[d$method == "Pike"], lim$bias[lim$method == "Pike"])
  expect_equal(d$value, d$bias)
})

test_that("The absolute bias is floored", {
  lim <- small_limit()
  d <- plot_data_LimitHR(lim, c("Pike", "MH_Halley"), "abs", floor = 1e-12)
  expect_equal(d$value, pmax(abs(d$bias), 1e-12))
  expect_true(all(d$value >= 1e-12))
  d <- plot_data_LimitHR(lim, "Peto", "abs", labels = "Peto (1972)")
  expect_identical(unique(d$label), "Peto (1972)")
})

test_that("Invalid input gives informative errors", {
  lim <- small_limit()
  expect_error(plot_data_LimitHR(lim, "PiHE", "abs"), "methods")
  expect_error(plot_data_LimitHR(lim, character(0), "abs"), "methods")
  expect_error(plot_data_LimitHR(lim, c("Pike", "Peto"), "abs",
                                 labels = "Pike"), "labels")
  expect_error(plot_data_LimitHR(data.frame(HR = 1), "Pike", "abs"),
               "columns")
  expect_error(plot_data_LimitHR(lim, "Pike", "abs", floor = -1), "floor")
})
