# Tests for the internal function plot_labels()

test_that("Underscores become hyphens unless labels are given", {
  expect_identical(plot_labels(c("MH_Halley", "Pike", "Cox_NR")),
                   c("MH-Halley", "Pike", "Cox-NR"))
  expect_identical(plot_labels(c("MH", "Pike"), c("a", "b")), c("a", "b"))
  expect_identical(plot_labels("MH", 1), "1")
  expect_error(plot_labels(c("MH", "Pike"), "a"), "one label per method")
})
