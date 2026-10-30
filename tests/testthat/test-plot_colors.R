# Tests for the internal function plot_colors()

test_that("The default palette and given colors are recycled", {
  expect_identical(plot_colors(2), c("#0072B2", "#D55E00"))
  expect_length(plot_colors(10), 10)
  expect_identical(plot_colors(10)[9], "#0072B2")
  expect_identical(plot_colors(3, "red"), rep("red", 3))
  expect_identical(plot_colors(3, c("red", "blue")), c("red", "blue", "red"))
  expect_error(plot_colors(2, character(0)), "col")
})
