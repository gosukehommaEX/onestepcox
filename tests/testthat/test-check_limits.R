# Tests for the internal function check_limits()

test_that("NULL and two increasing numbers are accepted", {
  expect_null(check_limits(NULL, "xlim"))
  expect_identical(check_limits(c(-1, 2), "ylim"), c(-1, 2))
  expect_identical(check_limits(c(0.5, 2), "xlim", positive = TRUE),
                   c(0.5, 2))
})

test_that("Other limits give an error that names the argument", {
  expect_error(check_limits(c(2, 1), "ylim"), "'ylim' must be two increasing")
  expect_error(check_limits(c(1, 1), "ylim"), "ylim")
  expect_error(check_limits(1, "ylim"), "ylim")
  expect_error(check_limits(c(1, NA), "ylim"), "ylim")
  expect_error(check_limits(c(1, Inf), "ylim"), "ylim")
  expect_error(check_limits(c("1", "2"), "ylim"), "ylim")
  expect_error(check_limits(c(0, 1), "xlim", positive = TRUE),
               "'xlim' must be two increasing positive numbers")
  expect_error(check_limits(c(-1, 1), "xlim", positive = TRUE), "xlim")
})
