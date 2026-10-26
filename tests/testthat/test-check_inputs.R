# Tests for the internal function check_inputs()

test_that("The treatment indicator is group != control", {
  expect_identical(check_inputs(c(1, 2, 3), c(1, 0, 1), c("C", "T", "T"),
                                "C", NULL),
                   c(FALSE, TRUE, TRUE))
  expect_identical(check_inputs(c(1, 2), c(1, 0), factor(c("a", "b")),
                                "b", c(1, 1)),
                   c(TRUE, FALSE))
})

test_that("Invalid data give informative errors", {
  expect_error(check_inputs(c(1, 2), c(1, 0, 1), c("C", "T"), "C", NULL),
               "same length")
  expect_error(check_inputs(c(1, 2), c(1, 0), c("C", "T"), "C", 1),
               "strata")
  expect_error(check_inputs(c(1, 2), c(1, 0), c("C", "T"), c("C", "T"),
                            NULL),
               "single value")
  expect_error(check_inputs(c("1", "2"), c(1, 0), c("C", "T"), "C", NULL),
               "numeric")
  expect_error(check_inputs(c(1, NA), c(1, 0), c("C", "T"), "C", NULL),
               "missing")
  expect_error(check_inputs(c(1, 2), c(1, 0), c("C", "T"), "C", c(1, NA)),
               "missing")
  expect_error(check_inputs(c(1, 2), c(1, 2), c("C", "T"), "C", NULL),
               "coded")
  expect_error(check_inputs(c(1, 2), c(1, 0), c("C", "T"), "X", NULL),
               "control")
})
