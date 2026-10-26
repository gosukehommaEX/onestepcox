# Tests for the internal function risk_sets()

test_that("Counts are doubles and stratum is NULL without strata", {
  rs <- risk_sets(c(2, 1, 3), c(1, 1, 0), c(TRUE, FALSE, TRUE))
  expect_null(rs$stratum)
  expect_type(rs$n_T, "double")
  expect_type(rs$O, "double")
  expect_equal(rs$time, c(1, 2))
  expect_equal(rs$n_T, c(2, 2))
  expect_equal(rs$n_C, c(1, 0))
  expect_equal(rs$O_T, c(0, 1))
  expect_equal(rs$O, c(1, 1))
})

test_that("Logical event indicators give the same summaries", {
  a <- risk_sets(c(2, 1, 3, 3), c(1, 1, 0, 1), c(TRUE, FALSE, TRUE, FALSE))
  b <- risk_sets(c(2, 1, 3, 3), c(TRUE, TRUE, FALSE, TRUE),
                 c(TRUE, FALSE, TRUE, FALSE))
  expect_identical(a, b)
})

test_that("A censored time equal to an event time stays in the risk set", {
  # The censored subject at time 2 precedes the event at time 2
  rs <- risk_sets(c(2, 2, 3), c(0, 1, 1), c(FALSE, FALSE, TRUE))
  expect_equal(rs$n_C, c(2, 0))
  expect_equal(rs$n_T, c(1, 1))
})

test_that("is_sorted = TRUE checks the order of the rows", {
  expect_error(risk_sets(c(1, 2, 3), c(1, 1, 1), c(TRUE, FALSE, TRUE),
                         strata = c("a", "b", "a"), is_sorted = TRUE),
               "contiguous")
  expect_error(risk_sets(c(2, 1, 3), c(1, 1, 1), c(TRUE, FALSE, TRUE),
                         strata = c("a", "a", "b"), is_sorted = TRUE),
               "within strata")
  # A decrease of time at the start of a new stratum is allowed
  rs <- risk_sets(c(2, 3, 1), c(1, 1, 1), c(TRUE, FALSE, TRUE),
                  strata = c("a", "a", "b"), is_sorted = TRUE)
  expect_equal(rs$stratum, c("a", "a", "b"))
  expect_equal(rs$time, c(2, 3, 1))
})

test_that("Large counts do not overflow", {
  n <- 100000L
  rs <- risk_sets(seq_len(n), rep(1, n), rep(c(TRUE, FALSE), n / 2))
  expect_equal(rs$n_T[1] * rs$n_C[1] * rs$O[1] * (rs$n_T[1] + rs$n_C[1]),
               5e4 * 5e4 * 1 * 1e5)
})
