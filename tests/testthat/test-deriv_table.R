# Tests for the internal function deriv_table()

rs_hand <- list(n_T = c(5, 4), n_C = c(5, 3), O_T = c(1, 2), O = c(2, 3))

test_that("Breslow keeps one term per event time with weight O_k", {
  dt <- deriv_table(rs_hand, "breslow")
  expect_equal(dt$n_T, c(5, 4))
  expect_equal(dt$n_C, c(5, 3))
  expect_equal(dt$w, c(2, 3))
  expect_equal(dt$O_T, 3)
})

test_that("Efron splits each tie block into O_k reduced risk sets", {
  dt <- deriv_table(rs_hand, "efron")
  # Block 1: O = 2, O_T = 1, O_C = 1; r = 0, 1
  # Block 2: O = 3, O_T = 2, O_C = 1; r = 0, 1, 2
  expect_equal(dt$n_T, c(5, 5 - 1 / 2, 4, 4 - 2 / 3, 4 - 4 / 3))
  expect_equal(dt$n_C, c(5, 5 - 1 / 2, 3, 3 - 1 / 3, 3 - 2 / 3))
  expect_equal(dt$w, rep(1, 5))
  expect_equal(dt$O_T, 3)
})

test_that("Without ties, Efron and Breslow give the same terms", {
  rs <- list(n_T = c(3, 2), n_C = c(2, 2), O_T = c(1, 0), O = c(1, 1))
  expect_identical(deriv_table(rs, "efron"), deriv_table(rs, "breslow"))
})
