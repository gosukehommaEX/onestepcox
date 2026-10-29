# Tests for print.LimitHR()

test_that("print() shows the design, the analysis times, and the bias", {
  lim <- small_limit()
  out <- utils::capture.output(v <- withVisible(print(lim)))
  expect_false(v$visible)
  expect_identical(v$value, lim)
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "n_T : n_C = 100 : 100, a = 12, E = 160", fixed = TRUE)
  expect_match(txt, "HR = 0.25", fixed = TRUE)
  expect_match(txt, "Analysis time:", fixed = TRUE)
  expect_match(txt, "Bias (limit - log HR):", fixed = TRUE)
  expect_equal(sum(grepl("^Pike_Chebyshev ", out)), 1)
})

test_that("print() shows the bias of each method and hazard ratio", {
  lim <- small_limit()
  out <- utils::capture.output(print(lim, digits = 6))
  b <- lim$bias[lim$method == "Pike" & lim$HR == 4]
  row <- out[grepl("^Pike ", out)]
  expect_length(row, 1)
  expect_match(row, formatC(b, digits = 6, format = "g"), fixed = TRUE)
})

test_that("Without its columns, a LimitHR object prints as a data frame", {
  lim <- small_limit()
  part <- lim[, c("HR", "limit")]
  out <- utils::capture.output(print(part))
  expect_match(out[1], "HR")
  expect_match(out[1], "limit")
})
