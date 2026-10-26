# Tests for the internal function one_step()

d0 <- c(U = 1.5, I = 20, J = 4, L = -3)

test_that("The three steps follow their formulas", {
  U <- 1.5
  I <- 20
  J <- 4
  L <- -3
  for (type in c("Newton", "Halley", "Chebyshev")) {
    delta <- switch(type,
                    Newton = U / I,
                    Halley = 2 * U * I / (2 * I^2 + U * J),
                    Chebyshev = U / I - J * U^2 / (2 * I^3))
    st <- one_step(0.3, d0, type)
    expect_equal(st[["eta"]], 0.3 + delta, tolerance = 1e-14, info = type)
    expect_equal(st[["se"]], 1 / sqrt(I + J * delta + L * delta^2 / 2),
                 tolerance = 1e-14, info = type)
    expect_named(st, c("eta", "se"))
  }
})

test_that("Halley and Chebyshev agree to second order in U / I", {
  # The difference is of order (U / I)^3
  for (U in c(1e-1, 1e-2, 1e-3)) {
    d <- c(U = U, I = 1, J = 0.5, L = 0)
    diff <- one_step(0, d, "Halley")[["eta"]] -
      one_step(0, d, "Chebyshev")[["eta"]]
    expect_lt(abs(diff), U^3)
    expect_gt(abs(diff), U^3 / 100)
  }
})

test_that("A non-positive Halley denominator or information gives NA", {
  st <- one_step(0, c(U = -10, I = 1, J = 1, L = 0), "Halley")
  expect_true(is.na(st[["eta"]]))
  expect_true(is.na(st[["se"]]))
  st <- one_step(0, c(U = 3, I = 1, J = 0, L = -10), "Newton")
  expect_equal(st[["eta"]], 3)
  expect_true(is.na(st[["se"]]))
})
