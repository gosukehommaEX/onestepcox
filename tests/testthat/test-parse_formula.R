# Tests for parse_formula()

toy <- data.frame(
  t = c(5, 3, 8, 2, 7, 4),
  e = c(1, 0, 1, 1, 0, 1),
  g = factor(c("b", "a", "b", "a", "b", "a"), levels = c("b", "a")),
  s = c(1, 1, 2, 2, 1, 2),
  u = c("x", "y", "x", "y", "y", "x"),
  stringsAsFactors = FALSE
)

test_that("The variables are taken from the data in their order", {
  f <- parse_formula(Surv(t, e) ~ g, toy)
  expect_identical(f$time, toy$t)
  expect_identical(f$event, toy$e)
  expect_identical(f$group, toy$g)
  expect_null(f$strata)
  # The first level of a factor is the control
  expect_identical(f$control, "b")
  expect_identical(f$n_removed, 0L)
})

test_that("The data may be given as the second argument or by name", {
  a <- parse_formula(Surv(t, e) ~ g, data = toy)
  b <- parse_formula(Surv(t, e) ~ g, event = toy)
  expect_identical(a, b)
})

test_that("Surv() and strata() may be named with the package", {
  a <- parse_formula(Surv(t, e) ~ g + strata(s), toy)
  b <- parse_formula(survival::Surv(t, e) ~ survival::strata(s) + g, toy)
  expect_identical(a, b)
  expect_identical(a$strata, toy$s)
})

test_that("Named arguments of Surv() are matched by name", {
  a <- parse_formula(Surv(t, e) ~ g, toy)
  b <- parse_formula(Surv(event = e, time = t) ~ g, toy)
  d <- parse_formula(Surv(event = e, t) ~ g, toy)
  expect_identical(a, b)
  expect_identical(a, d)
})

test_that("Several strata variables define strata by their combinations", {
  a <- parse_formula(Surv(t, e) ~ g + strata(s, u), toy)
  b <- parse_formula(Surv(t, e) ~ g + strata(s) + strata(u), toy)
  expect_identical(a$strata, paste(toy$s, toy$u, sep = ":"))
  expect_identical(a, b)
})

test_that("The control of other group variables is the smallest value", {
  expect_identical(parse_formula(Surv(t, e) ~ u, toy)$control, "x")
  expect_identical(parse_formula(Surv(t, e) ~ s, toy)$control, 1)
  # Expressions on the right-hand side; I() is removed from the values
  f <- parse_formula(Surv(t, e) ~ I(t > 4), toy)
  expect_identical(f$group, toy$t > 4)
  expect_identical(f$control, FALSE)
  f <- parse_formula(Surv(t, e) ~ factor(u, levels = c("y", "x")), toy)
  expect_identical(f$control, "y")
})

test_that("An event indicator coded 1 and 2 is recoded to 0 and 1", {
  d <- toy
  d$e2 <- d$e + 1
  f <- parse_formula(Surv(t, e2) ~ g, d)
  expect_equal(f$event, toy$e)
})

test_that("Rows with missing values are removed", {
  d <- toy
  d$t[2] <- NA
  d$s[5] <- NA
  f <- parse_formula(Surv(t, e) ~ g + strata(s), d)
  expect_identical(f$n_removed, 2L)
  expect_identical(f$time, toy$t[-c(2, 5)])
  expect_identical(f$strata, toy$s[-c(2, 5)])
})

test_that("The variables may come from the environment of the formula", {
  tt <- toy$t
  ee <- toy$e
  gg <- toy$u
  f <- parse_formula(Surv(tt, ee) ~ gg)
  expect_identical(f$time, toy$t)
  expect_identical(f$group, toy$u)
})

test_that("Invalid formulas and arguments give informative errors", {
  pf <- parse_formula
  expect_error(pf(t ~ g, toy), "Surv")
  expect_error(pf(~ g, toy), "must have the form")
  expect_error(pf(Surv(t, t, e) ~ g, toy), "right-censored")
  expect_error(pf(Surv(t, e, type = "right") ~ g, toy), "right-censored")
  expect_error(pf(Surv(t, e) ~ g + u, toy), "one group")
  expect_error(pf(Surv(t, e) ~ g * u, toy), "one group")
  expect_error(pf(Surv(t, e) ~ g:u, toy), "one group")
  expect_error(pf(Surv(t, e) ~ ., toy), "one group")
  expect_error(pf(Surv(t, e) ~ 1, toy), "one group")
  expect_error(pf(Surv(t, e) ~ strata(s), toy), "group variable")
  expect_error(pf(Surv(t, e) ~ g + strata(s, na.group = TRUE), toy),
               "variables only")
  expect_error(pf(Surv(t, e) ~ t, toy), "exactly two")
  expect_error(pf(Surv(t, e) ~ g, toy, group = toy$g), "by name")
  expect_error(pf(Surv(t, e) ~ g, toy, strata = toy$s), "strata")
  expect_error(pf(Surv(t, e) ~ g, data = toy, event = toy), "second")
  expect_error(pf(Surv(t, e) ~ g, data = 1:3), "data frame")
})
