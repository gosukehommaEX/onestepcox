# Tests for the internal function method_label()

test_that("Every method of EstHR() has a description", {
  lab <- method_label(all_methods)
  expect_named(lab, all_methods)
  expect_false(anyNA(lab))
  expect_match(lab[["MH_Halley"]], "proposed")
  expect_identical(method_label("PiHE"), c(PiHE = NA_character_))
})
