# Tests for print.CompareOC()

test_that("print() shows the settings and a table for each scenario", {
  res <- small_compare()
  out <- utils::capture.output(v <- withVisible(print(res)))
  expect_false(v$visible)
  expect_identical(v$value, res)
  txt <- paste(out, collapse = "\n")
  expect_match(txt, "2 scenarios; 8 replicates each; reference: Cox",
               fixed = TRUE)
  expect_match(txt, "Scenario S1 (true HR = 0.6)", fixed = TRUE)
  expect_match(txt, "Scenario S2 (true HR = 2)", fixed = TRUE)
  expect_match(txt, "se_ratio", fixed = TRUE)
  # One row per scenario (more if the table is wrapped)
  expect_gte(sum(grepl("^MH_Halley ", out)), 2)
})

test_that("print() shows the list when the summary is missing", {
  res <- small_compare()
  res$summary <- NULL
  out <- utils::capture.output(print(res))
  expect_match(out[1], "estimates", fixed = TRUE)
})
