# Tests for the internal function plot_data_CompareOC()

test_that("Deviations are matched by replicate and floored", {
  res <- small_compare()
  e <- res$estimates
  d <- plot_data_CompareOC(res, "deviation")
  expect_identical(attr(d, "reference"), "Cox")
  expect_identical(attr(d, "methods"), c("Pike", "MH_Halley"))
  expect_identical(attr(d, "labels"), c("Pike", "MH-Halley"))
  expect_identical(attr(d, "scenarios"), c("S1", "S2"))
  ref <- e$coef[e$label == "S1" & e$method == "Cox"]
  pk <- e$coef[e$label == "S1" & e$method == "Pike"]
  expect_equal(d$value[d$scenario == "S1" & d$method == "Pike"],
               pmax(abs(pk - ref), 1e-16))
  expect_identical(d$label[d$method == "MH_Halley"][1L], "MH-Halley")
  d <- plot_data_CompareOC(res, "deviation", floor = 1e-3)
  expect_true(all(d$value >= 1e-3))
})

test_that("Summaries are taken by scenario and method", {
  res <- small_compare()
  s <- res$summary
  d <- plot_data_CompareOC(res, "bias")
  expect_identical(attr(d, "methods"), c("Cox", "Pike", "MH_Halley"))
  expect_identical(d$scenario, rep(c("S1", "S2"), 3))
  expect_identical(d$method, rep(c("Cox", "Pike", "MH_Halley"), each = 2))
  expect_equal(d$value[d$scenario == "S2" & d$method == "Pike"],
               s$bias[s$label == "S2" & s$method == "Pike"])
  d <- plot_data_CompareOC(res, "se_ratio")
  expect_identical(attr(d, "methods"), c("Pike", "MH_Halley"))
  key <- match(paste(d$scenario, d$method), paste(s$label, s$method))
  expect_equal(d$value, s$se_ratio_median[key])
  d <- plot_data_CompareOC(res, "time", methods = "MH_Halley",
                           labels = "Proposed")
  expect_identical(d$label, c("Proposed", "Proposed"))
  expect_equal(d$value, s$time_per_sim[s$method == "MH_Halley"])
})

test_that("Invalid input gives informative errors", {
  res <- small_compare()
  expect_error(plot_data_CompareOC(res, "bias", methods = "Peto"), "methods")
  expect_error(plot_data_CompareOC(res, "bias", methods = character(0)),
               "no method")
  expect_error(plot_data_CompareOC(res, "deviation", labels = "Pike"),
               "labels")
  expect_error(plot_data_CompareOC(list(summary = 1), "bias"),
               "output of CompareOC")
  expect_error(plot_data_CompareOC(res, "deviation", floor = 0), "floor")
  expect_error(plot_data_CompareOC(res, "deviation", floor = c(1, 2)),
               "floor")
})
