# Tests for CompareOC()

scenarios_demo <- function(extra = NULL) {
  base <- data.frame(
    label    = c("S1", "S2"),
    n_T      = c(80, 100),
    n_C      = c(80, 100),
    hazard_C = c(log(2) / 18, log(2) / 12),
    HR       = c(0.7, 0.8),
    a        = c(12, 12),
    E        = c(60, 90),
    stringsAsFactors = FALSE
  )
  if (!is.null(extra)) base <- cbind(base, extra)
  base
}

test_that("Output has the documented structure", {
  res <- CompareOC(
    scenarios = scenarios_demo(),
    methods   = c("Cox", "Peto", "Pike", "MH_Halley"),
    nsim      = 30,
    seed      = 1
  )
  expect_named(res, c("estimates", "summary"), ignore.order = TRUE)
  expect_named(res$estimates,
               c("label", "simID", "method", "HR"),
               ignore.order = TRUE)
  expect_named(res$summary,
               c("label", "HR_true", "method", "mean", "bias",
                 "rmse", "time_per_sim"),
               ignore.order = TRUE)
  # 2 scenarios * 4 methods = 8 summary rows
  expect_equal(nrow(res$summary), 2 * 4)
})

test_that("Reference method has bias = 0 and rmse = 0 by construction", {
  res <- CompareOC(
    scenarios = scenarios_demo(),
    methods   = c("Cox", "Peto", "Pike", "MH_Halley"),
    nsim      = 30,
    seed      = 1,
    reference = "Cox"
  )
  ref_rows <- res$summary[res$summary$method == "Cox", ]
  expect_true(all(ref_rows$bias == 0))
  expect_true(all(ref_rows$rmse == 0))
})

test_that("Reference must be one of the supplied methods", {
  expect_error(
    CompareOC(
      scenarios = scenarios_demo(),
      methods   = c("Peto", "Pike", "MH_Halley"),
      nsim      = 10,
      seed      = 1,
      reference = "Cox"
    ),
    "reference"
  )
})

test_that("Missing required scenario columns triggers an error", {
  bad <- scenarios_demo()
  bad$E <- NULL
  expect_error(
    CompareOC(scenarios = bad, methods = c("Cox", "MH_Halley"),
              nsim = 10, seed = 1),
    "missing required columns"
  )
})

test_that("dropout_hazard column is optional and defaults to 0", {
  sc_no <- scenarios_demo()
  sc_with <- scenarios_demo(extra = data.frame(dropout_hazard = c(0, 0)))
  res_no   <- CompareOC(sc_no,   methods = c("Cox", "MH_Halley"),
                        nsim = 20, seed = 7)
  res_with <- CompareOC(sc_with, methods = c("Cox", "MH_Halley"),
                        nsim = 20, seed = 7)
  expect_equal(res_no$summary$mean, res_with$summary$mean,
               tolerance = 1e-10)
})

test_that("MH_Halley RMSE vs Cox is much smaller than Pike's", {
  # Sanity check of the article's main empirical claim, on a small grid.
  res <- CompareOC(
    scenarios = scenarios_demo(),
    methods   = c("Cox", "Pike", "MH_Halley"),
    nsim      = 100,
    seed      = 11,
    reference = "Cox"
  )
  for (lab in unique(res$summary$label)) {
    rmse_pike <- res$summary$rmse[res$summary$label == lab &
                                    res$summary$method == "Pike"]
    rmse_pihe <- res$summary$rmse[res$summary$label == lab &
                                    res$summary$method == "MH_Halley"]
    expect_lt(rmse_pihe, rmse_pike / 10)
  }
})
