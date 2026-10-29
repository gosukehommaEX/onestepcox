# Tests for plot.CompareOC()

test_that("Every quantity gives a ggplot that is drawn silently", {
  res <- small_compare()
  on_null_device({
    for (w in c("deviation", "bias", "rmse", "coverage", "se_ratio",
                "time")) {
      p <- plot(res, what = w)
      expect_true(inherits(p, "ggplot"))
      expect_silent(print(p))
    }
  })
})

test_that("The points are at the summaries of the methods", {
  res <- small_compare()
  s <- res$summary
  # No reference line for the RMSE: the layers are the lines and the points
  pts <- ggplot2::layer_data(plot(res, what = "rmse"), 2L)
  expect_equal(sort(pts$y), sort(s$rmse))
  # The bias has the reference line at 0 first
  p <- plot(res, what = "bias", methods = "Pike")
  expect_equal(ggplot2::layer_data(p, 1L)$yintercept, 0)
  expect_equal(sort(ggplot2::layer_data(p, 3L)$y),
               sort(s$bias[s$method == "Pike"]))
})

test_that("The legend is below the panel by default", {
  res <- small_compare()
  on_null_device({
    expect_identical(legend_boxes(plot(res)), "guide-box-bottom")
    expect_identical(legend_boxes(plot(res, what = "bias")),
                     "guide-box-bottom")
    expect_identical(legend_boxes(plot(res, legend_pos = "none")),
                     character(0))
    expect_identical(legend_boxes(plot(res, what = "rmse",
                                       legend_pos = c(0.98, 0.98))),
                     "guide-box-inside")
  })
})

test_that("plot() accepts the graphical options", {
  res <- small_compare()
  on_null_device({
    expect_silent(print(plot(res, methods = c("MH_Halley", "Pike"),
                             col = c("black", "red"), base_size = 13,
                             legend_pos = "right", ylim = c(1e-18, 1),
                             main = "Deviation", shape = c(1, 2),
                             size = 1)))
    expect_silent(print(plot(res, what = "coverage", methods = "Cox",
                             linetype = 2, linewidth = 0.5, shape = 1,
                             size = 3, legend_pos = "none",
                             ylim = c(0, 1))))
    expect_silent(print(plot(res, labels = c("Pike anchor", "Proposed"))))
    expect_silent(print(plot(res, what = "rmse", legend_pos = "top",
                             legend_ncol = 1, xlab = "", ylab = "RMSE")))
  })
})

test_that("Invalid input gives informative errors", {
  res <- small_compare()
  expect_error(plot(res, methods = "Peto"), "methods")
  expect_error(plot(res, what = "power"), "should be one of")
  expect_error(plot(res, legend_pos = "middle"), "legend_pos")
  expect_error(plot(res, labels = "Pike"), "labels")
  # The deviations are on a logarithmic axis, the bias on a linear one
  expect_error(plot(res, ylim = c(0, 1)), "ylim")
  expect_error(plot(res, what = "bias", ylim = c(1, -1)), "ylim")
  expect_true(inherits(plot(res, what = "bias", ylim = c(-1, 1)), "ggplot"))
})
