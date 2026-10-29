# Tests for plot.LimitHR()

test_that("Both scales give a ggplot that is drawn silently", {
  lim <- small_limit()
  on_null_device({
    p <- plot(lim)
    expect_true(inherits(p, "ggplot"))
    expect_silent(print(p))
    expect_silent(print(plot(lim, scale = "signed")))
  })
})

test_that("The points are at the bias of the methods", {
  lim <- small_limit()
  # Signed scale: the reference line at 0, the lines, and the points
  p <- plot(lim, methods = c("Pike", "Peto"), scale = "signed")
  pts <- ggplot2::layer_data(p, 3L)
  expect_equal(sort(pts$y),
               sort(lim$bias[lim$method == "Pike" | lim$method == "Peto"]))
  expect_equal(sort(unique(10^pts$x)), c(0.25, 0.5, 1, 2, 4))
  # Absolute scale: the lines and the points, on a logarithmic axis
  p <- plot(lim, methods = "Pike")
  pts <- ggplot2::layer_data(p, 2L)
  expect_equal(10^pts$y[order(pts$x)],
               pmax(abs(lim$bias[lim$method == "Pike"]), 1e-16))
})

test_that("The legend is below the panel by default", {
  lim <- small_limit()
  on_null_device({
    expect_identical(legend_boxes(plot(lim)), "guide-box-bottom")
    expect_identical(legend_boxes(plot(lim, legend_pos = "right")),
                     "guide-box-right")
    expect_identical(legend_boxes(plot(lim, legend_pos = c(0.02, 0.98))),
                     "guide-box-inside")
  })
})

test_that("plot() accepts the graphical options", {
  lim <- small_limit()
  on_null_device({
    expect_silent(print(plot(lim, methods = "Peto", col = "black",
                             linetype = 3, linewidth = 0.5, shape = 2,
                             size = 2, base_size = 9, legend_pos = "none",
                             main = "Peto", xlab = "HR", ylab = "bias")))
    expect_silent(print(plot(lim, methods = c("Pike", "Pike_Halley"),
                             labels = c("Pike", "Pike, Halley step"))))
    expect_silent(print(plot(lim, legend_pos = "bottom", legend_ncol = 2,
                             xlim = c(0.2, 5), ylim = c(1e-10, 1))))
    expect_silent(print(plot(lim, scale = "signed", floor = 1e-12,
                             legend_pos = c(0.02, 0.02))))
  })
})

test_that("Invalid input gives informative errors", {
  lim <- small_limit()
  expect_error(plot(lim, methods = "PiHE"), "methods")
  expect_error(plot(lim, scale = "log"), "should be one of")
  expect_error(plot(lim, legend_pos = 3), "legend_pos")
  expect_error(plot(lim, labels = "Pike"), "labels")
  expect_error(plot(lim, xlim = c(0, 4)), "xlim")
  expect_error(plot(lim, ylim = c(0, 1)), "ylim")
  expect_true(inherits(plot(lim, scale = "signed", ylim = c(-1, 1)),
                       "ggplot"))
})

test_that("The methods have different shapes and line types by default", {
  p <- plot(small_limit())
  expect_length(unique(ggplot2::layer_data(p, 1L)$linetype), 5L)
  expect_length(unique(ggplot2::layer_data(p, 2L)$shape), 5L)
})
