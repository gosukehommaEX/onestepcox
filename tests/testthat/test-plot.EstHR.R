# Tests for plot.EstHR()

test_that("plot() returns a ggplot that is drawn silently", {
  p <- plot(small_est())
  expect_true(inherits(p, "ggplot"))
  on_null_device(expect_silent(print(p)))
})

test_that("The points and intervals are at the estimates, top to bottom", {
  est <- small_est()
  # Without reference lines the layers are the intervals and the points
  p <- plot(est, ref_lines = FALSE)
  pts <- ggplot2::layer_data(p, 2L)
  pts <- pts[order(-pts$y), ]
  expect_equal(10^pts$x, unname(est[, "exp(coef)"]))
  expect_equal(pts$y, 5:1)
  seg <- ggplot2::layer_data(p, 1L)
  seg <- seg[order(-seg$y), ]
  expect_equal(10^seg$x, unname(est[, "lower .95"]))
  expect_equal(10^seg$xend, unname(est[, "upper .95"]))
  # With them, the reference lines at 1 and at the Cox estimate come first
  p <- plot(est)
  expect_equal(10^ggplot2::layer_data(p, 1L)$xintercept, 1)
  expect_equal(10^ggplot2::layer_data(p, 2L)$xintercept,
               est["Cox", "exp(coef)"])
})

test_that("The legend is not drawn by default, and below the panel if asked", {
  est <- small_est()
  on_null_device({
    expect_identical(legend_boxes(plot(est)), character(0))
    expect_identical(legend_boxes(plot(est, legend_pos = "bottom")),
                     "guide-box-bottom")
    expect_identical(legend_boxes(plot(est, legend_pos = c(1, 1))),
                     "guide-box-inside")
  })
})

test_that("plot() accepts the graphical options", {
  est <- small_est()
  single <- EstHR(Surv(rfstime, status) ~ hormon, data = survival::gbsg,
                  method = "MH_Halley")
  cox <- EstHR(Surv(rfstime, status) ~ hormon, data = survival::gbsg,
               method = "Cox")
  on_null_device({
    expect_silent(print(plot(est, col = "black", shape = 19, size = 2,
                             linewidth = 1, linetype = 2, base_size = 14,
                             main = "GBSG", xlab = "HR")))
    expect_silent(print(plot(est, show_values = FALSE, ref_lines = FALSE,
                             xlim = c(0.4, 1.2), value_digits = 2)))
    expect_silent(print(plot(est, legend_pos = "bottom", legend_ncol = 2,
                             ref_col = "red")))
    expect_silent(print(plot(est, labels = c("A", "B", "C", "D", "E"))))
    expect_silent(print(plot(single)))
    expect_silent(print(plot(cox)))
    # Intervals beyond the limits are cut at the edge of the panel
    expect_silent(print(plot(est, xlim = c(0.6, 0.8))))
  })
})

test_that("Invalid input gives informative errors", {
  na_est <- EstHR(c(1, 2, 3, 4), c(0, 0, 0, 0), c("C", "T", "C", "T"), "C",
                  "MH")
  est <- small_est()
  expect_error(plot(na_est), "no estimate")
  expect_error(plot(est, legend_pos = "middle"), "legend_pos")
  expect_error(plot(est, labels = c("A", "B")), "labels")
  expect_error(plot(est, xlim = c(1, 0.5)), "xlim")
  expect_error(plot(est, xlim = c(0, 1)), "xlim")
  expect_error(plot(est, base_size = -1), "base_size")
})
