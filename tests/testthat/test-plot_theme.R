# Tests for the internal function plot_theme()

test_that("The theme is added to a plot that is drawn silently", {
  df <- data.frame(x = 1:3, y = c(2, 1, 3))
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y)) +
    ggplot2::geom_point() + plot_theme(14)
  on_null_device(expect_silent(print(p)))
})

test_that("An invalid base size gives an informative error", {
  expect_error(plot_theme(0), "base_size")
  expect_error(plot_theme(c(10, 12)), "base_size")
  expect_error(plot_theme("11"), "base_size")
  expect_error(plot_theme(NA_real_), "base_size")
})
