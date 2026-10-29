# Tests for the internal function plot_legend()

test_that("The legend is drawn outside or inside the panel, or not at all", {
  df <- data.frame(x = 1:3, y = c(2, 1, 3), g = c("a", "b", "c"))
  base <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, colour = g)) +
    ggplot2::geom_point() + plot_theme()
  add <- function(pos) {
    base + plot_legend(pos, k = 3, aesthetics = "colour")
  }
  on_null_device({
    expect_identical(legend_boxes(add("bottom")), "guide-box-bottom")
    expect_identical(legend_boxes(add("top")), "guide-box-top")
    expect_identical(legend_boxes(add("right")), "guide-box-right")
    expect_identical(legend_boxes(add("left")), "guide-box-left")
    expect_identical(legend_boxes(add("none")), character(0))
    expect_identical(legend_boxes(add(c(1, 1))), "guide-box-inside")
    expect_silent(print(add(c(0, 0))))
    expect_silent(print(base + plot_legend("bottom", legend_ncol = 1,
                                           aesthetics = "colour")))
  })
})

test_that("The result has two parts, the theme and the guides", {
  expect_length(plot_legend("bottom", k = 5), 2L)
  expect_length(plot_legend(c(0.5, 0.5), legend_ncol = 2), 2L)
})

test_that("Invalid positions and numbers of columns give errors", {
  expect_error(plot_legend("middle"), "legend_pos")
  expect_error(plot_legend(NA_character_), "legend_pos")
  expect_error(plot_legend(c("bottom", "top")), "legend_pos")
  expect_error(plot_legend(1), "legend_pos")
  expect_error(plot_legend(c(0.5, 2)), "legend_pos")
  expect_error(plot_legend(c(0.5, NA)), "legend_pos")
  expect_error(plot_legend("bottom", legend_ncol = 0), "legend_ncol")
  expect_error(plot_legend("bottom", legend_ncol = c(1, 2)), "legend_ncol")
  expect_error(plot_legend("bottom", legend_ncol = "2"), "legend_ncol")
})
