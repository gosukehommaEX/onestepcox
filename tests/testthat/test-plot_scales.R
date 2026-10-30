# Tests for the internal function plot_scales()

test_that("Each method is mapped to its color and shape", {
  df <- data.frame(x = 1:3, y = c(2, 1, 3),
                   m = factor(c("B", "A", "C"), levels = c("A", "B", "C")))
  sc <- plot_scales(c("A", "B", "C"), c("a", "b", "c"),
                    col = c("red", "green", "blue"), shape = c(1, 2, 3),
                    linetype = c(1, 2, 3))
  expect_length(sc, 3L)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, colour = m,
                                        shape = m)) +
    ggplot2::geom_point() + sc[1:2]
  ld <- ggplot2::layer_data(p, 1L)
  expect_identical(ld$colour, c("green", "red", "blue"))
  expect_equal(ld$shape, c(2, 1, 3))
  # All scales have all methods as limits, also when a layer lacks some
  for (s in sc) expect_identical(s$get_limits(), c("A", "B", "C"))
})

test_that("The default colors are those of plot_colors(), recycled", {
  df <- data.frame(x = 1:3, y = c(2, 1, 3), m = c("A", "B", "C"))
  sc <- plot_scales(c("A", "B", "C"), c("a", "b", "c"),
                    aesthetics = "colour")
  expect_length(sc, 1L)
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, colour = m)) +
    ggplot2::geom_point() + sc
  expect_identical(ggplot2::layer_data(p, 1L)$colour, plot_colors(3))
  sc <- plot_scales(c("A", "B", "C"), c("a", "b", "c"), col = "black",
                    aesthetics = "colour")
  p <- ggplot2::ggplot(df, ggplot2::aes(x = x, y = y, colour = m)) +
    ggplot2::geom_point() + sc
  expect_identical(ggplot2::layer_data(p, 1L)$colour, rep("black", 3))
})

test_that("Empty values give informative errors", {
  expect_error(plot_scales("A", "a", shape = numeric(0)), "shape")
  expect_error(plot_scales("A", "a", linetype = character(0)), "linetype")
  expect_error(plot_scales("A", "a", col = character(0)), "col")
})
