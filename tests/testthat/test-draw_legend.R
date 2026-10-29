# Tests for the internal function draw_legend()

test_that("The legend is drawn at a keyword or at coordinates, or not", {
  on_null_device({
    graphics::plot(1:2)
    expect_null(draw_legend("none", c("a", "b")))
    out <- draw_legend("topright", c("a", "b"), col = 1:2, lty = 1)
    expect_named(out, c("rect", "text"))
    out <- draw_legend(c(1.2, 1.8), c("a", "b"), col = 1:2, lty = 1)
    expect_equal(out$rect$left, 1.2)
    expect_equal(out$rect$top, 1.8)
  })
})

test_that("An invalid position gives an informative error", {
  on_null_device({
    graphics::plot(1:2)
    expect_error(draw_legend("middle", "a"), "legend_pos")
    expect_error(draw_legend(c(1, NA), "a"), "legend_pos")
    expect_error(draw_legend(1, "a"), "legend_pos")
  })
})
