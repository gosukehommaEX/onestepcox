# Tests for plot.LimitHR()

test_that("Both scales are plotted silently", {
  lim <- small_limit()
  on_null_device({
    expect_silent(plot(lim))
    expect_silent(plot(lim, scale = "signed"))
  })
})

test_that("plot() returns the bias of the methods invisibly", {
  lim <- small_limit()
  on_null_device({
    v <- withVisible(plot(lim, methods = c("Pike", "MH_Halley")))
    expect_false(v$visible)
    expect_identical(colnames(v$value), c("Pike", "MH_Halley"))
    expect_equal(v$value[, "Pike"], lim$bias[lim$method == "Pike"],
                 ignore_attr = TRUE)
  })
})

test_that("plot() accepts the graphical options", {
  lim <- small_limit()
  on_null_device({
    expect_silent(plot(lim, methods = "Peto", col = "black", lty = 3,
                       lwd = 1, pch = 2, cex = 1, base_size = 0.8,
                       legend_pos = "none", main = "Peto", xlab = "HR",
                       ylab = "bias"))
    expect_silent(plot(lim, legend_pos = "bottomright", legend_ncol = 2,
                       xlim = c(0.2, 5), ylim = c(1e-10, 1)))
    expect_silent(plot(lim, scale = "signed", legend_pos = c(0.3, 0),
                       floor = 1e-12, las = 1))
  })
})

test_that("plot() restores the graphical parameters", {
  lim <- small_limit()
  on_null_device({
    old <- graphics::par("cex.lab")
    plot(lim, base_size = 2)
    expect_identical(graphics::par("cex.lab"), old)
  })
})

test_that("Invalid input gives informative errors", {
  lim <- small_limit()
  on_null_device({
    expect_error(plot(lim, methods = "PiHE"), "methods")
    expect_error(plot(lim, scale = "log"), "should be one of")
    expect_error(plot(lim, legend_pos = 3), "legend_pos")
  })
})
