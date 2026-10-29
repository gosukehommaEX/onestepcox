# Tests for plot.CompareOC()

test_that("Every quantity is plotted silently", {
  res <- small_compare()
  on_null_device({
    for (w in c("deviation", "bias", "rmse", "coverage", "se_ratio",
                "time")) {
      expect_silent(plot(res, what = w))
    }
  })
})

test_that("plot() returns the plotted values invisibly", {
  res <- small_compare()
  s <- res$summary
  e <- res$estimates
  on_null_device({
    v <- withVisible(plot(res, what = "bias"))
    expect_false(v$visible)
    expect_identical(dimnames(v$value),
                     list(c("S1", "S2"), c("Cox", "Pike", "MH_Halley")))
    expect_equal(v$value["S2", "Pike"],
                 s$bias[s$label == "S2" & s$method == "Pike"])
    d <- plot(res, what = "deviation")
    expect_named(d, c("S1:Pike", "S1:MH_Halley", "S2:Pike", "S2:MH_Halley"))
    ref <- e$coef[e$label == "S1" & e$method == "Cox"]
    pk <- e$coef[e$label == "S1" & e$method == "Pike"]
    expect_equal(d[["S1:Pike"]], pmax(abs(pk - ref), 1e-16))
  })
})

test_that("plot() accepts the graphical options", {
  res <- small_compare()
  on_null_device({
    expect_silent(plot(res, methods = c("MH_Halley", "Pike"),
                       col = c("black", "red"), base_size = 1.2,
                       legend_pos = "bottomleft", ylim = c(1e-18, 1),
                       main = "Deviation", boxwex = 0.5))
    expect_silent(plot(res, what = "coverage", methods = "Cox", lty = 2,
                       lwd = 1, pch = 1, cex = 1.5, legend_pos = "none",
                       ylim = c(0, 1)))
    expect_silent(plot(res, labels = c("Pike anchor", "Proposed")))
    expect_silent(plot(res, what = "rmse", legend_pos = c(1, 0.5),
                       legend_ncol = 3, xlab = "", ylab = "RMSE", las = 1))
  })
})

test_that("plot() restores the graphical parameters", {
  res <- small_compare()
  on_null_device({
    old <- graphics::par("cex.axis")
    plot(res, base_size = 1.5)
    expect_identical(graphics::par("cex.axis"), old)
  })
})

test_that("Invalid input gives informative errors", {
  res <- small_compare()
  on_null_device({
    expect_error(plot(res, methods = "Peto"), "methods")
    expect_error(plot(res, what = "power"), "should be one of")
    expect_error(plot(res, legend_pos = "middle"), "legend_pos")
    expect_error(plot(res, labels = "Pike"), "labels")
  })
})
