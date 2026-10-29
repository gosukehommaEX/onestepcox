# Tests for plot.EstHR()

est_for_plot <- function() {
  EstHR(Surv(rfstime, status) ~ hormon + strata(grade),
        data = survival::gbsg,
        method = c("MH_Halley", "MH", "Pike", "Peto", "Cox"))
}

test_that("plot() draws silently and returns the plotted values invisibly", {
  est <- est_for_plot()
  on_null_device({
    expect_silent(v <- plot(est))
    expect_identical(v$method, attr(est, "method"))
    expect_equal(v$HR, unname(est[, "exp(coef)"]))
    expect_equal(v$lower, unname(est[, "lower .95"]))
    expect_equal(v$y, 5:1)
    expect_false(withVisible(plot(est))$visible)
  })
})

test_that("plot() accepts the graphical options", {
  est <- est_for_plot()
  single <- EstHR(Surv(rfstime, status) ~ hormon, data = survival::gbsg,
                  method = "MH_Halley")
  on_null_device({
    expect_silent(plot(est, col = "black", pch = 19, cex = 1, lwd = 3,
                       lty = 2, base_size = 1.3, main = "GBSG", xlab = "HR"))
    expect_silent(plot(est, show_values = FALSE, ref_lines = FALSE,
                       xlim = c(0.4, 1.2), value_digits = 2))
    expect_silent(plot(est, legend_pos = "topright", legend_ncol = 2,
                       legend_inset = 0.05))
    expect_silent(plot(est, legend_pos = c(0.5, 4)))
    expect_silent(plot(est, mar = c(4, 8, 2, 2), ref_col = "red", las = 1))
    expect_silent(plot(single))
    # One Cox estimate, which is also the reference line
    cox <- EstHR(Surv(rfstime, status) ~ hormon, data = survival::gbsg,
                 method = "Cox")
    expect_silent(v <- plot(cox))
    expect_equal(v$HR, cox[["exp(coef)"]])
  })
})

test_that("plot() restores the graphical parameters", {
  est <- est_for_plot()
  on_null_device({
    old <- graphics::par("mar")
    plot(est)
    expect_identical(graphics::par("mar"), old)
  })
})

test_that("Invalid input gives informative errors", {
  na_est <- EstHR(c(1, 2, 3, 4), c(0, 0, 0, 0), c("C", "T", "C", "T"), "C",
                  "MH")
  on_null_device({
    expect_error(plot(na_est), "no estimate")
    expect_error(plot(est_for_plot(), legend_pos = "middle"), "legend_pos")
  })
})
