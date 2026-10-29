# Evaluate plotting code on a null PDF device, which writes no file, and
# close the device afterwards
on_null_device <- function(code) {
  grDevices::pdf(NULL)
  on.exit(grDevices::dev.off())
  code
}

# Small results of CompareOC() and LimitHR() shared by the tests of the
# print() and plot() methods
small_compare <- function() {
  scenarios <- data.frame(label = c("S1", "S2"), n_T = 60, n_C = 60,
                          hazard_C = log(2) / 18, HR = c(0.6, 2), a = 12,
                          E = 60, stringsAsFactors = FALSE)
  CompareOC(scenarios, c("Cox", "Pike", "MH_Halley"), nsim = 8, seed = 1)
}

small_limit <- function() {
  LimitHR(HR = c(0.25, 0.5, 1, 2, 4), hazard_C = log(2) / 18, n_T = 100,
          n_C = 100, a = 12, E = 160)
}

# Result of EstHR() with five methods shared by the tests of plot.EstHR()
# and plot_data_EstHR()
small_est <- function() {
  EstHR(Surv(rfstime, status) ~ hormon + strata(grade),
        data = survival::gbsg,
        method = c("MH_Halley", "MH", "Pike", "Peto", "Cox"))
}

# Names of the legend boxes of a ggplot that are drawn, such as
# "guide-box-bottom" for a legend below the panel; empty without a legend.
# Building the table of graphical objects needs an open device.
legend_boxes <- function(p) {
  g <- ggplot2::ggplotGrob(p)
  drawn <- !vapply(g$grobs, inherits, logical(1), what = "zeroGrob")
  g$layout$name[grepl("^guide-box", g$layout$name) & drawn]
}
