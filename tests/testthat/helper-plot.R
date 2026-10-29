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
