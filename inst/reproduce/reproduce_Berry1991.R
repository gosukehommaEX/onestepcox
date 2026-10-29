# Reproduction of Figure 1 of Berry, Kitchin and Mock (1991), "A comparison
# of two simple hazard ratio estimators based on the logrank test",
# Statistics in Medicine 10: 749-755.
#
# Design of the article (section 'Simulation study'): 100 subjects per group,
# all at risk from the start, no censoring before the end of the study. The
# group of the next death is drawn from a Bernoulli distribution with
# probability n_iA theta / (n_iA theta + n_iB) for group A. The study stops
# after 100, 80, or 60 deaths (50, 60, or 70 per cent censoring); if there
# is no death in group B by then, deaths are generated until one occurs in
# group B. Figure 1 shows the mean of the Peto and Pike estimates of the
# hazard ratio of group A to group B over 1000 simulated trials. With equal
# group sizes, the labels of the groups do not affect the comparison.
#
# The figure has no table, so its values were read from the PDF
# (Berry1991_figure1_digitized.csv; the method is described in
# Handover/checks_20261028/digitize_Berry1991_fig1.py of the project folder).
# Each value is compared with the mean of 10,000 simulated trials. The
# tolerance is 3 times the Monte Carlo standard error of the difference
# (1000 trials in the article and 10,000 here) plus the reading error.
#
# Run from the root of the package after devtools::load_all() (or with the
# installed package). Output: inst/reproduce/results/Berry1991_figure1.csv.
# Running time: a few minutes.

if (!"onestepcox" %in% loadedNamespaces()) library(onestepcox)

nsim <- 10000
n_article <- 1000
deaths <- c("50" = 100, "60" = 80, "70" = 60)

# Simulate the death sequences of nsim trials. Returns a logical matrix with
# one row per trial: TRUE if the k-th death is in group A, NA after the end.
sim_deaths <- function(theta, D, nsim) {
  nA <- rep(100, nsim)
  nB <- rep(100, nsim)
  OB <- rep(0, nsim)
  active <- rep(TRUE, nsim)
  out <- matrix(NA, nsim, 200)
  k <- 0L
  while (any(active)) {
    k <- k + 1L
    isA <- stats::runif(nsim) < nA * theta / (nA * theta + nB)
    out[active, k] <- isA[active]
    nA <- nA - (active & isA)
    nB <- nB - (active & !isA)
    OB <- OB + (active & !isA)
    if (k >= D) active <- active & OB == 0
  }
  out[, seq_len(k), drop = FALSE]
}

# Peto and Pike estimates of one trial: the k-th death occurs at time k and
# the survivors are censored after the last death
estimate_trial <- function(isA) {
  isA <- isA[!is.na(isA)]
  K <- length(isA)
  dA <- sum(isA)
  time <- c(seq_len(K), rep(K + 1, 200 - K))
  event <- c(rep(1, K), rep(0, 200 - K))
  group <- c(ifelse(isA, "A", "B"), rep("A", 100 - dA), rep("B", 100 - (K - dA)))
  c(Peto = EstHR(time, event, group, "B", "Peto", is_sorted = TRUE)[["exp(coef)"]],
    Pike = EstHR(time, event, group, "B", "Pike", is_sorted = TRUE)[["exp(coef)"]])
}

judge <- function(difference, tolerance) {
  ifelse(abs(difference) <= tolerance, "PASS", "FAIL")
}
# The judgement must fail for a value moved by twice the tolerance
stopifnot(judge(2, 1) == "FAIL", judge(0.5, 1) == "PASS")

digitized <- utils::read.csv("inst/reproduce/Berry1991_figure1_digitized.csv",
                             stringsAsFactors = FALSE)
set.seed(1991)
settings <- unique(digitized[, c("censoring", "theta")])
results <- list()
for (i in seq_len(nrow(settings))) {
  cens <- settings$censoring[i]
  theta <- settings$theta[i]
  dmat <- sim_deaths(theta, deaths[[as.character(cens)]], nsim)
  est <- t(apply(dmat, 1L, estimate_trial))
  for (m in c("Peto", "Pike")) {
    row <- digitized[digitized$censoring == cens & digitized$theta == theta &
                       digitized$estimator == m, ]
    x <- est[, m]
    sdx <- stats::sd(x)
    mc_se <- sdx * sqrt(1 / n_article + 1 / nsim)
    tol <- 3 * mc_se + row$reading_error
    results[[length(results) + 1L]] <- data.frame(
      censoring = cens, theta = theta, estimator = m,
      article = row$value, reading_error = row$reading_error,
      mean = mean(x), sd = sdx, mc_se = mc_se, tolerance = tol,
      difference = mean(x) - row$value,
      status = judge(mean(x) - row$value, tol),
      stringsAsFactors = FALSE
    )
  }
}
results <- do.call(rbind, results)
dir.create("inst/reproduce/results", showWarnings = FALSE)
utils::write.csv(results, "inst/reproduce/results/Berry1991_figure1.csv",
                 row.names = FALSE)
print(results[, c("censoring", "theta", "estimator", "article", "mean",
                  "tolerance", "status")], digits = 4)
cat("Berry et al. (1991), Figure 1:", sum(results$status == "PASS"), "PASS,",
    sum(results$status == "FAIL"), "FAIL of", nrow(results), "\n")
