# Validation of LimitHR() against large-sample simulations.
#
# For each hazard ratio (0.2, 0.5, 2, 5, 10) and sample size n (1,000,
# 10,000, 100,000), GenSimData() generates trials with control hazard
# log(2) / 18, 1:1 allocation, uniform accrual over 12 months, and the
# analysis at 0.8 n events. For each estimator m of EstHR(), the difference
# d = (estimate of m) - (Cox estimate by "Cox_NR") converges to the
# large-sample bias of m given by LimitHR(), because the Cox estimate
# converges to the true log hazard ratio. At n = 100,000 (50 trials), the
# check is PASS if |mean(d) - bias| <= 4 SE + 0.01 |bias|, where SE is the
# standard error of mean(d); the term 0.01 |bias| allows for the finite-sample
# part of the deviation. The results for n = 1,000 and 10,000 are reported
# as INFO, to show the convergence.
#
# The one-step corrections from the MH anchor are reported as INFO at every
# n. Their limit is the true log hazard ratio, as for MH itself, so their
# limiting bias is 0 by the same check as MH. Their deviation from the Cox
# estimate at a finite n is a higher-order term of the anchor error e (of
# order e^2 for Newton and e^3 for Halley and Chebyshev); for Newton it has a
# fixed sign and a mean of order 1/n, larger than the standard error of the
# mean, so a test of a zero mean at a fixed n does not apply to them.
#
# Run from the root of the package after devtools::load_all() (or with the
# installed package). Output: inst/validation/results/validate_LimitHR.csv.
# Running time: a few minutes.

if (!"onestepcox" %in% loadedNamespaces()) library(onestepcox)

hazard_C <- log(2) / 18
hrs <- c(0.2, 0.5, 2, 5, 10)
sizes <- data.frame(n = c(1000, 10000, 100000), nsim = c(200, 50, 50),
                    chunk = c(200, 50, 10))
methods <- c("MH_Halley", "MH_Newton", "Pike_Halley", "Pike_Newton",
             "MH_Chebyshev", "Pike_Chebyshev", "MH", "Pike", "Peto")
mh_steps <- c("MH_Halley", "MH_Newton", "MH_Chebyshev")

judge <- function(mean_d, bias, se) {
  if (abs(mean_d - bias) <= 4 * se + 0.01 * abs(bias)) "PASS" else "FAIL"
}
# The judgement must fail for a mean moved by twice the tolerance and pass
# for one moved by half of it (bias 0.1, standard error 1e-4)
tol_check <- 4 * 1e-4 + 0.01 * 0.1
stopifnot(judge(0.1 + 2 * tol_check, 0.1, 1e-4) == "FAIL",
          judge(0.1 + 0.5 * tol_check, 0.1, 1e-4) == "PASS")

# Limits for the design (only the ratios n_T : n_C : E matter)
lim <- LimitHR(hrs, hazard_C, n_T = 1, n_C = 1, a = 12, E = 1.6)

# Differences from the Cox estimate in each simulated trial
sim_diff <- function(hr, n, nsim, chunk, seed) {
  out <- list()
  for (k in seq_len(nsim / chunk)) {
    dat <- GenSimData(nsim = chunk, n_T = n / 2, n_C = n / 2,
                      hazard_C = hazard_C, HR = hr, a = 12, E = 0.8 * n,
                      seed = seed + k)
    out[[k]] <- dat[, {
      cox <- EstHR(time, event, group, "C", "Cox_NR",
                   is_sorted = TRUE)[["coef"]]
      as.list(vapply(methods, function(m) {
        EstHR(time, event, group, "C", m, is_sorted = TRUE)[["coef"]] - cox
      }, numeric(1)))
    }, by = simID]
  }
  do.call(rbind, out)
}

rows <- list()
for (i in seq_along(hrs)) {
  for (j in seq_len(nrow(sizes))) {
    n <- sizes$n[j]
    d <- sim_diff(hrs[i], n, sizes$nsim[j], sizes$chunk[j],
                  seed = 1000 * i + 100 * j)
    for (m in methods) {
      x <- d[[m]]
      bias <- lim$bias[lim$HR == hrs[i] & lim$method == m]
      se <- stats::sd(x) / sqrt(length(x))
      status <- if (n == max(sizes$n) && !m %in% mh_steps) {
        judge(mean(x), bias, se)
      } else {
        "INFO"
      }
      rows[[length(rows) + 1L]] <- data.frame(
        check_id = paste(hrs[i], n, m, sep = ":"), HR = hrs[i], n = n,
        nsim = length(x), method = m, mean_diff = mean(x), se = se,
        limit_bias = bias, tolerance = 4 * se + 0.01 * abs(bias),
        status = status, stringsAsFactors = FALSE
      )
    }
  }
}
results <- do.call(rbind, rows)
dir.create("inst/validation/results", showWarnings = FALSE)
utils::write.csv(results, "inst/validation/results/validate_LimitHR.csv",
                 row.names = FALSE)
print(results[results$n == max(sizes$n),
              c("HR", "method", "mean_diff", "limit_bias", "tolerance",
                "status")], digits = 4)
cat("validate_LimitHR:", paste(names(table(results$status)),
                               table(results$status), collapse = ", "), "\n")
