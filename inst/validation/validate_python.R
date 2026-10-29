# Validation of EstHR() against an independent Python implementation.
#
# inst/validation/python/make_reference.py generates 203 datasets
# (data/python_datasets.csv; settings in data/python_datasets_info.csv) and
# computes every estimator with inst/validation/python/reference.py
# (data/python_results.csv), which counts the risk sets by brute force and
# finds the Cox estimate by Brent's method. This script computes the same
# estimates with EstHR() for the Efron and the Breslow handling of ties and
# compares coef and se:
#   - methods other than "Cox": PASS if both are missing, or if
#     |R - Python| <= 1e-10 * max(1, |Python|) (1e-9 for "Cox_NR");
#   - "Cox" (survival::coxph()) against the Python Cox estimate with
#     tolerance 1e-8; INFO when the Python Cox estimate does not exist
#     (coxph() then returns a large finite value with a warning).
# Stratified estimates are used for the datasets with more than one stratum.
# A failure listed in explained.csv is reported as EXPLAINED.
#
# Run from the root of the package after devtools::load_all() (or with the
# installed package). Output: inst/validation/results/validate_python.csv.
# Running time: about a minute.

if (!"onestepcox" %in% loadedNamespaces()) library(onestepcox)

dat <- utils::read.csv("inst/validation/data/python_datasets.csv",
                       stringsAsFactors = FALSE)
info <- utils::read.csv("inst/validation/data/python_datasets_info.csv",
                        stringsAsFactors = FALSE)
ref <- utils::read.csv("inst/validation/data/python_results.csv",
                       stringsAsFactors = FALSE)
explained <- utils::read.csv("inst/validation/explained.csv",
                             stringsAsFactors = FALSE)

judge <- function(x, y, tol) {
  if (is.na(x) && is.na(y)) return("PASS")
  if (is.na(x) || is.na(y)) return("FAIL")
  if (abs(x - y) <= tol * max(1, abs(y))) "PASS" else "FAIL"
}
# The judgement must fail for a value moved beyond the tolerance and for a
# missing value on one side only
stopifnot(judge(1 + 1e-8, 1, 1e-10) == "FAIL", judge(NA, 1, 1e-10) == "FAIL",
          judge(1 + 1e-12, 1, 1e-10) == "PASS", judge(NA, NA, 1e-10) == "PASS")

methods <- c("MH_Halley", "MH_Newton", "Pike_Halley", "Pike_Newton",
             "MH_Chebyshev", "Pike_Chebyshev", "MH", "Pike", "Peto", "Cox_NR")
rows <- vector("list", nrow(info) * 2L * (length(methods) + 1L) * 2L)
idx <- 0L
for (i in info$dataset) {
  d <- dat[dat$dataset == i, ]
  strata <- if (info$n_strata[info$dataset == i] > 1) d$stratum else NULL
  for (ti in c("efron", "breslow")) {
    for (m in c(methods, "Cox")) {
      warn <- character(0)
      est <- withCallingHandlers(
        EstHR(d$time, d$event, d$group, "C", m, ties = ti, strata = strata),
        warning = function(w) {
          warn <<- c(warn, conditionMessage(w))
          invokeRestart("muffleWarning")
        }
      )
      r <- ref[ref$dataset == i & ref$ties == ti &
                 ref$method == (if (m == "Cox") "Cox_NR" else m), ]
      stopifnot(nrow(r) == 1L)
      tol <- switch(m, Cox = 1e-8, Cox_NR = 1e-9, 1e-10)
      for (q in c("coef", "se")) {
        x <- unname(est[[if (q == "coef") "coef" else "se(coef)"]])
        y <- r[[q]]
        status <- judge(x, y, tol)
        if (m == "Cox" && is.na(y)) status <- "INFO"
        check_id <- paste(i, ti, m, q, sep = ":")
        if (status == "FAIL" && check_id %in% explained$check_id) {
          status <- "EXPLAINED"
        }
        idx <- idx + 1L
        rows[[idx]] <- data.frame(
          check_id = check_id, dataset = i, ties = ti, method = m,
          quantity = q, R = x, Python = y, tolerance = tol, status = status,
          note = paste(unique(warn), collapse = "; "),
          stringsAsFactors = FALSE
        )
      }
    }
  }
}
results <- do.call(rbind, rows[seq_len(idx)])

# A deliberately perturbed reference value must give a failure
x0 <- results$Python[which(!is.na(results$Python))[1]]
stopifnot(judge(x0 * (1 + 1e-8) + 1e-8, x0, 1e-10) == "FAIL")

dir.create("inst/validation/results", showWarnings = FALSE)
utils::write.csv(results, "inst/validation/results/validate_python.csv",
                 row.names = FALSE)
print(table(results$method, results$status))
cat("validate_python:", paste(names(table(results$status)),
                              table(results$status), collapse = ", "), "\n")
