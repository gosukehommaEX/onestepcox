# Reproduction of the analytical results for test data 1 of Appendix E of
# Therneau and Grambsch (2000), Modeling Survival Data: Extending the Cox
# Model (Springer). The data have one binary covariate, so the Cox score,
# the information, and the estimate are available in closed form for the
# Breslow and the Efron handling of ties. With r = exp(beta):
#   Breslow: U = (6 + 3 r - r^2) / ((r + 1) (r + 3)),
#            I = r / (r + 1)^2 + 6 r / (r + 3)^2,
#            estimate r = 3 / 2 + sqrt(33) / 2;
#   Efron:   U = (30 + 23 r - r^3) / ((r + 1) (r + 3) (r + 5)),
#            I = g(r / (r + 1)) + g(r / (r + 3)) + g(r / (r + 5)),
#            g(x) = x - x^2, estimate r = the positive root of
#            r^3 - 23 r - 30 = 0.
# The formulas are those of the tests book1.R and book2.R of the 'survival'
# package, which take them from the appendix; they are to be checked against
# the book (see README.md).
#
# Run from the root of the package after devtools::load_all() (or with the
# installed package). Output:
# inst/reproduce/results/TherneauGrambsch2000.csv.

if (!"onestepcox" %in% loadedNamespaces()) library(onestepcox)

# Test data 1 without the subject whose status is missing
time   <- c(9, 1, 1, 6, 6, 8)
status <- c(1, 1, 0, 1, 1, 0)
x      <- c(0, 1, 1, 1, 0, 0)

truth <- list(
  breslow = list(
    U = function(r) (6 + 3 * r - r^2) / ((r + 1) * (r + 3)),
    I = function(r) r / (r + 1)^2 + 6 * r / (r + 3)^2,
    r_hat = 3 / 2 + sqrt(33) / 2
  ),
  efron = list(
    U = function(r) (30 + 23 * r - r^3) / ((r + 1) * (r + 3) * (r + 5)),
    I = function(r) {
      g <- function(p) p - p^2
      g(r / (r + 1)) + g(r / (r + 3)) + g(r / (r + 5))
    },
    # Largest real root of the cubic by the trigonometric formula
    r_hat = 2 * sqrt(23 / 3) * cos(acos(3 * (-30) / (2 * (-23)) *
                                          sqrt(3 / 23)) / 3)
  )
)
stopifnot(abs(truth$efron$r_hat^3 - 23 * truth$efron$r_hat - 30) < 1e-10)

judge <- function(value, reference, tol) {
  ifelse(abs(value - reference) <= tol * max(1, abs(reference)), "PASS",
         "FAIL")
}
# The judgement must fail for a value moved by more than the tolerance
stopifnot(judge(1 + 1e-6, 1, 1e-9) == "FAIL", judge(1 + 1e-12, 1, 1e-9) == "PASS")

rows <- list()
add <- function(ties, quantity, beta, value, reference, tol) {
  rows[[length(rows) + 1L]] <<- data.frame(
    ties = ties, quantity = quantity, beta = beta, value = value,
    reference = reference, difference = value - reference, tolerance = tol,
    status = judge(value, reference, tol), stringsAsFactors = FALSE
  )
}

rs <- onestepcox:::risk_sets(time, status, x == 1)
for (ti in c("breslow", "efron")) {
  tr <- truth[[ti]]
  dt <- onestepcox:::deriv_table(rs, ti)
  # Score and information at several values of beta
  for (beta in c(0, 0.5, log(tr$r_hat), 2.5)) {
    d <- onestepcox:::score_derivs(beta, dt)
    add(ti, "score U", beta, d[["U"]], tr$U(exp(beta)), 1e-12)
    add(ti, "information I", beta, d[["I"]], tr$I(exp(beta)), 1e-12)
  }
  # Cox estimate and its standard error
  beta_hat <- log(tr$r_hat)
  se_hat <- 1 / sqrt(tr$I(tr$r_hat))
  for (m in c("Cox", "Cox_NR")) {
    est <- EstHR(time, status, x, control = 0, method = m, ties = ti)
    add(ti, paste(m, "coef"), beta_hat, est[["coef"]], beta_hat, 1e-9)
    add(ti, paste(m, "se"), beta_hat, est[["se(coef)"]], se_hat, 1e-9)
  }
}
results <- do.call(rbind, rows)
dir.create("inst/reproduce/results", showWarnings = FALSE)
utils::write.csv(results, "inst/reproduce/results/TherneauGrambsch2000.csv",
                 row.names = FALSE)
print(results, digits = 6)
cat("Therneau and Grambsch (2000), test data 1:",
    sum(results$status == "PASS"), "PASS,", sum(results$status == "FAIL"),
    "FAIL of", nrow(results), "\n")
