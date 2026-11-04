# Efficient Mantel-Haenszel (EMH) estimator of Anderson and Bernstein (1985)
# under the fixed alternatives of the main text: Section S10 and Table S9 of
# the Supporting Information.
#
# With the numbers at risk Y_T, Y_C and the numbers of events O_T, O_C at the
# event times, the Cox estimate (Breslow approximation) solves
#   theta = sum(O_T Y_C / (Y_C + theta Y_T)) / sum(O_C Y_T / (Y_C + theta Y_T))
# (their Equation (5)). The Mantel-Haenszel estimate is the right-hand side at
# theta = 1, and the EMH estimate is the right-hand side at the
# Mantel-Haenszel estimate, which Anderson and Bernstein identify with the
# two-step estimator of Begun and Reid (1983) started from the
# Mantel-Haenszel estimate.
#
# The datasets are those of the fixed alternatives of
# inst/scripts/data_generate_main.R (fixed.rds): the seeds of the section
# "fixed" are used again on purpose, so that the datasets are the same.
# Scenario i of growth_scenarios("fixed") uses the seed base + 100 (i - 1),
# and chunk k of a scenario run in chunks of 10 uses base + 100 (i - 1) +
# (k - 1), as in run_growth() of settings.R. The Mantel-Haenszel estimate
# recomputed here is compared with that of fixed.rds in each scenario, which
# confirms that the datasets are the same, and the EMH estimate is compared
# with the reference Cox_NR of fixed.rds (the root of the Cox score by
# Newton-Raphson iteration). With continuous event times, the Breslow and
# Efron approximations coincide.
#
# Run from the root of the package (open onestepcox.Rproj) after
# devtools::load_all(), and after inst/scripts/data_generate_main.R:
#   source("inst/sup_info/data_generate_emh.R")
# For a quick test run, first run options(onestepcox.quick = TRUE); the
# results are then written to inst/sup_info/results/quick/ and need the
# quick results of data_generate_main.R. This script only computes and saves
# the results; table_and_figure_supplement.R makes Table S9 from them.

source(file.path("inst", "scripts", "settings.R"))
out <- results_dir("sup_info")
main <- results_dir("scripts")
fixed <- readRDS(file.path(main, "data", "fixed.rds"))
sc <- growth_scenarios("fixed")
if (!isTRUE(all.equal(sc, fixed$scenarios, check.attributes = FALSE))) {
  stop("the scenarios differ from those of fixed.rds; rerun ",
       "inst/scripts/data_generate_main.R first")
}
run_log <- list()

# Mantel-Haenszel and EMH estimates of the log hazard ratio from a risk table
# of RiskTable() (columns n_T, n_C, O_T, O). Rows with no subject at risk in
# one group add zero to both sums. NA when a sum is not positive.
mh_emh <- function(rs) {
  y <- rs$n_T + rs$n_C
  o_c <- rs$O - rs$O_T
  R <- sum(rs$O_T * rs$n_C / y)
  S <- sum(o_c * rs$n_T / y)
  if (!(R > 0 && S > 0)) return(c(mh = NA_real_, emh = NA_real_))
  theta0 <- R / S
  w <- 1 / (rs$n_C + theta0 * rs$n_T)
  c(mh = log(R / S), emh = log(sum(rs$O_T * rs$n_C * w) /
                                 sum(o_c * rs$n_T * w)))
}

# Estimates of fixed.rds in wide form, one row per scenario and replicate
fe <- fixed$estimates[fixed$estimates$method %in%
                        c("Cox_NR", "MH", "MH_Newton", "MH_Halley"),
                      c("label", "rep", "method", "coef")]
fw <- stats::reshape(fe, idvar = c("label", "rep"), timevar = "method",
                     direction = "wide")

# Risk tables of the first replicate of the scenarios with n up to 1,000,
# kept for the independent check in Python
n_keep_table <- 1000

t0 <- proc.time()[["elapsed"]]
est_rows <- vector("list", nrow(sc))
tables_kept <- list()
for (i in seq_len(nrow(sc))) {
  nsim <- n_rep(sc$nsim[i])
  chunk <- if (sc$n[i] >= 3e4) min(10L, nsim) else nsim
  starts <- seq(1L, nsim, by = chunk)
  message(sprintf("[emh] scenario %d of %d: %s (%d replicates)", i, nrow(sc),
                  sc$label[i], nsim))
  pieces <- vector("list", length(starts))
  for (k in seq_along(starts)) {
    size <- min(chunk, nsim - starts[k] + 1L)
    dat <- onestepcox::GenSimData(
      nsim = size, n_T = as.integer(sc$n_T[i]), n_C = as.integer(sc$n_C[i]),
      hazard_C = sc$hazard_C[i], HR = sc$HR[i], a = sc$a[i],
      E = as.integer(sc$E[i]), dropout_hazard = 0, strata_prop = NULL,
      strata_hazard = NULL, time_unit = 0,
      seed = seeds[["fixed"]] + 100L * (i - 1L) + (k - 1L)
    )
    # Row indices of each replicate; the rows are sorted by simID and time
    idx <- split(seq_len(nrow(dat)), dat$simID)
    est <- lapply(names(idx), function(s) {
      d <- dat[idx[[s]], ]
      rs <- onestepcox::RiskTable(d$time, d$event, d$group, control = "C",
                                  is_sorted = TRUE)
      rep_id <- starts[k] - 1L + as.integer(s)
      if (rep_id == 1L && sc$n[i] <= n_keep_table) {
        tables_kept[[sc$label[i]]] <<- rs
      }
      v <- mh_emh(rs)
      data.frame(label = sc$label[i], rep = rep_id, mh = v[["mh"]],
                 emh = v[["emh"]], stringsAsFactors = FALSE)
    })
    pieces[[k]] <- do.call(rbind, est)
  }
  e <- do.call(rbind, pieces)
  # The recomputed Mantel-Haenszel estimates must agree with fixed.rds up to
  # rounding (and be NA in the same replicates); otherwise the datasets are
  # not the same, and the script stops at the first such scenario
  f <- fw[fw$label == sc$label[i], ]
  f <- f[match(e$rep, f$rep), ]
  if (nrow(e) != sum(fw$label == sc$label[i]) || anyNA(f$rep) ||
      !identical(is.na(e$mh), is.na(f$coef.MH)) ||
      !isTRUE(max(c(0, abs(e$mh - f$coef.MH)), na.rm = TRUE) <= 1e-10)) {
    stop("the Mantel-Haenszel estimates of ", sc$label[i],
         " differ from those of fixed.rds")
  }
  est_rows[[i]] <- e
}
estimates <- do.call(rbind, est_rows)
rownames(estimates) <- NULL
run_log[[length(run_log) + 1L]] <-
  log_entry("emh", nrow(sc), sum(vapply(sc$nsim, n_rep, integer(1))), t0)

m <- merge(estimates, fw, by = c("label", "rep"), sort = FALSE)
if (nrow(m) != nrow(estimates) || nrow(m) != nrow(fw)) {
  stop("the replicates differ from those of fixed.rds")
}
m <- m[order(match(m$label, sc$label), m$rep), ]
rownames(m) <- NULL

# Summaries per scenario. Replicates without a finite Cox estimate are left
# out, as in summarize_methods() of settings.R (n_undefined); a Newton or
# Halley step that is not defined in the other replicates is counted in
# n_newton_na or n_halley_na. mh_check is the largest difference between the
# Mantel-Haenszel estimate recomputed here and that of fixed.rds.
summary_rows <- lapply(seq_len(nrow(sc)), function(i) {
  x <- m[m$label == sc$label[i], ]
  ok <- !is.na(x$coef.Cox_NR)
  d <- list(MH = x$coef.MH - x$coef.Cox_NR,
            EMH = x$emh - x$coef.Cox_NR,
            MH_Newton = x$coef.MH_Newton - x$coef.Cox_NR,
            MH_Halley = x$coef.MH_Halley - x$coef.Cox_NR)
  row <- data.frame(label = sc$label[i], n = sc$n[i], HR = sc$HR[i],
                    nsim = nrow(x), n_undefined = sum(!ok),
                    n_emh_na = sum(ok & is.na(x$emh)),
                    n_newton_na = sum(ok & is.na(x$coef.MH_Newton)),
                    n_halley_na = sum(ok & is.na(x$coef.MH_Halley)),
                    mh_check = max(c(0, abs(x$mh[ok] - x$coef.MH[ok]))),
                    stringsAsFactors = FALSE)
  for (nm in names(d)) {
    a <- abs(d[[nm]][ok])
    a <- a[!is.na(a)]
    row[[paste0(nm, "_median")]] <- stats::median(a)
    row[[paste0(nm, "_q99")]] <- stats::quantile(a, 0.99, names = FALSE)
    row[[paste0(nm, "_max")]] <- max(a)
  }
  row
})
summary_emh <- do.call(rbind, summary_rows)
rownames(summary_emh) <- NULL

# The summaries of MH, MH_Newton, and MH_Halley must reproduce those of
# fixed.rds, which are computed over the same replicates
for (meth in c("MH", "MH_Newton", "MH_Halley")) {
  fs <- fixed$summary[fixed$summary$method == meth, ]
  fs <- fs[match(summary_emh$label, fs$label), ]
  for (st in c("median", "q99", "max")) {
    if (!isTRUE(all.equal(summary_emh[[paste0(meth, "_", st)]],
                          fs[[paste0("dev_", st)]]))) {
      stop("the summary ", st, " of ", meth, " differs from fixed.rds")
    }
  }
}

saveRDS(list(scenarios = sc, summary = summary_emh,
             estimates = m[, c("label", "rep", "mh", "emh", "coef.Cox_NR",
                               "coef.MH", "coef.MH_Newton",
                               "coef.MH_Halley")],
             tables = tables_kept),
        file.path(out, "data", "emh.rds"))
utils::write.csv(summary_emh, file.path(out, "data", "emh_summary.csv"),
                 row.names = FALSE)
# The kept risk tables and the estimates of their replicates (17 significant
# digits), for the independent check in Python
kept <- do.call(rbind, lapply(names(tables_kept), function(l) {
  cbind(data.frame(label = l, stringsAsFactors = FALSE),
        tables_kept[[l]][, c("time", "n_T", "n_C", "O_T", "O")])
}))
kept$time <- sprintf("%.17g", kept$time)
utils::write.csv(kept, file.path(out, "data", "emh_check_tables.csv"),
                 row.names = FALSE)
chk <- m[m$rep == 1L & m$label %in% names(tables_kept),
         c("label", "rep", "mh", "emh", "coef.Cox_NR", "coef.MH",
           "coef.MH_Newton", "coef.MH_Halley")]
for (v in names(chk)[-(1:2)]) chk[[v]] <- sprintf("%.17g", chk[[v]])
utils::write.csv(chk, file.path(out, "data", "emh_check_estimates.csv"),
                 row.names = FALSE)
write_run_log(run_log, out, "run_log_emh")
message("Finished. Results in ", file.path(out, "data", "emh.rds"))
