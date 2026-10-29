# Settings shared by the scripts of inst/scripts/ (main text) and
# inst/sup_info/ (Supporting Information): the design of the simulations,
# the methods, the scenario tables, the seeds, the number of replicates
# (reduced in a quick test run), and the folders of the results. Each script
# sources this file first. Run the scripts from the root of the package.

if (!file.exists("DESCRIPTION") ||
    !identical(unname(read.dcf("DESCRIPTION", fields = "Package")[1L, 1L]),
               "onestepcox")) {
  stop("Run the scripts from the root of the onestepcox package")
}
if (!isNamespaceLoaded("onestepcox")) library(onestepcox)

# Quick test run, with options(onestepcox.quick = TRUE): 1/100 of the
# replicates (at least 2), fewer hazard ratios for the large-sample limits,
# fewer timing repetitions, and the results in results/quick/
quick <- isTRUE(getOption("onestepcox.quick", FALSE))
n_rep <- function(n) {
  if (quick) max(2L, as.integer(ceiling(n / 100))) else as.integer(n)
}

# Folder of the results of "scripts" (main text) or "sup_info"
# (Supporting Information), with the subfolders data, tables, and figures
results_dir <- function(folder) {
  d <- file.path("inst", folder, "results")
  if (quick) d <- file.path(d, "quick")
  for (s in c("data", "tables", "figures")) {
    dir.create(file.path(d, s), recursive = TRUE, showWarnings = FALSE)
  }
  d
}

# Design of GenSimData(): exponential survival times with control hazard
# log(2) / 18 per month, uniform accrual over 12 months, and the analysis at
# the E-th event with E = 0.4 (n_T + n_C). E = 0.4 n keeps the dropout
# scenarios far from the number of subjects who ever have an event
# (Handover/checks_20261029/check_dropout_capacity.py).
hazard_C <- log(2) / 18
accrual <- 12
event_frac <- 0.4
E_of <- function(n_T, n_C) as.integer(round(event_frac * (n_T + n_C)))
dropout_10 <- -log(0.9) / 12  # dropout of 10 per cent per year

# Methods of the simulations and the applications
methods <- c("Cox", "MH_Halley", "MH_Newton", "Pike_Halley", "Pike_Newton",
             "MH", "Pike", "Peto", "Pike_Chebyshev")

# Reference of the deviations in the simulations: the Cox estimate by
# Newton-Raphson iteration on the risk-set summaries ("Cox_NR"), which stops
# when the step is below 1e-12 and is therefore accurate to rounding. The
# estimate of survival::coxph() ("Cox") stops on the relative change of the
# log partial likelihood, so that it can differ from the root of the score by
# about 1e-12 to 1e-11; it is kept as the usual Cox estimate, and its
# deviation from "Cox_NR" is recorded as a check.
reference <- "Cox_NR"
methods_all <- c(methods, reference)

# Base seeds of the sections. Scenario i of a section uses
# base + 100 (i - 1), as in CompareOC(), and chunk k of a scenario run in
# chunks uses base + 100 (i - 1) + (k - 1). The bases are at least 50,000
# apart, far more than the seeds used in a section (at most 4,200), so that
# the seeds of different sections never coincide.
seeds <- c(fixed = 100001L, local = 200001L, gridA = 300001L,
           gridB = 400001L, gridB_breslow = 450001L, gridC = 500001L,
           timing = 600001L)

# Growth of n under fixed and local alternatives (Figures 2 to 4), 1:1
n_growth <- c(100, 300, 1000, 3000, 1e4, 3e4, 1e5)
nsim_growth <- function(n) {
  if (n <= 3000) 1000L else if (n < 1e5) 200L else 100L
}
theta_fixed <- c(0.1, 0.2, 0.5, 2, 5, 10)
c_local <- c(-1, 1) * 10 * log(2)  # hazard ratio 0.5 and 2 at n = 100

growth_scenarios <- function(type = c("fixed", "local")) {
  type <- match.arg(type)
  if (type == "fixed") {
    # n varies fastest, then the hazard ratio
    g <- expand.grid(n = n_growth, theta = theta_fixed)
    setting <- paste0("HR ", g$theta)
    hr <- g$theta
    tag <- paste0("HR", g$theta)
  } else {
    g <- expand.grid(n = n_growth, c = c_local)
    setting <- ifelse(g$c < 0, "Local, HR 0.5 at n = 100",
                      "Local, HR 2 at n = 100")
    hr <- exp(g$c / sqrt(g$n))
    tag <- ifelse(g$c < 0, "cminus", "cplus")
  }
  data.frame(label = sprintf("%s_%s_n%.0f", type, tag, g$n),
             setting = setting, n = g$n, n_T = g$n / 2, n_C = g$n / 2,
             hazard_C = hazard_C, HR = hr, a = accrual,
             E = E_of(g$n / 2, g$n / 2),
             nsim = vapply(g$n, nsim_growth, integer(1)),
             stringsAsFactors = FALSE)
}

# Grid A: continuous times, no strata, no dropout (Table 1, Table S2)
gridA_scenarios <- function() {
  # The hazard ratio varies fastest, then the allocation, then n_C
  g <- expand.grid(theta = c(0.2, 0.5, 0.7, 1, 2, 5), r = c(1, 2),
                   n_C = c(30, 100, 300))
  data.frame(label = sprintf("A_nC%.0f_r%.0f_HR%s", g$n_C, g$r, g$theta),
             setting = "Continuous", n_T = g$r * g$n_C, n_C = g$n_C,
             r = g$r, hazard_C = hazard_C, HR = g$theta, a = accrual,
             E = E_of(g$r * g$n_C, g$n_C), stringsAsFactors = FALSE)
}

# Grid B: dropout, ties, and strata, with n_T = n_C = 100 (Table 1,
# Table S3). B5 combines the ties of B2 with the strata of B4.
gridB_settings <- list(
  B1 = list(name = "Dropout", dropout_hazard = dropout_10, time_unit = 0,
            strata_prop = NULL, strata_hazard = NULL),
  B2 = list(name = "Ties", dropout_hazard = 0, time_unit = 1,
            strata_prop = NULL, strata_hazard = NULL),
  B3 = list(name = "2 strata", dropout_hazard = 0, time_unit = 0,
            strata_prop = c(0.5, 0.5), strata_hazard = c(1, 3)),
  B4 = list(name = "20 strata", dropout_hazard = 0, time_unit = 0,
            strata_prop = rep(1 / 20, 20),
            strata_hazard = exp(seq(-1, 1, length.out = 20))),
  B5 = list(name = "Ties, 20 strata", dropout_hazard = 0, time_unit = 1,
            strata_prop = rep(1 / 20, 20),
            strata_hazard = exp(seq(-1, 1, length.out = 20)))
)

gridB_scenarios <- function(codes = names(gridB_settings)) {
  # The hazard ratio varies fastest, then the setting
  g <- expand.grid(theta = c(0.2, 0.5, 1, 2, 5), code = codes,
                   stringsAsFactors = FALSE)
  # unname(): the repeated names would otherwise become duplicate row names
  def <- unname(gridB_settings[g$code])
  out <- data.frame(
    label = sprintf("%s_HR%s", g$code, g$theta), code = g$code,
    setting = vapply(def, function(d) d$name, character(1)),
    n_T = 100, n_C = 100, r = 1, hazard_C = hazard_C, HR = g$theta,
    a = accrual, E = E_of(100, 100),
    dropout_hazard = vapply(def, function(d) d$dropout_hazard, numeric(1)),
    time_unit = vapply(def, function(d) d$time_unit, numeric(1)),
    stringsAsFactors = FALSE
  )
  # List columns passed to GenSimData() by CompareOC(); NULL means no strata
  out$strata_prop <- lapply(def, function(d) d$strata_prop)
  out$strata_hazard <- lapply(def, function(d) d$strata_hazard)
  out
}

# Grid C: small samples and extreme hazard ratios (Table S4)
gridC_scenarios <- function() {
  g <- expand.grid(theta = c(0.1, 0.2, 0.5, 2, 5, 10), r = c(1, 2),
                   n = c(30, 60))
  n_C <- g$n / (1 + g$r)
  data.frame(label = sprintf("C_n%.0f_r%.0f_HR%s", g$n, g$r, g$theta),
             setting = "Small sample", n = g$n, n_T = g$n - n_C,
             n_C = n_C, r = g$r, hazard_C = hazard_C, HR = g$theta,
             a = accrual, E = E_of(g$n - n_C, n_C),
             stringsAsFactors = FALSE)
}

# Run CompareOC() on one scenario with all methods and the reference. The
# warning of survival::coxph() about a possibly infinite coefficient
# ("Loglik converged before variable ...; coefficient may be infinite"),
# given when the Cox estimate is not finite (see summarize_methods()) and
# possibly in replicates with a finite estimate (its tolerance toler.inf is
# the square root of eps = 1e-12, the default tol of EstHR(), so 1e-6), is
# counted and not printed; other warnings are printed. The result has the
# estimates, the total computing time of each method, and the number of
# such warnings.
run_compare <- function(sc_row, nsim, seed, ties = "efron") {
  n_warn <- 0L
  res <- withCallingHandlers(
    onestepcox::CompareOC(sc_row, methods_all, nsim = nsim, seed = seed,
                          reference = reference, ties = ties),
    warning = function(w) {
      if (grepl("converged before variable", conditionMessage(w),
                fixed = TRUE)) {
        n_warn <<- n_warn + 1L
        invokeRestart("muffleWarning")
      }
    }
  )
  tp <- res$summary$time_per_sim[match(methods_all, res$summary$method)]
  list(estimates = res$estimates,
       time = stats::setNames(tp * nsim, methods_all), n_warn = n_warn)
}

# Summaries of all methods in one scenario against the reference. In the
# replicates in which the reference has no estimate, all methods are set to
# NA, so that every method is summarized over the same replicates; their
# number is n_undefined. The reference has no estimate when the Cox estimate
# is not finite: when a group has no event (survival::coxph() then gives a
# very large estimate with a warning), or when every event of a group occurs
# at a time when no subject of the other group is at risk in the stratum, so
# that the Mantel-Haenszel anchor is not finite (R = 0 or S = 0). It has none
# either, with a warning, when the Newton-Raphson iteration does not
# converge. id is the column that identifies the replicate.
summarize_methods <- function(e, id, HR_true, time_total, nsim) {
  undefined <- unique(e[[id]][e$method == reference & is.na(e$coef)])
  e[e[[id]] %in% undefined, c("coef", "se", "lower", "upper")] <- NA_real_
  ref <- e[e$method == reference, , drop = FALSE]
  ref <- ref[order(ref[[id]]), , drop = FALSE]
  s <- do.call(rbind, lapply(methods_all, function(m) {
    x <- e[e$method == m, , drop = FALSE]
    x <- x[order(x[[id]]), , drop = FALSE]
    if (!identical(x[[id]], ref[[id]])) {
      stop("the replicates of ", m, " and ", reference, " do not match")
    }
    onestepcox:::summarize_oc(x, ref, HR_true, time_total[[m]] / nsim)
  }))
  data.frame(method = methods_all, nsim = nsim,
             n_undefined = length(undefined), s, stringsAsFactors = FALSE)
}

# Run the scenarios of a grid one by one with 10,000 replicates and keep the
# summaries (the estimates of the grids are not kept). The attribute nsim of
# the result is the number of replicates.
run_grid <- function(sc, base_seed, section, ties = "efron") {
  nsim <- n_rep(10000L)
  out <- vector("list", nrow(sc))
  for (i in seq_len(nrow(sc))) {
    message(sprintf("[%s] scenario %d of %d: %s", section, i, nrow(sc),
                    sc$label[i]))
    r <- run_compare(sc[i, , drop = FALSE], nsim,
                     base_seed + 100L * (i - 1L), ties)
    s <- summarize_methods(r$estimates, "simID", sc$HR[i], r$time, nsim)
    out[[i]] <- data.frame(label = sc$label[i], HR_true = sc$HR[i], s,
                           n_warning = r$n_warn, stringsAsFactors = FALSE)
  }
  s <- do.call(rbind, out)
  rownames(s) <- NULL
  attr(s, "nsim") <- nsim
  s
}

# Run the growth scenarios and keep the estimates of every replicate. For
# n of 30,000 or more, the replicates are run in chunks of 10 to limit the
# memory; the replicate index rep runs over all chunks, and the summaries
# are computed over all replicates.
run_growth <- function(sc, base_seed, section) {
  est <- vector("list", nrow(sc))
  summ <- vector("list", nrow(sc))
  for (i in seq_len(nrow(sc))) {
    nsim <- n_rep(sc$nsim[i])
    chunk <- if (sc$n[i] >= 3e4) min(10L, nsim) else nsim
    starts <- seq(1L, nsim, by = chunk)
    message(sprintf("[%s] scenario %d of %d: %s (%d replicates)", section,
                    i, nrow(sc), sc$label[i], nsim))
    pieces <- vector("list", length(starts))
    time_m <- stats::setNames(numeric(length(methods_all)), methods_all)
    n_warn <- 0L
    for (k in seq_along(starts)) {
      size <- min(chunk, nsim - starts[k] + 1L)
      r <- run_compare(sc[i, , drop = FALSE], size,
                       base_seed + 100L * (i - 1L) + (k - 1L))
      e <- r$estimates
      e$rep <- starts[k] - 1L + e$simID
      e$simID <- NULL
      pieces[[k]] <- e
      time_m <- time_m + r$time
      n_warn <- n_warn + r$n_warn
    }
    e <- do.call(rbind, pieces)
    s <- summarize_methods(e, "rep", sc$HR[i], time_m, nsim)
    summ[[i]] <- data.frame(label = sc$label[i], s, n_warning = n_warn,
                            stringsAsFactors = FALSE)
    est[[i]] <- e
  }
  e <- do.call(rbind, est)
  s <- do.call(rbind, summ)
  rownames(e) <- NULL
  rownames(s) <- NULL
  list(scenarios = sc, estimates = e, summary = s)
}

# One row of the run log: elapsed minutes since t0
log_entry <- function(section, scenarios, replicates, t0) {
  data.frame(section = section, scenarios = scenarios,
             replicates = replicates,
             elapsed_min = round((proc.time()[["elapsed"]] - t0) / 60, 2),
             stringsAsFactors = FALSE)
}

# Write the run log with the versions, and the session information
write_run_log <- function(entries, dir, name) {
  lg <- do.call(rbind, entries)
  lg$quick <- quick
  lg$date <- format(Sys.time(), "%Y-%m-%d %H:%M")
  lg$R <- R.version.string
  lg$onestepcox <- as.character(utils::packageVersion("onestepcox"))
  lg$survival <- as.character(utils::packageVersion("survival"))
  utils::write.csv(lg, file.path(dir, "data", paste0(name, ".csv")),
                   row.names = FALSE)
  writeLines(utils::capture.output(utils::sessionInfo()),
             file.path(dir, "data", paste0(name, "_session.txt")))
  invisible(lg)
}
