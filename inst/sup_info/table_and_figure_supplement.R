# Tables and figures of the Supporting Information, made from the results
# saved by inst/sup_info/data_generate_supplement.R,
# inst/sup_info/data_generate_grouped.R,
# inst/sup_info/data_generate_emh.R, and
# inst/scripts/data_generate_main.R (no simulation is run here), and from the
# summary of the reproduction and validation scripts
# (inst/validation/results/summary.csv).
#
# Run from the root of the package after devtools::load_all():
#   source("inst/sup_info/table_and_figure_supplement.R")
# With options(onestepcox.quick = TRUE), the results of the quick test run
# are used and the output goes to inst/sup_info/results/quick/.
#
# The long tables use the LaTeX package longtable, Tables S2 to S4 are set
# on landscape pages (package pdflscape), and the paragraph columns of
# Tables S1 and S7 are ragged right (package array). The deviations are
# taken from the Cox estimate "Cox_NR" (see inst/scripts/settings.R).

source(file.path("inst", "scripts", "settings.R"))
source(file.path("inst", "scripts", "figure_helpers.R"))
out <- results_dir("sup_info")
main <- results_dir("scripts")
fig_dir <- file.path(out, "figures")
tab_dir <- file.path(out, "tables")
read_main <- function(f) readRDS(file.path(main, "data", f))
read_sup <- function(f) readRDS(file.path(out, "data", f))

limit <- read_main("limit.rds")
fixed <- read_main("fixed.rds")
local_alt <- read_main("local.rds")
gridA <- read_main("gridA.rds")
gridB <- read_main("gridB.rds")
application <- read_main("application.rds")
limit_sup <- read_sup("limit_sup.rds")
gridB_breslow <- read_sup("gridB_breslow.rds")
gridC <- read_sup("gridC.rds")
timing <- read_sup("timing.rds")
grouped <- read_sup("grouped.rds")
emh <- read_sup("emh.rds")
validation <- utils::read.csv(file.path("inst", "validation", "results",
                                        "summary.csv"),
                              stringsAsFactors = FALSE)

# formatC() can pad its result with blanks, which trimws() removes
fmt_list <- function(v) {
  paste(trimws(formatC(v, format = "fg", big.mark = ",")), collapse = ", ")
}
fmt_pct <- function(v, digits) fmt_fixed(100 * v, digits)
# Paragraph column of width w with the text ragged right (package array)
p_col <- function(w) paste0(">{\\raggedright\\arraybackslash}p{", w, "}")

# Figure S1: large-sample bias with 2:1 allocation and with dropout -------
ps1 <- limit_plot(list("2:1 allocation" = limit_sup$r2,
                       "1:1, dropout 10% per year" = limit_sup$dropout))
save_figure(ps1, fig_dir, "FigS1", height_mm = 85)

# Figure S2: the growth scenarios with HR 0.2 and 5 -----------------------
gs <- rbind(growth_summary(fixed), growth_summary(local_alt))
s2_settings <- c("HR 0.2", "HR 5")
quantities <- c("|log HR - log HR of Cox|",
                "|log HR - log HR of Cox| / SE of Cox")
d2 <- rbind(
  cbind(growth_long(gs, c("dev_median", "dev_q99"), s2_settings),
        quantity = quantities[1L]),
  cbind(growth_long(gs, c("devse_median", "devse_q99"), s2_settings),
        quantity = quantities[2L])
)
d2$quantity <- factor(d2$quantity, levels = quantities)
pl2 <- plateau_data(limit$fixed, s2_settings)
pl2$quantity <- factor(quantities[1L], levels = quantities)
ps2 <- growth_plot(d2, NULL, s2_settings, plateau = pl2,
                   facet = ggplot2::facet_grid(quantity ~ setting,
                                               scales = "free_y"))
save_figure(ps2, fig_dir, "FigS2", height_mm = 130)

# Table S1: the scenarios -------------------------------------------------
nsim_grid <- attr(gridA$summary, "nsim")
growth_nsim <- unique(gs[, c("n", "nsim")])
growth_nsim <- growth_nsim[order(growth_nsim$n), ]
rep_text <- paste(vapply(unique(growth_nsim$nsim), function(k) {
  nn <- growth_nsim$n[growth_nsim$nsim == k]
  # n outside math mode, where the comma of 1,000 would be spaced as
  # punctuation, and tied to its first value
  paste0(fmt_int(k), " ($n$~=~", fmt_list(nn), ")")
}, character(1)), collapse = "; ")
scA <- gridA$scenarios
scB <- gridB$scenarios
scC <- gridC$scenarios
b_rows <- vapply(names(gridB_settings), function(code) {
  d <- gridB_settings[[code]]
  other <- c(
    if (d$dropout_hazard > 0) "dropout 10\\% per year",
    if (d$time_unit > 0) "times rounded up to months",
    if (!is.null(d$strata_prop)) {
      paste0(length(d$strata_prop), " strata, hazard multipliers ",
             if (length(d$strata_prop) == 2L) {
               paste(d$strata_hazard, collapse = " and ")
             } else {
               "$e^{-1}$ to $e^{1}$"
             })
    }
  )
  paste0(code, " (", d$name, ") & 1:1 & ",
         fmt_list(unique(scB$n_T + scB$n_C)), " & ",
         fmt_list(unique(scB$HR[scB$code == code])), " & ",
         fmt_int(nsim_grid), " & ", paste(other, collapse = "; "), " \\\\")
}, character(1))
tabS1 <- c(
  "\\begin{table}",
  paste0("\\caption{Scenarios of the simulations. In all scenarios the ",
         "survival times are exponential with control hazard ",
         "$\\log(2)/18$ per month, accrual is uniform over ", accrual,
         " months, and the analysis takes place at the $E$-th event with ",
         "$E = ", event_frac, "(n_T + n_C)$, rounded}"),
  "\\label{tab:scenarios}",
  "\\centering",
  "\\footnotesize",
  "\\setlength{\\tabcolsep}{4pt}",
  paste0("\\begin{tabular}{", p_col("2.3cm"), "l", p_col("2.4cm"),
         p_col("2.5cm"), p_col("2.6cm"), p_col("2.8cm"), "}"),
  "\\toprule",
  paste0("Scenarios & $n_T{:}n_C$ & $n_T + n_C$ & True HR & Replicates & ",
         "Dropout, ties, strata \\\\"),
  "\\midrule",
  paste0("Fixed alternatives & 1:1 & ", fmt_list(n_growth), " & ",
         fmt_list(theta_fixed), " & ", rep_text, " & none \\\\"),
  # Two formulas, so that the cell can be broken between them
  paste0("Local alternatives & 1:1 & ", fmt_list(n_growth), " & ",
         "$\\exp(c / \\sqrt{n})$ with $c = \\pm 10 \\log(2)$ & ", rep_text,
         " & none \\\\"),
  paste0("A & ", paste0(unique(scA$r), ":1", collapse = ", "), " & ",
         "$n_C$~=~", fmt_list(unique(scA$n_C)), " & ",
         fmt_list(unique(scA$HR)), " & ", fmt_int(nsim_grid),
         " & none \\\\"),
  b_rows,
  paste0("C & ", paste0(unique(scC$r), ":1", collapse = ", "), " & ",
         fmt_list(unique(scC$n)), " & ", fmt_list(unique(scC$HR)), " & ",
         fmt_int(nsim_grid), " & none \\\\"),
  "\\bottomrule",
  "\\end{tabular}",
  "\\par\\smallskip",
  paste0("{\\footnotesize The deviations from the Cox estimate are taken ",
         "from the root of the score found by Newton-Raphson iteration on ",
         "the risk-set summaries (tolerance $10^{-12}$ on the step). ",
         "Replicates in which the Cox estimate is infinite (for example, ",
         "those in which a group has no event, where ",
         "\\texttt{survival::coxph()} returns a large value with a warning) ",
         "are left out of the summaries of all methods.\\par}"),
  "\\end{table}"
)
writeLines(tabS1, file.path(tab_dir, "tableS1_scenarios.tex"))

# Tables S2 and S3: full summaries of grids A and B ------------------------
oc_methods <- methods  # survival::coxph() and the closed-form estimators
oc_header <- c(
  paste0("Method & Bias & SD & RMSE & $\\overline{\\mathrm{SE}}/\\mathrm{SD}$ ",
         "& Coverage & \\multicolumn{3}{c}{",
         "$|\\hat{\\eta} - \\hat{\\eta}_{\\mathrm{Cox}}|$} & ",
         "Mismatch \\\\"),
  "\\cmidrule(lr){7-9}",
  " & & & & & (\\%) & Median & 99th pct. & Max & (\\%) \\\\"
)
oc_body <- function(s, sc, describe) {
  unlist(lapply(seq_len(nrow(sc)), function(i) {
    r <- s[s$label == sc$label[i], ]
    r <- r[match(oc_methods, r$method), ]
    c(paste0("\\multicolumn{10}{l}{", describe(sc[i, ]), "; replicates ",
             "without a finite Cox estimate: ", fmt_int(r$n_undefined[1L]),
             "} \\\\"),
      paste0(label_of(r$method), " & ", fmt_fixed(r$bias, 4), " & ",
             fmt_fixed(r$sd, 4), " & ", fmt_fixed(r$rmse, 4), " & ",
             fmt_fixed(r$se_mean / r$sd, 3), " & ",
             fmt_pct(r$coverage, 1), " & ", fmt_sci(r$dev_median), " & ",
             fmt_sci(r$dev_q99), " & ", fmt_sci(r$dev_max), " & ",
             fmt_pct(r$mismatch_2dp, 2), " \\\\"))
  }))
}
dev_note <- paste0(
  "$|\\hat{\\eta} - \\hat{\\eta}_{\\mathrm{Cox}}|$: absolute ",
  "difference from the root of the score (for Cox, the difference of ",
  "\\texttt{survival::coxph()} from it)"
)
mismatch_note <- paste0(
  "mismatch: percentage of replicates in which the hazard ratio written ",
  "with two decimals differs from that of the Cox estimate"
)
pike_chebyshev_note <- paste0(
  "Pike-Chebyshev: the step of Chebyshev's method from the Pike anchor, ",
  "which the function \\texttt{coxph\\_fast()} of the \\texttt{R} package ",
  "\\texttt{FastSurvival} computes under the name Pike-Halley estimator."
)
oc_note <- paste0(
  "Bias, SD, and RMSE of the estimate of the log hazard ratio; ",
  "$\\overline{\\mathrm{SE}}/\\mathrm{SD}$: mean standard error divided ",
  "by the standard deviation of the estimates; coverage of the 95\\% Wald ",
  "interval; ", dev_note, "; ", mismatch_note, ". ", pike_chebyshev_note
)
describe_A <- function(x) {
  paste0("$n_T = ", x$n_T, "$, $n_C = ", x$n_C, "$, $E = ", x$E,
         "$, $\\theta = ", x$HR, "$")
}
tabS2 <- longtable_lines(
  "lrrrrrrrrr", 10,
  paste0("Operating characteristics in grid A (continuous times, no ",
         "dropout, no strata), ", fmt_int(nsim_grid),
         " replicates per scenario"),
  "tab:gridA", oc_header,
  oc_body(gridA$summary, gridA$scenarios, describe_A), oc_note,
  landscape = TRUE
)
writeLines(tabS2, file.path(tab_dir, "tableS2_gridA.tex"))

describe_B <- function(ties) {
  function(x) {
    paste0(x$code, " (", x$setting, "), $\\theta = ", x$HR, "$, ",
           ties, " ties")
  }
}
tabS3 <- longtable_lines(
  "lrrrrrrrrr", 10,
  paste0("Operating characteristics in grid B ($n_T = n_C = 100$, $E = ",
         unique(gridB$scenarios$E), "$), ", fmt_int(nsim_grid),
         " replicates per scenario; B2 and B5 also with the Breslow ",
         "handling of ties for all methods"),
  "tab:gridB", oc_header,
  c(oc_body(gridB$summary, gridB$scenarios, describe_B("Efron")),
    oc_body(gridB_breslow$summary, gridB_breslow$scenarios,
            describe_B("Breslow"))),
  oc_note, landscape = TRUE
)
writeLines(tabS3, file.path(tab_dir, "tableS3_gridB.tex"))

# Table S4: small samples and extreme hazard ratios (grid C) --------------
c_body <- unlist(lapply(seq_len(nrow(scC)), function(i) {
  r <- gridC$summary[gridC$summary$label == scC$label[i], ]
  r <- r[match(oc_methods, r$method), ]
  c(paste0("\\multicolumn{8}{l}{$n_T = ", scC$n_T[i], "$, $n_C = ",
           scC$n_C[i], "$, $E = ", scC$E[i], "$, $\\theta = ", scC$HR[i],
           "$; replicates without a finite Cox estimate: ",
           fmt_int(r$n_undefined[1L]), "} \\\\"),
    paste0(label_of(r$method), " & ", fmt_int(r$n_NA - r$n_undefined), " & ",
           fmt_fixed(r$bias, 3), " & ", fmt_pct(r$coverage, 1), " & ",
           fmt_sci(r$dev_median), " & ", fmt_sci(r$dev_q99), " & ",
           fmt_sci(r$dev_max), " & ", fmt_pct(r$mismatch_2dp, 2), " \\\\"))
}))
tabS4 <- longtable_lines(
  "lrrrrrrr", 8,
  paste0("Small samples and extreme hazard ratios (grid C), ",
         fmt_int(nsim_grid), " replicates per scenario"),
  "tab:gridC",
  c(paste0("Method & No estimate & Bias & Coverage & ",
           "\\multicolumn{3}{c}{$|\\hat{\\eta} - \\hat{\\eta}_{\\mathrm{Cox}}|$}",
           " & Mismatch \\\\"),
    "\\cmidrule(lr){5-7}",
    " & & & (\\%) & Median & 99th pct. & Max & (\\%) \\\\"),
  c_body,
  paste0("No estimate: replicates without an estimate of the method among ",
         "those with a finite Cox estimate; bias of the estimate of the log ",
         "hazard ratio; coverage of the 95\\% Wald interval; ", dev_note,
         "; ", mismatch_note, ". ", pike_chebyshev_note),
  landscape = TRUE
)
writeLines(tabS4, file.path(tab_dir, "tableS4_gridC.tex"))

# Table S5: computing time --------------------------------------------------
tm <- timing$times
it <- timing$iterations
n_time <- sort(unique(tm$n))
components <- c(
  aggregation = "Risk-set summaries (aggregation)",
  mh_halley_step = "MH anchor and Halley step",
  newton_1e8 = "Newton-Raphson from MH, tolerance $10^{-8}$",
  newton_1e12 = "Newton-Raphson from MH, tolerance $10^{-12}$",
  EstHR_MH_Halley = "\\texttt{EstHR()}, MH-Halley",
  EstHR_Cox_NR = "\\texttt{EstHR()}, Cox-NR",
  EstHR_Cox = "\\texttt{EstHR()}, Cox",
  coxph_default = "\\texttt{survival::coxph()}, default settings"
)
s5_rows <- vapply(names(components), function(cmp) {
  v <- tm$median_us[tm$expr == cmp][match(n_time, tm$n[tm$expr == cmp])]
  paste0(components[[cmp]], " & ",
         paste(trimws(formatC(v, format = "fg", digits = 3,
                              big.mark = ",")), collapse = " & "),
         " \\\\")
}, character(1))
iter_rows <- c(
  paste0("Iterations, tolerance $10^{-8}$ & ",
         paste(it$iter_1e8[match(n_time, it$n)], collapse = " & "), " \\\\"),
  paste0("Iterations, tolerance $10^{-12}$ & ",
         paste(it$iter_1e12[match(n_time, it$n)], collapse = " & "),
         " \\\\")
)
ev <- tm[match(n_time, tm$n), ]
tabS5 <- c(
  "\\begin{table}",
  paste0("\\caption{Computing time in microseconds (median of ",
         fmt_int(unique(tm$neval)), " repetitions) for one dataset with ",
         "hazard ratio 0.5 and 1:1 allocation}"),
  "\\label{tab:timing}",
  "\\centering",
  "\\footnotesize",
  paste0("\\begin{tabular}{l", strrep("r", length(n_time)), "}"),
  "\\toprule",
  paste0("$n_T + n_C$ & ", paste(fmt_int(n_time), collapse = " & "),
         " \\\\"),
  paste0("Events & ", paste(fmt_int(ev$events), collapse = " & "), " \\\\"),
  paste0("Distinct event times & ",
         paste(fmt_int(ev$event_times), collapse = " & "), " \\\\"),
  "\\midrule",
  s5_rows,
  "\\midrule",
  iter_rows,
  "\\bottomrule",
  "\\end{tabular}",
  "\\par\\smallskip",
  paste0("{\\footnotesize The first row times the computation of the ",
         "risk-set summaries; the next three rows start from summaries ",
         "computed beforehand, so that the closed-form step and the ",
         "iteration are compared on the same summaries. The times of ",
         "\\texttt{EstHR()} include the checks of the input and the ",
         "risk-set summaries. Machine: ",
         tex_escape(paste(timing$machine, collapse = ", ")), ".\\par}"),
  "\\end{table}"
)
writeLines(tabS5, file.path(tab_dir, "tableS5_timing.tex"))

# Table S6: applications, all methods, Efron and Breslow -------------------
s6_analyses <- c(
  gbsg = "\\texttt{gbsg}",
  gbsg_grade = "\\texttt{gbsg}, stratified by tumor grade",
  pharmacoSmoking = "\\texttt{pharmacoSmoking}",
  pharmacoSmoking_ageGroup2 = paste0("\\texttt{pharmacoSmoking}, ",
                                     "stratified by age group")
)
s6_body <- unlist(lapply(names(s6_analyses), function(a) {
  unlist(lapply(c("efron", "breslow"), function(ties) {
    r <- application[application$analysis == a & application$ties == ties, ]
    r <- r[match(methods, r$method), ]
    cox_coef <- r$coef[r$method == "Cox"]
    # The description in two lines, so that it does not widen the table
    c(paste0("\\multicolumn{4}{l}{", s6_analyses[[a]], ", ",
             if (ties == "efron") "Efron" else "Breslow", " ties} \\\\"),
      paste0("\\multicolumn{4}{l}{", r$event_times[1L],
             " distinct event times",
             if (r$n_strata[1L] > 1L) " within strata", "; ",
             r$tied_events[1L], " of ", r$events_T[1L] + r$events_C[1L],
             " events tied} \\\\"),
      paste0(label_of(r$method), " & ",
             formatC(r$HR, format = "f", digits = 4), " (",
             formatC(r$lower, format = "f", digits = 4), ", ",
             formatC(r$upper, format = "f", digits = 4), ") & ",
             formatC(r$se, format = "f", digits = 4), " & ",
             fmt_sci(abs(r$coef - cox_coef)), " \\\\"))
  }))
}))
tabS6 <- longtable_lines(
  "llrr", 4,
  paste0("Hazard ratio estimates in the applications with the Efron and ",
         "Breslow handling of ties"),
  "tab:application_all",
  "Method & HR (95\\% CI) & SE (log HR) & $|\\hat{\\eta} - \\hat{\\eta}_{\\mathrm{Cox}}|$ \\\\",
  s6_body,
  paste0("Cox: \\texttt{survival::coxph()} with the same handling of ties. ",
         "Tied events: events at a time with more than one event.")
)
writeLines(tabS6, file.path(tab_dir, "tableS6_application.tex"))

# Table S7: reproduction and validation -------------------------------------
checks <- c(
  Berry1991_figure1.csv = paste0("Reproduction of Figure 1 of Berry, ",
                                 "Kitchin and Mock (1991): Peto and Pike ",
                                 "estimates in simulated trials"),
  TherneauGrambsch2000.csv = paste0("Reproduction of test data 1 of ",
                                    "Therneau and Grambsch (2000): score, ",
                                    "information, and Cox estimate"),
  validate_python.csv = paste0("All estimators against an independent ",
                               "implementation in Python"),
  validate_LimitHR.csv = paste0("\\texttt{LimitHR()} against simulations ",
                                "with up to 100,000 subjects")
)
v <- validation
v$what <- checks[basename(v$file)]
v$what[is.na(v$what)] <- tex_escape(basename(v$file)[is.na(v$what)])
tabS7 <- c(
  "\\begin{table}",
  paste0("\\caption{Results of the reproduction and validation scripts of ",
         "the package (folders \\texttt{inst/reproduce} and ",
         "\\texttt{inst/validation})}"),
  "\\label{tab:validation}",
  "\\centering",
  "\\footnotesize",
  "\\setlength{\\tabcolsep}{4pt}",
  paste0("\\begin{tabular}{", p_col("5.2cm"), "rrrrr}"),
  "\\toprule",
  "Check & Checks & PASS & EXPLAINED & FAIL & INFO \\\\",
  "\\midrule",
  paste0(v$what, " & ", fmt_int(v$checks), " & ", fmt_int(v$PASS), " & ",
         fmt_int(v$EXPLAINED), " & ", fmt_int(v$FAIL), " & ",
         fmt_int(v$INFO), " \\\\"),
  "\\bottomrule",
  "\\end{tabular}",
  "\\par\\smallskip",
  paste0("{\\footnotesize INFO: results reported without a pass or fail ",
         "criterion; their reasons are given in ",
         "\\texttt{inst/validation/README.md}.\\par}"),
  "\\end{table}"
)
writeLines(tabS7, file.path(tab_dir, "tableS7_validation.tex"))

# Table S8: grouped event times (design of grid B2) -------------------------
gl <- grouped$limits
gs8 <- grouped$simulation
# Differences below 1e-14 in absolute value are at the level of rounding
fmt_diff <- function(v) ifelse(abs(v) < 1e-14, "$< 10^{-14}$", fmt_sci(v))
s8_lim <- paste0(
  "$", gl$theta, "$ & ", fmt_sci(gl$mh - log(gl$theta)), " & ",
  fmt_sci(gl$cox_efron - log(gl$theta)), " & ",
  fmt_sci(gl$cox_breslow - log(gl$theta)), " & ", fmt_diff(gl$diff_efron),
  " & ", fmt_diff(gl$diff_breslow), " \\\\"
)
s8_n <- sort(unique(gs8$n))
s8_sim <- vapply(s8_n, function(n) {
  e <- gs8[gs8$n == n & gs8$ties == "efron", ]
  b <- gs8[gs8$n == n & gs8$ties == "breslow", ]
  paste0(fmt_int(n), " & ", fmt_int(e$nsim), " & ", fmt_sci(e$mh_mean),
         " & ", fmt_sci(e$halley_median_abs), " & ", fmt_sci(b$mh_mean),
         " & ", fmt_sci(b$halley_median_abs), " \\\\")
}, character(1))
# The same datasets for both approximations: count each replicate once
s8_undefined <- sum(gs8$n_undefined[gs8$ties == "efron"])
s8_undefined_text <- if (s8_undefined == 0L) {
  "every replicate had a finite Cox estimate"
} else {
  paste0(fmt_int(s8_undefined), " replicates without a finite Cox estimate ",
         "are left out")
}
tabS8 <- c(
  "\\begin{table}",
  paste0("\\caption{Event times rounded up to whole months in the design ",
         "of grid B2 (1:1 allocation, no dropout, no strata): (a) ",
         "large-sample limits and (b) simulations with the hazard ratio ",
         fmt_list(unique(gs8$HR)), "}"),
  "\\label{tab:grouped}",
  "\\centering",
  "\\footnotesize",
  "\\begin{tabular}{rrrrrr}",
  "\\toprule",
  "\\multicolumn{6}{l}{(a) Large-sample limits} \\\\",
  "\\midrule",
  paste0("$\\theta$ & MH & Cox, Efron & Cox, Breslow & ",
         "\\multicolumn{2}{c}{MH-Halley $-$ Cox} \\\\"),
  "\\cmidrule(lr){5-6}",
  " & & & & Efron & Breslow \\\\",
  "\\midrule",
  s8_lim,
  "\\midrule",
  "\\multicolumn{6}{l}{(b) Simulations} \\\\",
  "\\midrule",
  paste0("$n_T + n_C$ & Replicates & \\multicolumn{2}{c}{Efron} & ",
         "\\multicolumn{2}{c}{Breslow} \\\\"),
  "\\cmidrule(lr){3-4}\\cmidrule(lr){5-6}",
  paste0(" & & Mean $d_{\\mathrm{MH}}$ & Median $|d_{\\mathrm{H}}|$ & ",
         "Mean $d_{\\mathrm{MH}}$ & Median $|d_{\\mathrm{H}}|$ \\\\"),
  "\\midrule",
  s8_sim,
  "\\bottomrule",
  "\\end{tabular}",
  "\\par\\smallskip",
  paste0("{\\footnotesize (a) Limits of the Mantel--Haenszel (MH) and Cox ",
         "estimators of the log hazard ratio minus $\\log \\theta$, ",
         "computed from the limits of the numbers at risk and of events in ",
         "each month, and the limit of the MH-Halley estimator minus that of ",
         "the Cox estimator with the same approximation for ties; ",
         "$< 10^{-14}$: below $10^{-14}$ in absolute value, the level of ",
         "rounding. (b) ",
         "$d_{\\mathrm{MH}}$ and $d_{\\mathrm{H}}$: the Mantel--Haenszel and ",
         "MH-Halley estimates minus the root of the Cox score with the same ",
         "approximation (Newton--Raphson iteration), on the same datasets ",
         "for both approximations; ", s8_undefined_text, ".\\par}"),
  "\\end{table}"
)
writeLines(tabS8, file.path(tab_dir, "tableS8_grouped.tex"))

# Table S9: the efficient Mantel-Haenszel estimator -------------------------
# Medians of |estimate - Cox estimate| under the fixed alternatives, for the
# sample sizes that are powers of ten (data_generate_emh.R)
e9 <- emh$summary
s9_n <- c(100, 1000, 1e4, 1e5)
e9 <- e9[e9$n %in% s9_n, ]
e9 <- e9[order(e9$HR, e9$n), ]
s9_rows <- unlist(lapply(sort(unique(e9$HR)), function(h) {
  x <- e9[e9$HR == h, ]
  first <- c(paste0("$", h, "$"), rep("", nrow(x) - 1L))
  c(paste0(first, " & ", fmt_int(x$n), " & ", fmt_int(x$nsim), " & ",
           fmt_sci(x$MH_median), " & ", fmt_sci(x$EMH_median), " & ",
           fmt_sci(x$MH_Newton_median), " & ", fmt_sci(x$MH_Halley_median),
           " \\\\"),
    if (h != max(e9$HR)) "\\midrule")
}))
s9_undefined <- sum(e9$n_undefined)
s9_undefined_text <- if (s9_undefined == 0L) {
  "every replicate had a finite Cox estimate"
} else {
  paste0(fmt_int(s9_undefined), " replicates without a finite Cox estimate ",
         "are left out")
}
s9_step_na <- sum(e9$n_emh_na + e9$n_newton_na + e9$n_halley_na)
s9_step_text <- if (s9_step_na == 0L) {
  "and all estimators were defined in the other replicates"
} else {
  paste0("and an estimator was not defined in ", fmt_int(s9_step_na),
         " cases of the other replicates, which are left out of its median")
}
tabS9 <- c(
  "\\begin{table}",
  paste0("\\caption{Median distance from the Cox estimate of the log hazard ",
         "ratio under the fixed alternatives: the Mantel--Haenszel ",
         "estimator, the efficient Mantel--Haenszel estimator (EMH), and ",
         "one Newton and one Halley step from the Mantel--Haenszel ",
         "estimate}"),
  "\\label{tab:emh}",
  "\\centering",
  "\\footnotesize",
  "\\begin{tabular}{rrrrrrr}",
  "\\toprule",
  paste0("$\\theta$ & $n_T + n_C$ & Replicates & ",
         "\\multicolumn{4}{c}{Median of ",
         "$|\\hat{\\eta} - \\hat{\\eta}_{\\mathrm{Cox}}|$} \\\\"),
  "\\cmidrule(lr){4-7}",
  " & & & MH & EMH & MH-Newton & MH-Halley \\\\",
  "\\midrule",
  s9_rows,
  "\\bottomrule",
  "\\end{tabular}",
  "\\par\\smallskip",
  paste0("{\\footnotesize The datasets of the fixed alternatives of ",
         "Figure~2 of the main article (Table~\\ref{tab:scenarios}), with ",
         "$n_T = n_C$ and continuous event times, for the sample sizes that ",
         "are powers of ten. MH: Mantel--Haenszel estimator. EMH: the ",
         "efficient Mantel--Haenszel estimator of Anderson and Bernstein ",
         "(1985), the solution of the Cox estimating equation with its ",
         "weights evaluated at the Mantel--Haenszel estimate. MH-Newton and ",
         "MH-Halley: one Newton and one Halley step from the ",
         "Mantel--Haenszel estimate (Section~2.3 of the main article). ",
         "$\\hat{\\eta}_{\\mathrm{Cox}}$: the root of the Cox score by ",
         "Newton--Raphson iteration on the same numbers at risk and of ",
         "events; ", s9_undefined_text, ", ", s9_step_text, ".\\par}"),
  "\\end{table}"
)
writeLines(tabS9, file.path(tab_dir, "tableS9_emh.tex"))

message("Finished. Tables in ", tab_dir, "; figures in ", fig_dir)
