# Tables and figures of the main text, made from the results saved by
# data_generate_main.R (no simulation is run here).
#
# Run from the root of the package after devtools::load_all():
#   source("inst/scripts/table_and_figure_manuscript.R")
# With options(onestepcox.quick = TRUE), the results of the quick test run
# in inst/scripts/results/quick/ are used and the output goes there.
#
# Figures: PDF and EPS (fonts embedded by the cairo devices), 119 mm wide,
# Arial. Tables: complete table environments with caption, label, and notes.
# The deviations are taken from the Cox estimate "Cox_NR" (see settings.R).

source(file.path("inst", "scripts", "settings.R"))
source(file.path("inst", "scripts", "figure_helpers.R"))
out <- results_dir("scripts")
data_dir <- file.path(out, "data")
fig_dir <- file.path(out, "figures")
tab_dir <- file.path(out, "tables")

limit <- readRDS(file.path(data_dir, "limit.rds"))
fixed <- readRDS(file.path(data_dir, "fixed.rds"))
local_alt <- readRDS(file.path(data_dir, "local.rds"))
gridA <- readRDS(file.path(data_dir, "gridA.rds"))
gridB <- readRDS(file.path(data_dir, "gridB.rds"))
application <- readRDS(file.path(data_dir, "application.rds"))

# Figure 1: large-sample bias against the hazard ratio --------------------
p1 <- limit_plot(list("1:1" = limit$grid))
save_figure(p1, fig_dir, "Fig1", height_mm = 85)

# Figures 2 and 3: deviation from the Cox estimate against n --------------
gs <- rbind(growth_summary(fixed), growth_summary(local_alt))
main_settings <- c("HR 0.1", "HR 0.5", "HR 2", "HR 10",
                   "Local, HR 0.5 at n = 100", "Local, HR 2 at n = 100")
p2 <- growth_plot(growth_long(gs, c("dev_median", "dev_q99"), main_settings),
                  "|log HR - log HR of Cox|", main_settings,
                  plateau = plateau_data(limit$fixed, main_settings))
save_figure(p2, fig_dir, "Fig2", height_mm = 165)
p3 <- growth_plot(growth_long(gs, c("devse_median", "devse_q99"),
                              main_settings),
                  "|log HR - log HR of Cox| / SE of Cox", main_settings)
save_figure(p3, fig_dir, "Fig3", height_mm = 165)

# Figure 4: error after one step against the error of the anchor ----------
ge <- rbind(growth_estimates(fixed), growth_estimates(local_alt))
key <- paste(ge$label, ge$rep)
cox <- ge[ge$method == reference, ]
wide <- data.frame(label = cox$label, rep = cox$rep, n = cox$n,
                   cox = cox$coef, stringsAsFactors = FALSE)
wide_key <- paste(wide$label, wide$rep)
for (m in c("MH", "Pike", "MH_Newton", "MH_Halley", "Pike_Newton",
            "Pike_Halley")) {
  sel <- ge$method == m
  wide[[m]] <- ge$coef[sel][match(wide_key, key[sel])]
}
# Points drawn: in each scenario, 200 replicates chosen at random and the 20
# with the largest anchor errors (all replicates are in the data)
set.seed(20261029)
keep <- unlist(lapply(split(seq_len(nrow(wide)), wide$label), function(idx) {
  worst <- idx[order(-pmax(abs(wide$MH[idx] - wide$cox[idx]),
                           abs(wide$Pike[idx] - wide$cox[idx]),
                           na.rm = TRUE))][seq_len(min(20L, length(idx)))]
  pick <- idx[sample.int(length(idx), min(200L, length(idx)))]
  union(worst, pick)
}))
w <- wide[sort(keep), ]
steps <- data.frame(anchor = c("MH", "MH", "Pike", "Pike"),
                    step = c("Newton", "Halley", "Newton", "Halley"),
                    method = c("MH_Newton", "MH_Halley", "Pike_Newton",
                               "Pike_Halley"),
                    stringsAsFactors = FALSE)
step_levels <- c("Newton step", "Halley step")
pts <- do.call(rbind, lapply(seq_len(nrow(steps)), function(j) {
  data.frame(anchor = steps$anchor[j],
             step = paste0(steps$step[j], " step"),
             e = abs(w[[steps$anchor[j]]] - w$cox),
             err = abs(w[[steps$method[j]]] - w$cox),
             stringsAsFactors = FALSE)
}))
pts <- pts[is.finite(pts$e) & is.finite(pts$err), ]
pts$e <- pmax(pts$e, floor_value)
pts$err <- pmax(pts$err, floor_value)
pts$step <- factor(pts$step, levels = step_levels)
# Dashed lines K |e|^p (p = 2 for Newton, 3 for Halley), with K the largest
# ratio |one-step - Cox| / |e|^p over all replicates with |e| >= 1e-3. With
# smaller anchor errors the error after one step is near the rounding error,
# and the ratio would reflect the rounding rather than the step.
e_min <- 1e-3
envelope <- do.call(rbind, lapply(c("Newton", "Halley"), function(s) {
  p <- if (s == "Newton") 2 else 3
  e_all <- c(abs(wide$MH - wide$cox), abs(wide$Pike - wide$cox))
  err_all <- c(abs(wide[[paste0("MH_", s)]] - wide$cox),
               abs(wide[[paste0("Pike_", s)]] - wide$cox))
  ok <- is.finite(e_all) & is.finite(err_all) & e_all >= e_min
  if (!any(ok)) return(NULL)
  ratio <- err_all[ok] / e_all[ok]^p
  data.frame(step = factor(paste0(s, " step"), levels = step_levels),
             power = p, K = max(ratio), n_used = sum(ok),
             slope = p, intercept = log10(max(ratio)))
}))
p4 <- ggplot2::ggplot(pts, ggplot2::aes(x = .data$e, y = .data$err,
                                        colour = .data$anchor,
                                        shape = .data$anchor)) +
  ggplot2::geom_point(size = 0.4)
if (!is.null(envelope)) {
  p4 <- p4 + ggplot2::geom_abline(
    data = envelope,
    mapping = ggplot2::aes(slope = .data$slope, intercept = .data$intercept),
    linetype = "dashed", colour = "grey30", linewidth = 0.4
  )
}
p4 <- p4 +
  ggplot2::facet_wrap(~ step, ncol = 2) +
  ggplot2::scale_x_log10() +
  ggplot2::scale_y_log10() +
  onestepcox:::plot_scales(c("MH", "Pike"), c("MH anchor", "Pike anchor"),
                           shape = c(16, 2),
                           aesthetics = c("colour", "shape")) +
  ggplot2::labs(x = "|anchor - Cox| (log HR)",
                y = "|one-step estimate - Cox| (log HR)") +
  onestepcox:::plot_theme(base_size) +
  onestepcox:::plot_legend("bottom", k = 2,
                           aesthetics = c("colour", "shape"))
save_figure(p4, fig_dir, "Fig4", height_mm = 80)
saveRDS(envelope, file.path(data_dir, "fig4_envelope.rds"))

# Table 1: coverage and bias, n_T = n_C = 100 ----------------------------
nsim_grid <- attr(gridA$summary, "nsim")
sA <- merge(gridA$summary,
            gridA$scenarios[, c("label", "setting", "n_C", "r", "E")],
            by = "label", sort = FALSE)
sA <- sA[sA$n_C == 100 & sA$r == 1, ]
sB <- merge(gridB$summary,
            gridB$scenarios[, c("label", "setting", "code", "E")],
            by = "label", sort = FALSE)
sB <- sB[sB$code %in% c("B1", "B2", "B3", "B4"), ]
cols <- c("label", "setting", "HR_true", "method", "coverage", "bias")
s1 <- rbind(sA[, cols], sB[, cols])
t1_settings <- c("Continuous", "Dropout", "Ties", "2 strata", "20 strata")
t1_hr <- c(0.2, 0.5, 1, 2, 5)
t1_methods <- c("Cox", "MH_Halley", "Pike", "Peto")
# Method labels split after the hyphen over two lines, aligned at the bottom
# and on the right, so that the eight columns fit the text width (119 mm).
# With fixed = TRUE the replacement is inserted as it is: "-\\\\" is the
# three characters -, \, and \.
stack_label <- function(x) {
  ifelse(grepl("-", x, fixed = TRUE),
         paste0("\\begin{tabular}[b]{@{}r@{}}",
                sub("-", "-\\\\", x, fixed = TRUE), "\\end{tabular}"),
         x)
}
lines <- character(0)
for (st in t1_settings) {
  for (k in seq_along(t1_hr)) {
    h <- t1_hr[k]
    r <- s1[s1$setting == st & abs(s1$HR_true - h) < 1e-9, ]
    cov <- r$coverage[match(t1_methods, r$method)]
    bias <- r$bias[match(t1_methods, r$method)]
    lines <- c(lines, paste0(
      if (k == 1L) st else "", " & $", h, "$ & ",
      paste(fmt_fixed(100 * cov, 1), collapse = " & "), " & ",
      paste(fmt_fixed(100 * bias, 2), collapse = " & "), " \\\\"
    ))
  }
  if (st != t1_settings[length(t1_settings)]) lines <- c(lines, "\\midrule")
}
E1 <- unique(c(sA$E, sB$E))
tab1 <- c(
  "\\begin{table}",
  paste0("\\caption{Coverage of the 95\\% Wald confidence intervals and ",
         "bias of the estimators of the log hazard ratio, from ",
         fmt_int(nsim_grid), " replicates of each scenario with ",
         "$n_T = n_C = 100$ and the analysis at ", E1,
         " events. The settings are described in Table~S1}"),
  "\\label{tab:operating}",
  "\\centering",
  "\\small",
  "\\setlength{\\tabcolsep}{4pt}",
  "\\begin{tabular}{llrrrrrrrr}",
  "\\toprule",
  paste0(" & & \\multicolumn{4}{c}{Coverage (\\%)} & ",
         "\\multicolumn{4}{c}{Bias $\\times 100$} \\\\"),
  "\\cmidrule(lr){3-6}\\cmidrule(lr){7-10}",
  paste0("Setting & $\\theta$ & ",
         paste(rep(stack_label(label_of(t1_methods)), 2), collapse = " & "),
         " \\\\"),
  "\\midrule",
  lines,
  "\\bottomrule",
  "\\end{tabular}",
  "\\par\\smallskip",
  paste0("{\\footnotesize $\\theta$: true hazard ratio. Cox: Cox estimate ",
         "of \\texttt{survival::coxph()} with the Efron handling of ties and ",
         "its Wald interval; MH-Halley: the proposed estimator. The Monte ",
         "Carlo standard error of a coverage of 95\\% is ",
         formatC(100 * sqrt(0.95 * 0.05 / nsim_grid), format = "f",
                 digits = 2),
         " percentage points.\\par}"),
  "\\end{table}"
)
writeLines(tab1, file.path(tab_dir, "table1_operating.tex"))

# Table 2: applications, Efron handling of ties ----------------------------
t2_methods <- c("Cox", "MH_Halley", "MH_Newton", "Pike_Halley",
                "Pike_Newton", "MH", "Pike", "Peto")
t2_analyses <- c(
  gbsg = "\\texttt{gbsg}",
  gbsg_grade = "\\texttt{gbsg}, stratified by tumor grade",
  pharmacoSmoking = "\\texttt{pharmacoSmoking}",
  pharmacoSmoking_ageGroup2 = paste0("\\texttt{pharmacoSmoking}, ",
                                     "stratified by age group")
)
app <- application[application$ties == "efron", ]
lines <- character(0)
for (a in names(t2_analyses)) {
  r <- app[app$analysis == a, ]
  r <- r[match(t2_methods, r$method), ]
  cox_coef <- r$coef[r$method == "Cox"]
  # The description in two lines, so that it does not widen the table
  head_line <- c(
    paste0("\\multicolumn{3}{l}{", t2_analyses[[a]], "} \\\\"),
    paste0("\\multicolumn{3}{l}{$n_T = ", r$n_T[1L], "$, $n_C = ",
           r$n_C[1L], "$, ", r$events_T[1L] + r$events_C[1L], " events",
           if (r$n_strata[1L] > 1L) paste0(", ", r$n_strata[1L], " strata"),
           "} \\\\")
  )
  body <- paste0(
    label_of(r$method), " & ",
    formatC(r$HR, format = "f", digits = 3), " (",
    formatC(r$lower, format = "f", digits = 3), ", ",
    formatC(r$upper, format = "f", digits = 3), ") & ",
    fmt_sci(abs(r$coef - cox_coef)), " \\\\"
  )
  lines <- c(lines, head_line, body)
  if (a != names(t2_analyses)[length(t2_analyses)]) {
    lines <- c(lines, "\\midrule")
  }
}
tab2 <- c(
  "\\begin{table}",
  paste0("\\caption{Hazard ratio estimates with 95\\% confidence intervals ",
         "in the two applications, with the Efron handling of ties, and the ",
         "absolute difference of the log hazard ratio from the Cox ",
         "estimate}"),
  "\\label{tab:application}",
  "\\centering",
  "\\small",
  "\\begin{tabular}{llr}",
  "\\toprule",
  "Method & HR (95\\% CI) & $|\\hat{\\eta} - \\hat{\\eta}_{\\mathrm{Cox}}|$ \\\\",
  "\\midrule",
  lines,
  "\\bottomrule",
  "\\end{tabular}",
  "\\par\\smallskip",
  paste0("{\\footnotesize Treatment groups: \\texttt{hormon} = ",
         tex_escape(app$treatment[app$analysis == "gbsg"][1L]),
         " in \\texttt{gbsg} and \\texttt{grp} = ",
         tex_escape(app$treatment[app$analysis == "pharmacoSmoking"][1L]),
         " in \\texttt{pharmacoSmoking}. Cox: \\texttt{survival::coxph()}. ",
         "The results with the Breslow handling of ties are in ",
         "Table~S6.\\par}"),
  "\\end{table}"
)
writeLines(tab2, file.path(tab_dir, "table2_application.tex"))

message("Finished. Tables in ", tab_dir, "; figures in ", fig_dir)
