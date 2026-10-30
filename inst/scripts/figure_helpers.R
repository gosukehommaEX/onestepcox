# Helpers shared by table_and_figure_manuscript.R and
# inst/sup_info/table_and_figure_supplement.R: saving figures, formatting
# numbers for LaTeX, and the figures of the growth scenarios. Sourced after
# settings.R.

base_size <- 8
floor_value <- 1e-16  # smaller values are drawn at this value on log axes
label_of <- function(m) onestepcox:::plot_labels(m)

# Save a figure as PDF and EPS, 119 mm wide (the width of the small-format
# journals of Springer), with the fonts embedded by the cairo devices
save_figure <- function(p, dir, name, height_mm, width_mm = 119) {
  w <- width_mm / 25.4
  h <- height_mm / 25.4
  grDevices::cairo_pdf(file.path(dir, paste0(name, ".pdf")), width = w,
                       height = h, family = "Arial")
  print(p)
  grDevices::dev.off()
  grDevices::cairo_ps(file.path(dir, paste0(name, ".eps")), width = w,
                      height = h, family = "Arial", onefile = FALSE)
  print(p)
  grDevices::dev.off()
  invisible(name)
}

# Numbers in LaTeX math mode. Rounding first and adding 0 turns a negative
# zero into a zero, so that "-0.00" is not written.
fmt_fixed <- function(v, digits) {
  x <- round(v, digits) + 0
  ifelse(is.finite(x), paste0("$", formatC(x, format = "f", digits = digits),
                              "$"), "NA")
}
fmt_sci <- function(v) {
  vapply(v, function(z) {
    if (!is.finite(z)) return("NA")
    # Exactly zero, as for the reference itself
    if (abs(z) < 1e-300) return("$0$")
    e <- floor(log10(abs(z)))
    m <- z / 10^e
    if (abs(round(m, 1)) >= 10) {
      m <- m / 10
      e <- e + 1
    }
    # Without the factor 10^0
    if (as.integer(e) == 0L) return(sprintf("$%.1f$", m))
    sprintf("$%.1f \\times 10^{%d}$", m, as.integer(e))
  }, character(1))
}
fmt_int <- function(v) formatC(v, format = "d", big.mark = ",")

# Escape the characters with a special meaning in LaTeX, in strings taken
# from the system or from file names. With fixed = TRUE the replacement is
# inserted as it is: "\\_" is the two characters \ and _.
tex_escape <- function(x) {
  x <- gsub("\\", "\\textbackslash ", x, fixed = TRUE)
  for (ch in c("&", "%", "$", "#", "_", "{", "}")) {
    x <- gsub(ch, paste0("\\", ch), x, fixed = TRUE)
  }
  x <- gsub("~", "\\textasciitilde ", x, fixed = TRUE)
  gsub("^", "\\textasciicircum ", x, fixed = TRUE)
}

# Large-sample limits for Figure 1 and Figure S1: absolute biases below
# 1e-10 (near the hazard ratio 1, where they are of the order of the
# accuracy of the numerical integration) are set to NA, so that they are not
# drawn
limit_for_plot <- function(lim) {
  lim$bias[abs(lim$bias) < 1e-10] <- NA_real_
  lim
}

# Figure of the absolute large-sample bias against the hazard ratio for one
# design, or one panel per design. lims is a named list of LimitHR() results.
# The lines are broken where the bias is below 1e-10 (see limit_for_plot()).
limit_methods <- c("Pike", "Peto", "Pike_Newton", "Pike_Halley")
limit_plot <- function(lims) {
  d <- do.call(rbind, lapply(names(lims), function(nm) {
    x <- onestepcox:::plot_data_LimitHR(limit_for_plot(lims[[nm]]),
                                        limit_methods, "abs")
    x$design <- nm
    x
  }))
  d$method <- factor(d$method, levels = limit_methods)
  d$design <- factor(d$design, levels = names(lims))
  p <- ggplot2::ggplot(d, ggplot2::aes(x = .data$HR, y = .data$value,
                                       colour = .data$method,
                                       shape = .data$method,
                                       linetype = .data$method,
                                       group = .data$method)) +
    ggplot2::geom_line(linewidth = 0.5, na.rm = TRUE) +
    ggplot2::geom_point(size = 0.6, na.rm = TRUE) +
    ggplot2::scale_x_log10(breaks = c(0.1, 0.2, 0.5, 1, 2, 5, 10)) +
    ggplot2::scale_y_log10() +
    onestepcox:::plot_scales(limit_methods, label_of(limit_methods),
                             shape = c(16, 17, 15, 18),
                             linetype = c("solid", "dashed", "dotted",
                                          "dotdash")) +
    ggplot2::labs(x = "Hazard ratio",
                  y = "|Large-sample bias| of the log hazard ratio") +
    onestepcox:::plot_theme(base_size) +
    # One row of four methods (the default would be three and one)
    onestepcox:::plot_legend("bottom", legend_ncol = length(limit_methods),
                             k = length(limit_methods))
  if (length(lims) > 1L) p <- p + ggplot2::facet_wrap(~ design, ncol = 2)
  p
}

# Summaries and estimates of the growth scenarios with their setting and n
growth_summary <- function(g) {
  sc <- g$scenarios[, c("label", "setting", "n", "HR")]
  merge(g$summary, sc, by = "label", sort = FALSE)
}
growth_estimates <- function(g) {
  merge(g$estimates, g$scenarios[, c("label", "setting", "n")],
        by = "label", sort = FALSE)
}

growth_methods <- c("MH_Halley", "MH_Newton", "Pike_Halley", "Pike_Newton",
                    "MH", "Pike")
growth_shapes <- c(16, 17, 15, 18, 1, 2)
stat_names <- c("Median", "99th percentile")

# Long format with one row per setting, n, method, and statistic
growth_long <- function(s, cols, settings) {
  s <- s[s$setting %in% settings & s$method %in% growth_methods, ]
  do.call(rbind, lapply(seq_along(cols), function(j) {
    data.frame(setting = s$setting, n = s$n, method = s$method,
               statistic = stat_names[j], value = s[[cols[j]]],
               stringsAsFactors = FALSE)
  }))
}

# Large-sample limits of |Pike-anchored estimator - Cox| at the fixed hazard
# ratios, for the dotted lines of Figure 2 (the Mantel-Haenszel-anchored
# estimators converge to the Cox limit)
plateau_data <- function(limit_fixed, settings) {
  lim <- as.data.frame(limit_fixed)
  pl <- lim[lim$method %in% c("Pike", "Pike_Newton", "Pike_Halley"), ]
  cox <- lim[lim$method == "Cox", ]
  pl$value <- pmax(abs(pl$limit - cox$limit[match(pl$HR, cox$HR)]),
                   floor_value)
  pl$setting <- paste0("HR ", pl$HR)
  pl[pl$setting %in% settings, c("setting", "method", "value")]
}

# Figure of a statistic of the deviation against n, one panel per setting
# (facet_wrap), or with the layer facet given (for example facet_grid())
growth_plot <- function(d, ylab, settings, plateau = NULL, facet = NULL) {
  d$value <- pmax(d$value, floor_value)
  d$method <- factor(d$method, levels = growth_methods)
  d$setting <- factor(d$setting, levels = settings)
  d$statistic <- factor(d$statistic, levels = stat_names)
  p <- ggplot2::ggplot(d, ggplot2::aes(
    x = .data$n, y = .data$value, colour = .data$method,
    shape = .data$method, linetype = .data$statistic,
    group = interaction(.data$method, .data$statistic)
  ))
  if (!is.null(plateau) && nrow(plateau) > 0L) {
    plateau$method <- factor(plateau$method, levels = growth_methods)
    plateau$setting <- factor(plateau$setting, levels = settings)
    p <- p + ggplot2::geom_hline(
      data = plateau,
      mapping = ggplot2::aes(yintercept = .data$value,
                             colour = .data$method),
      linetype = "dotted", linewidth = 0.3, inherit.aes = FALSE,
      show.legend = FALSE
    )
  }
  if (is.null(facet)) facet <- ggplot2::facet_wrap(~ setting, ncol = 2)
  p +
    ggplot2::geom_line(linewidth = 0.4) +
    ggplot2::geom_point(size = 1) +
    facet +
    # Twice the default expansion, so that the label 100,000 of one panel
    # does not run into the label 100 of the next
    ggplot2::scale_x_log10(breaks = c(1e2, 1e3, 1e4, 1e5),
                           labels = c("100", "1,000", "10,000",
                                      "100,000"),
                           expand = ggplot2::expansion(mult = 0.1)) +
    ggplot2::scale_y_log10() +
    onestepcox:::plot_scales(growth_methods, label_of(growth_methods),
                             shape = growth_shapes,
                             aesthetics = c("colour", "shape")) +
    ggplot2::scale_linetype_manual(values = c("solid", "dashed"),
                                   breaks = stat_names, name = NULL) +
    ggplot2::labs(x = "Total sample size n", y = ylab) +
    onestepcox:::plot_theme(base_size) +
    onestepcox:::plot_legend("bottom", k = length(growth_methods),
                             aesthetics = c("colour", "shape")) +
    ggplot2::theme(legend.box = "vertical")
}

# A long table for the Supporting Information (package longtable): the
# header is repeated on every page. body is a character vector of rows. The
# table and the note are set in \footnotesize; with landscape = TRUE, on
# landscape pages (package pdflscape), for the tables wider than the text.
longtable_lines <- function(colspec, ncol, caption, label, header, body,
                            note = NULL, landscape = FALSE) {
  out <- c("\\begingroup",
           "\\footnotesize",
           paste0("\\begin{longtable}{", colspec, "}"),
           paste0("\\caption{", caption, "}\\label{", label, "}\\\\"),
           "\\toprule", header, "\\midrule", "\\endfirsthead",
           paste0("\\multicolumn{", ncol, "}{l}{Table~\\thetable{} ",
                  "(continued)}\\\\"),
           "\\toprule", header, "\\midrule", "\\endhead",
           "\\bottomrule", "\\endfoot",
           body,
           "\\end{longtable}",
           if (!is.null(note)) paste0(note, "\\par"),
           "\\endgroup")
  if (landscape) out <- c("\\begin{landscape}", out, "\\end{landscape}")
  out
}
