# Check the numbers of the article against the results of the scripts.
#
# The ledger inst/scripts/number_locations.csv (UTF-8) has one row per number
# written in the main text, the Supporting Information, the cover letter, or
# the title page,
# with the columns
#   id           unique identifier
#   file         main.tex, supporting_information.tex, cover_letter.tex, or
#                title_page.tex
#   anchor       text of the .tex around the number, found exactly once
#   text_value   the number as written in the text
#   source_file  .rds or .csv file of the results, relative to the root of
#                the package (required for type computed)
#   source_expr  one R expression that gives the number: for type computed
#                it uses x, the content of source_file; for type input it
#                uses an object of settings.R (such as n_rep() or event_frac)
#   rounding     rule that writes the number as in the text (see below)
#   type         computed, input, or literature:<bibkey>
# For each row of type computed or input, the expression is evaluated, the
# value is written with the rounding rule, and the result is compared with
# text_value. Rows of type literature are listed without a check here; the
# citation is checked by Manuscript/tools/check_numbers.py.
#
# Rounding rules (the same rules are implemented in check_numbers.py):
#   fixed:d   d decimals, with a comma every three digits   (0.083, 12,345.7)
#   pct:d     100 times the value with d decimals            (95.1)
#   int       integer with a comma every three digits        (10,438)
#   signif:d  d (>= 1) significant digits, no exponent and no trailing
#             zeros                                          (0.083, 1,200)
#   approx:d  the same as signif:d, for "about x times"     (30)
#   order     the power of ten nearest on the log scale      (10^{-6})
#   sci:d     mantissa with d decimals times a power of ten
#                                                          (2.5 \times 10^{-6})
#   text      a character value, compared as it is           (0.1.0)
# The values are rounded with round() before they are written, as the tables
# of the article are (fmt_fixed() in figure_helpers.R), or with signif() for
# the rules signif and approx. A negative number is written -0.012 (in the text $-0.012$). The comparison
# ignores white space and reads "{,}" as ",". The rules are checked first
# against inst/scripts/number_rules_test.csv.
#
# Run from the root of the package, in a fresh R session, after the data
# scripts, and again after any change of the results, the settings, or the
# ledger:
#   devtools::load_all()
#   source("inst/scripts/number_check.R")
# (settings.R calls library(onestepcox) unless the package is loaded.)
# The result is written to inst/scripts/results/number_check.csv, with the
# MD5 of each source file (of settings.R for an input without a source file),
# so that check_numbers.py can tell whether the results or the settings
# changed after this script was run.

if (!file.exists("DESCRIPTION") ||
    !identical(unname(read.dcf("DESCRIPTION", fields = "Package")[1L, 1L]),
               "onestepcox")) {
  stop("Run the script from the root of the onestepcox package")
}
out_dir <- file.path("inst", "scripts", "results")
out_file <- file.path(out_dir, "number_check.csv")
# An old result must not outlive a run that stops
unlink(out_file)

# Settings of the full run (not the quick test run), for type input
settings_env <- new.env()
local({
  old <- options(onestepcox.quick = FALSE)
  on.exit(options(old))
  sys.source(file.path("inst", "scripts", "settings.R"), envir = settings_env)
})

normalize_value <- function(s) {
  gsub("{,}", ",", gsub("[[:space:]]+", "", s), fixed = TRUE)
}

# Write the value v with the rounding rule
format_rule <- function(v, rule) {
  parts <- strsplit(rule, ":", fixed = TRUE)[[1L]]
  kind <- parts[1L]
  d <- if (length(parts) == 2L) suppressWarnings(as.integer(parts[2L])) else NA
  valid <- if (kind %in% c("fixed", "pct", "sci")) {
    length(parts) == 2L && !is.na(d) && d >= 0L && d <= 7L
  } else if (kind %in% c("signif", "approx")) {
    length(parts) == 2L && !is.na(d) && d >= 1L && d <= 7L
  } else {
    kind %in% c("int", "order", "text") && length(parts) == 1L
  }
  if (!valid) stop("unknown rounding rule '", rule, "'")
  if (kind == "text") {
    if (!is.character(v)) stop("the rule text needs a character value")
    return(v)
  }
  if (!is.numeric(v) || !is.finite(v)) {
    stop("the value is not a finite number")
  }
  switch(
    kind,
    fixed = formatC(round(v, d) + 0, format = "f", digits = d,
                    big.mark = ",", decimal.mark = "."),
    pct = formatC(round(100 * v, d) + 0, format = "f", digits = d,
                  big.mark = ",", decimal.mark = "."),
    # format "f" with 0 digits rather than "d", which fails beyond the range
    # of integers
    int = formatC(round(v) + 0, format = "f", digits = 0L, big.mark = ",",
                  decimal.mark = "."),
    signif = ,
    approx = format(signif(v, d) + 0, digits = 15L, scientific = FALSE,
                    trim = TRUE, big.mark = ",", decimal.mark = "."),
    order = paste0("10^{", round(log10(abs(v))), "}"),
    sci = {
      if (v == 0) stop("sci needs a non-zero value")
      e <- floor(log10(abs(v)))
      m <- v / 10^e
      if (abs(round(m, d)) >= 10) {
        m <- m / 10
        e <- e + 1
      }
      paste0(formatC(m, format = "f", digits = d, decimal.mark = "."),
             " \\times 10^{", e, "}")
    }
  )
}

# Check of the rules against the expected strings
rule_tests <- utils::read.csv(file.path("inst", "scripts",
                                        "number_rules_test.csv"),
                              colClasses = "character",
                              na.strings = character(0),
                              fileEncoding = "UTF-8-BOM")
rule_ok <- vapply(seq_len(nrow(rule_tests)), function(i) {
  v <- rule_tests$value[i]
  if (rule_tests$rule[i] != "text") v <- as.numeric(v)
  identical(normalize_value(format_rule(v, rule_tests$rule[i])),
            normalize_value(rule_tests$expected[i]))
}, logical(1))
if (!all(rule_ok)) {
  stop("the rounding rules fail on rows ",
       paste(which(!rule_ok), collapse = ", "),
       " of inst/scripts/number_rules_test.csv")
}

ledger_path <- file.path("inst", "scripts", "number_locations.csv")
ledger <- utils::read.csv(ledger_path, colClasses = "character",
                          na.strings = character(0),
                          fileEncoding = "UTF-8-BOM")
columns <- c("id", "file", "anchor", "text_value", "source_file",
             "source_expr", "rounding", "type")
if (!identical(names(ledger), columns)) {
  stop("the columns of ", ledger_path, " must be ",
       paste(columns, collapse = ", "))
}
if (anyDuplicated(ledger$id) > 0L) {
  stop("duplicated id in ", ledger_path, ": ",
       paste(unique(ledger$id[duplicated(ledger$id)]), collapse = ", "))
}

cache <- list()
read_source <- function(f) {
  if (grepl("/quick/", gsub("\\\\", "/", f))) {
    stop("the results of the quick test run are not used")
  }
  if (is.null(cache[[f]])) {
    if (!file.exists(f)) stop("file not found: ", f)
    cache[[f]] <<- if (grepl("\\.rds$", f)) {
      readRDS(f)
    } else if (grepl("\\.csv$", f)) {
      utils::read.csv(f, stringsAsFactors = FALSE)
    } else {
      stop("source_file must be an .rds or .csv file")
    }
  }
  cache[[f]]
}

check_row <- function(r) {
  out <- c(formatted = "", recomputed = "", result = "", message = "")
  if (startsWith(r$type, "literature:")) {
    out[["result"]] <- "LITERATURE"
    return(out)
  }
  if (!(r$type %in% c("computed", "input"))) {
    out[["result"]] <- "FAIL"
    out[["message"]] <- "type must be computed, input, or literature:<bibkey>"
    return(out)
  }
  res <- tryCatch({
    e <- parse(text = r$source_expr)
    if (length(e) != 1L) stop("source_expr must be one expression")
    nm <- all.names(e)
    if (r$type == "computed" &&
        (!nzchar(r$source_file) || !("x" %in% nm))) {
      stop("a computed value must be read from source_file through x")
    }
    if (r$type == "input" && !any(nm %in% ls(settings_env))) {
      stop("an input value must come from an object of settings.R")
    }
    env <- new.env(parent = settings_env)
    # x is always defined, so that no object x of the workspace is used
    assign("x", if (nzchar(r$source_file)) read_source(r$source_file),
           envir = env)
    v <- eval(e[[1L]], envir = env)
    if (length(v) != 1L) stop("the expression must give one value")
    v <- unname(v)
    if (is.factor(v)) v <- as.character(v)
    v <- as.vector(v)
    list(v = v, formatted = format_rule(v, r$rounding))
  }, error = function(e) e)
  if (inherits(res, "error")) {
    out[["result"]] <- "FAIL"
    out[["message"]] <- conditionMessage(res)
    return(out)
  }
  out[["formatted"]] <- res$formatted
  out[["recomputed"]] <- if (is.numeric(res$v)) {
    sprintf("%.17g", res$v)
  } else {
    as.character(res$v)
  }
  same <- identical(normalize_value(res$formatted),
                    normalize_value(r$text_value))
  out[["result"]] <- if (same) "PASS" else "FAIL"
  if (!same) out[["message"]] <- "the text value differs"
  out
}

checked <- lapply(seq_len(nrow(ledger)), function(i) check_row(ledger[i, ]))
# MD5 of the source file of each row; for an input without a source file,
# that of settings.R
settings_path <- file.path("inst", "scripts", "settings.R")
source_md5 <- vapply(seq_len(nrow(ledger)), function(i) {
  f <- ledger$source_file[i]
  if (!nzchar(f) && ledger$type[i] == "input") f <- settings_path
  if (nzchar(f) && file.exists(f)) unname(tools::md5sum(f)) else ""
}, character(1))
result <- data.frame(
  id = ledger$id, file = ledger$file, type = ledger$type,
  text_value = ledger$text_value,
  formatted = vapply(checked, `[[`, character(1), "formatted"),
  recomputed = vapply(checked, `[[`, character(1), "recomputed"),
  rounding = ledger$rounding,
  result = vapply(checked, `[[`, character(1), "result"),
  message = vapply(checked, `[[`, character(1), "message"),
  source_file = ledger$source_file, source_expr = ledger$source_expr,
  source_md5 = source_md5,
  stringsAsFactors = FALSE
)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(result, out_file, row.names = FALSE, fileEncoding = "UTF-8")
n_fail <- sum(result$result == "FAIL")
message(sprintf(paste0("number_check: %d rows; PASS %d, FAIL %d, ",
                       "LITERATURE %d; rounding rules checked on %d cases. ",
                       "Result in %s"),
                nrow(result), sum(result$result == "PASS"), n_fail,
                sum(result$result == "LITERATURE"), nrow(rule_tests),
                out_file))
if (n_fail > 0L) {
  print(result[result$result == "FAIL",
               c("id", "text_value", "formatted", "message")],
        row.names = FALSE)
}
