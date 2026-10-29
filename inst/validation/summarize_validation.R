# Summary of the reproduction and validation results.
#
# Reads every results table written by the scripts in inst/reproduce/ and
# inst/validation/ (each with a column 'status' taking the values PASS,
# EXPLAINED, FAIL, and INFO), and writes the number of checks with each
# status per table to inst/validation/results/summary.csv. The rows with
# FAIL or EXPLAINED are printed. The summary is generated from the results;
# it is not written by hand.
#
# Run from the root of the package after the other scripts.

files <- c(list.files("inst/reproduce/results", pattern = "[.]csv$",
                      full.names = TRUE),
           list.files("inst/validation/results", pattern = "[.]csv$",
                      full.names = TRUE))
files <- files[basename(files) != "summary.csv"]
if (length(files) == 0L) stop("no results found; run the scripts first")

levels_status <- c("PASS", "EXPLAINED", "FAIL", "INFO")
summary_rows <- lapply(files, function(f) {
  r <- utils::read.csv(f, stringsAsFactors = FALSE)
  if (!"status" %in% names(r)) stop(f, " has no column 'status'")
  if (!all(r$status %in% levels_status)) stop(f, " has an unknown status")
  counts <- table(factor(r$status, levels = levels_status))
  bad <- r[r$status %in% c("FAIL", "EXPLAINED"), , drop = FALSE]
  if (nrow(bad) > 0L) {
    cat("\n", f, ": rows with FAIL or EXPLAINED\n", sep = "")
    print(bad)
  }
  data.frame(file = f, checks = nrow(r), PASS = counts[["PASS"]],
             EXPLAINED = counts[["EXPLAINED"]], FAIL = counts[["FAIL"]],
             INFO = counts[["INFO"]], stringsAsFactors = FALSE)
})
summary_tab <- do.call(rbind, summary_rows)
dir.create("inst/validation/results", showWarnings = FALSE)
utils::write.csv(summary_tab, "inst/validation/results/summary.csv",
                 row.names = FALSE)
print(summary_tab, row.names = FALSE)
if (any(summary_tab$FAIL > 0)) {
  cat("\nThere are FAIL results.\n")
} else {
  cat("\nNo FAIL results.\n")
}
