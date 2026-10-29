# Copy the tables and figures of the article from the results of the full run
# to the folders of the manuscript, without subfolders, and record the MD5 of
# every source and copy:
#   inst/scripts/results/tables/*.tex, figures/*.pdf, *.eps
#     -> Manuscript/Main_body/
#   inst/sup_info/results/tables/*.tex, figures/*.pdf, *.eps
#     -> Manuscript/Supp_info/
# The results of the quick test run (results/quick/) are not copied, and only
# .tex tables and .pdf and .eps figures are copied. Copies of earlier runs are
# not removed: verify_copies.py reports a table or figure of the manuscript
# folders that has no source any more, which is then moved by hand to
# _to_delete/.
# The record is written to Manuscript/tools/logs/copy_results.csv and checked
# by Manuscript/tools/verify_copies.py.
#
# Run from the root of the package after the table and figure scripts:
#   source("inst/scripts/copy_results.R")
# The manuscript folder is ../../Manuscript (the package is in
# <project>/Rpackage/onestepcox); another folder can be given with
# options(onestepcox.manuscript = "<path>").

if (!file.exists("DESCRIPTION") ||
    !identical(unname(read.dcf("DESCRIPTION", fields = "Package")[1L, 1L]),
               "onestepcox")) {
  stop("Run the script from the root of the onestepcox package")
}
manuscript <- getOption("onestepcox.manuscript",
                        file.path("..", "..", "Manuscript"))
targets <- list(
  list(from = file.path("inst", "scripts", "results"),
       to = file.path(manuscript, "Main_body")),
  list(from = file.path("inst", "sup_info", "results"),
       to = file.path(manuscript, "Supp_info"))
)
for (tg in targets) {
  if (!dir.exists(tg$to)) stop("folder not found: ", tg$to)
}

copied <- lapply(targets, function(tg) {
  files <- c(list.files(file.path(tg$from, "tables"), pattern = "\\.tex$",
                        full.names = TRUE),
             list.files(file.path(tg$from, "figures"),
                        pattern = "\\.(pdf|eps)$", full.names = TRUE))
  if (length(files) == 0L) {
    stop("no tables or figures in ", tg$from,
         "; run the table and figure scripts first")
  }
  dest <- file.path(tg$to, basename(files))
  ok <- file.copy(files, dest, overwrite = TRUE)
  if (!all(ok)) {
    stop("the copy failed: ", paste(files[!ok], collapse = ", "))
  }
  md5_source <- unname(tools::md5sum(files))
  md5_copy <- unname(tools::md5sum(dest))
  data.frame(file = basename(files), source = files, destination = dest,
             md5_source = md5_source, md5_copy = md5_copy,
             result = ifelse(md5_source == md5_copy, "PASS", "FAIL"),
             stringsAsFactors = FALSE)
})
copied <- do.call(rbind, copied)
copied$date <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
log_dir <- file.path(manuscript, "tools", "logs")
dir.create(log_dir, recursive = TRUE, showWarnings = FALSE)
utils::write.csv(copied, file.path(log_dir, "copy_results.csv"),
                 row.names = FALSE)
message(sprintf(paste0("copy_results: %d files copied; MD5 identical: %d, ",
                       "different: %d. Record in %s"),
                nrow(copied), sum(copied$result == "PASS"),
                sum(copied$result != "PASS"),
                file.path(log_dir, "copy_results.csv")))
if (!isTRUE(all(copied$result == "PASS"))) {
  stop("the MD5 of a copy differs from its source: ",
       paste(copied$file[copied$result != "PASS"], collapse = ", "))
}
