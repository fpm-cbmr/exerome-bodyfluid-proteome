# Shared helpers for repository runners (run_public.R and run_full.R)

suppressPackageStartupMessages(library(here))

assert_files_exist <- function(paths, context = "required inputs") {
  missing <- paths[!file.exists(paths)]
  if (length(missing) > 0) {
    msg <- paste(
      "Missing", context, "(", length(missing), "):",
      paste(sprintf(" - %s", missing), collapse = "\n"),
      sep = "\n"
    )
    stop(msg, call. = FALSE)
  }
}

ensure_dirs <- function(paths) {
  for (p in paths) dir.create(p, recursive = TRUE, showWarnings = FALSE)
}

run_script <- function(path) {
  message("\n=== ", path, " ===")
  status <- system2("Rscript", shQuote(here::here(path)), stdout = "", stderr = "")
  ok <- identical(status, 0L)
  if (!ok) message("  FAILED (exit ", status, ")")
  invisible(ok)
}

run_many <- function(paths) {
  stats <- vapply(paths, run_script, logical(1), USE.NAMES = TRUE)
  message("\n=== summary ===")
  message(sprintf("%d/%d scripts completed", sum(stats), length(stats)))
  if (any(!stats)) {
    failed <- names(stats)[!stats]
    message("failed: ", paste(failed, collapse = ", "))
  }
  invisible(stats)
}
