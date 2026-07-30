#!/usr/bin/env Rscript
# =============================================================================
# Launch the MedMethods Explorer locally.
#
#   Rscript app/run_local.R              # serves on http://127.0.0.1:7800
#   Rscript app/run_local.R 8080         # a different port
#   Rscript app/run_local.R --check      # check prerequisites and exit
#
# Run it from the repository root (the directory containing app/).
# =============================================================================

args <- commandArgs(trailingOnly = TRUE)
check_only <- "--check" %in% args
port <- suppressWarnings(as.integer(args[!grepl("^--", args)][1]))
if (is.na(port)) port <- 7800L

APP_DIR <- if (dir.exists("app")) "app" else
  if (file.exists("app.R")) "." else
    stop("Could not find the app directory. Run this from the repository root.",
         call. = FALSE)

UI_PKGS <- c("shiny", "bslib", "bsicons", "DT", "plotly", "shinycssloaders",
             "markdown")

cat("MedMethods Explorer\n")
cat(strrep("-", 60), "\n")
cat("R:        ", R.version.string, "\n")
cat("app dir:  ", normalizePath(APP_DIR), "\n\n")

# ---- the engine --------------------------------------------------------------
if (!requireNamespace("MedMethods", quietly = TRUE)) {
  cat("MedMethods is not installed. Installing it now...\n")
  ok <- FALSE
  # a local clone (the package sits at the repository root) comes first
  if (file.exists("DESCRIPTION")) {
    ok <- tryCatch({
      utils::install.packages(".", repos = NULL, type = "source")
      requireNamespace("MedMethods", quietly = TRUE)
    }, error = function(e) FALSE)
  }
  if (!ok) {
    if (!requireNamespace("remotes", quietly = TRUE))
      utils::install.packages("remotes", repos = "https://cloud.r-project.org")
    remotes::install_github("zhaoyi1026/MedMethods", upgrade = "never")
  }
  if (!requireNamespace("MedMethods", quietly = TRUE))
    stop("Could not install MedMethods. Install it manually, then re-run.",
         call. = FALSE)
}
cat("engine:    MedMethods", as.character(utils::packageVersion("MedMethods")),
    sprintf("(%d methods)\n", length(MedMethods::med_methods())))

# ---- the UI packages --------------------------------------------------------
missing <- UI_PKGS[!vapply(UI_PKGS, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing)) {
  cat("installing UI packages:", paste(missing, collapse = ", "), "\n")
  utils::install.packages(missing, repos = "https://cloud.r-project.org")
  still <- missing[!vapply(missing, requireNamespace, logical(1), quietly = TRUE)]
  if (length(still))
    stop("Could not install: ", paste(still, collapse = ", "), call. = FALSE)
}
cat("UI:       ", paste(UI_PKGS, collapse = ", "), "\n")

# `optimx` is only needed for macc()'s default optimizer; the app passes
# optimizer = "bobyqa" instead, so it is genuinely optional.
n_methods <- length(list.files(file.path(APP_DIR, "methods"),
                               pattern = "_method\\.R$", recursive = TRUE))
cat("plugins:  ", n_methods, "method pages\n")

if (check_only) {
  cat("\nAll prerequisites satisfied. Start the app with:\n")
  cat("  Rscript", file.path(APP_DIR, "run_local.R"), "\n")
  quit(status = 0)
}

cat("\nStarting on http://127.0.0.1:", port, "\n", sep = "")
cat("Press Ctrl+C to stop. Your data stays on this machine.\n\n")
shiny::runApp(APP_DIR, host = "127.0.0.1", port = port, launch.browser = TRUE)
