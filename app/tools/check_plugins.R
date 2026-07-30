##############################################################################
# Headless check of every app plugin.
#
#   Rscript app/tools/check_plugins.R          (from the repository root)
#
# For each registered method this runs the whole pipeline the Shiny page would:
#   example() -> run(params from the declared defaults) -> summarize() -> plots()
# and also round-trips the example through export_example() -> parse() so the
# upload path is exercised with data of exactly the shape the app advertises.
#
# It renders no UI, so it catches method/plumbing errors without a browser.
##############################################################################
suppressPackageStartupMessages({
  library(shiny); library(bslib); library(bsicons)
  library(DT); library(plotly); library(markdown)
})

args <- commandArgs(trailingOnly = TRUE)
APP_DIR <- if (length(args) >= 1) normalizePath(args[1]) else {
  cand <- c("app", file.path("MedMethods-Rpkg", "260729", "MedMethods", "app"), ".")
  hit <- cand[file.exists(file.path(cand, "app.R"))]
  if (!length(hit)) stop("Could not locate the app directory; pass it as an argument.")
  normalizePath(hit[1])
}
cat("app dir:", APP_DIR, "\n\n")

source(file.path(APP_DIR, "R", "registry.R"), local = TRUE)
for (f in c("io_helpers.R", "ui_helpers.R", "mod_method.R", "pkg_methods.R"))
  source(file.path(APP_DIR, "R", f), local = TRUE)
for (mf in list.files(file.path(APP_DIR, "methods"), pattern = "_method\\.R$",
                      recursive = TRUE, full.names = TRUE))
  source(mf, local = TRUE)

METHODS <- list_methods()
cat(sprintf("%d plugin(s) registered: %s\n\n", length(METHODS),
            paste(names(METHODS), collapse = ", ")))

quiet <- function(expr) {
  lf <- tempfile(); con <- file(lf, open = "wt")
  sink(con); sink(con, type = "message")
  v <- tryCatch(suppressWarnings(expr),
                error = function(e) structure(conditionMessage(e), class = "cfail"))
  sink(type = "message"); sink(); close(con); unlink(lf)
  v
}
failed <- function(x) inherits(x, "cfail")

# defaults exactly as the sidebar would supply them
default_params <- function(spec) {
  out <- list()
  for (p in spec$params) out[[p$id]] <- p$default
  out
}

results <- list()
for (id in names(METHODS)) {
  spec <- METHODS[[id]]
  cat(sprintf("=== %-12s %s\n", id, spec$full_name))
  rec <- list(id = id, ok = TRUE, notes = character(0))
  t0 <- proc.time()[["elapsed"]]

  d <- quiet(spec$example())
  if (failed(d)) {
    cat("    example()   FAIL:", as.character(d), "\n"); rec$ok <- FALSE
    results[[id]] <- rec; next
  }
  cat("    example()   OK  ", spec$describe_data(d), "\n")

  # round-trip through the documented upload format
  if (is.function(spec$export_example) && is.function(spec$parse)) {
    ex <- quiet(spec$export_example(d))
    if (failed(ex)) {
      cat("    export()    FAIL:", as.character(ex), "\n"); rec$ok <- FALSE
    } else {
      d2 <- quiet(spec$parse(ex, list(add_intercept = TRUE)))
      if (failed(d2)) {
        cat("    parse()     FAIL:", as.character(d2), "\n"); rec$ok <- FALSE
      } else {
        cat("    parse()     OK   (round-tripped the exported example)\n")
      }
    }
  }

  res <- quiet(spec$run(d, default_params(spec)))
  if (failed(res)) {
    cat("    run()       FAIL:", as.character(res), "\n"); rec$ok <- FALSE
    results[[id]] <- rec; next
  }
  cat(sprintf("    run()       OK   (%.1fs)\n", proc.time()[["elapsed"]] - t0))

  tabs <- if (is.function(spec$summarize)) quiet(spec$summarize(res)) else list()
  if (failed(tabs)) {
    cat("    summarize() FAIL:", as.character(tabs), "\n"); rec$ok <- FALSE
  } else {
    bad <- names(tabs)[!vapply(tabs, is.data.frame, logical(1))]
    cat(sprintf("    summarize() OK   %d table(s): %s\n", length(tabs),
                paste(names(tabs), collapse = " | ")))
    if (length(bad)) {
      cat("      NOT data.frames:", paste(bad, collapse = ", "), "\n"); rec$ok <- FALSE
    }
    tv <- tabs[["Truth vs estimate"]]
    if (!is.null(tv)) {
      cat("      truth vs estimate:\n")
      print(utils::head(tv, 8), row.names = FALSE, digits = 4)
    }
  }

  pls <- if (is.function(spec$plots)) quiet(spec$plots(res)) else list()
  if (failed(pls)) {
    cat("    plots()     FAIL:", as.character(pls), "\n"); rec$ok <- FALSE
  } else {
    bad <- names(pls)[!vapply(pls, function(p) !is.null(p$plot), logical(1))]
    cat(sprintf("    plots()     OK   %d plot(s): %s\n", length(pls),
                paste(names(pls), collapse = " | ")))
    if (length(bad)) {
      cat("      missing $plot:", paste(bad, collapse = ", "), "\n"); rec$ok <- FALSE
    }
  }
  cat("\n")
  results[[id]] <- rec
}

nfail <- sum(!vapply(results, function(r) r$ok, logical(1)))
cat(strrep("=", 70), "\n")
cat(sprintf("%d of %d plugin(s) fully OK\n", length(results) - nfail, length(results)))
if (nfail) {
  cat("failing:", paste(names(Filter(function(r) !r$ok, results)), collapse = ", "), "\n")
  quit(status = 1)
}
