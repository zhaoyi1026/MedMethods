# =============================================================================
# Bridge to the installed MedMethods package.
# -----------------------------------------------------------------------------
# The app uses the installed package as its engine: nothing here sources the
# original method folders. Unlike the CAP app, no private-environment cloning is
# needed, because MedMethods exports every method wrapper publicly -- a plugin
# just calls MedMethods::gmed(), MedMethods::pathlasso(), and so on.
#
# med_pkg() checks the package is available and returns its namespace, so a
# missing install produces one clear message instead of ten "object not found"
# errors. med_internal_fn() reaches an unexported helper of a method module for
# the rare case a plugin needs one.
# =============================================================================

med_pkg <- function() {
  if (!requireNamespace("MedMethods", quietly = TRUE))
    stop("The 'MedMethods' package is required but not installed. Install it with ",
         "remotes::install_github(\"zhaoyi1026/MedMethods\"), or from a local ",
         "clone with install.packages(\".\", repos = NULL, type = \"source\").",
         call. = FALSE)
  asNamespace("MedMethods")
}

#' Fetch an exported MedMethods function by name.
med_fn <- function(name) {
  ns <- med_pkg()
  f <- get0(name, envir = ns, inherits = FALSE)
  if (!is.function(f))
    stop(sprintf("MedMethods does not export '%s'.", name), call. = FALSE)
  f
}

#' Fetch an internal function of one of the package's method modules.
#' Thin wrapper over MedMethods::med_internal(method, fn).
med_internal_fn <- function(method, fn) {
  med_fn("med_internal")(method, fn)
}

#' Which method modules the installed package carries.
med_available <- function() {
  tryCatch(med_fn("med_methods")(), error = function(e) character(0))
}

#' Installed package version, for the footer.
med_version <- function() {
  tryCatch(as.character(utils::packageVersion("MedMethods")),
           error = function(e) NA_character_)
}

# The app and the package live in one repository but are installed separately, so
# they can drift: pulling new app code while an older MedMethods is still
# installed produced "unused arguments" errors with no hint of the cause. This
# check names the problem instead. Plugins that can work either way (see the gma
# page) degrade gracefully rather than relying on it.
MED_MIN_VERSION <- "0.2.0"

#' Where the engine is installed, so a stale copy can be located.
med_libpath <- function() {
  p <- tryCatch(dirname(system.file(package = "MedMethods")), error = function(e) NA)
  if (length(p) != 1 || is.na(p) || !nzchar(p)) NA_character_ else p
}

#' Path to the package source, if the app is running from inside the repository.
#'
#' NB under shiny::runApp("app") the working directory is the APP directory, not
#' the repository root -- so looking for DESCRIPTION in getwd() always fails and
#' wrongly sends people to install_github. Look one level up as well.
med_source_dir <- function() {
  cands <- c(".", "..", if (exists("APP_DIR", inherits = TRUE))
    file.path(get("APP_DIR", inherits = TRUE), ".."))
  for (d in cands) {
    if (file.exists(file.path(d, "DESCRIPTION")) &&
        any(grepl("^Package:\\s*MedMethods",
                  readLines(file.path(d, "DESCRIPTION"), warn = FALSE))))
      return(normalizePath(d))
  }
  NA_character_
}

#' The one command that fixes a stale install, given where the app is running.
med_fix_command <- function() {
  src <- med_source_dir()
  if (!is.na(src))
    sprintf("install.packages('%s', repos = NULL, type = 'source')", src)
  else
    "remotes::install_github('zhaoyi1026/MedMethods')"
}

#' Report any mismatch between the app and the installed engine.
#'
#' Checks the version string AND the actual capabilities, because a package can
#' be reinstalled without the version changing. Every problem is phrased so the
#' reader knows what to do about it.
#' @return NULL if all good, otherwise a character vector of problems.
med_engine_problems <- function() {
  v <- med_version()
  if (is.na(v)) return("MedMethods is not installed.")
  lp <- med_libpath()

  # --- Case 1: reinstalled, but this R session still holds the OLD namespace ---
  # packageVersion() reads the DESCRIPTION on disk; the functions the app calls
  # come from the loaded namespace. After reinstalling without restarting R, the
  # two disagree -- and no amount of reinstalling fixes it. Detect that first, or
  # the advice sends people round in circles.
  loaded <- tryCatch(as.character(getNamespaceVersion("MedMethods")),
                     error = function(e) NA_character_)
  if (!is.na(loaded) && !identical(loaded, v))
    return(c(
      sprintf("This R session is still using MedMethods %s, but %s is installed on disk.",
              loaded, v),
      "The package was reinstalled without restarting R, so the session kept the old copy. Reinstalling again will NOT help.",
      "Restart R (in RStudio: Session > Restart R) and start the app again.",
      sprintf("Installed at: %s", if (is.na(lp)) "(unknown)" else lp)))

  probs <- character(0)
  if (utils::compareVersion(v, MED_MIN_VERSION) < 0)
    probs <- c(probs, sprintf(
      "Version %s is installed; this app expects %s or newer.", v, MED_MIN_VERSION))

  # --- Case 2: genuinely old install (capability check, version-independent) ---
  gen <- tryCatch(med_fn("gma_example"), error = function(e) NULL)
  if (!is.null(gen) && !("model.type" %in% names(formals(gen))))
    probs <- c(probs, paste(
      "The installed gma_example() has no 'model.type' argument, so it predates",
      "the two-level example. The gma page still works -- it falls back to its own",
      "generator -- but the installed package is out of date."))

  if (!length(probs)) return(NULL)
  c(probs,
    sprintf("Installed at: %s  (R %s.%s)", if (is.na(lp)) "(unknown)" else lp,
            R.version$major, R.version$minor),
    "Note that packages install per R version -- installing under a different R will not update this one.",
    sprintf("To fix, run this in R, then RESTART R and start the app again:  %s",
            med_fix_command()))
}
