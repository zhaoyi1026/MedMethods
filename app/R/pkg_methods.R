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
