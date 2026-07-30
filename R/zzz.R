##############################################################################
# Method-module loader.
#
# Each method's original R code lives in inst/method/<id>.R and is sourced into
# its OWN environment at load time, with the package namespace as parent. This
# keeps identically-named internals from different methods apart (for example
# `obj.func` is defined by gmed, hdexposure and multimodal; `soft.thred` by
# hdexposure and multimodal; `BC.CI` by spcma and cfma) while still letting each
# module resolve its imports (nlme::gls, MASS::ginv, ...) through the namespace.
##############################################################################

# id -> environment holding that method's functions
.med_envs <- new.env(parent = emptyenv())

# Load order is irrelevant (R resolves function bodies at call time), but the
# vector also defines the public method catalogue.
.med_module_ids <- c(
  "macc", "gma", "spcma", "cfma", "pathwaylasso",
  "multimodal", "hdexposure", "pcma", "gmed", "heteromed"
)

.onLoad <- function(libname, pkgname) {
  ns <- asNamespace(pkgname)
  for (id in .med_module_ids) {
    path <- system.file("method", paste0(id, ".R"),
                        package = pkgname, lib.loc = libname)
    if (!nzchar(path) || !file.exists(path)) {
      warning("MedMethods: method module '", id, "' not found", call. = FALSE)
      next
    }
    env <- new.env(parent = ns)
    sys.source(path, envir = env, keep.source = FALSE)
    assign(id, env, envir = .med_envs)
  }
  invisible(NULL)
}

##############################################################################
# Internal accessors
##############################################################################

.med_env <- function(method) {
  env <- get0(method, envir = .med_envs, inherits = FALSE)
  if (is.null(env))
    stop("MedMethods: unknown method module '", method, "'. Available: ",
         paste(.med_module_ids, collapse = ", "), call. = FALSE)
  env
}

.med_get <- function(method, fn) {
  env <- .med_env(method)
  f <- get0(fn, envir = env, inherits = FALSE)
  if (!is.function(f))
    stop("MedMethods: '", fn, "' is not a function in method module '",
         method, "'", call. = FALSE)
  f
}

# Forward a wrapper call to the private implementation.
#
# The exported wrappers in R/api.R carry the implementation's formals verbatim
# (they are generated from them by tools/build_medpkg.R), so re-dispatching the
# matched call leaves argument matching, missingness and defaults untouched: any
# argument the caller did not supply is simply absent from `.call` and the
# implementation applies its own default.
.med_call <- function(method, fn, .call, .env) {
  .call[[1L]] <- .med_get(method, fn)
  eval(.call, .env)
}

##############################################################################
# Public escape hatches
##############################################################################

#' Method modules bundled in MedMethods
#'
#' @return A character vector of method-module ids.
#' @export
med_methods <- function() .med_module_ids

#' Reach an internal function of a method module
#'
#' Every method ships as its original R code, sourced into a private
#' environment. `med_internal()` returns one of those functions so that
#' intermediate helpers, which have no exported wrapper, remain usable.
#'
#' @param method Method-module id; see [med_methods()].
#' @param fn Name of a function inside that module. If missing, the names of all
#'   functions in the module are returned.
#' @return The requested function, or a character vector of available names when
#'   `fn` is missing.
#' @export
med_internal <- function(method, fn) {
  if (missing(fn)) return(sort(ls(.med_env(method))))
  .med_get(method, fn)
}
