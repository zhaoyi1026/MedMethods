##############################################################################
# Assemble the MedMethods R package from the per-method source folders.
#
# Reproducible: re-running regenerates inst/method/*.R from the original
# method code in the sibling folders listed in ../setting.md.
#
# Usage:  Rscript tools/build_medpkg.R [project_root] [pkg_dir]
##############################################################################

args <- commandArgs(trailingOnly = TRUE)
ROOT <- if (length(args) >= 1) args[1] else
  "/Users/yizhao/Dropbox/MyFolder/Biostat-IU/Projects/Yi/mediation/code"
PKG  <- if (length(args) >= 2) args[2] else
  file.path(ROOT, "MedMethods-Rpkg", "260729", "MedMethods")

cat("project root :", ROOT, "\n")
cat("package dir  :", PKG, "\n\n")

##############################################################################
# Method manifest: id -> source files (relative to ROOT), in source order.
# Each method's code is sourced into its OWN private environment at .onLoad,
# so identically-named internals across methods (obj.func, soft.thred, BC.CI,
# eigen.solve, ...) never collide and no renaming is needed.
##############################################################################
METHODS <- list(
  macc = list(
    title = "Multilevel mediation analysis with structured unmeasured confounding",
    dir   = "macc/macc/R"
  ),
  gma = list(
    title = "Granger mediation analysis",
    dir   = "gma/gma/R"
  ),
  spcma = list(
    title = "Sparse principal component based mediation analysis",
    dir   = "spcma/spcma-master/R"
  ),
  cfma = list(
    title = "Causal functional mediation analysis",
    dir   = "cfma/cfma-master/R"
  ),
  pathwaylasso = list(
    title = "Pathway Lasso",
    files = c("pathway_lasso/V1/functions.R",
              "pathway_lasso/V1/ADMM_adp_functions.R",
              "pathway_lasso/V1/VSS_adp_functions.R")
  ),
  multimodal = list(
    title = "Multimodal mediation analysis with two blocks of mediators",
    files = "multimodal/multimodal_integration-master/PathLasso.R"
  ),
  hdexposure = list(
    title = "Mediation with high-dimensional exposures and mediators",
    files = "HDExposureMediator/HDExposureMediator-main/HD-CauseMediation.R"
  ),
  pcma = list(
    title = "Principal component mediation analysis",
    files = "PCMA/PCMA-main/functions.R"
  ),
  gmed = list(
    title = "Mediation analysis with a graph mediator",
    files = "GMed/GMed-main/CAPMediation.R"
  ),
  heteromed = list(
    title = "Heterogeneous mediation effects",
    files = "hetero_mediation/V1/HeterMed.R"
  )
)

##############################################################################
# Strip top-level library()/require()/source()/load() calls. The package
# resolves these through NAMESPACE imports instead.
##############################################################################
strip_loaders <- function(lines) {
  pat <- "^[[:space:]]*(library|require|suppressMessages\\(library|source|load)[[:space:]]*\\("
  drop <- grepl(pat, lines)
  if (any(drop)) {
    for (i in which(drop))
      lines[i] <- paste0("# [MedMethods] removed at assembly: ", trimws(lines[i]))
  }
  lines
}

collect <- function(spec) {
  files <- if (!is.null(spec$files)) file.path(ROOT, spec$files) else
    sort(list.files(file.path(ROOT, spec$dir), pattern = "\\.R$", full.names = TRUE))
  missing <- files[!file.exists(files)]
  if (length(missing)) stop("missing source file(s): ", paste(missing, collapse = ", "))
  out <- character(0)
  for (f in files) {
    out <- c(out,
             sprintf("### ---- from %s ----", sub(paste0("^", ROOT, "/"), "", f)),
             strip_loaders(readLines(f, warn = FALSE)),
             "")
  }
  list(lines = out, files = files)
}

##############################################################################
# gma correction: the asymptotic-variance branch built a companion matrix with
# mismatched blocks, so cma.uni.delta.ts.arp.error(var.asmp=TRUE) errored for
# any VAR lag p >= 2 ("number of columns of matrices must match"). t(W.hat) is
# 2 x 2p while the lower block was 2(p-1) x (2p-1); it must be 2(p-1) x 2p.
# Also thread var.asmp through the sensitivity wrapper so gma()'s
# single.var.asmp choice reaches it (mirrors gma/gma_functions.R).
##############################################################################
patch_gma <- function(lines) {
  n_before <- length(lines)
  hit <- rep(FALSE, 3)

  # (1) companion matrix
  i <- grep("Fm<-rbind(cbind(t(W.hat)),cbind(diag(rep(1,2*(p-1))),matrix(0,2*(p-1),1)))",
            lines, fixed = TRUE)
  if (length(i) == 1L) {
    lines[i] <- sub("matrix(0,2*(p-1),1)", "matrix(0,2*(p-1),2)", lines[i], fixed = TRUE)
    lines[i] <- sub("rbind(cbind(t(W.hat)),", "rbind(t(W.hat),", lines[i], fixed = TRUE)
    hit[1] <- TRUE
  }

  # (2) var.asmp argument on the sensitivity wrapper
  sig <- "function(dat,delta=seq(-1,1,by=0.01),p=1,conf.level=0.95)"
  i <- which(trimws(lines) == sig)
  if (length(i) == 1L) {
    lines[i] <- sub(",conf.level=0.95)", ",conf.level=0.95,var.asmp=TRUE)", lines[i], fixed = TRUE)
    hit[2] <- TRUE
  }
  i <- grep("re<-cma.uni.delta.ts.arp.error(dat,delta=delta[i],p=p,conf.level=conf.level)",
            lines, fixed = TRUE)
  if (length(i) == 1L)
    lines[i] <- sub("conf.level=conf.level)", "conf.level=conf.level,var.asmp=var.asmp)",
                    lines[i], fixed = TRUE)

  # (3) gma() passes its single.var.asmp choice into the sensitivity analysis
  i <- grep("re.cma.sens<-cma.uni.sens.ts.arp.error(dat,delta=sens.delta,p=p,conf.level=conf.level)",
            lines, fixed = TRUE)
  if (length(i) == 1L) {
    lines[i] <- sub("conf.level=conf.level)",
                    "conf.level=conf.level,var.asmp=single.var.asmp)", lines[i], fixed = TRUE)
    hit[3] <- TRUE
  }

  if (!all(hit))
    stop("gma patch failed to apply cleanly (companion/sens-arg/gma-call = ",
         paste(hit, collapse = "/"), ")")
  stopifnot(length(lines) == n_before)
  lines
}

##############################################################################
# gmed correction: CAPMediation_coef() and CAPMediation_refit() build
#   M.cov <- array(NA, c(p, p, n))
# but never define n or p -- both were leaking in from the calling script's
# workspace. Inside a package the lookup escapes to the user's globalenv, so the
# functions fail ("object 'p' not found") or, worse, silently pick up unrelated
# values. Derive both from M, as every other function in the module does.
# Only these two sites are affected; the other six M.cov constructions set n and
# p from M beforehand.
##############################################################################
patch_gmed <- function(lines) {
  targets <- c("CAPMediation_coef<-function(X,M,Y,theta)",
               "CAPMediation_refit<-function(X,M,Y,Theta)")
  alloc <- "  M.cov<-array(NA,c(p,p,n))"
  for (tg in targets) {
    start <- grep(tg, lines, fixed = TRUE)
    if (length(start) != 1L) stop("gmed patch: cannot locate ", tg)
    j <- which(lines == alloc & seq_along(lines) > start)
    if (!length(j)) stop("gmed patch: no M.cov allocation after ", tg)
    j <- j[1L]
    lines <- append(lines,
                    c("  n<-length(M)", "  p<-ncol(M[[1]])"), after = j - 1L)
  }
  lines
}

##############################################################################
# pathwaylasso correction: sim.data_dep() takes the exposure -> mediator
# coefficients as `a`, but the j == 1 branch reads `A[1,1]`. In the original
# example script a global `A` held the same values, so the slip was invisible;
# in a package the lookup escapes to the user's workspace. M.coef <- rbind(a,
# Delta) makes M.coef[1,1] == a[1,1], so the first mediator's coefficient is
# a[1,1] -- consistent with the j > 1 branch.
##############################################################################
patch_pathwaylasso <- function(lines) {
  bad <- "      M[,j]<-Z*A[1,1]+eps1[,j]"
  i <- which(lines == bad)
  if (length(i) != 1L) stop("pathwaylasso patch: cannot locate the A[1,1] line")
  lines[i] <- "      M[,j]<-Z*a[1,1]+eps1[,j]"
  lines
}

##############################################################################
# heteromed correction: fit.inf.OLS() applies sqrt() TWICE to the coefficient
# standard errors. It computes
#     vecTheta.se <- sqrt(diag(cov.vecTheta))   # already standard errors
#     Theta.se    <- matrix(vecTheta.se, ncol = 2)
# and then builds every coefficient table with `SE = sqrt(Theta.se[...])`, i.e.
# the square root of a standard error -- dimensionally incoherent, and it
# inflates every SE, z-value, p-value and CI in the alpha/beta/gamma tables.
#
# Confirmed empirically two ways on hetermed_example():
#   (a) the reported SE scales as n^(-1/4), not n^(-1/2): the ratio
#       SE(n=600)/SE(n=2400) is 1.42, whereas a correct SE gives 2.00;
#   (b) the SQUARE of the reported SE matches the Monte Carlo SD of alpha1-hat
#       over 300 replicates, while the reported SE itself is ~7x too large.
# The NIE/NDE tables use cov.vecTheta directly (t(h) %*% cov.vecTheta %*% h) and
# are unaffected.
#
# Fix: drop the six spurious sqrt() calls so the tables report Theta.se itself.
##############################################################################
patch_heteromed <- function(lines) {
  idx <- grep("^\\s*(alpha0|alpha1|gamma0|gamma1|beta0|beta1)\\.out<-data\\.frame\\(Estimate=.*SE=sqrt\\(Theta\\.se\\[",
              lines)
  if (length(idx) != 6L)
    stop("heteromed patch: expected 6 coefficient-table lines, found ", length(idx))
  for (i in idx) {
    before <- lines[i]
    lines[i] <- sub("SE=sqrt(Theta.se[", "SE=(Theta.se[", lines[i], fixed = TRUE)
    if (identical(before, lines[i]))
      stop("heteromed patch: substitution failed on line ", i)
  }
  lines
}

##############################################################################
# Write inst/method/<id>.R
##############################################################################
dir.create(file.path(PKG, "inst", "method"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(PKG, "R"),              recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(PKG, "tools"),          recursive = TRUE, showWarnings = FALSE)

for (id in names(METHODS)) {
  spec <- METHODS[[id]]
  got  <- collect(spec)
  body <- got$lines
  body <- switch(id,
                 gma          = patch_gma(body),
                 gmed         = patch_gmed(body),
                 pathwaylasso = patch_pathwaylasso(body),
                 heteromed    = patch_heteromed(body),
                 body)

  header <- c(
    "##############################################################################",
    sprintf("# MedMethods method module: %s", id),
    sprintf("# %s", spec$title),
    "#",
    "# Assembled by tools/build_medpkg.R from the original method sources.",
    "# Sourced into a private environment at .onLoad (see R/zzz.R), so internal",
    "# helper names may safely collide with those of other method modules.",
    "# Do not edit by hand -- edit the source files and re-run the build script.",
    "##############################################################################",
    ""
  )
  out     <- file.path(PKG, "inst", "method", paste0(id, ".R"))
  content <- c(header, body)
  writeLines(content, out)

  # parse check + count of top-level function definitions
  exprs <- tryCatch(parse(out), error = function(e) {
    cat("  PARSE ERROR:", conditionMessage(e), "\n"); NULL
  })
  nfun <- if (is.null(exprs)) NA_integer_ else
    sum(vapply(exprs, function(e)
      is.call(e) && as.character(e[[1]]) %in% c("<-", "=") &&
        is.call(e[[3]]) && identical(as.character(e[[3]][[1]]), "function"),
      logical(1)))
  cat(sprintf("%-14s %2d file(s) -> %-16s %5d lines, %2s functions, parse %s\n",
              id, length(got$files), basename(out), length(content), nfun,
              if (is.null(exprs)) "FAIL" else "OK"))
}

cat("\nDone.\n")

##############################################################################
# Generate R/api.R: exported wrappers whose formals are copied verbatim from
# the private implementations, so signatures can never drift from the method
# code. Each wrapper re-dispatches its matched call via .med_call().
##############################################################################

API <- list(
  # macc -----------------------------------------------------------------
  list("macc",              "macc",         "macc",
       "Multilevel mediation analysis under structured unmeasured confounding"),
  list("macc_sim_single",   "macc",         "sim.data.single",
       "Simulate single-level data for macc()"),
  list("macc_sim_multi",    "macc",         "sim.data.multi",
       "Simulate multilevel data for macc()"),
  # gma ------------------------------------------------------------------
  list("gma",               "gma",          "gma",
       "Granger mediation analysis of time series"),
  list("gma_sim_single",    "gma",          "sim.data.ts.single",
       "Simulate a single time series for gma()"),
  list("gma_sim_two",       "gma",          "sim.data.ts.two",
       "Simulate two-level time series for gma()"),
  # spcma ----------------------------------------------------------------
  list("spcma",             "spcma",        "spcma",
       "Sparse principal component based high-dimensional mediation analysis"),
  list("mcma_bk",           "spcma",        "mcma_BK",
       "Baron-Kenny style multivariate mediation analysis"),
  list("mcma_pca",          "spcma",        "mcma_PCA",
       "Principal component based mediation analysis"),
  list("plot_spcma",        "spcma",        "plot_spcma",
       "Plot spcma coefficient estimates"),
  # cfma -----------------------------------------------------------------
  list("cfma_concurrent",      "cfma", "FMA.concurrent",
       "Concurrent functional mediation model"),
  list("cfma_concurrent_cv",   "cfma", "FMA.concurrent.CV",
       "Cross-validation for the concurrent functional mediation model"),
  list("cfma_concurrent_boot", "cfma", "FMA.concurrent.boot",
       "Bootstrap inference for the concurrent functional mediation model"),
  list("cfma_historical",      "cfma", "FMA.historical",
       "Historical-influence functional mediation model"),
  list("cfma_historical_cv",   "cfma", "FMA.historical.CV",
       "Cross-validation for the historical-influence functional mediation model"),
  list("cfma_historical_boot", "cfma", "FMA.historical.boot",
       "Bootstrap inference for the historical-influence functional mediation model"),
  # pathway lasso --------------------------------------------------------
  list("pathlasso",          "pathwaylasso", "pathlasso",
       "Pathway Lasso estimation and selection with high-dimensional mediators"),
  list("pathlasso_ksc",      "pathwaylasso", "mediation_net_ADMM_NC_KSC",
       "Pathway Lasso tuning by Kappa selection criterion"),
  list("pathlasso_vss",      "pathwaylasso", "mediation_net_ADMM_NC_VSS",
       "Variable selection stability for a single Pathway Lasso split"),
  list("pathlasso_vss_rep",  "pathwaylasso", "mediation_net_ADMM_NC_VSS_rep",
       "Variable selection stability over a Pathway Lasso tuning path"),
  list("pathlasso_sim",      "pathwaylasso", "sim.data_dep",
       "Simulate correlated-error data for pathlasso()"),
  # multimodal -----------------------------------------------------------
  list("pathlasso2b",       "multimodal",   "pathlasso.2b",
       "Multimodal pathway analysis with two blocks of high-dimensional mediators"),
  # HD exposures + mediators ---------------------------------------------
  list("hdmediation",       "hdexposure",   "HDCauseMediation",
       "Mediation analysis with high-dimensional exposures and mediators"),
  list("hdmediation_pca",   "hdexposure",   "HDCauseMediationPCA",
       "Principal-component variant of hdmediation()"),
  # PCMA -----------------------------------------------------------------
  list("pcma",              "pcma",         "HDEM.loglike.opt",
       "Principal component mediation analysis for multiple exposures and mediators"),
  list("pcma_inf",          "pcma",         "HDEM.inf",
       "Bootstrap inference for pcma()"),
  list("pcma_inf_asmp",     "pcma",         "HDEM.inf.asmp",
       "Asymptotic inference for pcma()"),
  list("pcma_coef",         "pcma",         "HDEM.coef",
       "Effect estimates at fixed pcma() projections"),
  # GMed -----------------------------------------------------------------
  list("gmed",              "gmed",         "CAPMediation",
       "Mediation analysis with a graph (covariance matrix) mediator"),
  list("gmed_boot",         "gmed",         "CAPMediation_boot",
       "Bootstrap inference for gmed()"),
  list("gmed_refit",        "gmed",         "CAPMediation_refit",
       "Refit gmed() effects at fixed projections"),
  list("gmed_refit_boot",   "gmed",         "CAPMediation_refit_boot",
       "Bootstrap inference for gmed_refit()"),
  list("gmed_coef",         "gmed",         "CAPMediation_coef",
       "Effect estimates at a fixed gmed() projection"),
  # heterogeneous mediation ----------------------------------------------
  list("hetermed",          "heteromed",    "med.inter",
       "Heterogeneous (moderated) causal mediation effect estimation"),
  list("hetermed_inf",      "heteromed",    "med.inter.inf",
       "Inference for hetermed()"),
  list("hetermed_ite",      "heteromed",    "med.inter.ITE",
       "Individual treatment/mediation effects from hetermed() parameters")
)

# --- load the assembled modules so we can read real formals ----------------
mod_env <- list()
for (id in names(METHODS)) {
  e <- new.env(parent = globalenv())
  suppressWarnings(sys.source(file.path(PKG, "inst", "method", paste0(id, ".R")),
                              envir = e, keep.source = FALSE))
  mod_env[[id]] <- e
}

dep1 <- function(x) paste(deparse(x, width.cutoff = 500L), collapse = " ")

# Render "a, b = 1, ... " wrapped to a sane width with hanging indent.
wrap_args <- function(parts, indent) {
  out <- character(0); cur <- ""
  for (k in seq_along(parts)) {
    piece <- paste0(parts[k], if (k < length(parts)) "," else "")
    cand <- if (nzchar(cur)) paste(cur, piece) else piece
    if (nchar(cand) + nchar(indent) > 92L && nzchar(cur)) {
      out <- c(out, cur); cur <- piece
    } else cur <- cand
  }
  c(out, cur)
}

api_lines <- c(
  "##############################################################################",
  "# Exported wrappers.",
  "#",
  "# GENERATED by tools/build_medpkg.R -- do not edit by hand.",
  "#",
  "# Every signature below is copied verbatim from the corresponding function in",
  "# inst/method/<module>.R, so the public API cannot drift from the method code.",
  "# Each body re-dispatches the matched call into the module's private",
  "# environment, leaving argument matching and defaults to the implementation.",
  "##############################################################################",
  ""
)

for (spec in API) {
  wname <- spec[[1]]; mod <- spec[[2]]; fname <- spec[[3]]; title <- spec[[4]]
  f <- get0(fname, envir = mod_env[[mod]], inherits = FALSE)
  if (!is.function(f))
    stop("API entry '", wname, "': ", mod, "::", fname, " not found")
  fmls <- formals(f)
  # NB: index fmls directly inside deparse(); binding an empty default to a
  # local variable would turn it into a missing-argument promise. deparse()
  # yields "" both for arguments without a default and for `...`.
  parts <- vapply(seq_along(fmls), function(i) {
    nm <- names(fmls)[i]
    d  <- dep1(fmls[[i]])
    if (!nzchar(d)) nm else paste0(nm, " = ", d)
  }, character(1), USE.NAMES = FALSE)

  head_txt <- paste0(wname, " <- function(")
  if (!length(parts)) {
    sig <- paste0(head_txt, ") {")
  } else {
    wrapped <- wrap_args(parts, strrep(" ", nchar(head_txt)))
    sig <- paste0(head_txt, wrapped[1])
    if (length(wrapped) > 1)
      sig <- c(sig, paste0(strrep(" ", nchar(head_txt)), wrapped[-1]))
    sig[length(sig)] <- paste0(sig[length(sig)], ") {")
  }

  api_lines <- c(api_lines,
    sprintf("#' %s", title),
    "#'",
    sprintf("#' Wrapper for `%s()` in the `%s` method module.", fname, mod),
    "#'",
    sprintf("#' @param ... Passed to the implementation; see `med_internal(\"%s\", \"%s\")`.", mod, fname),
    sprintf("#' @return The value of `%s::%s()`.", mod, fname),
    "#' @export",
    sig,
    sprintf('  .med_call("%s", "%s", match.call(), parent.frame())', mod, fname),
    "}",
    ""
  )
}

writeLines(api_lines, file.path(PKG, "R", "api.R"))
ok <- tryCatch({ parse(file.path(PKG, "R", "api.R")); TRUE }, error = function(e) {
  cat("api.R PARSE ERROR:", conditionMessage(e), "\n"); FALSE })
cat(sprintf("\nR/api.R: %d wrappers, %d lines, parse %s\n",
            length(API), length(api_lines), if (ok) "OK" else "FAIL"))
