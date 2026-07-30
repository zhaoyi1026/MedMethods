# =============================================================================
# Data I/O helpers
# -----------------------------------------------------------------------------
# Translate what a user uploads (CSV / .rds / .RData) into the structures the
# mediation methods expect. The methods here use five shapes:
#
#   matrix   n x p numeric              spcma, pathlasso, hdmediation, pcma,
#                                       pathlasso2b, cfma (N x ntp), hetermed Z
#   vector   length n numeric           the outcome Y, and hetermed's X / M
#   list     n matrices of T_i x p      gmed's covariance mediator
#   dat      data.frame Z, M, R         macc / gma single-level trial data
#   dat list list of such data.frames   macc two-/three-level data
#
# Everything funnels through read_upload(); the coercers below are shared by the
# method plugins so each one only declares which shape it wants.
# =============================================================================

#' Read any uploaded file into an R object.
#'
#' - .rds               -> the stored object
#' - .RData / .rda      -> the single stored object, or a named list if several
#' - .csv / .tsv / .txt -> a data.frame
read_upload <- function(path, name = path) {
  ext <- tolower(tools::file_ext(name))
  if (ext == "rds") return(readRDS(path))
  if (ext %in% c("rdata", "rda")) {
    e <- new.env()
    load(path, envir = e)
    objs <- ls(e)
    if (length(objs) == 1) return(get(objs[1], envir = e))
    return(mget(objs, envir = e))
  }
  read_table_any(path, name)
}

#' Read an uploaded delimited file into a data.frame, guessing the separator.
read_table_any <- function(path, name = path) {
  ext <- tolower(tools::file_ext(name))
  sep <- if (ext %in% c("tsv", "txt")) "\t" else ","
  df <- utils::read.csv(path, sep = sep, header = TRUE,
                        stringsAsFactors = FALSE, check.names = FALSE)
  if (ncol(df) < 2) {
    alt <- if (sep == ",") "\t" else ","
    df2 <- utils::read.csv(path, sep = alt, header = TRUE,
                           stringsAsFactors = FALSE, check.names = FALSE)
    if (ncol(df2) > ncol(df)) df <- df2
  }
  df
}

#' Unwrap a named .RData container, looking for one of `keys`.
unwrap_named <- function(obj, keys, what) {
  if (is.list(obj) && !is.data.frame(obj) && !is.matrix(obj)) {
    looks_like_matrix_list <- length(obj) > 0 &&
      all(vapply(obj, function(m) is.matrix(m) || is.data.frame(m), logical(1)))
    if (!looks_like_matrix_list) {
      hit <- intersect(keys, names(obj))
      if (length(hit) == 0)
        stop(sprintf("Could not find %s in the uploaded file (looked for: %s).",
                     what, paste(keys, collapse = ", ")))
      return(obj[[hit[1]]])
    }
  }
  obj
}

#' Is this column plausibly a subject/row identifier rather than data?
.looks_like_id <- function(x) {
  if (is.character(x) || is.factor(x)) return(TRUE)
  v <- suppressWarnings(as.numeric(x))
  if (anyNA(v)) return(TRUE)
  # a 1..n sequence is an index, not a covariate
  all(v == seq_along(v))
}

#' Coerce an upload into an `n x p` numeric matrix.
#'
#' Accepts a matrix, a data.frame (an obvious leading id column is dropped), or
#' a named .RData container holding one under a recognised key.
#' @param n expected number of rows; NULL to skip the check.
#' @param drop_id drop a leading identifier column if one is detected.
as_num_matrix <- function(obj, name = "data", n = NULL, drop_id = TRUE,
                          keys = NULL) {
  obj <- unwrap_named(obj, keys %||% c(name, toupper(name), tolower(name)), name)
  if (is.data.frame(obj)) {
    if (drop_id && ncol(obj) > 1 && .looks_like_id(obj[[1]])) obj <- obj[, -1, drop = FALSE]
    nm <- names(obj)
    obj <- vapply(obj, function(x) suppressWarnings(as.numeric(x)), numeric(nrow(obj)))
    obj <- matrix(obj, ncol = length(nm), dimnames = list(NULL, nm))
  }
  if (is.vector(obj) && is.numeric(obj)) obj <- matrix(obj, ncol = 1)
  M <- as.matrix(obj)
  storage.mode(M) <- "double"
  if (anyNA(M))
    stop(sprintf("'%s' contains non-numeric or missing values.", name))
  if (!is.null(n) && nrow(M) != n)
    stop(sprintf("'%s' has %d rows but %d were expected.", name, nrow(M), n))
  if (is.null(colnames(M)))
    colnames(M) <- if (ncol(M) == 1) name else paste0(name, seq_len(ncol(M)))
  M
}

#' Coerce an upload into a length-`n` numeric vector.
as_num_vector <- function(obj, name = "Y", n = NULL, keys = NULL) {
  obj <- unwrap_named(obj, keys %||% c(name, toupper(name), tolower(name)), name)
  if (is.data.frame(obj)) {
    # if there are two columns and the first looks like an id, take the second
    obj <- if (ncol(obj) >= 2 && .looks_like_id(obj[[1]])) obj[[2]] else obj[[ncol(obj)]]
  }
  v <- suppressWarnings(as.numeric(unlist(obj, use.names = FALSE)))
  if (anyNA(v))
    stop(sprintf("'%s' contains non-numeric or missing values.", name))
  if (!is.null(n) && length(v) != n)
    stop(sprintf("'%s' has length %d but %d values were expected.",
                 name, length(v), n))
  v
}

#' Coerce an upload into a list of `T_i x p` matrices (gmed's mediator).
#'
#' Accepts a list of matrices, or a long data.frame whose first column is the
#' subject id and whose remaining columns are the p variables.
as_matrix_list <- function(obj, name = "M") {
  obj <- unwrap_named(obj, c("M", "M_list", "Mlist", "mediator"), name)

  if (is.data.frame(obj)) {
    if (ncol(obj) < 2)
      stop(sprintf("'%s' needs a subject-id column plus at least one variable.", name))
    ids_raw <- obj[[1]]
    vals <- obj[, -1, drop = FALSE]
    nmv <- names(vals)
    vals <- vapply(vals, function(x) suppressWarnings(as.numeric(x)),
                   numeric(nrow(vals)))
    vals <- matrix(vals, ncol = length(nmv), dimnames = list(NULL, nmv))
    if (anyNA(vals))
      stop(sprintf("'%s' contains non-numeric or missing values.", name))
    ids <- unique(ids_raw)
    Ml <- lapply(ids, function(i) vals[ids_raw == i, , drop = FALSE])
    names(Ml) <- as.character(ids)
  } else {
    if (!is.list(obj) ||
        !all(vapply(obj, function(m) is.matrix(m) || is.data.frame(m), logical(1))))
      stop(sprintf("'%s' must be a list of T_i x p matrices, or a long CSV (id + p columns).",
                   name))
    Ml <- lapply(obj, function(m) {
      m <- as.matrix(m); storage.mode(m) <- "double"; m
    })
    if (is.null(names(Ml)) || any(names(Ml) == ""))
      names(Ml) <- paste0("S", seq_along(Ml))
  }

  if (length(unique(vapply(Ml, ncol, integer(1)))) != 1L)
    stop(sprintf("All '%s' matrices must have the same number of columns (p).", name))
  Tv <- vapply(Ml, nrow, integer(1))
  if (any(Tv < 2))
    stop(sprintf("Each '%s' matrix needs at least 2 rows to estimate a covariance.", name))
  vn <- colnames(Ml[[1]])
  if (is.null(vn)) vn <- paste0("V", seq_len(ncol(Ml[[1]])))
  list(M = Ml, ids = names(Ml), Tvec = Tv, p = ncol(Ml[[1]]), var_names = vn)
}

#' Coerce an upload into a `Z, M, R` trial-level data.frame (macc / gma).
as_dat_df <- function(obj, name = "data") {
  obj <- unwrap_named(obj, c("dat", "data", "data1"), name)
  if (is.matrix(obj)) obj <- as.data.frame(obj)
  if (!is.data.frame(obj))
    stop(sprintf("'%s' must be a table with columns Z, M and R.", name))
  nm <- names(obj)
  # accept Y as a synonym for R (the outcome)
  if (!"R" %in% nm && "Y" %in% nm) names(obj)[nm == "Y"] <- "R"
  # accept X as a synonym for Z (the treatment)
  nm <- names(obj)
  if (!"Z" %in% nm && "X" %in% nm) names(obj)[nm == "X"] <- "Z"
  need <- c("Z", "M", "R")
  miss <- setdiff(need, names(obj))
  if (length(miss))
    stop(sprintf("'%s' is missing column(s): %s. Expected Z (treatment), M (mediator), R or Y (outcome).",
                 name, paste(miss, collapse = ", ")))
  out <- obj[, need, drop = FALSE]
  out[] <- lapply(out, function(x) suppressWarnings(as.numeric(x)))
  if (anyNA(out))
    stop(sprintf("'%s' contains non-numeric or missing values in Z / M / R.", name))
  out
}

#' Coerce an upload into a list of `Z, M, R` data.frames, one per subject.
#'
#' Accepts a list of such tables, or one long table with a leading subject-id
#' column (`id` / `subject` / `Sub`, else the first column).
as_dat_list <- function(obj, name = "data") {
  obj <- unwrap_named(obj, c("dat", "data", "data2", "data3"), name)
  if (is.list(obj) && !is.data.frame(obj) &&
      all(vapply(obj, function(d) is.data.frame(d) || is.matrix(d), logical(1)))) {
    return(lapply(obj, as_dat_df, name = name))
  }
  if (!is.data.frame(obj))
    stop(sprintf("'%s' must be a list of per-subject tables, or one long table with a subject-id column.",
                 name))
  idcol <- intersect(c("id", "ID", "subject", "Subject", "Sub", "sub"), names(obj))
  idcol <- if (length(idcol)) idcol[1] else names(obj)[1]
  ids <- obj[[idcol]]
  rest <- obj[, setdiff(names(obj), idcol), drop = FALSE]
  out <- lapply(unique(ids), function(i) as_dat_df(rest[ids == i, , drop = FALSE], name))
  names(out) <- as.character(unique(ids))
  if (length(out) < 2)
    stop(sprintf("'%s' resolved to %d subject(s); the multilevel models need several.",
                 name, length(out)))
  out
}

# ---------------------------------------------------------------------------
# Shared preview builders
# ---------------------------------------------------------------------------

#' A value-box row summarising an X / M / Y dataset.
xmy_preview <- function(d, extra = NULL, note = NULL) {
  boxes <- list(
    bslib::value_box("Subjects", nrow(as.matrix(d$M)), theme = "primary"),
    bslib::value_box("Exposures", if (is.null(dim(d$X))) 1L else ncol(d$X),
                     theme = "secondary"),
    bslib::value_box("Mediators", ncol(as.matrix(d$M)), theme = "secondary")
  )
  if (!is.null(extra)) boxes <- c(boxes, extra)
  w <- rep(floor(12 / length(boxes)), length(boxes))
  tagList(
    do.call(bslib::layout_columns, c(list(col_widths = w), boxes)),
    if (!is.null(note)) tags$p(class = "small text-muted mt-2", note),
    if (!is.null(d$truth))
      tags$p(class = "small text-success mt-1",
             bsicons::bs_icon("check-circle"),
             " Simulated data with known truth — the results include a truth-vs-estimate table.")
  )
}

#' Turn a list of subject matrices into a long data.frame (for CSV export).
matrix_list_to_long <- function(Ml, id_name = "id", var_names = NULL) {
  if (is.null(var_names)) var_names <- colnames(Ml[[1]]) %||%
    paste0("V", seq_len(ncol(Ml[[1]])))
  out <- do.call(rbind, lapply(seq_along(Ml), function(i) {
    m <- as.data.frame(Ml[[i]])
    names(m) <- var_names
    cbind(stats::setNames(data.frame(names(Ml)[i] %||% i), id_name), m)
  }))
  rownames(out) <- NULL
  out
}

#' Turn a matrix into a data.frame with an id column (for CSV export).
matrix_to_df <- function(M, id = NULL, id_name = "id") {
  M <- as.matrix(M)
  df <- as.data.frame(M)
  if (!is.null(id)) df <- cbind(stats::setNames(data.frame(id), id_name), df)
  rownames(df) <- NULL
  df
}

#' Cosine similarity up to sign (projections are identified only up to sign).
abs_cos <- function(a, b) {
  a <- as.numeric(a); b <- as.numeric(b)
  abs(sum(a * b)) / sqrt(sum(a^2) * sum(b^2))
}
