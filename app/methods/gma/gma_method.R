# =============================================================================
# gma -- Granger mediation analysis for time series
# -----------------------------------------------------------------------------
# Mediation for a time series whose errors follow a VAR(p) process, with a
# correlation parameter absorbing unmeasured confounding:
#   M_t = Z_t A + E_1t,   R_t = Z_t C + M_t B + E_2t,   (E_1, E_2) VAR(p)
# As in macc, the SINGLE-level model cannot identify the error correlation delta;
# supply it, or use the two-level model across subjects.
# =============================================================================

# Two-level by default: the error correlation delta is only identifiable across
# several series, so a single-series example could never demonstrate estimating
# it. 120 series of about 150 time points each -- enough for delta and all three
# path coefficients to be recovered (about 15s to fit).
#
# The two-level branch of gma_example() was added after the package's first
# release, so an older installed MedMethods has a gma_example() that does not
# accept `model.type` and would fail with "unused arguments". Rather than depend
# on that, build the data here from gma_sim_two(), which has been exported all
# along; med_gen_two_level() below is the same construction the package uses.
gma_app_example <- function(N = 120L, n.time = 150L) {
  gen <- med_fn("gma_example")
  if ("model.type" %in% names(formals(gen))) {
    d <- gen("twolevel", N = N, n.time = n.time)
    out <- list(dat = d$dat, truth = d$truth)
  } else {
    out <- gma_two_level_fallback(N, n.time)
  }
  out$preview_ui <- gma_preview(out)
  out
}

# Self-contained two-level generator, used when the installed package predates
# gma_example("twolevel"). Identical construction: N series, each its own VAR(1)
# process, series-level coefficients drawn around (A, B, C) with covariance Lambda.
gma_two_level_fallback <- function(N = 120L, n.time = 150L, delta = 0.5) {
  A <- 0.5; B <- -1; C <- 0.5
  Sigma <- matrix(c(1, 2 * delta, 2 * delta, 4), 2, 2)
  Delta <- matrix(c(2, delta * sqrt(2 * 8), delta * sqrt(2 * 8), 8), 2, 2)
  W     <- matrix(c(-0.809, 0.154, -0.618, -0.5), 2, 2)
  Lambda <- diag(0.5, 3)
  set.seed(2000L)
  ni <- matrix(stats::rpois(N, n.time), N, 1)
  set.seed(1000L)
  Z.list <- lapply(seq_len(N), function(i)
    matrix(stats::rbinom(ni[i, 1], size = 1L, prob = 0.5), ni[i, 1], 1L))
  set.seed(1000L)
  sim <- med_fn("gma_sim_two")(Z.list, N, theta = c(A, B, C), Sigma, W,
                               Delta = Delta, p = 1L, Lambda = Lambda,
                               nburn = 500L)
  list(dat = sim$data,
       truth = list(A = A, B = B, C = C, C2 = C + A * B, ABp = A * B,
                    delta = delta, p = 1L, W = W, Lambda = Lambda))
}

gma_app_parse <- function(files, opts) {
  obj <- files$dat
  is_multi <- (is.list(obj) && !is.data.frame(obj) &&
               all(vapply(obj, function(d) is.data.frame(d) || is.matrix(d),
                          logical(1)))) ||
              (is.data.frame(obj) &&
               any(c("id", "ID", "subject", "Subject", "Sub", "sub") %in% names(obj)))
  dat <- if (is_multi) as_dat_list(obj, "data") else as_dat_df(obj, "data")
  out <- list(dat = dat, truth = NULL)
  out$preview_ui <- gma_preview(out)
  out
}

gma_is_multi <- function(d) is.list(d$dat) && !is.data.frame(d$dat)

gma_preview <- function(d) {
  multi <- gma_is_multi(d)
  ntp <- if (multi) vapply(d$dat, nrow, integer(1)) else nrow(d$dat)
  tagList(
    bslib::layout_columns(
      col_widths = c(4, 4, 4),
      bslib::value_box("Structure", if (multi) "Two-level" else "Single series",
                       theme = "primary"),
      bslib::value_box("Series", if (multi) length(d$dat) else 1L,
                       theme = "secondary"),
      bslib::value_box("Time points", if (multi)
        sprintf("%d-%d", min(ntp), max(ntp)) else ntp, theme = "secondary")
    ),
    if (!multi)
      div(class = "alert alert-warning p-2 small mb-0 mt-2",
          bsicons::bs_icon("exclamation-triangle"),
          " With a single series the error correlation is NOT identifiable: supply delta below. At delta = 0 the mediator-to-outcome path is biased toward zero -- on simulated single-series data with a true B = -1, delta = 0 gives about -0.003. Add an id column to fit several series and have delta estimated instead."),
    if (!is.null(d$truth))
      tags$p(class = "small text-success mt-2", bsicons::bs_icon("check-circle"),
             sprintf(" Simulated VAR(1) series: A = %.1f, B = %.0f, C = %.1f, delta = %.1f%s",
                     d$truth$A, d$truth$B, d$truth$C, d$truth$delta,
                     if (multi)
                       " -- and with several series delta can be estimated rather than supplied."
                     else "."))
  )
}

gma_describe <- function(d) {
  if (gma_is_multi(d))
    sprintf("two-level: %d series, %d time points in total.",
            length(d$dat), sum(vapply(d$dat, nrow, integer(1))))
  else sprintf("single series of %d time points.", nrow(d$dat))
}

gma_run <- function(d, params) {
  multi <- gma_is_multi(d)
  mt <- if (multi) "twolevel" else "single"
  est_delta <- isTRUE(params$estimate_delta) && multi
  fit <- med_fn("gma")(d$dat, model.type = mt,
                       method = params$method %||% "HL",
                       delta = if (est_delta) NULL else (params$delta %||% 0),
                       p = params$p %||% 1,
                       single.var.asmp = isTRUE(params$var.asmp),
                       conf.level = params$conf.level %||% 0.95)
  list(fit = fit, d = d, model.type = mt, estimated_delta = est_delta,
       params = params)
}

gma_summarize <- function(res) {
  fit <- res$fit; out <- list()
  cf <- as.data.frame(fit$Coefficients)
  out[["Coefficients"]] <- cbind(Term = rownames(fit$Coefficients), cf)

  out[["Model"]] <- data.frame(
    Quantity = c("model type", "VAR lag p", "delta", "delta source",
                 "variance for the single-level fit"),
    Value = c(res$model.type, res$params$p %||% 1,
              format(round(as.numeric(fit$delta %||% res$params$delta %||% NA), 4)),
              if (res$estimated_delta) "estimated across series"
              else "supplied by you",
              if (isTRUE(res$params$var.asmp)) "asymptotic" else "empirical"),
    check.names = FALSE)

  if (!is.null(fit$W)) {
    W <- as.matrix(fit$W)
    out[["VAR transition matrix (W)"]] <- cbind(
      Row = paste0("r", seq_len(nrow(W))), as.data.frame(W))
  }

  if (!is.null(res$d$truth)) {
    tr <- res$d$truth
    rn <- rownames(fit$Coefficients)
    ab_row <- intersect(c("AB.p", "ABp", "AB.prod"), rn)[1]
    pick <- function(r) if (r %in% rn) fit$Coefficients[r, "Estimate"] else NA_real_
    out[["Truth vs estimate"]] <- data.frame(
      Quantity = c("A", "B", "C", "AB (indirect)", "delta"),
      True = c(tr$A, tr$B, tr$C, tr$ABp, tr$delta),
      Estimate = c(pick("A"), pick("B"), pick("C"),
                   if (!is.na(ab_row)) fit$Coefficients[ab_row, "Estimate"] else NA_real_,
                   as.numeric(fit$delta %||% res$params$delta %||% NA)),
      check.names = FALSE)
  }
  out
}

gma_plots <- function(res) {
  fit <- res$fit; plots <- list()
  cf <- as.data.frame(fit$Coefficients)
  cf$Term <- factor(rownames(fit$Coefficients), levels = rownames(fit$Coefficients))
  has_ci <- all(c("LB", "UB") %in% names(cf))
  p <- plotly::plot_ly(cf, x = ~Term, y = ~Estimate, type = "bar",
                       marker = list(color = "#2c7fb8"),
                       error_y = if (has_ci)
                         list(type = "data", symmetric = FALSE,
                              array = cf$UB - cf$Estimate,
                              arrayminus = cf$Estimate - cf$LB) else NULL)
  truth <- res$d$truth
  if (!is.null(truth)) {
    tv <- setNames(rep(NA_real_, nrow(cf)), rownames(fit$Coefficients))
    for (nm in names(tv))
      tv[nm] <- switch(nm, A = truth$A, B = truth$B, C = truth$C,
                       C2 = truth$C2, AB.p = truth$ABp, AB.d = truth$ABp,
                       NA_real_)
    p <- plotly::add_trace(p, x = cf$Term, y = as.numeric(tv), type = "scatter",
                           mode = "markers", name = "truth",
                           marker = list(color = "grey35", size = 9,
                                         symbol = "circle-open"))
  }
  plots[["Coefficients"]] <- list(
    plot = p |> plotly::layout(title = "Mediation coefficients",
                               yaxis = list(title = "estimate")),
    data = cbind(Term = rownames(fit$Coefficients),
                 cf[, setdiff(names(cf), "Term")]))

  dd <- if (gma_is_multi(res$d)) res$d$dat[[1]] else res$d$dat
  n <- nrow(dd)
  keep <- seq_len(min(n, 400L))
  tdf <- data.frame(t = keep, M = dd$M[keep], R = dd$R[keep], Z = dd$Z[keep])
  plots[["Series"]] <- list(
    plot = plotly::plot_ly(tdf, x = ~t) |>
      plotly::add_lines(y = ~M, name = "mediator M") |>
      plotly::add_lines(y = ~R, name = "outcome R") |>
      plotly::add_lines(y = ~Z, name = "treatment Z",
                        line = list(dash = "dot")) |>
      plotly::layout(title = sprintf("First %d time points", length(keep)),
                     xaxis = list(title = "time"),
                     yaxis = list(title = "value")),
    data = tdf)
  plots
}

register_method(list(
  id = "gma",
  name = "gma",
  full_name = "Granger Mediation Analysis of Time Series",
  short = "Mediation for time series whose errors follow a VAR(p) process, with a correlation parameter absorbing unmeasured confounding between the mediator and outcome equations.",
  status = "ready",
  tags = c("time series", "VAR errors", "unmeasured confounding"),
  paper = list(
    citation = "Zhao, Y., & Luo, X. (2019). Granger mediation analysis of multiple time series with an application to functional magnetic resonance imaging. Biometrics, 75(3), 788-798.",
    url = "https://doi.org/10.1111/biom.13056"),
  explain = file.path(APP_DIR, "methods", "gma", "explain.md"),
  example_note = paste("Two-level: 120 VAR(1) series of about 150 time points each,",
                       "with A = 0.5, B = -1, C = 0.5 and error correlation",
                       "delta = 0.5. Because there are several series, delta is",
                       "identifiable and is estimated rather than supplied."),
  data_inputs = list(
    list(id = "dat", label = "Time series (Z, M, R)",
         help = "A table with columns Z (treatment at each time point), M (mediator) and R or Y (outcome), one row per time point in order. Add an 'id' column -- one block of rows per series -- or upload an .rds list of tables, to fit the two-level model across several series. Only then is the error correlation delta identifiable.")
  ),
  params = list(
    list(id = "p", label = "VAR lag order (p)", type = "integer",
         default = 1, min = 1, max = 6,
         help = "Order of the autoregressive error process: how many previous time points the errors depend on. The first p observations are used to condition on, so a higher order costs degrees of freedom."),
    list(id = "estimate_delta", label = "Estimate the error correlation delta",
         type = "checkbox", default = TRUE,
         help = "Only possible with two-level data (several series), which is what the built-in example provides. For a single series delta is not identifiable and must be supplied below."),
    list(id = "delta", label = "delta, if supplied", type = "numeric",
         default = 0.5, min = -0.95, max = 0.95, step = 0.05,
         help = "For a single series this is required. Note that 0 assumes no confounding and biases the mediator-to-outcome path toward zero."),
    list(id = "method", label = "Two-level estimator", type = "select",
         choices = c("Two-stage (TS)" = "TS", "Hierarchical likelihood (HL)" = "HL",
                     "HL then TS" = "HL-TS"),
         default = "TS",
         help = "Ignored for a single series. TS is the quicker of the two; HL maximises the hierarchical likelihood and is slower."),
    list(id = "var.asmp", label = "Asymptotic variances (single series)",
         type = "checkbox", default = TRUE,
         help = "Uses the VAR companion-matrix asymptotic variance rather than the empirical one."),
    list(id = "conf.level", label = "Confidence level", type = "numeric",
         default = 0.95, min = 0.5, max = 0.999, step = 0.01)
  ),
  example = gma_app_example,
  export_example = function(d) {
    dd <- d$dat
    if (is.data.frame(dd)) return(list(dat = dd))
    long <- do.call(rbind, lapply(seq_along(dd), function(i)
      cbind(id = names(dd)[i] %||% i, dd[[i]])))
    list(dat = long)
  },
  parse = gma_app_parse,
  describe_data = gma_describe,
  run = gma_run,
  summarize = gma_summarize,
  plots = gma_plots
))
