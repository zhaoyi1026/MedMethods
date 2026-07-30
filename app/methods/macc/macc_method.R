# =============================================================================
# macc -- multilevel mediation under structured unmeasured confounding
# -----------------------------------------------------------------------------
#   M_t = Z_t A + E_1t,      R_t = Z_t C + M_t B + E_2t
# The two error terms are correlated (correlation delta), which is what
# unmeasured mediator-outcome confounding looks like. In the SINGLE-level model
# delta is not identifiable and must be supplied; the multilevel models borrow
# strength across subjects and ESTIMATE it.
# =============================================================================

macc_app_example <- function() {
  d <- med_fn("macc_example")("twolevel", N = 50L, n.trial = 100L)
  out <- list(dat = d$dat, model.type = "twolevel", truth = d$truth)
  out$preview_ui <- macc_preview(out)
  out
}

macc_app_parse <- function(files, opts) {
  obj <- files$dat
  # a list of per-subject tables -> multilevel; a single table -> single level
  is_multi <- (is.list(obj) && !is.data.frame(obj) &&
               all(vapply(obj, function(d) is.data.frame(d) || is.matrix(d),
                          logical(1)))) ||
              (is.data.frame(obj) &&
               any(c("id", "ID", "subject", "Subject", "Sub", "sub") %in% names(obj)))
  dat <- if (is_multi) as_dat_list(obj, "data") else as_dat_df(obj, "data")
  out <- list(dat = dat, model.type = if (is_multi) "twolevel" else "single",
              truth = NULL)
  out$preview_ui <- macc_preview(out)
  out
}

macc_is_multi <- function(d) is.list(d$dat) && !is.data.frame(d$dat)

macc_preview <- function(d) {
  multi <- macc_is_multi(d)
  ntr <- if (multi) vapply(d$dat, nrow, integer(1)) else nrow(d$dat)
  tagList(
    bslib::layout_columns(
      col_widths = c(4, 4, 4),
      bslib::value_box("Structure", if (multi) "Multilevel" else "Single level",
                       theme = "primary"),
      bslib::value_box("Subjects", if (multi) length(d$dat) else 1L,
                       theme = "secondary"),
      bslib::value_box("Trials", if (multi)
        sprintf("%d-%d per subject", min(ntr), max(ntr)) else ntr,
        theme = "secondary")
    ),
    if (!multi)
      div(class = "alert alert-warning p-2 small mb-0 mt-2",
          bsicons::bs_icon("exclamation-triangle"),
          " With a single level the error correlation is NOT identifiable: you must supply delta below. Leaving it at 0 biases the mediator-to-outcome path toward zero."),
    if (!is.null(d$truth))
      tags$p(class = "small text-success mt-2", bsicons::bs_icon("check-circle"),
             sprintf(" Simulated: A = %.1f, B = %.0f, C = %.1f, delta = %.1f, so AB = %.1f.",
                     d$truth$A, d$truth$B, d$truth$C, d$truth$delta, d$truth$ABp))
  )
}

macc_describe <- function(d) {
  if (macc_is_multi(d))
    sprintf("multilevel: %d subjects, %d trials in total.",
            length(d$dat), sum(vapply(d$dat, nrow, integer(1))))
  else sprintf("single level: %d trials.", nrow(d$dat))
}

macc_run <- function(d, params) {
  multi <- macc_is_multi(d)
  mt <- if (multi) (params$model.type %||% "twolevel") else "single"
  est_delta <- isTRUE(params$estimate_delta) && multi
  delta <- if (est_delta) NULL else (params$delta %||% 0)
  fit <- med_fn("macc")(d$dat, model.type = mt,
                        method = params$method %||% "HL",
                        delta = delta,
                        optimizer = "bobyqa",     # avoids the optimx dependency
                        conf.level = params$conf.level %||% 0.95)
  list(fit = fit, d = d, model.type = mt, estimated_delta = est_delta,
       params = params)
}

macc_summarize <- function(res) {
  fit <- res$fit; out <- list()
  cf <- as.data.frame(fit$Coefficients)
  out[["Coefficients"]] <- cbind(Term = rownames(fit$Coefficients), cf)

  out[["Model"]] <- data.frame(
    Quantity = c("model type", "delta (error correlation)",
                 "delta source", "method"),
    Value = c(res$model.type,
              format(round(as.numeric(fit$delta %||% res$params$delta %||% NA), 4)),
              if (res$estimated_delta) "estimated from the multilevel structure"
              else "supplied by you",
              res$params$method %||% "HL"),
    check.names = FALSE)

  if (!is.null(res$d$truth)) {
    tr <- res$d$truth
    rn <- rownames(fit$Coefficients)
    ab_row <- intersect(c("AB.prod", "ABp"), rn)[1]
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

macc_plots <- function(res) {
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
    for (nm in names(tv)) {
      tv[nm] <- switch(nm, A = truth$A, B = truth$B, C = truth$C,
                       C2 = truth$C2, ABp = truth$ABp, ABd = truth$ABd,
                       AB.prod = truth$ABp, AB.diff = truth$ABd, NA_real_)
    }
    p <- plotly::add_trace(p, x = cf$Term, y = as.numeric(tv), type = "scatter",
                           mode = "markers", name = "truth",
                           marker = list(color = "grey35", size = 9,
                                         symbol = "circle-open"))
  }
  plots[["Coefficients"]] <- list(
    plot = p |> plotly::layout(title = "Mediation coefficients",
                               yaxis = list(title = "estimate")),
    data = cbind(Term = rownames(fit$Coefficients), cf[, setdiff(names(cf), "Term")]))

  # trial-level scatter of the mediator against the outcome
  dd <- if (macc_is_multi(res$d)) do.call(rbind, res$d$dat) else res$d$dat
  sdf <- data.frame(M = dd$M, R = dd$R, Z = factor(dd$Z))
  if (nrow(sdf) > 4000) sdf <- sdf[sample.int(nrow(sdf), 4000), ]
  plots[["Mediator vs outcome"]] <- list(
    plot = plotly::plot_ly(sdf, x = ~M, y = ~R, color = ~Z, type = "scatter",
                           mode = "markers",
                           marker = list(size = 5, opacity = 0.5)) |>
      plotly::layout(title = "Trial-level mediator vs outcome, by treatment",
                     xaxis = list(title = "mediator M"),
                     yaxis = list(title = "outcome R")),
    data = sdf)
  plots
}

register_method(list(
  id = "macc",
  name = "macc",
  full_name = "Multilevel Mediation Analysis under Structured Unmeasured Confounding",
  short = "The mediator and outcome errors are correlated, which is what unmeasured mediator-outcome confounding looks like. With several subjects that correlation becomes identifiable and is estimated rather than assumed.",
  status = "ready",
  tags = c("multilevel", "unmeasured confounding", "correlated errors"),
  paper = list(
    citation = "Zhao, Y., & Luo, X. (2023). Multilevel mediation analysis with structured unmeasured mediator-outcome confounding. Computational Statistics & Data Analysis, 179, 107623.",
    url = "https://doi.org/10.1016/j.csda.2022.107623"),
  explain = file.path(APP_DIR, "methods", "macc", "explain.md"),
  example_note = paste("Two-level data: 50 subjects with about 100 trials each,",
                       "true A = 0.5, B = -1, C = 0.5 and error correlation",
                       "delta = 0.5 (so the true indirect effect is -0.5)."),
  data_inputs = list(
    list(id = "dat", label = "Trial-level data (Z, M, R)",
         help = "A table with columns Z (treatment), M (mediator) and R or Y (outcome). Add an 'id' column to fit the multilevel model -- one block of trials per subject. An .rds holding a list of per-subject tables also works.")
  ),
  params = list(
    list(id = "model.type", label = "Model", type = "select",
         choices = c("Two-level" = "twolevel", "Three-level (sessions)" = "multilevel"),
         default = "twolevel",
         help = "Ignored for single-level data. Three-level expects nested sessions within subjects."),
    list(id = "method", label = "Estimator", type = "select",
         choices = c("Hierarchical likelihood (HL)" = "HL", "Two-stage (TS)" = "TS",
                     "HL then TS" = "HL-TS"),
         default = "HL"),
    list(id = "estimate_delta", label = "Estimate the error correlation delta",
         type = "checkbox", default = TRUE,
         help = "The point of the multilevel model. Uncheck to fix delta at the value below instead (required for single-level data)."),
    list(id = "delta", label = "delta, if supplied", type = "numeric",
         default = 0.5, min = -0.95, max = 0.95, step = 0.05,
         help = "Only used when the box above is unchecked. Note that 0 assumes no confounding and biases the mediator-to-outcome path toward zero."),
    list(id = "conf.level", label = "Confidence level", type = "numeric",
         default = 0.95, min = 0.5, max = 0.999, step = 0.01)
  ),
  example = macc_app_example,
  export_example = function(d) {
    dd <- d$dat
    long <- do.call(rbind, lapply(seq_along(dd), function(i)
      cbind(id = names(dd)[i] %||% i, dd[[i]])))
    list(dat = long)
  },
  parse = macc_app_parse,
  describe_data = macc_describe,
  run = macc_run,
  summarize = macc_summarize,
  plots = macc_plots
))
