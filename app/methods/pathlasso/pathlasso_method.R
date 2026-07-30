# =============================================================================
# Pathway Lasso -- pathway estimation and selection with high-dimensional mediators
# -----------------------------------------------------------------------------
#   M = Z A + E,        Y = Z C + M B + e
# The penalty acts on the PRODUCTS A_j B_j (the pathway effects) rather than on
# A and B separately, so whole pathways are selected or dropped together.
# =============================================================================

pathlasso_app_example <- function() {
  d <- med_fn("pathlasso_example")()        # n = 200, k = 50, signal 1:4
  out <- list(X = as.matrix(d$X), M = d$M, Y = as.numeric(d$Y), truth = d$truth,
              m_names = colnames(d$M) %||% paste0("M", seq_len(ncol(d$M))))
  out$preview_ui <- pathlasso_preview(out)
  out
}

pathlasso_app_parse <- function(files, opts) {
  X <- as_num_matrix(files$X, "X")
  if (ncol(X) != 1L)
    stop("Pathway Lasso takes a single exposure: 'X' has ", ncol(X), " columns.")
  n <- nrow(X)
  M <- as_num_matrix(files$M, "M", n = n)
  Y <- as_num_vector(files$Y, "Y", n = n)
  out <- list(X = X, M = M, Y = Y, truth = NULL, m_names = colnames(M))
  out$preview_ui <- pathlasso_preview(out)
  out
}

pathlasso_preview <- function(d) {
  tagList(
    bslib::layout_columns(
      col_widths = c(4, 4, 4),
      bslib::value_box("Subjects", nrow(d$M), theme = "primary"),
      bslib::value_box("Candidate mediators", ncol(d$M), theme = "secondary"),
      bslib::value_box("Exposure", "1", theme = "secondary")
    ),
    tags$p(class = "small text-muted mt-2",
           "Each mediator gets an a-path (exposure to mediator) and a b-path (mediator to outcome); the penalty selects their product."),
    if (!is.null(d$truth))
      tags$p(class = "small text-success", bsicons::bs_icon("check-circle"),
             sprintf(" Simulated: mediators %s carry signal (true products +-4); the rest are noise.",
                     paste(d$truth$signal, collapse = ", ")))
  )
}

pathlasso_describe <- function(d)
  sprintf("%d subjects; %d candidate mediators; single exposure.",
          nrow(d$M), ncol(d$M))

pathlasso_run <- function(d, params) {
  lam <- params$lambda %||% 0.001
  tuned <- NULL
  if (isTRUE(params$tune)) {
    grid <- 10^seq(-4, -1, length.out = 4)
    tuned <- tryCatch(
      med_fn("pathlasso_ksc")(d$X, d$M, d$Y, lambda = grid,
                              n.rep = params$n.rep %||% 3,
                              omega = params$omega %||% 0,
                              phi = params$phi %||% 1,
                              max.itr = params$max.itr %||% 3000, tol = 1e-8),
      error = function(e) NULL)
    if (!is.null(tuned) && !is.null(tuned$lambda.est)) lam <- tuned$lambda.est
    if (!is.null(tuned)) tuned$grid <- grid      # the return carries no grid itself
  }
  fit <- med_fn("pathlasso")(d$X, d$M, d$Y, lambda = lam,
                             omega = params$omega %||% 0,
                             phi = params$phi %||% 1,
                             max.itr = params$max.itr %||% 3000, tol = 1e-8)
  list(fit = fit, tuned = tuned, lambda = lam, d = d, params = params)
}

pathlasso_summarize <- function(res) {
  fit <- res$fit; out <- list()
  A <- as.numeric(fit$A); B <- as.numeric(fit$B); AB <- as.numeric(fit$AB)
  thr <- 0.05 * max(abs(AB), na.rm = TRUE)
  sel <- which(abs(AB) > thr)

  out[["Pathway effects"]] <- data.frame(
    Mediator = res$d$m_names, a_path = A, b_path = B, product_AB = AB,
    Selected = ifelse(abs(AB) > thr, "yes", ""), check.names = FALSE)

  out[["Selected pathways"]] <- if (length(sel))
    data.frame(Mediator = res$d$m_names[sel], a_path = A[sel], b_path = B[sel],
               product_AB = AB[sel], check.names = FALSE)
  else data.frame(Note = "No pathway exceeded the selection threshold at this lambda.")

  out[["Fit summary"]] <- data.frame(
    Quantity = c("lambda used", "direct effect C", "non-zero pathways",
                 "converged", "max |AB|"),
    Value = c(format(res$lambda), format(round(as.numeric(fit$C)[1], 4)),
              length(sel), as.character(isTRUE(fit$converge)),
              format(round(max(abs(AB)), 4))),
    check.names = FALSE)

  if (!is.null(fit$constraint)) {
    cst <- fit$constraint
    out[["ADMM constraints"]] <- if (is.data.frame(cst))
      cbind(Constraint = rownames(cst), cst) else
      data.frame(Constraint = names(cst), Value = unlist(cst, use.names = FALSE))
  }

  if (!is.null(res$tuned)) {
    g <- res$tuned$grid
    k <- if (!is.null(res$tuned$vss)) as.numeric(res$tuned$vss) else rep(NA_real_, length(g))
    out[["Selection stability (kappa)"]] <- data.frame(
      lambda = g, mean_kappa = k[seq_along(g)], check.names = FALSE)
  }

  if (!is.null(res$d$truth)) {
    tr <- res$d$truth
    tp <- intersect(sel, tr$signal)
    out[["Truth vs estimate"]] <- data.frame(
      Quantity = c("true signal mediators", "selected mediators",
                   "true positives", "false positives",
                   "true AB on signal", "estimated AB on signal"),
      Value = c(paste(tr$signal, collapse = ", "),
                if (length(sel)) paste(sel, collapse = ", ") else "(none)",
                length(tp), length(setdiff(sel, tr$signal)),
                paste(round(tr$AB[tr$signal], 2), collapse = ", "),
                paste(round(AB[tr$signal], 2), collapse = ", ")),
      check.names = FALSE)
  }
  out
}

pathlasso_plots <- function(res) {
  fit <- res$fit; plots <- list()
  A <- as.numeric(fit$A); B <- as.numeric(fit$B); AB <- as.numeric(fit$AB)
  idx <- seq_along(AB)
  truth_ab <- if (!is.null(res$d$truth)) res$d$truth$AB else NULL

  pdf_ <- data.frame(Mediator = idx, product_AB = AB)
  p <- plotly::plot_ly(pdf_, x = ~Mediator, y = ~product_AB, type = "bar",
                       name = "estimate", marker = list(color = "#c0392b"))
  if (!is.null(truth_ab))
    p <- plotly::add_trace(p, x = idx, y = truth_ab, type = "scatter",
                           mode = "markers", name = "truth",
                           marker = list(color = "grey35", size = 7,
                                         symbol = "circle-open"))
  plots[["Pathway effects (AB)"]] <- list(
    plot = p |> plotly::layout(title = "Pathway effect A_j x B_j per mediator",
                               xaxis = list(title = "mediator"),
                               yaxis = list(title = "A x B")),
    data = if (!is.null(truth_ab))
      cbind(pdf_, true_AB = truth_ab) else pdf_)

  abdf <- data.frame(Mediator = idx, a_path = A, b_path = B)
  plots[["a-path and b-path"]] <- list(
    plot = plotly::plot_ly(abdf, x = ~Mediator, y = ~a_path, type = "bar",
                           name = "a-path") |>
      plotly::add_trace(y = ~b_path, name = "b-path") |>
      plotly::layout(title = "Estimated a- and b-paths", barmode = "group",
                     xaxis = list(title = "mediator"),
                     yaxis = list(title = "coefficient")),
    data = abdf)

  if (!is.null(res$tuned) && !is.null(res$tuned$vss)) {
    g <- res$tuned$grid
    kdf <- data.frame(lambda = g,
                      mean_kappa = as.numeric(res$tuned$vss)[seq_along(g)])
    plots[["Selection stability"]] <- list(
      plot = plotly::plot_ly(kdf, x = ~lambda, y = ~mean_kappa,
                             type = "scatter", mode = "lines+markers") |>
        plotly::layout(title = "Selection stability (Cohen's kappa) across lambda",
                       xaxis = list(title = "lambda", type = "log"),
                       yaxis = list(title = "mean kappa")),
      data = kdf)
  }
  plots
}

register_method(list(
  id = "pathlasso",
  name = "Pathway Lasso",
  full_name = "Pathway Lasso: Pathway Estimation and Selection with High-Dimensional Mediators",
  short = "Penalises the PRODUCT of each mediator's a-path and b-path, so whole exposure -> mediator -> outcome pathways are selected or dropped as units.",
  status = "ready",
  tags = c("high-dimensional mediators", "pathway selection", "ADMM"),
  paper = list(
    citation = "Zhao, Y., & Luo, X. (2022). Pathway Lasso: pathway estimation and selection with high-dimensional mediators. Statistics and Its Interface, 15(1), 39-50.",
    url = NULL),
  explain = file.path(APP_DIR, "methods", "pathlasso", "explain.md"),
  example_note = paste("200 subjects and 50 candidate mediators, of which four",
                       "carry signal with products of mixed sign (true +-4)."),
  data_inputs = list(
    list(id = "X", label = "Exposure (n values)",
         help = "A single exposure, one value per subject."),
    list(id = "M", label = "Mediators (n x k)",
         help = "An n x k numeric matrix of candidate mediators, same subject order as the exposure."),
    list(id = "Y", label = "Outcome (n values)",
         help = "One numeric value per subject, in the same order.")
  ),
  params = list(
    list(id = "lambda", label = "Penalty (lambda)", type = "numeric",
         default = 0.001, min = 0, max = 10, step = 0.0005,
         help = "Weight on the pathway term sum_j (|A_j B_j| + phi(A_j^2 + B_j^2)) plus |C|. On standardized data the useful range is small -- around 1e-3 here. The jump from dense to empty solutions is abrupt, so scan rather than guess."),
    list(id = "tune", label = "Pick lambda by selection stability",
         type = "checkbox", default = FALSE,
         help = "Scans 1e-4 to 1e-1 and keeps the most parsimonious stable value. Slower (it refits on split halves). Check the kappa table -- near-zero kappa means the criterion is uninformative for this sample size."),
    list(id = "n.rep", label = "Stability splits", type = "integer",
         default = 3, min = 2, max = 10,
         help = "Only used when tuning by selection stability."),
    list(id = "omega", label = "Extra l1 on the paths (omega)", type = "numeric",
         default = 0, min = 0, max = 5, step = 0.05,
         help = "Weight on sum_j (|A_j| + |B_j|) -- separate shrinkage on the individual paths, in the spirit of the elastic net. Note that even small values can zero out the whole solution."),
    list(id = "phi", label = "Ridge weight (phi)", type = "numeric",
         default = 1, min = 0.5, max = 10, step = 0.1,
         help = "Weight on the phi(A_j^2 + B_j^2) term. |A_j B_j| alone is not convex; |ab| + phi(a^2 + b^2) is convex if and only if phi >= 1/2, which is why the minimum here is 0.5."),
    list(id = "max.itr", label = "Max ADMM iterations", type = "integer",
         default = 3000, min = 100, max = 20000)
  ),
  example = pathlasso_app_example,
  export_example = function(d) list(X = data.frame(X = d$X[, 1]),
                                    M = matrix_to_df(d$M),
                                    Y = data.frame(Y = d$Y)),
  parse = pathlasso_app_parse,
  describe_data = pathlasso_describe,
  run = pathlasso_run,
  summarize = pathlasso_summarize,
  plots = pathlasso_plots
))
