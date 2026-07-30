# =============================================================================
# HDMediation -- mediation with high-dimensional exposures AND mediators
# -----------------------------------------------------------------------------
#   M = X alpha + E,        Y = X gamma + M beta + e
# A pathway-lasso style penalty on mu = alpha diag(beta) selects individual
# exposure -> mediator -> outcome paths out of the q x p grid.
# =============================================================================

hdmed_app_example <- function() {
  d <- med_fn("hdmediation_example")()      # n = 100, r = 20, p = 20
  out <- list(X = d$X, M = d$M, Y = as.numeric(d$Y), truth = d$truth,
              x_names = colnames(d$X) %||% paste0("X", seq_len(ncol(d$X))),
              m_names = colnames(d$M) %||% paste0("M", seq_len(ncol(d$M))))
  out$preview_ui <- hdmed_preview(out)
  out
}

hdmed_app_parse <- function(files, opts) {
  X <- as_num_matrix(files$X, "X")
  n <- nrow(X)
  M <- as_num_matrix(files$M, "M", n = n)
  Y <- as_num_vector(files$Y, "Y", n = n)
  out <- list(X = X, M = M, Y = Y, truth = NULL,
              x_names = colnames(X), m_names = colnames(M))
  out$preview_ui <- hdmed_preview(out)
  out
}

hdmed_preview <- function(d) {
  n <- nrow(d$X); q <- ncol(d$X)
  tagList(
    bslib::layout_columns(
      col_widths = c(3, 3, 3, 3),
      bslib::value_box("Subjects", n, theme = "primary"),
      bslib::value_box("Exposures", q, theme = "secondary"),
      bslib::value_box("Mediators", ncol(d$M), theme = "secondary"),
      bslib::value_box("Candidate paths", q * ncol(d$M), theme = "secondary")
    ),
    if (n <= q)
      div(class = "alert alert-warning p-2 small mb-0 mt-2",
          bsicons::bs_icon("exclamation-triangle"),
          sprintf(" There are %d exposures but only %d subjects. The direct estimator regresses on all exposures at once and needs n > number of exposures -- switch on the principal-component option below.", q, n)),
    tags$p(class = "small text-muted mt-2",
           "Each (exposure, mediator) pair is one candidate pathway; the penalty selects among them."),
    if (!is.null(d$truth))
      tags$p(class = "small text-success", bsicons::bs_icon("check-circle"),
             " Simulated: three true pathways (exposure 1 -> mediator 1, 2 -> 2, 3 -> 3).")
  )
}

hdmed_describe <- function(d)
  sprintf("%d subjects; %d exposures; %d mediators (%d candidate pathways).",
          nrow(d$X), ncol(d$X), ncol(d$M), ncol(d$X) * ncol(d$M))

hdmed_run <- function(d, params) {
  use_pca <- isTRUE(params$use_pca)
  args <- list(d$X, d$M, d$Y,
               lambda = params$lambda %||% 2.5,
               pi = params$pi %||% 0.5,
               phi = params$phi %||% 2,
               delta = params$delta %||% 0.5,
               max.itr = params$max.itr %||% 5000, tol = 1e-8)
  fit <- if (use_pca)
    do.call(med_fn("hdmediation_pca"),
            c(args, list(adaptive = TRUE, var.prop = params$var.prop %||% 0.9)))
  else
    do.call(med_fn("hdmediation"), args)
  list(fit = fit, d = d, use_pca = use_pca, params = params)
}

hdmed_summarize <- function(res) {
  fit <- res$fit; out <- list()
  IE <- as.matrix(fit$IE)
  xn <- if (res$use_pca) paste0("PC", seq_len(nrow(IE))) else res$d$x_names
  mn <- res$d$m_names
  thr <- 1e-3

  sel <- which(abs(IE) > thr, arr.ind = TRUE)
  out[["Selected pathways"]] <- if (nrow(sel))
    data.frame(Exposure = xn[sel[, 1]], Mediator = mn[sel[, 2]],
               IE = IE[sel], a_path = as.matrix(fit$alpha)[sel],
               b_path = as.numeric(fit$beta)[sel[, 2]], check.names = FALSE)
  else data.frame(Note = "No pathway survived at this lambda -- try a smaller value.")

  out[["Fit summary"]] <- data.frame(
    Quantity = c("lambda", "non-zero pathways", "candidate pathways",
                 "converged", "objective"),
    Value = c(format(fit$lambda), nrow(sel), length(IE),
              as.character(isTRUE(fit$converge)),
              format(round(as.numeric(fit$logLik$obj %||% NA), 3))),
    check.names = FALSE)

  out[["b-path (mediator -> outcome)"]] <- data.frame(
    Mediator = mn, beta = as.numeric(fit$beta), check.names = FALSE)
  out[["direct effect (gamma)"]] <- data.frame(
    Exposure = xn, gamma = as.numeric(fit$gamma), check.names = FALSE)

  if (res$use_pca)
    out[["Principal components"]] <- data.frame(
      Quantity = c("components kept", "variance explained"),
      Value = c(fit$n_pc %||% fit$n.pc %||% NA,
                round(as.numeric(fit$var_explained %||% NA), 3)),
      check.names = FALSE)

  if (!is.null(res$d$truth) && !res$use_pca) {
    tr <- res$d$truth
    tp <- tr$signal
    found <- vapply(seq_len(nrow(tp)), function(i)
      abs(IE[tp[i, 1], tp[i, 2]]) > thr, logical(1))
    out[["Truth vs estimate"]] <- data.frame(
      Pathway = paste0(xn[tp[, 1]], " -> ", mn[tp[, 2]]),
      True_IE = tr$IE[tp], Estimated_IE = IE[tp],
      Selected = ifelse(found, "yes", "no"), check.names = FALSE)
    out[["Selection summary"]] <- data.frame(
      Quantity = c("true pathways", "true pathways selected", "total selected",
                   "candidate pathways", "note"),
      Value = c(nrow(tp), sum(found), nrow(sel), length(IE),
                "This penalty shrinks the effect MAGNITUDES hard: read the table as pathway selection, not as effect-size estimation. Lower lambda for less shrinkage (and less sparsity)."),
      check.names = FALSE)
  }
  out
}

hdmed_plots <- function(res) {
  fit <- res$fit; plots <- list()
  IE <- as.matrix(fit$IE)
  xn <- if (res$use_pca) paste0("PC", seq_len(nrow(IE))) else res$d$x_names
  mn <- res$d$m_names

  plots[["Pathway effect heatmap"]] <- list(
    plot = plotly::plot_ly(z = IE, x = mn, y = xn, type = "heatmap",
                           colorscale = "RdBu", zmid = 0,
                           colorbar = list(title = "IE")) |>
      plotly::layout(title = "Exposure x mediator pathway effects",
                     xaxis = list(title = "mediator"),
                     yaxis = list(title = "exposure")),
    data = data.frame(Exposure = rep(xn, times = ncol(IE)),
                      Mediator = rep(mn, each = nrow(IE)),
                      IE = as.numeric(IE)))

  bdf <- data.frame(Mediator = factor(mn, levels = mn),
                    beta = as.numeric(fit$beta))
  plots[["b-path"]] <- list(
    plot = plotly::plot_ly(bdf, x = ~Mediator, y = ~beta, type = "bar",
                           marker = list(color = "#2c7fb8")) |>
      plotly::layout(title = "Mediator -> outcome coefficients",
                     yaxis = list(title = "beta")),
    data = bdf)

  gdf <- data.frame(Exposure = factor(xn, levels = xn),
                    gamma = as.numeric(fit$gamma))
  plots[["Direct effects"]] <- list(
    plot = plotly::plot_ly(gdf, x = ~Exposure, y = ~gamma, type = "bar",
                           marker = list(color = "#c0392b")) |>
      plotly::layout(title = "Direct effect of each exposure",
                     yaxis = list(title = "gamma")),
    data = gdf)
  plots
}

register_method(list(
  id = "hdmediation",
  name = "HDMediation",
  full_name = "Mediation Analysis with High-Dimensional Exposures and Mediators",
  short = "Both sides are high-dimensional: many exposures AND many mediators. Three penalty terms select individual mediators, individual exposures, and direct effects out of the full q x p grid of candidate pathways.",
  status = "ready",
  tags = c("high-dimensional exposures", "high-dimensional mediators",
           "sparse group lasso", "ADMM"),
  paper = list(
    citation = "Zhao, Y., Li, L., & Alzheimer's Disease Neuroimaging Initiative (2022). Multimodal data integration via mediation analysis with high-dimensional exposures and mediators. Human Brain Mapping, 43(8), 2519-2533.",
    url = "https://doi.org/10.1002/hbm.25800"),
  explain = file.path(APP_DIR, "methods", "hdmediation", "explain.md"),
  example_note = paste("100 subjects, 20 exposures and 20 mediators (400",
                       "candidate pathways) of which three carry signal."),
  data_inputs = list(
    list(id = "X", label = "Exposures (n x q)",
         help = "An n x q numeric matrix. The direct estimator needs n > q; for wider data switch on the principal-component option."),
    list(id = "M", label = "Mediators (n x p)",
         help = "An n x p numeric matrix in the same subject order."),
    list(id = "Y", label = "Outcome (n values)",
         help = "One numeric value per subject, in the same order.")
  ),
  params = list(
    list(id = "lambda", label = "Penalty (lambda)", type = "numeric",
         default = 2.5, min = 0, max = 50, step = 0.25,
         help = "Overall strength. It is split across the three penalty terms as lambda1 = pi*lambda (mediator selection), lambda2 = (1-pi)*lambda (exposure selection) and lambda3 = lambda (direct effects). The dense-to-empty transition is abrupt: on the built-in example lambda = 2 keeps 50 paths, 2.5 keeps 22, and 5 keeps none. Scan and select by BIC on real data."),
    list(id = "pi", label = "Mediator vs exposure selection (pi)", type = "numeric",
         default = 0.5, min = 0, max = 1, step = 0.05,
         help = "Splits lambda between R1 and R2. pi = 1 puts all the weight on R1, which selects individual MEDIATORS by shrinking all paths through a mediator together; pi = 0 puts it all on the group term R2, which selects individual EXPOSURES by shrinking all paths out of an exposure together."),
    list(id = "phi", label = "Convexity weight (c0)", type = "numeric",
         default = 2, min = 0.5, max = 10, step = 0.5,
         help = "Weight on the c0(alpha^2 + beta^2) term inside R1. |alpha*beta| alone is not convex; the sum is convex when c0 >= 1/2, which is why the minimum here is 0.5. The paper fixes c0 = 2."),
    list(id = "delta", label = "Extra l1 on alpha / beta (c1)", type = "numeric",
         default = 0.5, min = 0, max = 5, step = 0.1,
         help = "Weight on the plain lasso term over the individual path coefficients alpha_jk and beta_k."),
    list(id = "use_pca", label = "Reduce exposures to principal components first",
         type = "checkbox", default = FALSE,
         help = "Required when the number of exposures approaches the sample size. Results are then reported per component."),
    list(id = "var.prop", label = "Variance retained by the components",
         type = "numeric", default = 0.9, min = 0.1, max = 0.999, step = 0.05,
         help = "Only used with the principal-component option."),
    list(id = "max.itr", label = "Max ADMM iterations", type = "integer",
         default = 5000, min = 100, max = 20000)
  ),
  example = hdmed_app_example,
  export_example = function(d) list(X = matrix_to_df(d$X), M = matrix_to_df(d$M),
                                    Y = data.frame(Y = d$Y)),
  parse = hdmed_app_parse,
  describe_data = hdmed_describe,
  run = hdmed_run,
  summarize = hdmed_summarize,
  plots = hdmed_plots
))
