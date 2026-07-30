# =============================================================================
# Multimodal -- two blocks of high-dimensional mediators
# -----------------------------------------------------------------------------
#   M1 = X beta  + e1
#   M2 = X zeta  + M1 Lambda + e2
#   Y  = X delta + M1 theta  + M2 pi + e
# Pathways run X -> M1 -> Y, X -> M2 -> Y, and X -> M1 -> M2 -> Y.
# =============================================================================

pl2b_app_example <- function() {
  d <- med_fn("pathlasso2b_example")()      # n = 200, p1 = 20, p2 = 30
  out <- list(X = as.numeric(d$X), M1 = d$M1, M2 = d$M2, Y = as.numeric(d$Y),
              truth = d$truth,
              m1_names = colnames(d$M1) %||% paste0("A", seq_len(ncol(d$M1))),
              m2_names = colnames(d$M2) %||% paste0("B", seq_len(ncol(d$M2))))
  out$preview_ui <- pl2b_preview(out)
  out
}

pl2b_app_parse <- function(files, opts) {
  X <- as_num_vector(files$X, "X")
  n <- length(X)
  M1 <- as_num_matrix(files$M1, "M1", n = n)
  M2 <- as_num_matrix(files$M2, "M2", n = n)
  Y <- as_num_vector(files$Y, "Y", n = n)
  out <- list(X = X, M1 = M1, M2 = M2, Y = Y, truth = NULL,
              m1_names = colnames(M1), m2_names = colnames(M2))
  out$preview_ui <- pl2b_preview(out)
  out
}

pl2b_preview <- function(d) {
  tagList(
    bslib::layout_columns(
      col_widths = c(3, 3, 3, 3),
      bslib::value_box("Subjects", length(d$X), theme = "primary"),
      bslib::value_box("Block 1 (M1)", ncol(d$M1), theme = "secondary"),
      bslib::value_box("Block 2 (M2)", ncol(d$M2), theme = "secondary"),
      bslib::value_box("Cross paths", ncol(d$M1) * ncol(d$M2), theme = "secondary")
    ),
    tags$p(class = "small text-muted mt-2",
           "M1 feeds into M2, so the model carries three pathway families: X -> M1 -> Y, X -> M2 -> Y and X -> M1 -> M2 -> Y."),
    if (!is.null(d$truth))
      tags$p(class = "small text-success", bsicons::bs_icon("check-circle"),
             " Simulated: mediators 1-4 of each block carry signal.")
  )
}

pl2b_describe <- function(d)
  sprintf("%d subjects; block 1 has %d mediators, block 2 has %d.",
          length(d$X), ncol(d$M1), ncol(d$M2))

pl2b_run <- function(d, params) {
  k <- params$kappa %||% 5
  mu <- params$mu %||% 2
  fit <- med_fn("pathlasso2b")(d$X, d$M1, d$M2, d$Y,
                               kappa1 = k, kappa2 = k, kappa3 = k, kappa4 = k,
                               nu1 = params$nu %||% 2, nu2 = params$nu %||% 2,
                               mu1 = mu, mu2 = mu,
                               max.itr = params$max.itr %||% 3000,
                               tol = 1e-6)
  list(fit = fit, d = d, params = params)
}

pl2b_summarize <- function(res) {
  fit <- res$fit; out <- list()
  thr <- function(v) 0.05 * max(abs(v), na.rm = TRUE)
  ie1 <- as.numeric(fit$IE.M1); ie2 <- as.numeric(fit$IE.M2)
  cross <- as.matrix(fit$IE.M1M2)

  out[["Block 1 pathways (X -> M1 -> Y)"]] <- data.frame(
    Mediator = res$d$m1_names, a_path = as.numeric(fit$beta),
    b_path = as.numeric(fit$theta), IE = ie1,
    Selected = ifelse(abs(ie1) > thr(ie1), "yes", ""), check.names = FALSE)

  out[["Block 2 pathways (X -> M2 -> Y)"]] <- data.frame(
    Mediator = res$d$m2_names, a_path = as.numeric(fit$zeta),
    b_path = as.numeric(fit$pi), IE = ie2,
    Selected = ifelse(abs(ie2) > thr(ie2), "yes", ""), check.names = FALSE)

  csel <- which(abs(cross) > thr(cross), arr.ind = TRUE)
  out[["Cross-block pathways (X -> M1 -> M2 -> Y)"]] <- if (nrow(csel))
    data.frame(M1 = res$d$m1_names[csel[, 1]], M2 = res$d$m2_names[csel[, 2]],
               IE = cross[csel], check.names = FALSE)
  else data.frame(Note = "No cross-block pathway exceeded the selection threshold.")

  out[["Fit summary"]] <- data.frame(
    Quantity = c("direct effect (delta)", "block 1 selected", "block 2 selected",
                 "cross-block selected", "converged"),
    Value = c(format(round(as.numeric(fit$delta)[1], 4)),
              sum(abs(ie1) > thr(ie1)), sum(abs(ie2) > thr(ie2)),
              nrow(csel), as.character(isTRUE(fit$converge))),
    check.names = FALSE)

  out[["M1 -> M2 coefficients (Lambda)"]] <- {
    L <- as.matrix(fit$Lambda)
    lsel <- which(abs(L) > 1e-3, arr.ind = TRUE)
    if (nrow(lsel))
      data.frame(M1 = res$d$m1_names[lsel[, 1]], M2 = res$d$m2_names[lsel[, 2]],
                 Lambda = L[lsel], check.names = FALSE)
    else data.frame(Note = "All M1 -> M2 coefficients shrank to zero.")
  }

  if (!is.null(res$d$truth)) {
    tr <- res$d$truth
    t1 <- which(abs(tr$IE.M1) > 0); t2 <- which(abs(tr$IE.M2) > 0)
    s1 <- which(abs(ie1) > thr(ie1)); s2 <- which(abs(ie2) > thr(ie2))
    out[["Truth vs estimate"]] <- data.frame(
      Quantity = c("block 1 true signal", "block 1 selected",
                   "block 1 true positives",
                   "block 2 true signal", "block 2 selected",
                   "block 2 true positives"),
      Value = c(paste(t1, collapse = ", "),
                if (length(s1)) paste(s1, collapse = ", ") else "(none)",
                length(intersect(s1, t1)),
                paste(t2, collapse = ", "),
                if (length(s2)) paste(s2, collapse = ", ") else "(none)",
                length(intersect(s2, t2))),
      check.names = FALSE)
  }
  out
}

pl2b_plots <- function(res) {
  fit <- res$fit; plots <- list()
  ie1 <- as.numeric(fit$IE.M1); ie2 <- as.numeric(fit$IE.M2)

  bar_ie <- function(v, labels, truth, title) {
    df <- data.frame(Mediator = seq_along(v), IE = v)
    p <- plotly::plot_ly(df, x = ~Mediator, y = ~IE, type = "bar",
                         name = "estimate", marker = list(color = "#2c7fb8"))
    if (!is.null(truth))
      p <- plotly::add_trace(p, x = seq_along(truth), y = truth,
                             type = "scatter", mode = "markers", name = "truth",
                             marker = list(color = "grey35", size = 7,
                                           symbol = "circle-open"))
    list(plot = p |> plotly::layout(title = title,
                                    xaxis = list(title = "mediator"),
                                    yaxis = list(title = "indirect effect")),
         data = if (!is.null(truth)) cbind(df, true_IE = truth) else df)
  }
  plots[["Block 1 indirect effects"]] <-
    bar_ie(ie1, res$d$m1_names,
           if (!is.null(res$d$truth)) as.numeric(res$d$truth$IE.M1) else NULL,
           "X -> M1 -> Y pathway effects")
  plots[["Block 2 indirect effects"]] <-
    bar_ie(ie2, res$d$m2_names,
           if (!is.null(res$d$truth)) as.numeric(res$d$truth$IE.M2) else NULL,
           "X -> M2 -> Y pathway effects")

  cross <- as.matrix(fit$IE.M1M2)
  plots[["Cross-block heatmap"]] <- list(
    plot = plotly::plot_ly(z = cross, x = res$d$m2_names, y = res$d$m1_names,
                           type = "heatmap", colorscale = "RdBu", zmid = 0,
                           colorbar = list(title = "IE")) |>
      plotly::layout(title = "X -> M1 -> M2 -> Y pathway effects",
                     xaxis = list(title = "block 2 mediator"),
                     yaxis = list(title = "block 1 mediator")),
    data = data.frame(M1 = rep(res$d$m1_names, times = ncol(cross)),
                      M2 = rep(res$d$m2_names, each = nrow(cross)),
                      IE = as.numeric(cross)))
  plots
}

register_method(list(
  id = "pathlasso2b",
  name = "Multimodal",
  full_name = "Multimodal Pathway Analysis with Two Blocks of High-Dimensional Mediators",
  short = "Two mediator blocks in sequence (X -> M1 -> M2 -> Y), so the model separates block-1 pathways, block-2 pathways, and pathways that run through both.",
  status = "ready",
  tags = c("two mediator blocks", "multimodal", "pathway selection", "ADMM"),
  paper = list(
    citation = "Zhao, Y., Li, L., & Caffo, B. S. (2021). Multimodal neuroimaging data integration and pathway analysis. Biometrics, 77(3), 879-889.",
    url = "https://doi.org/10.1111/biom.13351"),
  explain = file.path(APP_DIR, "methods", "pathlasso2b", "explain.md"),
  example_note = paste("200 subjects, 20 mediators in block 1 and 30 in block 2;",
                       "the first four of each block carry signal."),
  data_inputs = list(
    list(id = "X", label = "Exposure (n values)",
         help = "A single exposure, one value per subject."),
    list(id = "M1", label = "Mediator block 1 (n x p1)",
         help = "The upstream mediator block, one row per subject."),
    list(id = "M2", label = "Mediator block 2 (n x p2)",
         help = "The downstream mediator block, which M1 may feed into."),
    list(id = "Y", label = "Outcome (n values)",
         help = "One numeric value per subject, in the same order.")
  ),
  params = list(
    list(id = "kappa", label = "Pathway penalty (kappa)", type = "numeric",
         default = 5, min = 0, max = 50, step = 0.5,
         help = "Applied to all four penalty terms. Larger is sparser."),
    list(id = "mu", label = "Extra l1 on the paths (mu)", type = "numeric",
         default = 2, min = 0, max = 20, step = 0.25,
         help = "Adds separate shrinkage on each path; raise it for a sparser solution."),
    list(id = "nu", label = "Ridge weight (nu)", type = "numeric",
         default = 2, min = 0, max = 10, step = 0.5),
    list(id = "max.itr", label = "Max ADMM iterations", type = "integer",
         default = 3000, min = 500, max = 20000,
         help = "This method needs iterations in the thousands; at a few hundred it reports non-convergence.")
  ),
  example = pl2b_app_example,
  export_example = function(d) list(X = data.frame(X = d$X),
                                    M1 = matrix_to_df(d$M1),
                                    M2 = matrix_to_df(d$M2),
                                    Y = data.frame(Y = d$Y)),
  parse = pl2b_app_parse,
  describe_data = pl2b_describe,
  run = pl2b_run,
  summarize = pl2b_summarize,
  plots = pl2b_plots
))
