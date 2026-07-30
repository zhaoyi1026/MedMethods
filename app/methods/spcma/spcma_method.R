# =============================================================================
# spcma -- sparse principal component based high-dimensional mediation analysis
# -----------------------------------------------------------------------------
# Mediation through the leading (sparse) principal components of a
# high-dimensional mediator:
#   Mtilde_j = alpha_j X + xi_j          (component score j)
#   Y        = gamma X + sum_j beta_j Mtilde_j + eta
# The loadings are regularised by a fused lasso, so they come out piecewise
# constant. `PC.run = TRUE` also fits the plain-PCA version for comparison.
# =============================================================================

spcma_app_example <- function() {
  d <- med_fn("spcma_example")()            # n = 200, p = 50, 3 blocks of 10
  out <- list(X = as.matrix(d$X), M = d$M, Y = as.numeric(d$Y), truth = d$truth,
              m_names = colnames(d$M) %||% paste0("M", seq_len(ncol(d$M))))
  out$preview_ui <- spcma_preview(out)
  out
}

spcma_app_parse <- function(files, opts) {
  X <- as_num_matrix(files$X, "X")
  if (ncol(X) != 1L)
    stop("spcma takes a single exposure: 'X' has ", ncol(X), " columns.")
  n <- nrow(X)
  M <- as_num_matrix(files$M, "M", n = n)
  Y <- as_num_vector(files$Y, "Y", n = n)
  out <- list(X = X, M = M, Y = Y, truth = NULL, m_names = colnames(M))
  out$preview_ui <- spcma_preview(out)
  out
}

spcma_preview <- function(d) {
  tagList(
    bslib::layout_columns(
      col_widths = c(4, 4, 4),
      bslib::value_box("Subjects", nrow(d$M), theme = "primary"),
      bslib::value_box("Mediators (p)", ncol(d$M), theme = "secondary"),
      bslib::value_box("Exposure", "1 (binary or continuous)", theme = "secondary")
    ),
    tags$p(class = "small text-muted mt-2",
           "The mediators are reduced to a few sparse principal components; mediation is assessed component by component."),
    if (!is.null(d$truth))
      tags$p(class = "small text-success", bsicons::bs_icon("check-circle"),
             sprintf(" Simulated: the leading components are blocks of %d consecutive mediators, with true IE = 4, -1.5, 0.",
                     d$truth$n.block %||% 10))
  )
}

spcma_describe <- function(d)
  sprintf("%d subjects; %d mediators; single exposure.", nrow(d$M), ncol(d$M))

spcma_run <- function(d, params) {
  sims <- params$sims %||% 0
  fit <- med_fn("spcma")(d$X, d$M, d$Y,
                         adaptive = isTRUE(params$adaptive),
                         var.per = params$var.per %||% 0.8,
                         gamma = params$gamma %||% 0,
                         per.jump = params$per.jump %||% 0.7,
                         boot = sims > 0, sims = max(sims, 1),
                         boot.ci.type = "bca",
                         conf.level = params$conf.level %||% 0.95,
                         p.adj.method = params$p.adj.method %||% "BH",
                         PC.run = isTRUE(params$PC.run))
  list(fit = fit, d = d, sims = sims, params = params)
}

.spcma_eff_table <- function(part, label) {
  if (is.null(part)) return(NULL)
  ie <- as.data.frame(part$IE)
  comp <- paste0("PC", seq_len(nrow(ie)))
  cbind(Model = label, Component = comp, ie)
}

spcma_summarize <- function(res) {
  fit <- res$fit; out <- list()
  sp <- fit$SPCA; pc <- fit$PCA

  ie <- rbind(.spcma_eff_table(sp, "Sparse PC"), .spcma_eff_table(pc, "Plain PC"))
  if (!is.null(ie)) out[["Indirect effect by component"]] <- ie

  parts <- list("Sparse PC" = sp, "Plain PC" = pc)
  for (nm in c("alpha", "beta", "gamma")) {
    rows <- NULL
    for (lab in names(parts)) {
      pt <- parts[[lab]]
      if (is.null(pt) || is.null(pt[[nm]])) next
      m <- as.data.frame(pt[[nm]])
      rows <- rbind(rows, cbind(Model = lab,
                                Component = paste0("PC", seq_len(nrow(m))), m))
    }
    if (!is.null(rows)) out[[paste0(nm, " estimates")]] <- rows
  }

  if (!is.null(sp$W)) {
    W <- as.matrix(sp$W)
    df <- as.data.frame(W); names(df) <- paste0("PC", seq_len(ncol(W)))
    out[["Sparse loadings (W)"]] <- cbind(Mediator = res$d$m_names, df)
  }
  if (!is.null(sp$var.per))
    out[["Variance explained"]] <- data.frame(
      Quantity = "cumulative proportion of mediator variance retained",
      Value = as.numeric(sp$var.per)[1], check.names = FALSE)

  if (!is.null(res$d$truth) && !is.null(sp$W)) {
    tr <- res$d$truth; W <- as.matrix(sp$W)
    k <- min(ncol(W), 3L)
    out[["Truth vs estimate"]] <- data.frame(
      Component = paste0("PC", seq_len(k)),
      Loading_cosine = vapply(seq_len(k),
                              function(j) abs_cos(W[, j], tr$Phi[, j]), numeric(1)),
      True_IE = tr$IE[seq_len(k)],
      Estimated_IE = as.numeric(sp$IE[seq_len(k), "Estimate"]),
      check.names = FALSE)
  }
  out
}

spcma_plots <- function(res) {
  fit <- res$fit; sp <- fit$SPCA; plots <- list()

  if (!is.null(sp$W)) {
    W <- as.matrix(sp$W); k <- ncol(W)
    long <- do.call(rbind, lapply(seq_len(k), function(j)
      data.frame(Mediator = res$d$m_names, Component = paste0("PC", j),
                 Loading = W[, j])))
    long$Mediator <- factor(long$Mediator, levels = res$d$m_names)
    plots[["Sparse loadings"]] <- list(
      plot = plotly::plot_ly(long, x = ~Mediator, y = ~Loading,
                             color = ~Component, type = "bar") |>
        plotly::layout(title = "Sparse principal component loadings (fused lasso)",
                       barmode = "group",
                       xaxis = list(title = "mediator")),
      data = long)
  }

  if (!is.null(sp$IE)) {
    ie <- as.data.frame(sp$IE)
    ie$Component <- factor(paste0("PC", seq_len(nrow(ie))),
                           levels = paste0("PC", seq_len(nrow(ie))))
    has_ci <- all(c("LB", "UB") %in% names(ie))
    p <- plotly::plot_ly(ie, x = ~Component, y = ~Estimate, type = "bar",
                         marker = list(color = "#2c7fb8"),
                         error_y = if (has_ci)
                           list(type = "data", symmetric = FALSE,
                                array = ie$UB - ie$Estimate,
                                arrayminus = ie$Estimate - ie$LB) else NULL)
    plots[["Indirect effect by component"]] <- list(
      plot = p |> plotly::layout(title = "Indirect effect per sparse component",
                                 yaxis = list(title = "IE")),
      data = ie)
  }

  # component score vs outcome for the leading component
  if (!is.null(sp$W)) {
    sc <- as.numeric(res$d$M %*% as.matrix(sp$W)[, 1])
    sdf <- data.frame(component_score = sc, Y = res$d$Y,
                      Exposure = factor(res$d$X[, 1]))
    plots[["Component 1 -> outcome"]] <- list(
      plot = plotly::plot_ly(sdf, x = ~component_score, y = ~Y, color = ~Exposure,
                             type = "scatter", mode = "markers",
                             marker = list(size = 6, opacity = 0.6)) |>
        plotly::layout(title = "Outcome vs leading component score, by exposure",
                       xaxis = list(title = "M W_1"),
                       yaxis = list(title = "outcome Y")),
      data = sdf)
  }
  plots
}

register_method(list(
  id = "spcma",
  name = "spcma",
  full_name = "Sparse Principal Component Based High-Dimensional Mediation Analysis",
  short = "Reduces a high-dimensional mediator to a few principal components whose loadings are regularised by a fused lasso (so they are piecewise constant), then assesses mediation component by component.",
  status = "ready",
  tags = c("high-dimensional mediators", "sparse PCA", "fused lasso", "bootstrap"),
  paper = list(
    citation = "Zhao, Y., Lindquist, M. A., & Caffo, B. S. (2020). Sparse principal component based high-dimensional mediation analysis. Computational Statistics & Data Analysis, 142, 106835.",
    url = "https://doi.org/10.1016/j.csda.2019.106835"),
  explain = file.path(APP_DIR, "methods", "spcma", "explain.md"),
  example_note = paste("200 subjects and 50 mediators whose leading components",
                       "are blocks of 10 consecutive variables; true component",
                       "indirect effects 4, -1.5 and 0."),
  data_inputs = list(
    list(id = "X", label = "Exposure (n values)",
         help = "A single exposure, one value per subject (binary or continuous)."),
    list(id = "M", label = "Mediators (n x p)",
         help = "An n x p numeric matrix of mediators, same subject order as the exposure."),
    list(id = "Y", label = "Outcome (n values)",
         help = "One numeric value per subject, in the same order.")
  ),
  params = list(
    list(id = "adaptive", label = "Choose the number of components adaptively",
         type = "checkbox", default = TRUE,
         help = "Keep components until the retained variance proportion below is reached."),
    list(id = "var.per", label = "Variance proportion to retain", type = "numeric",
         default = 0.8, min = 0.1, max = 0.999, step = 0.05),
    list(id = "gamma", label = "Fused-lasso sparsity (gamma)", type = "numeric",
         default = 0, min = 0, max = 5, step = 0.1,
         help = "0 penalises only differences between neighbouring loadings; larger values also shrink the loadings themselves toward zero."),
    list(id = "per.jump", label = "Loading-path jump proportion", type = "numeric",
         default = 0.7, min = 0.1, max = 1, step = 0.05,
         help = "Where along the fused-lasso path the loadings are read off."),
    list(id = "PC.run", label = "Also fit the plain-PCA version",
         type = "checkbox", default = TRUE,
         help = "Provides an unregularised comparison in the results tables."),
    list(id = "sims", label = "Bootstrap replicates", type = "integer",
         default = 200, min = 0, max = 2000,
         help = "BCa bootstrap for the component effects; 0 skips inference."),
    list(id = "p.adj.method", label = "Multiplicity adjustment", type = "select",
         choices = c("BH" = "BH", "Bonferroni" = "bonferroni", "BY" = "BY"),
         default = "BH"),
    list(id = "conf.level", label = "Confidence level", type = "numeric",
         default = 0.95, min = 0.5, max = 0.999, step = 0.01)
  ),
  example = spcma_app_example,
  export_example = function(d) list(X = data.frame(X = d$X[, 1]),
                                    M = matrix_to_df(d$M),
                                    Y = data.frame(Y = d$Y)),
  parse = spcma_app_parse,
  describe_data = spcma_describe,
  run = spcma_run,
  summarize = spcma_summarize,
  plots = spcma_plots
))
