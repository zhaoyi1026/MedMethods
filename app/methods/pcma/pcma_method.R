# =============================================================================
# PCMA -- principal component mediation analysis
# -----------------------------------------------------------------------------
# Orthogonal projections of the exposures (Phi) and mediators (Psi) define
# component scores carrying parallel mediation mechanisms:
#   M Psi = (X Phi) alpha + e,   Y = (X Phi) gamma + (M Psi) beta + u
#   indirect effect per component  IE = alpha * beta
# =============================================================================

pcma_app_example <- function() {
  d <- med_fn("pcma_example")()             # n = 400, p = 5, q = 10
  # NB: keep Y as an n x 1 MATRIX -- pcma's bootstrap resamples with Y[idx, ],
  # which fails on a plain vector ("incorrect number of dimensions").
  out <- list(X = d$X, M = d$M, Y = as.matrix(d$Y), truth = d$truth,
              x_names = colnames(d$X) %||% paste0("X", seq_len(ncol(d$X))),
              m_names = colnames(d$M) %||% paste0("M", seq_len(ncol(d$M))))
  out$preview_ui <- pcma_preview(out)
  out
}

pcma_app_parse <- function(files, opts) {
  X <- as_num_matrix(files$X, "X")
  n <- nrow(X)
  M <- as_num_matrix(files$M, "M", n = n)
  Y <- matrix(as_num_vector(files$Y, "Y", n = n), ncol = 1)   # matrix, see above
  out <- list(X = X, M = M, Y = Y, truth = NULL,
              x_names = colnames(X), m_names = colnames(M))
  out$preview_ui <- pcma_preview(out)
  out
}

pcma_preview <- function(d) {
  xmy_preview(d, extra = list(bslib::value_box("Outcome", "1 per subject",
                                               theme = "secondary")),
              note = sprintf("%d exposures projected to components; %d mediators projected to components.",
                             ncol(d$X), ncol(d$M)))
}

pcma_describe <- function(d)
  sprintf("%d subjects; %d exposures; %d mediators.",
          nrow(d$X), ncol(d$X), ncol(d$M))

pcma_run <- function(d, params) {
  nD <- params$nD %||% 1
  sims <- params$sims %||% 0
  fit <- med_fn("pcma")(d$X, d$M, d$Y, stop.crt = "nD", nD = nD,
                        boot = sims > 0, sims = max(sims, 1),
                        boot.ci.type = "perc",
                        conf.level = params$conf.level %||% 0.95,
                        ninitial = params$ninitial %||% 3,
                        seed = 100, verbose = FALSE)
  list(fit = fit, d = d, nD = nD, sims = sims, params = params)
}

pcma_summarize <- function(res) {
  fit <- res$fit
  Phi <- as.matrix(fit$Phi); Psi <- as.matrix(fit$Psi)
  nD <- ncol(Phi); dirs <- paste0("D", seq_len(nD)); out <- list()

  if (res$sims > 0 && length(fit$coef.inference)) {
    eff <- do.call(rbind, lapply(seq_along(fit$coef.inference), function(j) {
      ci <- fit$coef.inference[[j]]
      if (is.null(ci) || !length(ci)) return(NULL)
      data.frame(Component = dirs[j], Effect = rownames(ci),
                 as.data.frame(ci), check.names = FALSE, row.names = NULL)
    }))
    if (!is.null(eff)) out[["Effects (bootstrap)"]] <- eff
  } else {
    # no bootstrap requested: report the point estimates at the fitted projections
    eff <- do.call(rbind, lapply(seq_len(nD), function(j) {
      cf <- med_fn("pcma_coef")(res$d$X, res$d$M, res$d$Y,
                                phi = Phi[, j, drop = FALSE],
                                psi = Psi[, j, drop = FALSE])
      data.frame(Component = dirs[j],
                 alpha = as.numeric(cf$alpha)[1], beta = as.numeric(cf$beta)[1],
                 gamma = as.numeric(cf$gamma)[1], IE = as.numeric(cf$IE)[1],
                 check.names = FALSE)
    }))
    out[["Effects (point estimates)"]] <- eff
  }

  ph <- as.data.frame(Phi); names(ph) <- dirs
  out[["Phi exposure loadings"]] <- cbind(Exposure = res$d$x_names, ph)
  ps <- as.data.frame(Psi); names(ps) <- dirs
  out[["Psi mediator loadings"]] <- cbind(Mediator = res$d$m_names, ps)

  if (nD > 1) {
    op <- as.data.frame(round(t(Phi) %*% Phi, 6)); names(op) <- dirs
    out[["Orthogonality Phi'Phi"]] <- cbind(" " = dirs, op)
    os <- as.data.frame(round(t(Psi) %*% Psi, 6)); names(os) <- dirs
    out[["Orthogonality Psi'Psi"]] <- cbind(" " = dirs, os)
  }

  if (!is.null(res$d$truth)) {
    tr <- res$d$truth
    k <- min(nD, ncol(tr$Phi))
    out[["Truth vs estimate"]] <- data.frame(
      Component = dirs[seq_len(k)],
      Phi_cosine = vapply(seq_len(k), function(j) abs_cos(Phi[, j], tr$Phi[, j]),
                          numeric(1)),
      Psi_cosine = vapply(seq_len(k), function(j) abs_cos(Psi[, j], tr$Psi[, j]),
                          numeric(1)),
      True_alpha = diag(tr$alpha)[seq_len(k)],
      True_beta = tr$beta[seq_len(k)],
      check.names = FALSE)
  }
  out
}

pcma_plots <- function(res) {
  fit <- res$fit
  Phi <- as.matrix(fit$Phi); Psi <- as.matrix(fit$Psi)
  nD <- ncol(Phi); dirs <- paste0("D", seq_len(nD)); plots <- list()

  loadings_plot <- function(W, labels, title, axis) {
    long <- do.call(rbind, lapply(seq_len(ncol(W)), function(j)
      data.frame(Term = labels, Component = dirs[j], Loading = W[, j])))
    long$Term <- factor(long$Term, levels = labels)
    list(plot = plotly::plot_ly(long, x = ~Term, y = ~Loading,
                                color = ~Component, type = "bar") |>
           plotly::layout(title = title, barmode = "group",
                          xaxis = list(title = axis)),
         data = long)
  }
  plots[["Phi exposure loadings"]] <-
    loadings_plot(Phi, res$d$x_names, "Exposure projection loadings (Phi)", "exposure")
  plots[["Psi mediator loadings"]] <-
    loadings_plot(Psi, res$d$m_names, "Mediator projection loadings (Psi)", "mediator")

  # component scores: mediator score against exposure score, the a-path
  xs <- as.numeric(res$d$X %*% Phi[, 1])
  ms <- as.numeric(res$d$M %*% Psi[, 1])
  sdf <- data.frame(exposure_score = xs, mediator_score = ms,
                    Y = as.numeric(res$d$Y))
  plots[["Component 1 a-path"]] <- list(
    plot = plotly::plot_ly(sdf, x = ~exposure_score, y = ~mediator_score,
                           type = "scatter", mode = "markers",
                           marker = list(size = 6, opacity = 0.6,
                                         color = "#2c7fb8")) |>
      plotly::layout(title = "Component 1: mediator score vs exposure score",
                     xaxis = list(title = "X Phi_1"),
                     yaxis = list(title = "M Psi_1")),
    data = sdf)
  plots[["Component 1 b-path"]] <- list(
    plot = plotly::plot_ly(sdf, x = ~mediator_score, y = ~Y,
                           type = "scatter", mode = "markers",
                           marker = list(size = 6, opacity = 0.6,
                                         color = "#c0392b")) |>
      plotly::layout(title = "Component 1: outcome vs mediator score",
                     xaxis = list(title = "M Psi_1"),
                     yaxis = list(title = "outcome Y")),
    data = sdf)

  if (res$sims > 0 && length(fit$coef.inference)) {
    ci <- fit$coef.inference[[1]]
    if (!is.null(ci) && length(ci)) {
      keep <- intersect(c("alpha", "beta", "gamma", "IE"), rownames(ci))
      edf <- data.frame(Effect = factor(keep, levels = keep),
                        Estimate = ci[keep, "Estimate"],
                        LB = ci[keep, "LB"], UB = ci[keep, "UB"])
      plots[["Effects with CI"]] <- list(
        plot = plotly::plot_ly(edf, x = ~Effect, y = ~Estimate, type = "bar",
                               marker = list(color = "#2c7fb8"),
                               error_y = list(type = "data", symmetric = FALSE,
                                              array = edf$UB - edf$Estimate,
                                              arrayminus = edf$Estimate - edf$LB)) |>
          plotly::layout(title = "Component 1 effects with bootstrap CI",
                         yaxis = list(title = "estimate")),
        data = edf)
    }
  }
  plots
}

register_method(list(
  id = "pcma",
  name = "PCMA",
  full_name = "Principal Component Mediation Analysis",
  short = "Finds orthogonal projections of the exposures and of the mediators whose component scores carry parallel mediation mechanisms, giving one indirect effect per component.",
  status = "ready",
  tags = c("multiple exposures", "multiple mediators", "projections", "bootstrap"),
  paper = list(
    citation = "Zhao, Y. (2024). Mediation analysis with multiple exposures and multiple mediators. Statistics in Medicine, 43(25), 4887-4898.",
    url = "https://doi.org/10.1002/sim.10215"),
  explain = file.path(APP_DIR, "methods", "pcma", "explain.md"),
  example_note = paste("400 subjects, 5 exposures and 10 mediators, with two",
                       "mediating components (true alpha = 2, 2 and beta = 2, 1)."),
  data_inputs = list(
    list(id = "X", label = "Exposures (n x p)",
         help = "An n x p numeric matrix, one row per subject. A CSV whose first column is an id is also accepted (the id column is dropped)."),
    list(id = "M", label = "Mediators (n x q)",
         help = "An n x q numeric matrix in the same subject order as the exposures."),
    list(id = "Y", label = "Outcome (n values)",
         help = "One numeric value per subject, in the same order.")
  ),
  params = list(
    list(id = "nD", label = "Number of components (nD)", type = "integer",
         default = 1, min = 1, max = 5,
         help = "Components are extracted one at a time, deflating the data after each."),
    list(id = "ninitial", label = "Random initialisations", type = "integer",
         default = 3, min = 1, max = 20,
         help = "Restarts of the projection search; the best log-likelihood wins."),
    list(id = "sims", label = "Bootstrap replicates", type = "integer",
         default = 200, min = 0, max = 2000,
         help = "Bootstrap inference for alpha / beta / gamma / IE; 0 reports point estimates only."),
    list(id = "conf.level", label = "Confidence level", type = "numeric",
         default = 0.95, min = 0.5, max = 0.999, step = 0.01)
  ),
  example = pcma_app_example,
  export_example = function(d) list(X = matrix_to_df(d$X), M = matrix_to_df(d$M),
                                    Y = data.frame(Y = as.numeric(d$Y))),
  parse = pcma_app_parse,
  describe_data = pcma_describe,
  run = pcma_run,
  summarize = pcma_summarize,
  plots = pcma_plots
))
