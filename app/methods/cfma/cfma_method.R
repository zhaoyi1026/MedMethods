# =============================================================================
# cfma -- causal functional mediation analysis
# -----------------------------------------------------------------------------
# Treatment, mediator and outcome are all functions of time. The concurrent
# model links them at each time point:
#   M(t) = Z(t) alpha(t) + e1(t)
#   Y(t) = Z(t) gamma(t) + M(t) beta(t) + e2(t)
#   indirect effect curve  IE(t) = alpha(t) beta(t)
# The historical model instead integrates the effect over the recent past.
# =============================================================================

cfma_app_example <- function() {
  d <- med_fn("cfma_example")()             # N = 200 subjects, 150 time points
  out <- list(Z = d$Z, M = d$M, Y = d$Y, timeinv = d$timeinv, truth = d$truth)
  out$preview_ui <- cfma_preview(out)
  out
}

cfma_app_parse <- function(files, opts) {
  Z <- as_num_matrix(files$Z, "Z", drop_id = FALSE)
  M <- as_num_matrix(files$M, "M", n = nrow(Z), drop_id = FALSE)
  Y <- as_num_matrix(files$Y, "Y", n = nrow(Z), drop_id = FALSE)
  if (ncol(M) != ncol(Z) || ncol(Y) != ncol(Z))
    stop(sprintf("Z, M and Y must share the same number of time points (got %d, %d, %d).",
                 ncol(Z), ncol(M), ncol(Y)))
  out <- list(Z = Z, M = M, Y = Y, timeinv = c(0, 1), truth = NULL)
  out$preview_ui <- cfma_preview(out)
  out
}

cfma_preview <- function(d) {
  tagList(
    bslib::layout_columns(
      col_widths = c(4, 4, 4),
      bslib::value_box("Subjects", nrow(d$Z), theme = "primary"),
      bslib::value_box("Time points", ncol(d$Z), theme = "secondary"),
      bslib::value_box("Time range",
                       paste(d$timeinv, collapse = " to "), theme = "secondary")
    ),
    tags$p(class = "small text-muted mt-2",
           "Each of the treatment, mediator and outcome is one curve per subject; the coefficients are curves too."),
    if (!is.null(d$truth))
      tags$p(class = "small text-success", bsicons::bs_icon("check-circle"),
             " Simulated: the true alpha(t), beta(t) and gamma(t) curves are known, so the results overlay them on the estimates.")
  )
}

cfma_describe <- function(d)
  sprintf("%d subjects, each with %d time points.", nrow(d$Z), ncol(d$Z))

cfma_run <- function(d, params) {
  model <- params$model %||% "concurrent"
  if (model == "concurrent") {
    fit <- med_fn("cfma_concurrent")(d$Z, d$M, d$Y,
                                     intercept = isTRUE(params$intercept),
                                     nbasis = params$nbasis %||% 3,
                                     timeinv = d$timeinv,
                                     lambda.m = params$lambda %||% 0.01,
                                     lambda.y = params$lambda %||% 0.01)
  } else {
    fit <- med_fn("cfma_historical")(d$Z, d$M, d$Y,
                                     intercept = isTRUE(params$intercept),
                                     nbasis1 = params$nbasis %||% 3,
                                     nbasis2 = params$nbasis %||% 3,
                                     timeinv = d$timeinv,
                                     lambda1.m = params$lambda %||% 0.01,
                                     lambda2.m = params$lambda %||% 0.01,
                                     lambda1.y = params$lambda %||% 0.01,
                                     lambda2.y = params$lambda %||% 0.01)
  }
  list(fit = fit, d = d, model = model, params = params)
}

# pull the curves out of the fit, guarding against the historical model's
# different shape (there the "curves" are surfaces, so only IE/DE are 1-D)
.cfma_curves <- function(fit) {
  g <- function(x) if (is.null(x)) NULL else as.numeric(x)
  alpha <- if (!is.null(fit$M$curve) && is.matrix(fit$M$curve)) g(fit$M$curve[1, ]) else NULL
  gamma <- if (!is.null(fit$Y$curve) && is.matrix(fit$Y$curve) && nrow(fit$Y$curve) >= 1)
    g(fit$Y$curve[1, ]) else NULL
  beta  <- if (!is.null(fit$Y$curve) && is.matrix(fit$Y$curve) && nrow(fit$Y$curve) >= 2)
    g(fit$Y$curve[2, ]) else NULL
  list(alpha = alpha, beta = beta, gamma = gamma,
       IE = g(fit$IE$curve), DE = g(fit$DE$curve))
}

cfma_summarize <- function(res) {
  cur <- .cfma_curves(res$fit); out <- list()
  L <- max(vapply(cur, function(v) length(v) %||% 0L, integer(1)))
  tg <- seq(res$d$timeinv[1], res$d$timeinv[2], length.out = max(L, 1L))
  pad <- function(v) if (is.null(v) || !length(v)) rep(NA_real_, L) else
    c(v, rep(NA_real_, L - length(v)))[seq_len(L)]

  out[["Estimated curves"]] <- data.frame(
    time = tg, alpha = pad(cur$alpha), beta = pad(cur$beta),
    gamma = pad(cur$gamma), IE = pad(cur$IE), DE = pad(cur$DE),
    check.names = FALSE)

  summ <- function(v) if (is.null(v) || !length(v)) c(NA, NA, NA) else
    c(mean(v), min(v), max(v))
  out[["Curve summaries"]] <- data.frame(
    Curve = c("alpha(t)  treatment -> mediator", "beta(t)   mediator -> outcome",
              "gamma(t)  direct", "IE(t)  indirect", "DE(t)  direct"),
    Mean = c(summ(cur$alpha)[1], summ(cur$beta)[1], summ(cur$gamma)[1],
             summ(cur$IE)[1], summ(cur$DE)[1]),
    Min = c(summ(cur$alpha)[2], summ(cur$beta)[2], summ(cur$gamma)[2],
            summ(cur$IE)[2], summ(cur$DE)[2]),
    Max = c(summ(cur$alpha)[3], summ(cur$beta)[3], summ(cur$gamma)[3],
            summ(cur$IE)[3], summ(cur$DE)[3]),
    check.names = FALSE)

  if (!is.null(res$d$truth) && res$model == "concurrent") {
    tr <- res$d$truth
    cmp <- function(est, tru) {
      if (is.null(est) || !length(est) || length(est) != length(tru))
        return(c(NA_real_, NA_real_))
      c(stats::cor(est, tru), sqrt(mean((est - tru)^2)))
    }
    a <- cmp(cur$alpha, tr$alpha); b <- cmp(cur$beta, tr$beta)
    g <- cmp(cur$gamma, tr$gamma); i <- cmp(cur$IE, tr$IE)
    out[["Truth vs estimate"]] <- data.frame(
      Curve = c("alpha(t)", "beta(t)", "gamma(t)", "IE(t)"),
      Correlation = c(a[1], b[1], g[1], i[1]),
      RMSE = c(a[2], b[2], g[2], i[2]),
      check.names = FALSE)
  }
  out
}

cfma_plots <- function(res) {
  cur <- .cfma_curves(res$fit); plots <- list()
  tr <- res$d$truth
  L <- length(cur$IE %||% cur$alpha %||% numeric(0))
  if (!L) return(plots)
  tg <- seq(res$d$timeinv[1], res$d$timeinv[2], length.out = L)

  one <- function(est, tru, title, ylab) {
    if (is.null(est) || !length(est)) return(NULL)
    n <- min(length(est), L)
    df <- data.frame(time = tg[seq_len(n)], estimate = est[seq_len(n)])
    p <- plotly::plot_ly(df, x = ~time, y = ~estimate, type = "scatter",
                         mode = "lines", name = "estimate",
                         line = list(color = "#c0392b", width = 2))
    if (!is.null(tru) && length(tru) == n) {
      df$truth <- tru
      p <- plotly::add_trace(p, y = ~truth, name = "truth",
                             line = list(color = "grey40", width = 5,
                                         dash = "dot"))
    }
    list(plot = p |> plotly::layout(title = title,
                                    xaxis = list(title = "time"),
                                    yaxis = list(title = ylab)),
         data = df)
  }
  add <- function(nm, obj) if (!is.null(obj)) plots[[nm]] <<- obj

  add("alpha(t): treatment -> mediator",
      one(cur$alpha, tr$alpha, "Treatment -> mediator curve", "alpha(t)"))
  add("beta(t): mediator -> outcome",
      one(cur$beta, tr$beta, "Mediator -> outcome curve", "beta(t)"))
  add("gamma(t): direct",
      one(cur$gamma, tr$gamma, "Direct-effect curve", "gamma(t)"))
  add("IE(t): indirect effect",
      one(cur$IE, tr$IE, "Indirect-effect curve", "IE(t)"))
  add("DE(t): direct effect",
      one(cur$DE, if (!is.null(tr)) tr$DE else NULL,
          "Direct-effect curve (DE)", "DE(t)"))

  # a few subject trajectories, so the input shape is visible
  ns <- min(nrow(res$d$Z), 6L)
  tdf <- do.call(rbind, lapply(seq_len(ns), function(i) data.frame(
    time = seq(res$d$timeinv[1], res$d$timeinv[2], length.out = ncol(res$d$Z)),
    value = as.numeric(res$d$M[i, ]), subject = paste0("S", i))))
  plots[["Mediator trajectories"]] <- list(
    plot = plotly::plot_ly(tdf, x = ~time, y = ~value, color = ~subject,
                           type = "scatter", mode = "lines") |>
      plotly::layout(title = sprintf("Mediator curves for the first %d subjects", ns),
                     xaxis = list(title = "time"),
                     yaxis = list(title = "M(t)")),
    data = tdf)
  plots
}

register_method(list(
  id = "cfma",
  name = "cfma",
  full_name = "Causal Functional Mediation Analysis",
  short = "Treatment, mediator and outcome are all curves over time, so the mediation coefficients are curves too. Offers a concurrent model (effects act at the same instant) and a historical model (effects integrate over the recent past).",
  status = "ready",
  tags = c("functional data", "time-varying effects", "fMRI"),
  paper = list(
    citation = "Zhao, Y., Luo, X., Sobel, M. E., Lindquist, M. A., & Caffo, B. S. (2025). Causal functional mediation analysis with an application to functional magnetic resonance imaging data. Biostatistics, 26(1), kxaf019.",
    url = "https://doi.org/10.1093/biostatistics/kxaf019"),
  explain = file.path(APP_DIR, "methods", "cfma", "explain.md"),
  example_note = paste("200 subjects observed at 150 time points, with known",
                       "alpha(t), beta(t) and gamma(t) curves overlaid on the",
                       "estimates in the results."),
  data_inputs = list(
    list(id = "Z", label = "Treatment curves (N x T)",
         help = "One row per subject, one column per time point. All three files must share the same dimensions."),
    list(id = "M", label = "Mediator curves (N x T)",
         help = "One row per subject, one column per time point, same order as the treatment."),
    list(id = "Y", label = "Outcome curves (N x T)",
         help = "One row per subject, one column per time point, same order as the treatment.")
  ),
  params = list(
    list(id = "model", label = "Model", type = "select",
         choices = c("Concurrent" = "concurrent",
                     "Historical influence" = "historical"),
         default = "concurrent",
         help = "Concurrent links the curves instant by instant; historical integrates the effect over the recent past (slower)."),
    list(id = "nbasis", label = "Fourier basis functions", type = "integer",
         default = 3, min = 1, max = 15,
         help = "More basis functions allow wigglier coefficient curves."),
    list(id = "lambda", label = "Smoothing penalty (lambda)", type = "numeric",
         default = 0.01, min = 0, max = 10, step = 0.005,
         help = "Penalises curvature in the coefficient curves; applied to both the mediator and outcome models."),
    list(id = "intercept", label = "Include an intercept curve",
         type = "checkbox", default = FALSE)
  ),
  example = cfma_app_example,
  export_example = function(d) list(Z = matrix_to_df(d$Z), M = matrix_to_df(d$M),
                                    Y = matrix_to_df(d$Y)),
  parse = cfma_app_parse,
  describe_data = cfma_describe,
  run = cfma_run,
  summarize = cfma_summarize,
  plots = cfma_plots
))
