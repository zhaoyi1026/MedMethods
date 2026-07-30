# =============================================================================
# GMed -- mediation analysis with a graph (covariance-matrix) mediator
# -----------------------------------------------------------------------------
# The mediator is a subject-level covariance matrix; a projection theta
# summarises it as a scalar log-variance:
#   mediator:  log(theta' Sigma_i theta) = alpha0 + x_i' alpha + b_i
#   outcome:   Y_i = gamma0 + x_i' gamma + beta * log(theta' Sigma_i theta) + e_i
#   indirect effect  IE = alpha[exposure] * beta
# =============================================================================

gmed_app_example <- function() {
  d <- med_fn("gmed_example")()          # n = 100, p = 10, Ti = 500
  ids <- paste0("S", seq_along(d$M))
  names(d$M) <- ids
  out <- list(M = d$M, X = as.matrix(d$X), Y = as.numeric(d$Y), ids = ids,
              m_names = paste0("V", seq_len(ncol(d$M[[1]]))),
              x_names = colnames(d$X) %||% "X",
              truth = d$truth)
  out$preview_ui <- gmed_preview(out)
  out
}

gmed_app_parse <- function(files, opts) {
  mp <- as_matrix_list(files$M, "M")
  n <- length(mp$M)
  X <- as_num_matrix(files$X, "X", n = n)      # NO intercept: added internally
  Y <- as_num_vector(files$Y, "Y", n = n)
  out <- list(M = mp$M, X = X, Y = Y, ids = mp$ids,
              m_names = mp$var_names, x_names = colnames(X), truth = NULL)
  out$preview_ui <- gmed_preview(out)
  out
}

gmed_preview <- function(d) {
  Tv <- vapply(d$M, nrow, integer(1))
  rng <- function(x) if (length(unique(x)) == 1) as.character(x[1])
                     else paste0(min(x), "-", max(x))
  tagList(
    bslib::layout_columns(
      col_widths = c(3, 3, 3, 3),
      bslib::value_box("Subjects", length(d$M), theme = "primary"),
      bslib::value_box("Mediator dim (p)", ncol(d$M[[1]]), theme = "secondary"),
      bslib::value_box("Rows / subject", rng(Tv), theme = "secondary"),
      bslib::value_box("Exposure", d$x_names[1], theme = "secondary")
    ),
    tags$p(class = "small text-muted mt-2",
           sprintf("Each subject's mediator is the covariance of their %s matrix.%s",
                   paste(rng(Tv), "x", ncol(d$M[[1]])),
                   if (length(d$x_names) > 1)
                     paste0(" Adjusted covariates: ",
                            paste(d$x_names[-1], collapse = ", "), ".") else "")),
    if (!is.null(d$truth))
      tags$p(class = "small text-success", bsicons::bs_icon("check-circle"),
             " Simulated: one mediating direction with alpha = beta = 1, so IE = 1.")
  )
}

gmed_describe <- function(d)
  sprintf("%d subjects; mediator dimension p = %d; exposure: %s.",
          length(d$M), ncol(d$M[[1]]), d$x_names[1])

gmed_run <- function(d, params) {
  stop_crt <- params$stop.crt %||% "nD"
  nD <- if (stop_crt == "nD") (params$nD %||% 1) else NULL
  fit <- med_fn("gmed")(d$X, d$M, d$Y, H = NULL,      # H = NULL -> average covariance
                        stop.crt = stop_crt, nD = nD,
                        DfD.thred = params$DfD.thred %||% 2,
                        ninitial = params$ninitial %||% 5,
                        seed = 100, verbose = FALSE)
  theta <- as.matrix(fit$theta)
  sims <- params$sims %||% 0
  inference <- NULL
  if (sims > 0) {
    inference <- lapply(seq_len(ncol(theta)), function(j)
      tryCatch(med_fn("gmed_boot")(d$X, d$M, d$Y, theta = theta[, j], H = NULL,
                                   boot = TRUE, sims = sims,
                                   boot.ci.type = "se", verbose = FALSE),
               error = function(e) NULL))
  }
  # projected log-variance per subject, for the mediator -> outcome plot
  score <- vapply(d$M, function(Mi)
    as.numeric(t(theta[, 1]) %*% stats::cov(Mi) %*% theta[, 1]), numeric(1))
  list(fit = fit, theta = theta, inference = inference, score = score,
       d = d, params = params)
}

gmed_summarize <- function(res) {
  fit <- res$fit; theta <- res$theta; nD <- ncol(theta)
  dirs <- paste0("C", seq_len(nD)); out <- list()
  labs <- c(alpha = "alpha  exposure -> mediator",
            beta  = "beta   mediator -> outcome",
            gamma = "gamma  direct effect",
            IE    = "IE     indirect (alpha x beta)",
            DE    = "DE     direct effect")

  eff <- do.call(rbind, lapply(seq_len(nD), function(j) {
    bi <- if (!is.null(res$inference)) res$inference[[j]] else NULL
    if (!is.null(bi)) {
      df <- as.data.frame(bi$coef)
      data.frame(Direction = dirs[j],
                 Effect = labs[rownames(bi$coef)] %||% rownames(bi$coef),
                 df, check.names = FALSE, row.names = NULL)
    } else {
      cf <- fit$coef[, j]
      data.frame(Direction = dirs[j], Effect = labs[rownames(fit$coef)],
                 Estimate = cf, SE = NA, statistics = NA, pvalue = NA,
                 LB = NA, UB = NA, check.names = FALSE, row.names = NULL)
    }
  }))
  out[["Mediation effects"]] <- eff

  if (!is.null(res$inference) && !is.null(res$inference[[1]]$coef.other)) {
    oth <- do.call(rbind, lapply(seq_len(nD), function(j) {
      bi <- res$inference[[j]]; if (is.null(bi)) return(NULL)
      data.frame(Direction = dirs[j], Term = rownames(bi$coef.other),
                 as.data.frame(bi$coef.other), check.names = FALSE,
                 row.names = NULL)
    }))
    if (!is.null(oth)) out[["Other coefficients"]] <- oth
  }

  th <- as.data.frame(theta); names(th) <- dirs
  out[["theta mediator loadings"]] <- cbind(Node = res$d$m_names, th)

  if (!is.null(fit$DfD)) {
    v <- if (is.list(fit$DfD)) as.numeric(fit$DfD$avg.level) else as.numeric(fit$DfD)
    if (length(v) < nD) v <- c(v, rep(NA, nD - length(v)))
    out[["DfD across dimensions"]] <- data.frame(Dimension = dirs,
                                                 DfD = v[seq_len(nD)])
  }
  if (nD > 1) {
    o <- as.data.frame(round(t(theta) %*% theta, 6)); names(o) <- dirs
    out[["Orthogonality theta'theta"]] <- cbind(" " = dirs, o)
  }

  if (!is.null(res$d$truth)) {
    tr <- res$d$truth
    j <- which.max(vapply(seq_len(nD),
                          function(k) abs_cos(theta[, k], tr$theta), numeric(1)))
    cf <- fit$coef[, j]
    out[["Truth vs estimate"]] <- data.frame(
      Quantity = c("theta cosine", "alpha (exp -> med)", "beta (med -> out)",
                   "gamma (direct)", "IE (alpha x beta)"),
      True = c(1, tr$alpha, tr$beta, tr$gamma, tr$IE),
      Estimate = c(abs_cos(theta[, j], tr$theta),
                   cf["alpha"], cf["beta"], cf["gamma"], cf["IE"]),
      check.names = FALSE)
  }
  out
}

gmed_plots <- function(res) {
  theta <- res$theta; nD <- ncol(theta)
  dirs <- paste0("C", seq_len(nD)); plots <- list()

  long <- do.call(rbind, lapply(seq_len(nD), function(j)
    data.frame(Node = res$d$m_names, Direction = dirs[j], Loading = theta[, j])))
  long$Node <- factor(long$Node, levels = res$d$m_names)
  plots[["theta loadings"]] <- list(
    plot = plotly::plot_ly(long, x = ~Node, y = ~Loading, color = ~Direction,
                           type = "bar") |>
      plotly::layout(title = "Mediator projection loadings (theta)",
                     barmode = "group"),
    data = long)

  if (!is.null(res$inference) && !is.null(res$inference[[1]])) {
    bi <- res$inference[[1]]$coef
    keep <- intersect(c("alpha", "beta", "gamma", "IE"), rownames(bi))
    edf <- data.frame(Effect = factor(keep, levels = keep),
                      Estimate = bi[keep, "Estimate"],
                      LB = bi[keep, "LB"], UB = bi[keep, "UB"])
    plots[["Effects with CI"]] <- list(
      plot = plotly::plot_ly(edf, x = ~Effect, y = ~Estimate, type = "bar",
                             marker = list(color = "#2c7fb8"),
                             error_y = list(type = "data", symmetric = FALSE,
                                            array = edf$UB - edf$Estimate,
                                            arrayminus = edf$Estimate - edf$LB)) |>
        plotly::layout(title = "Mediation effects (C1) with bootstrap 95% CI",
                       yaxis = list(title = "estimate")),
      data = edf)
  }

  sdf <- data.frame(id = res$d$ids, log_variance = log(pmax(res$score, 1e-12)),
                    Y = res$d$Y, Exposure = factor(res$d$X[, 1]))
  plots[["Mediator -> outcome"]] <- list(
    plot = plotly::plot_ly(sdf, x = ~log_variance, y = ~Y, color = ~Exposure,
                           type = "scatter", mode = "markers",
                           marker = list(size = 8, opacity = 0.7),
                           text = ~id, hovertemplate = "%{text}<extra></extra>") |>
      plotly::layout(title = "Outcome vs projected log-variance, by exposure",
                     xaxis = list(title = "log(theta' Sigma theta)"),
                     yaxis = list(title = "outcome Y")),
    data = sdf)
  plots
}

register_method(list(
  id = "gmed",
  name = "GMed",
  full_name = "Mediation Analysis with a Graph (Covariance) Mediator",
  short = "The mediator is a subject-level covariance matrix (for example a connectivity matrix). A learned projection theta reduces it to a scalar log-variance that carries the exposure -> mediator -> outcome effect.",
  status = "ready",
  tags = c("covariance mediator", "graph mediator", "bootstrap"),
  paper = list(
    citation = "Xu, Y., & Zhao, Y. (2025). Mediation analysis with graph mediator. Biostatistics, 26(1), kxaf004.",
    url = "https://doi.org/10.1093/biostatistics/kxaf004"),
  explain = file.path(APP_DIR, "methods", "gmed", "explain.md"),
  example_note = paste("100 subjects, mediator dimension p = 10, 500 rows per",
                       "subject, one binary exposure. One mediating direction",
                       "with alpha = beta = 1, so the true IE = 1."),
  data_inputs = list(
    list(id = "M", label = "Mediator (list of T x p matrices)",
         help = "An .rds / .RData holding a list of length n whose i-th element is a T_i x p matrix; subject i's mediator is its covariance. A long CSV (subject-id column + p variable columns) also works."),
    list(id = "X", label = "Exposure / covariates (n x nX)",
         help = "One row per subject, in the same order as the mediator list. The FIRST column is the exposure of interest; any further columns are adjusted for. Do NOT include an intercept -- it is added internally."),
    list(id = "Y", label = "Outcome (n values)",
         help = "One numeric value per subject, in the same order as the mediator list.")
  ),
  params = list(
    list(id = "stop.crt", label = "Number of directions chosen by", type = "select",
         choices = c("Fixed number (nD)" = "nD", "DfD criterion" = "DfD"),
         default = "nD",
         help = "Fix the count, or keep adding directions while deviation-from-diagonality stays below the threshold."),
    list(id = "nD", label = "Number of directions (nD)", type = "integer",
         default = 1, min = 1, max = 5,
         help = "Used with the fixed-number criterion. Each direction adds a bootstrap pass."),
    list(id = "DfD.thred", label = "DfD threshold", type = "numeric",
         default = 2, min = 1, max = 100, step = 0.5),
    list(id = "ninitial", label = "Random initialisations", type = "integer",
         default = 5, min = 1, max = 20,
         help = "Restarts for the theta search; more is safer but slower."),
    list(id = "sims", label = "Bootstrap replicates", type = "integer",
         default = 200, min = 0, max = 2000,
         help = "Subject-level bootstrap for alpha / beta / IE inference; 0 skips it.")
  ),
  example = gmed_app_example,
  export_example = function(d) list(
    M = d$M, X = matrix_to_df(d$X, id = d$ids), Y = data.frame(id = d$ids, Y = d$Y)),
  parse = gmed_app_parse,
  describe_data = gmed_describe,
  run = gmed_run,
  summarize = gmed_summarize,
  plots = gmed_plots
))
