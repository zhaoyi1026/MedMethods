# =============================================================================
# HeterMed -- heterogeneous (moderated) mediation effects
# -----------------------------------------------------------------------------
#   M_i = Z_i'(alpha0 + X_i alpha1) + e_i
#   Y_i = Z_i'(gamma0 + X_i gamma1) + (beta0 + X_i beta1) M_i + u_i
# Both paths are moderated by Z, so every subject has their own indirect effect.
# =============================================================================

hetermed_app_example <- function() {
  d <- med_fn("hetermed_example")(n = 600L)
  out <- list(X = d$X, M = d$M, Y = d$Y, Z = d$Z, truth = d$truth,
              z_names = colnames(d$Z) %||% paste0("Z", seq_len(ncol(d$Z))))
  out$preview_ui <- hetermed_preview(out)
  out
}

hetermed_app_parse <- function(files, opts) {
  X <- as_num_vector(files$X, "X (treatment)")
  n <- length(X)
  M <- as_num_vector(files$M, "M (mediator)", n = n)
  Y <- as_num_vector(files$Y, "Y (outcome)", n = n)
  Z <- as_num_matrix(files$Z, "Z", n = n)
  if (!all(sort(unique(X)) %in% c(-1, 1)))
    stop("The treatment X must be coded -1 / +1. Found: ",
         paste(utils::head(sort(unique(X)), 5), collapse = ", "), ".")
  # the model needs a leading intercept column in Z; add one if absent
  if (!all(Z[, 1] == 1)) {
    Z <- cbind(Intercept = 1, Z)
  } else {
    colnames(Z)[1] <- "Intercept"
  }
  out <- list(X = X, M = M, Y = Y, Z = Z, truth = NULL, z_names = colnames(Z))
  out$preview_ui <- hetermed_preview(out)
  out
}

hetermed_preview <- function(d) {
  tagList(
    bslib::layout_columns(
      col_widths = c(3, 3, 3, 3),
      bslib::value_box("Subjects", length(d$X), theme = "primary"),
      bslib::value_box("Moderators", ncol(d$Z) - 1L, theme = "secondary"),
      bslib::value_box("Treated (+1)", sum(d$X == 1), theme = "secondary"),
      bslib::value_box("Control (-1)", sum(d$X == -1), theme = "secondary")
    ),
    tags$p(class = "small text-muted mt-2",
           sprintf("Moderators: %s (the first column is the intercept).",
                   paste(setdiff(d$z_names, "Intercept"), collapse = ", "))),
    if (!is.null(d$truth))
      tags$p(class = "small text-success",
             bsicons::bs_icon("check-circle"),
             " Simulated: true beta0 = 0.8, beta1 = 0.2, alpha1 = (0.5, 0.4, 0).")
  )
}

hetermed_describe <- function(d)
  sprintf("%d subjects; %d moderator(s); treatment coded -1/+1.",
          length(d$X), ncol(d$Z) - 1L)

hetermed_run <- function(d, params) {
  meth <- params$method %||% "OLS"
  fit <- med_fn("hetermed")(d$X, d$M, d$Y, d$Z, method = meth)
  inf <- NULL
  if (isTRUE(params$inference)) {
    inf <- tryCatch(
      med_fn("hetermed_inf")(d$X, d$M, d$Y, d$Z, fit, method = meth,
                             conf.level = params$conf.level %||% 0.95),
      error = function(e) structure(conditionMessage(e), class = "inf_failed"))
  }
  list(fit = fit, inf = inf, d = d, params = params)
}

hetermed_summarize <- function(res) {
  fit <- res$fit; inf <- res$inf; zn <- res$d$z_names; out <- list()
  ok <- !is.null(inf) && !inherits(inf, "inf_failed")

  tab <- function(nm, est, lab) {
    if (ok && !is.null(inf[[nm]])) {
      df <- as.data.frame(inf[[nm]])
      df <- cbind(Term = rownames(inf[[nm]]) %||% lab, df)
      rownames(df) <- NULL
      df
    } else {
      data.frame(Term = lab, Estimate = as.numeric(est), check.names = FALSE)
    }
  }
  out[["a-path: alpha0 (main)"]]      <- tab("alpha0", fit$alpha0, zn)
  out[["a-path: alpha1 (moderated)"]] <- tab("alpha1", fit$alpha1, zn)
  out[["b-path: beta0 / beta1"]] <- rbind(
    tab("beta0", fit$beta0, "M"), tab("beta1", fit$beta1, "M x X"))
  out[["direct: gamma0 (main)"]]      <- tab("gamma0", fit$gamma0, zn)
  out[["direct: gamma1 (moderated)"]] <- tab("gamma1", fit$gamma1, zn)

  ite <- as.data.frame(fit$ITE)
  ite <- cbind(Subject = seq_len(nrow(ite)), ite)
  if (ok && !is.null(inf$NIE)) {
    ite$NIE_SE <- inf$NIE$SE; ite$NIE_LB <- inf$NIE$LB; ite$NIE_UB <- inf$NIE$UB
    ite$NDE_SE <- inf$NDE$SE; ite$NDE_LB <- inf$NDE$LB; ite$NDE_UB <- inf$NDE$UB
  }
  out[["Per-subject effects"]] <- ite

  out[["Effect summary"]] <- data.frame(
    Quantity = c("mean NIE", "SD of NIE", "mean NDE", "SD of NDE"),
    Value = c(mean(fit$ITE[, "NIE"]), stats::sd(fit$ITE[, "NIE"]),
              mean(fit$ITE[, "NDE"]), stats::sd(fit$ITE[, "NDE"])),
    check.names = FALSE)

  if (!is.null(res$d$truth)) {
    tr <- res$d$truth
    out[["Truth vs estimate"]] <- data.frame(
      Quantity = c("beta0", "beta1",
                   paste0("alpha1[", zn, "]"), "mean NIE",
                   "correlation of per-subject NIE with truth"),
      True = c(tr$beta0, tr$beta1, tr$alpha1, mean(tr$NIE), 1),
      Estimate = c(fit$beta0, fit$beta1, fit$alpha1,
                   mean(fit$ITE[, "NIE"]),
                   stats::cor(fit$ITE[, "NIE"], tr$NIE)),
      check.names = FALSE)
  }
  if (!is.null(inf) && inherits(inf, "inf_failed"))
    out[["Inference note"]] <- data.frame(Message = as.character(inf))
  out
}

hetermed_plots <- function(res) {
  fit <- res$fit; inf <- res$inf; plots <- list()
  ok <- !is.null(inf) && !inherits(inf, "inf_failed")
  ite <- as.data.frame(fit$ITE)

  # (a) distribution of the per-subject indirect effect
  hdf <- data.frame(NIE = ite$NIE, Treatment = factor(ite$Treatment))
  plots[["Indirect effect spread"]] <- list(
    plot = plotly::plot_ly(hdf, x = ~NIE, color = ~Treatment, type = "histogram",
                           nbinsx = 40, opacity = 0.75) |>
      plotly::layout(title = "Per-subject natural indirect effect",
                     barmode = "overlay",
                     xaxis = list(title = "NIE"), yaxis = list(title = "subjects"),
                     shapes = list(list(type = "line", y0 = 0, y1 = 1,
                                        yref = "paper",
                                        x0 = mean(ite$NIE), x1 = mean(ite$NIE),
                                        line = list(dash = "dot")))),
    data = hdf)

  # (b) estimated vs true per-subject NIE (example data only)
  if (!is.null(res$d$truth)) {
    sdf <- data.frame(True = res$d$truth$NIE, Estimate = ite$NIE,
                      Treatment = factor(ite$Treatment))
    r <- stats::cor(sdf$True, sdf$Estimate)
    rng <- range(c(sdf$True, sdf$Estimate))
    plots[["Recovery of per-subject NIE"]] <- list(
      plot = plotly::plot_ly(sdf, x = ~True, y = ~Estimate, color = ~Treatment,
                             type = "scatter", mode = "markers",
                             marker = list(size = 6, opacity = 0.7)) |>
        plotly::add_trace(x = rng, y = rng, inherit = FALSE, type = "scatter",
                          mode = "lines", line = list(dash = "dash", color = "grey"),
                          showlegend = FALSE) |>
        plotly::layout(title = sprintf("Estimated vs true individual NIE (r = %.3f)", r),
                       xaxis = list(title = "true NIE"),
                       yaxis = list(title = "estimated NIE")),
      data = sdf)
  }

  # (c) moderated a-path coefficients with CIs
  if (ok && !is.null(inf$alpha1)) {
    a <- as.data.frame(inf$alpha1)
    a$Term <- rownames(inf$alpha1) %||% res$d$z_names
    plots[["Moderated a-path (alpha1)"]] <- list(
      plot = plotly::plot_ly(a, x = ~Term, y = ~Estimate, type = "bar",
                             marker = list(color = "#2c7fb8"),
                             error_y = list(type = "data", symmetric = FALSE,
                                            array = a$UB - a$Estimate,
                                            arrayminus = a$Estimate - a$LB)) |>
        plotly::layout(title = "alpha1: moderation of the treatment -> mediator path",
                       yaxis = list(title = "estimate")),
      data = a)
  }

  # (d) the moderated b-path: outcome vs mediator by treatment arm
  bdf <- data.frame(M = res$d$M, Y = res$d$Y, Treatment = factor(res$d$X))
  plots[["Moderated b-path"]] <- list(
    plot = plotly::plot_ly(bdf, x = ~M, y = ~Y, color = ~Treatment,
                           type = "scatter", mode = "markers",
                           marker = list(size = 6, opacity = 0.6)) |>
      plotly::layout(title = sprintf("Outcome vs mediator by arm (slope %.3f vs %.3f)",
                                     fit$beta0 + fit$beta1, fit$beta0 - fit$beta1),
                     xaxis = list(title = "mediator M"),
                     yaxis = list(title = "outcome Y")),
    data = bdf)
  plots
}

register_method(list(
  id = "hetermed",
  name = "HeterMed",
  full_name = "Heterogeneous Causal Mediation Effects",
  short = "Both the treatment -> mediator and mediator -> outcome paths are moderated by covariates, so each subject has their own natural indirect and direct effect.",
  status = "ready",
  tags = c("moderated mediation", "heterogeneous effects", "asymptotic inference"),
  paper = list(
    citation = "Zhao, Y., Li, C., & Tu, W. (2025). Estimation of Heterogeneous Causal Mediation Effects in a Hypertension Treatment Trial. arXiv preprint arXiv:2512.12043.",
    url = "https://doi.org/10.48550/arXiv.2512.12043"),
  explain = file.path(APP_DIR, "methods", "hetermed", "explain.md"),
  example_note = paste("600 subjects, treatment coded -1/+1, two moderators.",
                       "True beta0 = 0.8, beta1 = 0.2, alpha1 = (0.5, 0.4, 0)."),
  data_inputs = list(
    list(id = "X", label = "Treatment (n values, coded -1 / +1)",
         help = "One value per subject, coded -1 and +1. A CSV column or an .rds numeric vector."),
    list(id = "M", label = "Mediator (n values)",
         help = "One numeric value per subject, in the same order as the treatment."),
    list(id = "Y", label = "Outcome (n values)",
         help = "One numeric value per subject, in the same order as the treatment."),
    list(id = "Z", label = "Moderators (n x p)",
         help = "An n x p numeric matrix of covariates that moderate both paths. A leading column of 1s is added automatically if you do not supply one.")
  ),
  params = list(
    list(id = "method", label = "Estimator", type = "select",
         choices = c("OLS" = "OLS", "Generalized lasso" = "genlasso"),
         default = "OLS",
         help = "OLS fits both arms unpenalised. The generalized-lasso option shrinks the moderation toward a common effect (slower; bootstrap inference)."),
    list(id = "inference", label = "Compute standard errors and CIs",
         type = "checkbox", default = TRUE,
         help = "Asymptotic inference for OLS; bootstrap for the generalized lasso."),
    list(id = "conf.level", label = "Confidence level", type = "numeric",
         default = 0.95, min = 0.5, max = 0.999, step = 0.01)
  ),
  example = hetermed_app_example,
  export_example = function(d) list(
    X = data.frame(X = d$X), M = data.frame(M = d$M), Y = data.frame(Y = d$Y),
    Z = matrix_to_df(d$Z)),
  parse = hetermed_app_parse,
  describe_data = hetermed_describe,
  run = hetermed_run,
  summarize = hetermed_summarize,
  plots = hetermed_plots
))
