##############################################################################
# Phase 1 acceptance test: every exported method wrapper runs end-to-end.
##############################################################################
suppressMessages(library(MedMethods))
ROOT <- "/Users/yizhao/Dropbox/MyFolder/Biostat-IU/Projects/Yi/mediation/code"

# Load only the DATA objects out of a workspace dump (these .RData files also
# contain the authors' functions, which must not shadow the package).
dat_from <- function(relpath, want) {
  e <- new.env()
  load(file.path(ROOT, relpath), envir = e)
  setNames(lapply(want, function(o) get(o, envir = e)), want)
}

RESULTS <- list()
run <- function(label, expr, note = function(v) "") {
  t0 <- proc.time()[["elapsed"]]
  v <- NULL
  # several methods print progress chatter on both streams; keep the report readable
  logf <- tempfile(); con <- file(logf, open = "wt")
  sink(con); sink(con, type = "message")
  v <- tryCatch(suppressWarnings(expr),
                error = function(e) structure(conditionMessage(e), class = "smokefail"))
  sink(type = "message"); sink(); close(con); unlink(logf)
  el <- proc.time()[["elapsed"]] - t0
  ok <- !inherits(v, "smokefail")
  RESULTS[[length(RESULTS) + 1L]] <<- list(label = label, ok = ok, secs = el,
    note = if (ok) tryCatch(note(v), error = function(e) "<note failed>") else as.character(v))
  cat(sprintf("%-34s %-4s %6.1fs  %s\n", label, if (ok) "OK" else "FAIL", el,
              RESULTS[[length(RESULTS)]]$note))
  invisible(v)
}
fmt <- function(x, d = 3) paste(format(round(as.numeric(x), d), nsmall = d), collapse = " ")

cat("\n== macc =====================================================\n")
d <- dat_from("macc/macc/data/env.single.rda", "env.single")$env.single
macc_single <- get("data1", d)
run("macc(single, delta=0.5)",
    macc(macc_single, model.type = "single", delta = 0.5),
    function(v) paste("ABp =", fmt(v$Coefficients["ABp", "Estimate"])))
d2 <- dat_from("macc/macc/data/env.two.rda", "env.two")$env.two
macc_two <- get("data2", d2)
run("macc(twolevel, HL, delta=0.5)",
    macc(macc_two, model.type = "twolevel", method = "HL", delta = 0.5),
    function(v) paste("AB.prod =", fmt(v$Coefficients["AB.prod", "Estimate"])))
run("macc_sim_single()", {
  Sigma <- matrix(c(1, 1, 1, 4), 2, 2); Theta <- matrix(c(0.5, 0, 0.5, -1), 2, 2)
  set.seed(100); Z <- matrix(rbinom(100, 1, 0.5), 100, 1)
  set.seed(5000); macc_sim_single(Z, Theta, Sigma)
}, function(v) paste("names =", paste(names(v), collapse = ",")))

cat("\n== gma ======================================================\n")
g <- dat_from("gma/gma/data/env.single.rda", "env.single")$env.single
gma_single <- get("data1", g)
gma_note <- function(v) paste0("A = ", fmt(v$Coefficients["A", "Estimate"]),
                               ", AB.p = ", fmt(v$Coefficients["AB.p", "Estimate"]),
                               " (SE ", fmt(v$Coefficients["AB.p", "SE"]), ")")
run("gma(single, p=1, var.asmp=TRUE)",
    gma(gma_single, model.type = "single", p = 1, single.var.asmp = TRUE), gma_note)
run("gma(single, p=2, var.asmp=TRUE) FIX",
    gma(gma_single, model.type = "single", p = 2, single.var.asmp = TRUE), gma_note)
run("gma(single, p=3, var.asmp=TRUE) FIX",
    gma(gma_single, model.type = "single", p = 3, single.var.asmp = TRUE), gma_note)
run("gma(single, p=4, var.asmp=TRUE) FIX",
    gma(gma_single, model.type = "single", p = 4, single.var.asmp = TRUE), gma_note)
run("gma(single, p=2, var.asmp=FALSE)",
    gma(gma_single, model.type = "single", p = 2, single.var.asmp = FALSE), gma_note)

cat("\n== spcma ====================================================\n")
s <- dat_from("spcma/spcma-master/data/env.example.rda", "env.example")$env.example
sX <- get("X", s); sM <- get("M", s); sY <- get("Y", s)
run("spcma(adaptive, no boot)",
    spcma(sX, sM, sY, adaptive = TRUE, var.per = 0.75, boot = FALSE, PC.run = FALSE),
    function(v) paste("IE rows =", nrow(v$SPCA$IE)))
run("mcma_pca(no boot)",
    mcma_pca(sX, sM, sY, adaptive = TRUE, var.per = 0.75, boot = FALSE),
    function(v) paste("IE rows =", nrow(v$IE)))
run("mcma_bk(sims=50)",
    mcma_bk(sX, sM, sY, sims = 50, boot = TRUE),
    function(v) paste("IE rows =", nrow(v$IE)))

cat("\n== cfma =====================================================\n")
cc <- dat_from("cfma/cfma-master/data/env.concurrent.rda", "env.concurrent")$env.concurrent
cZ <- get("Z", cc); cM <- get("M", cc); cY <- get("Y", cc)
run("cfma_concurrent()",
    cfma_concurrent(cZ, cM, cY, intercept = FALSE, timeinv = c(0, 300)),
    function(v) paste("IE curve length =", length(v$IE$curve)))
ch <- dat_from("cfma/cfma-master/data/env.historical.rda", "env.historical")$env.historical
run("cfma_historical()",
    cfma_historical(get("Z", ch), get("M", ch), get("Y", ch), intercept = FALSE,
                    timeinv = c(0, 300)),
    function(v) paste("components =", paste(names(v), collapse = ",")))

cat("\n== pathway lasso ============================================\n")
pl <- dat_from("pathway_lasso/V1/eg.RData", c("X", "M", "Y"))
run("pathlasso(lambda=1)",
    pathlasso(pl$X, pl$M, pl$Y, lambda = 1, omega = 0, phi = 1, max.itr = 300),
    function(v) paste("nonzero AB =", sum(abs(v$AB) > 1e-3), "of", length(v$AB),
                      "| max|AB| =", fmt(max(abs(v$AB)))))
run("pathlasso(lambda=0.5, omega=0.1)",
    pathlasso(pl$X, pl$M, pl$Y, lambda = 0.5, omega = 0.1, phi = 1, max.itr = 300),
    function(v) paste("nonzero AB =", sum(abs(v$AB) > 1e-3), "of", length(v$AB),
                      "| max|AB| =", fmt(max(abs(v$AB)))))
run("pathlasso_ksc(3 lambdas)",
    pathlasso_ksc(pl$X, pl$M, pl$Y, lambda = c(0.5, 1, 2), n.rep = 2, max.itr = 200),
    function(v) paste("components =", paste(head(names(v), 5), collapse = ",")))
run("pathlasso_sim()", {
  set.seed(1); k <- 8; n <- 40
  A <- matrix(c(1, -1, 1, rep(0, k - 3)), nrow = 1)
  B <- matrix(c(1, 1, -1, rep(0, k - 3)), ncol = 1)
  pathlasso_sim(n, matrix(rbinom(n, 1, 0.5), n, 1), A, B, 1,
                Delta = diag(k), Xi1 = diag(k), Sigma2 = matrix(1, 1, 1))
}, function(v) paste("names =", paste(names(v), collapse = ",")))

cat("\n== multimodal ===============================================\n")
mm <- dat_from("multimodal/multimodal_integration-master/example.RData",
               c("X", "M1", "M2", "Y"))
run("pathlasso2b(kappa=0.1)",
    pathlasso2b(mm$X, mm$M1, mm$M2, mm$Y, kappa1 = 0.1, kappa2 = 0.1,
                kappa3 = 0.1, kappa4 = 0.1, nu1 = 2, nu2 = 2, mu1 = 0.05, mu2 = 0.05,
                max.itr = 300),
    function(v) paste0("nonzero IE.M1 = ", sum(abs(v$IE.M1) > 1e-3), "/", length(v$IE.M1),
                       ", IE.M2 = ", sum(abs(v$IE.M2) > 1e-3), "/", length(v$IE.M2),
                       ", IE.M1M2 = ", sum(abs(v$IE.M1M2) > 1e-3)))

cat("\n== HD exposures + mediators =================================\n")
hd <- dat_from("HDExposureMediator/HDExposureMediator-main/example.RData",
               c("X", "M", "Y"))
# The published data is r = p = n = 100. The non-PCA estimator regresses on all
# r exposures directly, so it needs n > r; the PCA variant is what the paper uses
# at r = n. Take a 20-exposure / 20-mediator subset for the direct estimator.
hdX <- hd$X[, 1:20, drop = FALSE]; hdM <- hd$M[, 1:20, drop = FALSE]
run("hdmediation(r=20, p=20)",
    hdmediation(hdX, hdM, hd$Y, lambda = 1, pi = 0.5, phi = 2, delta = 0.5,
                max.itr = 200),
    function(v) paste0("nonzero IE = ", sum(abs(v$IE) > 1e-3), "/", length(v$IE),
                       ", nonzero gamma = ", sum(abs(v$gamma) > 1e-3)))
run("hdmediation_pca(r=100, adaptive)",
    hdmediation_pca(hd$X, hd$M, hd$Y, adaptive = TRUE, var.prop = 0.9,
                    lambda = 1, pi = 0.5, phi = 2, delta = 0.5, max.itr = 200),
    function(v) paste0("nonzero IE = ", sum(abs(v$IE) > 1e-3), "/", length(v$IE),
                       ", n.pc = ", if (!is.null(v$n.pc)) v$n.pc else NA))

cat("\n== PCMA =====================================================\n")
# Generated from the PCMA model itself: orthonormal exposure/mediator
# projections Phi and Psi, one mediating component with a-path alpha, b-path
# beta and direct effect gamma on the projected scores.
set.seed(20260729)
n <- 400; p <- 3; q <- 5
Phi_full <- qr.Q(qr(matrix(rnorm(p * p), p, p)))
Psi_full <- qr.Q(qr(matrix(rnorm(q * q), q, q)))
truth_phi <- Phi_full[, 1]; truth_psi <- Psi_full[, 1]
alpha_t <- 2; beta_t <- 1.5; gamma_t <- 0.5
pX <- matrix(rnorm(n * p), n, p)
u <- as.vector(pX %*% truth_phi)                       # exposure score
S <- matrix(rnorm(n * q, sd = 0.3), n, q)              # mediator scores
S[, 1] <- alpha_t * u + rnorm(n, sd = 0.3)             # mediating component
pM <- S %*% t(Psi_full)                                # so pM %*% Psi_full == S
pY <- matrix(beta_t * S[, 1] + gamma_t * u + rnorm(n, sd = 0.3), n, 1)
absco <- function(a, b) abs(sum(a * b) / sqrt(sum(a^2) * sum(b^2)))
pcma_fit <- run("pcma(nD=1, no boot)",
    pcma(pX, pM, pY, stop.crt = "nD", nD = 1, boot = FALSE, ninitial = 2),
    function(v) paste0("Phi cos = ", fmt(absco(v$Phi[, 1], truth_phi)),
                       ", Psi cos = ", fmt(absco(v$Psi[, 1], truth_psi))))
run("pcma(nD=1, boot, sims=50)",
    pcma(pX, pM, pY, stop.crt = "nD", nD = 1, boot = TRUE, sims = 50,
         ninitial = 2, verbose = FALSE),
    function(v) paste("inference rows =", nrow(v$coef.inference[[1]])))
run("pcma_coef(at truth)",
    pcma_coef(pX, pM, pY, phi = matrix(truth_phi, ncol = 1),
              psi = matrix(truth_psi, ncol = 1)),
    function(v) paste0("alpha = ", fmt(v$alpha), " (true ", alpha_t, "), beta = ",
                       fmt(v$beta), " (true ", beta_t, "), IE = ",
                       fmt(v$alpha * v$beta), " (true ", alpha_t * beta_t, ")"))
run("pcma_inf_asmp(at fit)",
    pcma_inf_asmp(pX, pM, pY, phi = pcma_fit$Phi, psi = pcma_fit$Psi),
    function(v) paste("components =", paste(head(names(v), 5), collapse = ",")))

cat("\n== GMed =====================================================\n")
gm <- dat_from("GMed/GMed-main/example.RData", c("X", "M", "Y", "H"))
gmed_fit <- run("gmed(nD=1)",
    gmed(gm$X, gm$M, gm$Y, H = gm$H, stop.crt = "nD", nD = 1, ninitial = 2,
         verbose = FALSE),
    function(v) paste0("alpha = ", fmt(v$coef["alpha", 1]),
                       ", beta = ", fmt(v$coef["beta", 1]),
                       ", IE = ", fmt(v$coef["IE", 1])))
# gmed_coef / gmed_refit are the two functions whose undefined n/p were patched.
# NB: despite its comment, CAPMediation_coef_Mcov() adds the intercept itself and
# names X's columns "X","W1",...; pass the exposure/covariate matrix without a 1s column.
run("gmed_coef(at fitted theta)",
    gmed_coef(gm$X, gm$M, gm$Y, theta = gmed_fit$theta),
    function(v) paste0("alpha = ", fmt(v$alpha), ", beta = ", fmt(v$beta),
                       ", IE = ", fmt(v$IE)))
run("gmed_refit(at fitted theta)",
    gmed_refit(gm$X, gm$M, gm$Y, Theta = as.matrix(gmed_fit$theta)),
    function(v) paste("components =", paste(head(names(v), 5), collapse = ",")))
run("gmed_boot(sims=50)",
    gmed_boot(gm$X, gm$M, gm$Y, theta = gmed_fit$theta, H = gm$H, sims = 50,
              verbose = FALSE),
    function(v) paste0("IE = ", fmt(v$coef["IE", "Estimate"]),
                       " (p = ", fmt(v$coef["IE", "pvalue"]), ")"))
run("gmed_refit_boot(sims=50)",
    gmed_refit_boot(gm$X, gm$M, gm$Y, Theta = as.matrix(gmed_fit$theta),
                    H = gm$H, sims = 50, verbose = FALSE),
    function(v) paste("components =", paste(head(names(v), 4), collapse = ",")))

cat("\n== heterogeneous mediation ==================================\n")
set.seed(20260729)
n <- 300
hX <- sample(c(-1, 1), n, replace = TRUE)
hZ <- cbind(Intercept = 1, Z1 = rnorm(n), Z2 = rbinom(n, 1, 0.5))
a0 <- c(0.2, 0.3, 0); a1 <- c(0.5, 0.4, 0)
g0 <- c(0.1, 0, 0);   g1 <- c(0.3, 0, 0)
b0 <- 0.8; b1 <- 0.2
hM <- as.vector(hZ %*% a0 + hX * (hZ %*% a1)) + rnorm(n, sd = 0.5)
hY <- as.vector(hZ %*% g0 + hX * (hZ %*% g1)) + (b0 + b1 * hX) * hM + rnorm(n, sd = 0.5)
run("hetermed(OLS)",
    hetermed(hX, hM, hY, hZ, method = "OLS"),
    function(v) paste("beta0 =", fmt(v$beta0), "beta1 =", fmt(v$beta1)))
run("hetermed_inf(OLS)", {
  fit <- hetermed(hX, hM, hY, hZ, method = "OLS")
  hetermed_inf(hX, hM, hY, hZ, fit, method = "OLS")
}, function(v) paste("components =", paste(head(names(v), 6), collapse = ",")))
run("hetermed(genlasso)",
    hetermed(hX, hM, hY, hZ, method = "genlasso"),
    function(v) paste("beta0 =", fmt(v$beta0), "beta1 =", fmt(v$beta1)))
run("hetermed_ite()", {
  fit <- hetermed(hX, hM, hY, hZ, method = "OLS")
  hetermed_ite(hX, hZ, fit$alpha0, fit$alpha1, fit$beta0, fit$beta1,
               fit$gamma0, fit$gamma1)
}, function(v) paste("mean NIE =", fmt(mean(v[, "NIE"])),
                     "mean NDE =", fmt(mean(v[, "NDE"]))))

cat("\n== med_methods() / med_internal() ===========================\n")
run("med_methods()", med_methods(), function(v) paste(length(v), "modules"))
run("med_internal('gmed')", med_internal("gmed"), function(v) paste(length(v), "functions"))

cat("\n============================================================\n")
nfail <- sum(!vapply(RESULTS, `[[`, logical(1), "ok"))
cat(sprintf("%d/%d checks passed", length(RESULTS) - nfail, length(RESULTS)), "\n")
if (nfail) {
  cat("\nFAILURES:\n")
  for (r in RESULTS) if (!r$ok) cat(sprintf(" - %-34s %s\n", r$label, r$note))
}
