##############################################################################
# Verify each *_example() generator feeds its wrapper and recovers the truth.
# These are the exact settings used in README.md.
##############################################################################
suppressMessages(library(MedMethods))
f <- function(x, d = 3) paste(format(round(as.numeric(x), d), nsmall = d), collapse = " ")
abscos <- function(a, b) abs(sum(a * b)) / sqrt(sum(a^2) * sum(b^2))
quiet <- function(expr) {
  lf <- tempfile(); con <- file(lf, open = "wt")
  sink(con); sink(con, type = "message")
  v <- tryCatch(suppressWarnings(expr),
                error = function(e) structure(conditionMessage(e), class = "efail"))
  sink(type = "message"); sink(); close(con); unlink(lf); v
}
rep2 <- function(label, txt) cat(sprintf("  %-28s %s\n", label, txt))
sec  <- function(x) cat(sprintf("\n=== %s %s\n", x, strrep("=", max(0, 50 - nchar(x)))))
# "selected" = clearly non-zero relative to the largest effect
nz <- function(v, frac = 0.05) which(abs(v) > frac * max(abs(v)))

sec("macc_example (single)")
d <- macc_example()
fit <- quiet(macc(d$dat, model.type = "single", delta = d$truth$delta))
rep2("A / B / C  (0.5,-1,0.5)", f(fit$Coefficients[c("A","B","C"), "Estimate"]))
rep2("ABp  (true -0.5)", f(fit$Coefficients["ABp", "Estimate"]))

sec("macc_example (twolevel)")
d <- macc_example("twolevel", N = 50L, n.trial = 100L)
fit <- quiet(macc(d$dat, model.type = "twolevel", method = "HL", delta = 0.5))
if (inherits(fit, "efail")) rep2("FAIL", fit) else {
  rep2("A / B / C  (0.5,-1,0.5)", f(fit$Coefficients[c("A","B","C"), "Estimate"]))
  rep2("AB.prod  (true -0.5)", f(fit$Coefficients["AB.prod", "Estimate"]))
}

sec("gma_example")
d <- gma_example()
for (pp in 1:3) {
  fit <- quiet(gma(d$dat, model.type = "single", p = pp,
                   delta = d$truth$delta, single.var.asmp = TRUE))
  if (inherits(fit, "efail")) rep2(paste0("p=", pp, " FAIL"), fit) else
    rep2(paste0("p=", pp, " A/B/ABp (.5,-1,-.5)"),
         f(fit$Coefficients[c("A","B","AB.p"), "Estimate"]))
}

sec("spcma_example")
d <- spcma_example()
fit <- quiet(spcma(d$X, d$M, d$Y, adaptive = TRUE, var.per = 0.8, boot = FALSE,
                   PC.run = TRUE))
if (inherits(fit, "efail")) rep2("FAIL", fit) else {
  W <- fit$SPCA$W
  rep2("n components", ncol(W))
  rep2("loading cos to Phi 1/2/3",
       f(vapply(1:min(3, ncol(W)), function(j) abscos(W[, j], d$truth$Phi[, j]), 0)))
  rep2("IE 1/2/3  (true 4,-1.5,0)", f(fit$SPCA$IE[1:3, "Estimate"]))
}

sec("cfma_example")
d <- cfma_example()
fit <- quiet(cfma_concurrent(d$Z, d$M, d$Y, intercept = FALSE, timeinv = d$timeinv))
if (inherits(fit, "efail")) rep2("FAIL", fit) else {
  rep2("cor: alpha / beta / IE",
       f(c(stats::cor(fit$M$curve[1, ], d$truth$alpha),
           stats::cor(fit$Y$curve[2, ], d$truth$beta),
           stats::cor(fit$IE$curve,     d$truth$IE))))
  rep2("gamma(t) RMSE (range 0-0.8)",
       f(sqrt(mean((fit$Y$curve[1, ] - d$truth$gamma)^2))))
}

sec("pathlasso_example")
d <- pathlasso_example()
fit <- quiet(pathlasso(d$X, d$M, d$Y, lambda = 0.001, omega = 0, phi = 1,
                       max.itr = 3000, tol = 1e-8))
if (inherits(fit, "efail")) rep2("FAIL", fit) else {
  sel <- nz(fit$AB)
  rep2("selected mediators", paste(sel, collapse = ","))
  rep2("true signal", paste(d$truth$signal, collapse = ","))
  rep2("AB 1:4  (true 4,-4,-4,4)", f(fit$AB[1:4]))
}

sec("pathlasso2b_example")
d <- pathlasso2b_example()
fit <- quiet(pathlasso2b(d$X, d$M1, d$M2, d$Y, kappa1 = 5, kappa2 = 5, kappa3 = 5,
                         kappa4 = 5, nu1 = 2, nu2 = 2, mu1 = 2, mu2 = 2,
                         max.itr = 3000))
if (inherits(fit, "efail")) rep2("FAIL", fit) else {
  t1 <- which(abs(d$truth$IE.M1) > 0); t2 <- which(abs(d$truth$IE.M2) > 0)
  s1 <- nz(fit$IE.M1); s2 <- nz(fit$IE.M2)
  rep2("block 1 selected / true", paste(paste(s1, collapse = ","), "/",
                                        paste(t1, collapse = ",")))
  rep2("block 2 selected / true", paste(paste(s2, collapse = ","), "/",
                                        paste(t2, collapse = ",")))
}

sec("hdmediation_example")
d <- hdmediation_example()
fit <- quiet(hdmediation(d$X, d$M, d$Y, lambda = 2, pi = 0.5, phi = 2,
                         delta = 0.5, max.itr = 1000))
if (inherits(fit, "efail")) rep2("FAIL", fit) else {
  tp  <- d$truth$signal
  sel <- which(abs(fit$IE) > 1e-3, arr.ind = TRUE)
  hit <- sum(apply(tp, 1, function(r) any(sel[, 1] == r[1] & sel[, 2] == r[2])))
  rep2("true exposure-mediator paths", paste(apply(tp, 1, paste, collapse = "->"),
                                             collapse = " "))
  rep2("recovered / total selected", paste(hit, "of", nrow(tp), "|", nrow(sel),
                                           "non-zero of", length(fit$IE)))
}

sec("pcma_example")
d <- pcma_example()
fit <- quiet(pcma(d$X, d$M, d$Y, stop.crt = "nD", nD = 2, boot = FALSE, ninitial = 3))
if (inherits(fit, "efail")) rep2("FAIL", fit) else {
  rep2("Phi cos (dir 1,2)", f(c(abscos(fit$Phi[, 1], d$truth$Phi[, 1]),
                                abscos(fit$Phi[, 2], d$truth$Phi[, 2]))))
  rep2("Psi cos (dir 1,2)", f(c(abscos(fit$Psi[, 1], d$truth$Psi[, 1]),
                                abscos(fit$Psi[, 2], d$truth$Psi[, 2]))))
  co <- quiet(pcma_coef(d$X, d$M, d$Y, phi = fit$Phi, psi = fit$Psi))
  if (!inherits(co, "efail")) {
    # projections are identified only up to sign, so compare magnitudes
    rep2("|alpha| diag (true 2,2)", f(abs(diag(as.matrix(co$alpha)))))
    rep2("|beta|  (true 2,1)", f(abs(as.numeric(co$beta))))
  }
}

sec("gmed_example")
d <- gmed_example()
fit <- quiet(gmed(d$X, d$M, d$Y, stop.crt = "nD", nD = 1, ninitial = 5,
                  verbose = FALSE))
if (inherits(fit, "efail")) rep2("FAIL", fit) else {
  rep2("theta cos to truth", f(abscos(fit$theta, d$truth$theta)))
  rep2("alpha / beta / IE (all 1)", f(fit$coef[c("alpha","beta","IE"), 1]))
  bt <- quiet(gmed_boot(d$X, d$M, d$Y, theta = fit$theta, sims = 200,
                        verbose = FALSE))
  if (!inherits(bt, "efail"))
    rep2("boot IE est (p-value)", paste0(f(bt$coef["IE", "Estimate"]), " (",
                                         f(bt$coef["IE", "pvalue"]), ")"))
}

sec("hetermed_example")
d <- hetermed_example()
fit <- quiet(hetermed(d$X, d$M, d$Y, d$Z, method = "OLS"))
if (inherits(fit, "efail")) rep2("FAIL", fit) else {
  rep2("beta0 / beta1 (0.8, 0.2)", f(c(fit$beta0, fit$beta1)))
  rep2("alpha1  (true .5,.4,0)", f(fit$alpha1))
  ite <- quiet(hetermed_ite(d$X, d$Z, fit$alpha0, fit$alpha1, fit$beta0,
                            fit$beta1, fit$gamma0, fit$gamma1))
  rep2("NIE cor with true NIE", f(stats::cor(ite[, "NIE"], d$truth$NIE)))
  rep2("mean NIE est / true", paste(f(mean(ite[, "NIE"])), "/",
                                    f(mean(d$truth$NIE))))
}
cat("\nDONE\n")
