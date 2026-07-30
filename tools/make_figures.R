##############################################################################
# Regenerate the recovery figures shown in README.md.
#
#   Rscript tools/make_figures.R [output_dir]
#
# Each panel fits a method to its own *_example() data and plots the estimates
# against the ground truth the generator returns.
##############################################################################
suppressMessages(library(MedMethods))

args <- commandArgs(trailingOnly = TRUE)
OUT <- if (length(args) >= 1) args[1] else "man/figures"
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)
cat("writing figures to", normalizePath(OUT), "\n\n")

quiet <- function(expr) {
  lf <- tempfile(); con <- file(lf, open = "wt")
  sink(con); sink(con, type = "message")
  v <- tryCatch(suppressWarnings(expr), error = function(e) e)
  sink(type = "message"); sink(); close(con); unlink(lf)
  v
}
abscos <- function(a, b) abs(sum(a * b)) / sqrt(sum(a^2) * sum(b^2))
# align an estimated direction to the truth's sign (projections are signless)
align <- function(est, truth) if (sum(est * truth) < 0) -est else est

COL_T <- "grey35"     # truth
COL_E <- "#C0392B"    # estimate
png_open <- function(f, w = 1800, h = 700) {
  grDevices::png(file.path(OUT, f), width = w, height = h, res = 150)
}

## ---------------------------------------------------------------- spcma -----
cat("spcma.png ... ")
d   <- spcma_example()
fit <- quiet(spcma(d$X, d$M, d$Y, adaptive = TRUE, var.per = 0.8,
                   boot = FALSE, PC.run = TRUE))
png_open("spcma.png", 1800, 620)
graphics::par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
for (j in 1:2) {
  tr <- d$truth$Phi[, j]
  es <- align(fit$SPCA$W[, j], tr)
  plot(tr, type = "h", lwd = 3, col = COL_T, xlab = "mediator",
       ylab = "loading",
       main = sprintf("Sparse component %d (cos = %.3f)", j, abscos(es, tr)))
  graphics::points(es, pch = 19, cex = 0.7, col = COL_E)
  graphics::abline(h = 0, col = "grey80")
  graphics::legend("topright", c("truth", "estimate"), bty = "n",
                   col = c(COL_T, COL_E), lwd = c(3, NA), pch = c(NA, 19))
}
ie  <- fit$SPCA$IE[, "Estimate"]
tie <- d$truth$IE[seq_along(ie)]
bp <- graphics::barplot(rbind(tie, ie), beside = TRUE, col = c(COL_T, COL_E),
                        names.arg = paste0("C", seq_along(ie)),
                        xlab = "component", ylab = "indirect effect",
                        main = "Indirect effect by component")
graphics::abline(h = 0, col = "grey60")
graphics::legend("topright", c("truth", "estimate"), fill = c(COL_T, COL_E),
                 bty = "n")
grDevices::dev.off()
cat("done\n")

## ----------------------------------------------------------------- cfma -----
cat("cfma.png ... ")
d   <- cfma_example()
fit <- quiet(cfma_concurrent(d$Z, d$M, d$Y, intercept = FALSE,
                             timeinv = d$timeinv))
tg  <- seq(d$timeinv[1], d$timeinv[2], length.out = length(fit$IE$curve))
png_open("cfma.png", 1800, 620)
graphics::par(mfrow = c(1, 3), mar = c(4, 4, 3, 1))
panels <- list(
  list(est = fit$M$curve[1, ], tru = d$truth$alpha, lab = expression(alpha(t)),
       main = "Treatment -> mediator"),
  list(est = fit$Y$curve[2, ], tru = d$truth$beta,  lab = expression(beta(t)),
       main = "Mediator -> outcome"),
  list(est = fit$IE$curve,     tru = d$truth$IE,    lab = "IE(t)",
       main = "Indirect effect")
)
for (pn in panels) {
  yl <- range(c(pn$est, pn$tru))
  plot(tg, pn$tru, type = "l", lwd = 6, col = COL_T, ylim = yl,
       xlab = "time", ylab = pn$lab, main = pn$main)
  graphics::lines(tg, pn$est, lwd = 2, col = COL_E)
  graphics::legend("bottomright", c("truth", "estimate"), bty = "n",
                   col = c(COL_T, COL_E), lwd = c(6, 2))
}
grDevices::dev.off()
cat("done\n")

## ------------------------------------------------------------ pathlasso -----
cat("pathlasso.png ... ")
d   <- pathlasso_example()
fit <- quiet(pathlasso(d$X, d$M, d$Y, lambda = 0.001, omega = 0, phi = 1,
                       max.itr = 3000, tol = 1e-8))
png_open("pathlasso.png", 1250, 620)
graphics::par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
for (w in list(list(e = fit$A, t = d$truth$A, m = "a path (X -> M)"),
               list(e = fit$B, t = d$truth$B, m = "b path (M -> Y)"))) {
  plot(w$t, type = "h", lwd = 3, col = COL_T, xlab = "mediator",
       ylab = "coefficient", main = w$m,
       ylim = range(c(w$e, w$t)))
  graphics::points(w$e, pch = 19, cex = 0.7, col = COL_E)
  graphics::abline(h = 0, col = "grey80")
  graphics::legend("topright", c("truth", "estimate"), bty = "n",
                   col = c(COL_T, COL_E), lwd = c(3, NA), pch = c(NA, 19))
}
grDevices::dev.off()
cat("done\n")

## ----------------------------------------------------------------- gmed -----
cat("gmed.png ... ")
d   <- gmed_example()
fit <- quiet(gmed(d$X, d$M, d$Y, stop.crt = "nD", nD = 1, ninitial = 5,
                  verbose = FALSE))
th  <- align(as.numeric(fit$theta), d$truth$theta)
score <- vapply(d$M, function(Mi) as.numeric(t(th) %*% stats::cov(Mi) %*% th), 0)
png_open("gmed.png", 1250, 620)
graphics::par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
plot(d$truth$theta, type = "h", lwd = 3, col = COL_T, xlab = "node",
     ylab = expression(theta), ylim = range(c(th, d$truth$theta)),
     main = sprintf("Mediating direction (cos = %.3f)",
                    abscos(th, d$truth$theta)))
graphics::points(th, pch = 19, cex = 0.9, col = COL_E)
graphics::abline(h = 0, col = "grey80")
graphics::legend("topright", c("truth", "estimate"), bty = "n",
                 col = c(COL_T, COL_E), lwd = c(3, NA), pch = c(NA, 19))
grp <- as.factor(drop(d$X))
plot(log(score), d$Y, pch = 19, cex = 0.8,
     col = ifelse(grp == "1", COL_E, "#2C6FA6"),
     xlab = expression(log(theta ~ "'" ~ Sigma[i] ~ theta)), ylab = "outcome Y",
     main = "Projected log-variance vs outcome")
graphics::abline(stats::lm(d$Y ~ log(score)), lwd = 2, col = "grey30")
graphics::legend("topleft", c("X = 0", "X = 1"), pch = 19, bty = "n",
                 col = c("#2C6FA6", COL_E))
grDevices::dev.off()
cat("done\n")

## ------------------------------------------------------------- hetermed -----
cat("hetermed.png ... ")
d   <- hetermed_example()
fit <- quiet(hetermed(d$X, d$M, d$Y, d$Z, method = "OLS"))
ite <- quiet(hetermed_ite(d$X, d$Z, fit$alpha0, fit$alpha1, fit$beta0,
                          fit$beta1, fit$gamma0, fit$gamma1))
png_open("hetermed.png", 1250, 620)
graphics::par(mfrow = c(1, 2), mar = c(4, 4, 3, 1))
plot(d$truth$NIE, ite[, "NIE"], pch = 19, cex = 0.6, col = COL_E,
     xlab = "true individual NIE", ylab = "estimated individual NIE",
     main = sprintf("Individual indirect effect (r = %.3f)",
                    stats::cor(ite[, "NIE"], d$truth$NIE)))
graphics::abline(0, 1, lwd = 2, col = COL_T)
graphics::hist(ite[, "NIE"], breaks = 30, col = "grey85", border = "white",
               xlab = "estimated individual NIE",
               main = "Heterogeneity of the indirect effect")
graphics::abline(v = mean(ite[, "NIE"]), lwd = 3, col = COL_E)
graphics::abline(v = mean(d$truth$NIE), lwd = 3, col = COL_T, lty = 2)
graphics::legend("topright", c("mean (est)", "mean (true)"), bty = "n",
                 col = c(COL_E, COL_T), lwd = 3, lty = c(1, 2))
grDevices::dev.off()
cat("done\n")

cat("\nFigures written:\n")
print(list.files(OUT, pattern = "[.]png$"))
