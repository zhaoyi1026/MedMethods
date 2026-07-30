##############################################################################
# Built-in example-data generators, one per method.
#
# Each returns a named list of inputs shaped exactly as the matching wrapper
# expects, plus a `truth` list of the data-generating parameters so estimates
# can be checked. Self-contained: base R only (no mvtnorm dependency).
##############################################################################

# Multivariate normal draws via an eigen root. Fills the normal matrix BY ROW so
# that a given seed reproduces mvtnorm::rmvnorm().
.med_rmvnorm <- function(n, mean, sigma) {
  sigma <- as.matrix(sigma)
  ev <- eigen(sigma, symmetric = TRUE)
  R  <- ev$vectors %*% (t(ev$vectors) * sqrt(pmax(ev$values, 0)))
  Z  <- matrix(stats::rnorm(n * ncol(sigma)), nrow = n, byrow = TRUE)
  sweep(Z %*% R, 2L, mean, "+")
}

# Orthonormal basis from a seed matrix, sign-normalised so the largest-magnitude
# entry of each column is positive (makes a fixed seed give a stable basis).
.med_orth <- function(seed_mat) {
  Q <- qr.Q(qr(seed_mat))
  for (j in seq_len(ncol(Q))) {
    if (Q[which.max(abs(Q[, j])), j] < 0) Q[, j] <- -Q[, j]
  }
  Q
}

# Cosine similarity up to sign (projections are identified only up to sign).
.med_abscos <- function(a, b) {
  a <- as.numeric(a); b <- as.numeric(b)
  abs(sum(a * b)) / sqrt(sum(a^2) * sum(b^2))
}

#' Example data for the mediation methods
#'
#' Self-contained synthetic-data generators, one per method, each returning the
#' inputs in exactly the shape the matching wrapper expects plus a `truth`
#' element holding the data-generating parameters. A full run is therefore two
#' lines, e.g. `d <- gmed_example(); gmed(d$X, d$M, d$Y, nD = 1)`.
#'
#' @param n number of subjects (or time points, for `gma_example`).
#' @param N number of subjects for the functional model (`cfma_example`), or
#'   number of series (`gma_example`, two-level).
#' @param n.time expected time points per series (`gma_example`, two-level).
#' @param K number of sessions per subject (`macc_example`, three-level data).
#' @param p mediator dimension. For `pcma_example` the number of exposures; for
#'   `hdmediation_example` the number of mediators.
#' @param q number of mediators (`pcma_example`).
#' @param r number of exposures (`hdmediation_example`).
#' @param p1,p2 sizes of the two mediator blocks (`pathlasso2b_example`).
#' @param k number of mediators (`pathlasso_example`).
#' @param Ti within-subject sample size, i.e. rows per mediator matrix
#'   (`gmed_example`).
#' @param ntp number of time points (`cfma_example`).
#' @param n.trial trials per subject (`macc_example`, multilevel data).
#' @param model.type `"single"`, `"twolevel"` or `"three"` (`macc_example`);
#'   `"single"` or `"twolevel"` (`gma_example`). Only the multilevel forms make
#'   the error correlation `delta` identifiable.
#' @param delta error correlation between the mediator and outcome models
#'   (`macc_example`, `gma_example`).
#' @param seed RNG seed.
#' @return A named list of inputs plus a `truth` list of data-generating
#'   parameters.
#' @examples
#' d <- hetermed_example()
#' str(d, max.level = 1)
#' @name med_examples
NULL


#' @rdname med_examples
#' @export
macc_example <- function(model.type = c("single", "twolevel", "three"),
                         n = 100L, N = 50L, K = 4L, n.trial = 100L,
                         delta = 0.5, seed = 5000L) {
  model.type <- match.arg(model.type)
  A <- 0.5; B <- -1; C <- 0.5
  # Theta[1,1] = A, Theta[1,2] = C, Theta[2,2] = B
  Theta <- matrix(c(A, 0, C, B), 2, 2)
  Sigma <- matrix(c(1, 2 * delta, 2 * delta, 4), 2, 2)
  truth <- list(A = A, B = B, C = C, C2 = C + A * B,
                ABp = A * B, ABd = A * B, delta = delta)

  if (model.type == "single") {
    set.seed(100L)
    Z <- matrix(stats::rbinom(n, size = 1L, prob = 0.5), n, 1L)
    set.seed(seed)
    dat <- .med_get("macc", "sim.data.single")(Z, Theta, Sigma)
    return(list(dat = dat, model.type = "single", truth = truth))
  }

  # NB: sim.data.multi() wraps the trial-level data in a `data` element, while
  # macc() expects that element itself.
  Lambda <- diag(rep(0.5^2, 3))
  Psi    <- diag(rep(0.5^2, 3))
  set.seed(100L)
  # trials per subject vary, as in the published multilevel example data
  ntr <- round(stats::runif(N, 0.75 * n.trial, 1.25 * n.trial))

  if (model.type == "twolevel") {
    Z.list <- lapply(seq_len(N), function(i)
      stats::rbinom(ntr[i], size = 1L, prob = 0.5))
    set.seed(seed)
    sim <- .med_get("macc", "sim.data.multi")(Z.list, N = N, K = 1L,
                                              Theta = Theta, Sigma = Sigma,
                                              Lambda = Lambda)
    return(list(dat = sim$data, model.type = "twolevel",
                truth = c(truth, list(A.subj = sim$A, B.subj = sim$B,
                                      C.subj = sim$C))))
  }

  Z.list <- lapply(seq_len(N), function(i)
    lapply(seq_len(K), function(k) stats::rbinom(ntr[i], size = 1L, prob = 0.5)))
  set.seed(seed)
  sim <- .med_get("macc", "sim.data.multi")(Z.list, N = N, K = K,
                                            Theta = Theta, Sigma = Sigma,
                                            Psi = Psi, Lambda = Lambda)
  list(dat = sim$data, model.type = "multilevel",
       truth = c(truth, list(A.subj = sim$A, B.subj = sim$B, C.subj = sim$C)))
}


#' @rdname med_examples
#' @export
gma_example <- function(model.type = c("single", "twolevel"), n = 500L,
                        N = 40L, n.time = 150L, delta = 0.5, seed = 1000L) {
  model.type <- match.arg(model.type)
  A <- 0.5; B <- -1; C <- 0.5
  Sigma <- matrix(c(1, 2 * delta, 2 * delta, 4), 2, 2)
  # VAR(1) initial-condition covariance and transition matrix (2p x 2, p = 1)
  Delta <- matrix(c(2, delta * sqrt(2 * 8), delta * sqrt(2 * 8), 8), 2, 2)
  W     <- matrix(c(-0.809, 0.154, -0.618, -0.5), 2, 2)
  truth <- list(A = A, B = B, C = C, C2 = C + A * B, ABp = A * B,
                delta = delta, p = 1L, W = W)

  if (model.type == "single") {
    set.seed(1000L)
    Z <- matrix(stats::rbinom(n, size = 1L, prob = 0.5), n, 1L)
    set.seed(seed)
    # NB: sim.data.ts.single() returns list(data, error); gma() wants `data`.
    sim <- .med_get("gma", "sim.data.ts.single")(n, Z, A, B, C, Sigma, W,
                                                 Delta = Delta, p = 1L,
                                                 nburn = 1000L)
    truth$error <- sim$error
    return(list(dat = sim$data, model.type = "single", truth = truth))
  }

  # Two-level: N series, each its own VAR(1) process, with the subject-level
  # coefficients drawn around (A, B, C) with covariance Lambda. Only here is the
  # error correlation `delta` identifiable, because it is separated from the
  # between-series variability.
  Lambda <- diag(0.5, 3)
  set.seed(2000L)
  ni <- matrix(stats::rpois(N, n.time), N, 1)
  set.seed(1000L)
  Z.list <- lapply(seq_len(N), function(i)
    matrix(stats::rbinom(ni[i, 1], size = 1L, prob = 0.5), ni[i, 1], 1L))
  set.seed(seed)
  sim <- .med_get("gma", "sim.data.ts.two")(Z.list, N, theta = c(A, B, C),
                                            Sigma, W, Delta = Delta, p = 1L,
                                            Lambda = Lambda, nburn = 500L)
  truth$Lambda <- Lambda
  truth$A.series <- sim$A; truth$B.series <- sim$B; truth$C.series <- sim$C
  list(dat = sim$data, model.type = "twolevel", truth = truth)
}


#' @rdname med_examples
#' @export
spcma_example <- function(n = 200L, p = 50L, seed = 2026L) {
  # Piecewise-constant leading loadings -- what the fused-lasso penalty in
  # spcma() is designed to recover: three blocks of consecutive mediators.
  # The block width adapts to p so that the three blocks always fit (at p = 50,
  # the default, this gives the intended width of 10).
  if (p < 6L) stop("spcma_example(): p must be at least 6.")
  nb  <- max(1L, p %/% 5L)
  if (3L * nb > p) nb <- p %/% 3L
  Phi0 <- matrix(0, p, p)
  Phi0[1:nb, 1]                  <- 1
  Phi0[(nb + 1):(2 * nb), 2]     <- 1
  Phi0[(2 * nb + 1):(3 * nb), 3] <- 1
  set.seed(seed)
  Phi0[, 4:p] <- matrix(stats::rnorm(p * (p - 3)), p, p - 3)
  Phi <- .med_orth(Phi0)

  alpha <- c(2, -1.5, 0, rep(0, p - 3))   # exposure -> latent component
  beta  <- c(2,  1,   0, rep(0, p - 3))   # latent component -> outcome
  gamma <- 1                              # direct effect

  # The three structured components must carry the LARGEST variance, otherwise
  # the (sparse) PCA step has no reason to find them before the noise directions.
  sds <- c(4, 3, 2, rep(0.3, p - 3))

  set.seed(seed)
  X    <- matrix(stats::rbinom(n, size = 1L, prob = 0.5), n, 1L)
  Mlat <- X %*% matrix(alpha, nrow = 1L) +
    sweep(matrix(stats::rnorm(n * p), n, p), 2L, sds, "*")  # component scores
  M    <- Mlat %*% t(Phi)                                   # so M %*% Phi == Mlat
  Y    <- matrix(gamma * X + Mlat %*% beta + stats::rnorm(n, sd = 1), n, 1L)

  list(X = X, M = M, Y = Y,
       truth = list(Phi = Phi, alpha = alpha, beta = beta, gamma = gamma,
                    IE = alpha * beta, n.block = nb, comp.sd = sds))
}


#' @rdname med_examples
#' @export
cfma_example <- function(N = 200L, ntp = 150L, seed = 2026L) {
  tg <- seq(0, 1, length.out = ntp)

  # Time-varying coefficient curves, all clearly non-constant
  alpha <- 1 + 0.5 * sin(2 * pi * tg)
  beta  <- 1 + 0.5 * cos(2 * pi * tg)
  gamma <- 0.8 * sin(pi * tg)

  # Treatment trajectory: boxcar event onsets convolved with a single-gamma HRF
  hrf <- function(u) ifelse(u < 0, 0, u^5 * exp(-u) / gamma(6))
  h   <- hrf(seq(0, 20, length.out = 21L))
  h   <- h / max(h)
  set.seed(seed)
  Z <- t(vapply(seq_len(N), function(i) {
    onset <- rep(0, ntp)
    starts <- seq(5L, ntp - 10L, by = 25L)
    for (s in starts) onset[s:(s + 4L)] <- 1
    z <- stats::filter(onset, h, sides = 1L)
    z <- as.numeric(z); z[is.na(z)] <- 0
    z + stats::rnorm(ntp, sd = 0.05)
  }, numeric(ntp)))

  # The mediator error must be substantial: with M nearly a deterministic
  # multiple of Z the two become collinear and the direct effect gamma(t) is
  # only weakly identified.
  M <- Z * matrix(alpha, N, ntp, byrow = TRUE) +
    matrix(stats::rnorm(N * ntp, sd = 1), N, ntp)
  Y <- Z * matrix(gamma, N, ntp, byrow = TRUE) +
    M * matrix(beta, N, ntp, byrow = TRUE) +
    matrix(stats::rnorm(N * ntp, sd = 0.3), N, ntp)

  list(Z = Z, M = M, Y = Y, timeinv = c(0, 1),
       truth = list(alpha = alpha, beta = beta, gamma = gamma,
                    IE = alpha * beta, DE = gamma, timegrids = tg))
}


#' @rdname med_examples
#' @export
pathlasso_example <- function(n = 200L, k = 50L, seed = 100L) {
  # Four signal mediators with a-path / b-path products of mixed sign; the
  # remaining k - 4 mediators are pure noise.
  #   M = X A + e1,  Y = X C + M B + e2
  amp <- 2
  A <- matrix(c(c(amp, -amp, amp, -amp), rep(0, k - 4L)), nrow = 1L)
  B <- matrix(c(c(amp,  amp, -amp, -amp), rep(0, k - 4L)), ncol = 1L)
  AB <- drop(t(A) * B)
  C  <- 1

  set.seed(seed)
  X <- matrix(stats::rbinom(n, size = 1L, prob = 0.5), ncol = 1L)
  M <- X %*% A + matrix(stats::rnorm(n * k, sd = 1), n, k)
  Y <- matrix(X * C + M %*% B + stats::rnorm(n, sd = 1), n, 1L)

  list(X = X, M = M, Y = Y,
       truth = list(A = drop(A), B = drop(B), C = C, AB = AB,
                    signal = which(abs(AB) > 0)))
}


#' @rdname med_examples
#' @export
pathlasso2b_example <- function(n = 200L, p1 = 20L, p2 = 30L, seed = 100L) {
  # X -> M1 -> M2 -> Y with two blocks of mediators:
  #   M1 = X beta  + e1
  #   M2 = X zeta  + M1 Lambda + e2
  #   Y  = X delta + M1 theta  + M2 pi + eY
  nz1 <- 4L; nz2 <- 4L
  beta  <- matrix(c(rep(c(1, -1), length.out = nz1), rep(0, p1 - nz1)), nrow = 1L)
  theta <- c(rep(c(1.5, 1.5), length.out = nz1), rep(0, p1 - nz1))
  zeta  <- matrix(c(rep(c(1, -1), length.out = nz2), rep(0, p2 - nz2)), nrow = 1L)
  pi    <- c(rep(c(1.5, -1.5), length.out = nz2), rep(0, p2 - nz2))
  Lambda <- matrix(0, p1, p2)
  Lambda[cbind(1:nz1, 1:nz1)] <- 1
  delta <- 1

  set.seed(seed)
  X  <- matrix(stats::rbinom(n, size = 1L, prob = 0.5), n, 1L)
  M1 <- X %*% beta + matrix(stats::rnorm(n * p1, sd = 1), n, p1)
  M2 <- X %*% zeta + M1 %*% Lambda + matrix(stats::rnorm(n * p2, sd = 1), n, p2)
  Y  <- matrix(X * delta + M1 %*% theta + M2 %*% pi + stats::rnorm(n, sd = 1),
               n, 1L)

  list(X = drop(X), M1 = M1, M2 = M2, Y = Y,
       truth = list(beta = drop(beta), theta = theta, zeta = drop(zeta),
                    pi = pi, Lambda = Lambda, delta = delta,
                    IE.M1 = drop(beta) * theta,
                    IE.M2 = drop(zeta) * pi,
                    IE.M1M2 = outer(drop(beta), pi) * Lambda))
}


#' @rdname med_examples
#' @export
hdmediation_example <- function(n = 100L, r = 20L, p = 20L, seed = 100L) {
  # High-dimensional exposures X (n x r) and mediators M (n x p):
  #   M = X alpha + eM,   Y = X gamma + M beta + eY
  # with sparse alpha / beta so only a few exposure-mediator paths carry signal.
  nz <- 3L
  alpha <- matrix(0, r, p)
  alpha[cbind(1:nz, 1:nz)] <- c(1.5, -1.5, 1.5)
  beta  <- c(c(1.5, 1.5, -1.5), rep(0, p - nz))
  gamma <- c(1, rep(0, r - 1L))

  set.seed(seed)
  X <- matrix(stats::rnorm(n * r), n, r)
  M <- X %*% alpha + matrix(stats::rnorm(n * p, sd = 0.5), n, p)
  Y <- matrix(X %*% gamma + M %*% beta + stats::rnorm(n, sd = 0.5), n, 1L)

  IE <- alpha %*% diag(beta)
  list(X = X, M = M, Y = Y,
       truth = list(alpha = alpha, beta = beta, gamma = gamma, IE = IE,
                    signal = which(abs(IE) > 0, arr.ind = TRUE)))
}


#' @rdname med_examples
#' @export
pcma_example <- function(n = 400L, p = 5L, q = 10L, seed = 100L) {
  # Orthogonal exposure (Phi) and mediator (Psi) projections with parallel
  # mediation mechanisms on the projected scores.
  set.seed(50L)
  Phi <- .med_orth(matrix(stats::rnorm(p * p), p, p))
  set.seed(100L)
  Psi <- .med_orth(matrix(stats::rnorm(q * q), q, q))

  sigma2.X <- sort(exp(seq(3, -4, length.out = p)), decreasing = TRUE)
  Sigma.X  <- Phi %*% diag(sigma2.X) %*% t(Phi)

  alpha0 <- c(2, 2)
  alpha  <- matrix(0, p, q)
  alpha[seq_along(alpha0), seq_along(alpha0)] <- diag(alpha0)
  beta   <- c(c(2, 1), rep(0, q - 2L))
  gamma  <- c(c(1, -1), rep(0, p - 2L))
  sigma.M <- rep(1, q); sigma.Y <- 1

  set.seed(seed)
  X   <- .med_rmvnorm(n, rep(0, p), Sigma.X)
  X.t <- X %*% Phi                                  # exposure component scores
  e.M <- .med_rmvnorm(n, rep(0, q), diag(sigma.M^2))
  M.t <- matrix(0, n, q)                            # mediator component scores
  kk  <- min(p, q)
  for (k in seq_len(kk)) M.t[, k] <- X.t[, k] * alpha[k, k] + e.M[, k]
  if (q > p) M.t[, (p + 1L):q] <- e.M[, (p + 1L):q]
  M <- M.t %*% t(Psi)
  Y <- matrix(X.t %*% gamma + M.t %*% beta + stats::rnorm(n, sd = sigma.Y),
              n, 1L)

  list(X = X, M = M, Y = Y,
       truth = list(Phi = Phi, Psi = Psi, alpha = alpha, beta = beta,
                    gamma = gamma, IE = alpha %*% diag(beta),
                    nD = length(alpha0)))
}


#' @rdname med_examples
#' @export
gmed_example <- function(n = 100L, p = 10L, Ti = 500L, seed = 100L) {
  # Covariance (graph) mediator: the exposure shifts the log-eigenvalue of the
  # subject covariance along ONE direction of a common eigenbasis, and the
  # outcome depends on that log-eigenvalue.
  #   log(theta' Sigma_i theta) = alpha0 + X alpha + eta_i
  #   Y = gamma0 + X gamma + beta * log(theta' Sigma_i theta) + eps
  set.seed(100L)
  Gamma <- .med_orth(matrix(stats::runif(p), nrow = p, ncol = p))

  sig.dim <- 2L                       # the single mediating direction
  alpha0  <- 1
  alpha   <- 1                        # exposure -> log-variance (a-path)
  beta    <- 1                        # log-variance -> outcome (b-path)
  gamma0  <- 1
  gamma   <- 1                        # direct effect
  base    <- seq(1.5, -3, length.out = p)   # covariate-free log-variances
  tau <- 0.1; sigma <- 0.1

  set.seed(seed)
  X <- matrix(stats::rbinom(n, size = 1L, prob = 0.5), n, 1L)
  colnames(X) <- "X"

  d <- matrix(NA_real_, n, p)         # subject eigenvalues
  for (j in seq_len(p)) {
    d[, j] <- if (j == sig.dim) {
      exp(alpha0 + X * alpha + stats::rnorm(n, mean = 0, sd = tau))
    } else {
      exp(base[j])
    }
  }
  M <- vector("list", n)
  for (i in seq_len(n)) {
    M[[i]] <- .med_rmvnorm(Ti, rep(0, p), Gamma %*% diag(d[i, ]) %*% t(Gamma))
  }
  set.seed(seed)
  Y <- matrix(gamma0 + X * gamma + log(d[, sig.dim]) * beta +
                stats::rnorm(n, mean = 0, sd = sigma), n, 1L)

  list(X = X, M = M, Y = Y,
       truth = list(Gamma = Gamma, theta = Gamma[, sig.dim], dim = sig.dim,
                    alpha0 = alpha0, alpha = alpha, beta = beta,
                    gamma = gamma, gamma0 = gamma0, IE = alpha * beta))
}


#' @rdname med_examples
#' @export
hetermed_example <- function(n = 600L, seed = 100L) {
  # Treatment coded +/-1; the a-path and b-path are both moderated by Z, so the
  # natural indirect effect varies across subjects.
  alpha0 <- c(0.2, 0.3, 0.0)
  alpha1 <- c(0.5, 0.4, 0.0)
  gamma0 <- c(0.1, 0.0, 0.0)
  gamma1 <- c(0.3, 0.0, 0.0)
  beta0  <- 0.8
  beta1  <- 0.2

  set.seed(seed)
  X <- sample(c(-1, 1), n, replace = TRUE)
  Z <- cbind(Intercept = 1,
             Z1 = stats::rnorm(n),
             Z2 = stats::rbinom(n, size = 1L, prob = 0.5))
  M <- as.vector(Z %*% alpha0 + X * (Z %*% alpha1)) + stats::rnorm(n, sd = 0.5)
  Y <- as.vector(Z %*% gamma0 + X * (Z %*% gamma1)) +
    (beta0 + beta1 * X) * M + stats::rnorm(n, sd = 0.5)

  list(X = X, M = M, Y = Y, Z = Z,
       truth = list(alpha0 = alpha0, alpha1 = alpha1,
                    beta0 = beta0, beta1 = beta1,
                    gamma0 = gamma0, gamma1 = gamma1,
                    NIE = 2 * (beta0 + beta1 * X) * as.vector(Z %*% alpha1),
                    NDE = 2 * as.vector(Z %*% gamma1) +
                      2 * beta1 * as.vector(Z %*% alpha0 +
                                              X * (Z %*% alpha1))))
}
