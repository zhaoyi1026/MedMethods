"""Heterogeneous causal mediation effects.

The treatment ``X`` is coded +1 / -1 and both the a-path and the b-path are
moderated by covariates ``Z``, so every subject has their own natural indirect
and direct effect:

    M model:  M_i = Z_i'(alpha0 + X_i alpha1) + e_i
    Y model:  Y_i = Z_i'(gamma0 + X_i gamma1) + (beta0 + X_i beta1) M_i + u_i

    NIE_i = 2 (beta0 + beta1 X_i) (Z_i' alpha1)
    NDE_i = 2 Z_i' gamma1 + 2 beta1 (Z_i' alpha0 + X_i Z_i' alpha1)

``med_inter`` fits the two models by stacking the arm-specific regressions;
``med_inter_inf`` adds asymptotic (sandwich-free, model-based) inference for the
coefficients and delta-method standard errors for the per-subject effects.

Only the OLS estimator is ported. The generalized-lasso variant needs an
equivalent of R's ``genlasso`` (a generalized-lasso path solver), which has no
drop-in Python counterpart; use ``MedMethods::hetermed(method = "genlasso")``.

Standard-error note
-------------------
The R implementation builds its coefficient tables with ``SE = sqrt(Theta.se)``
where ``Theta.se`` already holds standard errors, i.e. it takes the square root
twice. This port reports the correct standard errors. Confirmed empirically:
the R value scales as ``n^(-1/4)`` rather than ``n^(-1/2)``, and its *square*
matches the Monte Carlo SD of the estimator. (R's NIE/NDE tables use the
covariance matrix directly and were always correct.)
"""
from __future__ import annotations

import numpy as np
from scipy.stats import norm

__all__ = ["med_inter", "med_inter_ite", "med_inter_inf"]


def _prep(X, M, Y, Z):
    X = np.asarray(X, dtype=float).ravel()
    M = np.asarray(M, dtype=float).ravel()
    Y = np.asarray(Y, dtype=float).ravel()
    Z = np.asarray(Z, dtype=float)
    if Z.ndim == 1:
        Z = Z[:, None]
    return X, M, Y, Z


def med_inter_ite(X, Z, alpha0, alpha1, beta0, beta1, gamma0, gamma1):
    """Per-subject natural indirect / direct effects (== R ``med.inter.ITE``)."""
    X = np.asarray(X, dtype=float).ravel()
    Z = np.asarray(Z, dtype=float)
    if Z.ndim == 1:
        Z = Z[:, None]
    alpha0 = np.asarray(alpha0, dtype=float).ravel()
    alpha1 = np.asarray(alpha1, dtype=float).ravel()
    gamma1 = np.asarray(gamma1, dtype=float).ravel()
    Za1 = Z @ alpha1
    NIE = 2.0 * (beta0 + beta1 * X) * Za1
    NDE = 2.0 * (Z @ gamma1) + 2.0 * beta1 * (Z @ alpha0 + X * Za1)
    return {"treatment": X, "NIE": NIE, "NDE": NDE}


def _stacked_ols(y, D, idx1, idx2, q):
    """Fit arm-specific coefficients by stacking the two arms, then split.

    Returns ``(phi0, phi1)`` = the average and half-difference of the two arms'
    coefficient vectors, which are the ``*0`` (main) and ``*1`` (moderated)
    parameters of the model.
    """
    yy = np.concatenate([y[idx1], y[idx2]])
    top = np.hstack([D[idx1], np.zeros((idx1.size, q))])
    bot = np.hstack([np.zeros((idx2.size, q)), D[idx2]])
    XX = np.vstack([top, bot])
    th = np.linalg.solve(XX.T @ XX, XX.T @ yy)
    th1, th2 = th[:q], th[q:]
    return (th1 + th2) / 2.0, (th1 - th2) / 2.0


def med_inter(X, M, Y, Z, method="OLS"):
    """Estimate the moderated mediation model (== R ``med.inter``).

    Parameters
    ----------
    X : (n,) treatment coded +1 / -1.
    M : (n,) mediator.
    Y : (n,) outcome.
    Z : (n, p) covariates whose FIRST column is a column of ones.
    method : only ``"OLS"`` is implemented.
    """
    if method != "OLS":
        raise NotImplementedError(
            "Only method='OLS' is ported. The generalized-lasso variant needs an "
            "equivalent of R's genlasso path solver; use "
            "MedMethods::hetermed(method = 'genlasso') in R."
        )
    X, M, Y, Z = _prep(X, M, Y, Z)
    n, p = Z.shape
    idx1 = np.flatnonzero(X == 1)
    idx2 = np.flatnonzero(X == -1)
    if idx1.size == 0 or idx2.size == 0:
        raise ValueError("X must be coded +1 / -1 with both arms present")

    alpha0, alpha1 = _stacked_ols(M, Z, idx1, idx2, p)

    W = np.hstack([Z, M[:, None]])
    phi0, phi1 = _stacked_ols(Y, W, idx1, idx2, p + 1)
    gamma0, beta0 = phi0[:p], phi0[p]
    gamma1, beta1 = phi1[:p], phi1[p]

    ite = med_inter_ite(X, Z, alpha0, alpha1, beta0, beta1, gamma0, gamma1)
    return {
        "ITE": ite, "alpha0": alpha0, "alpha1": alpha1,
        "beta0": float(beta0), "beta1": float(beta1),
        "gamma0": gamma0, "gamma1": gamma1,
    }


def med_inter_inf(X, M, Y, Z, fit, conf_level=0.95):
    """Asymptotic inference for :func:`med_inter` (== R ``fit.inf.OLS``).

    Returns per-coefficient tables (estimate, SE, z, p, CI) and per-subject
    NIE/NDE tables with delta-method standard errors.
    """
    X, M, Y, Z = _prep(X, M, Y, Z)
    n, p = Z.shape
    a0 = np.asarray(fit["alpha0"], dtype=float).ravel()
    a1 = np.asarray(fit["alpha1"], dtype=float).ravel()
    g0 = np.asarray(fit["gamma0"], dtype=float).ravel()
    g1 = np.asarray(fit["gamma1"], dtype=float).ravel()
    b0 = float(fit["beta0"])
    b1 = float(fit["beta1"])

    U = X[:, None] * Z          # diag(X) %*% Z
    V = X * M                   # diag(X) %*% M

    s2m = np.mean((M - Z @ a0 - U @ a1) ** 2)
    s2y = np.mean((Y - Z @ g0 - U @ g1 - M * b0 - V * b1) ** 2)

    Q = (Z.T @ Z) / n
    mm = np.array([np.mean(X), np.mean(X ** 2), np.mean(X ** 3), np.mean(X ** 4)])
    Qa0, Qa1 = Q @ a0, Q @ a1

    Qzm = Qa0 + mm[0] * Qa1
    Qzv = mm[0] * Qa0 + mm[1] * Qa1
    Qum = mm[0] * Qa0 + mm[1] * Qa1
    Quv = mm[1] * Qa0 + mm[2] * Qa1
    a0Qa0 = a0 @ Q @ a0
    a0Qa1 = a0 @ Q @ a1
    a1Qa1 = a1 @ Q @ a1
    Qm = a0Qa0 + 2 * mm[0] * a0Qa1 + mm[1] * a1Qa1 + s2m
    Qmv = mm[0] * a0Qa0 + 2 * mm[1] * a0Qa1 + mm[2] * a1Qa1 + mm[0] * s2m
    Qv = mm[1] * a0Qa0 + 2 * mm[2] * a0Qa1 + mm[3] * a1Qa1 + mm[1] * s2m

    # assembled explicitly: np.block() mis-parses the scalar corner cells
    d = 2 * p + 2
    Qx = np.zeros((d, d))
    Qx[:p, :p] = Q
    Qx[:p, p:2 * p] = mm[0] * Q
    Qx[p:2 * p, :p] = mm[0] * Q
    Qx[p:2 * p, p:2 * p] = mm[1] * Q
    Qx[:p, 2 * p] = Qzm
    Qx[2 * p, :p] = Qzm
    Qx[:p, 2 * p + 1] = Qzv
    Qx[2 * p + 1, :p] = Qzv
    Qx[p:2 * p, 2 * p] = Qum
    Qx[2 * p, p:2 * p] = Qum
    Qx[p:2 * p, 2 * p + 1] = Quv
    Qx[2 * p + 1, p:2 * p] = Quv
    Qx[2 * p, 2 * p] = Qm
    Qx[2 * p, 2 * p + 1] = Qmv
    Qx[2 * p + 1, 2 * p] = Qmv
    Qx[2 * p + 1, 2 * p + 1] = Qv
    Qx_inv = np.linalg.pinv(Qx)
    Sigma = np.diag([s2m, s2y])
    cov_vec = np.kron(Sigma, Qx_inv) / n

    # NB: this is already the standard error. The R code takes sqrt() of it
    # again when building the tables; see the module docstring.
    se_mat = np.sqrt(np.diag(cov_vec)).reshape(2, -1).T   # (2p+2) x 2, column-major

    zv = norm.ppf(1 - (1 - conf_level) / 2)

    def _tab(est, se, names):
        est = np.atleast_1d(np.asarray(est, dtype=float))
        se = np.atleast_1d(np.asarray(se, dtype=float))
        z = est / se
        return {
            "names": list(names), "estimate": est, "SE": se, "zvalue": z,
            "pvalue": 2 * (1 - norm.cdf(np.abs(z))),
            "LB": est - zv * se, "UB": est + zv * se,
        }

    znames = ["Intercept"] + ["Z%d" % (j + 1) for j in range(p - 1)]
    out = {
        "alpha0": _tab(a0, se_mat[:p, 0], znames),
        "alpha1": _tab(a1, se_mat[p:2 * p, 0], znames),
        "gamma0": _tab(g0, se_mat[:p, 1], znames),
        "gamma1": _tab(g1, se_mat[p:2 * p, 1], znames),
        "beta0": _tab([b0], se_mat[2 * p, 1], ["M"]),
        "beta1": _tab([b1], se_mat[2 * p + 1, 1], ["M"]),
    }

    # delta-method variances for the per-subject effects
    zeros_p = np.zeros(p)
    ie_var = np.empty(n)
    de_var = np.empty(n)
    for i in range(n):
        z_i = Z[i]
        x_i = X[i]
        a1z = a1 @ z_i
        h1 = np.concatenate([
            zeros_p, 2 * (b0 + b1 * x_i) * z_i, np.zeros(2 + 2 * p),
            [2 * a1z], [2 * x_i * a1z],
        ])
        h2 = np.concatenate([
            2 * b1 * z_i, 2 * b1 * x_i * z_i, np.zeros(2 + p), 2 * z_i, [0.0],
            [2 * (a0 @ z_i + x_i * a1z)],
        ])
        ie_var[i] = h1 @ cov_vec @ h1
        de_var[i] = h2 @ cov_vec @ h2

    ite = fit["ITE"]
    ie_se = np.sqrt(ie_var)
    de_se = np.sqrt(de_var)
    out["NIE"] = {"treatment": X, "estimate": ite["NIE"], "SE": ie_se,
                  "LB": ite["NIE"] - zv * ie_se, "UB": ite["NIE"] + zv * ie_se}
    out["NDE"] = {"treatment": X, "estimate": ite["NDE"], "SE": de_se,
                  "LB": ite["NDE"] - zv * de_se, "UB": ite["NDE"] + zv * de_se}
    out["vecTheta_cov"] = cov_vec
    out["sigma2"] = {"M": float(s2m), "Y": float(s2y)}
    return out
