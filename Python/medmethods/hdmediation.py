"""Mediation analysis with high-dimensional exposures and high-dimensional mediators.

    M = X alpha + E,        Y = X gamma + M beta + e

with a pathway-lasso style penalty on the exposure-mediator effects
``mu = alpha diag(beta)``, combining an l1 term, a sparse-group term across each
exposure's row of ``mu``, and ridge terms on ``alpha`` and ``beta``. Solved by
ADMM on the splitting ``mu = alpha diag(beta)``.

Tuning: ``lambda`` (overall strength), ``pi`` (l1 vs group split), ``phi``
(ridge), ``delta`` (extra l1 on alpha / beta).

Ported from the R implementation behind ``MedMethods::hdmediation()``.
"""
from __future__ import annotations

import warnings

import numpy as np

from ._core import soft_threshold

__all__ = ["hdmediation", "hdmediation_std", "hdmediation_pca", "hd_obj"]


def _sg_lasso_iden(z, lam1, lam2):
    """Sparse-group-lasso proximal step (== R ``sg.lasso.iden``)."""
    s = soft_threshold(z, lam2)
    nrm = np.sqrt(np.sum(s ** 2))
    if nrm == 0:
        return np.zeros_like(z)
    return max(nrm - lam1, 0.0) * s / nrm


def hd_obj(X, M, Y, alpha, beta, gamma, lam=0.0, pi=1.0, phi=2.0, delta=0.1):
    """Objective value with its penalty decomposition (== R ``obj.func``)."""
    X = np.asarray(X, dtype=float)
    M = np.asarray(M, dtype=float)
    Y = np.asarray(Y, dtype=float).ravel()
    alpha = np.atleast_2d(np.asarray(alpha, dtype=float))
    beta = np.asarray(beta, dtype=float).ravel()
    gamma = np.asarray(gamma, dtype=float).ravel()
    q, p = alpha.shape

    R = M - X @ alpha
    ll = (np.sum(R * R) + np.sum((Y - X @ gamma - M @ beta) ** 2)) / 2.0

    mu = alpha * beta[None, :]
    # R1: pathway-lasso penalty
    r1 = (np.sum(np.abs(mu)) + phi * np.sum(alpha ** 2) + phi * q * np.sum(beta ** 2)
          + delta * np.sum(np.abs(alpha)) + delta * np.sum(np.abs(beta))
          + np.sum(np.abs(gamma)))
    # R2: group penalty across each exposure's row of mu
    r2 = np.sqrt(p) * np.sum(np.sqrt(np.sum(mu ** 2, axis=1))) + np.sum(np.abs(gamma))

    pen1 = pi * lam * r1
    pen2 = (1.0 - pi) * lam * r2
    return {"logLik": float(ll), "R1": float(pen1), "R2": float(pen2),
            "obj": float(ll + pen1 + pen2)}


def hdmediation_std(X, M, Y, lam=0.0, pi=1.0, phi=2.0, delta=0.1, rho=1.0,
                    max_itr=10000, tol=1e-6, rho_increase=False,
                    alpha0=None, beta0=None, gamma0=None):
    """ADMM on already-standardized data (== R ``HDCauseMediation.std``)."""
    X = np.asarray(X, dtype=float)
    M = np.asarray(M, dtype=float)
    Y = np.asarray(Y, dtype=float).ravel()
    n, q = X.shape
    p = M.shape[1]

    rho0 = rho if rho_increase else 0.0
    alpha = np.zeros((q, p)) if alpha0 is None else \
        np.array(alpha0, dtype=float).reshape(q, p)
    beta = np.zeros(p) if beta0 is None else np.asarray(beta0, dtype=float).ravel()
    gamma = np.zeros(q) if gamma0 is None else np.asarray(gamma0, dtype=float).ravel()
    tau = np.zeros((q, p))

    MtM = M.T @ M
    XtX = X.T @ X
    xj_sq = np.einsum("ij,ij->j", X, X)      # diag(X'X)
    Ip = np.eye(p)

    s, diff = 0, 100.0
    mu = np.zeros((q, p))
    while s <= max_itr and diff > tol:
        s += 1

        # --- mu: sparse-group-lasso prox on alpha diag(beta) - tau/rho ---
        z = alpha * beta[None, :] - tau / rho
        lam1 = (1.0 - pi) * lam * np.sqrt(p) / rho
        lam2 = pi * lam / rho
        mu = np.vstack([_sg_lasso_iden(z[j], lam1, lam2) for j in range(q)])

        # --- alpha, row by row.  V is diagonal, so solve() is elementwise ---
        alpha_new = np.empty((q, p))
        for j in range(q):
            if q > 1:
                others = np.delete(np.arange(q), j)
                resid = M - X[:, others] @ alpha[others, :]
            else:
                resid = M
            w = resid.T @ X[:, j] + beta * tau[j] + rho * beta * mu[j]
            dV = rho * beta ** 2 + (2 * pi * lam * phi + xj_sq[j])
            alpha_new[j] = soft_threshold(w, pi * lam * delta) / dV

        # --- beta (full system: M'M is not diagonal) ---
        Vb = MtM + 2 * pi * lam * phi * q * Ip
        wb = M.T @ (Y - X @ gamma)
        for j in range(q):
            Vb = Vb + rho * np.diag(alpha_new[j] ** 2)
            wb = wb + alpha_new[j] * tau[j] + rho * alpha_new[j] * mu[j]
        beta_new = np.linalg.solve(Vb, soft_threshold(wb, pi * lam * delta))

        # --- gamma ---
        gamma_new = np.linalg.solve(XtX, soft_threshold(X.T @ (Y - M @ beta_new), lam))

        # --- dual update ---
        tau = tau + rho * (mu - alpha_new * beta_new[None, :])

        diff = max(np.max(np.abs(alpha_new - alpha)),
                   np.max(np.abs(beta_new - beta)),
                   np.max(np.abs(gamma_new - gamma)))
        alpha, beta, gamma = alpha_new, beta_new, gamma_new
        rho = rho + rho0

    converged = s <= max_itr
    if not converged:
        warnings.warn("hdmediation: method does not converge")

    return {"lambda": lam, "IE": mu, "alpha": alpha, "beta": beta, "gamma": gamma,
            "logLik": hd_obj(X, M, Y, alpha, beta, gamma, lam, pi, phi, delta),
            "converge": converged, "nitr": s}


def hdmediation(X, M, Y, lam=0.0, pi=1.0, phi=2.0, delta=0.1, rho=1.0,
                standardize=True, max_itr=10000, tol=1e-6, rho_increase=False,
                alpha0=None, beta0=None, gamma0=None):
    """High-dimensional exposures and mediators (== R ``hdmediation``).

    With ``standardize=True`` the ADMM runs on standardized data and the
    estimates are mapped back to the original scale; the standardized fit is
    returned under ``out_scaled``.

    ``lambda`` is spelled ``lam`` here (``lambda`` is a Python keyword).

    Note: this estimator regresses on all ``q`` exposures directly, so it needs
    ``n > q``. For ``q`` comparable to ``n`` use :func:`hdmediation_pca`.
    """
    X = np.asarray(X, dtype=float)
    M = np.asarray(M, dtype=float)
    Y = np.asarray(Y, dtype=float).ravel()
    n, q = X.shape
    p = M.shape[1]
    if not standardize:
        return hdmediation_std(X, M, Y, lam, pi, phi, delta, rho, max_itr, tol,
                               rho_increase, alpha0, beta0, gamma0)

    x_sd = X.std(axis=0, ddof=1)
    m_sd = M.std(axis=0, ddof=1)
    y_sd = Y.std(ddof=1)
    Xs = (X - X.mean(axis=0)) / x_sd
    Ms = (M - M.mean(axis=0)) / m_sd
    Ys = (Y - Y.mean()) / y_sd

    r = hdmediation_std(Xs, Ms, Ys, lam, pi, phi, delta, rho, max_itr, tol,
                        rho_increase, alpha0, beta0, gamma0)

    alpha = r["alpha"] * (m_sd[None, :] / x_sd[:, None])
    IE = r["IE"] * (y_sd / x_sd[:, None])
    beta = r["beta"] * (y_sd / m_sd)
    gamma = r["gamma"] * (y_sd / x_sd)
    return {"lambda": lam, "IE": IE, "alpha": alpha, "beta": beta, "gamma": gamma,
            "logLik": hd_obj(X, M, Y, alpha, beta, gamma, lam, pi, phi, delta),
            "converge": r["converge"], "out_scaled": r}


def hdmediation_pca(X, M, Y, adaptive=False, var_prop=0.9, n_pc=None, **kwargs):
    """Principal-component variant (== R ``hdmediation_pca``).

    Replaces the exposures by their leading principal components before fitting,
    which is the appropriate route when the number of exposures is comparable to
    the sample size. Returns the fit on the component scale plus the loadings
    (``rotation``) needed to map back.
    """
    X = np.asarray(X, dtype=float)
    Xc = X - X.mean(axis=0)
    U, sv, Vt = np.linalg.svd(Xc, full_matrices=False)
    var = sv ** 2
    prop = np.cumsum(var) / np.sum(var)
    if n_pc is None or adaptive:
        k = int(np.searchsorted(prop, var_prop) + 1)
        if n_pc is not None:
            k = min(k, int(n_pc)) if adaptive else int(n_pc)
    else:
        k = int(n_pc)
    k = max(1, min(k, Xc.shape[1]))
    scores_ = Xc @ Vt[:k].T
    out = hdmediation(scores_, M, Y, **kwargs)
    out["n_pc"] = k
    out["rotation"] = Vt[:k].T
    out["var_explained"] = float(prop[k - 1])
    return out
