"""Principal component mediation analysis for multiple exposures and mediators.

Orthogonal projections ``phi`` (exposures) and ``psi`` (mediators) define scalar
component scores that carry parallel mediation mechanisms:

    M psi = (X phi) alpha + e,      e ~ N(0, sigma2)
    Y     = (X phi) gamma + (M psi) beta + u,   u ~ N(0, tau2)

with indirect effect ``IE = alpha * beta``. ``pcma`` minimizes the joint negative
log-likelihood over ``(phi, psi, alpha, beta, gamma, sigma2, tau2)`` by block
updates, restarting from several initial projections, then deflates for higher
components.

Ported from the R ``HDEM.*`` functions behind ``MedMethods::pcma()``.
"""
from __future__ import annotations

import numpy as np
from scipy.optimize import brentq
from scipy.stats import norm

from ._core import ginv_solve, norm_col

__all__ = [
    "pcma_coef", "pcma_loglike", "pcma_d1_base", "pcma_d1", "pcma", "pcma_inf",
]


def _center(A):
    A = np.asarray(A, dtype=float)
    return A - A.mean(axis=0)


def _resid(y, x):
    """Residual of ``y`` after projecting out the single column ``x``."""
    b = (x @ y) / (x @ x)
    return y - b * x


def pcma_coef(X, M, Y, phi, psi):
    """Path coefficients at fixed projections (== R ``HDEM.coef``).

    ``alpha`` regresses the mediator score on the exposure score; ``beta`` and
    ``gamma`` are the partial slopes in the outcome model. Implemented by
    residualisation rather than the R code's explicit ``n x n`` projectors --
    algebraically identical, but O(n) instead of O(n^2) memory.
    """
    X = np.asarray(X, dtype=float)
    M = np.asarray(M, dtype=float)
    Y = np.asarray(Y, dtype=float).ravel()
    phi = np.asarray(phi, dtype=float).ravel()
    psi = np.asarray(psi, dtype=float).ravel()

    xt = _center(X @ phi)
    mt = _center(M @ psi)
    y = Y - Y.mean()

    alpha = (xt @ mt) / (xt @ xt)
    # beta: regress y on mt after removing xt from both  (Hx = I - P_xt)
    mt_x, y_x = _resid(mt, xt), _resid(y, xt)
    beta = (mt_x @ y_x) / (mt_x @ mt_x)
    # gamma: regress y on xt after removing mt from both (Hm = I - P_mt)
    xt_m, y_m = _resid(xt, mt), _resid(y, mt)
    gamma = (xt_m @ y_m) / (xt_m @ xt_m)

    return {"alpha": float(alpha), "beta": float(beta), "gamma": float(gamma),
            "IE": float(alpha * beta)}


def pcma_loglike(X, M, Y, phi, psi, alpha, beta, gamma, sigma2, tau2):
    """Joint negative log-likelihood (== R ``HDEM.loglike``)."""
    X = np.asarray(X, dtype=float)
    M = np.asarray(M, dtype=float)
    Y = np.asarray(Y, dtype=float).ravel()
    n = Y.size
    xt = X @ np.asarray(phi, dtype=float).ravel()
    mt = M @ np.asarray(psi, dtype=float).ravel()
    s1 = np.sum((mt - xt * alpha) ** 2) / sigma2
    s2 = np.sum((Y - xt * gamma - mt * beta) ** 2) / tau2
    s3 = n * np.log(sigma2) + n * np.log(tau2)
    return float(s1 + s2 + s3)


def _unit_ridge(A, U, upper=1e6):
    """Solve ``(A + lam I)^{-1} U`` with ``lam >= 0`` chosen so the result is unit norm.

    Mirrors the R ``uniroot(..., lower = 0, upper = 1e6)`` on
    ``|| (A + lam I)^{-1} U ||^2 - 1``; if no sign change is bracketed, R falls
    back to ``lam = 0``.
    """
    d = A.shape[0]
    I = np.eye(d)

    def f(lam):
        v = ginv_solve(A + lam * I, U)
        return float(v @ v) - 1.0

    lam = 0.0
    try:
        f0, f1 = f(0.0), f(upper)
        if np.isfinite(f0) and np.isfinite(f1) and f0 * f1 < 0:
            lam = brentq(f, 0.0, upper)
    except Exception:
        lam = 0.0
    v = ginv_solve(A + lam * I, U)
    nrm = np.sqrt(np.sum(v ** 2))
    if nrm != 1.0 and nrm > 0:
        v = v / nrm
    return v, lam


def pcma_d1_base(X, M, Y, max_itr=1000, tol=1e-4, phi0=None, psi0=None):
    """One component from one starting point (== R ``HDEM.loglike.D1.base``)."""
    X = np.asarray(X, dtype=float)
    M = np.asarray(M, dtype=float)
    Y = np.asarray(Y, dtype=float).ravel()
    n, p = X.shape
    q = M.shape[1]

    phi = np.full(p, 1.0 / np.sqrt(p)) if phi0 is None else \
        np.asarray(phi0, dtype=float).ravel().copy()
    psi = np.full(q, 1.0 / np.sqrt(q)) if psi0 is None else \
        np.asarray(psi0, dtype=float).ravel().copy()

    a_old = b_old = g_old = 0.0
    XtX = X.T @ X
    MtM = M.T @ M
    s = 0
    diff = 100.0
    alpha = beta = gamma = 0.0
    sigma2 = tau2 = 1.0
    while s <= max_itr and diff > tol:
        s += 1
        c = pcma_coef(X, M, Y, phi, psi)
        alpha, beta, gamma = c["alpha"], c["beta"], c["gamma"]

        xt = X @ phi
        mt = M @ psi
        sigma2 = np.sum((mt - xt * alpha) ** 2) / n
        tau2 = np.sum((Y - xt * gamma - mt * beta) ** 2) / n

        # phi update
        A = (alpha ** 2 / sigma2 + gamma ** 2 / tau2) * XtX
        U = X.T @ ((alpha / sigma2) * mt + (gamma / tau2) * Y
                   - (beta * gamma / tau2) * mt)
        phi, _ = _unit_ridge(A, U)

        # psi update (uses the NEW phi, as in R)
        xt_new = X @ phi
        B = (1.0 / sigma2 + beta ** 2 / tau2) * MtM
        V = M.T @ ((alpha / sigma2 - beta * gamma / tau2) * xt_new
                   + (beta / tau2) * Y)
        psi, _ = _unit_ridge(B, V)

        diff = max(abs(a_old - alpha), abs(b_old - beta), abs(g_old - gamma))
        a_old, b_old, g_old = alpha, beta, gamma

    if phi[np.argmax(np.abs(phi))] < 0:
        phi = -phi
    if psi[np.argmax(np.abs(psi))] < 0:
        psi = -psi

    return {
        "alpha": alpha, "beta": beta, "gamma": gamma, "phi": phi, "psi": psi,
        "sigma2": float(sigma2), "tau2": float(tau2),
        "logLik": pcma_loglike(X, M, Y, phi, psi, alpha, beta, gamma, sigma2, tau2),
        "convergence": s < max_itr, "nitr": s,
    }


def _init_mat(d, ncol, seed):
    rng = np.random.RandomState(seed)
    T = rng.normal(size=(d, ncol))
    return T / np.sqrt((T ** 2).sum(axis=0, keepdims=True))


def pcma_d1(X, M, Y, max_itr=1000, tol=1e-4, phi0_mat=None, psi0_mat=None,
            ninitial=None, seed=100):
    """Leading component, best of several starts (== R ``HDEM.loglike.D1``).

    The starting projections come from numpy's RNG, so the individual starts
    differ from R's; the selected optimum (minimum negative log-likelihood) is
    what matters. Pass ``phi0_mat`` / ``psi0_mat`` to use identical starts.
    """
    X = np.asarray(X, dtype=float)
    M = np.asarray(M, dtype=float)
    p, q = X.shape[1], M.shape[1]
    ncol = max(p, q) + 1 + 5
    if phi0_mat is None:
        phi0_mat = _init_mat(p, ncol, seed)
    else:
        phi0_mat = np.atleast_2d(np.asarray(phi0_mat, dtype=float))
    if psi0_mat is None:
        psi0_mat = _init_mat(q, ncol, seed)
    else:
        psi0_mat = np.atleast_2d(np.asarray(psi0_mat, dtype=float))

    navail = min(phi0_mat.shape[1], psi0_mat.shape[1])
    if ninitial is None:
        ninitial = min(navail, 10)
    ninitial = min(ninitial, navail)

    best, best_obj = None, np.inf
    for k in range(ninitial):
        try:
            r = pcma_d1_base(X, M, Y, max_itr=max_itr, tol=tol,
                             phi0=phi0_mat[:, k], psi0=psi0_mat[:, k])
        except Exception:
            continue
        if np.isfinite(r["logLik"]) and r["logLik"] < best_obj:
            best_obj, best = r["logLik"], r
    if best is None:
        raise RuntimeError("pcma: every starting point failed")
    return best


def pcma_inf(X, M, Y, phi, psi, sims=1000, conf_level=0.95,
             boot_ci_type="perc", seed=100):
    """Bootstrap inference at fixed projections (== R ``HDEM.inf``).

    Resamples subjects and re-estimates the path coefficients with ``phi`` and
    ``psi`` held fixed. R seeds each replicate with ``set.seed(100 + b)`` and
    uses R's sampler, so individual replicates differ from this port's; the
    resulting standard errors agree up to Monte Carlo error.
    """
    X = np.asarray(X, dtype=float)
    M = np.asarray(M, dtype=float)
    Y = np.asarray(Y, dtype=float).ravel()
    n = Y.size
    rng = np.random.RandomState(seed)
    out = np.full((sims, 5), np.nan)
    for b in range(sims):
        idx = rng.randint(0, n, size=n)
        try:
            c = pcma_coef(X[idx], M[idx], Y[idx], phi, psi)
        except Exception:
            continue
        out[b] = [c["alpha"], c["beta"], c["gamma"], c["IE"], c["gamma"]]

    names = ["alpha", "beta", "gamma", "IE", "DE"]
    est = np.nanmean(out, axis=0)
    se = np.nanstd(out, axis=0, ddof=1)
    stat = est / se
    pval = 2 * (1 - norm.cdf(np.abs(stat)))
    if boot_ci_type == "perc":
        lo = np.nanpercentile(out, 100 * (1 - conf_level) / 2, axis=0)
        hi = np.nanpercentile(out, 100 * (1 - (1 - conf_level) / 2), axis=0)
    else:
        z = norm.ppf(1 - (1 - conf_level) / 2)
        lo, hi = est - z * se, est + z * se
    return {"names": names, "estimate": est, "SE": se, "statistic": stat,
            "pvalue": pval, "LB": lo, "UB": hi, "boot": out}


def pcma(X, M, Y, nD=1, max_itr=1000, tol=1e-4, ninitial=None, seed=100,
         boot=False, sims=1000, conf_level=0.95, boot_ci_type="perc",
         verbose=False):
    """Principal component mediation analysis (== R ``pcma`` with ``stop.crt="nD"``).

    Returns the exposure projections ``Phi`` (p x nD), mediator projections
    ``Psi`` (q x nD), the per-component path coefficients, and -- when
    ``boot=True`` -- bootstrap inference per component.
    """
    X = np.asarray(X, dtype=float)
    M = np.asarray(M, dtype=float)
    Y = np.asarray(Y, dtype=float).ravel()

    r1 = pcma_d1(X, M, Y, max_itr=max_itr, tol=tol, ninitial=ninitial, seed=seed)
    Phi, Psi = [r1["phi"]], [r1["psi"]]
    res = [r1]
    inference = []
    if boot:
        inference.append(pcma_inf(X, M, Y, r1["phi"], r1["psi"], sims=sims,
                                  conf_level=conf_level,
                                  boot_ci_type=boot_ci_type, seed=seed))
    if verbose:
        print("Component 1")

    for j in range(2, nD + 1):
        P0 = np.column_stack(Phi)
        S0 = np.column_stack(Psi)
        b0 = np.array([r["beta"] for r in res])
        g0 = np.array([r["gamma"] for r in res])
        # deflate: strip identified components from X, M and their effect from Y
        Xt = X - X @ (P0 @ P0.T)
        Mt = M - M @ (S0 @ S0.T)
        Yt = Y - (X @ P0) @ g0 - (M @ S0) @ b0
        try:
            rk = pcma_d1(Xt, Mt, Yt, max_itr=max_itr, tol=tol,
                         ninitial=ninitial, seed=seed)
        except Exception:
            break
        rk["orthogonality"] = {"phi": rk["phi"] @ P0, "psi": rk["psi"] @ S0}
        Phi.append(rk["phi"])
        Psi.append(rk["psi"])
        res.append(rk)
        if boot:
            inference.append(pcma_inf(Xt, Mt, Yt, rk["phi"], rk["psi"], sims=sims,
                                      conf_level=conf_level,
                                      boot_ci_type=boot_ci_type, seed=seed))
        if verbose:
            print("Component %d" % len(Phi))

    out = {
        "Phi": np.column_stack(Phi), "Psi": np.column_stack(Psi),
        "alpha": np.array([r["alpha"] for r in res]),
        "beta": np.array([r["beta"] for r in res]),
        "gamma": np.array([r["gamma"] for r in res]),
        "IE": np.array([r["alpha"] * r["beta"] for r in res]),
        "sigma2": np.array([r["sigma2"] for r in res]),
        "tau2": np.array([r["tau2"] for r in res]),
        "logLik": np.array([r["logLik"] for r in res]),
        "components": ["D%d" % (i + 1) for i in range(len(Phi))],
    }
    if boot:
        out["inference"] = inference
    if len(Phi) > 1:
        out["Phi_orthogonality"] = np.column_stack(Phi).T @ np.column_stack(Phi)
        out["Psi_orthogonality"] = np.column_stack(Psi).T @ np.column_stack(Psi)
    return out
