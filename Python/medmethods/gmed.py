"""Mediation analysis with a graph (covariance-matrix) mediator.

The exposure shifts the log-variance of a subject-level covariance mediator along
a projection ``theta``, and the outcome depends on that log-variance:

    M model:  log(theta' Sigma_i theta) = alpha0 + X_i' alpha + b_i
    Y model:  Y_i = gamma0 + X_i' gamma + beta * log(theta' Sigma_i theta) + e_i

with indirect effect ``IE = alpha_x * beta`` and direct effect ``DE = gamma_x``.

Ported from the R implementation behind ``MedMethods::gmed()``.

Random-effect note
------------------
R fits the M model with ``nlme::lme(score ~ X, random = ~1|ID)`` where every
subject contributes a *single* score, so the random-intercept groups are
singletons and the random-intercept / residual variance split is not
identifiable -- the marginal likelihood depends only on their sum, and R's
optimiser returns an arbitrary split. This port uses a deterministic
random-intercept fit: the fixed effects are GLS = OLS (identifiable) and
``blup_shrink`` controls the BLUP. The identifiable quantities -- ``theta``,
``alpha``, ``beta``, ``gamma`` and ``IE`` -- agree with R; only the
non-identifiable ``alpha0_rnd`` / ``tau2`` can differ.
"""
from __future__ import annotations

import numpy as np
from scipy.stats import norm

from ._core import cov_cube, scores, accum, gamma_solve, ginv_solve, norm_col, diag_level

__all__ = ["gmed_coef", "gmed_d1", "gmed", "gmed_boot", "gmed_refit"]


def _as2d(X):
    X = np.asarray(X, dtype=float)
    return X[:, None] if X.ndim == 1 else X


def gmed_coef(X, M, Y, theta, M_cov=None, nT=None, blup_shrink=1.0):
    """Path coefficients at a fixed projection (== R ``CAPMediation_coef``).

    ``X`` is the exposure/covariate matrix WITHOUT an intercept column (the
    intercept is added internally), its first column being the exposure of
    interest.
    """
    X = _as2d(X)
    Y = np.asarray(Y, dtype=float).ravel()
    n = Y.size
    if M_cov is None:
        M_cov, nT = cov_cube(M, ddof=1)
    theta = np.asarray(theta, dtype=float).ravel()
    score = scores(M_cov, theta)
    logs = np.log(score)

    # Y model: Y ~ [1, X, log(score)]
    Z = np.column_stack([np.ones(n), X, logs])
    mu = ginv_solve(Z.T @ Z, Z.T @ Y)
    gamma0 = float(mu[0])
    gamma = mu[1:-1]
    beta = float(mu[-1])
    sigma2 = float(np.mean((Y - Z @ mu) ** 2))

    # M model fixed effects (GLS == OLS with singleton groups)
    Zm = np.column_stack([np.ones(n), X])
    am = ginv_solve(Zm.T @ Zm, Zm.T @ logs)
    alpha0 = float(am[0])
    alpha = am[1:]
    fitted = Zm @ am
    # per-subject intercept only: the X-slope enters separately via X @ alpha
    alpha0_rnd = alpha0 + blup_shrink * (logs - fitted)
    tau2 = float(np.mean((alpha0_rnd - alpha0) ** 2))

    return {
        "theta": theta, "alpha": alpha, "beta": beta, "gamma": gamma,
        "IE": float(alpha[0] * beta), "DE": float(gamma[0]),
        "alpha0": alpha0, "alpha0_rnd": alpha0_rnd, "gamma0": gamma0,
        "tau2": tau2, "sigma2": sigma2, "score": score,
    }


def _obj(X, M_cov, Y, nT, theta, c):
    """Objective (== R ``obj.func``)."""
    X = _as2d(X)
    Y = np.asarray(Y, dtype=float).ravel()
    score = scores(M_cov, theta)
    lp = c["alpha0_rnd"] + X @ c["alpha"]
    l1 = np.sum((lp + score * np.exp(-lp)) * nT) / 2.0
    l2 = np.sum((Y - c["gamma0"] - X @ c["gamma"] - c["beta"] * np.log(score)) ** 2
                / c["sigma2"] + np.log(c["sigma2"])) / 2.0
    l3 = np.sum((c["alpha0_rnd"] - c["alpha0"]) ** 2 / c["tau2"]
                + np.log(c["tau2"])) / 2.0
    return float(l1 + l2 + l3)


def gmed_d1(X, M, Y, H=None, theta0=None, max_itr=1000, tol=1e-4,
            blup_shrink=1.0, M_cov=None, nT=None):
    """One mediating direction from one start (== R ``CAPMediation_D1``)."""
    X = _as2d(X)
    Y = np.asarray(Y, dtype=float).ravel()
    if M_cov is None:
        M_cov, nT = cov_cube(M, ddof=1)
    p = M_cov.shape[0]
    if H is None:
        H = M_cov.mean(axis=2)
    theta = np.full(p, 1.0 / np.sqrt(p)) if theta0 is None else \
        np.asarray(theta0, dtype=float).ravel().copy()

    c = gmed_coef(X, M, Y, theta, M_cov=M_cov, nT=nT, blup_shrink=blup_shrink)
    s, diff = 0, 100.0
    while s <= max_itr and diff > tol:
        s += 1
        score = scores(M_cov, theta)
        U = np.exp(-c["alpha0_rnd"] - X @ c["alpha"])
        V = Y - c["gamma0"] - X @ c["gamma"]
        beta = c["beta"]
        w = nT * U - (2 * beta * (V - beta * np.log(score))) / (c["sigma2"] * score)
        theta_new = gamma_solve(accum(M_cov, w), H)
        c_new = gmed_coef(X, M, Y, theta_new, M_cov=M_cov, nT=nT,
                          blup_shrink=blup_shrink)
        diff = np.max(np.abs(np.concatenate([
            c_new["alpha"] - c["alpha"], [c_new["beta"] - c["beta"]],
            c_new["gamma"] - c["gamma"]])))
        c, theta = c_new, theta_new

    theta = norm_col(theta)
    out = gmed_coef(X, M, Y, theta, M_cov=M_cov, nT=nT, blup_shrink=blup_shrink)
    out["obj"] = _obj(X, M_cov, Y, nT, theta, out)
    out["convergence"] = s < max_itr
    out["nitr"] = s
    return out


def _theta0_mat(p, ninitial, seed):
    rng = np.random.RandomState(seed)
    T = rng.normal(size=(p, p + 1 + 5))
    T = T / np.sqrt((T ** 2).sum(axis=0, keepdims=True))
    return T[:, :ninitial]


def _d1_opt(X, M, Y, M_cov, nT, H, max_itr, tol, ninitial, seed, blup_shrink):
    p = M_cov.shape[0]
    if ninitial is None:
        ninitial = min(p, 10)
    T = _theta0_mat(p, ninitial, seed)
    best, best_obj = None, np.inf
    for k in range(T.shape[1]):
        try:
            rk = gmed_d1(X, M, Y, H=H, theta0=T[:, k], max_itr=max_itr, tol=tol,
                         blup_shrink=blup_shrink, M_cov=M_cov, nT=nT)
        except Exception:
            continue
        # Select on the objective at the H-UNSCALED theta (as R does): the
        # 1/sqrt(theta'H theta) rescaling penalises high-variance background
        # directions so the mediating direction wins.
        try:
            g = rk["theta"]
            g_un = g / np.sqrt(g @ H @ g)
            c_un = gmed_coef(X, M, Y, g_un, M_cov=M_cov, nT=nT,
                             blup_shrink=blup_shrink)
            obj = _obj(X, M_cov, Y, nT, g_un, c_un)
        except Exception:
            obj = np.inf
        if obj < best_obj:
            best_obj, best = obj, rk
    if best is None:
        raise RuntimeError("gmed: every starting point failed")
    return best


def gmed(X, M, Y, H=None, stop_crt="nD", nD=1, DfD_thred=2.0, max_itr=1000,
         tol=1e-4, ninitial=None, seed=100, blup_shrink=1.0, verbose=False):
    """Graph-mediator mediation analysis (== R ``gmed``).

    Parameters
    ----------
    X : (n,) or (n, nX) exposure and covariates, no intercept column; the first
        column is the exposure of interest.
    M : length-n sequence of ``(T_i, p)`` mediator data matrices.
    Y : (n,) outcome.
    H : ``(p, p)`` constraint matrix; ``None`` uses the average covariance, which
        is what R's ``H = NULL`` computes and what recovers the mediating
        direction (an identity ``H`` collapses onto a background direction).
    """
    X = _as2d(X)
    Y = np.asarray(Y, dtype=float).ravel()
    M_cov, nT = cov_cube(M, ddof=1)
    if H is None:
        H = M_cov.mean(axis=2)

    r1 = _d1_opt(X, M, Y, M_cov, nT, H, max_itr, tol, ninitial, seed, blup_shrink)
    Theta, res = [r1["theta"]], [r1]
    if verbose:
        print("Component 1")

    def _deflate(T0):
        P = T0 @ T0.T
        Mt = [np.asarray(Mi, dtype=float) - np.asarray(Mi, dtype=float) @ P
              for Mi in M]
        Sc, nt = cov_cube(Mt, ddof=1)
        return Mt, Sc, nt

    if stop_crt == "nD":
        for _ in range(2, nD + 1):
            Mt, Sc, nt = _deflate(np.column_stack(Theta))
            try:
                rk = _d1_opt(X, Mt, Y, Sc, nt, Sc.mean(axis=2), max_itr, tol,
                             ninitial, seed, blup_shrink)
            except Exception:
                break
            Theta.append(rk["theta"])
            res.append(rk)
            if verbose:
                print("Component %d" % len(Theta))
    elif stop_crt == "DfD":
        DfD = 1.0
        while DfD < DfD_thred:
            Mt, Sc, nt = _deflate(np.column_stack(Theta))
            try:
                rk = _d1_opt(X, Mt, Y, Sc, nt, Sc.mean(axis=2), max_itr, tol,
                             ninitial, seed, blup_shrink)
            except Exception:
                break
            cand = np.column_stack(Theta + [rk["theta"]])
            DfD = diag_level(M_cov, cand, nT)["avg_level"][-1]
            if np.isfinite(DfD) and DfD < DfD_thred:
                Theta.append(rk["theta"])
                res.append(rk)
                if verbose:
                    print("Component %d" % len(Theta))
            else:
                break
    else:
        raise ValueError("stop_crt must be 'nD' or 'DfD'")

    Tm = np.column_stack(Theta)
    out = {
        "theta": Tm,
        "alpha": np.column_stack([r["alpha"] for r in res]),
        "beta": np.array([r["beta"] for r in res]),
        "gamma": np.column_stack([r["gamma"] for r in res]),
        "IE": np.array([r["IE"] for r in res]),
        "DE": np.array([r["DE"] for r in res]),
        "components": ["C%d" % (i + 1) for i in range(Tm.shape[1])],
    }
    if Tm.shape[1] >= 2:
        out["DfD"] = diag_level(M_cov, Tm, nT)
        out["orthogonality"] = Tm.T @ Tm
    return out


def gmed_refit(X, M, Y, Theta, blup_shrink=1.0):
    """Path coefficients at fixed projections (== R ``gmed_refit``).

    Note: R's ``CAPMediation_coef`` / ``CAPMediation_refit`` referenced undefined
    ``n`` and ``p`` (they leaked in from the caller's workspace); the R package
    derives them from ``M``, as this port does.
    """
    Theta = np.asarray(Theta, dtype=float)
    if Theta.ndim == 1:
        Theta = Theta[:, None]
    M_cov, nT = cov_cube(M, ddof=1)
    out = []
    for k in range(Theta.shape[1]):
        out.append(gmed_coef(X, M, Y, Theta[:, k], M_cov=M_cov, nT=nT,
                             blup_shrink=blup_shrink))
    return {
        "theta": Theta,
        "alpha": np.column_stack([o["alpha"] for o in out]),
        "beta": np.array([o["beta"] for o in out]),
        "gamma": np.column_stack([o["gamma"] for o in out]),
        "IE": np.array([o["IE"] for o in out]),
        "per_direction": out,
    }


def gmed_boot(X, M, Y, theta, sims=1000, conf_level=0.95, boot_ci_type="se",
              seed=100, blup_shrink=1.0):
    """Bootstrap inference at a fixed projection (== R ``gmed_boot``).

    Resamples subjects and re-estimates the path coefficients with ``theta`` held
    fixed. R seeds each replicate with ``set.seed(seed.boot + b)`` and uses R's
    sampler, so individual replicates differ from this port's; standard errors
    agree up to Monte Carlo error.
    """
    X = _as2d(X)
    Y = np.asarray(Y, dtype=float).ravel()
    n = Y.size
    theta = np.asarray(theta, dtype=float).ravel()
    rng = np.random.RandomState(seed)
    names = ["alpha", "beta", "gamma", "IE", "DE"]
    draws = np.full((sims, 5), np.nan)
    for b in range(sims):
        idx = rng.randint(0, n, size=n)
        try:
            c = gmed_coef(X[idx], [M[i] for i in idx], Y[idx], theta,
                          blup_shrink=blup_shrink)
        except Exception:
            continue
        draws[b] = [c["alpha"][0], c["beta"], c["gamma"][0], c["IE"], c["gamma"][0]]

    est = np.nanmean(draws, axis=0)
    se = np.nanstd(draws, axis=0, ddof=1)
    stat = est / se
    pval = 2 * (1 - norm.cdf(np.abs(stat)))
    if boot_ci_type == "perc":
        lo = np.nanpercentile(draws, 100 * (1 - conf_level) / 2, axis=0)
        hi = np.nanpercentile(draws, 100 * (1 - (1 - conf_level) / 2), axis=0)
    else:
        z = norm.ppf(1 - (1 - conf_level) / 2)
        lo, hi = est - z * se, est + z * se
    return {"names": names, "estimate": est, "SE": se, "statistic": stat,
            "pvalue": pval, "LB": lo, "UB": hi, "boot": draws}
