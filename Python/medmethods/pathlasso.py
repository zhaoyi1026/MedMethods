"""Pathway Lasso: pathway estimation and selection with high-dimensional mediators.

    M = Z A + E,        Y = Z C + M B + e

The penalty acts on the *products* ``A_j B_j`` (the pathway effects) rather than
on ``A`` and ``B`` separately, plus ridge (``phi``) and l1 (``omega``) terms.
Solved by ADMM with the splittings ``Theta = alpha`` and ``D = beta``, where the
``(alpha_j, beta_j)`` block has a closed-form proximal solution.

Ported from the R implementation behind ``MedMethods::pathlasso()``.

The ``lambda`` scale matters: the dense-to-empty transition is abrupt, and on
standardized data the useful range is small (the R package's example uses
``lambda = 1e-3``). Select over a path -- see :func:`pathlasso_ksc`.
"""
from __future__ import annotations

import warnings

import numpy as np

__all__ = ["pathlasso", "pathlasso_admm", "prox_pair", "pathlasso_ksc"]


def _soft(mu, lam):
    return max(abs(mu) - lam, 0.0) * np.sign(mu)


def prox_pair(lam, phi1, phi2, mu1, mu2, omega1=0.0, omega2=0.0):
    """Closed-form minimiser of the penalised pair problem (== R ``solution``).

    Minimises
    ``phi1 a^2/2 + phi2 b^2/2 - mu1 a - mu2 b + lam |a||b| + omega1|a| + omega2|b|``
    by checking the four smooth branches and then the axis / origin solutions.
    """
    if lam == 0:
        return _soft(mu1, omega1) / phi1, _soft(mu2, omega2) / phi2

    de = phi1 * phi2 - lam ** 2
    x1 = phi2 * (mu1 - omega1) - lam * (mu2 - omega2)
    x2 = phi2 * (mu1 - omega1) + lam * (mu2 + omega2)
    x3 = phi2 * (mu1 + omega1) + lam * (mu2 - omega2)
    x4 = phi2 * (mu1 + omega1) - lam * (mu2 + omega2)
    y1 = phi1 * (mu2 - omega2) - lam * (mu1 - omega1)
    y2 = phi1 * (mu2 + omega2) + lam * (mu1 - omega1)
    y3 = phi1 * (mu2 - omega2) + lam * (mu1 + omega1)
    y4 = phi1 * (mu2 + omega2) - lam * (mu1 + omega1)

    if x1 > 0 and y1 > 0:
        return x1 / de, y1 / de
    if x2 > 0 and y2 < 0:
        return x2 / de, y2 / de
    if x3 < 0 and y3 > 0:
        return x3 / de, y3 / de
    if x4 < 0 and y4 < 0:
        return x4 / de, y4 / de
    if abs(mu1) > omega1 and (phi1 * abs(mu2) - lam * abs(mu1)) <= (-lam * omega1 + phi1 * omega2):
        return _soft(mu1, omega1) / phi1, 0.0
    if abs(mu2) > omega2 and (phi2 * abs(mu1) - lam * abs(mu2)) <= (-lam * omega2 + phi2 * omega1):
        return 0.0, _soft(mu2, omega2) / phi2
    return 0.0, 0.0


def pathlasso_admm(Z, M, R, lam=1.0, omega=0.0, phi=1.0, Phi1=None, Phi2=None,
                   rho=1.0, rho_increase=False, tol=1e-10, max_itr=10000,
                   thred=1e-10, est_thred=False, Sigma1=None, Sigma2=None,
                   Theta0=None, D0=None, alpha0=None, beta0=None):
    """ADMM on (already standardized) data (== R ``mediation_net_ADMM_NC``)."""
    Z = np.asarray(Z, dtype=float).reshape(-1, 1)
    M = np.asarray(M, dtype=float)
    R = np.asarray(R, dtype=float).reshape(-1, 1)
    n, k = M.shape

    if Phi1 is None:
        Phi1 = np.diag(np.concatenate([[0.0], np.ones(k)]))
    else:
        Phi1 = np.asarray(Phi1, dtype=float)
    if Phi2 is None:
        Phi2 = np.eye(k + 1)
    else:
        Phi2 = np.asarray(Phi2, dtype=float)
    ph1d, ph2d = np.diag(Phi1).copy(), np.diag(Phi2).copy()

    X = np.hstack([Z, M])
    e1 = np.concatenate([[1.0], np.zeros(k)])
    ZtZ = float(Z.T @ Z)

    Theta = np.zeros((1, k + 1)) if Theta0 is None else \
        np.array(Theta0, dtype=float).reshape(1, k + 1)
    alpha = np.zeros((1, k + 1)) if alpha0 is None else \
        np.array(alpha0, dtype=float).reshape(1, k + 1)
    D = np.zeros((k + 1, 1)) if D0 is None else \
        np.array(D0, dtype=float).reshape(k + 1, 1)
    beta = np.zeros((k + 1, 1)) if beta0 is None else \
        np.array(beta0, dtype=float).reshape(k + 1, 1)

    def _S1(Th):
        Rm = M - Z @ Th[:, 1:]
        return np.diag(np.diag(Rm.T @ Rm / n))

    def _S2(Dm):
        Rr = R - X @ Dm
        return float(Rr.T @ Rr)

    S1 = _S1(Theta) if Sigma1 is None else np.asarray(Sigma1, dtype=float)
    S2 = _S2(D) if Sigma2 is None else float(np.asarray(Sigma2))

    nu1 = np.zeros(k + 1)
    nu2 = np.zeros((k + 1, 1))
    nu3 = 0.0
    rho0 = rho if rho_increase else 0.0

    J = np.concatenate([[0.0], np.ones(k)])
    if lam == 0:
        J[0] = 1.0

    Ik1 = np.eye(k + 1)
    e1e1 = np.outer(e1, e1)
    diff, s = 100.0, 0
    Theta_prev = Theta.copy()
    D_prev = D.copy()
    while diff >= tol and s <= max_itr:
        s += 1
        Om1 = np.zeros((k + 1, k + 1))
        Om1[1:, 1:] = np.linalg.inv(S1) / n
        w2 = 1.0 / (S2 * n)

        # Theta
        de_T = np.linalg.inv(ZtZ * Om1 + 2 * rho * Ik1 + 2 * rho * e1e1)
        nu_T = (Z.T @ X) @ Om1 - nu1[None, :] + 2 * rho * alpha + (2 * rho - nu3) * e1[None, :]
        Theta = nu_T @ de_T
        if est_thred:
            Theta[np.abs(Theta) < thred] = 0.0

        # D
        de_D = np.linalg.inv(w2 * (X.T @ X) + 2 * rho * Ik1)
        nu_D = w2 * (X.T @ R) - nu2 + 2 * rho * beta
        D = de_D @ nu_D
        if est_thred:
            D[np.abs(D) < thred] = 0.0

        # (alpha_j, beta_j) blocks
        a_new = np.empty(k + 1)
        b_new = np.empty(k + 1)
        for j in range(k + 1):
            pp1 = 2 * lam * phi * ph1d[j] ** 2 + 2 * rho
            pp2 = 2 * lam * phi * ph2d[j] ** 2 + 2 * rho
            mu1 = nu1[j] + 2 * rho * Theta[0, j]
            mu2 = nu2[j, 0] + 2 * rho * D[j, 0]
            a_new[j], b_new[j] = prox_pair(
                lam * ph1d[j] * ph2d[j], pp1, pp2, mu1, mu2,
                omega1=omega * ph1d[j] * J[j], omega2=omega * ph2d[j] * J[j])
        alpha_new = a_new[None, :]
        beta_new = b_new[:, None]

        nu1 = nu1 + 2 * rho * (Theta - alpha_new).ravel()
        nu2 = nu2 + 2 * rho * (D - beta_new)
        nu3 = nu3 + 2 * rho * float(Theta @ e1 - 1.0)

        # convergence on the change in Theta and D since the previous iterate
        diff = max(np.max(np.abs(Theta - Theta_prev)), np.max(np.abs(D - D_prev)))
        Theta_prev, D_prev = Theta.copy(), D.copy()
        alpha, beta = alpha_new, beta_new

        if Sigma1 is None:
            S1 = _S1(Theta)
        if Sigma2 is None:
            S2 = _S2(D)

    converged = s <= max_itr
    if not converged:
        warnings.warn("pathlasso: method does not converge")

    A = alpha[0, 1:][None, :]
    C = np.array([[beta[0, 0]]])
    B = beta[1:, :]
    constraint = {
        "Theta=alpha": float(np.max(np.abs(Theta - alpha))),
        "D=beta": float(np.max(np.abs(D - beta))),
        "Theta[1]": float(Theta[0, 0]),
    }
    return {"lambda": lam, "omega": omega, "phi": phi, "rho": rho,
            "A": A, "B": B, "C": C, "AB": (A.ravel() * B.ravel()),
            "Theta": Theta, "D": D, "alpha": alpha, "beta": beta,
            "converge": converged, "nitr": s, "constraint": constraint}


def pathlasso(X, M, Y, lam=1.0, omega=0.0, phi=1.0, Phi1=None, Phi2=None,
              rho=1.0, rho_increase=False, standardize=True, tol=1e-10,
              max_itr=10000, thred=1e-10, est_thred=False, Sigma1=None,
              Sigma2=None, Theta0=None, D0=None, alpha0=None, beta0=None):
    """Pathway Lasso (== R ``pathlasso``).

    ``lambda`` is spelled ``lam`` here (``lambda`` is a Python keyword). With
    ``standardize=True`` the ADMM runs on standardized data and the estimates are
    mapped back to the original scale; the standardized fit is under
    ``out_scaled``.

    Returns the a-path ``A`` (1 x k), b-path ``B`` (k x 1), direct effect ``C``
    and the pathway effects ``AB`` (length k).
    """
    Xv = np.asarray(X, dtype=float).ravel()
    M = np.asarray(M, dtype=float)
    Yv = np.asarray(Y, dtype=float).ravel()
    if not standardize:
        return pathlasso_admm(Xv, M, Yv, lam, omega, phi, Phi1, Phi2, rho,
                              rho_increase, tol, max_itr, thred, est_thred,
                              Sigma1, Sigma2, Theta0, D0, alpha0, beta0)

    x_sd = Xv.std(ddof=1)
    m_sd = M.std(axis=0, ddof=1)
    y_sd = Yv.std(ddof=1)
    Xs = (Xv - Xv.mean()) / x_sd
    Ms = (M - M.mean(axis=0)) / m_sd
    Ys = (Yv - Yv.mean()) / y_sd

    r = pathlasso_admm(Xs, Ms, Ys, lam, omega, phi, Phi1, Phi2, rho,
                       rho_increase, tol, max_itr, thred, est_thred,
                       Sigma1, Sigma2, Theta0, D0, alpha0, beta0)
    A = r["A"] * (m_sd[None, :] / x_sd)
    B = r["B"] * (y_sd / m_sd)[:, None]
    C = r["C"] * (y_sd / x_sd)
    return {"lambda": lam, "omega": omega, "phi": phi, "rho": rho,
            "A": A, "B": B, "C": C, "AB": A.ravel() * B.ravel(),
            "converge": r["converge"], "constraint": r["constraint"],
            "out_scaled": r}


def _cohen_kappa(s1, s2):
    """Agreement between two binary selection vectors (== R ``cohen.kappa``)."""
    s1 = np.asarray(s1, dtype=bool)
    s2 = np.asarray(s2, dtype=bool)
    n = s1.size
    if n == 0:
        return np.nan
    po = np.mean(s1 == s2)
    p1 = np.mean(s1) * np.mean(s2) + (1 - np.mean(s1)) * (1 - np.mean(s2))
    if p1 == 1.0:
        return 1.0
    return (po - p1) / (1 - p1)


def pathlasso_ksc(X, M, Y, lam_grid, n_rep=5, zero_cutoff=1e-3, vss_cut=0.1,
                  seed=100, **kwargs):
    """Tuning by selection stability (== R ``pathlasso_ksc``, kappa criterion).

    For each ``lambda`` the sample is split in half ``n_rep`` times; the two
    halves are fitted independently and the agreement of their selected pathway
    sets is measured by Cohen's kappa. The chosen ``lambda`` is the largest one
    whose mean kappa is within ``vss_cut`` of the maximum, matching the R
    routine's "most parsimonious among the stable" rule.
    """
    Xv = np.asarray(X, dtype=float).ravel()
    M = np.asarray(M, dtype=float)
    Yv = np.asarray(Y, dtype=float).ravel()
    n = Yv.size
    lam_grid = np.atleast_1d(np.asarray(lam_grid, dtype=float))
    rng = np.random.RandomState(seed)
    kappa = np.full((len(lam_grid), n_rep), np.nan)
    for r in range(n_rep):
        perm = rng.permutation(n)
        h1, h2 = perm[: n // 2], perm[n // 2:]
        for li, lam in enumerate(lam_grid):
            try:
                f1 = pathlasso(Xv[h1], M[h1], Yv[h1], lam=lam, **kwargs)
                f2 = pathlasso(Xv[h2], M[h2], Yv[h2], lam=lam, **kwargs)
            except Exception:
                continue
            kappa[li, r] = _cohen_kappa(np.abs(f1["AB"]) > zero_cutoff,
                                        np.abs(f2["AB"]) > zero_cutoff)
    mean_kappa = np.nanmean(kappa, axis=1)
    if np.all(np.isnan(mean_kappa)):
        raise RuntimeError("pathlasso_ksc: every fit failed")
    best = np.nanmax(mean_kappa)
    ok = np.flatnonzero(mean_kappa >= best - vss_cut)
    idx = int(ok[np.argmax(lam_grid[ok])])
    return {"lambda_grid": lam_grid, "kappa": kappa, "mean_kappa": mean_kappa,
            "lambda_index": idx, "lambda_est": float(lam_grid[idx])}
