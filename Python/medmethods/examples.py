"""Built-in example-data generators, one per implemented method.

Each returns the inputs in exactly the shape the matching function expects, plus
a ``truth`` dict of the data-generating parameters so estimates can be checked.
Self-contained: numpy only.

These mirror the settings of the R package's ``*_example()`` generators, so the
structure and truth match; the realised data differ because the RNGs differ.
"""
from __future__ import annotations

import numpy as np

__all__ = [
    "hetermed_example", "pcma_example", "gmed_example", "hdmediation_example",
    "pathlasso_example",
]


def _orth(A):
    """Orthonormal basis, sign-fixed so each column's largest entry is positive."""
    Q, _ = np.linalg.qr(A)
    for j in range(Q.shape[1]):
        if Q[np.argmax(np.abs(Q[:, j])), j] < 0:
            Q[:, j] = -Q[:, j]
    return Q


def _rmvnorm(rng, n, mean, sigma):
    ev, V = np.linalg.eigh(sigma)
    R = V @ (V.T * np.sqrt(np.maximum(ev, 0))[:, None])
    return rng.normal(size=(n, sigma.shape[0])) @ R + mean


def hetermed_example(n=600, seed=100):
    """Moderated mediation: both the a-path and b-path vary with covariates."""
    rng = np.random.RandomState(seed)
    alpha0 = np.array([0.2, 0.3, 0.0])
    alpha1 = np.array([0.5, 0.4, 0.0])
    gamma0 = np.array([0.1, 0.0, 0.0])
    gamma1 = np.array([0.3, 0.0, 0.0])
    beta0, beta1 = 0.8, 0.2

    X = rng.choice([-1.0, 1.0], size=n)
    Z = np.column_stack([np.ones(n), rng.normal(size=n),
                         rng.binomial(1, 0.5, size=n).astype(float)])
    M = Z @ alpha0 + X * (Z @ alpha1) + rng.normal(scale=0.5, size=n)
    Y = (Z @ gamma0 + X * (Z @ gamma1) + (beta0 + beta1 * X) * M
         + rng.normal(scale=0.5, size=n))
    Za1 = Z @ alpha1
    return {
        "X": X, "M": M, "Y": Y, "Z": Z,
        "truth": {
            "alpha0": alpha0, "alpha1": alpha1, "beta0": beta0, "beta1": beta1,
            "gamma0": gamma0, "gamma1": gamma1,
            "NIE": 2 * (beta0 + beta1 * X) * Za1,
            "NDE": 2 * (Z @ gamma1) + 2 * beta1 * (Z @ alpha0 + X * Za1),
        },
    }


def pcma_example(n=400, p=5, q=10, seed=100):
    """Orthogonal exposure/mediator projections with parallel mediation."""
    rng = np.random.RandomState(50)
    Phi = _orth(rng.normal(size=(p, p)))
    rng = np.random.RandomState(100)
    Psi = _orth(rng.normal(size=(q, q)))

    sigma2_x = np.sort(np.exp(np.linspace(3, -4, p)))[::-1]
    Sigma_X = Phi @ np.diag(sigma2_x) @ Phi.T

    alpha0 = np.array([2.0, 2.0])
    alpha = np.zeros((p, q))
    alpha[:2, :2] = np.diag(alpha0)
    beta = np.concatenate([[2.0, 1.0], np.zeros(q - 2)])
    gamma = np.concatenate([[1.0, -1.0], np.zeros(p - 2)])

    rng = np.random.RandomState(seed)
    X = _rmvnorm(rng, n, np.zeros(p), Sigma_X)
    Xt = X @ Phi
    eM = rng.normal(size=(n, q))
    Mt = np.zeros((n, q))
    kk = min(p, q)
    for k in range(kk):
        Mt[:, k] = Xt[:, k] * alpha[k, k] + eM[:, k]
    if q > p:
        Mt[:, p:] = eM[:, p:]
    M = Mt @ Psi.T
    Y = Xt @ gamma + Mt @ beta + rng.normal(size=n)
    return {"X": X, "M": M, "Y": Y,
            "truth": {"Phi": Phi, "Psi": Psi, "alpha": alpha, "beta": beta,
                      "gamma": gamma, "IE": alpha @ np.diag(beta), "nD": 2}}


def gmed_example(n=100, p=10, Ti=500, seed=100):
    """Covariance (graph) mediator with a single mediating direction."""
    rng = np.random.RandomState(100)
    Gamma = _orth(rng.uniform(size=(p, p)))
    sig_dim = 1                      # 0-based index of the mediating direction
    alpha0 = alpha = beta = gamma0 = gamma = 1.0
    base = np.linspace(1.5, -3.0, p)
    tau, sigma = 0.1, 0.1

    rng = np.random.RandomState(seed)
    X = rng.binomial(1, 0.5, size=n).astype(float)
    d = np.empty((n, p))
    for j in range(p):
        if j == sig_dim:
            d[:, j] = np.exp(alpha0 + X * alpha + rng.normal(scale=tau, size=n))
        else:
            d[:, j] = np.exp(base[j])
    M = [_rmvnorm(rng, Ti, np.zeros(p), Gamma @ np.diag(d[i]) @ Gamma.T)
         for i in range(n)]
    Y = (gamma0 + X * gamma + np.log(d[:, sig_dim]) * beta
         + rng.normal(scale=sigma, size=n))
    return {"X": X, "M": M, "Y": Y,
            "truth": {"Gamma": Gamma, "theta": Gamma[:, sig_dim], "dim": sig_dim,
                      "alpha": alpha, "beta": beta, "gamma": gamma,
                      "IE": alpha * beta}}


def hdmediation_example(n=100, r=20, p=20, seed=100):
    """High-dimensional exposures and mediators with three signal pathways."""
    nz = 3
    alpha = np.zeros((r, p))
    alpha[np.arange(nz), np.arange(nz)] = [1.5, -1.5, 1.5]
    beta = np.concatenate([[1.5, 1.5, -1.5], np.zeros(p - nz)])
    gamma = np.concatenate([[1.0], np.zeros(r - 1)])

    rng = np.random.RandomState(seed)
    X = rng.normal(size=(n, r))
    M = X @ alpha + rng.normal(scale=0.5, size=(n, p))
    Y = X @ gamma + M @ beta + rng.normal(scale=0.5, size=n)
    IE = alpha @ np.diag(beta)
    return {"X": X, "M": M, "Y": Y,
            "truth": {"alpha": alpha, "beta": beta, "gamma": gamma, "IE": IE,
                      "signal": np.argwhere(np.abs(IE) > 0)}}


def pathlasso_example(n=200, k=50, seed=100):
    """Four signal mediators of mixed sign among ``k``."""
    amp = 2.0
    A = np.concatenate([[amp, -amp, amp, -amp], np.zeros(k - 4)])
    B = np.concatenate([[amp, amp, -amp, -amp], np.zeros(k - 4)])
    C = 1.0
    rng = np.random.RandomState(seed)
    X = rng.binomial(1, 0.5, size=n).astype(float)
    M = np.outer(X, A) + rng.normal(size=(n, k))
    Y = X * C + M @ B + rng.normal(size=n)
    AB = A * B
    return {"X": X, "M": M, "Y": Y,
            "truth": {"A": A, "B": B, "C": C, "AB": AB,
                      "signal": np.flatnonzero(np.abs(AB) > 0)}}
