"""Shared machinery for the mediation methods.

Mirrors the helpers used across the R implementations: generalized-inverse
solves, soft thresholding, per-subject covariance cubes, projected scores,
weighted accumulation, the generalized smallest-eigenvector solve, and the
deviation-from-diagonality (DfD).

Kept deliberately py3.7 / numpy>=1.16 friendly: no ``np.random.default_rng``
(added in numpy 1.17), so ``RandomState`` is used throughout.
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import eigh

__all__ = [
    "as_matrix_list", "ginv_solve", "soft_threshold", "cov_cube", "scores",
    "accum", "gamma_solve", "diag_level", "standardize", "norm_col",
]


def as_matrix_list(M):
    """Coerce a sequence of array-likes to a list of 2-D float arrays."""
    out = []
    for Mi in M:
        A = np.asarray(Mi, dtype=float)
        out.append(A.reshape(A.shape[0], -1))
    return out


def ginv_solve(A, b):
    """``MASS::ginv(A) %*% b`` -- Moore-Penrose solve, as the R code uses."""
    return np.linalg.pinv(np.asarray(A, dtype=float)) @ np.asarray(b, dtype=float)


def soft_threshold(mu, lam):
    """Elementwise soft threshold ``sign(mu) * max(|mu| - lam, 0)``.

    Matches the R ``soft.thred`` / ``soft_thred``.
    """
    mu = np.asarray(mu, dtype=float)
    lam = np.asarray(lam, dtype=float)
    return np.sign(mu) * np.maximum(np.abs(mu) - lam, 0.0)


def cov_cube(M, ddof=1):
    """Per-subject covariance cube of shape ``(p, p, n)`` plus sample sizes.

    ``ddof=1`` reproduces R's ``cov()`` (denominator ``T_i - 1``), which is what
    the graph-mediator code uses.
    """
    M = as_matrix_list(M)
    n = len(M)
    p = M[0].shape[1]
    cube = np.empty((p, p, n))
    nT = np.empty(n)
    for i, Mi in enumerate(M):
        Ti = Mi.shape[0]
        nT[i] = Ti
        Mc = Mi - Mi.mean(axis=0, keepdims=True)
        cube[:, :, i] = (Mc.T @ Mc) / (Ti - ddof)
    return cube, nT


def scores(Sigma, v):
    """Projected variances ``v' Sigma_i v`` for all ``i`` (Sigma is p x p x n)."""
    v = np.asarray(v, dtype=float).ravel()
    return np.einsum("i,ijk,j->k", v, Sigma, v)


def accum(Sigma, w):
    """Weighted accumulation ``sum_i w_i Sigma_i`` (Sigma is p x p x n)."""
    return Sigma @ np.asarray(w, dtype=float).ravel()


def gamma_solve(A, H):
    """Smallest generalized eigenvector: minimize ``g'Ag`` s.t. ``g'Hg = 1``.

    Mirrors the R ``eigen.solve`` / ``gamma.solve``. ``scipy.linalg.eigh(A, H)``
    returns ascending generalized eigenpairs with ``V'HV = I``.
    """
    A = np.asarray(A, dtype=float)
    H = np.asarray(H, dtype=float)
    A = 0.5 * (A + A.T)
    H = 0.5 * (H + H.T)
    _, vecs = eigh(A, H)
    return vecs[:, 0]


def norm_col(v):
    """Unit-normalise and sign-fix so the largest-magnitude entry is positive."""
    v = np.asarray(v, dtype=float).ravel()
    v = v / np.sqrt(np.sum(v ** 2))
    if v[np.argmax(np.abs(v))] < 0:
        v = -v
    return v


def diag_level(Sigma, Theta, nT):
    """Deviation-from-diagonality across the leading directions.

    Matches the R ``diag.level``: for each ``k``, the size-weighted geometric
    mean over subjects of ``prod(diag(Mk)) / det(Mk)`` where
    ``Mk = Theta[:, :k]' Sigma_i Theta[:, :k]``.
    """
    Theta = np.asarray(Theta, dtype=float)
    if Theta.ndim == 1:
        Theta = Theta[:, None]
    n = Sigma.shape[2]
    r = Theta.shape[1]
    nT = np.asarray(nT, dtype=float)
    dl = np.ones((n, r))
    for i in range(n):
        for j in range(1, r):
            P = Theta[:, : j + 1]
            Mt = P.T @ Sigma[:, :, i] @ P
            dl[i, j] = np.prod(np.diag(Mt)) / np.linalg.det(Mt)
    w = nT / nT.sum()
    avg = np.array([np.prod(dl[:, j] ** w) for j in range(r)])
    return {"avg_level": avg, "sub_level": dl}


def standardize(*arrays):
    """Centre and scale each array by its (n-1) SD, returning values and scales.

    Reproduces the ``standardize = TRUE`` preprocessing in the ADMM methods:
    each column is centred and divided by its sample SD.
    """
    out = []
    for A in arrays:
        A = np.asarray(A, dtype=float)
        one_d = A.ndim == 1
        Am = A.reshape(A.shape[0], -1) if one_d else A
        mu = Am.mean(axis=0)
        sd = Am.std(axis=0, ddof=1)
        sd = np.where(sd == 0, 1.0, sd)
        Z = (Am - mu) / sd
        out.append((Z.ravel() if one_d else Z, mu, sd))
    return out
