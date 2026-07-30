"""Multimodal pathway analysis with two blocks of high-dimensional mediators.

    M1 = X beta + e1
    M2 = X zeta + M1 Lambda + e2
    Y  = X delta + M1 theta + M2 pi + e

Not yet ported. The ADMM itself is straightforward -- it has the same shape as
:mod:`medmethods.pathlasso`, reusing the same closed-form pair proximal step --
but the ``Lambda`` block is updated with a lasso fit per column of ``M2``, which
the R code delegates to ``glmnet``. Reproducing it needs a self-contained
coordinate-descent lasso matching glmnet's ``(1/2n)||y - b0 - Xb||^2 + lam|b|_1``
objective. This is the most tractable of the remaining ports.
"""
from __future__ import annotations

from ._stub import make_stub

__all__ = ["pathlasso2b"]

pathlasso2b = make_stub(
    "pathlasso2b", "two blocks of high-dimensional mediators", "pathlasso2b",
    "Its inner Lambda update is a per-column lasso that R delegates to glmnet, "
    "so it needs a self-contained coordinate-descent lasso.")
