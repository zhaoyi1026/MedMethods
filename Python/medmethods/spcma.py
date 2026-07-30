"""Sparse principal component based high-dimensional mediation analysis.

Mediation through the leading (sparse) principal components of a
high-dimensional mediator, with the loadings regularised by a fused lasso so
they are piecewise constant.

Not yet ported. Two pieces have no drop-in Python counterpart:

* the sparse-PCA loadings are fitted with R's ``genlasso`` (fused lasso over a
  penalty matrix / graph), which needs a generalized-lasso path solver; and
* ``mcma_bk`` delegates its bootstrap to the CRAN ``mediation`` package.

The plain principal-component variant (``mcma_pca``) needs neither and is the
easiest next port: PCA of the mediator, then a per-component single-mediator
mediation fit with a bootstrap.
"""
from __future__ import annotations

from ._stub import make_stub

__all__ = ["spcma", "mcma_pca", "mcma_bk"]

spcma = make_stub("spcma", "sparse principal component mediation analysis", "spcma",
                  "Its sparse loadings need a generalized-lasso (fused lasso) path "
                  "solver equivalent to R's genlasso.")
mcma_pca = make_stub("mcma_pca", "principal component mediation analysis", "mcma_pca",
                     "Planned as the next port -- it needs only PCA plus a "
                     "per-component bootstrap.")
mcma_bk = make_stub("mcma_bk", "Baron-Kenny multivariate mediation", "mcma_bk",
                    "R delegates its bootstrap to the CRAN mediation package.")
