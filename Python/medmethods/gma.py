"""Granger mediation analysis of time series.

A mediation model for time series with VAR(p) autoregressive errors and a
correlation parameter capturing unmeasured confounding, in single-level and
two-level (multi-subject) forms.

Not yet ported. The single-level model is closed-form linear algebra -- OLS
projections plus an asymptotic variance built from the VAR companion matrix --
and is the tractable half. The two-level model, where the correlation becomes
identifiable, needs ``nlme::gls`` / ``lme``-style generalized least squares and
a profile likelihood over the correlation.

As with ``macc``, the single-level fit requires ``delta`` to be supplied; the
default of zero biases the b-path toward zero.
"""
from __future__ import annotations

from ._stub import make_stub

__all__ = ["gma"]

gma = make_stub("gma", "Granger mediation analysis", "gma",
                "Its two-level model needs nlme-style generalized least squares "
                "and a profile likelihood over the error correlation.")
