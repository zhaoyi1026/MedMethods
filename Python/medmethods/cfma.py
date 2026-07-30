"""Causal functional mediation analysis.

Functional treatment, mediator and outcome; a concurrent model

    M(t) = Z(t) alpha(t) + e1(t),   Y(t) = Z(t) gamma(t) + M(t) beta(t) + e2(t)

and a historical-influence model in which the effect at time t integrates over
the past.

Not yet ported. The estimation is deterministic penalised least squares in a
Fourier basis, so it is a mechanical port, but it needs the basis machinery
(``fourier.basis``, the second-derivative penalty ``Ld2.fourier``, and the
concurrent / historical design-matrix builders) plus the cross-validation and
bootstrap wrappers.
"""
from __future__ import annotations

from ._stub import make_stub

__all__ = ["cfma_concurrent", "cfma_historical"]

_why = ("It needs the Fourier-basis and design-matrix machinery (basis, "
        "second-derivative penalty, concurrent/historical designs).")
cfma_concurrent = make_stub("cfma_concurrent", "concurrent functional mediation",
                            "cfma_concurrent", _why)
cfma_historical = make_stub("cfma_historical", "historical functional mediation",
                            "cfma_historical", _why)
