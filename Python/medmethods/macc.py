"""Multilevel mediation analysis with structured unmeasured confounding.

Single-, two- and three-level mediation models in which the mediator and outcome
errors are correlated; the correlation ``delta`` becomes identifiable from the
multilevel structure and is estimated rather than assumed.

Not yet ported. The estimator maximises a hierarchical / marginal likelihood
over the correlation and the variance components, which the R code delegates to
``nlme::lme`` and ``lme4::lmer`` (with ``optimx`` as the inner optimiser). A
Python port needs an equivalent linear-mixed-model fit -- either a
self-contained REML implementation or a dependency such as statsmodels, which
the CAP project deliberately avoided.

Note that the *single-level* model is a closed-form calculation, but it requires
``delta`` to be supplied: it is not identifiable from one level, and leaving it
at zero biases the b-path toward zero.
"""
from __future__ import annotations

from ._stub import make_stub

__all__ = ["macc"]

macc = make_stub("macc", "multilevel mediation under structured confounding", "macc",
                 "It maximises a hierarchical likelihood over the variance "
                 "components, which R delegates to nlme/lme4.")
