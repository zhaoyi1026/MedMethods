"""medmethods -- causal mediation analysis methods, ported from the R package
`MedMethods <https://github.com/zhaoyi1026/MedMethods>`_.

Implemented and cross-checked against R:

    hetermed      heterogeneous (moderated) mediation effects   [bit-identical]
    pcma          principal component mediation analysis        [verified]
    gmed          mediation with a graph (covariance) mediator  [verified]
    hdmediation   high-dimensional exposures and mediators      [bit-identical]
    pathlasso     Pathway Lasso                                 [bit-identical]

Planned (these raise ``NotImplementedError`` with a pointer to R):

    pathlasso2b   two blocks of high-dimensional mediators
    cfma          causal functional mediation analysis
    spcma         sparse principal component mediation analysis
    macc          multilevel mediation, structured confounding
    gma           Granger mediation analysis

See the package README for verification details and the porting roadmap.
"""
from __future__ import annotations

from . import _core
from . import examples
from . import references as _references

from .hetermed import med_inter, med_inter_inf, med_inter_ite
from .pcma import pcma, pcma_coef, pcma_d1, pcma_inf, pcma_loglike
from .gmed import gmed, gmed_boot, gmed_coef, gmed_d1, gmed_refit
from .hdmediation import hdmediation, hdmediation_pca, hdmediation_std
from .pathlasso import pathlasso, pathlasso_ksc, prox_pair

# not yet ported -- these raise NotImplementedError
from .pathlasso2b import pathlasso2b
from .cfma import cfma_concurrent, cfma_historical
from .spcma import spcma, mcma_pca, mcma_bk
from .macc import macc
from .gma import gma

from .references import REFERENCES, references

# `hetermed` mirrors the R wrapper name for med_inter
hetermed = med_inter
hetermed_inf = med_inter_inf
hetermed_ite = med_inter_ite

__all__ = [
    # implemented
    "med_inter", "med_inter_inf", "med_inter_ite",
    "hetermed", "hetermed_inf", "hetermed_ite",
    "pcma", "pcma_coef", "pcma_d1", "pcma_inf", "pcma_loglike",
    "gmed", "gmed_boot", "gmed_coef", "gmed_d1", "gmed_refit",
    "hdmediation", "hdmediation_pca", "hdmediation_std",
    "pathlasso", "pathlasso_ksc", "prox_pair",
    # planned
    "pathlasso2b", "cfma_concurrent", "cfma_historical",
    "spcma", "mcma_pca", "mcma_bk", "macc", "gma",
    # utilities
    "examples", "references", "REFERENCES", "_core",
]
__version__ = "0.1.0"
