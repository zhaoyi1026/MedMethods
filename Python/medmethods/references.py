"""Publication references for the mediation methods.

Source of truth: ``setting.md`` in the project root. Each entry gives the
citation, DOI and URL for that method's paper.
"""
from __future__ import annotations

__all__ = ["REFERENCES", "references"]


def _r(citation, doi):
    return {"citation": citation, "doi": doi, "url": "https://doi.org/" + doi}


REFERENCES = {
    "macc": [_r("Zhao, Y., & Luo, X. (2023). Multilevel mediation analysis with "
                "structured unmeasured mediator-outcome confounding. Computational "
                "Statistics & Data Analysis, 179, 107623.",
                "10.1016/j.csda.2022.107623")],
    "gma": [_r("Zhao, Y., & Luo, X. (2019). Granger mediation analysis of multiple "
               "time series with an application to functional magnetic resonance "
               "imaging. Biometrics, 75(3), 788-798.", "10.1111/biom.13056")],
    "spcma": [_r("Zhao, Y., Lindquist, M. A., & Caffo, B. S. (2020). Sparse principal "
                 "component based high-dimensional mediation analysis. Computational "
                 "Statistics & Data Analysis, 142, 106835.",
                 "10.1016/j.csda.2019.106835")],
    "pathlasso": [{"citation": ("Zhao, Y., & Luo, X. (2022). Pathway Lasso: pathway "
                                "estimation and selection with high-dimensional "
                                "mediators. Statistics and Its Interface, 15(1), "
                                "39-50."),
                   "doi": None, "url": None}],
    "pathlasso2b": [_r("Zhao, Y., Li, L., & Caffo, B. S. (2021). Multimodal "
                       "neuroimaging data integration and pathway analysis. "
                       "Biometrics, 77(3), 879-889.", "10.1111/biom.13351")],
    "hdmediation": [_r("Zhao, Y., Li, L., & Alzheimer's Disease Neuroimaging "
                       "Initiative (2022). Multimodal data integration via mediation "
                       "analysis with high-dimensional exposures and mediators. "
                       "Human Brain Mapping, 43(8), 2519-2533.",
                       "10.1002/hbm.25800")],
    "pcma": [_r("Zhao, Y. (2024). Mediation analysis with multiple exposures and "
                "multiple mediators. Statistics in Medicine, 43(25), 4887-4898.",
                "10.1002/sim.10215")],
    "cfma": [_r("Zhao, Y., Luo, X., Sobel, M. E., Lindquist, M. A., & Caffo, B. S. "
                "(2025). Causal functional mediation analysis with an application to "
                "functional magnetic resonance imaging data. Biostatistics, 26(1), "
                "kxaf019.", "10.1093/biostatistics/kxaf019")],
    "gmed": [_r("Xu, Y., & Zhao, Y. (2025). Mediation analysis with graph mediator. "
                "Biostatistics, 26(1), kxaf004.", "10.1093/biostatistics/kxaf004")],
    "hetermed": [_r("Zhao, Y., Li, C., & Tu, W. (2025). Estimation of heterogeneous "
                    "causal mediation effects in a hypertension treatment trial. "
                    "arXiv preprint arXiv:2512.12043.",
                    "10.48550/arXiv.2512.12043")],
}


def references(method=None):
    """Print the reference(s) for ``method``, or for every method."""
    keys = sorted(REFERENCES) if method is None else [method]
    for k in keys:
        if k not in REFERENCES:
            raise KeyError("unknown method %r; known: %s"
                           % (k, ", ".join(sorted(REFERENCES))))
        print(k)
        for ref in REFERENCES[k]:
            print("  " + ref["citation"])
            if ref["url"]:
                print("  " + ref["url"])
