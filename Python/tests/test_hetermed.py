"""Tests for the heterogeneous-mediation port."""
import numpy as np
import pytest

from medmethods import examples, hetermed, hetermed_inf, hetermed_ite


def test_recovers_truth():
    d = examples.hetermed_example(n=800, seed=3)
    f = hetermed(d["X"], d["M"], d["Y"], d["Z"])
    assert abs(f["beta0"] - 0.8) < 0.1
    assert abs(f["beta1"] - 0.2) < 0.1
    assert np.max(np.abs(f["alpha1"] - d["truth"]["alpha1"])) < 0.15


def test_ite_matches_closed_form():
    d = examples.hetermed_example(n=300, seed=4)
    f = hetermed(d["X"], d["M"], d["Y"], d["Z"])
    ite = hetermed_ite(d["X"], d["Z"], f["alpha0"], f["alpha1"], f["beta0"],
                       f["beta1"], f["gamma0"], f["gamma1"])
    assert np.allclose(ite["NIE"], f["ITE"]["NIE"])
    assert np.corrcoef(ite["NIE"], d["truth"]["NIE"])[0, 1] > 0.9


def test_se_scales_as_root_n():
    """SE must scale as n^(-1/2).

    R's coefficient tables applied sqrt() to an already-computed standard error,
    which gives n^(-1/4); this port reports the correct value.
    """
    se = {}
    for n in (500, 2000):
        d = examples.hetermed_example(n=n, seed=7)
        f = hetermed(d["X"], d["M"], d["Y"], d["Z"])
        se[n] = hetermed_inf(d["X"], d["M"], d["Y"], d["Z"], f)["alpha1"]["SE"]
    ratio = se[500] / se[2000]
    assert np.all(ratio > 1.7) and np.all(ratio < 2.3), ratio


def test_inference_tables_are_consistent():
    d = examples.hetermed_example(n=400, seed=8)
    f = hetermed(d["X"], d["M"], d["Y"], d["Z"])
    inf = hetermed_inf(d["X"], d["M"], d["Y"], d["Z"], f)
    for key in ("alpha0", "alpha1", "gamma0", "gamma1", "beta0", "beta1"):
        t = inf[key]
        assert np.allclose(t["zvalue"], t["estimate"] / t["SE"])
        assert np.all(t["LB"] < t["estimate"]) and np.all(t["estimate"] < t["UB"])
    # the moderated a-path is real and should be detected at n = 400
    assert inf["alpha1"]["pvalue"][0] < 0.01


def test_requires_pm1_coding():
    d = examples.hetermed_example(n=100, seed=5)
    with pytest.raises(ValueError):
        hetermed(np.zeros_like(d["X"]), d["M"], d["Y"], d["Z"])


def test_genlasso_not_ported():
    d = examples.hetermed_example(n=100, seed=5)
    with pytest.raises(NotImplementedError):
        hetermed(d["X"], d["M"], d["Y"], d["Z"], method="genlasso")
