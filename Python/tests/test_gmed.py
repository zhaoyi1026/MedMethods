"""Tests for the graph-mediator port."""
import numpy as np
from medmethods import examples, gmed, gmed_coef, gmed_boot, gmed_refit


def _cos(a, b):
    return abs(a @ b) / np.sqrt((a @ a) * (b @ b))


def test_coef_at_true_theta():
    """The a-path is recovered directly; the b-path is attenuated at finite Ti."""
    d = examples.gmed_example(n=60, p=8, Ti=300, seed=2)
    c = gmed_coef(d["X"], d["M"], d["Y"], d["truth"]["theta"])
    assert abs(c["alpha"][0] - 1.0) < 0.2
    assert 0.4 < c["beta"] < 1.2          # attenuated but the right sign/scale


def test_b_path_attenuation_shrinks_with_Ti():
    """The b-path is attenuated by errors-in-variables and improves with Ti.

    ``log(theta' Sigma_i theta)`` estimates the true log-variance with error of
    order ``sqrt(2/Ti)``, while the generator gives it only ``tau = 0.1`` of
    within-arm signal, so the attenuation shrinks slowly. The reliable property
    is that it moves the right way, not that it vanishes at any given Ti. This is
    a property of the estimator, not of the port.
    """
    betas = []
    for Ti in (150, 1500):
        d = examples.gmed_example(n=60, p=8, Ti=Ti, seed=2)
        betas.append(gmed_coef(d["X"], d["M"], d["Y"], d["truth"]["theta"])["beta"])
    assert betas[0] < betas[1]
    assert betas[1] > 0.6


def test_recovers_direction_and_effects():
    d = examples.gmed_example(n=80, p=8, Ti=400, seed=3)
    f = gmed(d["X"], d["M"], d["Y"], nD=1, ninitial=4)
    assert f["theta"].shape == (8, 1)
    assert _cos(f["theta"][:, 0], d["truth"]["theta"]) > 0.9
    assert abs(f["alpha"][0, 0] - 1.0) < 0.25
    assert f["IE"][0] > 0.3


def test_identity_H_finds_a_different_direction():
    """H must be the average covariance; an identity H collapses onto a
    covariate-free background direction."""
    d = examples.gmed_example(n=60, p=8, Ti=300, seed=4)
    good = gmed(d["X"], d["M"], d["Y"], nD=1, ninitial=4)
    bad = gmed(d["X"], d["M"], d["Y"], H=np.eye(8), nD=1, ninitial=4)
    c_good = _cos(good["theta"][:, 0], d["truth"]["theta"])
    c_bad = _cos(bad["theta"][:, 0], d["truth"]["theta"])
    assert c_good > c_bad


def test_refit_matches_main_fit():
    d = examples.gmed_example(n=60, p=8, Ti=300, seed=5)
    f = gmed(d["X"], d["M"], d["Y"], nD=1, ninitial=3)
    r = gmed_refit(d["X"], d["M"], d["Y"], f["theta"])
    assert np.allclose(r["beta"], f["beta"])
    assert np.allclose(r["IE"], f["IE"])


def test_boot_gives_finite_se():
    d = examples.gmed_example(n=50, p=6, Ti=200, seed=6)
    f = gmed(d["X"], d["M"], d["Y"], nD=1, ninitial=2)
    b = gmed_boot(d["X"], d["M"], d["Y"], theta=f["theta"][:, 0], sims=60)
    assert np.all(np.isfinite(b["SE"]))
    assert b["names"][3] == "IE"
