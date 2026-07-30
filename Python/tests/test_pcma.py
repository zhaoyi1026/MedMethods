"""Tests for the principal-component mediation port."""
import numpy as np
from medmethods import examples, pcma, pcma_coef, pcma_loglike


def _cos(a, b):
    return abs(a @ b) / np.sqrt((a @ a) * (b @ b))


def test_coef_at_true_projections():
    d = examples.pcma_example(n=600, seed=2)
    t = d["truth"]
    c = pcma_coef(d["X"], d["M"], d["Y"], t["Phi"][:, 0], t["Psi"][:, 0])
    # alpha[0,0] = 2 and beta[0] = 2 by construction; projections are sign-free
    assert abs(abs(c["alpha"]) - 2.0) < 0.3
    assert abs(abs(c["beta"]) - 2.0) < 0.4


def test_recovers_projections():
    d = examples.pcma_example(n=600, seed=2)
    f = pcma(d["X"], d["M"], d["Y"], nD=2, ninitial=3)
    assert f["Phi"].shape == (5, 2)
    assert f["Psi"].shape == (10, 2)
    assert _cos(f["Phi"][:, 0], d["truth"]["Phi"][:, 0]) > 0.9
    assert _cos(f["Psi"][:, 0], d["truth"]["Psi"][:, 0]) > 0.9


def test_loglike_decreases_from_random_start():
    d = examples.pcma_example(n=300, seed=6)
    rng = np.random.RandomState(0)
    p, q = d["X"].shape[1], d["M"].shape[1]
    phi = rng.normal(size=p); phi /= np.sqrt(phi @ phi)
    psi = rng.normal(size=q); psi /= np.sqrt(psi @ psi)
    c = pcma_coef(d["X"], d["M"], d["Y"], phi, psi)
    start = pcma_loglike(d["X"], d["M"], d["Y"], phi, psi, c["alpha"], c["beta"],
                         c["gamma"], 1.0, 1.0)
    f = pcma(d["X"], d["M"], d["Y"], nD=1, ninitial=3)
    assert f["logLik"][0] < start
