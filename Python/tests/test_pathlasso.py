"""Tests for the Pathway Lasso port."""
import numpy as np
from medmethods import examples, pathlasso, pathlasso_ksc
from medmethods.pathlasso import prox_pair


def test_recovers_signal_and_zeroes_noise():
    d = examples.pathlasso_example(n=200, k=12, seed=2)
    f = pathlasso(d["X"], d["M"], d["Y"], lam=0.001, omega=0, phi=1,
                  max_itr=3000, tol=1e-8)
    AB = f["AB"]
    sig = d["truth"]["signal"]
    assert np.all(np.abs(AB[sig]) > 1.0)
    noise = np.setdiff1d(np.arange(12), sig)
    assert np.max(np.abs(AB[noise])) < np.min(np.abs(AB[sig]))


def test_large_lambda_shrinks_to_zero():
    d = examples.pathlasso_example(n=120, k=10, seed=3)
    f = pathlasso(d["X"], d["M"], d["Y"], lam=1.0, omega=0, phi=1,
                  max_itr=2000, tol=1e-8)
    assert np.max(np.abs(f["AB"])) < 1e-2


def test_prox_pair_zero_lambda_is_soft_threshold():
    a, b = prox_pair(0.0, 2.0, 3.0, 1.5, -0.9, omega1=0.5, omega2=0.4)
    assert abs(a - (1.5 - 0.5) / 2.0) < 1e-12
    assert abs(b - (-0.9 + 0.4) / 3.0) < 1e-12


def test_prox_pair_origin_when_signal_small():
    a, b = prox_pair(1.0, 2.0, 2.0, 0.01, 0.01, omega1=5.0, omega2=5.0)
    assert a == 0.0 and b == 0.0


def test_ksc_returns_a_grid_choice():
    d = examples.pathlasso_example(n=100, k=8, seed=4)
    out = pathlasso_ksc(d["X"], d["M"], d["Y"], lam_grid=[0.001, 0.01, 0.1],
                        n_rep=2, max_itr=800, tol=1e-6)
    assert out["lambda_est"] in (0.001, 0.01, 0.1)
    assert out["kappa"].shape == (3, 2)
