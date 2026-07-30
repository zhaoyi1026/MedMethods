"""Tests for the high-dimensional exposures/mediators port."""
import numpy as np
from medmethods import examples, hdmediation, hdmediation_pca
from medmethods.hdmediation import hd_obj


def test_sparsity_increases_with_lambda():
    d = examples.hdmediation_example(n=100, r=8, p=8, seed=2)
    counts = []
    for lam in (0.5, 2.0, 4.0):
        f = hdmediation(d["X"], d["M"], d["Y"], lam=lam, pi=0.5, phi=2, delta=0.5,
                        max_itr=2000, tol=1e-8)
        counts.append(int(np.sum(np.abs(f["IE"]) > 1e-3)))
    assert counts[0] >= counts[1] >= counts[2]


def test_recovers_signal_paths():
    d = examples.hdmediation_example(n=120, r=8, p=8, seed=3)
    f = hdmediation(d["X"], d["M"], d["Y"], lam=0.5, pi=0.5, phi=2, delta=0.5,
                    max_itr=2000, tol=1e-8)
    truth = set(map(tuple, d["truth"]["signal"]))
    sel = set(map(tuple, np.argwhere(np.abs(f["IE"]) > 1e-3)))
    assert len(truth & sel) >= 2


def test_unpenalised_matches_ols():
    """With lambda = 0 the M model reduces to column-wise OLS."""
    d = examples.hdmediation_example(n=200, r=4, p=4, seed=4)
    f = hdmediation(d["X"], d["M"], d["Y"], lam=0.0, standardize=False,
                    max_itr=4000, tol=1e-10)
    ols = np.linalg.lstsq(d["X"], d["M"], rcond=None)[0]
    assert np.max(np.abs(f["alpha"] - ols)) < 1e-3


def test_obj_decomposition():
    d = examples.hdmediation_example(n=60, r=4, p=4, seed=5)
    f = hdmediation(d["X"], d["M"], d["Y"], lam=1.0, pi=0.5, phi=2, delta=0.5,
                    max_itr=1000)
    o = hd_obj(d["X"], d["M"], d["Y"], f["alpha"], f["beta"], f["gamma"],
               1.0, 0.5, 2, 0.5)
    assert abs(o["obj"] - (o["logLik"] + o["R1"] + o["R2"])) < 1e-8


def test_pca_variant_runs():
    d = examples.hdmediation_example(n=60, r=20, p=6, seed=6)
    f = hdmediation_pca(d["X"], d["M"], d["Y"], adaptive=True, var_prop=0.9,
                        lam=1.0, pi=0.5, phi=2, delta=0.5, max_itr=1000)
    assert 1 <= f["n_pc"] <= 20
    assert f["IE"].shape[0] == f["n_pc"]
