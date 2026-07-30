"""Shape / truth-consistency checks for the example generators."""
import numpy as np
from medmethods import examples, REFERENCES


def test_all_generators_have_truth():
    for name in ["hetermed_example", "pcma_example", "gmed_example",
                 "hdmediation_example", "pathlasso_example"]:
        d = getattr(examples, name)()
        assert "truth" in d, name


def test_shapes():
    d = examples.hetermed_example(n=50)
    assert d["X"].shape == (50,) and d["Z"].shape == (50, 3)
    d = examples.pcma_example(n=40, p=5, q=10)
    assert d["X"].shape == (40, 5) and d["M"].shape == (40, 10)
    d = examples.gmed_example(n=12, p=6, Ti=30)
    assert len(d["M"]) == 12 and d["M"][0].shape == (30, 6)
    d = examples.hdmediation_example(n=30, r=6, p=6)
    assert d["X"].shape == (30, 6) and d["truth"]["IE"].shape == (6, 6)
    d = examples.pathlasso_example(n=30, k=9)
    assert d["M"].shape == (30, 9) and d["truth"]["AB"].shape == (9,)


def test_gmed_projection_is_unit_norm():
    d = examples.gmed_example(n=8, p=6, Ti=20)
    th = d["truth"]["theta"]
    assert abs(np.sqrt(th @ th) - 1.0) < 1e-12


def test_references_cover_all_ten_methods():
    assert len(REFERENCES) == 10
    for k, v in REFERENCES.items():
        assert v and v[0]["citation"], k
