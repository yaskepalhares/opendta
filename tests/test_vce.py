"""Robust variance (estimation/vce.py): frequency weights count each observation w times, so a
model fitted with [fweight=w] has the variance of the same model fitted on the expanded data
([U] 20.23 Weighted estimation; [U] 20.22 Obtaining robust variance estimates)."""

import numpy as np
import pytest

from opendta.estimation.vce import bread, sandwich

rng = np.random.default_rng(3)


@pytest.mark.parametrize("clustered", [False, True])
def test_sandwich_fweight_equals_expanded_data(clustered):
    S = rng.normal(size=(30, 3))
    w = rng.integers(1, 4, 30).astype(float)
    g = rng.integers(0, 6, 30)
    V = bread(-np.eye(3) * 5 - 0.1)
    rep = np.repeat(np.arange(30), w.astype(int))
    a, na = sandwich(V, S, w, fweight=True, n=w.sum(), clusters=g if clustered else None)
    b, nb = sandwich(V, S[rep], np.ones(len(rep)), fweight=False, n=len(rep),
                     clusters=g[rep] if clustered else None)
    assert np.allclose(a, b, rtol=1e-12)
    assert na == nb


def test_bread_handles_singular_and_empty():
    assert bread(np.zeros((0, 0))).shape == (0, 0)
    V = bread(-np.array([[1.0, 1.0], [1.0, 1.0]]))
    assert np.allclose(V, V.T)


SETUP = """clear
set seed 5
set obs 120
gen x = runiform()
gen g = mod(_n, 10)
gen t = ceil(-ln(runiform())/exp(-1 + x) * 10)
gen d = runiform() < .8
gen fw = 1 + mod(_n, 3)
gen k = mod(_n, 3)
"""


@pytest.mark.parametrize("weighted,expanded", [
    ("stset t [fweight=fw], failure(d)\nstcox x, vce(cluster g)",
     "expand fw\nstset t, failure(d)\nstcox x, vce(cluster g)"),
    ("stset t [fweight=fw], failure(d)\nstreg x, distribution(weibull) vce(robust)",
     "expand fw\nstset t, failure(d)\nstreg x, distribution(weibull) vce(robust)"),
    ("ologit k x [fweight=fw], vce(cluster g)", "expand fw\nologit k x, vce(cluster g)"),
])
def test_fweight_robust_variance_equals_expanded_data(run, weighted, expanded):
    run(SETUP + "quietly {\n" + weighted + "\n}")
    run("matrix A = e(V)")
    run(SETUP + "quietly {\n" + expanded + "\n}")
    run("matrix B = e(V)\nscalar r = mreldif(A, B)")
    assert run.session.scalars["r"] < 1e-6
