"""sem: path models and confirmatory factor analysis by ML, checked against closed forms."""

import numpy as np
import pytest


@pytest.fixture
def semdata(run):
    rng = np.random.default_rng(5)
    N = 300
    F = rng.normal(size=N)
    G = .5 * F + rng.normal(scale=.8, size=N)
    cols = {"x1": F + rng.normal(scale=.6, size=N), "x2": .8 * F + rng.normal(scale=.7, size=N),
            "x3": 1.2 * F + rng.normal(scale=.5, size=N), "y1": G + rng.normal(scale=.6, size=N),
            "y2": .9 * G + rng.normal(scale=.6, size=N), "y3": 1.1 * G + rng.normal(scale=.6, size=N),
            "z": rng.normal(size=N)}
    cols["w"] = .5 * cols["z"] + .3 * cols["x1"] + rng.normal(size=N)
    run("clear\ninput " + " ".join("double " + k for k in cols) + "\n"
        + "\n".join(" ".join(repr(float(cols[k][i])) for k in cols) for i in range(N)) + "\nend")
    run.cols = cols
    return run


def _b(run, eq, name):
    e = run.session.e["b"]
    for j, (q, c) in enumerate(zip(e.coleq, e.colnames)):
        if q == eq and c == name:
            return float(e.data[0, j])
    raise KeyError((eq, name))


def test_path_model_is_ols(semdata):
    run = semdata
    out = run("sem (w <- z x1)")
    assert "Structural equation model" in out and "chi2(0)" in out
    X = np.column_stack([run.cols["z"], run.cols["x1"], np.ones(300)])
    b = np.linalg.lstsq(X, run.cols["w"], rcond=None)[0]
    assert _b(run, "w", "z") == pytest.approx(b[0], rel=1e-6)
    assert _b(run, "w", "x1") == pytest.approx(b[1], rel=1e-6)
    e = run.cols["w"] - X @ b
    assert _b(run, "var(e.w)", "_cons") == pytest.approx(e @ e / 300, rel=1e-5)


def test_cfa_just_identified(semdata):
    run = semdata
    run("sem (F -> x1 x2 x3)")
    X = np.column_stack([run.cols[k] for k in ("x1", "x2", "x3")])
    S = np.cov(X.T, bias=True)
    assert _b(run, "x2", "F") == pytest.approx(S[1, 2] / S[0, 2], rel=1e-5)
    assert _b(run, "var(F)", "_cons") == pytest.approx(S[0, 1] * S[0, 2] / S[1, 2], rel=1e-5)
    assert run.eval("e(df_ms)") == 0


def test_full_model_gof_predict(semdata):
    run = semdata
    out = run("sem (F -> x1 x2 x3) (G -> y1 y2 y3) (G <- F)")
    # endogenous latent variables first, as in Stata 14 (compat 1302)
    assert "( 1)  [y1]G = 1" in out and "( 2)  [x1]F = 1" in out
    assert run.eval("e(df_ms)") == 8
    out = run("estat gof, stats(all)")
    for word in ("RMSEA", "CFI", "TLI", "SRMR", "AIC"):
        assert word in out
    assert 0 <= run.eval("r(cfi)") <= 1
    out = run("sem (F -> x1 x2 x3) (G -> y1 y2 y3) (G <- F), standardized")
    assert "Standardized" in out
    run("predict f, latent(F)")
    f = np.asarray(run.session.data.get("f").data)
    assert np.corrcoef(f, run.cols["x3"])[0, 1] > .8


def test_errors(semdata):
    run = semdata
    run("capture sem (F -> x1 x2)")      # factor with two indicators: not identified
    assert run.eval("_rc") != 0
    run("capture sem w z")
    assert run.eval("_rc") == 198
