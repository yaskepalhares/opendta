"""Fase 5a: regress, variáveis fatoriais, séries temporais e pós-estimação."""

import numpy as np
import pytest

DATA = """clear
set seed 2024
set obs 60
gen g = ceil(_n/20)
gen h = mod(_n, 2)
gen x1 = round(runiform()*10, .1)
gen x2 = round(runiform()*5, .01)
gen y = 3 + 2*x1 - 1.5*x2 + g + round(runiform()*4, .01)
"""


def _arrays(run, *names):
    ds = run.session.data
    return [ds.get(n).data.astype(float) for n in names]


def test_regress_matches_least_squares(run):
    run(DATA)
    out = run("regress y x1 x2")
    y, x1, x2 = _arrays(run, "y", "x1", "x2")
    X = np.column_stack([x1, x2, np.ones_like(x1)])
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    assert run.eval("_b[x1]") == pytest.approx(b[0], rel=1e-10)
    assert run.eval("_b[_cons]") == pytest.approx(b[2], rel=1e-10)
    e = y - X @ b
    s2 = e @ e / (60 - 3)
    se = np.sqrt(np.diag(s2 * np.linalg.inv(X.T @ X)))
    assert run.eval("_se[x2]") == pytest.approx(se[1], rel=1e-9)
    assert "Number of obs   =        60" in out
    assert "      Source |       SS           df       MS      Number of obs" in out
    assert run.eval("e(N)") == 60 and run.eval("e(df_r)") == 57


def test_robust_and_cluster(run):
    run(DATA)
    run("regress y x1 x2, vce(robust)")
    y, x1, x2, g = _arrays(run, "y", "x1", "x2", "g")
    X = np.column_stack([x1, x2, np.ones_like(x1)])
    A = np.linalg.inv(X.T @ X)
    b = A @ X.T @ y
    e = y - X @ b
    V = 60 / 57 * A @ (X * (e ** 2)[:, None]).T @ X @ A
    assert run.eval("_se[x1]") == pytest.approx(np.sqrt(V[0, 0]), rel=1e-9)
    out = run("regress y x1 x2, vce(cluster g)")
    assert "(Std. Err. adjusted for 3 clusters in g)" in out
    assert run.eval("e(df_r)") == 2


def test_factor_variables_names(run):
    run(DATA)
    run("regress y i.g##i.h x1")
    b = run.session.e["b"]
    assert b.colnames == ["1b.g", "2.g", "3.g", "0b.h", "1.h", "1b.g#0b.h", "1b.g#1o.h",
                          "2o.g#0b.h", "2.g#1.h", "3o.g#0b.h", "3.g#1.h", "x1", "_cons"]
    out = run("regress y ib2.g x1")
    assert "          1  |" in out and "          3  |" in out


def test_collinearity_note(run):
    run(DATA + "gen x3 = 2*x1\n")
    out = run("regress y x1 x3 x2")
    # como o Stata (compat 0501): com x3 = 2*x1, sai x1
    assert "note: x1 omitted because of collinearity" in out
    assert "x1 |          0  (omitted)" in out


def test_tsset_and_lags(run):
    out = run("clear\nset obs 6\ngen t = _n\ngen y = t^2\ntsset t\ngen ly = L.y\ngen dy = D.y")
    assert "time variable:  t, 1 to 6" in out
    ly, dy = _arrays(run, "ly", "dy")
    assert ly[0] >= 8.98e307 and ly[3] == 9
    assert dy[2] == 9 - 4
    run("xtset, clear")


def test_lags_respect_panels_and_gaps(run):
    run("clear\nset obs 8\ngen id = ceil(_n/4)\ngen t = mod(_n-1,4)+1\nreplace t = 9 in 4\n"
        "gen y = _n\nxtset id t\ngen ly = L.y")
    (ly,) = _arrays(run, "ly")
    assert ly[0] >= 8.98e307 and ly[1] == 1 and ly[3] >= 8.98e307 and ly[4] >= 8.98e307


def test_test_and_lincom(run):
    run(DATA)
    run("regress y x1 x2")
    t = run.eval("_b[x1]/_se[x1]")
    run("test x1")
    assert run.eval("r(F)") == pytest.approx(t * t, rel=1e-9)
    out = run("lincom x1 + x2")
    assert " ( 1)  x1 + x2 = 0" in out
    assert run.eval("r(estimate)") == pytest.approx(run.eval("_b[x1] + _b[x2]"), rel=1e-12)


def test_predict_and_esample(run):
    run(DATA + "replace x2 = . in 5\n")
    run("regress y x1 x2 if g < 3")
    out = run("predict p")
    assert "(option xb assumed; fitted values)" in out
    run("predict r, residuals")
    assert run.eval("e(N)") == 39
    out = run("count if e(sample)")
    assert out.strip() == "39"
    run("sort x1")
    assert run("count if e(sample)").strip() == "39"
    p, r, y = _arrays(run, "p", "r", "y")
    ok = r < 8.98e307
    assert np.allclose(p[ok] + r[ok], y[ok])


def test_estimates_store_restore(run):
    run(DATA)
    run("regress y x1\nestimates store a\nregress y x2")
    assert run.eval("e(df_m)") == 1
    run("estimates restore a")
    assert run.eval("_b[x1]") > 1
    out = run("estimates table a")
    assert "Variable |" in out
