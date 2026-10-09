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


ML_DATA = """clear
set seed 611
set obs 80
gen x1 = round(runiform()*10, .1)
gen x2 = round(runiform()*5, .01)
gen y = runiform() < invlogit(-2 + .3*x1 - .2*x2)
gen c = floor(-ln(runiform())*exp(.1 + .1*x1))
"""


def test_logit_score_is_zero_at_estimate(run):
    run(ML_DATA)
    out = run("logit y x1 x2")
    assert "Iteration 0:   log likelihood =" in out and "Logistic regression" in out
    y, x1, x2 = _arrays(run, "y", "x1", "x2")
    b = run.session.e["b"].data.ravel()
    X = np.column_stack([x1, x2, np.ones_like(x1)])
    p = 1 / (1 + np.exp(-X @ b))
    assert np.allclose(X.T @ (y - p), 0, atol=1e-3)
    H = (X * (p * (1 - p))[:, None]).T @ X
    assert run.eval("_se[x1]") == pytest.approx(np.sqrt(np.linalg.inv(H)[0, 0]), rel=1e-6)
    assert run.eval("e(chi2)") == pytest.approx(2 * (run.eval("e(ll)") - run.eval("e(ll_0)")))


def test_logistic_odds_ratios_and_predict(run):
    run(ML_DATA)
    out = run("logistic y x1 x2")
    assert "Odds Ratio" in out and "Iteration" not in out
    out = run("predict p")
    assert "(option pr assumed; Pr(y))" in out
    (p,) = _arrays(run, "p")
    assert np.all((p > 0) & (p < 1))


def test_poisson_matches_score_equations(run):
    run(ML_DATA)
    run("poisson c x1 x2")
    c, x1, x2 = _arrays(run, "c", "x1", "x2")
    b = run.session.e["b"].data.ravel()
    X = np.column_stack([x1, x2, np.ones_like(x1)])
    assert np.allclose(X.T @ (c - np.exp(X @ b)), 0, atol=1e-3)
    out = run("poisson c x1 x2, irr vce(robust)")
    assert "IRR" in out and "Wald chi2(2)" in out and "log pseudolikelihood" in out


def test_glm_binomial_equals_logit(run):
    run(ML_DATA)
    run("quietly logit y x1 x2")
    b1 = run.session.e["b"].data.ravel().copy()
    run("quietly glm y x1 x2, family(binomial)")
    b2 = run.session.e["b"].data.ravel()
    assert np.allclose(b1, b2, atol=1e-6)


def test_fixed_effects_match_dummies(run):
    run("clear\nset seed 9\nset obs 40\ngen id = ceil(_n/4)\ngen t = mod(_n-1,4)+1\n"
        "gen x = runiform()*10 + id\ngen y = x + id + runiform()\nqui xtset id t")
    run("quietly regress y x i.id")
    bx = run.eval("_b[x]")
    run("quietly xtreg y x, fe")
    assert run.eval("_b[x]") == pytest.approx(bx, rel=1e-10)
    run("quietly areg y x, absorb(id)")
    assert run.eval("_b[x]") == pytest.approx(bx, rel=1e-10)


def test_margins_linear_model(run):
    run(DATA)
    run("quietly regress y x1 x2")
    out = run("margins, dydx(x1)")
    assert "Average marginal effects" in out
    assert run.session.r["b"].data[0, 0] == pytest.approx(run.eval("_b[x1]"), rel=1e-6)


def test_ordered_and_multinomial_run(run):
    run(ML_DATA + "gen r = 1 + (x1 > 3) + (x1 > 7)\n")
    out = run("ologit r x2")
    assert "/cut1" in out and "/cut2" in out
    out = run("mlogit r x2")
    assert "(base outcome)" in out
