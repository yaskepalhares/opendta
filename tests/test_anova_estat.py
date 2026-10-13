"""Fase 8: oneway, anova, tabulate summarize() e estat, conferidos com contas
independentes (scipy e numpy)."""

import math

import numpy as np
import pytest
from scipy import stats as st

DATA = """clear
set seed 802
set obs 90
gen g = mod(_n, 3) + 1
gen h = _n > 40
gen x = rnormal(5, 2)
gen y = 2 + .5*x + g + h*(g == 2) + rnormal()
gen d = y > 6
"""


def _col(run, name):
    return run.session.data.get(name).data.astype(float)


def test_oneway_matches_f_oneway(run):
    run(DATA)
    out = run("oneway y g")
    assert "Analysis of Variance" in out
    assert "Bartlett's test for equal variances" in out
    y, g = _col(run, "y"), _col(run, "g")
    F, p = st.f_oneway(*[y[g == k] for k in (1, 2, 3)])
    assert run.session.r["F"] == pytest.approx(F, rel=1e-10)
    stat, _ = st.bartlett(*[y[g == k] for k in (1, 2, 3)])
    assert run.session.r["chi2bart"] == pytest.approx(stat, rel=1e-10)


def test_tabulate_summarize(run):
    run(DATA)
    out = run("tabulate g, summarize(y)")
    assert "Summary of y" in out
    y, g = _col(run, "y"), _col(run, "g")
    assert f"{y[g == 1].mean():.6f}"[:8] in out.replace(" ", "")


def _partial_ss(y, X_full, X_red):
    def rss(X):
        b, *_ = np.linalg.lstsq(X, y, rcond=None)
        r = y - X @ b
        return r @ r
    return rss(X_red) - rss(X_full)


def test_anova_partial_ss_interaction(run):
    run(DATA)
    run("anova y g##h")
    e = run.session.e
    y, g, h = _col(run, "y"), _col(run, "g"), _col(run, "h")
    # codificação por desvios, independente da implementação
    G = np.column_stack([(g == 1) * 1.0 - (g == 3), (g == 2) * 1.0 - (g == 3)])
    H = ((h == 0) * 1.0 - (h == 1) * 1.0)[:, None]
    GH = G * H
    one = np.ones((len(y), 1))
    full = np.column_stack([one, G, H, GH])
    assert e["ss_1"] == pytest.approx(_partial_ss(y, full, np.column_stack([one, H, GH])), rel=1e-9)
    assert e["ss_3"] == pytest.approx(_partial_ss(y, full, np.column_stack([one, G, H])), rel=1e-9)
    assert e["df_1"] == 2 and e["df_3"] == 2
    assert e["cmd"] == "anova"
    # o modelo é o mesmo do regress com fatores
    b, *_ = np.linalg.lstsq(full, y, rcond=None)
    r = y - full @ b
    assert e["rss"] == pytest.approx(r @ r, rel=1e-9)


def test_anova_continuous(run):
    run(DATA)
    run("anova y g c.x")
    e = run.session.e
    assert e["term_2"] == "x"
    assert e["df_2"] == 1


def test_estat_ic_and_vce(run):
    run(DATA)
    run("quietly regress y x i.g")
    out = run("estat ic")
    assert "Akaike's information criterion" in out
    S = run.session.r["S"]
    ll = run.session.e["ll"]
    k = run.session.e["rank"]
    assert S.data[0, 4] == pytest.approx(-2 * ll + 2 * k)
    assert S.data[0, 5] == pytest.approx(-2 * ll + math.log(90) * k)
    out = run("estat vce, correlation")
    assert "Correlation matrix of coefficients of regress model" in out
    V = run.session.r["V"].data
    assert V[0, 0] == pytest.approx(1.0)


def test_estat_hettest_default_matches_formula(run):
    run(DATA)
    run("quietly regress y x")
    run("estat hettest")
    y, x = _col(run, "y"), _col(run, "x")
    X = np.column_stack([x, np.ones_like(x)])
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    fit = X @ b
    e = y - fit
    u = e * e / (e @ e / len(y))
    Z = np.column_stack([fit, np.ones_like(x)])
    g, *_ = np.linalg.lstsq(Z, u, rcond=None)
    ess = ((Z @ g - u.mean()) ** 2).sum()
    assert run.session.r["chi2"] == pytest.approx(ess / 2, rel=1e-9)


def test_estat_classification_and_gof(run):
    run(DATA)
    run("quietly logit d x")
    out = run("estat classification")
    assert "Correctly classified" in out
    r = run.session.r
    assert 0 <= r["P_corr"] <= 100
    run("estat gof, group(10)")
    assert run.session.r["df"] == run.session.r["m"] - 2


def test_estat_requires_valid_model(run):
    run(DATA)
    run("quietly logit d x")
    run("capture estat hettest")
    assert run.session.rc == 321
