"""Estatística descritiva: summarize, ci, ttest, prtest, correlate, centile,
pctile/xtile, tabulate e table. Valores conferidos com numpy/scipy."""

import numpy as np
import pytest
from scipy import stats as st


def data(run):
    run("clear\nset obs 20\ngen id = _n\ngen grupo = mod(_n, 2)\n"
        'label define g 0 "Controle" 1 "Tratado"\nlabel values grupo g\n'
        "gen x = id * 1.5 + mod(id, 3)\ngen y = 100 - id * 2 + mod(id * 7, 5)\nreplace x = . in 3")
    ids = np.arange(1, 21)
    x = ids * 1.5 + ids % 3
    x = np.delete(x, 2)
    return ids, x


def test_summarize(run):
    ids, x = data(run)
    out = run("summarize x y")
    assert out.splitlines()[1] == "    Variable |        Obs        Mean    Std. Dev.       Min        Max"
    assert "           x |         19    17.44737    8.693699        2.5         32" in out
    r = run.session.r
    assert r["N"] == 20 and r["mean"] == 81
    run("summarize x, detail")
    r = run.session.r
    assert r["N"] == 19 and abs(r["sd"] - np.std(x, ddof=1)) < 1e-10
    assert abs(r["skewness"] - st.skew(x)) < 1e-10 and abs(r["kurtosis"] - st.kurtosis(x, fisher=False)) < 1e-10
    assert r["p50"] == 18 and r["p25"] == 9.5
    run("summarize x [fw=grupo+1], meanonly")
    w = (ids % 2 + 1)
    xx = ids * 1.5 + ids % 3
    keep = ids != 3
    assert abs(run.session.r["mean"] - np.average(xx[keep], weights=w[keep])) < 1e-10
    assert run.session.r["N"] == w[keep].sum()
    out = run("by grupo, sort: summarize x")
    assert "-> grupo = Controle" in out and "-> grupo = Tratado" in out


def test_ci_ttest_prtest(run):
    ids, x = data(run)
    run("ci x")
    lo, hi = st.t.interval(0.95, len(x) - 1, loc=x.mean(), scale=st.sem(x))
    assert abs(run.session.r["lb"] - lo) < 1e-9 and abs(run.session.r["ub"] - hi) < 1e-9
    out = run("ttest x == 15")
    res = st.ttest_1samp(x, 15)
    assert abs(run.session.r["t"] - res.statistic) < 1e-9 and abs(run.session.r["p"] - res.pvalue) < 1e-9
    assert "One-sample t test" in out and "Pr(|T| > |t|)" in out
    run("ttest x, by(grupo)")
    g = ids % 2
    a = (ids * 1.5 + ids % 3)
    m = ids != 3
    res = st.ttest_ind(a[m & (g == 0)], a[m & (g == 1)])
    assert abs(run.session.r["t"] - res.statistic) < 1e-9 and abs(run.session.r["p"] - res.pvalue) < 1e-9
    run("ttest x, by(grupo) unequal")
    res = st.ttest_ind(a[m & (g == 0)], a[m & (g == 1)], equal_var=False)
    assert abs(run.session.r["p"] - res.pvalue) < 1e-9
    run("ttest x == y")
    y = 100 - ids * 2 + (ids * 7) % 5
    res = st.ttest_rel(a[m], y[m])
    assert abs(run.session.r["t"] - res.statistic) < 1e-9
    run("gen b = x > 15 if x < .\nprtest b == 0.5")
    assert abs(run.session.r["P_1"] - (x > 15).mean()) < 1e-12


def test_correlate_pwcorr(run):
    ids, x = data(run)
    out = run("correlate x y")
    assert out.startswith("(obs=19)")
    a = ids * 1.5 + ids % 3
    y = 100 - ids * 2 + (ids * 7) % 5
    m = ids != 3
    assert abs(run.session.r["rho"] - np.corrcoef(a[m], y[m])[0, 1]) < 1e-12
    out = run("pwcorr x y id, sig obs star(0.05)")
    assert "-0.9929*" in out and "       20" in out


def test_centile_pctile_xtile(run):
    ids, x = data(run)
    run("centile x, centile(50)")
    assert run.session.r["c_1"] == 18
    run("xtile q = x, nq(4)\npctile p = x, nq(4)")
    d = run.session.data
    assert list(d.get("p").data[:3]) == [9.5, 18, 25]
    q = d.get("q").data
    assert sorted(set(q[q < 1e300].tolist())) == [1, 2, 3, 4]
    run("_pctile x, p(10 90)")
    assert run.session.r == {"r1": 5.0, "r2": 29.5}


def test_tabulate(run):
    data(run)
    run("gen c3 = mod(id, 3)")
    out = run("tabulate grupo")
    assert "  Controle |         10       50.00       50.00" in out
    assert "     Total |         20      100.00" in out
    out = run("tabulate c3 grupo, chi2 lrchi2 V exact")
    T = np.array([[3, 3], [3, 4], [4, 3]])
    chi2, p, dof, _ = st.chi2_contingency(T, correction=False)
    assert abs(run.session.r["chi2"] - chi2) < 1e-9 and abs(run.session.r["p"] - p) < 1e-9
    assert "Pearson chi2(2) =   0.2857   Pr = 0.867" in out and "Fisher's exact" in out
    run("tabulate grupo c3 in 1/4, exact")
    out = run("tabulate c3 grupo, row col")
    assert "| Key" in out and "row percentage" in out
    run("tabulate grupo, generate(gg) matcell(F)")
    assert run.session.data.has("gg1") and run.session.matrices["F"].data.ravel().tolist() == [10, 10]
    out = run("table grupo, contents(freq mean x)")
    # colunas de largura única (o maior rótulo ou valor; Stata 14)
    assert " Controle |       10      17.6" in out


def test_fisher_rxc_matches_2x2():
    from scipy.stats import fisher_exact
    from opendta.commands.tabulate import fisher_rxc
    for T in ([[3, 1], [1, 3]], [[5, 0], [1, 4]], [[2, 7], [8, 2]]):
        assert fisher_rxc(np.array(T)) == pytest.approx(fisher_exact(T)[1], abs=1e-9)


def test_tabstat(run):
    data(run)
    out = run("tabstat x y, stats(mean sd n)")
    assert "    mean |  17.44737        81" in out and "       N |        19        20" in out
    out = run("tabstat x, by(grupo) stats(mean sd)")
    assert "Controle |      17.6" in out and "   Total |  17.44737" in out
    out = run("tabstat x y, s(p50 iqr) c(s)")
    assert "variable |       p50       iqr" in out
