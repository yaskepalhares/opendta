"""Pieces of existing commands completed for 1.0: sort ... in, display _request(), anova
(sequential, nested terms, error terms, repeated measures), the pa estimator of xt commands,
two-way tabulate, summarize(), misstable tree/nested, ci/cii variances, two-sample prtest,
reshape's advanced syntax, the old merge syntax, export excel to .xls and set rng kiss32."""

import numpy as np
import pytest


def test_sort_in_range(run):
    run("clear\ninput x\n5\n4\n3\n2\n1\nend\nsort x in 2/4")
    assert run("list x, noobs clean").split()[1:] == ["5", "2", "3", "4", "1"]
    assert run('display "`: sortedby\'"').strip() == ""


def test_display_request(run):
    answers = iter(["Maria", "42"])
    run.session.input_provider = lambda: next(answers)
    out = run('display "Name? " _request(nome)')
    assert out == "Name? Maria\n"
    assert run.session.macros.globals["nome"] == "Maria"
    run('display "Age? " _request(_idade)\ndisplay `idade\' + 1')
    assert "43" in run.cap.text


def test_display_request_cancelled(run):
    run.session.input_provider = lambda: None
    run('capture display "x" _request(a)')
    assert run.eval("_rc") == 1


def test_merge_old_syntax(run, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    run("clear\ninput id b\n1 10\n2 20\n4 40\nend\nsave u")
    run("clear\ninput id a\n1 1\n2 2\n3 3\nend\nsort id")
    out = run("merge id using u")
    assert "old merge syntax" in out
    assert run("list id a b _merge, noobs clean").split()[4:] == [
        "1", "1", "10", "3", "2", "2", "20", "3", "3", "3", ".", "1", "4", ".", "40", "2"]
    run("clear\ninput id a\n3 3\n1 1\nend")
    run("capture merge id using u")
    assert run.eval("_rc") == 5
    run("merge id using u, sort unique nokeep _merge(m)")
    assert run("count if m == 2").strip() == "0"
    run("clear\ninput z\n7\n8\n9\n9\nend\nmerge using u")       # by position
    assert run("count if _merge == 3").strip() == "3"


def test_reshape_advanced_syntax(run):
    run("clear\ninput id inc80 inc81 ue80 ue81\n1 10 11 0 1\n2 20 21 1 0\nend")
    run("reshape i id\nreshape j year\nreshape xij inc ue")
    out = run("reshape query")
    assert "id" in out and "year" in out and "inc ue" in out
    run("reshape long")
    assert run("count").strip() == "4" and run.session.data.has("year")
    run("reshape wide")
    assert run.session.data.has("inc81")
    run("reshape clear")
    run("capture reshape long")
    assert run.eval("_rc") == 198


def test_reshape_error(run):
    run("clear\ninput id year inc\n1 80 5\n1 80 6\n2 80 7\nend")
    run("capture reshape wide inc, i(id) j(year)")
    assert run.eval("_rc") == 9
    out = run("reshape error")
    assert "2 of 3 observations have repeated year values" in out
    run("clear\ninput id inc80 inc81\n1 1 2\n1 3 4\nend")
    run("capture reshape long inc, i(id) j(year)")
    out = run("reshape error")
    assert "2 of 2 observations have duplicate id values" in out


def test_ci_variances_and_proportions(run):
    run("clear\nset obs 30\ngen x = mod(_n * 7, 11)\ngen b = mod(_n, 3) == 0")
    out = run("ci variances x")
    assert "Variance" in out
    var = run.eval("r(Var)")
    assert var == pytest.approx(np.var([(k * 7) % 11 for k in range(1, 31)], ddof=1))
    from scipy import stats
    assert run.eval("r(lb)") == pytest.approx(29 * var / stats.chi2.ppf(.975, 29))
    run("ci variances x, sd bonett")
    assert run.eval("r(sd)") == pytest.approx(var ** .5)
    assert run.eval("r(lb)") < var ** .5 < run.eval("r(ub)")
    run("cii variances 30 " + repr(var))
    assert run.eval("r(lb)") == pytest.approx(29 * var / stats.chi2.ppf(.975, 29))
    out = run("ci proportions b, wilson")
    assert "Wilson" in out and run.eval("r(proportion)") == pytest.approx(1 / 3)


def test_two_sample_prtest(run):
    run("clear\nset obs 40\ngen g = _n > 15\ngen y = mod(_n, 3) == 0\ngen y2 = mod(_n, 2) == 0")
    out = run("prtest y, by(g)")
    assert "Two-sample test of proportions" in out and "under Ho:" in out
    p1, p2 = 5 / 15, 8 / 25
    pp = 13 / 40
    z = (p1 - p2) / (pp * (1 - pp) * (1 / 15 + 1 / 25)) ** .5
    assert run.eval("r(z)") == pytest.approx(z)
    run("prtest y == y2")
    assert run.eval("r(P_2)") == pytest.approx(.5)
    run("prtesti 15 " + repr(p1) + " 25 " + repr(p2))
    assert run.eval("r(z)") == pytest.approx(z)
    run("prtesti 15 5 25 8, count")
    assert run.eval("r(z)") == pytest.approx(z)


def test_rng_kiss32(run):
    run("clear\nset obs 5\nset rng kiss32\nset seed 42\ngen u = runiform()")
    assert run.eval("c(rng_current)") == "kiss32"
    u = run.session.data.get("u").data.copy()
    assert np.all((u > 0) & (u < 1)) and len(set(u)) == 5
    run("set seed 42\ngen v = runiform()")
    assert np.array_equal(run.session.data.get("v").data, u)
    run("set rng default\nset seed 42\ngen w = runiform()")
    assert run.eval("c(rng_current)") == "mt64"
    assert not np.array_equal(run.session.data.get("w").data, u)


def test_export_excel_xls(run, tmp_path, monkeypatch):
    pytest.importorskip("xlwt")
    pytest.importorskip("xlrd")
    monkeypatch.chdir(tmp_path)
    run('clear\ninput x str5 s\n1 "a"\n2.5 "bb"\nend\nformat x %9.2f')
    run("export excel using t.xls, firstrow(variables)")
    run("clear\nimport excel using t.xls, firstrow")
    assert run("list, noobs clean").split() == ["x", "s", "1", "a", "2.5", "bb"]
    run('clear\ninput y\n9\nend\nexport excel using t.xls, sheet("two") sheetmodify')
    run('clear\nimport excel using t.xls, sheet("Sheet1") firstrow')
    assert run.session.data.nobs == 2


DRUG = """clear
input person drug score
1 1 30
1 2 28
1 3 16
1 4 34
2 1 14
2 2 18
2 3 10
2 4 22
3 1 24
3 2 20
3 3 18
3 4 30
4 1 38
4 2 34
4 3 20
4 4 44
5 1 26
5 2 28
5 3 14
5 4 30
end"""


def test_anova_repeated(run):
    run(DRUG)
    out = run("anova score person drug, repeated(drug)")
    # [R] anova, repeated-measures example (same data and corrections)
    assert "Huynh-Feldt epsilon        =  1.0789" in out
    assert "*Huynh-Feldt epsilon reset to 1.0000" in out
    assert "Greenhouse-Geisser epsilon =  0.6049" in out
    assert "drug |      3    24.76   0.0000   0.0000   0.0006   0.0076" in out
    assert run.eval("e(F_2)") == pytest.approx(232.733333 / 9.4)


def test_anova_sequential_and_error_terms(run):
    import statsmodels.formula.api as smf
    from statsmodels.stats.anova import anova_lm
    import pandas as pd
    rng = np.random.default_rng(3)
    a = np.repeat([1, 2, 3], 7)
    b = np.tile([1, 2], 11)[:21]
    y = a + 2 * b + rng.normal(size=21)
    run("clear\ninput a b double y\n" + "\n".join(f"{p} {q} {float(r)!r}" for p, q, r in zip(a, b, y))
        + "\nend")
    out = run("anova y a b, sequential")
    assert "Seq. SS" in out
    ref = anova_lm(smf.ols("y ~ C(a) + C(b)", pd.DataFrame({"a": a, "b": b, "y": y})).fit(), typ=1)
    assert run.eval("e(ss_1)") == pytest.approx(ref["sum_sq"].iloc[0])
    assert run.eval("e(ss_2)") == pytest.approx(ref["sum_sq"].iloc[1])
    assert run.session.e["sstype"] == "sequential"
    run(DRUG + "\ngen g = person > 3")
    out = run("anova score g / person|g drug")
    ss_g, ss_err = run.eval("e(ss_1)"), run.eval("e(ss_2)")
    assert run.eval("e(F_1)") == pytest.approx((ss_g / 1) / (ss_err / 3))
    assert "person|g |" in out


def test_tabulate_twoway_summarize(run):
    run("clear\nset obs 30\ngen s = mod(_n,2)\ngen g = mod(_n,3)\ngen y = _n")
    out = run("tab s g, summarize(y)")
    assert "Means, Standard Deviations and Frequencies of y" in out
    flat = " ".join(out.split())
    # mean of y where s == 0 and g == 0: _n in 6 12 18 24 30
    assert "18 " in flat


def test_misstable_tree_nested(run):
    run("clear\nset obs 100\ngen a = _n\ngen b = _n\ngen c = _n\nreplace a = . if _n <= 10\n"
        "replace b = . if _n <= 20\nreplace c = . if _n > 95")
    out = run("misstable tree a b c")
    assert "Nested pattern of missing values" in out and "20%" in out
    out = run("misstable nested a b c")
    assert "a(10) -> b(20)" in out and "c(5)" in out


@pytest.fixture
def panel(run):
    rng = np.random.default_rng(1)
    G, T = 40, 5
    rows = []
    for i in range(G):
        u = rng.normal()
        for t in range(T):
            x = rng.normal()
            y = 1 + .5 * x + u + rng.normal()
            yb = int(rng.uniform() < 1 / (1 + np.exp(-(-.3 + .8 * x + u))))
            yc = int(rng.poisson(np.exp(.2 + .3 * x + .5 * u)))
            rows.append(f"{i} {t} {x!r} {y!r} {yb} {yc}")
    run("clear\ninput id t double x double y yb yc\n" + "\n".join(rows) + "\nend\nxtset id t")
    return run


def _sm_gee(run, dep, family):
    sm = pytest.importorskip("statsmodels.api")
    ds = run.session.data
    cols = {n: np.asarray(ds.get(n).data, dtype=float) for n in ("id", "t", "x", dep)}
    X = sm.add_constant(cols["x"])
    m = sm.GEE(cols[dep], X, groups=cols["id"], time=cols["t"], family=family,
               cov_struct=sm.cov_struct.Exchangeable())
    return m.fit(ddof_scale=0)


def test_xtgee_and_pa_match_statsmodels(panel):
    sm = pytest.importorskip("statsmodels.api")
    out = panel("xtreg y x, pa vce(robust)")
    assert "GEE population-averaged model" in out and "exchangeable" in out
    ref = _sm_gee(panel, "y", sm.families.Gaussian())
    b = panel.session.e["b"].data.ravel()
    assert b == pytest.approx(ref.params[[1, 0]], rel=1e-6)
    assert np.sqrt(np.diag(panel.session.e["V"].data)) == pytest.approx(ref.bse[[1, 0]], rel=1e-5)
    panel("xtpoisson yc x, pa vce(robust)")
    ref = _sm_gee(panel, "yc", sm.families.Poisson())
    assert panel.session.e["b"].data.ravel() == pytest.approx(ref.params[[1, 0]], rel=1e-6)
    assert np.sqrt(np.diag(panel.session.e["V"].data)) == pytest.approx(ref.bse[[1, 0]], rel=1e-5)
    panel("xtlogit yb x, pa")
    assert panel.session.e["cmd2"] == "xtlogit" and panel.session.e["family"] == "binomial"
    panel("xtgee y x, corr(ar 1)")
    assert panel.session.e["corr"] == "AR(1)"
    panel("xtgee y x, corr(unstructured)\npredict mu")
    assert panel.session.data.has("mu")


def test_xtreg_be_and_mle(panel):
    sm = pytest.importorskip("statsmodels.api")
    out = panel("xtreg y x, be")
    assert "Between regression (regression on group means)" in out
    ds = panel.session.data
    ids = np.asarray(ds.get("id").data)
    means = {n: np.array([np.asarray(ds.get(n).data)[ids == g].mean() for g in range(40)]) for n in ("x", "y")}
    ref = sm.OLS(means["y"], sm.add_constant(means["x"])).fit()
    assert panel.session.e["b"].data.ravel() == pytest.approx(ref.params[[1, 0]])
    out = panel("xtreg y x, mle")
    assert "Random-effects ML regression" in out and "LR test of sigma_u=0" in out
    smf = pytest.importorskip("statsmodels.formula.api")
    import pandas as pd
    df = pd.DataFrame({n: np.asarray(ds.get(n).data, dtype=float) for n in ("id", "x", "y")})
    md = smf.mixedlm("y ~ x", df, groups=df["id"]).fit(reml=False)
    assert panel.eval("e(ll)") == pytest.approx(md.llf, rel=1e-8)
    assert panel.session.e["sigma_u"] == pytest.approx(np.sqrt(md.cov_re.iloc[0, 0]), rel=1e-4)


def test_xtreg_predict_u_e(panel):
    panel("xtreg y x, fe\npredict u, u\npredict e, e\npredict ue, ue\ngen d = u + e - ue")
    d = np.asarray(panel.session.data.get("d").data)
    assert np.allclose(d, 0, atol=1e-5)
