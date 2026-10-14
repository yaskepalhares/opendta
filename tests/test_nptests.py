"""Fase 8: testes não paramétricos, de normalidade e comandos imediatos."""

import numpy as np
import pytest
from scipy import stats as st

from opendta.commands.nptests import _royston, _sk_z

DATA = """clear
set seed 801
set obs 60
gen y = round(rnormal(10, 3), .5)
gen g = mod(_n, 3)
gen b = _n > 25
gen x = y + rnormal()
"""


def _col(run, name):
    return run.session.data.get(name).data.astype(float)


def test_ranksum_matches_scipy(run):
    run(DATA)
    out = run("ranksum y, by(b)")
    assert "Two-sample Wilcoxon rank-sum (Mann-Whitney) test" in out
    y, b = _col(run, "y"), _col(run, "b")
    res = st.mannwhitneyu(y[b == 0], y[b == 1], use_continuity=False, method="asymptotic")
    z = run.session.r["z"]
    assert 2 * st.norm.sf(abs(z)) == pytest.approx(res.pvalue, rel=1e-9)
    assert run.session.r["N_1"] == 25


def test_kwallis_matches_scipy(run):
    run(DATA)
    run("kwallis y, by(g)")
    y, g = _col(run, "y"), _col(run, "g")
    H = st.kruskal(*[y[g == k] for k in range(3)]).statistic
    assert run.session.r["chi2_adj"] == pytest.approx(H, rel=1e-10)
    assert run.session.r["df"] == 2


def test_spearman_matches_scipy(run):
    run(DATA)
    out = run("spearman y x")
    r = st.spearmanr(_col(run, "y"), _col(run, "x"))
    assert run.session.r["rho"] == pytest.approx(r.statistic, rel=1e-12)
    assert run.session.r["p"] == pytest.approx(r.pvalue, rel=1e-8)
    assert "Test of Ho: y and x are independent" in out


def test_signrank_small_example(run):
    run("clear\ninput d\n1\n-2\n3\n0\n4\n-5\n6\nend")
    run("signrank d = 0")
    r = run.session.r
    # |d| = 1 2 3 0 4 5 6 -> postos 2 3 4 1 5 6 7; positivos: 2+4+5+7
    assert r["sum_pos"] == 18 and r["sum_neg"] == 9 and r["N_tie"] == 1
    E = (28 - 1) / 2
    var = 7 * 8 * 15 / 24 - 1 * 2 * 3 / 24
    assert r["z"] == pytest.approx((18 - E) / np.sqrt(var))


def test_swilk_royston_manual_example():
    # exemplo do manual [R] swilk (auto.dta, mpg): W = 0.94821, n = 74
    W, V, z, p = _royston(0.94821, 74)
    assert round(V, 3) == pytest.approx(3.335, abs=0.002)
    assert round(z, 3) == pytest.approx(2.627, abs=0.002)
    assert p == pytest.approx(0.0043, abs=0.0001)


def test_sktest_components_match_scipy():
    rng = np.random.default_rng(3)
    x = rng.gamma(2.0, size=80)
    z1, z2 = _sk_z(x)
    assert z1 == pytest.approx(st.skewtest(x).statistic, rel=1e-10)
    assert z2 == pytest.approx(st.kurtosistest(x).statistic, rel=1e-10)


def test_sdtest_by(run):
    run(DATA)
    out = run("sdtest y, by(b)")
    assert "Variance ratio test" in out and "Ho: ratio = 1" in out
    y, b = _col(run, "y"), _col(run, "b")
    f = y[b == 0].var(ddof=1) / y[b == 1].var(ddof=1)
    assert run.session.r["F"] == pytest.approx(f, rel=1e-12)


def test_immediate_commands(run):
    out = run("cii means 166 19509 4379")
    assert "18837.93    20180.07" in out
    run("ttesti 24 62.6 15.8 75")
    assert run.session.r["t"] == pytest.approx((62.6 - 75) / (15.8 / np.sqrt(24)))
    run("ttesti 20 20 5 32 15 4")
    assert run.session.r["mu_1"] == pytest.approx(20) and run.session.r["sd_2"] == pytest.approx(4)
    run("prtesti 50 .52 .7")
    assert run.session.r["z"] == pytest.approx((0.52 - 0.7) / np.sqrt(0.7 * 0.3 / 50))
    run("clear\nset obs 3\ngen keep = 1")
    out = run("tabi 30 18 \\ 38 14, chi2")
    assert "Pearson chi2(1)" in out
    assert run.session.data.nobs == 3          # os dados em memória voltam


def test_two_sample_ttest_has_no_full_rule_before_diff(run):
    out = run("ttesti 20 20 5 32 15 4")
    lines = out.splitlines()
    k = next(i for i, ln in enumerate(lines) if ln.startswith("combined"))
    assert lines[k + 1].startswith("---------+") and lines[k + 2].startswith("    diff")
