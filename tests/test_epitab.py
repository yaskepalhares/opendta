"""Epidemiology commands: epitab tables, odds and Mantel-Haenszel, symmetry,
standardization, kappa and ROC. Reference values come from closed formulas,
from SciPy and from the published examples quoted in the comments."""

import math

import numpy as np
import pytest
from scipy import stats as st
from scipy.stats.contingency import odds_ratio

from opendta.commands.kappa import cohen, fleiss
from opendta.core import epi as E

HANLEY = """clear
input disease rating pop
0 1 33
0 2 6
0 3 6
0 4 11
0 5 2
1 1 3
1 2 2
1 3 2
1 4 11
1 5 33
end
"""

STRATA = """clear
input case exp sex pop
1 1 0 10
1 0 0 20
0 1 0 90
0 0 0 80
1 1 1 25
1 0 1 15
0 1 1 75
0 0 1 85
end
"""


# ---------------------------------------------------------------------------
# 2x2 formulas
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("t", [(10, 20, 90, 80), (3, 1, 7, 12), (25, 15, 75, 85)])
def test_exact_or_limits_match_scipy(t):
    a, b, c, d = t
    ref = odds_ratio(np.array([[a, b], [c, d]]), kind="conditional").confidence_interval(0.95)
    lo, hi = E.or_ci_exact(a, b, c, d, 95)
    assert lo == pytest.approx(ref.low, rel=1e-6)
    assert hi == pytest.approx(ref.high, rel=1e-6)
    assert E.fisher_2x2(a, b, c, d)[1] == pytest.approx(st.fisher_exact([[a, b], [c, d]])[1], rel=1e-9)


def test_mh_single_stratum_is_woolf():
    est, lo, hi = E.mh_or([(10, 20, 90, 80)], 95)
    assert est == pytest.approx(10 * 80 / (20 * 90))
    assert (lo, hi) == pytest.approx(E.or_ci_woolf(10, 20, 90, 80, 95))


# ---------------------------------------------------------------------------
# cs, cc, ir, mcc
# ---------------------------------------------------------------------------

def test_csi_measures(run):
    out = run("csi 7 12 9 2")
    assert "Risk ratio" in out and "Prev. frac. ex." in out
    r = run.session.r
    r1, r0 = 7 / 16, 12 / 14
    assert r["rd"] == pytest.approx(r1 - r0)
    assert r["rr"] == pytest.approx(r1 / r0)
    se = math.sqrt(1 / 7 - 1 / 16 + 1 / 12 - 1 / 14)
    assert r["lb_rr"] == pytest.approx(r1 / r0 * math.exp(-st.norm.ppf(.975) * se))
    chi2 = st.chi2_contingency([[7, 12], [9, 2]], correction=False)[0]
    assert r["chi2"] == pytest.approx(chi2)


def test_cci_exact_and_fisher(run):
    out = run("cci 4 386 4 1250, exact")
    assert "(exact)" in out and "2-sided Fisher's exact P" in out
    r = run.session.r
    assert r["or"] == pytest.approx(4 * 1250 / (386 * 4))
    assert r["p_exact"] == pytest.approx(st.fisher_exact([[4, 386], [4, 1250]])[1], rel=1e-9)


def test_cc_by_mantel_haenszel(run):
    run(STRATA)
    out = run("cc case exp [fw=pop], by(sex)")
    assert "M-H combined" in out and "Test that combined OR = 1" in out
    r = run.session.r
    # stratum sums of ad/n and bc/n
    num = 10 * 80 / 200 + 25 * 85 / 200
    den = 20 * 90 / 200 + 15 * 75 / 200
    assert r["or"] == pytest.approx(num / den)
    assert r["crude"] == pytest.approx(1.0)


def test_cs_by_mh_risk_ratio(run):
    run(STRATA)
    run("cs case exp [fw=pop], by(sex)")
    num = 10 * 100 / 200 + 25 * 100 / 200
    den = 20 * 100 / 200 + 15 * 100 / 200
    assert run.session.r["rr"] == pytest.approx(num / den)


def test_iri_exact_ratio(run):
    out = run("iri 41 15 28010 19017")
    assert "(midp)" in out
    r = run.session.r
    assert r["irr"] == pytest.approx((41 / 28010) / (15 / 19017))
    # the exact limits come from the binomial limits of 41/56
    lo_p = st.beta.ppf(.025, 41, 16)
    assert r["lb_irr"] == pytest.approx(lo_p / (1 - lo_p) * 19017 / 28010, rel=1e-9)


def test_mcci_matches_manual_difference(run):
    out = run("mcci 8 8 3 8")
    assert "McNemar's chi2(1) =      2.27" in out
    r = run.session.r
    # [ST] epitab, mcc example: difference .1851852 with limits -.0822542 and .4526246
    assert r["D_f"] == pytest.approx(.1851852, abs=1e-7)
    assert r["lb_D_f"] == pytest.approx(-.0822542, abs=1e-7)
    assert r["ub_D_f"] == pytest.approx(.4526246, abs=1e-7)
    assert r["p_exact"] == pytest.approx(2 * st.binom.cdf(3, 11, .5))


# ---------------------------------------------------------------------------
# tabodds, mhodds, symmetry
# ---------------------------------------------------------------------------

ALC = """clear
input case alc age pop
1 1 1 2
0 1 1 60
1 2 1 4
0 2 1 40
1 3 1 8
0 3 1 20
1 1 2 10
0 1 2 50
1 2 2 15
0 2 2 30
1 3 2 25
0 3 2 15
end
"""


def test_tabodds_homogeneity_is_pearson(run):
    run(ALC)
    run("tabodds case alc [fw=pop]")
    chi2 = st.chi2_contingency([[12, 19, 33], [110, 70, 35]], correction=False)[0]
    assert run.session.r["chi2_hom"] == pytest.approx(chi2)
    assert run.session.r["df_hom"] == 2


def test_mhodds_compare_equals_tabodds_or(run):
    run(ALC)
    run("mhodds case alc age [fw=pop], compare(3,1)")
    a = run.session.r["or"]
    num = 8 * 60 / 90 + 25 * 50 / 100
    den = 2 * 20 / 90 + 10 * 15 / 100
    assert a == pytest.approx(num / den)


def test_symmetry_bowker_and_stuart_maxwell(run):
    run("""clear
input r c n
1 1 20
1 2 5
1 3 2
2 1 3
2 2 30
2 3 6
3 1 1
3 2 4
3 3 25
end""")
    run("symmetry r c [fw=n]")
    r = run.session.r
    assert r["chi2"] == pytest.approx(4 / 8 + 1 / 3 + 4 / 10)
    assert r["df"] == 3
    T = np.array([[20, 5, 2], [3, 30, 6], [1, 4, 25]], float)
    d = (T.sum(1) - T.sum(0))[:2]
    V = np.array([[T[0].sum() + T[:, 0].sum() - 2 * T[0, 0], -(T[0, 1] + T[1, 0])],
                  [-(T[0, 1] + T[1, 0]), T[1].sum() + T[:, 1].sum() - 2 * T[1, 1]]])
    assert r["chi2_sm"] == pytest.approx(float(d @ np.linalg.solve(V, d)))


# ---------------------------------------------------------------------------
# dstdize, istdize
# ---------------------------------------------------------------------------

POPS = """clear
input str6 state agecat pop deaths
Maine 1 100000 20
Maine 2 80000 60
Maine 3 40000 200
Iowa 1 150000 25
Iowa 2 90000 70
Iowa 3 30000 160
end
"""


def test_dstdize_adjusted_rate(run):
    run(POPS)
    run("dstdize deaths pop agecat, by(state)")
    std = np.array([250000, 170000, 70000], float)
    P = std / std.sum()
    maine = np.array([20 / 100000, 60 / 80000, 200 / 40000])
    adj = run.session.r["adj"].data[0]
    assert adj[1] == pytest.approx(float((maine * P).sum()))      # groups in sort order: Iowa, Maine


def test_istdize_smr(run, tmp_path):
    std = tmp_path / "std.dta"
    run(f"""clear
input agecat stdpop stddeaths
1 300000 50
2 200000 140
3 80000 400
end
save "{std}", replace""")
    run(POPS)
    run(f'istdize deaths pop agecat using "{std}", popvars(stddeaths stdpop) by(state)')
    R = np.array([50 / 300000, 140 / 200000, 400 / 80000])
    exp_iowa = float((np.array([150000, 90000, 30000]) * R).sum())
    assert run.session.r["smr"].data[0][0] == pytest.approx(255 / exp_iowa)


# ---------------------------------------------------------------------------
# kappa
# ---------------------------------------------------------------------------

def test_cohen_kappa_and_null_se():
    po, pe, k, se = cohen(np.array([[20, 5], [10, 15]], float), np.eye(2))
    assert (po, pe, k) == pytest.approx((.7, .5, .4))
    # Fleiss, Cohen and Everitt: sqrt(pe + pe^2 - sum r_i c_i (r_i + c_i)) / ((1 - pe) sqrt(n))
    assert se == pytest.approx(math.sqrt(.5 + .25 - .51) / (.5 * math.sqrt(50)))


def test_fleiss_kappa_classic_example():
    # Fleiss (1971): 10 subjects, 14 raters, 5 categories -> kappa = 0.210
    F = np.array([[0, 0, 0, 0, 14], [0, 2, 6, 4, 2], [0, 0, 3, 5, 6], [0, 3, 9, 2, 0], [2, 2, 8, 1, 1],
                  [7, 7, 0, 0, 0], [3, 2, 6, 3, 0], [2, 5, 3, 2, 2], [6, 5, 2, 1, 0], [0, 2, 2, 3, 7]], float)
    assert fleiss(F)[2] == pytest.approx(0.20993, abs=1e-5)


def test_kap_command(run):
    run("""clear
input a b n
1 1 20
1 2 5
2 1 10
2 2 15
end""")
    out = run("kap a b [fw=n]")
    assert "Kappa" in out
    assert run.session.r["kappa"] == pytest.approx(.4)


# ---------------------------------------------------------------------------
# ROC
# ---------------------------------------------------------------------------

def test_roctab_hanley_mcneil_manual_values(run):
    run(HANLEY)
    out = run("roctab disease rating [fw=pop]")
    # [R] roc, Hanley and McNeil (1982) data: 109 0.8932 0.0307 0.83295 0.95339
    assert "109     0.8932       0.0307        0.83295     0.95339" in out
    assert run.session.r["area"] == pytest.approx(0.89317, abs=1e-5)


def test_roccomp_two_classifiers(run):
    run(HANLEY + "expand pop\nset seed 3\ngen r2 = rating + (disease==1)*(runiform()<.3)")
    run("roccomp disease rating r2")
    r = run.session.r
    assert r["df"] == 1 and r["chi2"] > 0
    assert r["V"].data.shape == (2, 2)


def test_rocfit_binormal(run):
    run(HANLEY)
    run("rocfit disease rating [fw=pop], nolog")
    e = run.session.e
    # [R] rocfit example: log likelihood -123.64855, goodness of fit chi2(2) = 0.21
    assert e["ll"] == pytest.approx(-123.64855, abs=1e-4)
    assert e["chi2_gf"] == pytest.approx(0.21, abs=0.005)


def test_lroc_equals_roctab_of_predictions(run):
    run("""clear
set seed 1402
set obs 200
gen x = rnormal()
gen y = runiform() < invlogit(-0.5 + 1.2*x)
quietly logit y x
quietly predict p
""")
    run("lroc, nograph")
    a = run.session.r["area"]
    run("roctab y p")
    assert a == pytest.approx(run.session.r["area"], rel=1e-12)
