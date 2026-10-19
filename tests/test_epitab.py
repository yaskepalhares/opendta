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
    # population attributable fraction of the crude table (crude OR = 1)
    assert r["afp"] == pytest.approx(0.0, abs=1e-12)


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


DOLL = """clear
input agecat smokes deaths pyears
1 1 32 52407
2 1 104 43248
3 1 206 28612
4 1 186 12663
5 1 102 5317
1 0 2 18790
2 0 12 10673
3 0 28 5710
4 0 28 2585
5 0 31 1462
end
"""


def test_ir_istandard_pool_manual_homogeneity(run):
    # [ST] epitab, ir example 4 (Doll and Hill): istandard pool, level(90):
    # Test of homogeneity (direct) chi2(4) = 10.20, Pr>chi2 = 0.0372
    run(DOLL)
    out = run("ir deaths smokes pyears, by(agecat) level(90) istandard pool")
    assert "Pooled (direct)" in out and "I. Standardized" in out and "M-H combined" not in out
    assert " Test of homogeneity (direct) chi2(4) =     10.20  Pr>chi2 = 0.0372" in out
    a1 = [32, 104, 206, 186, 102]
    t1 = [52407, 43248, 28612, 12663, 5317]
    a0 = [2, 12, 28, 28, 31]
    t0 = [18790, 10673, 5710, 2585, 1462]
    # internal standardization: observed exposed deaths over those expected at the unexposed rates
    smr = sum(a1) / sum(t * b / u for t, b, u in zip(t1, a0, t0))
    assert run.session.r["irr"] == pytest.approx(smr)


def test_ir_istandard_ird(run):
    run(DOLL)
    out = run("ir deaths smokes pyears, by(agecat) istandard ird")
    assert "IRD" in out and "Test of homogeneity" not in out
    r = run.session.r
    d = [32 / 52407 - 2 / 18790, 104 / 43248 - 12 / 10673, 206 / 28612 - 28 / 5710,
         186 / 12663 - 28 / 2585, 102 / 5317 - 31 / 1462]
    w = [52407, 43248, 28612, 12663, 5317]
    assert r["ird"] == pytest.approx(sum(x * y for x, y in zip(w, d)) / sum(w))


def test_cs_standard_weights_and_rd(run):
    run(STRATA + "gen wgt = 1")
    run("cs case exp [fw=pop], by(sex) standard(wgt)")
    # equal weights: the ratio of the mean risks
    assert run.session.r["rr"] == pytest.approx((.1 + .25) / (.2 + .15))
    run("cs case exp [fw=pop], by(sex) standard(wgt) rd")
    assert run.session.r["rd"] == pytest.approx(((.1 - .2) + (.25 - .15)) / 2)
    run("cs case exp [fw=pop], by(sex) estandard")
    # external weights: the unexposed of each stratum (100 and 100)
    assert run.session.r["rr"] == pytest.approx((.1 + .25) / (.2 + .15))


def test_cs_pool_is_inverse_variance(run):
    run(STRATA)
    out = run("cs case exp [fw=pop], by(sex) pool")
    # pool adds the directly pooled row and its test; the M-H row and test stay (compat 1601)
    assert "Pooled (direct)" in out and "Test of homogeneity (direct)" in out and "M-H combined" in out
    lr = [math.log(.5), math.log(25 / 15)]
    v = [1 / 10 - 1 / 100 + 1 / 20 - 1 / 100, 1 / 25 - 1 / 100 + 1 / 15 - 1 / 100]
    m = sum(x / y for x, y in zip(lr, v)) / sum(1 / y for y in v)
    assert f"{math.exp(m):.6f}"[1:8] in out
    # r(chi2_p) in single precision, as Stata stores it
    assert run.session.r["chi2_p"] == pytest.approx(sum((x - m) ** 2 / y for x, y in zip(lr, v)), rel=1e-6)


def test_cs_cc_binomial_records(run):
    run(STRATA)
    run("cs case exp [fw=pop], by(sex)")
    rr = run.session.r["rr"]
    run("cc case exp [fw=pop], by(sex)")
    orr = run.session.r["or"]
    run("""clear
input cases n exp sex
10 100 1 0
20 100 0 0
25 100 1 1
15 100 0 1
end""")
    run("cs cases exp, by(sex) binomial(n)")
    assert run.session.r["rr"] == pytest.approx(rr)
    run("cc cases exp, binomial(n) by(sex)")
    assert run.session.r["or"] == pytest.approx(orr)


def test_cc_tarone_and_breslow_day(run):
    run(STRATA)
    out = run("cc case exp [fw=pop], by(sex) bd tarone")
    assert "Test of homogeneity (B-D)" in out and "Test of homogeneity (Tarone)" in out
    tables = [(10, 20, 90, 80), (25, 15, 75, 85)]
    psi = E.mh_or(tables, 95)[0]
    r = run.session.r
    assert r["chi2_bd"] == pytest.approx(E.breslow_day(tables, psi)[0])
    # Tarone's correction is never larger than the Breslow-Day statistic
    assert r["chi2_t"] <= r["chi2_bd"]


def test_cc_istandard_odds(run):
    run(STRATA)
    out = run("cc case exp [fw=pop], by(sex) istandard")
    assert "I. Standardized" in out and "Test that combined OR = 1" not in out
    # weights: exposed controls; sum a / sum c (b/d)
    expected = (10 + 25) / (90 * 20 / 80 + 75 * 15 / 85)
    assert run.session.r["or"] == pytest.approx(expected)


def test_epitab_option_checks(run):
    run(STRATA)
    # without by(), pool is ignored, as in Stata (compat 1601)
    run("cs case exp [fw=pop], pool")
    assert run.rc == 0
    run("cs case exp [fw=pop], by(sex) rd")
    assert run.rc == 198
    run("csi 10 20 90 80, istandard")
    assert run.rc == 198


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


def test_mhodds_by_strata_and_homogeneity(run):
    run(ALC)
    out = run("mhodds case alc [fw=pop], compare(3,1) by(age)")
    assert "by age" in out and "Test of homogeneity of ORs" in out
    r = run.session.r
    # overall: Mantel-Haenszel controlling for the by() variable
    num = 8 * 60 / 90 + 25 * 50 / 100
    den = 2 * 20 / 90 + 10 * 15 / 100
    assert r["or"] == pytest.approx(num / den)
    # homogeneity of the stratum log odds ratios (Woolf variances within a single table)
    b = math.log(num / den)
    lor = [math.log(8 * 60 / (2 * 20)), math.log(25 * 50 / (10 * 15))]
    v = [1 / 8 + 1 / 60 + 1 / 2 + 1 / 20, 1 / 25 + 1 / 50 + 1 / 10 + 1 / 15]
    assert r["chi2_hom"] == pytest.approx(sum((x - b) ** 2 / y for x, y in zip(lor, v)), rel=1e-6)
    assert r["df_hom"] == 1


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


def test_roctab_bamber_and_binomial_manual_values(run):
    run(HANLEY)
    out = run("roctab disease rating [fw=pop], bamber")
    # [R] roctab, Hanley and McNeil data: Bamber 109 0.8932 0.0306 0.83317 0.95317
    assert "109     0.8932       0.0306        0.83317     0.95317" in out
    run("roctab disease rating [fw=pop], binomial")
    # exact binomial limits: 0.81559 0.94180
    assert run.session.r["lb"] == pytest.approx(0.81559, abs=5e-6)
    assert run.session.r["ub"] == pytest.approx(0.94180, abs=5e-6)


def test_roctab_lorenz_gini_pietra(run):
    run(HANLEY)
    run("roctab disease rating [fw=pop]")
    area = run.session.r["area"]
    run("roctab disease rating [fw=pop], lorenz")
    r = run.session.r
    # Gini: twice the area between the Lorenz curve and the diagonal, 2 AUC - 1
    assert r["gini"] == pytest.approx(2 * area - 1, rel=1e-6)
    # Pietra: the largest sensitivity + specificity - 1 over the cutpoints
    sens = np.array([51, 48, 46, 44, 33]) / 51
    spec = np.array([0, 33, 39, 45, 56]) / 58
    assert r["pietra"] == pytest.approx(float((sens + spec - 1).max()))


def test_roccomp_binormal_matches_rocfit(run):
    run(HANLEY + "expand pop\nset seed 3\ngen r2 = rating + (disease==1)*(runiform()<.3)")
    run("rocfit disease rating, nolog")
    area, se = run.session.e["area"], run.session.e["se_area"]
    out = run("roccomp disease rating r2, binormal")
    # the binormal area and standard error of rocfit
    assert "Fitting binormal model for: rating" in out
    assert f"{area:.4f}       {se:.4f}" in out


def test_rocreg_nonparametric_auc_equals_roctab(run):
    run(HANLEY + "expand pop")
    run("roctab disease rating")
    auc = run.session.r["area"]
    run("rocreg disease rating, tiecorrected breps(50) bseed(1) nodots")
    e = run.session.e
    assert e["b"].data[0, 0] == pytest.approx(auc)
    assert e["V"].data[0, 0] > 0 and e["breps"] == "50"
    # without the tie correction, tied controls do not count below the case
    run("rocreg disease rating, nobootstrap")
    neg = np.repeat([1, 2, 3, 4, 5], [33, 6, 6, 11, 2])
    pos = np.repeat([1, 2, 3, 4, 5], [3, 2, 2, 11, 33])
    assert run.session.e["b"].data[0, 0] == pytest.approx(float((neg[None, :] < pos[:, None]).mean()))


def test_rocreg_summaries_and_strata(run):
    from opendta.estimation import rocreg as RR
    pv = np.array([.1, .5, .9, .95])
    assert RR.roc_at(pv, .5) == pytest.approx(.75)
    assert RR.pauc_at(pv, .2) == pytest.approx(np.mean(np.clip(.2 - (1 - pv), 0, None)))
    assert RR.invroc_at(pv, .5) == pytest.approx(.1)
    run(HANLEY + "expand pop\ngen g = mod(_n, 2)")
    run("rocreg disease rating, ctrlcov(g) tiecorrected auc roc(.3) nobootstrap")
    b = run.session.e["b"].data[0]
    assert run.session.e["ctrlmodel"] == "strata"
    # stratified controls: each case is placed among the controls of its own g
    d = np.array([run.eval(f"disease[{i}]") for i in range(1, 110)]) == 1
    y = np.array([run.eval(f"rating[{i}]") for i in range(1, 110)])
    g = np.arange(1, 110) % 2
    pv = np.array([RR.pv_empirical(y[i:i + 1], y[~d & (g == g[i])], True)[0] for i in np.flatnonzero(d)])
    # Stata's order: ROC(t) before the AUC (compat 1601)
    assert b[0] == pytest.approx(RR.roc_at(pv, .3))
    assert b[1] == pytest.approx(pv.mean())


def test_rocreg_probit_glm_binormal_data(run):
    run("""clear
set seed 1601
set obs 4000
gen d = _n > 2000
gen y = rnormal(cond(d, 1.5, 0), cond(d, 1.25, 1))
""")
    run("rocreg d y, probit nobootstrap fprpts(20)")
    b = run.session.e["b"].data[0]
    # binormal ROC: intercept mu/sigma1 = 1.2, slope sigma0/sigma1 = 0.8
    assert b[0] == pytest.approx(1.2, abs=0.1)
    assert b[1] == pytest.approx(0.8, abs=0.1)
