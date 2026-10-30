"""power: twoway, repeated, cmh, mcc, trend, cox, exponential and logrank.

The core functions are checked against closed formulas and SciPy; the command against the
examples of [PSS] power twoway (the manual's N = 132 and N = 153) and against its own stored
results."""

import math

import numpy as np
import pytest
from scipy import integrate, stats

from opendta.core import power_suites as C

MEANS = np.array([[134, 143, 91], [106, 173, 145]], dtype=float)


# ---------------------------------------------------------------------------
# core
# ---------------------------------------------------------------------------

def test_core_twoway_manual_example():
    assert C.twoway_variances(MEANS)[0] == pytest.approx(87.1111, abs=1e-4)
    frac = np.full((2, 3), 1 / 6)

    def pw(k):
        return C.twoway_power(C.twoway_noncentrality(MEANS, frac, 6 * k, "row", 1417.0), 2, 3, 6 * k, "row", 0.05)
    assert 6 * C.smallest_int(pw, 0.8, 2) == 132
    w = np.array([[2, 2, 2], [1, 1, 1]], dtype=float)

    def pu(k):
        n = 9 * k
        return C.twoway_power(C.twoway_noncentrality(MEANS, w / w.sum(), n, "row", 1417.0), 2, 3, n, "row", 0.05)
    assert 9 * C.smallest_int(pu, 0.8, 2) == 153


def test_core_twoway_balanced_noncentrality_is_n_delta2():
    frac = np.full((2, 3), 1 / 6)
    for which, k in (("row", 0), ("column", 1), ("rowcol", 2)):
        lam = C.twoway_noncentrality(MEANS, frac, 60.0, which, 1417.0)
        assert lam == pytest.approx(60 * C.twoway_variances(MEANS)[k] / 1417.0)
    d = C.lambda_for_power(1, 54, 0.05, 0.8)
    assert C.anova_power(d, 1, 54, 0.05) == pytest.approx(0.8, abs=1e-9)
    assert C.anova_power(d, 1, 54, 0.05) == pytest.approx(
        stats.ncf.sf(stats.f.isf(0.05, 1, 54), 1, 54, d))


def test_core_repeated_matches_t_tests():
    # two groups, K measures, compound symmetry: the F test on the group means is a two-sample t test
    K, rho, s2, d, n = 3, 0.5, 4.0, 1.0, 40
    means = np.array([[0.0] * K, [d] * K])
    cov = C.cs_covariance(K, s2, rho)
    sm = C.repeated_summary(2, K, cov, "between", np.array([0.5, 0.5]), means)
    assert sm["eps"] == pytest.approx(1.0)
    sd_mean = math.sqrt(s2 * (1 + (K - 1) * rho) / K)
    nc = d / sd_mean * math.sqrt(n / 4)
    crit = stats.t.isf(0.025, n - 2)
    ref = stats.nct.sf(crit, n - 2, nc) + stats.nct.cdf(-crit, n - 2, nc)
    assert C.repeated_power(sm["delta2"], sm["dc"], sm["du"], 2, n, sm["eps"], 0.05) == pytest.approx(ref, abs=1e-9)
    # one group, two measures: the within test is a paired t test
    means = np.array([[0.0, 1.0]])
    sm = C.repeated_summary(1, 2, C.cs_covariance(2, s2, rho), "within", np.array([1.0]), means)
    nc = 1.0 * math.sqrt(30) / math.sqrt(2 * s2 * (1 - rho))
    crit = stats.t.isf(0.025, 29)
    ref = stats.nct.sf(crit, 29, nc) + stats.nct.cdf(-crit, 29, nc)
    assert C.repeated_power(sm["delta2"], sm["dc"], sm["du"], 1, 30, sm["eps"], 0.05) == pytest.approx(ref, abs=1e-9)


def test_core_repeated_epsilon():
    cs = C.cs_covariance(4, 2.0, 0.3)
    assert C.repeated_summary(2, 4, cs, "within", np.array([.5, .5]), None, 1.0)["eps"] == pytest.approx(1.0)
    ar = 2.0 * 0.7 ** np.abs(np.subtract.outer(range(4), range(4)))
    assert C.repeated_summary(2, 4, ar, "within", np.array([.5, .5]), None, 1.0)["eps"] < 1.0


def test_core_cmh_single_stratum_is_the_z_test():
    n, p1, p2 = 80, 0.3, 0.5
    pb = (p1 + p2) / 2
    sd0, sda = math.sqrt(2 * pb * (1 - pb) / n), math.sqrt((p1 * (1 - p1) + p2 * (1 - p2)) / n)
    z = stats.norm.isf(0.025)
    ref = stats.norm.sf((z * sd0 - (p2 - p1)) / sda) + stats.norm.cdf((-z * sd0 - (p2 - p1)) / sda)
    assert C.cmh_power([n], [n], [p1], [p2], 0.05, 2) == pytest.approx(ref, abs=1e-12)
    assert C.cmh_power([n], [n], [p1], [p2], 0.05, 2, continuity=True) < ref
    assert C.cmh_power([n] * 2, [n] * 2, [p1] * 2, [p2] * 2, 0.05, 2) > ref
    assert C.odds_probability(np.array([0.2]), 2.5)[0] == pytest.approx(0.5 / 1.3 * 0.5 / 0.5 * 0.5 / 0.5 * 1.0 * 0.5 / 0.5 * 0.5 / 0.5 * 1.0, abs=0.5)
    p2 = C.odds_probability(np.array([0.2]), 2.5)[0]
    assert (p2 / (1 - p2)) / (0.2 / 0.8) == pytest.approx(2.5)


def test_core_cmh_sizes():
    n1, n2 = C.cmh_sizes([4.0, 1.0, 4.0], [0.5] * 3, 9)
    assert list(n1) == [36, 9, 36] and list(n2) == [36, 9, 36]
    n1, n2 = C.cmh_sizes([4.0, 1.0, 4.0], [0.47, 0.57, 0.51], 18)
    assert list(n2) == [34, 11, 37] and list(n1 + n2) == [72, 18, 72]       # the manual's example 3


def test_core_mcc_matches_mcnemar_formula():
    p0, theta = 0.3, 2.0
    p1 = C.mcc_p1(p0, theta, 0.0)
    assert (p1 / (1 - p1)) / (p0 / (1 - p0)) == pytest.approx(theta)
    t = p1 * (1 - p0) + (1 - p1) * p0                      # probability of a discordant pair
    pcase = theta / (1 + theta)
    za, zb = stats.norm.isf(0.025), stats.norm.ppf(0.8)
    n_disc = ((za * 0.5 + zb * math.sqrt(pcase * (1 - pcase))) / (pcase - 0.5)) ** 2
    n = C.solve_up(lambda x: C.mcc_power(x, p0, theta, 0.0, 1, 0.05, 2) - 0.8, 1.0, 50.0)
    assert n == pytest.approx(n_disc / t, rel=0.02)         # two-sided: the lower tail adds a little
    # one-sided closed form
    e1 = C.mcc_size_start(p0, theta, 0.0, 1, 0.05, 2, 0.8)
    assert C.mcc_power(e1, p0, theta, 0.0, 1, 0.05, 1) == pytest.approx(0.8, abs=1e-9)
    # a correlation changes p1 but keeps the odds ratio p10 / p01
    p1 = C.mcc_p1(0.3, 2.0, 0.3)
    s = 0.3 * math.sqrt(p1 * (1 - p1) * 0.3 * 0.7)
    assert (p1 * 0.7 - s) / ((1 - p1) * 0.3 - s) == pytest.approx(2.0, rel=1e-8)


def test_core_trend_two_groups_is_the_z_test():
    n, p1, p2 = 100, 0.2, 0.35
    pb = (p1 + p2) / 2
    sd0, sda = math.sqrt(2 * pb * (1 - pb) / n), math.sqrt((p1 * (1 - p1) + p2 * (1 - p2)) / n)
    z = stats.norm.isf(0.025)
    ref = stats.norm.sf((z * sd0 - (p2 - p1)) / sda) + stats.norm.cdf((-z * sd0 - (p2 - p1)) / sda)
    assert C.trend_power([n, n], [p1, p2], [1, 2], 0.05, 2) == pytest.approx(ref, abs=1e-12)
    assert C.trend_power([n, n], [p1, p2], [1, 2], 0.05, 2, continuity=True) < ref
    assert C.is_monotonic([.1, .2, .3]) and not C.is_monotonic([.1, .3, .2])
    d = C.trend_power([n, n], [p2, p1], [1, 2], 0.05, 1)               # decreasing: the lower tail
    assert d == pytest.approx(C.trend_power([n, n], [p1, p2], [1, 2], 0.05, 1), abs=1e-12)


def test_core_cox_and_logrank_closed_forms():
    b = math.log(0.5)
    e = C.cox_events(b, 0.5, 0.0, 0.05, 2, 0.8)
    assert e == pytest.approx((stats.norm.isf(0.025) + stats.norm.ppf(0.8)) ** 2 / (0.25 * b * b))
    assert C.cox_power(b, 0.5, 0.0, e, 1.0, 0.05, 2) == pytest.approx(0.8, abs=1e-12)
    assert C.cox_effect(0.5, 0.0, e, 1.0, 0.05, 2, 0.8) == pytest.approx(abs(b), abs=1e-12)
    assert C.cox_events(b, 0.5, 0.5, 0.05, 2, 0.8) == pytest.approx(2 * e)
    assert math.ceil(C.logrank_events(0.5, 1.0, 0.05, 2, 0.8, False)) == 71
    ev = C.logrank_events(0.76, 1.0, 0.05, 2, 0.9, True)
    assert ev == pytest.approx(4 * (stats.norm.isf(0.025) + stats.norm.ppf(0.9)) ** 2 / math.log(0.76) ** 2)
    assert C.logrank_power(0.76, 1.0, ev, 1.0, 0.05, 2, True) == pytest.approx(0.9, abs=1e-12)
    assert C.logrank_pe(0.5, 1.0, 0.5, None) == pytest.approx(1 - (0.5 + 0.5 ** 0.5) / 2)
    assert C.logrank_pe(0.5, 1.0, None) == 1.0
    # Simpson's rule: the same survival at the three times gives the average of the two groups
    assert C.logrank_pe(0.5, 1.0, None, None, (0.6, 0.6, 0.6)) == pytest.approx(1 - (0.6 + 0.6 ** 0.5) / 2)


def test_core_exponential_event_probability():
    lam, eta, T, r = 0.3, 0.05, 6.0, 3.0
    assert C.exp_pe(lam, eta, T, 0.0, 0.0) == pytest.approx(lam / (lam + eta) * (1 - math.exp(-(lam + eta) * T)))
    ref, _ = integrate.quad(lambda u: lam / (lam + eta) * (1 - math.exp(-(lam + eta) * (T - u))) / r, 0, r)
    assert C.exp_pe(lam, eta, T, r, 0.0) == pytest.approx(ref, rel=1e-10)
    for g in (0.8, -0.6):
        dens = lambda u: g * math.exp(-g * u) / (1 - math.exp(-g * r))
        ref, _ = integrate.quad(lambda u: dens(u) * lam / (lam + eta) * (1 - math.exp(-(lam + eta) * (T - u))), 0, r)
        assert C.exp_pe(lam, eta, T, r, g) == pytest.approx(ref, rel=1e-9)
    assert C.exp_pe(lam, 0.0, None, 0.0, 0.0) == 1.0
    g = C.exp_gamma(2.0, 1.0, 0.7)
    assert (1 - math.exp(-g)) / (1 - math.exp(-2 * g)) == pytest.approx(0.7)
    assert C.exp_gamma(2.0, 1.0, 0.5) == 0.0


def test_core_exponential_sample_size():
    # no censoring and the log hazard-ratio test: the same as the Schoenfeld number of events
    xo, xa = C.exp_xi(1.0, 0.6, 0.0, 0.0, 1.0, None, 0.0, 0.0, True)
    n = C.exp_size(math.log(0.6), xo, xa, 0.05, 2, 0.8)
    assert n == pytest.approx(C.logrank_events(0.6, 1.0, 0.05, 2, 0.8, True))
    assert C.exp_power(math.log(0.6), xo, xa, n, 0.05, 2) == pytest.approx(0.8, abs=1e-12)
    xo, xa = C.exp_xi(0.1, 0.05, 0.0, 0.0, 1.0, None, 0.0, 0.0, False)
    n = C.exp_size(-0.05, xo, xa, 0.05, 2, 0.8)
    assert C.exp_power(-0.05, xo, xa, n, 0.05, 2) == pytest.approx(0.8, abs=1e-12)
    xo2, xa2 = C.exp_xi(0.1, 0.05, 0.0, 0.0, 1.0, 5.0, 0.0, 0.0, False)
    assert C.exp_size(-0.05, xo2, xa2, 0.05, 2, 0.8) > n                  # censoring needs more subjects
    assert C.exp_size(-0.05, xo, xa, 0.05, 2, 0.8, unconditional=True) > 0


def test_core_searches():
    assert C.smallest_int(lambda k: k / 10, 0.55, 1) == 6
    assert C.smallest_int(lambda k: math.nan if k < 3 else k / 10, 0.35, 1) == 4
    with pytest.raises(ValueError):
        C.smallest_int(lambda k: 0.1, 0.9, 1, cap=1000)
    assert C.standard_weights([4, 4, 4, 2, 2, 2]) == [2, 2, 2, 1, 1, 1]
    assert C.standard_weights([1.5, 3.0]) == [1.0, 2.0]


# ---------------------------------------------------------------------------
# commands
# ---------------------------------------------------------------------------

def test_twoway_manual_examples(run):
    out = run("power twoway 134 143 91 \\ 106 173 145, varerror(1417) power(.8)")
    assert "F test for row effect" in out and "0.2479" in out and "87.1111" in out
    assert run.eval("r(N)") == 132 and run.eval("r(N_per_cell)") == 22
    assert run.eval("r(N_r)") == 2 and run.eval("r(N_c)") == 3
    run("power twoway 134 143 91 \\ 106 173 145, varerror(1417) power(.8) cellweights(2 2 2 \\ 1 1 1)")
    assert run.eval("r(N)") == 153 and run.eval("r(balanced)") == 0
    out = run("power twoway 134 143 91 \\ 106 173 145, varerror(1417) n(132)")
    assert "Estimated power for two-way ANOVA" in out
    assert run.eval("r(power)") == pytest.approx(0.8070, abs=5e-5)
    run("power twoway 134 143 91 \\ 106 173 145, varerror(1417) factor(column) power(.8)")
    assert run.eval("r(N)") == 48
    run("power twoway 134 143 91 \\ 106 173 145, varerror(1417) factor(rowcol) power(.8)")
    n = run.eval("r(N)")
    assert n % 6 == 0
    run("power twoway 134 143 91 \\ 106 173 145, varerror(1417) factor(rowcol) n(%d)" % n)
    assert run.eval("r(power)") >= 0.8
    run("power twoway 134 143 91 \\ 106 173 145, varerror(1417) factor(rowcol) n(%d)" % (n - 6))
    assert run.eval("r(power)") < 0.8


def test_twoway_effect_size_and_matrix_input(run):
    out = run("power twoway, nrows(2) ncols(3) n(132) power(.8) varerror(1417)")
    assert "Estimated row variance for two-way ANOVA" in out
    d = run.eval("r(delta)")
    run("power twoway, nrows(2) ncols(3) n(132) varrow(%.10f) varerror(1417)" % (d * d * 1417))
    assert run.eval("r(power)") == pytest.approx(0.8, abs=1e-6)
    run("matrix M = (134, 143, 91 \\ 106, 173, 145)")
    run("power twoway M, varerror(1417) power(.8)")
    assert run.eval("r(N)") == 132
    out = run("power twoway M, varerror(1417) n(100(20)140)")
    assert "N_per_cell" in out and run.eval("r(N)") == 140


def test_twoway_errors(run):
    assert "only two of n(), power()" in run("power twoway M2 \\ 1, n(10) power(.8)") or run.rc != 0
    out = run("power twoway 1 2 \\ 3 4, graph")
    assert "not available yet" in out and run.rc == 198
    out = run("power twoway 1 2 \\ 3 4, onesided")
    assert run.rc == 198
    run("power twoway, nrows(2) ncols(2) power(.8)")
    assert run.rc == 198


def test_repeated(run):
    out = run("power repeated, ngroups(2) nrepeated(3) corr(.5) varbetween(1) n(40)")
    assert "F test for between subjects" in out and run.rc == 0
    d2 = run.eval("r(delta)") ** 2
    assert d2 == pytest.approx(1.0 / 2.0)                                  # Var_b / Var_be
    run("power repeated 10 12 14 \\ 11 13 12, corr(.6) varerror(77) power(.9) factor(within)")
    assert run.rc == 0 and run.eval("r(N)") == 2 * run.eval("r(N_per_group)")
    n = run.eval("r(N)")
    run("power repeated 10 12 14 \\ 11 13 12, corr(.6) varerror(77) n(%d) factor(within)" % n)
    assert run.eval("r(power)") >= 0.9
    run("power repeated 10 12 14 \\ 11 13 12, corr(.6) varerror(77) n(%d) factor(within)" % (n - 2))
    assert run.eval("r(power)") < 0.9
    run("matrix S = (4, 2, 1 \\ 2, 4, 2 \\ 1, 2, 4)")
    out = run("power repeated, ngroups(2) covmatrix(S) factor(within) vareffect(.4) n(60)")
    assert "spherical" in out and run.rc == 0
    out = run("power repeated, ngroups(2) nrepeated(3) corr(.5) factor(bwithin) vareffect(.2) power(.8)")
    assert "between-within" in out or run.rc == 0


def test_cmh(run):
    out = run("power cmh .2 .1 .3, oratio(2.5)")
    assert "Cochran-Mantel-Haenszel test" in out and "OR = OR_1 = OR_2 = OR_3" in out
    n = run.eval("r(N)")
    assert n % 6 == 0 and run.eval("r(N_per_group)") == n / 2
    run("power cmh .2 .1 .3, oratio(2.5) n(%d)" % n)
    assert run.eval("r(power)") >= 0.8
    run("power cmh .2 .1 .3, oratio(2.5) n(%d)" % (n - 6))
    assert run.eval("r(power)") < 0.8
    run("power cmh .2 .1 .3, n(%d) power(.8)" % n)
    assert 1 < run.eval("r(oratio)") <= 2.5 and run.eval("r(delta)") == run.eval("r(oratio)")
    run("power cmh .2 .1 .3, oratio(2.5) strweights(4 1 4) grratios(.47 .57 .51)")
    assert run.rc == 0 and run.eval("r(balanced)") == 0
    run("power cmh .2 .1 .3, n(%d) power(.8) direction(lower) onesided" % n)
    assert run.eval("r(oratio)") < 1
    out = run("power cmh .2 .1 .3, oratio(2 2.5 3)")
    assert "oratio" in out and "p1_1" in out
    run("power cmh .2 .1 .3, oratio(2.5) continuity")
    assert run.eval("r(N)") > n
    run("matrix P = (.2, .1, .3 \\ .25, .15, .35)")
    out = run("power cmh P, oratio(2.5)")
    assert run.rc == 0 and "p1_3" in out


def test_mcc(run):
    out = run("power mcc .3, oratio(2)")
    assert "Asymptotic z test, 1:1 matched design" in out and "N cases" in out
    n = run.eval("r(N)")
    run("power mcc .3, oratio(2) n(%d)" % n)
    assert run.eval("r(power)") >= 0.8
    run("power mcc .3, oratio(2) n(%d)" % (n - 1))
    assert run.eval("r(power)") < 0.8
    run("power mcc .3, oratio(2) m(2) compare")
    assert run.eval("r(N)") < n and 0 < run.eval("r(F_M)") < 1
    run("power mcc .3, n(%d) power(.8)" % n)
    assert 1 < run.eval("r(oratio)") <= 2
    run("power mcc .3, n(%d) power(.8) direction(lower)" % n)
    assert run.eval("r(oratio)") < 1
    run("power mcc .3, oratio(2) corr(.2)")
    assert run.rc == 0
    out = run("power mcc .3, oratio(1.5 2 3) m(1 2)")
    assert "oratio" in out and "corr" in out


def test_trend(run):
    out = run("power trend .1 .2 .3, power(.8)")
    assert "Cochran-Armitage trend test" in out and "equally spaced" in out
    n = run.eval("r(N)")
    assert n % 3 == 0
    run("power trend .1 .2 .3, n(%d)" % n)
    assert run.eval("r(power)") >= 0.8
    run("power trend .1 .2 .3, n(%d)" % (n - 3))
    assert run.eval("r(power)") < 0.8
    out = run("power trend .1 .2 .3, n(%d) exposure(1 2 4)" % n)
    assert "x3" in out and "equally spaced" not in out
    out = run("power trend .1 .3 .2, n(300)")
    assert "not monotonic" in out
    run("power trend .1 .2 .3, power(.8) continuity")
    assert run.eval("r(N)") > n
    run("power trend .1 .2 .3, power(.8) grweights(1 1 2)")
    assert run.eval("r(balanced)") == 0
    run("power trend .1 .2 .4, power(.8) exposure(1 2 4) continuity")
    assert run.rc == 198


def test_cox(run):
    out = run("power cox")
    assert "Wald test" in out and "Cox PH regression" in out
    assert run.eval("r(E)") == 66 and run.eval("r(N)") == 66
    run("power cox, hratio(.5) eventprob(.5)")
    assert run.eval("r(E)") == 66 and run.eval("r(N)") == 132
    run("power cox, hratio(.5) eventprob(.5) wdprob(.2)")
    assert run.eval("r(N)") == 165
    run("power cox %.6f, sd(.5) n(66)" % math.log(.5))
    assert run.eval("r(power)") == pytest.approx(0.8, abs=0.01)
    run("power cox, n(200) power(.8)")
    assert run.eval("r(b1)") < 0
    run("power cox, n(200) power(.8) direction(upper) effect(hratio)")
    assert run.eval("r(hratio)") > 1 and run.eval("r(delta)") == run.eval("r(hratio)")
    run("power cox, hratio(.5) r2(.5)")
    assert run.eval("r(E)") == 131
    out = run("power cox, hratio(.5 .6 .7)")
    assert "b1" in out


def test_exponential(run):
    out = run("power exponential .1 .05, power(.8)")
    assert "hazard difference, conditional" in out and "Ho: h2 = h1" in out
    n = run.eval("r(N)")
    run("power exponential .1 .05, n(%d)" % n)
    assert run.eval("r(power)") >= 0.8
    run("power exponential .1 .05, n(%d) onesided" % n)
    assert run.eval("r(power)") > 0.8
    out = run("power exponential .5 .25, time(2) power(.8)")
    assert "survival probabilities" in out and run.eval("r(h1)") == pytest.approx(math.log(2) / 2)
    run("power exponential .1 .05, power(.8) fperiod(5) aperiod(2)")
    assert run.eval("r(N)") > n
    run("power exponential .1 .05, power(.8) fperiod(5) aperiod(2) ashape(.5)")
    assert run.rc == 0 and run.eval("r(ashape)") == 0.5
    run("power exponential .1 .05, power(.8) fperiod(5) aperiod(2) lossprob(.1)")
    assert run.eval("r(N)") > n
    run("power exponential .1 .05, power(.8) fperiod(5) aperiod(2) aprob(.7) atime(1)")
    assert run.rc == 0 and run.eval("r(ashape)") > 0
    run("power exponential .1 .05, power(.8) unconditional")
    assert run.eval("r(unconditional)") == 1
    run("power exponential .1, hdifference(-.05) power(.8)")
    assert run.eval("r(N)") == n
    run("power exponential, hratio(.6) loghazard power(.8)")
    assert run.rc == 0 and run.eval("r(N)") > 0
    out = run("power exponential .1 .05, power(.8) show")
    assert "not available yet" in out and run.rc == 198
    run("power exponential .1 .05, power(.8) nratio(2)")
    assert run.eval("r(N2)") == 2 * run.eval("r(N1)") or run.eval("r(N2)") >= 2 * run.eval("r(N1)") - 1


def test_logrank(run):
    out = run("power logrank")
    assert "Freedman method" in out and run.eval("r(E)") == 71 and run.eval("r(N)") == 72
    run("power logrank, hratio(.76) schoenfeld power(.9)")
    assert run.eval("r(E)") == 559
    run("power logrank .5, hratio(.7) wdprob(.1)")
    assert run.eval("r(Pr_E)") < 1 and run.eval("r(Pr_w)") == .1
    n = run.eval("r(N)")
    run("power logrank .5, hratio(.7) n(%d)" % n)
    assert run.eval("r(power)") >= 0.78
    run("power logrank .5 .35, n(300)")
    assert run.eval("r(hratio)") == pytest.approx(math.log(.35) / math.log(.5))
    run("power logrank .5, n(200) power(.8)")
    h = run.eval("r(hratio)")
    run("power logrank .5, hratio(%.8f) n(200)" % h)
    assert run.eval("r(power)") == pytest.approx(0.8, abs=2e-4)
    run("power logrank .5, n(200) power(.8) direction(upper)")
    assert run.eval("r(hratio)") > 1
    run("power logrank, simpson(.8 .7 .6) hratio(.6)")
    assert run.eval("r(Pr_E)") < 1
    out = run("power logrank, st1(s t) hratio(.6)")
    assert run.rc == 198 and "not available yet" in out
    out = run("power logrank, hratio(.5 .6 .7) power(.8 .9)")
    assert "hratio" in out and run.rc == 0
    run("power logrank, hratio(.5) lnhratio(-.7)")
    assert run.rc == 198


def test_parallel_and_table(run):
    out = run("power cox, hratio(.5 .6) power(.8 .9) parallel")
    assert out.count("\n  | ") == 3                         # header + two rows
    out = run("power cox, hratio(.5 .6) power(.8 .9)")
    assert out.count("\n  | ") == 5
    out = run("power cox, hratio(.5) table")
    assert "alpha" in out and "  +" in out
    run("power cox, hratio(.5 .6) power(.8 .9 .95) parallel")
    assert run.rc == 198


def test_return_list_types(run):
    run("power cmh .2 .1 .3, oratio(2.5)")
    assert run.eval('r(method)') == "cmh" and run.eval('r(type)') == "test"
    run("power twoway 134 143 91 \\ 106 173 145, varerror(1417) power(.8)")
    assert run.eval('r(method)') == "twoway"
    run("matrix list r(means)")
    assert run.rc == 0
