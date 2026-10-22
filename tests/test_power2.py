"""power: onevariance, twovariances, pairedproportions, oneway, the lrchi2 and fisher tests of
twoproportions and the effect options. The reference values are those of the examples of [PSS]
power (the manual) and independent computations."""

import math

import numpy as np
import pytest
from scipy import stats


def test_onevariance(run):
    run("power onevariance 4 9")
    assert run.eval("r(N)") == 24
    run("power onevariance 4 9, n(30)")
    assert run.eval("r(power)") == pytest.approx(0.8827, abs=5e-5)
    run("power onevariance 2 3, sd")
    assert run.eval("r(N)") == 24
    out = run("power onevariance 4, n(30) power(.8)")
    assert "8.1371" in out and "2.0343" in out
    out = run("power onevariance 4, n(30) power(.8) direction(lower)")
    assert "1.8267" in out and "0.4567" in out


def test_twovariances(run):
    run("power twovariances 4 2.25")
    assert run.eval("r(N)") == 194 and run.eval("r(N1)") == 97
    run("power twovariances 2 1.5, sd")
    assert run.eval("r(N)") == 194
    run("power twovariances 4, ratio(.5625)")
    assert run.eval("r(N)") == 194
    run("power twovariances 4 2.25, nratio(2)")
    assert (run.eval("r(N1)"), run.eval("r(N2)")) == (75, 150)
    run("power twovariances 4 2.25, n(250)")
    assert run.eval("r(power)") == pytest.approx(0.8908, abs=5e-5)
    out = run("power twovariances 4, n(250) power(.8)")
    assert "6.6291" in out and "1.6573" in out
    out = run("power twovariances 4 (1.5(.5)3), n(250)")
    assert ".9997" in out and ".3572" in out


def test_pairedproportions(run):
    run("power pairedproportions .105 .004")
    assert run.eval("r(N)") == 82
    run("power pairedproportions .105 .004, n(100)")
    assert run.eval("r(power)") == pytest.approx(0.8759, abs=5e-5)
    run("power pairedproportions .105, prdiscordant(.109) n(100)")
    assert run.eval("r(power)") == pytest.approx(0.8759, abs=5e-5)
    run("power pairedproportions .53 .4293, corr(.8)")
    assert run.eval("r(N)") == 82
    assert run.eval("r(delta)") == pytest.approx(-0.1007, abs=5e-5)
    run("power pairedproportions .53, oratio(.667) corr(.8)")
    assert run.eval("r(N)") == 82
    # McNemar's test against simulation of the trinomial counts
    rng = np.random.default_rng(1)
    n, p12, p21 = 60, .15, .05
    x = rng.multinomial(n, [p12, p21, 1 - p12 - p21], size=40000)
    stat = (x[:, 0] - x[:, 1]) ** 2 / np.maximum(x[:, 0] + x[:, 1], 1)
    run(f"power pairedproportions {p12} {p21}, n({n})")
    assert run.eval("r(power)") == pytest.approx(np.mean(stat > stats.chi2.isf(.05, 1)), abs=0.03)


def test_twoproportions_tests(run):
    run("power twoproportions .015 .001")
    assert run.eval("r(N)") == 1270
    run("power twoproportions .015 .001, test(lrchi2)")
    assert run.eval("r(N)") == 1062
    run("power twoproportions .015 .001, n(1100)")
    assert run.eval("r(power)") == pytest.approx(0.7416, abs=5e-5)
    run("power twoproportions .015 .001, test(lrchi2) n(1100)")
    assert 0.7 < run.eval("r(power)") < 0.9
    run("power twoproportions .1 .2, test(fisher)")
    assert run.rc == 198
    run("power twoproportions .1 .2, test(chi2) continuity test(lrchi2) n(100)")
    assert run.rc == 198


def test_fisher_power_brute_force(run):
    n1 = n2 = 12
    p1, p2, a = .2, .6, .05
    # the rejection region from the conditional distribution, built table by table
    pw = 0.0
    for x1 in range(n1 + 1):
        for x2 in range(n2 + 1):
            m = x1 + x2
            h = stats.hypergeom(n1 + n2, n2, m)
            lo = max(0, m - n1)
            hi = min(m, n2)
            p_up = h.sf(x2 - 1)
            p_lo = h.cdf(x2)
            if p_up <= a / 2 + 1e-12 or p_lo <= a / 2 + 1e-12:
                pw += stats.binom.pmf(x1, n1, p1) * stats.binom.pmf(x2, n2, p2)
    run("power twoproportions .2 .6, test(fisher) n(24)")
    assert run.eval("r(power)") == pytest.approx(pw, abs=1e-9)
    # Fisher's exact test is conservative: the size is below alpha
    assert run.eval("r(power)") < stats.norm.cdf(1.0)


def test_effect_options(run):
    run("power twoproportions .1, oratio(2)")
    p2 = 2 * .1 / (1 - .1 + 2 * .1)
    run(f"power twoproportions .1 {p2}")
    n_direct = run.eval("r(N)")
    run("power twoproportions .1, oratio(2)")
    assert run.eval("r(N)") == n_direct
    run("power twoproportions .1, ratio(2) n(100)")
    assert run.rc == 0
    out = run("power twoproportions .1, rrisk(2) effect(rrisk) n(100) test(lrchi2)")
    assert "(relative risk)" in out and "2.0000" in out


def test_oneway(run):
    f2 = (2 / 3) / 4
    run("power oneway 1 2 3, varerror(4)")
    n = run.eval("r(N_per_group)")
    g = 3

    def pw(m):
        d2 = g * m - g
        return stats.ncf.sf(stats.f.isf(.05, g - 1, d2), g - 1, d2, g * m * f2)
    assert pw(n) >= .8 > pw(n - 1)
    run("power oneway 1 2 3, varerror(4) n(60)")
    assert run.eval("r(power)") == pytest.approx(pw(20), abs=1e-9)
    rng = np.random.default_rng(3)
    y = rng.normal(size=(20000, 3, 20)) * 2 + np.array([1, 2, 3])[None, :, None]
    msb = 20 * y.mean(2).var(1, ddof=1)
    mse = y.var(2, ddof=1).mean(1)
    sim = np.mean(msb / mse > stats.f.isf(.05, 2, 57))
    assert pw(20) == pytest.approx(sim, abs=0.02)


def test_graph(run):
    run("power twovariances 4 (1.5 2 3), n(250) graph")
    assert run.rc == 0
