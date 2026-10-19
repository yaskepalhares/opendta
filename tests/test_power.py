"""power ([PSS]): sample size, power and effect size of the essential methods, against the
closed-form results of the manuals and textbooks."""

import math

import pytest
from scipy import stats


def test_onemean(run):
    out = run("power onemean 100 102, sd(5)")
    assert "Estimated sample size for a one-sample mean test" in out and "t test" in out
    n = run.eval("r(N)")
    # smallest n with t-test power >= .8
    def pw(k):
        nc = .4 * math.sqrt(k)
        c = stats.t.isf(.025, k - 1)
        return stats.nct.sf(c, k - 1, nc) + stats.nct.cdf(-c, k - 1, nc)
    assert pw(n) >= .8 > pw(n - 1)
    run("power onemean 100 102, sd(5) knownsd")
    z = (stats.norm.isf(.025) + stats.norm.isf(.2)) / .4
    assert run.eval("r(N)") == math.ceil(z ** 2)
    run("power onemean 100 102, sd(5) n(40)")
    assert run.eval("r(power)") == pytest.approx(pw(40), rel=1e-6)
    # Stata 14: with n() and no power(), the alternative mean is required (compat 1301)
    run("power onemean 100, sd(5) n(40)")
    assert run.rc == 198
    out = run("power onemean 100, sd(5) n(40) power(.8)")
    assert "Estimated effect size and target mean" in out
    assert pw(40) < .8 and run.eval("r(ma)") > 102


def test_twomeans_and_paired(run):
    run("power twomeans 12 14, sd(4)")
    assert run.eval("r(N_per_group)") == 64 and run.eval("r(N)") == 128
    run("power twomeans 12 14, sd1(4) sd2(6) nratio(2)")
    assert run.eval("r(N2)") == 2 * run.eval("r(N1)")
    run("power pairedmeans 10 12, sddiff(5)")
    assert run.eval("r(N)") == 52


def test_proportions_and_correlations(run):
    run("power twoproportions .3 .45")
    assert run.eval("r(N_per_group)") == 163
    run("power oneproportion .3 .4")
    z1, z2 = stats.norm.isf(.025), stats.norm.isf(.2)
    n = ((z1 * math.sqrt(.21) + z2 * math.sqrt(.24)) / .1) ** 2
    assert run.eval("r(N)") == math.ceil(n)
    run("power onecorrelation 0 .3")
    assert run.eval("r(N)") == 85
    run("power twocorrelations .3 .5, n(400)")
    assert 0 < run.eval("r(power)") < 1


def test_table_of_values(run):
    out = run("power onemean 100 (102 104 106), sd(5)")
    assert out.count("|     .05") == 3


def test_errors(run):
    run("capture power onemean 100 102, n(30) power(.9)")
    assert run.eval("_rc") == 198
    run("capture power nothing 1 2")
    assert run.eval("_rc") == 198
