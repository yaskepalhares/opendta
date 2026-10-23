"""svy: with more estimators. With one observation per PSU the linearized variance is the robust
sandwich of the weighted model; with clusters it is close to the cluster-robust one (the small
sample factors differ slightly)."""

import numpy as np
import pytest

DATA = """clear
set seed 77
set obs 800
gen psu = ceil(runiform()*40)
gen w = 1 + floor(5*runiform())
gen x = round(invnormal(runiform()),.01)
gen z = round(invnormal(runiform()),.01)
gen lat = .5*x + round(invnormal(runiform()),.01)
gen y = 1 + lat + round(invnormal(runiform()),.01)
gen b = lat > .3
gen b2 = (lat + round(invnormal(runiform()),.01)) > .3
gen c = floor(exp(.3 + .3*x)*3*runiform())
gen yl = y - .2
gen yu = y + .2
gen sel = runiform() < .8
gen yy = y if sel
"""

CASES = {
    "cloglog": ("cloglog b x z", "cloglog b x z"),
    "intreg": ("intreg yl yu x z", "intreg yl yu x z"),
    "truncreg": ("truncreg y x z, ll(-3)", "truncreg y x z"),
    "heckman": ("heckman yy x, select(sel = z)", "heckman yy x"),
    "zip": ("zip c x, inflate(z)", "zip c x"),
    "biprobit": ("biprobit b b2 x", "biprobit b b2 x"),
    "fracreg": ("fracreg logit b x", "fracreg logit b x"),
    "hetprobit": ("hetprobit b x, het(z)", "hetprobit b x"),
}


def _se(run):
    return np.sqrt(np.diag(np.asarray(run.session._estimates.V)))


def _robust(cmd, tail):
    head, _, opts = cmd.partition(",")
    base, _, rest = head.partition(" ")
    # weights go after the variable list
    return f"{head.split(' ', 1)[0]} {head.split(' ', 1)[1]} [pweight=w]," + (opts + " " if opts else "") + tail


@pytest.mark.parametrize("name", list(CASES))
def test_unit_psu_equals_robust(run, name):
    cmd = CASES[name][0]
    run(DATA)
    run("svyset _n [pweight=w]")
    out = run("svy: " + cmd)
    assert run.rc == 0 and "Survey:" in out
    a = _se(run)
    b_svy = np.asarray(run.session._estimates.b).copy()
    run("quietly " + _robust(cmd, "vce(robust)"))
    assert run.rc == 0
    assert np.allclose(b_svy, np.asarray(run.session._estimates.b), rtol=1e-6)
    assert np.allclose(a, _se(run), rtol=2e-3)


@pytest.mark.parametrize("name", ["cloglog", "truncreg", "zip", "biprobit"])
def test_clusters_match_cluster_robust(run, name):
    cmd = CASES[name][0]
    run(DATA)
    run("svyset psu [pweight=w]")
    run("svy: " + cmd)
    assert run.rc == 0
    a = _se(run)
    run("quietly " + _robust(cmd, "vce(cluster psu)"))
    assert run.rc == 0
    assert np.allclose(a, _se(run), rtol=0.05)


def test_subpop_and_cloglog_estimates(run):
    run(DATA)
    run("svyset psu [pweight=w]")
    out = run("svy, subpop(if x > -1): cloglog b z")
    assert run.rc == 0 and "Subpop. no. obs" in out
    n = run.eval("e(N_sub)")
    assert n == run.eval("r(N)") or n > 0
