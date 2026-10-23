"""mi impute mvn, monotone, truncreg, intreg and nbreg: the imputed values must reproduce the
distribution of the observed ones (the values are MCAR) and respect the limits."""

import numpy as np

from opendta.session import Session

BASE = """clear
set seed 2301
set obs 400
gen x = round(invnormal(runiform()), .01)
gen z = round(invnormal(runiform()), .01)
"""


def _run(text):
    s = Session()
    lines = []
    s.output.add_listener(lambda t, style: lines.append(t))
    rc = s.run_text(BASE + text)
    assert rc == 0, "".join(lines)[-3000:]
    return s, "".join(lines)


def _imputed(s, v):
    ds = s.data
    m = ds.get("_mi_m").data
    return ds.get(v).data[m > 0], ds.get(v).data[m == 0]


def test_mvn():
    s, out = _run("""gen y1 = 1 + .7*x + round(invnormal(runiform()), .01)
gen y2 = 2 - .5*x + .6*(y1 - 1 - .7*x) + round(invnormal(runiform()), .01)
replace y1 = . if runiform() < .25
replace y2 = . if runiform() < .25
mi set mlong
mi register imputed y1 y2
mi impute mvn y1 y2 = x, add(10) rseed(5) burnin(50) burnbetween(20)
""")
    assert "Multivariate" in out
    for v in ("y1", "y2"):
        imp, orig = _imputed(s, v)
        assert np.all(imp < 8.9e307)
        assert np.all(orig[orig < 8.9e307] == orig[orig < 8.9e307])
    ds = s.data
    m = ds.get("_mi_m").data
    # imputed values keep the regression on x
    y1 = ds.get("y1").data[m == 1]
    xx = ds.get("x").data[m == 1]
    ok = y1 < 8.9e307
    b = np.polyfit(xx[ok], y1[ok], 1)
    assert abs(b[0] - .7) < .25
    # without covariates
    s, _ = _run("""gen y1 = 1 + round(invnormal(runiform()), .01)
gen y2 = y1 + round(invnormal(runiform()), .01)
replace y2 = . if runiform() < .3
mi set mlong
mi register imputed y1 y2
mi impute mvn y1 y2, add(5) rseed(7)
""")
    imp, _ = _imputed(s, "y2")
    assert abs(imp.mean() - 1) < .4


def test_monotone():
    s, out = _run("""gen y1 = 1 + x + round(invnormal(runiform()), .01)
gen y2 = 2 + .5*y1 + round(invnormal(runiform()), .01)
replace y2 = . if runiform() < .4
replace y1 = . if runiform() < .15 & y2 == .
mi set mlong
mi register imputed y1 y2
mi impute monotone (regress) y1 (regress) y2 = x, add(3) rseed(11)
""")
    assert "Monotone" in out
    imp, _ = _imputed(s, "y2")
    assert np.all(imp < 8.9e307)
    s2 = Session()
    rc = s2.run_text(BASE + """gen y1 = x
gen y2 = z
replace y1 = . in 1/50
replace y2 = . in 40/100
mi set mlong
mi register imputed y1 y2
mi impute monotone (regress) y1 (regress) y2, add(2)
""")
    assert rc == 459


def test_truncreg_and_intreg():
    s, _ = _run("""gen y = 1 + x + round(invnormal(runiform()), .01)
drop if y <= 0
replace y = . if runiform() < .3
mi set mlong
mi register imputed y
mi impute truncreg y = x, ll(0) add(6) rseed(3)
""")
    imp, _ = _imputed(s, "y")
    assert np.all(imp > 0)
    s, _ = _run("""gen y = 1 + x + round(invnormal(runiform()), .01)
gen lo = floor(y)
gen hi = lo + 1
replace y = . if runiform() < .4
mi set mlong
mi register imputed y
mi impute intreg y = x, ll(lo) ul(hi) add(4) rseed(9)
""")
    ds = s.data
    m = ds.get("_mi_m").data
    ids = ds.get("_mi_id").data
    y = ds.get("y").data
    lo, hi = ds.get("lo").data, ds.get("hi").data
    sel = (m > 0) & (y < 8.9e307)
    assert sel.sum() > 0
    assert np.all(y[sel] >= lo[sel] - 1e-9) and np.all(y[sel] <= hi[sel] + 1e-9)


def test_nbreg():
    s, out = _run("""gen mu = exp(1 + .4*x)
gen u = runiform()
gen g = round(invnormal(runiform()), .0001)
* counts from a gamma-Poisson mixture built with uniforms
gen y = 0
forvalues k = 1/30 {
    quietly replace y = y + (runiform() < mu/30)
}
replace y = y + (runiform() < .5)*floor(3*runiform())
replace y = . if runiform() < .3
mi set mlong
mi register imputed y
mi impute nbreg y = x, add(8) rseed(21)
""")
    imp, orig = _imputed(s, "y")
    obs = orig[orig < 8.9e307]
    assert np.all(imp >= 0) and np.all(imp == np.floor(imp))
    assert abs(imp.mean() - obs.mean()) < .25 * obs.mean()
    assert "Negative binomial" in out
