"""Block 20h: eteffects by the manual's control function (residual t - Phi(z'pi)), with linear,
probit, fractional and exponential outcome models."""

import numpy as np
from scipy import stats

from opendta.session import Session

DATA = """clear
set obs 2000
set seed 1927
gen x = round(invnormal(runiform()), .01)
gen z = round(invnormal(runiform()), .01)
gen u = round(invnormal(runiform()), .01)
gen t = (.2 + .7*z + .5*x + u) > 0
gen yl = round(1 + x + 1.5*t + .6*u + .7*invnormal(runiform()), .01)
gen yb = (-.2 + .5*x + .8*t + .5*u + invnormal(runiform()) > 0)
gen yf = invlogit(-.2 + .5*x + .8*t + .5*u + invnormal(runiform()))
gen yc = round(exp(.2 + .3*x + .5*t + .3*u)*runiform()*2)
"""


def _run(text, s=None):
    s = s or Session()
    lines = []
    s.output.add_listener(lambda t, style: lines.append(t))
    rc = s.run_text(text)
    assert rc == 0, "".join(lines)[-3000:]
    return s, "".join(lines)


def _b(s):
    return np.asarray(s._estimates.b, dtype=float)


def _col(s, nm):
    return s.data.get(nm).data.astype(float)


def _nu(s):
    _run("quietly probit t x z\nquietly predict double ph, pr\n", s)
    return _col(s, "t") - _col(s, "ph")


def test_linear_pomeans_by_hand():
    s, _ = _run(DATA)
    nu = _nu(s)
    y, x, t = _col(s, "yl"), _col(s, "x"), _col(s, "t")
    W = np.column_stack([x, np.ones_like(x), nu])
    pom = []
    for lv in (0, 1):
        a = np.linalg.lstsq(W[t == lv], y[t == lv], rcond=None)[0]
        pom.append(float((W @ a).mean()))
    _run("eteffects (yl x) (t x z), pomeans\n", s)
    assert np.allclose(_b(s)[:2], pom, atol=1e-6)


def test_exponential_pomeans_by_hand():
    s, _ = _run(DATA)
    nu = _nu(s)
    y, x, t = _col(s, "yc"), _col(s, "x"), _col(s, "t")
    W = np.column_stack([x, np.ones_like(x), nu])
    pom = []
    for lv in (0, 1):
        m = t == lv
        a = np.zeros(3)
        for _ in range(50):
            mu = np.exp(W[m] @ a)
            a = a + np.linalg.solve((W[m] * mu[:, None]).T @ W[m], W[m].T @ (y[m] - mu))
        pom.append(float(np.exp(W @ a).mean()))
    _, out = _run("eteffects (yc x, exponential) (t x z), pomeans\n", s)
    assert np.allclose(_b(s)[:2], pom, atol=1e-6)
    assert "Outcome model  : exponential" in out


def test_probit_and_fractional_outcomes():
    s, _ = _run(DATA)
    nu = _nu(s)
    y, x, t = _col(s, "yb"), _col(s, "x"), _col(s, "t")
    W = np.column_stack([x, np.ones_like(x), nu])
    pom = []
    for lv in (0, 1):
        m = t == lv
        a = np.zeros(3)
        for _ in range(50):
            e = W[m] @ a
            F = stats.norm.cdf(e)
            f = stats.norm.pdf(e)
            g = W[m].T @ ((y[m] - F) * f / (F * (1 - F)))
            H = (W[m] * (f * f / (F * (1 - F)))[:, None]).T @ W[m]
            a = a + np.linalg.solve(H, g)
        pom.append(float(stats.norm.cdf(W @ a).mean()))
    _run("eteffects (yb x, probit) (t x z)\n", s)
    assert np.allclose(_b(s)[:2], [pom[1] - pom[0], pom[0]], atol=1e-6)
    _, out = _run("eteffects (yf x, fractional) (t x z)\n", s)
    assert "Outcome model  : fractional" in out


def test_outcome_checks():
    s = Session()
    _run(DATA, s)
    assert s.run_text("eteffects (yl x, probit) (t x z)\n") != 0
    assert s.run_text("eteffects (yl x) (t x z, logit)\n") == 198
