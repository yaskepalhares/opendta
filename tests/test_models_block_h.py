"""Block 20f: teffects ra, ipw, ipwra and aipw with multivalued treatments."""

import numpy as np

from opendta.core import treatfx as T
from opendta.session import Session

DATA = """clear
set obs 1200
set seed 1925
gen x = round(invnormal(runiform()), .01)
gen z = round(invnormal(runiform()), .01)
gen double u = runiform()
gen double p1 = exp(.3*x + .2*z)
gen double p2 = exp(-.2 + .5*x - .3*z)
gen t = 0
replace t = 1 if u < p1/(1+p1+p2)
replace t = 2 if u >= p1/(1+p1+p2) & u < (p1+p2)/(1+p1+p2)
gen y = round(1 + .5*x + t + .3*t*x + invnormal(runiform()), .01)
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


def _data(s, *names):
    return [s.data.get(nm).data.astype(float) for nm in names]


def test_ra_pomeans_are_arm_regression_means():
    s, out = _run(DATA + "teffects ra (y x) (t), pomeans\n")
    y, x, t = _data(s, "y", "x", "t")
    X = np.column_stack([x, np.ones_like(x)])
    ref = []
    for lv in (0, 1, 2):
        b = np.linalg.lstsq(X[t == lv], y[t == lv], rcond=None)[0]
        ref.append(float((X @ b).mean()))
    assert np.allclose(_b(s), ref, atol=1e-8)
    assert list(s._estimates.names) == ["0.t", "1.t", "2.t"]


def test_ate_rows_and_control():
    s, out = _run(DATA + "teffects ipwra (y x) (t x z)\n")
    assert "(1 vs 0)" in out and "(2 vs 0)" in out and "(multinomial) logit" in out
    b = _b(s)
    assert abs(b[0] - 1.0) < .2 and abs(b[1] - 2.0) < .2
    _, out = _run("teffects ra (y x) (t), control(1)\n", s)
    assert "(0 vs 1)" in out and "(2 vs 1)" in out
    assert s.eval("e(n_levels)") == 3 and s.eval("e(control)") == 1


def test_multi_with_two_levels_matches_binary():
    s, _ = _run(DATA + "keep if t < 2\n")
    y, x, z, t = _data(s, "y", "x", "z", "t")
    n = len(y)
    Xo = np.column_stack([x, np.ones(n)])
    Z = np.column_stack([x, z, np.ones(n)])
    for est in ("ra", "ipw", "ipwra", "aipw"):
        th2, V2, _ = T.teffects(est, y, t, Xo if est != "ipw" else np.ones((n, 1)), Z, "linear", "logit", "ate")
        thm, Vm = T.teffects_multi(est, y, t.astype(int), 2, Xo if est != "ipw" else np.ones((n, 1)), Z,
                                   "linear", "ate", 1)
        assert np.allclose(th2[:2], thm[:2], atol=1e-7)
        assert np.allclose(V2[:2, :2], Vm[:2, :2], rtol=1e-4, atol=1e-10)


def test_multivalued_errors():
    s = Session()
    _run(DATA, s)
    assert s.run_text("teffects ipw (y) (t x z, probit)\n") == 198
    assert s.run_text("teffects nnmatch (y x) (t)\n") == 198
