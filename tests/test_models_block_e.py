"""Block 20c: biprobit with partial and constraints(), frontier with constraints()."""

import numpy as np
from scipy.stats import multivariate_normal

from opendta.session import Session

DATA = """clear
set obs 2000
set seed 3520
gen x1 = round(invnormal(runiform()), .01)
gen x2 = round(invnormal(runiform()), .01)
gen x3 = round(invnormal(runiform()), .01)
gen u1 = invnormal(runiform())
gen u2 = .3*u1 + sqrt(.91)*invnormal(runiform())
gen y1 = (.8 + x2 + .3*x1 + u1 > 0)
gen y2 = (.5 - x3 + .3*x1 + u2 > 0)
gen e = round(invnormal(runiform())*.3, .01) - abs(round(invnormal(runiform())*.4, .01))
gen yf = 1 + .5*x1 + .3*x2 + e
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


def test_biprobit_constraints_hold_and_are_listed():
    s, out = _run(DATA + "constraint 1 [y1]x1 = [y2]x1\nbiprobit (y1 = x1 x2) (y2 = x1 x3), constraints(1)\n")
    b = _b(s)
    assert abs(b[0] - b[3]) < 1e-8
    assert " ( 1)  [y1]x1 - [y2]x1 = 0" in out
    assert "Wald chi2(3)" in out
    _, out2 = _run("matrix list e(Cns)\n", s)
    assert "_Cns" in out2


def test_biprobit_constraints_match_free_fit_when_not_binding():
    # a constraint that the free fit already satisfies leaves the estimates where they are
    s, _ = _run(DATA + "biprobit (y1 = x1 x2) (y2 = x1 x3)\n")
    b_free = _b(s)
    _run(f"constraint 2 [y2]x3 = {float(b_free[4])!r}\n", s)
    _run("biprobit (y1 = x1 x2) (y2 = x1 x3), constraints(2)\n", s)
    assert np.allclose(_b(s), b_free, atol=1e-5)


def test_biprobit_partial_likelihood():
    s, out = _run(DATA + "biprobit (y1 = x1 x2) (y2 = x1 x3), partial\n")
    assert "Partial observability bivariate probit" in out
    b = _b(s)
    _run("predict double e1, xb1\npredict double e2, xb2\n", s)
    e1 = s.data.get("e1").data.astype(float)
    e2 = s.data.get("e2").data.astype(float)
    z = s.data.get("y1").data.astype(float) * s.data.get("y2").data.astype(float)
    rho = np.tanh(b[-1])
    mvn = multivariate_normal(mean=[0, 0], cov=[[1, rho], [rho, 1]])
    P = mvn.cdf(np.column_stack([e1, e2]))
    ll = float(np.sum(np.where(z > 0, np.log(P), np.log1p(-P))))
    assert abs(ll - float(s._estimates.extra["ll"])) < 1e-4
    # the maximum: the score is zero
    assert abs(b[1] - 1.0) < 0.25 and abs(b[4] + 1.0) < 0.25


def test_frontier_constraints():
    s, out = _run(DATA + "constraint 3 x1 + x2 = .8\nfrontier yf x1 x2, constraints(3)\n")
    b = _b(s)
    assert abs(b[0] + b[1] - .8) < 1e-8
    assert " ( 1)  [yf]x1 + [yf]x2 = .8" in out


def test_frontier_constraint_on_variance():
    s, out = _run(DATA + "constraint 4 [lnsig2v]_cons = -2\nfrontier yf x1 x2, constraints(4)\n")
    assert abs(_b(s)[3] + 2) < 1e-8
    assert "[lnsig2v]_cons = -2" in out
