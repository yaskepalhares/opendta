"""Block 20g: ivprobit and ivtobit with twostep (Newey's minimum chi-squared estimator), also with
more than one endogenous regressor."""

import numpy as np

from opendta.session import Session

DATA = """clear
set obs 3000
set seed 1926
gen x = round(invnormal(runiform()), .01)
gen z1 = round(invnormal(runiform()), .01)
gen z2 = round(invnormal(runiform()), .01)
gen v = round(invnormal(runiform()), .01)
gen v2 = round(invnormal(runiform()), .01)
gen u = round(.5*v + sqrt(.75)*invnormal(runiform()), .01)
gen y2 = round(.5*x + z1 + .5*z2 + v, .01)
gen y3 = round(-.3*x + .4*z1 + z2 + v2, .01)
gen ys = .5 + .4*y2 - .6*x + u
gen yb = ys > 0
gen yt = max(ys, 0)
gen ys3 = .5 + .4*y2 - .3*y3 - .6*x + u
gen yb3 = ys3 > 0
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


def test_ivprobit_twostep_estimates_scaled_coefficients():
    # two-step probit estimates delta / sqrt(1 - rho^2) (manual, Methods and formulas)
    s, _ = _run(DATA + "ivprobit yb x (y2 = z1 z2), nolog\n")
    b_ml = _b(s)[:3]
    rho = s.eval("e(rho)")
    _, out = _run("ivprobit yb x (y2 = z1 z2), twostep\n", s)
    assert np.allclose(_b(s), b_ml / np.sqrt(1 - rho ** 2), rtol=0.03)
    assert "Two-step probit with endogenous regressors" in out
    assert "Checking reduced-form model..." in out
    assert "Wald test of exogeneity:     chi2(1)" in out
    assert s.eval("e(df_exog)") == 1


def test_ivtobit_twostep_close_to_ml():
    s, _ = _run(DATA + "ivtobit yt x (y2 = z1 z2), ll(0) nolog\n")
    b_ml = _b(s)[:3]
    _, out = _run("ivtobit yt x (y2 = z1 z2), ll(0) twostep\n", s)
    assert np.allclose(_b(s), b_ml, rtol=0.01)
    assert "left-censored observations at yt <= 0" in out


def test_twostep_several_endogenous():
    s, out = _run(DATA + "ivprobit yb3 x (y2 y3 = z1 z2), twostep\n")
    assert list(s._estimates.names) == ["y2", "y3", "x", "_cons"]
    assert "chi2(2)" in out
    b = _b(s)
    assert b[0] > 0 and b[1] < 0 and b[2] < 0


def test_twostep_exogenous_regressor_matches_probit():
    # with an endogenous regressor that is in fact exogenous, the estimates are near probit
    s, _ = _run(DATA + "gen w = round(z1 + .5*z2 + invnormal(runiform()), .01)\ngen yw = (.3 + .5*w + u > 0)\n"
                "probit yw w\n")
    b_p = _b(s)
    _run("ivprobit yw (w = z1 z2), twostep\n", s)
    assert np.allclose(_b(s), b_p, atol=0.15)


def test_twostep_errors():
    s = Session()
    _run(DATA, s)
    assert s.run_text("ivprobit yb x (y2 y3 = z1), twostep\n") == 481
    assert s.run_text("ivprobit yb x (y2 = z1 z2), twostep vce(robust)\n") == 198
    assert s.run_text("ivprobit yb x (y2 y3 = z1 z2)\n") == 198
