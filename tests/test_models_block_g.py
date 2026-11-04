"""Block 20e: slogit with dimension() above 1; exlogistic with binomial(), group(), condvars() and
estconstant; expoisson with exposure(), offset() and condvars()."""

import numpy as np

from opendta.session import Session

SLOGIT = """clear
set obs 1500
set seed 1924
gen x1 = round(invnormal(runiform()), .01)
gen x2 = round(invnormal(runiform()), .01)
gen x3 = round(invnormal(runiform()), .01)
gen double u1 = runiform()
gen double e1 = .5*x1 + .3*x2
gen double e2 = .4*x2 - .6*x3
gen double s1 = exp(.5 - e1)
gen double s2 = exp(.2 - e2)
gen double s3 = exp(-.1 - .5*e1 - .5*e2)
gen double tot = 1 + s1 + s2 + s3
gen y = 4
replace y = 1 if u1 < s1/tot
replace y = 2 if u1 >= s1/tot & u1 < (s1+s2)/tot
replace y = 3 if u1 >= (s1+s2)/tot & u1 < (s1+s2+s3)/tot
"""

EXACT = """clear
set obs 40
set seed 1924
gen g = ceil(_n/8)
gen x1 = round(runiform())
gen x2 = round(runiform()*2)
gen n = 1 + round(runiform()*3)
gen y = (runiform() < .3 + .3*x1)
gen yb = round(n*(.2 + .4*x1*runiform()))
gen t = 1 + round(runiform()*2)
gen c = round(t*exp(.2*x1)*runiform()*2)
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


def test_slogit_dim2_constraints_and_fit():
    s, out = _run(SLOGIT + "slogit y x1 x2 x3, dim(2)\n")
    for line in (" ( 1)  [phi1_1]_cons = 1", " ( 2)  [phi1_2]_cons = 0", " ( 3)  [phi2_1]_cons = 0",
                 " ( 4)  [phi2_2]_cons = 1"):
        assert line in out
    assert "dim1" in out and "dim2" in out and "/phi2_3" in out
    assert s.eval("e(k_dim)") == 2.0
    ll2 = s._estimates.extra["ll"]
    # the maximum: zero score, and dim(2) fits at least as well as dim(1)
    _run("slogit y x1 x2 x3\n", s)
    assert s._estimates.extra["ll"] <= ll2 + 1e-6
    # dim(3) with four outcomes is the multinomial logit
    _run("slogit y x1 x2 x3, dim(3)\n", s)
    ll3 = s._estimates.extra["ll"]
    _run("mlogit y x1 x2 x3\n", s)
    assert abs(ll3 - s._estimates.extra["ll"]) < 1e-4


def test_slogit_dim_too_large():
    s = Session()
    _run(SLOGIT, s)
    assert s.run_text("slogit y x1 x2 x3, dim(4)\n") == 198


def test_expoisson_exposure_brute_force():
    # conditional likelihood of b by enumeration: Pr(y | sum y) is proportional to
    # prod (t_i exp(b x_i))^y_i / y_i!
    import itertools
    from math import factorial
    from scipy import optimize
    x = np.array([0, 1, 0, 1, 1, 0.])
    tt = np.array([1, 2, 3, 1, 2, 2.])
    c = np.array([1, 2, 0, 1, 2, 0])
    total = int(c.sum())
    vecs = [v for v in itertools.product(range(total + 1), repeat=6) if sum(v) == total]
    T = np.array([np.dot(v, x) for v in vecs])
    lw = np.array([sum(vi * np.log(ti) - np.log(factorial(vi)) for vi, ti in zip(v, tt)) for v in vecs])
    t_obs = float(c @ x)

    def nll(b):
        return -(b * t_obs - np.log(np.sum(np.exp(lw + b * T))))
    b_ref = optimize.minimize_scalar(nll, bounds=(-10, 10), method="bounded", options={"xatol": 1e-10}).x
    data = "clear\ninput x t c\n" + "".join(f"{a:g} {b:g} {d:g}\n" for a, b, d in zip(x, tt, c)) + "end\n"
    s, _ = _run(data + "expoisson c x, exposure(t)\n")
    assert abs(_b(s)[0] - b_ref) < 1e-5
    _run("gen lt = ln(t)\nexpoisson c x, offset(lt)\n", s)
    assert abs(_b(s)[0] - b_ref) < 1e-5


def test_expoisson_condvars_absorbs_offset():
    s, _ = _run(EXACT + "expoisson c x1 x2\n")
    b = _b(s)
    _, out = _run("expoisson c x1, condvars(x2) offset(x2)\n", s)
    assert abs(_b(s)[0] - b[0]) < 1e-8
    # x2 only appears as the offset row (compat 1924)
    assert "x2 |          1  (offset)" in out and out.count("x2 |") == 1


def test_exlogistic_binomial_equals_expanded_data():
    s, _ = _run(EXACT + "exlogistic yb x1 x2, binomial(n)\n")
    b = _b(s)
    _run("gen id = _n\nexpand n\nbysort id: gen yy = (_n <= yb)\nexlogistic yy x1 x2\n", s)
    assert np.allclose(_b(s), b, atol=1e-6)


def test_exlogistic_group_condvars_estconstant():
    s, out = _run(EXACT + "exlogistic y x1 x2, group(g)\n")
    assert list(s._estimates.names) == ["x1", "x2"]
    # one group: the same as no group
    s2, _ = _run(EXACT + "replace g = 1\nexlogistic y x1 x2, group(g)\n")
    b_g = _b(s2)
    _run("exlogistic y x1 x2\n", s2)
    assert np.allclose(_b(s2), b_g, atol=1e-8)
    _, out = _run("exlogistic y x1 x2, estconstant coef\n", s2)
    assert "_cons" in out
    # with odds ratios the constant has no row (compat 1924)
    _, out = _run("exlogistic y x1 x2, estconstant\n", s2)
    assert "_cons" not in out
    _run("exlogistic y x1, condvars(x2)\n", s2)
    assert np.allclose(_b(s2), b_g[:1], atol=1e-8)


def test_exact_option_errors():
    s = Session()
    _run(EXACT, s)
    assert s.run_text("expoisson c x1, binomial(n)\n") == 198
    assert s.run_text("exlogistic y x1, group(g) estconstant\n") == 198
    assert s.run_text("exlogistic y x1, saving(dist)\n") == 198
