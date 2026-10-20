"""Block 5: meglm families, mecloglog, menbreg, meologit/meoprobit, meqr*, mestreg and crossed
effects in mixed. Each model is checked against the likelihood integrated with a 60-node
Gauss-Hermite rule (independent code) at the estimates, and the estimates against a direct
maximization of that likelihood."""

import numpy as np
from scipy import optimize, special, stats

from opendta.session import Session

DATA = """clear
set seed 2024
set obs 40
gen g = _n
gen u = rnormal()*.7
expand 12
gen x = rnormal()
gen mu = exp(.3 + .5*x + u*.6)
gen yg = .3 + .5*x + u + rnormal()*1.2
gen yb = runiform() < 1 - exp(-exp(-.2 + .6*x + u))
gen yn = rpoisson(mu*rgamma(2, .5))
gen lat = .8*x + u + logit(runiform())
gen yo = 1 + (lat > -.7) + (lat > .8)
"""


def _run(text, s=None):
    s = s or Session()
    lines = []
    s.output.add_listener(lambda t, style: lines.append(t))
    rc = s.run_text(text)
    assert rc == 0, "".join(lines)[-3000:]
    return s, "".join(lines)


def _col(s, name):
    return s.data.get(name).data.astype(float)


def _exact(s, lf_of_eta, npar_aux, yname="yg"):
    """ll(par), par = [_cons, x, aux..., ln sd]; lf_of_eta(eta N x K, aux, y)."""
    g = _col(s, "g").astype(int) - 1
    x = _col(s, "x")
    y = _col(s, yname)
    xn, w = np.polynomial.hermite_e.hermegauss(60)
    w = w / np.sqrt(2 * np.pi)

    def ll(par):
        eta = par[0] + par[1] * x[:, None] + np.exp(par[-1]) * xn[None, :]
        lf = lf_of_eta(eta, par[2:2 + npar_aux], y[:, None])
        S = np.zeros((g.max() + 1, len(xn)))
        np.add.at(S, g, lf)
        return float(np.sum(special.logsumexp(S + np.log(w)[None, :], axis=1)))
    return ll


def _check(s, f, tol=2e-4):
    b = s.e["b"].data.ravel()
    par = np.concatenate([[b[1], b[0]], b[2:]])
    assert abs(s.e["ll"] - f(par)) < tol
    r = optimize.minimize(lambda p: -f(p), par, method="BFGS", options={"gtol": 1e-6})
    assert np.max(np.abs(r.x - par)) < 5e-3


def test_meglm_gaussian_matches_mixed():
    s, _ = _run(DATA + "mixed yg x || g:")
    ll_mixed = s.e["ll"]
    s, out = _run("meglm yg x || g:", s)
    assert abs(s.e["ll"] - ll_mixed) < 1e-5
    assert "Gaussian" in out and "Identity" in out
    f = _exact(s, lambda eta, a, y: stats.norm.logpdf(y, eta, np.exp(a[0] / 2)), 1)
    _check(s, f)


def test_mecloglog():
    s, out = _run(DATA + "mecloglog yb x || g:")
    f = _exact(s, lambda eta, a, y: y * np.log(-np.expm1(-np.exp(eta))) + (1 - y) * (-np.exp(eta)), 0, "yb")
    _check(s, f)


def test_menbreg_mean_and_constant():
    s, _ = _run(DATA + "menbreg yn x || g:")
    f = _exact(s, lambda eta, a, y: stats.nbinom.logpmf(
        y, 1 / np.exp(a[0]), 1 / (1 + np.exp(a[0]) * np.exp(eta))), 1, "yn")
    _check(s, f)
    s, _ = _run(DATA + "menbreg yn x || g:, dispersion(constant)")
    f = _exact(s, lambda eta, a, y: stats.nbinom.logpmf(
        y, np.exp(eta) / np.exp(a[0]), 1 / (1 + np.exp(a[0]))), 1, "yn")
    _check(s, f)


def test_meologit_meoprobit():
    for cmd, F in (("meologit", special.expit), ("meoprobit", stats.norm.cdf)):
        s, _ = _run(DATA + f"{cmd} yo x || g:")
        b = s.e["b"].data.ravel()               # x, /cut1, /cut2, ln sd
        g = _col(s, "g").astype(int) - 1
        x, yi = _col(s, "x"), _col(s, "yo").astype(int) - 1
        xn, w = np.polynomial.hermite_e.hermegauss(60)
        w = w / np.sqrt(2 * np.pi)
        eta = b[0] * x[:, None] + np.exp(b[3]) * xn[None, :]
        cuts = np.concatenate([[-np.inf], b[1:3], [np.inf]])
        p = F(cuts[yi + 1][:, None] - eta) - F(cuts[yi][:, None] - eta)
        S = np.zeros((g.max() + 1, len(xn)))
        np.add.at(S, g, np.log(p))
        assert abs(s.e["ll"] - np.sum(special.logsumexp(S + np.log(w)[None, :], axis=1))) < 2e-4


def test_meqr_equals_me():
    s, _ = _run(DATA + "melogit yb x || g:")
    ll1 = s.e["ll"]
    s, _ = _run("meqrlogit yb x || g:", s)
    assert abs(s.e["ll"] - ll1) < 1e-5
    s, _ = _run("meqrpoisson yn x || g:", s)
    s2, _ = _run("mepoisson yn x || g:", s)
    assert abs(s.e["ll"] - s2.e["ll"]) < 1e-5


def test_mixed_crossed_effects_match_dense_ml():
    text = """clear
set seed 77
set obs 150
gen a = ceil(runiform()*8)
gen b = ceil(runiform()*6)
gen x = rnormal()
gen ua = .8*sin(a*7)
gen ub = .5*cos(b*5)
gen y = 1 + .6*x + ua + ub + rnormal()
"""
    s, out = _run(text + "mixed y x || _all: R.a || _all: R.b")
    y, x = _col(s, "y"), _col(s, "x")
    a, b = _col(s, "a").astype(int) - 1, _col(s, "b").astype(int) - 1
    n = len(y)
    X = np.column_stack([np.ones(n), x])
    Za, Zb = np.eye(a.max() + 1)[a], np.eye(b.max() + 1)[b]

    def nll(p):
        v = np.exp(2 * p)
        V = v[0] * Za @ Za.T + v[1] * Zb @ Zb.T + v[2] * np.eye(n)
        Vi = np.linalg.inv(V)
        be = np.linalg.solve(X.T @ Vi @ X, X.T @ Vi @ y)
        e = y - X @ be
        return .5 * (n * np.log(2 * np.pi) + np.linalg.slogdet(V)[1] + e @ Vi @ e)
    o = optimize.minimize(nll, [0, 0, 0], method="Nelder-Mead",
                          options=dict(xatol=1e-8, fatol=1e-10, maxiter=4000))
    assert abs(s.e["ll"] + o.fun) < 1e-4
    assert "Group variable: _all" in out and "var(R.a)" in out and "var(R.b)" in out
    s, _ = _run("mixed y x || _all: R.a || b:", s)
    assert abs(s.e["ll"] + o.fun) < 1e-4


def test_mestreg_weibull_matches_exact_likelihood():
    text = """clear
set seed 31
set obs 40
gen g = _n
gen u = rnormal()*.6
expand 10
gen x = rnormal()
gen t = (-ln(runiform())/exp(-1 + .5*x + u))^(1/1.4)
gen d = 1
stset t, failure(d)
"""
    s, out = _run(text + "mestreg x || g:, distribution(weibull) nohr")
    b = s.e["b"].data.ravel()               # x, _cons, /ln_p, ln sd
    g = _col(s, "g").astype(int) - 1
    x, t = _col(s, "x"), _col(s, "t")
    xn, w = np.polynomial.hermite_e.hermegauss(60)
    w = w / np.sqrt(2 * np.pi)
    p = np.exp(b[2])
    eta = b[1] + b[0] * x[:, None] + np.exp(b[3]) * xn[None, :]
    lf = eta + np.log(p) + (p - 1) * np.log(t)[:, None] - np.exp(eta) * (t ** p)[:, None]
    S = np.zeros((g.max() + 1, len(xn)))
    np.add.at(S, g, lf)
    assert abs(s.e["ll"] - np.sum(special.logsumexp(S + np.log(w)[None, :], axis=1))) < 2e-4
    assert "Mixed-effects Weibull regression" in out


def test_nonadaptive_quadrature_converges_to_the_adaptive_value():
    s, _ = _run(DATA + "mecloglog yb x || g:, intpoints(30)")
    ll_ad = s.e["ll"]
    s, out = _run("mecloglog yb x || g:, intmethod(ghermite) intpoints(30)", s)
    assert abs(s.e["ll"] - ll_ad) < 1e-4
    assert out.count("Iteration") < 40


def test_stddeviations_and_conditional_fixedonly():
    s, out = _run(DATA + "meglm yg x || g:, stddeviations")
    assert "sd(_cons)" in out and "sd(e.yg)" in out
    s, _ = _run("predict m1, mu conditional(fixedonly)\npredict m2, fixedonly\npredict m3, xb", s)
    assert np.allclose(_col(s, "m1"), _col(s, "m2")) and np.allclose(_col(s, "m1"), _col(s, "m3"))
