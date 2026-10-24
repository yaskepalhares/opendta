"""Block 20a: tpoisson/tnbreg ul() and variable limits, tnbreg and nbreg dispersion(constant),
glogit or, offset() in intreg, truncreg, heckman and biprobit, heckman mills(), sureg and reg3
constraints(), ivregress gmm wmatrix(hac), center and igmm, frontier uhet(), vhet() and cm(),
glm and binreg mu()."""

import math

import numpy as np
from scipy import special, stats

from opendta.session import Session

DATA = """clear
set seed 2018
set obs 400
gen t = _n
gen x = rnormal()
gen z = rnormal()
gen o = .3*runiform()
gen double lam = exp(.8 + .4*x)
gen y = rpoisson(lam)
gen y1 = 1 + .5*x + .5*z + rnormal()
gen y2 = 2 + .5*x - z + rnormal()
"""


def _run(text, s=None):
    s = s or Session()
    lines = []
    s.output.add_listener(lambda t, style: lines.append(t))
    rc = s.run_text(text)
    assert rc == 0, "".join(lines)[-3000:]
    return s, "".join(lines)


def _col(s, name, mask=None):
    v = s.data.get(name).data.astype(float)
    return v if mask is None else v[mask]


def test_tpoisson_upper_limit_likelihood():
    s, out = _run(DATA + "tpoisson y x if y > 0 & y < 6, ul(6)\n")
    assert "Truncation points: 0 - 6" in out
    b = s.e["b"].data.ravel()
    m = (_col(s, "y") > 0) & (_col(s, "y") < 6)
    y, x = _col(s, "y", m), _col(s, "x", m)
    mu = np.exp(b[0] * x + b[1])
    ll = stats.poisson.logpmf(y, mu) - np.log(stats.poisson.cdf(5, mu) - stats.poisson.cdf(0, mu))
    assert abs(ll.sum() - s.e["ll"]) < 1e-6
    # the mean of the truncated distribution
    s, _ = _run("predict n\n", s)
    n = _col(s, "n", m)
    j = np.arange(1, 6)
    pm = stats.poisson.pmf(j[None, :], mu[:, None])
    assert np.allclose(n, (pm * j).sum(1) / pm.sum(1))


def test_nb1_nbreg_and_tnbreg():
    s, out = _run(DATA + "nbreg y x, dispersion(constant)\n")
    assert "Dispersion     = constant" in out and "/lndelta" in out
    b = s.e["b"].data.ravel()
    y, x = _col(s, "y"), _col(s, "x")
    mu = np.exp(b[0] * x + b[1])
    d = math.exp(b[2])
    r = mu / d
    ll = special.gammaln(y + r) - special.gammaln(r) - special.gammaln(y + 1) - r * np.log1p(d) + y * (np.log(d) - np.log1p(d))
    assert abs(ll.sum() - s.e["ll"]) < 1e-6
    s, out = _run(DATA + "tnbreg y x if y > 0, dispersion(constant)\n")
    b = s.e["b"].data.ravel()
    m = _col(s, "y") > 0
    y, x = _col(s, "y", m), _col(s, "x", m)
    mu = np.exp(b[0] * x + b[1])
    d = math.exp(b[2])
    r = mu / d
    lf = special.gammaln(y + r) - special.gammaln(r) - special.gammaln(y + 1) - r * np.log1p(d) + y * (np.log(d) - np.log1p(d))
    p0 = np.exp(-r * np.log1p(d))
    assert abs((lf - np.log(1 - p0)).sum() - s.e["ll"]) < 1e-6


def test_variable_lower_limit():
    s, out = _run(DATA + "gen l = cond(x > 0, 1, 0)\ntpoisson y x if y > l, ll(l)\n")
    assert "Truncation point: l" in out
    b = s.e["b"].data.ravel()
    m = _col(s, "y") > _col(s, "l")
    y, x, lv = _col(s, "y", m), _col(s, "x", m), _col(s, "l", m)
    mu = np.exp(b[0] * x + b[1])
    ll = stats.poisson.logpmf(y, mu) - stats.poisson.logsf(lv, mu)
    assert abs(ll.sum() - s.e["ll"]) < 1e-6


def test_glogit_or():
    s, out = _run(DATA + "gen pop = y + 3\nglogit y pop x in 1/80, or\n")
    assert "Odds Ratio" in out
    b = s.e["b"].data.ravel()
    assert f"{math.exp(b[0]):.6g}"[:5] in out.replace(" ", "")


def test_offsets():
    s, out = _run(DATA + "gen lo = floor(y1)\ngen hi = lo + 1\nquietly intreg lo hi x z, offset(o)\n")
    b_off = s.e["b"].data.ravel().copy()
    s, out = _run("gen lo2 = lo - o\ngen hi2 = hi - o\nquietly intreg lo2 hi2 x z\nintreg lo hi x z, offset(o)\n", s)
    assert np.allclose(s.e["b"].data.ravel(), b_off, atol=1e-6)
    assert "(offset)" in out
    s, out = _run(DATA + "gen s = (x + rnormal()) > 0\ngen yh = y1 if s\n"
                  "heckman yh x, select(s = x z) twostep mills(imr)\nprobit s x z\npredict xbs, xb\n")
    imr = _col(s, "imr")
    xbs = _col(s, "xbs")
    assert np.allclose(imr, stats.norm.pdf(xbs) / stats.norm.cdf(xbs), atol=1e-6)
    s, out = _run("gen b1 = y1 > 1\ngen b2 = y2 > 2\nbiprobit b1 b2 x z, offset1(o)\n", s)
    assert "(offset)" in out


def test_system_constraints():
    s, out = _run(DATA + "constraint define 1 [y1]x = [y2]x\nsureg (y1 x z) (y2 x z), constraints(1)\n")
    b = s.e["b"].data.ravel()
    assert abs(b[0] - b[3]) < 1e-10
    assert "( 1)  [y1]x - [y2]x = 0" in out
    # the constrained GLS equals the unconstrained GLS of the model rewritten with a common slope
    s2, _ = _run("quietly reg3 (y1 x z) (y2 x z), constraints(1)\n", s)
    assert np.allclose(s2.e["b"].data.ravel(), b, atol=1e-8)


def test_gmm_options():
    base = DATA + """gen e = rnormal()
gen w = rnormal()
gen yy = 1 + x + e
gen xe = x + .5*e + z + .3*w
tsset t
"""
    s, _ = _run(base + "quietly ivregress gmm yy (xe = z w)\n")
    b_rob, V_rob = s.e["b"].data.ravel().copy(), s.e["V"].data.copy()
    s, out = _run("ivregress gmm yy (xe = z w), wmatrix(hac bartlett 0)\n", s)
    assert np.allclose(s.e["b"].data.ravel(), b_rob) and np.allclose(s.e["V"].data, V_rob)
    s, out = _run("ivregress gmm yy (xe = z w), wmatrix(hac bartlett 3)\n", s)
    assert "HAC Bartlett 3" in out
    s, _ = _run("quietly ivregress gmm yy (xe = z w), igmm\n", s)
    b_it = s.e["b"].data.ravel()
    # the iterated estimate is a fixed point of the weight matrix
    X = np.column_stack([_col(s, "xe"), np.ones(400)])
    Z = np.column_stack([_col(s, "z"), _col(s, "w"), np.ones(400)])
    y = _col(s, "yy")
    e = y - X @ np.array([b_it[0], b_it[1]])
    W = np.linalg.inv((Z * e[:, None] ** 2).T @ Z)
    ZX = Z.T @ X
    b_fp = np.linalg.solve(ZX.T @ W @ ZX, ZX.T @ W @ (Z.T @ y))
    assert np.allclose(b_fp, b_it, atol=1e-5)
    s, _ = _run("quietly ivregress gmm yy (xe = z w), center\n", s)
    assert not np.array_equal(s.e["V"].data, V_rob)


def test_frontier_heteroskedastic():
    base = DATA + "gen w = runiform()\ngen u = abs(rnormal())*exp(.5*w)\ngen yf = 1 + .5*x + .3*rnormal() - .5*u\n"
    s, _ = _run(base + "quietly frontier yf x\n")
    ll0 = s.e["ll"]
    s, out = _run("frontier yf x, uhet(w)\n", s)
    ll1 = s.e["ll"]
    assert ll1 >= ll0 - 1e-6 and "lnsig2u" in out and "sigma_v" in out
    s, _ = _run("quietly frontier yf x, uhet(w) vhet(w)\n", s)
    assert s.e["ll"] >= ll1 - 1e-6
    s, out = _run("frontier yf x, distribution(tnormal) cm(w)\n", s)
    assert "mu" in out
    assert Session().run_text(base + "frontier yf x, cm(w)\n") == 198


def test_glm_mu_start():
    s, _ = _run(DATA + "gen yb = y1 > 1\nquietly binreg yb x\n")
    b0 = s.e["b"].data.ravel().copy()
    s, _ = _run("gen m0 = .5\nquietly binreg yb x, mu(m0)\n", s)
    assert np.allclose(s.e["b"].data.ravel(), b0, atol=1e-6)
