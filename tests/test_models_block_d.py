"""Block 20b: gmm with igmm and derivatives(), ivpoisson cfunction with several endogenous
regressors, ivregress gmm with wmatrix(hac kernel opt)."""

import numpy as np

from opendta.session import Session

DATA = """clear
set seed 2020
set obs 400
gen t = _n
tsset t
gen x = rnormal()
gen z = rnormal()
gen w = rnormal()
gen u = rnormal()
gen e = rnormal()
gen e2 = rnormal()
gen e3 = rnormal()
gen yy = 1 + x + e
gen xe = x + .5*e + z + .3*w
gen w2 = x + z + .5*e2 + .3*w
gen w3 = x - z + .5*e3 + .3*u
gen yc = rpoisson(exp(.2 + .3*x + .2*w2 - .1*w3))
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


def test_gmm_igmm_matches_ivregress_igmm():
    s, out = _run(DATA + "gmm (yy - {b1}*xe - {b0}), instruments(z w) igmm\n")
    assert "Step 3" in out
    b = s.e["b"].data.ravel()
    s2, _ = _run("quietly ivregress gmm yy (xe = z w), igmm\n", s)
    assert np.allclose(b, s2.e["b"].data.ravel(), atol=1e-5)
    s3, out3 = _run("gmm (yy - {b1}*xe - {b0}), instruments(z w) igmm igmmiterate(1)\n", s)
    assert "Step 3" in out3 and "Step 4" not in out3


def test_gmm_derivatives_equal_numeric():
    s, _ = _run(DATA + "quietly gmm (yy - {b1}*xe - {b0}), instruments(z w)\n")
    b, V = s.e["b"].data.ravel().copy(), s.e["V"].data.copy()
    s, _ = _run("quietly gmm (yy - {b1}*xe - {b0}), instruments(z w) derivatives(/b1 = -xe) "
                "derivatives(/b0 = -1)\n", s)
    assert np.allclose(s.e["b"].data.ravel(), b, atol=1e-8) and np.allclose(s.e["V"].data, V, atol=1e-8)
    # a parameter without its derivative is an error
    rc = s.run_text("gmm (yy - {b1}*xe - {b0}), instruments(z w) derivatives(/b1 = -xe)\n")
    assert rc == 198


def test_ivpoisson_cfunction_two_endogenous():
    s, out = _run(DATA + "ivpoisson cfunction yc x (w2 w3 = z w u)\n")
    assert "/c_w2" in out and "/c_w3" in out and "Instrumented:  w2 w3" in out
    b = s.e["b"].data.ravel()
    assert len(b) == 16 and s.e["k"] == 16
    # exactly identified: the first-stage equations are OLS on the instruments
    Z = np.column_stack([_col(s, "x"), _col(s, "z"), _col(s, "w"), _col(s, "u"), np.ones(400)])
    pi = np.linalg.lstsq(Z, _col(s, "w2"), rcond=None)[0]
    assert np.allclose(b[4:9], pi, atol=1e-6)
    # the mean equation: the moments of the exponential mean with both residuals vanish
    v = [_col(s, n) - Z @ np.linalg.lstsq(Z, _col(s, n), rcond=None)[0] for n in ("w2", "w3")]
    X2 = np.column_stack([_col(s, "x"), _col(s, "w2"), _col(s, "w3"), np.ones(400), v[0], v[1]])
    be = np.array([b[0], b[1], b[2], b[3], b[14], b[15]])
    mu = np.exp(X2 @ be)
    assert np.abs((X2 * (_col(s, "yc") / mu - 1)[:, None]).mean(0)).max() < 1e-6


def test_hac_optimal_lag():
    s, _ = _run(DATA)
    for kern in ("bartlett", "parzen", "quadraticspectral"):
        s, out = _run(f"ivregress gmm yy (xe = z w), wmatrix(hac {kern} opt)\n", s)
        lag = int(s.e["hac_lag"])
        assert lag >= 0 and f"{lag} lags" in out
        s2, _ = _run(f"quietly ivregress gmm yy (xe = z w), wmatrix(hac {kern} {lag})\n", s)
        assert np.allclose(s2.e["b"].data, s.e["b"].data)
