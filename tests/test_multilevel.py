"""Multilevel models: mixed (checked against statsmodels MixedLM) and
melogit/mepoisson (checked against the likelihood integrated with
scipy.integrate.quad and maximized directly)."""

import numpy as np

from opendta.session import Session

LMM = """clear
set seed 1001
set obs 30
gen id = _n
gen u0 = rnormal()*2
gen u1 = rnormal()*.5
expand 8
bysort id: gen week = _n
gen x = runiform()
gen y = 10 + u0 + (1 + u1)*week + 2*x + rnormal()*1.5
"""

GLMM = """clear
set seed 1003
set obs 40
gen g = _n
gen u = rnormal()*.9
gen v = rnormal()*.3
expand 15
gen x = rnormal()
gen yb = runiform() < invlogit(-.4 + .8*x + u + v*x)
gen yc = rpoisson(exp(.2 + .5*x + u/2))
"""


def _run(text):
    s = Session()
    lines = []
    s.output.add_listener(lambda t, style: lines.append(t))
    rc = s.run_text(text)
    assert rc == 0, "".join(lines)[-3000:]
    return s, "".join(lines)


def _b(s):
    return s.e["b"].data.ravel()


def test_mixed_random_intercept_ml_and_reml():
    s, out = _run(LMM + "mixed y week x || id:")
    b = _b(s)
    assert abs(s.e["ll"] - (-534.7523279263447)) < 1e-6
    assert np.allclose(b[:3], [1.02968492, 2.34329639, 10.12384973], atol=1e-5)
    assert np.exp(2 * b[3]) == np.exp(2 * b[3]) and abs(np.exp(2 * b[3]) - 7.8228) < 2e-3
    assert "Mixed-effects ML regression" in out and "chibar2(01)" in out
    s.run_command("mixed y week x || id:, reml")
    assert abs(s.e["ll"] - (-536.3993448294261)) < 1e-6


def test_mixed_random_slope_unstructured():
    s, out = _run(LMM + "mixed y week x || id: week, covariance(unstructured)")
    assert abs(s.e["ll"] - (-503.69090055172876)) < 1e-6
    assert "cov(week,_cons)" in out
    s.run_command("mixed y week x || id: week")
    assert abs(s.e["ll"] - (-503.74239476706305)) < 1e-6


def test_mixed_nested_and_postestimation():
    s, out = _run("""clear
set seed 1002
set obs 20
gen school = _n
gen us = rnormal()*1.5
expand 4
bysort school: gen class = _n
gen uc = rnormal()*.8
expand 6
gen x = rnormal()
gen y = 5 + us + uc + 0.7*x + rnormal()
mixed y x || school: || class:
estat icc
predict yf, fitted
predict r, residuals""")
    assert abs(s.e["ll"] - (-753.654834940601)) < 1e-6
    assert "class|school" in out          # innermost level first, as Stata
    r = s.data.get("r").data
    assert abs(r.mean()) < 1e-6


def test_melogit_and_mepoisson_match_exact_integration():
    s, out = _run(GLMM + "melogit yb x || g:")
    b = _b(s)
    assert np.allclose(b[:2], [0.89783905, -0.55790145], atol=2e-5)
    assert abs(s.e["ll"] - (-355.94336106019165)) < 1e-5
    assert "var(_cons)|" in out
    s.run_command("mepoisson yc x || g:")
    b = _b(s)
    assert np.allclose(b[:2], [0.51818027, 0.17411692], atol=2e-5)
    assert abs(s.e["ll"] - (-881.8990902521404)) < 1e-4      # 7-point quadrature error


def test_melogit_icc_and_predict():
    s, out = _run(GLMM + "melogit yb x || g:\nestat icc\npredict p\npredict re, reffects")
    assert "Residual intraclass correlation" in out
    p = s.data.get("p").data
    assert np.all((p > 0) & (p < 1))


def test_reffects_one_variable_per_effect():
    s, out = _run(LMM + "mixed y week x || id: week\ncapture noisily predict re, reffects\n"
                  "predict re1 re0, reffects")
    assert "too few variables specified\nyou must specify 2 new variable(s)" in out
    ds = s.data
    # one value per group, in the order of e(revars): week, then _cons
    assert s.e["revars"] == "week _cons"
    ids = ds.get("id").data
    for v in ("re1", "re0"):
        vals = ds.get(v).data
        assert all(len(set(vals[ids == g].tolist())) == 1 for g in np.unique(ids))
    assert np.std(ds.get("re0").data) > np.std(ds.get("re1").data)


def test_icc_title_with_covariates():
    _s, out = _run(LMM + "mixed y week x || id:\nestat icc")
    assert "Residual intraclass correlation" in out
    _s, out = _run(LMM + "mixed y || id:\nestat icc")
    assert "\nIntraclass correlation" in out
