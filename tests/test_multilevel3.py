"""Block 18: the me commands with nested levels, weights and vce(robust); mestreg with nested
levels, weights and the median, mean and hazard predictions."""

import math

import numpy as np
from opendta.session import Session

DATA = """clear
set seed 1916
set obs 25
gen s = _n
gen us = rnormal()*.8
expand 4
bysort s: gen c = _n
gen uc = rnormal()*.5
expand 6
gen x = rnormal()
gen yb = runiform() < invlogit(-.2 + .7*x + us + uc)
gen yg = 1 + .5*x + us + uc + rnormal()
gen fw = 1 + mod(_n, 3)
"""


def _run(text, s=None):
    s = s or Session()
    lines = []
    s.output.add_listener(lambda t, style: lines.append(t))
    rc = s.run_text(text)
    assert rc == 0, "".join(lines)[-3000:]
    return s, "".join(lines)


def test_two_levels_gaussian_matches_mixed():
    s, _ = _run(DATA + "quietly mixed yg x || s: || c:\n")
    ll_mixed, b_mixed = s.e["ll"], s.e["b"].data.ravel()
    s, out = _run("meglm yg x || s: || c:, intpoints(12)\n", s)
    assert abs(s.e["ll"] - ll_mixed) < 1e-4
    assert np.allclose(s.e["b"].data.ravel()[:2], b_mixed[:2], atol=1e-4)
    assert "s>c" in out and "Group Variable" in out
    assert s.e["k_r"] == 2 and s.e["ivars"] == "s c"


def test_two_levels_logit_reffects_and_icc():
    s, out = _run(DATA + "melogit yb x || s: || c:\nestat icc\npredict r1 r2, reffects\n")
    assert "c|s" in out
    r1, r2 = s.data.get("r1").data, s.data.get("r2").data
    # the effect of s is constant within s, the effect of c within s and c
    sv, cv = s.data.get("s").data, s.data.get("c").data
    for g in np.unique(sv)[:5]:
        assert np.ptp(r1[sv == g]) < 1e-12
        assert np.ptp(r2[(sv == g) & (cv == 1)]) < 1e-12
    assert np.ptp(r2[sv == 1]) > 0
    # more quadrature points change the likelihood little
    ll7 = s.e["ll"]
    s, _ = _run("quietly melogit yb x || s: || c:, intpoints(11)\n", s)
    assert abs(s.e["ll"] - ll7) < 2e-2


def test_fweights_equal_expanded_data():
    s, _ = _run(DATA + "quietly melogit yb x [fw=fw] || s:\n")
    ll_w, b_w = s.e["ll"], s.e["b"].data.ravel()
    assert s.e["N"] == float(np.sum(s.data.get("fw").data))
    s, _ = _run("expand fw\nquietly melogit yb x || s:\n", s)
    assert abs(s.e["ll"] - ll_w) < 1e-4
    assert np.allclose(s.e["b"].data.ravel(), b_w, atol=1e-4)


def test_level_fweight_equals_duplicated_groups():
    # fweight(w) of the highest level counts each group w times
    s, _ = _run(DATA + "gen ws = 1 + mod(s, 2)\nquietly meglm yg x || s:, fweight(ws)\n")
    ll_w = s.e["ll"]
    s2, _ = _run(DATA + "gen ws = 1 + mod(s, 2)\nexpand ws, generate(dup)\nreplace s = s + 100*dup\n"
                 "quietly meglm yg x || s:\n")
    assert abs(s2.e["ll"] - ll_w) < 1e-4


def test_robust_sandwich():
    s, out = _run(DATA + "melogit yb x || s:, vce(robust)\n")
    assert "Robust" in out and "adjusted for 25 clusters in s" in out
    assert "LR test" not in out and "Log pseudolikelihood" in out
    assert s.e["vce"] == "robust" and s.e["N_clust"] == 25
    V = s.e["V"].data
    s3, _ = _run(DATA + "gen cl = ceil(s/5)\nquietly melogit yb x || s:, vce(cluster cl)\n")
    assert s3.e["N_clust"] == 5 and s3.e["clustvar"] == "cl"
    assert not np.allclose(s3.e["V"].data, V)
    # clusters must contain the highest-level groups
    assert Session().run_text(DATA + "gen cl = mod(_n, 3)\nmelogit yb x || s:, vce(cluster cl)\n") == 198


def test_pweights_imply_robust_and_level_weight_constant():
    s, out = _run(DATA + "melogit yb x [pw=fw] || s:\n")
    assert s.e["vce"] == "robust" and "Robust" in out
    s = Session()
    rc = s.run_text(DATA + "melogit yb x || s:, pweight(fw)\n")
    assert rc == 198


def test_mestreg_levels_and_predictions():
    st = DATA + """gen double tw = (-ln(runiform())/exp(-1 + .5*x + us + uc))^(1/1.4)
gen dw = tw < 2
replace tw = min(tw, 2)
stset tw, failure(dw)
"""
    s, out = _run(st + "mestreg x || s: || c:, distribution(weibull)\npredict md, median\npredict mn, mean\n"
                  "predict h, hazard\npredict eta, eta\n")
    assert "s>c" in out
    p = math.exp(float(s.e["b"].data.ravel()[2]))
    eta = s.data.get("eta").data
    lam = np.exp(eta)
    assert np.allclose(s.data.get("md").data, (math.log(2) / lam) ** (1 / p))
    assert np.allclose(s.data.get("mn").data, math.gamma(1 + 1 / p) * lam ** (-1 / p))
    t = s.data.get("_t").data
    assert np.allclose(s.data.get("h").data, p * t ** (p - 1) * lam)
    s, _ = _run(st + "quietly mestreg x || s:, distribution(lognormal)\npredict md, median\npredict mn, mean\n"
                "predict eta, eta\n")
    sig = math.exp(float(s.e["b"].data.ravel()[2]))
    eta = s.data.get("eta").data
    assert np.allclose(s.data.get("md").data, np.exp(eta))
    assert np.allclose(s.data.get("mn").data, np.exp(eta + sig * sig / 2))
    s, out = _run(st.replace("stset tw,", "stset tw [pw=fw],") + "mestreg x || s:, distribution(exponential)\n")
    assert "Robust" in out and s.e["vce"] == "robust"
