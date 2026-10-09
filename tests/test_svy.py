"""Fase 6: svyset e linearização de Taylor."""

import numpy as np
import pytest

DATA = """clear
set seed 600
set obs 200
gen estrato = ceil(_n/50)
gen upa = ceil(_n/10)
gen peso = 50 + mod(_n, 7)*10
gen sexo = mod(_n, 2)
gen imc = 22 + 2*sexo + round(runiform()*4, .01)
svyset upa [pweight=peso], strata(estrato)
"""


def _manual_mean_se(y, w, h, u):
    W = w.sum()
    m = (w * y).sum() / W
    z = w * (y - m) / W
    V = 0.0
    for hv in np.unique(h):
        sel = h == hv
        psus = np.unique(u[sel])
        zh = np.array([z[sel & (u == p)].sum() for p in psus])
        nh = len(psus)
        V += nh / (nh - 1) * ((zh - zh.mean()) ** 2).sum()
    return m, np.sqrt(V)


def test_svy_mean_matches_manual_linearization(run):
    run(DATA)
    out = run("svy: mean imc")
    ds = run.session.data
    y, w, h, u = (ds.get(n).data.astype(float) for n in ("imc", "peso", "estrato", "upa"))
    m, se = _manual_mean_se(y, w, h, u)
    b = run.session.e["b"].data[0, 0]
    V = run.session.e["V"].data[0, 0]
    assert b == pytest.approx(m, rel=1e-12)
    assert np.sqrt(V) == pytest.approx(se, rel=1e-10)
    assert "Linearized" in out and "Design df         =          16" in out


def test_subpop_keeps_design(run):
    run(DATA)
    run("svy, subpop(sexo): mean imc")
    assert run.eval("e(N_psu)") == 20 and run.eval("e(df_r)") == 16


def test_svy_regress_point_estimates_are_weighted_ls(run):
    run(DATA)
    run("quietly svy: regress imc sexo")
    ds = run.session.data
    y, w, x = (ds.get(n).data.astype(float) for n in ("imc", "peso", "sexo"))
    X = np.column_stack([x, np.ones_like(x)])
    b = np.linalg.solve((X * w[:, None]).T @ X, X.T @ (w * y))
    assert run.eval("_b[sexo]") == pytest.approx(b[0], rel=1e-10)


def test_requires_svyset(run):
    run("clear\nset obs 5\ngen y = _n")
    out = run("svy: mean y")
    assert "svyset" in out and run.rc == 119
