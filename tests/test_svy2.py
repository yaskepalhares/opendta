"""Fase 6b: svy com réplicas, vários estágios, estat effects e mais modelos."""

import numpy as np
import pytest

DATA = """clear
set seed 602
set obs 200
gen estrato = ceil(_n/50)
gen upa = ceil(_n/10)
gen upa2 = mod(_n, 2)
gen peso = 50 + mod(_n, 7)*10
gen sexo = mod(_n, 2)
gen idade = 18 + floor(runiform()*60)
gen imc = 22 + .05*idade + 2*sexo + round(runiform()*4, .01)
gen nivel = 1 + (imc > 26) + (imc > 28)
gen ativo = runiform() < invlogit(-1 + .02*idade - .5*sexo)
"""


def _data(run):
    ds = run.session.data
    return {n: ds.get(n).data.astype(float) for n in ("estrato", "upa", "upa2", "peso", "imc", "idade", "nivel")}


def _wmean(y, w):
    return float(w @ y / w.sum())


def test_jackknife_mean_matches_manual(run):
    run(DATA)
    run("svyset upa [pweight=peso], strata(estrato)")
    run("svy, vce(jackknife): mean imc")
    d = _data(run)
    y, w, h, u = d["imc"], d["peso"], d["estrato"], d["upa"]
    theta = _wmean(y, w)
    V = 0.0
    for hh in np.unique(h):
        units = np.unique(u[h == hh])
        nh = len(units)
        reps = []
        for uu in units:
            ww = w.copy()
            ww[h == hh] *= nh / (nh - 1)
            ww[(h == hh) & (u == uu)] = 0
            reps.append(_wmean(y, ww))
        reps = np.array(reps).astype(np.float32).astype(float)   # réplicas em float, como no Stata
        V += (nh - 1) / nh * ((reps - reps.mean()) ** 2).sum()
    e = run.session.e
    assert e["b"].data[0, 0] == pytest.approx(theta, rel=1e-12)
    assert e["V"].data[0, 0] == pytest.approx(V, rel=1e-10)
    assert e["vce"] == "jackknife"


def test_brr_mean_with_hadamard(run):
    run(DATA)
    run("svyset upa2 [pweight=peso], strata(estrato) vce(brr)")
    out = run("capture noisily svy: mean imc")
    assert "option hadamard() required when brrweights are not svyset" in out
    run("matrix h4 = (1, 1, 1, 1 \\ 1, -1, 1, -1 \\ 1, 1, -1, -1 \\ 1, -1, -1, 1)")
    out = run("svy, hadamard(h4): mean imc")
    assert "BRR" in out
    d = _data(run)
    y, w, h, u = d["imc"], d["peso"], d["estrato"], d["upa2"]
    H = np.array([[1, 1, 1, 1], [1, -1, 1, -1], [1, 1, -1, -1], [1, -1, -1, 1]], float)
    reps = []
    for r in range(4):
        ww = w.copy()
        for k, hh in enumerate(np.unique(h)):
            first = (h == hh) & (u == 0)
            second = (h == hh) & (u == 1)
            a, b = (2, 0) if H[r, k] > 0 else (0, 2)
            ww[first] *= a
            ww[second] *= b
        reps.append(_wmean(y, ww))
    reps = np.array(reps).astype(np.float32).astype(float)   # réplicas em float, como no Stata
    V = ((reps - reps.mean()) ** 2).sum() / 4
    assert run.session.e["V"].data[0, 0] == pytest.approx(V, rel=1e-10)


def test_bootstrap_and_sdr_replicate_weights(run):
    run(DATA)
    run("forvalues r = 1/8 {\n gen bw`r' = peso*(1 + .1*sin(_n*`r'))\n}")
    run("svyset [pweight=peso], bsrweight(bw*) vce(bootstrap) bsn(2)")
    run("svy: mean imc")
    ds = run.session.data
    y, w = ds.get("imc").data.astype(float), ds.get("peso").data.astype(float)
    reps = np.array([_wmean(y, ds.get(f"bw{r}").data.astype(float)) for r in range(1, 9)]).astype(np.float32).astype(float)
    V = 2 / 8 * ((reps - reps.mean()) ** 2).sum()
    assert run.session.e["V"].data[0, 0] == pytest.approx(V, rel=1e-10)
    run("svyset [pweight=peso], sdrweight(bw*) vce(sdr) mse")
    run("svy: mean imc")
    V = 4 / 8 * ((reps - _wmean(y, w)) ** 2).sum()
    assert run.session.e["V"].data[0, 0] == pytest.approx(V, rel=1e-10)


def test_multistage_display_and_second_stage(run):
    run(DATA)
    run("gen npop = 20\ngen n2 = 20")
    out = run("svyset upa [pweight=peso], strata(estrato) fpc(npop) || _n, fpc(n2)")
    assert "SU 2: <observations>" in out and "FPC 2: n2" in out
    run("svy: mean imc")
    v2 = run.session.e["V"].data[0, 0]
    run("svyset upa [pweight=peso], strata(estrato) fpc(npop)")
    run("svy: mean imc")
    v1 = run.session.e["V"].data[0, 0]
    # f1 = 5/20 > 0 e o 2º estágio com reposição parcial (10 de 20) soma variância
    assert v2 > v1


def test_estat_effects_srs_formula(run):
    run(DATA)
    run("svyset upa [pweight=peso], strata(estrato)")
    run("svy: mean imc")
    out = run("estat effects")
    assert "DEFF" in out and "DEFT" in out
    d = _data(run)
    y, w = d["imc"], d["peso"]
    M = w.sum()
    m = len(y)
    ybar = _wmean(y, w)
    Vsrs = M / (m - 1) * (w * ((y - ybar) / M) ** 2).sum()
    V = run.session.e["V"].data[0, 0]
    assert run.session.r["deff"].data[0, 0] == pytest.approx(V / Vsrs, rel=1e-10)


def test_svy_generic_models_point_estimates_match_weighted_fit(run):
    run(DATA)
    run("svyset upa [pweight=peso], strata(estrato)")
    run("svy: ologit nivel idade")
    b_svy = run.session.e["b"].data.copy()
    run("quietly ologit nivel idade [pweight=peso]")
    b_w = run.session.e["b"].data
    assert np.allclose(b_svy, b_w, rtol=1e-10)
    run("svy: ologit nivel idade")
    assert run.session.e["prefix"] == "svy"
    assert run.session.e["df_r"] == 16
    out = run("svy: mlogit nivel idade, rrr")
    assert "RRR" in out
