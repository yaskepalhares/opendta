"""irf after var: simple and orthogonalized impulse responses, cumulative responses and the
variance decomposition, against statsmodels' VAR results."""

import numpy as np
import pytest


@pytest.fixture
def varfit(run, tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    rng = np.random.default_rng(3)
    T = 200
    y = np.zeros((T, 2))
    for t in range(1, T):
        e = rng.normal(size=2)
        y[t] = [.5 * y[t - 1, 0] + .2 * y[t - 1, 1] + e[0], .3 * y[t - 1, 0] + .4 * y[t - 1, 1] + .5 * e[0] + e[1]]
    run("clear\ninput t double y1 double y2\n"
        + "\n".join(f"{i + 1} {float(y[i, 0])!r} {float(y[i, 1])!r}" for i in range(T)) + "\nend\ntsset t")
    run("quietly var y1 y2, lags(1 2)")
    run.y = y
    return run


def test_irf_create_and_values(varfit):
    run = varfit
    out = run("irf create a, set(res) step(5)")
    assert "res.irf" in out
    from opendta.estimation.irf import _load
    from pathlib import Path
    d = _load(Path("res.irf"))["a"]
    VAR = pytest.importorskip("statsmodels.tsa.api").VAR
    r = VAR(run.y).fit(2, trend="c")
    ir = r.irf(5)
    assert np.allclose(d["stats"]["irf"], ir.irfs, atol=1e-10)
    Sml = r.sigma_u * r.df_resid / r.nobs          # var keeps the ML Sigma, like Stata
    assert np.allclose(d["stats"]["oirf"], ir.irfs @ np.linalg.cholesky(Sml), atol=1e-10)
    assert np.allclose(d["stats"]["fevd"][1:6].transpose(1, 0, 2), r.fevd(6).decomp[:, :5, :], atol=1e-10)
    # standard errors: same formula, with the ML instead of the df-adjusted covariance
    ratio = d["se"]["irf"][1:] / ir.stderr(orth=False)[1:]
    assert np.allclose(ratio, np.sqrt(r.df_resid / r.nobs), rtol=1e-4)


def test_irf_table_graph_describe(varfit):
    run = varfit
    run("irf set res\nirf create a, step(4)\nirf create b, order(y2 y1) step(4)")
    out = run("irf table oirf, impulse(y1) response(y2)")
    assert "Results from a b" in out and "Lower" in out and "irfname = b" in out
    out = run("irf table coirf, irf(a) noci stderror")
    assert "S.E." in out and "Lower" not in out
    out = run("irf describe")
    assert "| var      y2 y1" in out and "(*) order is relevant" in out
    run("irf graph oirf, impulse(y1) response(y2) nodraw")
    assert run.rc == 0
    run("irf rename b c\nirf drop c")
    run("capture irf table irf, irf(c)")
    assert run.eval("_rc") == 111


def test_irf_errors(run):
    run("clear\nset obs 5\ngen y = _n\nregress y")
    run("capture irf create x, set(z)")
    assert run.eval("_rc") == 321
