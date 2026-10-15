"""Time series: arima, var, dfuller, corrgram and newey, checked against
statsmodels 0.15 (ARIMA state-space ML, VAR, adfuller, acf/pacf,
acorr_ljungbox, OLS with HAC covariance)."""

import numpy as np

from opendta.session import Session

DATA = """clear
set seed 1101
set obs 150
gen t = _n
tsset t
gen e = rnormal()
gen u = e
replace u = 0.6*u[_n-1] + e + 0.3*e[_n-1] in 2/l
gen x = rnormal()
gen y = 2 + 0.5*x + u
gen w = sum(rnormal() + .1)
gen z = 0.4*L.y + rnormal() in 2/l
"""


def _run(text):
    s = Session()
    lines = []
    s.output.add_listener(lambda t, style: lines.append(t))
    rc = s.run_text(DATA + text)
    assert rc == 0, "".join(lines)[-3000:]
    return s, "".join(lines)


def test_arima_exact_ml():
    s, out = _run("arima y x, ar(1) ma(1)\npredict r, residuals")
    b = s.e["b"].data.ravel()
    assert abs(s.e["ll"] - (-202.65585934555065)) < 1e-6
    assert np.allclose(b, [0.46993946, 2.08973926, 0.35086521, 0.32230723, 0.9329323811308793], atol=1e-4)
    assert "OPG" in out
    r = s.data.get("r").data
    assert abs(np.mean(r)) < 0.2
    s.run_command("arima w, arima(1,1,1)")
    assert abs(s.e["ll"] - (-193.4257728542805)) < 1e-5


def test_var():
    s, out = _run("var y z, lags(1/2)")
    assert abs(s.e["ll"] - (-415.18333998332474)) < 1e-8
    assert abs(s.e["fpe"] - 1.1152802499920909) < 1e-8
    assert abs(s.e["detsig_ml"] - 0.9733600658258537) < 1e-10
    s.run_command("vargranger")
    s.run_command("varsoc y z")


def test_dfuller_and_mackinnon():
    s, out = _run("dfuller w")
    assert abs(s.r["Zt"] - (-1.4940455677954572)) < 1e-9
    assert abs(s.r["p"] - 0.5364548279914909) < 1e-6
    s.run_command("dfuller w, trend lags(2)")
    assert abs(s.r["Zt"] - (-1.8964210254268725)) < 1e-9
    assert abs(s.r["p"] - 0.6565217400546705) < 1e-6
    s.run_command("dfuller y, noconstant")
    # no MacKinnon p-value without a constant (as Stata: neither displayed nor stored)
    assert "p" not in s.r


def test_dfuller_interpolated_critical_values():
    from opendta.estimation.tsmodels import _fuller_crit
    assert np.allclose(_fuller_crit("c", 79), [-3.5394, -2.9068, -2.5884], atol=1e-4)


def test_corrgram_and_newey():
    s, out = _run("corrgram y, lags(8)")
    assert np.allclose([s.r["ac1"], s.r["ac2"], s.r["ac3"]], [0.45739324, 0.18552095, 0.13820912], atol=1e-7)
    assert np.allclose([s.r["pac1"], s.r["pac2"], s.r["pac3"]], [0.46002521, -0.02807392, 0.07998498], atol=1e-7)
    assert abs(s.r["q8"] - 46.4545313) < 1e-5
    s.run_command("newey y x, lag(3)")
    se = np.sqrt(np.diag(s.e["V"].data))
    assert np.allclose(se, [0.10461293, 0.13942029], rtol=1e-6)
