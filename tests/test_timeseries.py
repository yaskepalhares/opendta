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


# ---------------------------------------------------------------------------
# Phase 18: checks against closed forms and between commands (the Stata reference is compat/1904)
# ---------------------------------------------------------------------------

DATA18 = """clear
set seed 1904
set obs 200
gen t = _n
tsset t
gen e1 = rnormal()
gen e2 = rnormal()
gen x = rnormal()
gen y = 0
gen z = 0
replace y = .5*L.y + .2*L.z + e1 in 2/l
replace z = .1*L.y + .4*L.z + e2 in 2/l
forvalues i = 2/200 {
    quietly replace y = .5*y[`i'-1] + .2*z[`i'-1] + e1 in `i'
    quietly replace z = .1*y[`i'-1] + .4*z[`i'-1] + e2 in `i'
}
gen w = sum(rnormal())
gen c1 = w + .5*e1
gen c2 = .5*w + .5*e2
gen h = 1
gen g = 0
forvalues i = 2/200 {
    quietly replace h = .2 + .2*g[`i'-1]^2 + .5*h[`i'-1] in `i'
    quietly replace g = sqrt(h)*e1 in `i'
}
"""


def _run18(text):
    s = Session()
    lines = []
    s.output.add_listener(lambda t, style: lines.append(t))
    rc = s.run_text(DATA18 + text)
    assert rc == 0, "".join(lines)[-3000:]
    return s, "".join(lines)


def test_sspace_ar1_matches_arima():
    s, _ = _run18("arima y, ar(1) nolog")
    ll_arima = s.e["ll"]
    s.run_command("sspace (z1 L.z1 e.z1, state noconstant) (y z1@1, noerror), nolog")
    assert abs(s.e["ll"] - ll_arima) < 1e-3


def test_hp_filter_solves_normal_equations():
    s, _ = _run18("tsfilter hp yc = y, smooth(1600) trend(yt)")
    y = s.data.get("y").data
    yt = s.data.get("yt").data
    n = len(y)
    K = np.zeros((n - 2, n))
    for i in range(n - 2):
        K[i, i:i + 3] = [1, -2, 1]
    tau = np.linalg.solve(np.eye(n) + 1600 * K.T @ K, y)
    assert np.allclose(yt, tau, atol=1e-8)
    assert np.allclose(s.data.get("yc").data, y - tau, atol=1e-8)


def test_prais_twostep_rho_is_residual_autocorrelation():
    s, _ = _run18("gen r = 1 + .5*x + 0*y\nreplace r = r + .6*L.r + e2 in 2/l\nregress r x\npredict ehat, residuals")
    e = s.data.get("ehat").data
    e0, e1 = e[1:], e[:-1]
    rho = float(e0 @ e1 / (e1 @ e1))
    s.run_command("prais r x, corc rhotype(regress) twostep")
    assert abs(s.e["rho"] - rho) < 1e-6


def test_arch_loglik_recursion():
    s, _ = _run18("arch g, arch(1) garch(1) nolog")
    b = s.e["b"].data.ravel()
    mu, a, gm, w = b[0], b[1], b[2], b[3]
    y = s.data.get("g").data - mu
    h = np.zeros(len(y))
    s0 = float(np.mean(y ** 2))
    ll = 0.0
    for t in range(len(y)):
        e2 = y[t - 1] ** 2 if t else s0
        hp = h[t - 1] if t else s0
        h[t] = w + a * e2 + gm * hp
        ll += -0.5 * (np.log(2 * np.pi) + np.log(h[t]) + y[t] ** 2 / h[t])
    assert abs(ll - s.e["ll"]) < 1e-4


def test_vecrank_vec_and_fcast():
    s, out = _run18("vecrank c1 c2, lags(2)")
    assert "trace" in out
    s.run_command("vec c1 c2, rank(1) lags(2)")
    assert np.isfinite(s.e["ll"])
    s.run_command("fcast compute v_, step(3)")
    assert np.isfinite(s.data.get("v_c1").data[-1])


def test_svar_recursive_matches_cholesky():
    s, _ = _run18("matrix A1 = (1, 0 \\ ., 1)\nmatrix B1 = (., 0 \\ 0, .)\n"
                  "svar y z, lags(1/2) aeq(A1) beq(B1)")
    A = s.e["A"].data
    B = s.e["B"].data
    P = np.linalg.solve(A, B)
    s.run_command("var y z, lags(1/2)")
    Sig = s.e["Sigma"].data
    assert np.allclose(P, np.linalg.cholesky(Sig), atol=1e-4)


def test_ucm_local_level_and_mgarch_run():
    s, out = _run18("ucm y, model(llevel) nolog")
    assert "Unobserved-components model" in out and np.isfinite(s.e["ll"])
    s.run_command("mgarch ccc (g y = , noconstant), arch(1) garch(1) nolog")
    assert np.isfinite(s.e["ll"])
    s.run_command("dfactor (y z = , noconstant) (f = , ar(1)), nolog")
    assert np.isfinite(s.e["ll"])


def test_pergram_values_and_bartlett_b():
    # compat 1904: the generated periodogram is |sum x e^(-iwt)|^2 of the demeaned series, and
    # B = sqrt(q) max|F_k - (k-1)/(q-1)| with q = n/2 + 1
    s, out = _run("pergram y, gen(pg) nograph\nwntestb y, table")
    y = s.data.get("y").data[:150]
    pg = s.data.get("pg").data[:76]
    f = np.fft.rfft(y - y.mean())
    assert np.allclose(pg[1:], np.abs(f[1:]) ** 2, rtol=1e-9)
    F = np.cumsum(pg) / pg.sum()
    q = 76
    B = np.sqrt(q) * np.abs(F - np.arange(q) / (q - 1)).max()
    assert f"{B:.4f}" in out


def test_program_echo_does_not_count_comments():
    # compat 1904: the prompt after a comment line of a program body is "." and comments do not advance the count
    _s, out = _run("program define zz\n  gen a1 = 1\n  * note\n  gen a2 = 2\n  gen a3 = 3\nend\n")
    lines = [ln for ln in out.splitlines() if "gen a" in ln or "note" in ln]
    assert lines[1].lstrip().startswith("2.") and lines[1].strip().endswith("* note")
    assert lines[2].startswith(".") and lines[3].lstrip().startswith("3.")


# -- block 17 (phase 29) -------------------------------------------------------------------------

def test_prais_robust_and_cluster(run):
    run("clear\nset seed 19\nset obs 120\ngen t = _n\ntsset t\ngen x = round(invnormal(runiform()),.01)\n"
        "gen e = round(invnormal(runiform()),.01)\nreplace e = .6*e[_n-1] + e in 2/l\ngen y = 1 + 2*x + e\n"
        "gen g = ceil(_n/10)")
    run("quietly prais y x\nmatrix b0 = e(b)")
    out = run("prais y x, vce(robust) nolog")
    assert "Semirobust" in out and run.eval("mreldif(e(b), b0)") < 1e-12
    run("quietly prais y x, vce(cluster g)")
    assert run.eval("e(N_clust)") == 12 and run.eval("e(df_r)") == 11


def test_arch_constraints_and_from(run):
    run("clear\nset seed 20\nset obs 300\ngen t = _n\ntsset t\ngen e = round(invnormal(runiform()),.001)\n"
        "gen x = round(invnormal(runiform()),.01)\ngen h = 1\ngen y = 0\n"
        "replace h = .2 + .3*(e[_n-1]*sqrt(h[_n-1]))^2 + .4*h[_n-1] in 2/l\nreplace y = 1 + .5*x + e*sqrt(h)")
    run("quietly arch y x, arch(1) garch(1)\nmatrix b0 = e(b)")
    run("constraint define 1 [ARCH]L.arch = .3")
    out = run("arch y x, arch(1) garch(1) constraints(1) nolog")
    assert "( 1)  [ARCH]L.arch = .3" in out and "(constrained)" in out
    assert abs(run.eval("el(e(b), 1, 3)") - .3) < 1e-10
    run("quietly arch y x, arch(1) garch(1) from(b0)")
    assert run.rc == 0 and run.eval("mreldif(e(b), b0)") < 1e-3


def test_arfima_robust_svar_constraints(run):
    run("clear\nset seed 21\nset obs 200\ngen t = _n\ntsset t\ngen e1 = round(invnormal(runiform()),.01)\n"
        "gen e2 = round(invnormal(runiform()),.01)\ngen y1 = e1\ngen y2 = e2\n"
        "replace y1 = .5*y1[_n-1] + e1 in 2/l\nreplace y2 = .3*y1[_n-1] + .4*y2[_n-1] + .5*e1 + e2 in 2/l")
    run("quietly arfima y1, ar(1) vce(robust)")
    assert run.rc == 0 and run.eval('e(vcetype) == "Robust"') == 1
    run("matrix A = (1,0\\.,1)\nmatrix B = (.,0\\0,.)\nquietly svar y1 y2, lags(1) aeq(A) beq(B)")
    ll0 = run.eval("e(ll)")
    run("constraint define 1 [a_1_1]_cons = 1\nconstraint define 2 [a_1_2]_cons = 0\n"
        "constraint define 3 [a_2_2]_cons = 1\nconstraint define 4 [b_1_2]_cons = 0\n"
        "constraint define 5 [b_2_1]_cons = 0")
    run("quietly svar y1 y2, lags(1) aconstraints(1/3) bconstraints(4 5)")
    assert abs(run.eval("e(ll)") - ll0) < 1e-6 and run.eval("e(N_cns)") == 5


def test_state_space_predict_and_psdensity(run):
    run("clear\nset seed 24\nset obs 150\ngen t = _n\ntsset t\ngen e = round(invnormal(runiform()),.01)\n"
        "gen y = e\nreplace y = .6*y[_n-1] + e in 2/l\ngen y2 = round(y + invnormal(runiform()), .01)")
    run("constraint define 1 [y]u = 1")
    run("quietly sspace (u L.u, state noconstant) (y u, noerror noconstant), constraints(1)")
    run("predict yh\npredict r, residuals")
    assert run.rc == 0 and abs(run.eval("y[10] - yh[10] - r[10]")) < 1e-6
    run("quietly dfactor (y y2 = , noconstant) (f = , ar(1))\npredict rs, rstandard equation(y2)")
    assert run.rc == 0
    run("quietly ucm y, model(rwalk) cycle(1, frequency(.5))\npredict u1\npsdensity sd w")
    assert run.rc == 0 and run.session.data.has("sd")
