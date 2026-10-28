"""Phase 19: panel data (quadrature models, linear panel estimators, frontier, unit roots, dynamic
GMM and description). The Stata reference is compat 1905; here the checks are identities between
commands and independent computations with numpy and scipy."""

import math

import numpy as np
from scipy import integrate, stats as st

PANEL = """
clear
set seed 1905
set obs 240
gen id = ceil(_n/8)
gen t = _n - (id-1)*8
xtset id t
gen u = rnormal()
bysort id: replace u = u[1]
gen x1 = rnormal()
gen z = rnormal()
gen x2 = .5*z + rnormal() + .5*u
gen y = 1 + x1 + .5*x2 + u + rnormal()
gen b = runiform() < invlogit(-.2 + .8*x1 + .7*u)
gen x3 = x2
"""


def _b(run, name="b"):
    return np.asarray(run.session.e[name].data, dtype=float).ravel()


def test_xtgls_iid_is_ols(run):
    run(PANEL + "quietly xtgls y x1 x2, panels(iid)")
    b1 = _b(run)
    run("quietly regress y x1 x2")
    assert np.allclose(b1, _b(run), atol=1e-9)


def test_xtivreg_with_own_instruments_matches_xtreg(run):
    run(PANEL + "quietly xtivreg y x1 (x2 = x3), fe")
    b1 = _b(run)
    run("quietly xtreg y x1 x2, fe")
    # xtivreg lists the instrumented variable first
    assert np.allclose(b1[[1, 0, 2]], _b(run), atol=1e-8)
    run("quietly xtivreg y x1 (x2 = x3), be")
    b2 = _b(run)
    run("quietly xtreg y x1 x2, be")
    assert np.allclose(b2[[1, 0, 2]], _b(run), atol=1e-8)
    run("quietly xtivreg y x1 (x2 = x3)")
    b3 = _b(run)
    run("quietly xtreg y x1 x2, re")
    assert np.allclose(b3[[1, 0, 2]], _b(run), atol=1e-6)


def test_xtivreg_fd_is_first_difference_2sls(run):
    run(PANEL + "quietly xtivreg y x1 (x2 = z), fd")
    b = _b(run)
    d = lambda v: np.where(np.r_[False, np.diff(np.arange(len(v)) // 8) == 0], np.r_[np.nan, np.diff(v)], np.nan)
    s = run.session
    y, x1, x2, z = (s.data.get(n).data.astype(float) for n in ("y", "x1", "x2", "z"))
    dy, dx1, dx2, dz = d(y), d(x1), d(x2), d(z)
    m = ~np.isnan(dy)
    one = np.ones(m.sum())
    X = np.column_stack([dx2[m], dx1[m], one])      # the differenced equation has a constant
    Z = np.column_stack([dz[m], dx1[m], one])
    Pz = Z @ np.linalg.solve(Z.T @ Z, Z.T)
    ref = np.linalg.solve(X.T @ Pz @ X, X.T @ Pz @ dy[m])
    assert np.allclose(b, ref, atol=1e-8)


def test_xtpcse_independent_has_ols_coefficients(run):
    run(PANEL + "quietly xtpcse y x1 x2")
    b1 = _b(run)
    run("quietly regress y x1 x2")
    assert np.allclose(b1, _b(run), atol=1e-9)


def test_xtregar_fe_rhof_zero_is_within(run):
    run(PANEL + "quietly xtregar y x1 x2, fe rhof(0)")
    b1 = _b(run)[:2]
    # the first observation of every panel is lost in the transformation
    run("quietly xtreg y x1 x2 if t > 1, fe")
    assert np.allclose(b1, _b(run)[:2], atol=1e-6)


def test_xtrc_runs_and_stores(run):
    run(PANEL + "xtrc y x1 x2")
    assert run.session.e["N_g"] == 30
    assert "Test of parameter constancy" in run.cap.text


def test_xtprobit_quadrature_is_stable(run):
    out = run(PANEL + "xtprobit b x1, re intpoints(12)")
    ll12 = run.session.e["ll"]
    run("quietly xtprobit b x1, re intpoints(40)")
    assert abs(ll12 - run.session.e["ll"]) < 1e-3
    assert "Random-effects probit regression" in out


def test_quadchk_restores_estimates(run):
    run(PANEL + "quietly xtprobit b x1, re")
    b0 = _b(run).copy()
    run("quadchk")
    assert np.allclose(_b(run), b0)


def test_frontier_group_loglik_matches_integration():
    from opendta.core import panelfrontier as pf
    from opendta.core.panelquad import Groups
    rng = np.random.default_rng(1)
    inv = np.repeat(np.arange(3), 4)
    e = rng.normal(size=12)
    d = np.exp(-0.2 * np.tile(np.arange(4)[::-1], 3))
    mu, su2, sv2 = 0.3, 0.8, 0.5
    got = pf.group_loglik(e, d, Groups(inv), mu, su2, sv2)
    for i in range(3):
        m = inv == i

        def f(u):
            return (np.prod(st.norm.pdf(e[m] + d[m] * u, 0, math.sqrt(sv2)))
                    * st.norm.pdf(u, mu, math.sqrt(su2)) / st.norm.cdf(mu / math.sqrt(su2)))
        assert abs(math.log(integrate.quad(f, 0, np.inf)[0]) - got[i]) < 1e-8


def test_xtfrontier_runs(run):
    run(PANEL + "gen ly = 1 + .5*x1 + rnormal()*.3 - abs(u)*.3\nxtfrontier ly x1, ti")
    assert 0 < run.session.e["gamma"] < 1


def test_fisher_statistics_follow_the_panel_pvalues():
    from opendta.core import panelroot as pr
    p = np.array([.2, .5, .01, .7])
    r = pr.fisher_combine(p)
    assert abs(r["P"] + 2 * np.log(p).sum()) < 1e-12
    assert abs(r["p_P"] - st.chi2.sf(r["P"], 8)) < 1e-12
    assert abs(r["Z"] - st.norm.ppf(p).sum() / 2) < 1e-12


def test_hadri_closed_form():
    from opendta.core import panelroot as pr
    rng = np.random.default_rng(3)
    Y = rng.normal(size=(12, 20))
    r = pr.hadri(Y)
    E = Y - Y.mean(axis=1, keepdims=True)
    S = np.cumsum(E, axis=1)
    lm = (S ** 2).sum() / (12 * 20 ** 2 * (E ** 2).sum() / (12 * (20 - 1)))
    assert abs(r["lm"] - lm) < 1e-12
    assert abs(r["z"] - math.sqrt(12) * (lm - 1 / 6) / math.sqrt(1 / 45)) < 1e-12


def test_unit_root_tests_separate_random_walk_and_stationary(run):
    run("""
clear
set seed 5
set obs 600
gen id = ceil(_n/30)
gen t = _n - (id-1)*30
xtset id t
gen e = rnormal()
gen rw = e
bysort id (t): replace rw = rw[_n-1] + e if _n>1
gen st = e
bysort id (t): replace st = .5*st[_n-1] + e if _n>1
""")
    run("quietly xtunitroot breitung rw")
    p_rw = run.session.r["p"]
    run("quietly xtunitroot breitung st")
    assert p_rw > 0.05 > run.session.r["p"]
    run("quietly xtunitroot hadri rw")
    assert run.session.r["p"] < 0.01
    assert "must specify one of dfuller or pperron" in run("capture noisily xtunitroot fisher rw, lags(1)")
    run("quietly xtunitroot fisher rw, dfuller lags(1)")
    assert run.session.r["p_Z"] > 0.05 > 0 and "P" in run.session.r


def test_xtabond_matches_independent_gmm(run):
    run("""
clear
set seed 3
set obs 400
gen id = ceil(_n/5)
gen t = _n - (id-1)*5
xtset id t
gen u = rnormal()
bysort id: replace u = u[1]
gen y = u + rnormal()
bysort id (t): replace y = .5*y[_n-1] + u + rnormal() if _n>1
quietly xtabond y, noconstant
""")
    b = _b(run)
    s = run.session
    Y = s.data.get("y").data.astype(float).reshape(80, 5)
    Xs, ys, Zs, Hs = [], [], [], []
    for i in range(80):
        rows_x, rows_y, rows_z, ts = [], [], [], []
        for t in range(2, 5):                       # periods 3..5 (0-based 2..4): Dy_t and Dy_{t-1}
            rows_x.append([Y[i, t - 1] - Y[i, t - 2]])
            rows_y.append(Y[i, t] - Y[i, t - 1])
            z = np.zeros(1 + 2 + 3)                 # instruments y_{1..t-2}: t=3 -> 1, t=4 -> 2, t=5 -> 3
            n0 = {2: 0, 3: 1, 4: 3}[t]
            z[n0:n0 + t - 1] = Y[i, :t - 1]
            rows_z.append(z)
            ts.append(t)
        Xs.append(np.array(rows_x))
        ys.append(np.array(rows_y))
        Zs.append(np.array(rows_z))
        H = 2 * np.eye(3) - np.eye(3, k=1) - np.eye(3, k=-1)
        Hs.append(H)
    from opendta.core import panelgmm as pg
    fit = pg.gmm_fit(Xs, ys, Zs, Hs)
    assert abs(fit["b"][0] - b[0]) < 1e-9


def test_estat_abond_and_sargan_run(run):
    run(PANEL + "quietly xtabond y x1, twostep vce(robust)\nestat abond")
    assert "Arellano-Bond test for zero autocorrelation" in run.cap.text
    run("quietly xtabond y x1\nestat sargan")
    assert "Sargan test of overidentifying restrictions" in run.cap.text


def test_xtdpdsys_runs(run):
    out = run(PANEL + "xtdpdsys y x1")
    assert "System dynamic panel-data estimation" in out
    assert "Instruments for level equation" in out


def test_xttab_between_and_within(run):
    out = run("""
clear
set obs 12
gen id = ceil(_n/4)
gen t = mod(_n-1, 4) + 1
xtset id t
gen g = (id == 1 & t <= 2) | id == 2
xttab g
""")
    assert "(n = 3)" in out
    # g = 1 in two panels (id 1 and id 2): between percent 66.67; within: (50 + 100)/2
    assert "66.67" in out and "75.00" in out


def test_xtdata_fe_and_be(run):
    run(PANEL + "xtdata x1 y, fe clear")
    s = run.session
    x = s.data.get("x1").data.astype(float).reshape(30, 8)
    assert np.allclose(x.mean(axis=1), x.mean())
    run(PANEL)
    means = run.session.data.get("x1").data.astype(float).reshape(30, 8).mean(axis=1)
    run(PANEL + "xtdata x1, be clear")
    x = run.session.data.get("x1").data.astype(float)
    assert len(x) == 30                                 # one observation per panel
    assert np.allclose(x, means)


# -- block 16 (phase 28) -------------------------------------------------------------------------

def _iv_panel(run):
    run("clear\nset seed 160\nset obs 300\ngen id = ceil(_n/10)\ngen t = mod(_n-1,10)+1\nxtset id t\n"
        "gen z = round(invnormal(runiform()),.01)\ngen x2 = round(invnormal(runiform()),.01)\n"
        "gen u = round(invnormal(runiform()),.01)\ngen x1 = z + .5*u + round(invnormal(runiform()),.01)\n"
        "gen y = 1 + x1 + x2 + u + id/10")


def test_xtivreg_regress_and_first(run):
    _iv_panel(run)
    run("quietly xtreg y x1 x2, re\nmatrix b0 = e(b)")
    run("quietly xtivreg y x2 (x1 = z), re regress")
    assert run.eval("mreldif(e(b), b0)") < 1e-8
    out = run("xtivreg y x2 (x1 = z), fe first")
    assert "First-stage within regression" in out and "Instrumented:   x1" in out
    assert run.eval('e(model) == "fe"') == 1


def test_xtdata_re_reproduces_xtreg_re(run):
    run("clear\nset seed 160\nset obs 300\ngen id = ceil(_n/10)\ngen t = mod(_n-1,10)+1\nxtset id t\n"
        "gen x1 = round(invnormal(runiform()),.01)\ngen y = 1 + x1 + round(invnormal(runiform()),.01) + mod(id,7)/3")
    run("quietly xtreg y x1, re")
    b1, b0 = run.eval("_b[x1]"), run.eval("_b[_cons]")
    r = run.eval("e(sigma_u)/e(sigma_e)")
    run(f"xtdata y x1, re ratio({r!r}) clear\nquietly regress y x1 constant, noconstant")
    assert abs(run.eval("_b[x1]") - b1) < 1e-6 and abs(run.eval("_b[constant]") - b0) < 1e-6


def test_xtrc_betas(run):
    run("clear\nset seed 17\nset obs 60\ngen id = ceil(_n/12)\ngen t = mod(_n-1,12)+1\nxtset id t\n"
        "gen x = round(invnormal(runiform()),.01)\ngen y = 1 + (1 + id/5)*x + round(invnormal(runiform()),.01)")
    out = run("xtrc y x, betas")
    assert "Group-specific coefficients" in out and "Group 5" in out
    assert run.eval("rowsof(e(beta_ps))") == 5


def test_xtunitroot_fisher_and_hadri_options(run):
    run("clear\nset seed 16\nset obs 400\ngen id = ceil(_n/20)\ngen t = mod(_n-1,20)+1\nxtset id t\n"
        "gen e = round(invnormal(runiform()), .01)\ngen y = e\nbysort id (t): replace y = .5*y[_n-1] + e if _n > 1")
    for opts in ("dfuller lags(1) drift", "dfuller lags(1) demean", "pperron lags(2)"):
        run(f"quietly xtunitroot fisher y, {opts}")
        assert run.rc == 0 and run.eval("r(p_P)") < .05
    run("quietly xtunitroot hadri y")
    z0 = run.eval("r(z)")
    run("quietly xtunitroot hadri y, kernel(bartlett 3)")
    assert run.rc == 0 and run.eval("r(z)") != z0


def test_xtdpd_artests_and_hascons(run):
    run("clear\nset seed 5\nset obs 300\ngen id = ceil(_n/10)\ngen t = mod(_n-1,10)+1\nxtset id t\n"
        "gen x = round(invnormal(runiform()),.01)\ngen y = round(invnormal(runiform()),.01)\n"
        "bysort id (t): replace y = .4*y[_n-1] + x + y if _n > 1\ngen one = 1")
    run("quietly xtabond y x, artests(3) vce(robust)")
    assert run.rc == 0 and run.eval("e(artests)") == 3
    out = run("estat abond")
    assert "|   3  |" in out
    run("quietly xtdpd y L.y x one, dgmmiv(y) lgmmiv(y) div(x) liv(one) hascons")
    # the constant regressor is dropped as collinear and the model has _cons (compat 1914)
    assert run.rc == 0 and run.session.e["b"].colnames == ["L.y", "x", "_cons"]
