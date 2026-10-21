"""Phase 22: systems, choice models, exact methods, nonlinear models and treatment effects. The
Stata reference is compat 1908; here the checks are identities between commands and independent
computations with numpy."""

import numpy as np

DATA = """
clear
set seed 2208
set obs 300
gen double x1 = round(invnormal(runiform()), .01)
gen double x2 = round(invnormal(runiform()), .01)
gen double z1 = round(invnormal(runiform()), .01)
gen double e1 = round(invnormal(runiform()), .01)
gen double y1 = round(1 + 2*x1 - x2 + e1, .01)
gen double y2 = round(-1 + z1 + .5*e1 + invnormal(runiform()), .01)
gen double sdv = round(.5 + runiform(), .01)
"""


def _col(run, name):
    return run.session.data.get(name).data.astype(float)


def _b(run):
    return np.asarray(run.session.e["b"].data, dtype=float).ravel()


# ---------------------------------------------------------------------------
# linear systems
# ---------------------------------------------------------------------------

def test_sureg_with_equal_regressors_is_ols(run):
    run(DATA + "quietly regress y1 x1 x2\nmatrix b1 = e(b)\nquietly regress y2 x1 x2\nmatrix b2 = e(b)\n"
        "quietly sureg (y1 x1 x2) (y2 x1 x2)")
    b = _b(run)
    b1 = np.asarray(run.session.matrices["b1"].data).ravel()
    b2 = np.asarray(run.session.matrices["b2"].data).ravel()
    assert np.allclose(b, np.concatenate([b1, b2]), atol=1e-8)


def test_sureg_matches_gls_computation(run):
    run(DATA + "quietly sureg (y1 x1 x2) (y2 z1)")
    n = 300
    one = np.ones(n)
    y1, y2 = _col(run, "y1"), _col(run, "y2")
    X1 = np.column_stack([_col(run, "x1"), _col(run, "x2"), one])
    X2 = np.column_stack([_col(run, "z1"), one])
    e = np.column_stack([y1 - X1 @ np.linalg.lstsq(X1, y1, rcond=None)[0],
                         y2 - X2 @ np.linalg.lstsq(X2, y2, rcond=None)[0]])
    Si = np.linalg.inv(e.T @ e / n)
    A = np.block([[Si[0, 0] * X1.T @ X1, Si[0, 1] * X1.T @ X2], [Si[1, 0] * X2.T @ X1, Si[1, 1] * X2.T @ X2]])
    c = np.concatenate([Si[0, 0] * X1.T @ y1 + Si[0, 1] * X1.T @ y2, Si[1, 0] * X2.T @ y1 + Si[1, 1] * X2.T @ y2])
    assert np.allclose(_b(run), np.linalg.solve(A, c), atol=1e-8)


def test_sureg_output_layout(run):
    out = run(DATA + "sureg (y1 x1 x2) (y2 z1), corr")
    assert "Seemingly unrelated regression" in out
    assert 'Equation          Obs Parms        RMSE    "R-sq"       chi2        P' in out
    assert "Breusch-Pagan test of independence: chi2(1)" in out


def test_cnsreg_imposes_the_constraint(run):
    run(DATA + "constraint define 1 x1 = -2*x2\nquietly cnsreg y1 x1 x2, constraints(1)")
    b = _b(run)
    assert abs(b[0] + 2 * b[1]) < 1e-9
    y, X = _col(run, "y1"), np.column_stack([_col(run, "x1"), _col(run, "x2"), np.ones(300)])
    Z = np.column_stack([X[:, 0] - 2 * X[:, 1], np.ones(300)])      # x1*b1 + x2*b2 with b1 = -2 b2: b2*(x2 - 2 x1)
    Z = np.column_stack([X[:, 1] - 2 * X[:, 0], np.ones(300)])
    g = np.linalg.lstsq(Z, y, rcond=None)[0]
    assert np.allclose(b, [-2 * g[0], g[0], g[1]], atol=1e-8)
    assert run.session.e["df_r"] == 300 - 3 + 1


def test_vwls_matches_weighted_least_squares(run):
    run(DATA + "quietly vwls y1 x1 x2, sd(sdv)")
    b = _b(run)
    y, X, w = _col(run, "y1"), np.column_stack([_col(run, "x1"), _col(run, "x2"), np.ones(300)]), 1 / _col(run, "sdv") ** 2
    ref = np.linalg.solve((X * w[:, None]).T @ X, (X * w[:, None]).T @ y)
    assert np.allclose(b, ref, atol=1e-8)
    se = np.sqrt(np.diag(np.linalg.inv((X * w[:, None]).T @ X)))
    V = np.asarray(run.session.e["V"].data)
    assert np.allclose(np.sqrt(np.diag(V)), se, atol=1e-8)


def test_eivreg_with_full_reliability_is_ols(run):
    run(DATA + "quietly regress y1 x1 x2\nmatrix b0 = e(b)\nquietly eivreg y1 x1 x2")
    b0 = np.asarray(run.session.matrices["b0"].data).ravel()
    assert np.allclose(_b(run), b0, atol=1e-8)


def test_eivreg_corrects_attenuation(run):
    run(DATA + "quietly regress y1 x1 x2\nmatrix b0 = e(b)\nquietly eivreg y1 x1 x2, reliab(x1 .8)")
    b0 = np.asarray(run.session.matrices["b0"].data).ravel()
    assert abs(_b(run)[0]) > abs(b0[0])


def test_rreg_resists_outliers_and_saves_weights(run):
    run(DATA + "replace y1 = y1 + 40 in 1/6\nquietly regress y1 x1 x2\nmatrix b0 = e(b)\n"
        "quietly rreg y1 x1 x2, nolog genwt(wt)")
    b0 = np.asarray(run.session.matrices["b0"].data).ravel()
    assert abs(_b(run)[0] - 2) < abs(b0[0] - 2)
    wt = _col(run, "wt")
    assert wt.min() >= 0 and wt.max() <= 1 and wt[:6].max() < .5


def test_reg3_2sls_equals_ivregress_when_exactly_identified(run):
    run(DATA + "quietly reg3 (y1 y2 x1) (y2 y1 z1), 2sls\nmatrix b3 = e(b)\n"
        "quietly ivregress 2sls y1 x1 (y2 = z1)\nmatrix bi = e(b)")
    b3 = np.asarray(run.session.matrices["b3"].data).ravel()
    bi = np.asarray(run.session.matrices["bi"].data).ravel()
    # reg3 lists y2 x1 _cons of the first equation; ivregress lists y2 x1 _cons as well
    assert np.allclose(b3[:3], bi[:3], atol=1e-6)


# ---------------------------------------------------------------------------
# choice models
# ---------------------------------------------------------------------------

CHOICE = """
clear
set seed 5
set obs 1200
gen id = ceil(_n/4)
gen alt = mod(_n-1,4)+1
gen double x1 = invnormal(runiform())
gen double x2 = invnormal(runiform())
gen double u = -ln(-ln(runiform()))
gen double v = .8*x1 - .5*x2 + u
bysort id (v): gen byte ch = (_n==4)
bysort id (v): gen byte rk = _n
"""


def test_clogit_matches_direct_maximization(run):
    from scipy import optimize as opt
    run(CHOICE + "quietly clogit ch x1 x2, group(id)")
    b = _b(run)
    y, X = _col(run, "ch"), np.column_stack([_col(run, "x1"), _col(run, "x2")])

    def nll(t):
        eta = (X @ t).reshape(-1, 4)
        return -np.sum((y.reshape(-1, 4) * eta).sum(1) - np.log(np.exp(eta).sum(1)))

    ref = opt.minimize(nll, np.zeros(2), method="BFGS", options={"gtol": 1e-10}).x
    assert np.allclose(b, ref, atol=1e-5)
    assert abs(run.session.e["ll"] + nll(ref)) < 1e-8
    assert abs(run.session.e["ll_0"] + 300 * np.log(4)) < 1e-9


def test_clogit_drops_uninformative_groups(run):
    out = run(CHOICE + "replace ch = 0 if id == 1\nclogit ch x1 x2, group(id)")
    assert "note: 1 group (4 obs) dropped because of all positive or all negative outcomes." in out
    assert run.session.e["N"] == 1196


def test_rologit_without_ties_is_an_exploded_logit(run):
    from scipy import optimize as opt
    run(CHOICE + "quietly rologit rk x1 x2, group(id)")
    b = _b(run)
    r, X = _col(run, "rk").reshape(-1, 4), np.column_stack([_col(run, "x1"), _col(run, "x2")]).reshape(-1, 4, 2)

    def nll(t):
        tot = 0.0
        for i in range(r.shape[0]):
            order = np.argsort(-r[i])
            e = X[i] @ t
            for m in range(3):
                rest = e[order[m:]]
                tot -= e[order[m]] - np.log(np.exp(rest).sum())
        return tot

    ref = opt.minimize(nll, np.zeros(2), method="BFGS", options={"gtol": 1e-9}).x
    assert np.allclose(b, ref, atol=1e-4)


def test_rologit_tie_methods_agree_without_ties(run):
    run(CHOICE + "quietly rologit rk x1 x2, group(id)\nmatrix b1 = e(b)\n"
        "quietly rologit rk x1 x2, group(id) ties(efron)\nmatrix b2 = e(b)\n"
        "quietly rologit rk x1 x2, group(id) ties(breslow)\nmatrix b3 = e(b)")
    m = run.session.matrices
    for k in ("b2", "b3"):
        assert np.allclose(np.asarray(m["b1"].data), np.asarray(m[k].data), atol=1e-7)


def test_asclogit_equals_clogit_on_the_expanded_design(run):
    run(CHOICE + "gen double inc = round(invnormal(runiform()), .01)\nbysort id (alt): replace inc = inc[1]\n"
        "quietly asclogit ch x1, case(id) alternatives(alt) casevars(inc) basealternative(1)\n"
        "matrix ba = e(b)\n"
        "forvalues j = 2/4 {\n gen byte d`j' = alt == `j'\n gen double i`j' = inc*d`j'\n}\n"
        "quietly clogit ch x1 d2 d3 d4 i2 i3 i4, group(id)\nmatrix bc = e(b)")
    ba = np.asarray(run.session.matrices["ba"].data).ravel()
    bc = np.asarray(run.session.matrices["bc"].data).ravel()
    # asclogit: x1, then [2]: inc, _cons, [3]: ..., [4]: ...; clogit: x1, d2..d4, i2..i4
    assert np.isclose(ba[0], bc[0], atol=1e-6)
    for k, j in enumerate((2, 3, 4)):
        assert np.isclose(ba[1 + 2 * k], bc[4 + k], atol=1e-6)          # inc in alternative j
        assert np.isclose(ba[2 + 2 * k], bc[1 + k], atol=1e-6)          # constant of alternative j


def test_ghk_matches_the_normal_cdf():
    from scipy import stats
    from opendta.core import choice as C
    S = np.array([[2.0, 1.0, .5], [1.0, 2.0, 1.0], [.5, 1.0, 2.0]])
    L = np.linalg.cholesky(S)
    mu = np.array([[.3, -.2, .5], [1.0, 1.2, .1]])
    U = C.halton(4000, 3)
    p = C.ghk_upper(mu, L, U)
    ref = [stats.multivariate_normal.cdf(m, mean=np.zeros(3), cov=S) for m in mu]
    assert np.allclose(p, ref, atol=2e-3)


def test_mprobit_with_two_outcomes_is_close_to_probit(run):
    run("""
clear
set seed 77
set obs 600
gen double x = round(invnormal(runiform()), .01)
gen byte y = (.2 + .7*x + invnormal(runiform())) > 0
quietly probit y x
matrix bp = e(b)
quietly mprobit y x, baseoutcome(0) intpoints(400)
matrix bm = e(b)
""")
    bp = np.asarray(run.session.matrices["bp"].data).ravel()
    bm = np.asarray(run.session.matrices["bm"].data).ravel()
    # the differences of two iid N(0,1) errors have variance 2: coefficients are sqrt(2) times probit's
    assert np.allclose(bm[[0, 1]], np.sqrt(2) * bp[[0, 1]], atol=1e-3)


def test_slogit_fixes_phi_and_theta_of_the_base(run):
    out = run(CHOICE.replace("rk", "rk") + "gen byte y4 = rk\nslogit y4 x1 x2, baseoutcome(4)")
    assert "/phi1_1" in out and "(constrained)" in out and "(base outcome)" in out
    assert run.session.e["converged"] == 1
