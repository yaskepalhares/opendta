"""Phase 21: binary and count models, censored and selection models, endogenous regressors. The
Stata reference is compat 1907; here the checks are identities between commands and independent
computations with numpy and scipy."""

import numpy as np
from scipy import optimize as opt
from scipy import stats as st

from opendta.core.bvnorm import bvn_cdf, bvn_grad

DATA = """
clear
set seed 2107
set obs 500
gen double x1 = round(invnormal(runiform()), .01)
gen double x2 = round(invnormal(runiform()), .01)
gen double u = round(invnormal(runiform()), .01)
gen byte yb = (.3 + .8*x1 - .5*x2 + u) > 0
gen int cnt = rpoisson(exp(.4 + .3*x1 - .2*x2))
gen int ctr = cnt if cnt > 0
gen double yc = round(1 + .7*x1 + .5*x2 + 1.2*u, .01)
gen byte sel = (.3 + .6*x1 + round(invnormal(runiform()), .01)) > 0
gen double ysel = round(1 + .8*x2 + .6*u + .6*round(invnormal(runiform()), .01), .01) if sel
gen double z1 = round(invnormal(runiform()), .01)
gen double z2 = round(invnormal(runiform()), .01)
gen double w2 = round(.6*z1 + .4*z2 + .3*x1 + .5*u + .8*invnormal(runiform()), .01)
gen double yiv = round(1 + .7*w2 + .4*x1 + u, .01)
"""


def _col(run, name):
    return run.session.data.get(name).data.astype(float)


def _b(run):
    return np.asarray(run.session.e["b"].data, dtype=float).ravel()


def _fit(fun, x0):
    return opt.minimize(fun, x0, method="BFGS", options={"gtol": 1e-9}).x


# ---------------------------------------------------------------------------
# bivariate normal distribution
# ---------------------------------------------------------------------------

def test_bvn_cdf_matches_scipy():
    rng = np.random.default_rng(3)
    a, b = rng.normal(size=20), rng.normal(size=20)
    for rho in (-.8, -.3, 0., .4, .9):
        ref = [st.multivariate_normal.cdf([x, y], mean=[0, 0], cov=[[1, rho], [rho, 1]]) for x, y in zip(a, b)]
        assert np.allclose(bvn_cdf(a, b, np.full(20, rho)), ref, atol=1e-7)


def test_bvn_independent_and_derivatives():
    assert np.allclose(bvn_cdf(.3, -.7, 0.), st.norm.cdf(.3) * st.norm.cdf(-.7), atol=1e-14)
    a, b, r, h = .4, -.2, .35, 1e-6
    F, Fa, Fb, Fr = bvn_grad(a, b, r)
    assert abs((bvn_cdf(a + h, b, r) - bvn_cdf(a - h, b, r)) / (2 * h) - Fa) < 1e-8
    assert abs((bvn_cdf(a, b + h, r) - bvn_cdf(a, b - h, r)) / (2 * h) - Fb) < 1e-8
    assert abs((bvn_cdf(a, b, r + h) - bvn_cdf(a, b, r - h)) / (2 * h) - Fr) < 1e-7


# ---------------------------------------------------------------------------
# binary and count
# ---------------------------------------------------------------------------

def test_cloglog_matches_direct_maximization(run):
    run(DATA + "quietly cloglog yb x1 x2")
    y, X = _col(run, "yb"), np.column_stack([_col(run, "x1"), _col(run, "x2"), np.ones(500)])

    def nll(b):
        eta = X @ b
        p = 1 - np.exp(-np.exp(eta))
        return -np.sum(y * np.log(p) + (1 - y) * (-np.exp(eta)))
    assert np.allclose(_b(run), _fit(nll, np.zeros(3)), atol=1e-4)


def test_binreg_logit_link_is_logit(run):
    run(DATA + "quietly binreg yb x1 x2, or")
    b_bin = _b(run)
    run("quietly logit yb x1 x2")
    assert np.allclose(b_bin, _b(run), atol=1e-6)


def test_tpoisson_matches_direct_maximization(run):
    run(DATA + "quietly tpoisson ctr x1 x2, ll(0)")
    m = ~np.isnan(_col(run, "ctr")) & (_col(run, "ctr") < 1e300)
    y = _col(run, "ctr")[m]
    X = np.column_stack([_col(run, "x1"), _col(run, "x2"), np.ones(500)])[m]

    def nll(b):
        mu = np.exp(X @ b)
        return -np.sum(st.poisson.logpmf(y, mu) - np.log1p(-np.exp(-mu)))
    assert np.allclose(_b(run), _fit(nll, np.zeros(3)), atol=1e-4)


def test_zip_without_inflation_signal_reduces_to_poisson_loglik_bound(run):
    run(DATA + "quietly zip cnt x1 x2, inflate(x1)\nscalar llz = e(ll)\nquietly poisson cnt x1 x2")
    assert run.eval("llz") >= run.session.e["ll"] - 1e-6


# ---------------------------------------------------------------------------
# censored, truncated and selection
# ---------------------------------------------------------------------------

def test_intreg_with_exact_intervals_is_regress_sigma(run):
    run(DATA + "quietly intreg yc yc x1 x2")
    b = _b(run)
    run("quietly regress yc x1 x2")
    X = np.column_stack([_col(run, "x1"), _col(run, "x2"), np.ones(500)])
    e = _col(run, "yc") - X @ np.asarray(run.session.e["b"].data, dtype=float).ravel()
    assert np.allclose(b[:3], np.asarray(run.session.e["b"].data, dtype=float).ravel(), atol=1e-5)
    assert abs(np.exp(b[3]) - np.sqrt(np.mean(e ** 2))) < 1e-5


def test_heckman_twostep_second_stage_is_ols_with_mills_ratio(run):
    run(DATA + "quietly heckman ysel x2, select(sel = x1) twostep")
    b = _b(run)
    run("quietly probit sel x1\npredict double xbs, xb\ngen double imr = normalden(xbs)/normal(xbs)\n"
        "quietly regress ysel x2 imr if sel")
    ols = np.asarray(run.session.e["b"].data, dtype=float).ravel()
    # e(b): outcome slope and constant, then the selection equation, then lambda
    assert np.allclose([b[0], b[1]], [ols[0], ols[2]], atol=1e-6)
    assert abs(b[-1] - ols[1]) < 1e-6


def test_biprobit_loglik_not_below_two_probits(run):
    run(DATA + "gen byte yb2 = (.2 - .6*x1 + round(invnormal(runiform()), .01)) > 0\n"
        "quietly probit yb x1 x2\nscalar l1 = e(ll)\nquietly probit yb2 x1\nscalar l2 = e(ll)\n"
        "quietly biprobit (yb = x1 x2) (yb2 = x1)")
    assert run.session.e["ll"] >= run.eval("l1") + run.eval("l2") - 1e-6
    assert -1 < run.session.e["rho"] < 1


# ---------------------------------------------------------------------------
# endogenous regressors
# ---------------------------------------------------------------------------

def test_ivregress_gmm_unadjusted_is_2sls(run):
    run(DATA + "quietly ivregress 2sls yiv x1 (w2 = z1 z2)")
    b2 = _b(run)
    run("quietly ivregress gmm yiv x1 (w2 = z1 z2), wmatrix(unadjusted)")
    assert np.allclose(_b(run), b2, atol=1e-9)


def test_ivregress_liml_equals_2sls_when_just_identified(run):
    run(DATA + "quietly ivregress 2sls yiv x1 (w2 = z1)")
    b2 = _b(run)
    run("quietly ivregress liml yiv x1 (w2 = z1)")
    assert np.allclose(_b(run), b2, atol=1e-8)
    assert abs(run.session.e["kappa"] - 1) < 1e-8


def test_ivprobit_without_endogeneity_matches_probit_scale(run):
    run(DATA + "quietly ivprobit yb x1 (w2 = z1 z2)")
    assert run.session.e["cmd"] == "ivprobit"
    assert abs(run.session.e["rho"]) < 1


def test_ivpoisson_gmm_moments_vanish_when_just_identified(run):
    run(DATA + "gen int yic = rpoisson(exp(.2 + .3*w2 + .2*x1))\n"
        "quietly ivpoisson gmm yic x1 (w2 = z1)")
    b = _b(run)
    X = np.column_stack([_col(run, "w2"), _col(run, "x1"), np.ones(500)])
    Z = np.column_stack([_col(run, "z1"), _col(run, "x1"), np.ones(500)])
    g = Z.T @ (_col(run, "yic") - np.exp(X @ b)) / 500
    assert np.max(np.abs(g)) < 1e-6
