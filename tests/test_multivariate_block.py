"""Phase 20: multivariate statistics (components, factors, rotation, reliability, multivariate tests
and models, clusters, discriminant analysis, correspondence analysis, scaling, procrustes). The
Stata reference is compat 1906; here the checks are identities between commands and independent
computations with numpy and scipy."""

import numpy as np
from scipy import stats as st
from scipy.cluster import hierarchy as H
from scipy.spatial.distance import pdist, squareform

MV = """
clear
set seed 1906
set obs 150
gen int id = _n
gen byte g = ceil(3*runiform())
gen byte g2 = 1 + (runiform() < .5)
gen double f1 = round(invnormal(runiform()), .01)
gen double f2 = round(invnormal(runiform()), .01)
gen double x1 = round(.8*f1 + .6*invnormal(runiform()) + .6*g, .01)
gen double x2 = round(.7*f1 + .7*invnormal(runiform()) - .4*g, .01)
gen double x3 = round(.6*f1 + .8*invnormal(runiform()) + .3*g, .01)
gen double x4 = round(.8*f2 + .6*invnormal(runiform()), .01)
gen double x5 = round(.7*f2 + .7*invnormal(runiform()), .01)
gen double x6 = round(.6*f2 + .8*invnormal(runiform()), .01)
gen double y1 = round(1 + .5*x1 + .3*x4 + invnormal(runiform()), .01)
gen double y2 = round(.5*x2 - .4*x5 + invnormal(runiform()), .01)
gen double y3 = round(.3*x3 + .2*g + invnormal(runiform()), .01)
gen byte a = ceil(3*runiform())
gen byte b = ceil(4*runiform())
replace b = a if runiform() < .4
"""


def _data(run, *names):
    d = run.session.data
    return np.column_stack([d.get(n).data.astype(float) for n in names])


def _mat(run, name):
    return np.asarray(run.session.e[name].data, dtype=float)


# ---------------------------------------------------------------------------
# components, factors, rotation, reliability
# ---------------------------------------------------------------------------

def test_pca_eigenvalues_are_those_of_the_correlation_matrix(run):
    run(MV + "quietly pca x1 x2 x3 x4 x5 x6")
    X = _data(run, "x1", "x2", "x3", "x4", "x5", "x6")
    ev = np.sort(np.linalg.eigvalsh(np.corrcoef(X.T)))[::-1]
    assert np.allclose(_mat(run, "Ev").ravel(), ev, atol=1e-9)
    assert run.session.e["trace"] == 6


def test_pca_covariance_and_loadings_are_orthonormal(run):
    run(MV + "quietly pca x1 x2 x3 x4 x5 x6, covariance")
    X = _data(run, "x1", "x2", "x3", "x4", "x5", "x6")
    ev = np.sort(np.linalg.eigvalsh(np.cov(X.T)))[::-1]
    assert np.allclose(_mat(run, "Ev").ravel(), ev, atol=1e-9)
    L = _mat(run, "L")
    assert np.allclose(L.T @ L, np.eye(6), atol=1e-9)


def test_pca_scores_are_uncorrelated_with_variance_equal_to_eigenvalue(run):
    run(MV + "quietly pca x1 x2 x3 x4 x5 x6, components(2)\nquietly predict double p1 p2")
    P = _data(run, "p1", "p2")
    ev = _mat(run, "Ev").ravel()
    assert np.allclose(P.var(0, ddof=1), ev[:2], atol=1e-9)
    assert abs(np.corrcoef(P.T)[0, 1]) < 1e-9


def test_factor_pcf_equals_pca_loadings(run):
    run(MV + "quietly pca x1 x2 x3 x4 x5 x6, components(2)")
    ev = _mat(run, "Ev").ravel()
    V = _mat(run, "L")[:, :2]
    out = run("factor x1 x2 x3 x4 x5 x6, pcf factors(2)")
    assert "Factor analysis/correlation" in out
    assert np.allclose(np.abs(_mat(run, "L")[:, :2]), np.abs(V * np.sqrt(ev[:2])), atol=1e-8)


def test_factor_ml_stores_the_likelihood_ratio_tests(run):
    run(MV + "quietly factor x1 x2 x3 x4 x5 x6, ml factors(2)")
    e = run.session.e
    assert e["df_m"] == 11 and e["df_r"] == 4 and e["df_1"] == 4
    assert e["chi2_i"] > 100 and e["p_i"] < 1e-10
    assert np.isclose(e["chi2_1"], 6.3557, atol=1e-3)


def test_rotation_preserves_communalities(run):
    run(MV + "quietly factor x1 x2 x3 x4 x5 x6, pf factors(2)")
    L0 = _mat(run, "L")[:, :2]
    run("quietly rotate")
    Lr = _mat(run, "r_L")
    T = _mat(run, "r_T")
    assert np.allclose((L0 ** 2).sum(1), (Lr ** 2).sum(1), atol=1e-8)
    assert np.allclose(T.T @ T, np.eye(2), atol=1e-8)
    assert np.allclose(L0 @ T, Lr, atol=1e-8)


def test_alpha_matches_the_covariance_formula(run):
    run(MV + "quietly alpha x1 x2 x3 x4")
    S = np.cov(_data(run, "x1", "x2", "x3", "x4").T)
    k = 4
    ref = k / (k - 1) * (1 - np.trace(S) / S.sum())
    assert np.isclose(run.session.r["alpha"], ref, atol=1e-9)
    run("quietly alpha x1 x2 x3 x4, std")
    R = np.corrcoef(_data(run, "x1", "x2", "x3", "x4").T)
    rbar = (R.sum() - k) / (k * (k - 1))
    assert np.isclose(run.session.r["alpha"], k * rbar / (1 + (k - 1) * rbar), atol=1e-9)


def test_core_varimax_is_orthogonal_and_improves_the_criterion():
    from opendta.core import multivariate as mv
    rng = np.random.default_rng(3)
    L = rng.normal(size=(8, 3))
    rotated, T, _ = mv.orthomax(L, 1.0)
    assert np.allclose(T.T @ T, np.eye(3), atol=1e-8)
    assert np.allclose(L @ T, rotated, atol=1e-8)
    # varimax criterion: sum over columns of the variance of the squared loadings
    crit = lambda A: ((A ** 2).var(0)).sum()
    assert crit(rotated) >= crit(L) - 1e-12


# ---------------------------------------------------------------------------
# multivariate tests and models
# ---------------------------------------------------------------------------

def test_hotelling_two_sample_matches_the_formula(run):
    run(MV + "quietly hotelling x1 x2 x3, by(g2)")
    X, g = _data(run, "x1", "x2", "x3"), _data(run, "g2").ravel()
    A, B = X[g == 1], X[g == 2]
    n1, n2 = len(A), len(B)
    S = ((n1 - 1) * np.cov(A.T) + (n2 - 1) * np.cov(B.T)) / (n1 + n2 - 2)
    d = A.mean(0) - B.mean(0)
    t2 = n1 * n2 / (n1 + n2) * d @ np.linalg.solve(S, d)
    assert np.isclose(run.session.r["T2"], t2, atol=1e-8)
    assert run.session.r["N"] == n1 + n2 and run.session.r["k"] == 3 and run.session.r["df"] == n1 + n2 - 3 - 1


def test_mvreg_equals_regress_equation_by_equation(run):
    run(MV + "quietly mvreg y1 y2 y3 = x1 x2")
    b = np.asarray(run.session.e["b"].data, dtype=float).ravel()
    X = np.column_stack([_data(run, "x1", "x2"), np.ones(150)])
    for k, y in enumerate(("y1", "y2", "y3")):
        beta = np.linalg.lstsq(X, _data(run, y).ravel(), rcond=None)[0]
        assert np.allclose(b[3 * k:3 * k + 3], beta, atol=1e-8)


def test_manova_wilks_matches_hotelling_for_two_groups(run):
    out = run(MV + "quietly hotelling y1 y2 y3, by(g2)\nlocal t2 = r(T2)\nmanova y1 y2 y3 = g2")
    assert "Wilks' lambda" in out
    y, g = _data(run, "y1", "y2", "y3"), _data(run, "g2").ravel()
    A, B = y[g == 1], y[g == 2]
    n1, n2 = len(A), len(B)
    S = ((n1 - 1) * np.cov(A.T) + (n2 - 1) * np.cov(B.T)) / (n1 + n2 - 2)
    d = A.mean(0) - B.mean(0)
    t2 = n1 * n2 / (n1 + n2) * d @ np.linalg.solve(S, d)
    wilks = 1.0 / (1.0 + t2 / (n1 + n2 - 2))
    assert f"{wilks:.4f}" in out


def test_canon_correlations_are_the_square_roots_of_the_eigenvalues(run):
    run(MV + "quietly canon (x1 x2 x3) (y1 y2 y3)")
    X, Y = _data(run, "x1", "x2", "x3"), _data(run, "y1", "y2", "y3")
    Z = np.column_stack([X, Y])
    R = np.corrcoef(Z.T)
    Rxx, Ryy, Rxy = R[:3, :3], R[3:, 3:], R[:3, 3:]
    M = np.linalg.solve(Rxx, Rxy) @ np.linalg.solve(Ryy, Rxy.T)
    rho = np.sqrt(np.sort(np.linalg.eigvals(M).real)[::-1])
    assert np.allclose(_mat(run, "ccorr").ravel(), rho, atol=1e-8)


def test_mvtest_normality_stores_the_p_values(run):
    run(MV + "quietly mvtest normality x1 x2 x3, stats(all)")
    r = run.session.r
    assert all(0 < r[k] < 1 for k in ("p_skew", "p_kurt", "p_hz", "p_dh"))


# ---------------------------------------------------------------------------
# clusters
# ---------------------------------------------------------------------------

def test_core_linkage_heights_match_scipy():
    from opendta.core import cluster as K
    rng = np.random.default_rng(5)
    X = rng.normal(size=(25, 3))
    D = squareform(pdist(X))
    for ours, theirs in (("single", "single"), ("complete", "complete"), ("average", "average"), ("waverage", "weighted")):
        h = [m[2] for m in K.linkage(D, ours)]
        ref = H.linkage(pdist(X), theirs)[:, 2]
        assert np.allclose(sorted(h), sorted(ref), atol=1e-9)
    D2 = squareform(pdist(X, "sqeuclidean"))
    h = [m[2] for m in K.linkage(D2, "ward")]
    ref = H.linkage(pdist(X), "ward")[:, 2] ** 2 / 2
    assert np.allclose(sorted(h), sorted(ref), atol=1e-9) or np.allclose(sorted(h), sorted(ref * 2), atol=1e-9)


def test_cluster_kmeans_generate_and_list(run):
    run(MV + "quietly cluster kmeans x1 x2 x3, k(3) name(km) start(firstk)")
    km = _data(run, "km").ravel()
    assert set(km) == {1.0, 2.0, 3.0}
    out = run("cluster dir")
    assert "km" in out
    run("quietly cluster averagelinkage x1 x2 x3, name(av)\nquietly cluster generate av3 = groups(3), name(av)")
    av3 = _data(run, "av3").ravel()
    X = _data(run, "x1", "x2", "x3")
    ref = H.fcluster(H.linkage(pdist(X), "average"), 3, "maxclust")
    # same partition up to label permutation
    pairs = set(zip(av3.astype(int), ref))
    assert len(pairs) == 3


def test_matrix_dissimilarity_l1(run):
    from opendta.commands.matrix import get_matrix
    run(MV + "quietly matrix dissimilarity D = x1 x2 x3 in 1/5, L1")
    X = _data(run, "x1", "x2", "x3")[:5]
    assert np.allclose(get_matrix(run.session, "D").data, squareform(pdist(X, "cityblock")), atol=1e-9)


# ---------------------------------------------------------------------------
# discriminant analysis
# ---------------------------------------------------------------------------

def test_discrim_lda_classifies_by_the_largest_posterior(run):
    out = run(MV + "discrim lda x1 x2 x3, group(g)\nquietly predict cls\nquietly predict p1 p2 p3, pr")
    assert "Linear discriminant analysis" in out and "Resubstitution classification summary" in out
    P = _data(run, "p1", "p2", "p3")
    cls = _data(run, "cls").ravel()
    assert np.allclose(P.sum(1), 1.0)
    assert np.array_equal(cls, P.argmax(1) + 1)
    X, g = _data(run, "x1", "x2", "x3"), _data(run, "g").ravel()
    means = np.array([X[g == k].mean(0) for k in (1, 2, 3)])
    S = sum((X[g == k] - means[i]).T @ (X[g == k] - means[i]) for i, k in enumerate((1, 2, 3))) / (150 - 3)
    d2 = np.array([[(x - m) @ np.linalg.solve(S, x - m) for m in means] for x in X])
    assert np.array_equal(cls, d2.argmin(1) + 1)


def test_discrim_lda_error_rate_matches_the_table(run):
    run(MV + "quietly discrim lda x1 x2 x3, group(g)\nquietly estat errorrate")
    rate = np.asarray(run.session.r["Rate"].data, dtype=float).ravel()
    run("quietly predict cls")
    cls, g = _data(run, "cls").ravel(), _data(run, "g").ravel()
    per = [np.mean(cls[g == k] != k) for k in (1, 2, 3)]
    assert np.allclose(rate[:3], per, atol=1e-12)
    assert np.isclose(rate[3], np.mean(per), atol=1e-12)


def test_candisc_eigenvalues_are_those_of_w_inverse_b(run):
    run(MV + "quietly candisc x1 x2 x3, group(g)")
    X, g = _data(run, "x1", "x2", "x3"), _data(run, "g").ravel()
    grand = X.mean(0)
    B = sum((g == k).sum() * np.outer(X[g == k].mean(0) - grand, X[g == k].mean(0) - grand) for k in (1, 2, 3))
    W = sum((X[g == k] - X[g == k].mean(0)).T @ (X[g == k] - X[g == k].mean(0)) for k in (1, 2, 3))
    ev = np.sort(np.linalg.eigvals(np.linalg.solve(W, B)).real)[::-1][:2]
    assert np.allclose(_mat(run, "Ev").ravel(), ev, atol=1e-8)


def test_discrim_knn_k1_resubstitution_is_perfect(run):
    run(MV + "quietly discrim knn x1 x2 x3, group(g) k(1)")
    assert run.session.e["k_nn"] == 1
    out = run("estat classtable")
    assert "Resubstitution classification table" in out
    # with k(1) every observation is its own nearest neighbor: perfect resubstitution
    run("quietly estat errorrate")
    assert np.allclose(np.asarray(run.session.r["Rate"].data, dtype=float), 0.0)


def test_discrim_logistic_classifies_better_than_chance(run):
    run(MV + "quietly discrim logistic x1 x2 x3, group(g)\nquietly predict cls")
    cls = _data(run, "cls").ravel()
    assert set(np.unique(cls)) <= {1.0, 2.0, 3.0}
    # the multinomial logit classifies at least as well as chance
    assert np.mean(cls == _data(run, "g").ravel()) > 0.5


# ---------------------------------------------------------------------------
# correspondence analysis, scaling, procrustes
# ---------------------------------------------------------------------------

def test_ca_total_inertia_is_chi2_over_n(run):
    out = run(MV + "ca a b")
    d = run.session.data
    tab = np.zeros((3, 4))
    for i, j in zip(d.get("a").data.astype(int), d.get("b").data.astype(int)):
        tab[i - 1, j - 1] += 1
    chi2 = st.chi2_contingency(tab, correction=False)[0]
    assert np.isclose(run.session.e["inertia"], chi2 / 150, atol=1e-10)
    assert np.isclose(run.session.e["X2"], chi2, atol=1e-8)
    assert "Correspondence analysis" in out


def test_ca_principal_inertias_sum_to_the_total():
    from opendta.core import ordination as O
    N = np.array([[10, 5, 2.0], [3, 12, 6], [1, 4, 15]])
    f = O.ca_fit(N)
    assert np.isclose(f["inertia"].sum(), f["total"])
    R, C = O.ca_coordinates(f, "principal")
    # barycentric relation: row principal coordinates are the weighted averages of column standard ones
    P = N / N.sum()
    prof = P / P.sum(1, keepdims=True)
    assert np.allclose(R, prof @ f["B"], atol=1e-9)


def test_mca_adjusted_inertias_follow_benzecri():
    from opendta.core import ordination as O
    rng = np.random.default_rng(2)
    codes = [rng.integers(0, 3, 80), rng.integers(0, 2, 80), rng.integers(0, 4, 80)]
    f = O.mca_fit(codes, [3, 2, 4])
    lam, Q = f["inertia"], 3
    keep = lam > 1 / Q
    assert np.allclose(f["adjusted"][keep], (Q / (Q - 1) * (lam[keep] - 1 / Q)) ** 2)
    assert np.all(f["adjusted"][~keep] == 0)


def test_mds_classical_reproduces_the_distances_with_all_dimensions():
    from opendta.core import ordination as O
    rng = np.random.default_rng(1)
    X = rng.normal(size=(9, 3))
    D = squareform(pdist(X))
    r = O.mds_classical(D, 3)
    assert np.allclose(squareform(pdist(r["Y"])), D, atol=1e-9)
    m1, m2 = O.mardia(r["ev"], 3)
    assert np.isclose(m1, 1.0) and np.isclose(m2, 1.0)


def test_mds_command_stores_the_configuration(run):
    out = run(MV + "mds x1 x2 x3 in 1/30, id(id) noplot")
    assert "Classical metric multidimensional scaling" in out
    e = run.session.e
    assert e["N"] == 30 and e["p"] == 2
    Y = _mat(run, "Y")
    X = _data(run, "x1", "x2", "x3")[:30]
    D = squareform(pdist(X))
    full = squareform(pdist(Y))
    assert np.all(full <= D + 1e-9)          # a projection never lengthens a distance


def test_procrustes_recovers_a_rotation():
    from opendta.core import ordination as O
    rng = np.random.default_rng(4)
    Y = rng.normal(size=(30, 2))
    th = 0.9
    Rm = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    X = 1.7 * Y @ Rm.T + np.array([2.0, -1.0])
    r = O.procrustes_orthogonal(Y, X)
    assert r["rss"] < 1e-20
    assert np.isclose(r["c"], 1 / 1.7)


def test_procrustes_command_dilation_is_one_with_norho(run):
    run(MV + "quietly procrustes (x1 x2) (x4 x5), norho")
    assert run.session.e["rho"] == 1
    T = _mat(run, "A")
    assert np.allclose(T.T @ T, np.eye(2), atol=1e-9)


def test_biplot_runs_and_core_reconstructs_the_data():
    from opendta.core import ordination as O
    rng = np.random.default_rng(6)
    X = rng.normal(size=(20, 4))
    b = O.biplot_coords(X, std=False, alpha=0.5)
    assert np.allclose(b["rows"] @ b["cols"].T, X - X.mean(0), atol=1e-9)


# ---------------------------------------------------------------------------
# discrim core
# ---------------------------------------------------------------------------

def test_core_discrim_lda_loo_matches_a_manual_refit():
    from opendta.core import discrim as DC
    rng = np.random.default_rng(7)
    X = np.vstack([rng.normal(size=(15, 2)), rng.normal(size=(15, 2)) + 1.5])
    codes = np.repeat([0, 1], 15)
    loo = DC.loo_scores("lda", X, codes, 2, None, "equal")
    i = 4
    keep = np.arange(30) != i
    st_ = DC.group_stats(X[keep], codes[keep], 2)
    sc, _ = DC.lda_scores(X[i:i + 1], st_, np.array([0.5, 0.5]))
    assert np.allclose(loo[i], sc[0])
