"""Block 23a: rotate and rotatemat with the gradient projection criteria, pca vce(normal) and
means, estat factors and rotatecompare, the estat subcommands of canon, mvtest covariances
(compound, block, equals) and correlations (equals, Jennrich), mvtest normality univariate and
bivariate, manova with error terms and manovatest test().

The numbers marked "manual" come from the examples of the [MV] manual whose data are printed
there (the 24 x 4 loading matrix of [MV] rotatemat, the correlation matrix of [MV] pca); the
other tests check identities between results."""

import numpy as np
import pytest

from opendta.core import multivariate as mv
from opendta.core import mvtests as T
from opendta.core import rotation as R

# [MV] rotatemat: loadings of 24 psychological tests on 4 factors (Harman)
L24 = np.array([601, 19, 388, 221, 372, -25, 252, 132, 413, -117, 388, 144, 487, -100, 254, 192, 691, -304, -279, 35,
                690, -409, -200, -76, 677, -409, -292, 84, 674, -189, -99, 122, 697, -454, -212, -80, 476, 534, -486, 92,
                558, 332, -142, -90, 472, 508, -139, 256, 602, 244, 28, 295, 423, 58, 15, -415, 394, 89, 97, -362,
                510, 95, 347, -249, 466, 197, -4, -381, 515, 312, 152, -147, 443, 89, 109, -150, 614, -118, 126, -38,
                589, 227, 57, 123, 608, -107, 127, -38, 687, -44, 138, 98, 651, 177, -212, -17], float).reshape(24, 4) / 1000

L24_DO = ("matrix L = (" + " \\ ".join(",".join(f"{v:g}" for v in row) for row in L24) + ")\n"
          "matrix colnames L = F1 F2 F3 F4\n")

MV = """clear
set seed 2329
set obs 200
gen f1 = rnormal()
gen f2 = rnormal()
gen g = 1 + mod(_n, 3)
gen h = 1 + mod(_n, 2)
forvalues i = 1/4 {
    gen x`i' = f1 + 0.4*f2 + rnormal()*.6
}
forvalues i = 5/8 {
    gen x`i' = f2 + 0.3*f1 + rnormal()*.6
}
gen y1 = x1 + g + rnormal()
gen y2 = x5 - h + rnormal()
"""


def _m(run, name):
    return run.session.r[name].data if name in run.session.r else run.session.e[name].data


# ---------------------------------------------------------------------------
# core/rotation
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("name,arg", [("quartimax", None), ("oblimin", 0.5), ("quartimin", None), ("cf", 0.3),
                                      ("entropy", None), ("tandem1", None), ("tandem2", None), ("bentler", None),
                                      ("oblimax", None), ("vgpf", None), ("parsimax", None), ("equamax", None)])
def test_criterion_gradients_match_finite_differences(name, arg):
    L = L24[:8, :3] + 0.01
    f, G = R.criterion(name, L, arg)
    num = np.zeros_like(L)
    h = 1e-6
    for i in range(L.shape[0]):
        for j in range(L.shape[1]):
            Lp, Lm = L.copy(), L.copy()
            Lp[i, j] += h
            Lm[i, j] -= h
            num[i, j] = (R.criterion(name, Lp, arg)[0] - R.criterion(name, Lm, arg)[0]) / (2 * h)
    assert np.allclose(G, num, atol=1e-5)


def test_target_and_partial_gradients():
    L = L24[:6, :2]
    H = np.ones_like(L) * 0.5
    W = (np.arange(12).reshape(6, 2) % 2).astype(float)
    for name, arg in (("target", H), ("partial", (H, W))):
        f, G = R.criterion(name, L, arg)
        Lp = L.copy()
        Lp[2, 1] += 1e-6
        assert G[2, 1] == pytest.approx((R.criterion(name, Lp, arg)[0] - f) / 1e-6, abs=1e-5)


def test_gpf_reproduces_the_rotatemat_examples():
    # manual: quartimax criterion -1.032898 after 35 iterations; oblique oblimin .1957363 after 78
    q = R.gpf(L24, "quartimax")
    assert q["f"] == pytest.approx(-1.032898, abs=5e-7) and q["iter"] == 35
    o = R.gpf(L24, "oblimin", 0.0, oblique=True)
    assert o["f"] == pytest.approx(0.1957363, abs=5e-8) and o["iter"] == 78
    L, Tm = R.orient(o["L"], o["T"])
    assert np.allclose(Tm[0], [0.823, 0.715, 0.584, 0.699], atol=6e-4)
    # manual: r(T)'r(T) has .4026978 below the diagonal in row 2
    assert (Tm.T @ Tm)[1, 0] == pytest.approx(0.4026978, abs=2e-6)
    assert np.allclose(L, L24 @ np.linalg.inv(Tm).T)


def test_orthogonal_criteria_that_the_manual_calls_equivalent():
    p, k = L24.shape
    v = R.orient(*[R.gpf(L24, "vgpf")[x] for x in ("L", "T")])[0]
    c = R.orient(*[R.gpf(L24, "cf", 1.0 / p)[x] for x in ("L", "T")])[0]
    q = R.orient(*[R.gpf(L24, "quartimax")[x] for x in ("L", "T")])[0]
    ox = R.orient(*[R.gpf(L24, "oblimax")[x] for x in ("L", "T")])[0]
    assert np.allclose(v, c, atol=1e-4)            # varimax = cf(1/p)
    assert np.allclose(q, ox, atol=1e-4)           # quartimax = oblimax for orthogonal rotations


# ---------------------------------------------------------------------------
# rotatemat and rotate
# ---------------------------------------------------------------------------

def test_rotatemat_varimax_promax_and_quartimax_output(run):
    out = run(L24_DO + "rotatemat L, format(%6.3f)")
    assert "Criterion               varimax" in out and "Rotation class          orthogonal" in out
    T_ = _m(run, "T")
    assert np.allclose(T_[0], [0.677, 0.438, 0.475, 0.352], atol=6e-4)        # manual
    assert np.allclose(_m(run, "AT")[0], [0.247, 0.151, 0.679, 0.128], atol=6e-4)
    out = run("rotatemat L, quartimax norotation blanks(0.3) format(%6.3f)")
    assert "Criterion value         -1.032898" in out and "Number of iterations    35" in out
    assert "Orthogonal rotation" not in out and "(blanks represent abs()<.3)" in out
    assert run.session.r["ctitle"] == "quartimax" and run.session.r["iter"] == 35
    run("rotatemat L, promax")
    assert np.allclose(_m(run, "T")[0], [0.841, 0.829, 0.663, 0.743], atol=6e-4)       # manual
    assert run.session.r["ctitle"] == "promax(3)" and run.session.r["class"] == "oblique"


def test_rotatemat_errors(run):
    run(L24_DO)
    run("rotatemat L, entropy oblique")
    assert run.rc == 198
    run("rotatemat L, varimax quartimax")
    assert run.rc == 198
    run("matrix Tg = J(3,4,0)\nrotatemat L, target(Tg)")
    assert run.rc == 503


def test_rotate_after_factor_with_new_criteria(run):
    run(MV + "quietly factor x1-x8, ipf factors(2)")
    out = run("rotate, oblimin")
    # oblimin without oblique is orthogonal ([MV] rotate: orthogonal is the default but for promax)
    assert "orthogonal oblimin (Kaiser off)" in out
    assert np.allclose(_m(run, "r_Phi"), np.eye(2))
    out = run("rotate, cf(0.5) oblique")
    assert "oblique cf (Kaiser off)" in out and "Rotated factors are correlated" in out
    assert run.session.e["r_ctitle"] == "Crawford-Ferguson(0.5)" and "r_fmin" in run.session.e
    run("matrix W = (1,0\\1,0\\1,0\\1,0\\0,1\\0,1\\0,1\\0,1)")
    out = run("rotate, target(W) normalize")
    assert "orthogonal target (Kaiser on)" in out
    Lr = _m(run, "r_L")
    assert (np.abs(Lr[:4, 0]) > np.abs(Lr[:4, 1])).all() and (np.abs(Lr[4:, 1]) > np.abs(Lr[4:, 0])).all()
    ev = _m(run, "r_Ev").ravel()
    assert ev[0] >= ev[1]                       # decreasing explained variance
    out = run("set seed 1\nrotate, oblimin oblique normalize protect(3)")
    assert out.count("Trial ") == 3


def test_rotate_entropy_is_orthogonal_only(run):
    run(MV + "quietly pca x1-x8, comp(3)")
    run("rotate, entropy oblique")
    assert run.rc == 198
    out = run("rotate, tandem2")
    assert "orthogonal tandem2" in out


# ---------------------------------------------------------------------------
# estat factors and estat rotatecompare
# ---------------------------------------------------------------------------

def test_estat_factors(run):
    run(MV + "quietly factor x1-x8, ml factors(2)")
    out = run("estat factors")
    assert "Factor analysis with different numbers of factors (maximum likelihood)" in out
    S = _m(run, "stats")
    p_ = 8
    assert S.shape == (4, 5)                    # 8 variables: up to 4 factors with df_r >= 0
    for k in range(1, 5):
        ll, dfm, dfr, aic, bic = S[k - 1]
        assert dfm == p_ * k - k * (k - 1) / 2 and dfr == ((p_ - k) ** 2 - (p_ + k)) / 2
        assert aic == pytest.approx(-2 * ll + 2 * dfm) and bic == pytest.approx(-2 * ll + dfm * np.log(200))
    assert S[1, 0] == pytest.approx(run.session.e["ll"], rel=1e-6)
    run("estat factors, factors(2)")
    assert _m(run, "stats").shape == (2, 5)


def test_estat_rotatecompare(run):
    run(MV + "quietly factor x1-x8, ipf factors(2)")
    run("estat rotatecompare")
    assert run.rc == 321
    run("quietly rotate, promax")
    out = run("estat rotatecompare")
    assert "Rotation matrix -- oblique promax (Kaiser off)" in out and "Unrotated" in out
    run(MV + "quietly pca x1-x8, comp(2)\nquietly rotate")
    out = run("estat rotatec")
    assert "Component loadings" in out


# ---------------------------------------------------------------------------
# pca vce(normal) and means
# ---------------------------------------------------------------------------

AUDIO = [[1], [.7775, 1], [.4012, .5366, 1], [.2554, .2749, .4250, 1], [.6963, .5515, .2391, .1790, 1],
         [.6416, .7070, .4460, .2632, .6634, 1], [.2372, .3597, .7011, .3165, .1589, .4142, 1],
         [.2041, .2169, .3262, .7097, .1321, .2201, .3746, 1]]


def test_pca_vce_against_the_audiometric_example():
    C = np.zeros((8, 8))
    for i, row in enumerate(AUDIO):
        for j, v in enumerate(row):
            C[i, j] = C[j, i] = v
    fit = mv.pca_fit(C)
    v = mv.pca_vce(fit["ev"], fit["V"], 100, 2)
    se = np.sqrt(np.diag(v["V"]))
    # manual (the correlations are printed with 4 decimals, hence the tolerances)
    assert se[0] == pytest.approx(0.5556453, abs=2e-4)
    assert v["b"][0] / se[0] == pytest.approx(np.sqrt(50))
    assert se[2] == pytest.approx(0.0429963, abs=2e-4)                 # lft500 in Comp1
    assert np.sqrt(v["v_rho"]) == pytest.approx(0.0273, abs=1e-4)
    assert v["se_prop"][0] == pytest.approx(0.0394, abs=1e-4) and v["se_prop"][1] == pytest.approx(0.0271, abs=1e-4)
    assert v["bias"][0] == pytest.approx(0.056663, abs=2e-4)
    lr = mv.pca_lr_tests(C, fit["ev"], 100)
    assert lr["chi2_i"] == pytest.approx(448.21, abs=0.2) and lr["df_i"] == 28
    assert lr["chi2_s"] == pytest.approx(451.11, abs=0.2) and lr["df_s"] == 35
    assert np.allclose(v["V"], v["V"].T)


def test_pca_vce_normal_and_means(run):
    out = run(MV + "pca x1-x4, comp(2) vce(normal) means")
    assert "Summary statistics of the variables" in out
    assert "(with PCA/correlation, SEs and tests are approximate)" in out
    assert "SEs assume multivariate normality" in out and "Explained variance by components" in out
    e = run.session.e
    assert e["properties"] == "b V eigen" and e["vce"] == "multivariate normality"
    assert e["b"].data.shape == (1, 2 + 2 * 4) and e["b"].coleq[0] == "Eigenvalues" and e["b"].coleq[2] == "Comp1"
    assert e["Ev_stats"].data.shape == (4, 5)
    out = run("pca x1-x4, comp(2) vce(normal) novce")
    assert "Principal components (eigenvectors)" in out and "Std. Err." not in out
    run("pca x1-x4, vce(robust)")
    assert run.rc == 198


# ---------------------------------------------------------------------------
# canon estat
# ---------------------------------------------------------------------------

def test_canon_estat_loadings_and_rotation(run):
    run(MV + "quietly canon (x1 x2 x3) (x5 x6 y2)")
    out = run("estat loadings")
    assert "Canonical loadings for variable list 1" in out
    assert "Correlation between variable list 2 and canonical variates from list 1" in out
    L11, L22 = _m(run, "canload11"), _m(run, "canload22")
    run("estat rotatecompare")
    assert run.rc == 321
    out = run("estat rotate, loadings")
    assert "Rotated canonical loadings" in out and "Rotation matrix" in out
    T_, AT = _m(run, "T"), _m(run, "AT")
    assert np.allclose(T_.T @ T_, np.eye(3), atol=1e-10)
    assert np.allclose(AT, np.vstack([L11, L22]) @ T_)
    assert run.session.r["coefficients"] == "canonical loadings" and run.session.r["criterion"] == "varimax"
    out = run("estat rotatecompare")
    assert "Unrotated canonical loadings" in out
    run("estat rotate, stdcoefs")
    assert run.session.r["coefficients"] == "standardized coefficients"
    out = run("estat correlations")
    assert "Correlations between variable lists 1 and 2" in out
    assert _m(run, "corr_var1").shape == (3, 3)


def test_canon_sign_rule_second_set_standardized(run):
    # the largest standardized coefficient of each second-set column is positive (the signs of
    # the auto example of [MV] canon)
    run(MV + "quietly canon (x1 x2 x3 y1) (x5 x6 x7 y2)")
    Z2 = _m(run, "stdcoef_var2")
    assert all(Z2[np.argmax(np.abs(Z2[:, j])), j] > 0 for j in range(Z2.shape[1]))


# ---------------------------------------------------------------------------
# mvtest
# ---------------------------------------------------------------------------

def test_block_test_with_single_variables_is_the_diagonal_test():
    rng = np.random.default_rng(1)
    X = rng.standard_normal((80, 4)) @ np.array([[1, .3, 0, 0], [0, 1, .2, 0], [0, 0, 1, .4], [0, 0, 0, 1]])
    S = np.cov(X, rowvar=False, bias=True)
    R_ = np.corrcoef(X.T)
    a = T.lr_block_cov(S, [[0], [1], [2], [3]], 80)
    b = T.lr_diagonal(R_, 80)
    assert a["chi2"] == pytest.approx(b["chi2"]) and a["df"] == b["df"]
    assert T.lr_equals_cov(S, S * 80 / 79, 80)["chi2"] == pytest.approx(0.0, abs=1e-10)
    assert T.jennrich_one(R_, R_, 80)["chi2"] == pytest.approx(0.0, abs=1e-10)
    assert T.jennrich_equal([R_, R_], [40, 40])["chi2"] == pytest.approx(0.0, abs=1e-10)


def test_mvtest_covariances_and_correlations(run):
    run(MV)
    out = run("mvtest covariances x1 x2 x5, compound")
    assert "Test that covariance matrix is compound symmetric" in out and run.session.r["df"] == 4
    out = run("mvtest covariances x1 x2 x5 x6, block(x1 x2 || x5 x6)")
    assert "block diagonal" in out and run.session.r["df"] == 4 and run.session.r["chi2type"] == "Adjusted LR"
    run("quietly correlate x1 x2 x5, covariance\nmatrix S0 = r(C)")
    run("mvtest covariances x1 x2 x5, equals(S0)")
    assert run.session.r["chi2"] == pytest.approx(0.0, abs=1e-8) and run.session.r["df"] == 6
    run("mvtest covariances x1 x2, by(g)")
    assert {"F_Box", "df_m_Box", "df_r_Box", "p_F_Box"} <= set(run.session.r)
    run("matrix C = I(3)")
    out = run("mvtest correlations x1 x2 x5, equals(C)")
    assert "Jennrich chi2(3)" in out and run.session.r["chi2type"] == "Jennrich"
    out = run("mvtest correlations x1 x2 x5, by(g h)")
    assert "Test of equality of correlation matrices across samples" in out and run.session.r["df"] == 15
    run("mvtest covariances x1 x2, compound spherical")
    assert run.rc == 198


def test_mvtest_normality_univariate_and_bivariate(run):
    run(MV)
    out = run("mvtest normality x1 x2 x3, univariate bivariate stats(all)")
    assert out.index("Test for univariate normality") < out.index("Doornik-Hansen test for bivariate normality") \
        < out.index("Mardia mSkewness") < out.index("Henze-Zirkler") < out.index("Doornik-Hansen  ")
    U, B = _m(run, "Utest"), _m(run, "Btest")
    assert U.shape == (3, 4) and B.shape == (3, 3)
    run("quietly sktest x2")
    assert U[1, 3] == pytest.approx(run.session.r["P_chi2"])
    X = np.column_stack([run.session.data.get(v).data for v in ("x1", "x3")])
    assert B[1, 0] == pytest.approx(T.doornik_hansen(X)["chi2"])
    run("mvtest normality x1 x2 x3, stats(hz kurt)")
    r = run.session.r
    assert {"hz", "E_hz", "V_hz", "z_hz", "p_hz", "mkurt", "z_mkurt"} <= set(r) and "p_dh" not in r
    assert r["p_hz"] == pytest.approx(2 * T.st.norm.sf(abs(r["z_hz"])))


# ---------------------------------------------------------------------------
# manova error terms and manovatest test()
# ---------------------------------------------------------------------------

NEST = """clear
set seed 77
set obs 48
gen a = 1 + mod(_n, 2)
gen b = 1 + mod(floor((_n - 1) / 2), 3)
gen c = 1 + mod(floor((_n - 1) / 6), 2)
gen u = rnormal() + a + 0.3*b
gen v = rnormal() - a + 0.2*c
"""


def test_manova_error_term_matches_manovatest(run):
    run(NEST)
    out = run("manova u v = a / b|a /")
    assert "b|a" in out and run.session.e["errorterm_1"] == "b|a"
    W1 = None
    for line in out.splitlines():
        if line.strip().startswith("a |W"):
            W1 = float(line.split("W")[1].split()[0])
    out2 = run("manovatest a / b|a")
    W2 = float(next(ln for ln in out2.splitlines() if ln.strip().startswith("a |W")).split("W")[1].split()[0])
    assert W1 == W2
    assert "errorterm_2" not in run.session.e


def test_manova_error_terms_follow_the_slash_rules(run):
    # T1 / T2 T3 /: T2 is the error of T1 and is not tested; T3 is tested against the residual
    run(NEST)
    out = run("manova u v = c / b|c a /")
    rows = [ln.split("|")[0].strip() for ln in out.splitlines() if "|W" in ln]
    assert rows == ["Model", "c", "a"]


def test_manovatest_test_matrix(run):
    run(NEST + "quietly manova u v = b")
    out = run("matrix C = (1,-1,0,0 \\ 1,0,-1,0)\nmanovatest, test(C)")
    assert "Test constraints" in out and "(1)    1.b - 2.b = 0" in out
    st1 = _m(run, "stat")
    run("quietly manovatest b")
    assert np.allclose(st1, _m(run, "stat"))
    run("matrix C = (1,-1,0)\nmanovatest, test(C)")
    assert run.rc == 503


def test_anova_keeps_multilevel_nesting_names(run):
    run(NEST)
    run("quietly manova u v = a / b|c|a /")
    assert run.session.e["term_2"] == "b|c|a"
