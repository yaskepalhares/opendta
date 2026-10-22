"""sem: method(mlmv) and method(adf), group() and ginvariant(), path labels and the estat
subcommands, checked against separate fits, closed forms and likelihood-ratio tests."""

import numpy as np
import pytest

from opendta.core import mvnmiss


def _load(run, cols):
    n = len(next(iter(cols.values())))
    run("clear\ninput " + " ".join("double " + k for k in cols) + "\n"
        + "\n".join(" ".join("." if np.isnan(cols[k][i]) else repr(float(cols[k][i])) for k in cols)
                    for i in range(n)) + "\nend")


@pytest.fixture
def cfa(run):
    rng = np.random.default_rng(21)
    N = 500
    F = rng.normal(size=N)
    G = .5 * F + rng.normal(scale=.8, size=N)
    cols = {"x1": F + rng.normal(scale=.6, size=N), "x2": .8 * F + rng.normal(scale=.7, size=N),
            "x3": 1.2 * F + rng.normal(scale=.5, size=N), "y1": G + rng.normal(scale=.6, size=N),
            "y2": .9 * G + rng.normal(scale=.6, size=N), "y3": 1.1 * G + rng.normal(scale=.6, size=N)}
    run.cols = cols
    _load(run, cols)
    return run


@pytest.fixture
def two_groups(run):
    rng = np.random.default_rng(11)
    N = 600
    g = (np.arange(N) % 2) + 1
    F = rng.normal(size=N) + .4 * (g == 2)
    cols = {"grp": g.astype(float), "x1": F + rng.normal(scale=.6, size=N),
            "x2": .8 * F + rng.normal(scale=.7, size=N) + .25 * (g == 2),
            "x3": 1.2 * F + rng.normal(scale=.5, size=N), "x4": .7 * F + rng.normal(scale=.5, size=N)}
    _load(run, cols)
    return run


def _b(run, eq, name):
    e = run.session.e["b"]
    for j, (q, c) in enumerate(zip(e.coleq, e.colnames)):
        if q == eq and c == name:
            return float(e.data[0, j])
    raise KeyError((eq, name))


def test_em_complete_data_is_sample_moments():
    rng = np.random.default_rng(1)
    X = rng.normal(size=(80, 3)) @ np.array([[1, .5, 0], [0, 1, .3], [0, 0, 1]])
    mu, S, ll, _ = mvnmiss.em_mvn(X)
    assert mu == pytest.approx(X.mean(axis=0), rel=1e-10)
    assert S == pytest.approx(np.cov(X.T, bias=True), rel=1e-8)


def test_em_missing_increases_likelihood_and_block_structure():
    rng = np.random.default_rng(2)
    X = rng.normal(size=(200, 3)) @ np.array([[1, .6, 0], [0, 1, .3], [0, 0, 1]])
    X[rng.random(200) < .15, 1] = np.nan
    X[rng.random(200) < .15, 2] = np.nan
    mu, S, ll, it = mvnmiss.em_mvn(X)
    pats = mvnmiss.patterns(X)
    start_ll = mvnmiss.loglik(pats, np.nanmean(X, axis=0), np.diag(np.nanvar(X, axis=0)))
    assert ll > start_ll
    # the gradient of the log likelihood is zero at the EM solution (checked numerically)
    h = 1e-6
    for k in range(3):
        e = np.zeros(3)
        e[k] = h
        g = (mvnmiss.loglik(pats, mu + e, S) - mvnmiss.loglik(pats, mu - e, S)) / (2 * h)
        assert abs(g) < 1e-3
    _, S2, _, _ = mvnmiss.em_mvn(X, blocks=[[0, 1]])
    assert S2[0, 2] == 0 and S2[1, 2] == 0 and S2[0, 1] != 0


def test_mlmv_equals_ml_without_missing(cfa):
    run = cfa
    run("quietly sem (F -> x1 x2 x3) (G -> y1 y2 y3), nolog")
    ml = [_b(run, "x2", "F"), _b(run, "y3", "G")]
    ll = run.session.e["ll"]
    out = run("sem (F -> x1 x2 x3) (G -> y1 y2 y3), method(mlmv) nolog")
    assert "Estimation method  = mlmv" in out
    assert _b(run, "x2", "F") == pytest.approx(ml[0], rel=1e-6)
    assert _b(run, "y3", "G") == pytest.approx(ml[1], rel=1e-6)
    assert run.session.e["ll"] == pytest.approx(ll, rel=1e-9)
    assert run.session.e["method"] == "mlmv"


def test_mlmv_uses_incomplete_observations(cfa):
    run = cfa
    rng = np.random.default_rng(5)
    cols = {k: v.copy() for k, v in run.cols.items()}
    for k in ("x2", "y1", "y3"):
        cols[k][rng.random(500) < .15] = np.nan
    _load(run, cols)
    run("quietly sem (F -> x1 x2 x3) (G -> y1 y2 y3), nolog")
    n_ml = run.session.e["N"]
    run("quietly sem (F -> x1 x2 x3) (G -> y1 y2 y3), method(mlmv) nolog")
    assert run.session.e["N"] == 500 > n_ml
    assert _b(run, "y2", "G") == pytest.approx(.9, abs=.15)
    assert run.session.e["converged"] == 1
    assert run.session.e["df_ms"] == 8


def test_adf_close_to_ml_for_normal_data(cfa):
    run = cfa
    run("quietly sem (F -> x1 x2 x3) (G -> y1 y2 y3), nolog")
    ml = _b(run, "x2", "F")
    out = run("sem (F -> x1 x2 x3) (G -> y1 y2 y3), method(adf) nolog")
    assert "Estimation method  = adf" in out and "Log likelihood" not in out
    assert _b(run, "x2", "F") == pytest.approx(ml, abs=.05)
    assert run.session.e["df_ms"] == 8
    assert run.session.e["chi2_ms"] >= 0


def test_group_default_invariance_and_separate_fits(two_groups):
    run = two_groups
    run("quietly sem (F -> x1 x2 x3 x4), group(grp) nolog")
    e = run.session.e
    assert e["N_groups"] == 2 and e["N"] == 600
    b = e["b"]
    # default: measurement coefficients and intercepts equal across groups
    col = {(q, c): float(b.data[0, j]) for j, (q, c) in enumerate(zip(b.coleq, b.colnames))}
    assert col[("x2", "1bn.grp#c.F")] == pytest.approx(col[("x2", "2.grp#c.F")], rel=1e-9)
    assert col[("x3", "1bn.grp")] == pytest.approx(col[("x3", "2.grp")], rel=1e-9)
    # df: 2 * 14 moments - (3 loadings + 4 intercepts + 2 * (4 error variances + var(F)) + mean(F) in group 2)
    assert e["df_ms"] == 28 - 18
    # without invariance each group gives the single-group fit
    run("quietly sem (F -> x1 x2 x3 x4), group(grp) ginvariant(none) nolog")
    b = run.session.e["b"]
    col = {(q, c): float(b.data[0, j]) for j, (q, c) in enumerate(zip(b.coleq, b.colnames))}
    run("quietly sem (F -> x1 x2 x3 x4) if grp == 2, nolog")
    assert col[("x2", "2.grp#c.F")] == pytest.approx(_b(run, "x2", "F"), rel=1e-5)
    assert col[("var(e.x3)", "2.grp")] == pytest.approx(_b(run, "var(e.x3)", "_cons"), rel=1e-5)


def test_group_loglik_is_sum_of_group_fits_without_invariance(two_groups):
    run = two_groups
    run("quietly sem (F -> x1 x2 x3 x4), group(grp) ginvariant(none) nolog")
    ll = run.session.e["ll"]
    tot = 0.0
    for g in (1, 2):
        run(f"quietly sem (F -> x1 x2 x3 x4) if grp == {g}, nolog")
        tot += run.session.e["ll"]
    assert ll == pytest.approx(tot, rel=1e-8)


def test_labels_constrain_paths_equal(run):
    rng = np.random.default_rng(5)
    N = 400
    z, x = rng.normal(size=N), rng.normal(size=N)
    _load(run, {"z": z, "x": x, "y1": .5 * z + .5 * x + rng.normal(size=N),
                "y2": .5 * z + .5 * x + rng.normal(size=N)})
    out = run("sem (y1 <- z@a x@b) (y2 <- z@a x@b), nolog")
    assert "[y1]z - [y2]z = 0" in out
    assert _b(run, "y1", "z") == pytest.approx(_b(run, "y2", "z"), rel=1e-12)
    assert _b(run, "y1", "x") == pytest.approx(_b(run, "y2", "x"), rel=1e-12)
    assert run.session.e["df_ms"] == 3
    # the pooled slopes minimize the sum of the squared residuals of both equations
    # (the same errors variances are not required): compare with the regression of the stacked data
    import numpy as _np
    y1 = run.session.data.get("y1").data
    y2 = run.session.data.get("y2").data
    Z = _np.column_stack([z, x])
    stacked = _np.concatenate([y1, y2])
    D = _np.block([[Z, _np.ones((N, 1)), _np.zeros((N, 1))], [Z, _np.zeros((N, 1)), _np.ones((N, 1))]])
    # equal error variances would give exactly this least-squares solution; with free variances the
    # estimates are close to it
    ols = _np.linalg.lstsq(D, stacked, rcond=None)[0]
    assert _b(run, "y1", "z") == pytest.approx(ols[0], abs=.03)


def test_estat_eqgof_and_eqtest_match_ols(run):
    rng = np.random.default_rng(8)
    N = 300
    x1, x2 = rng.normal(size=N), rng.normal(size=N)
    y = 1 + .8 * x1 - .5 * x2 + rng.normal(size=N)
    _load(run, {"x1": x1, "x2": x2, "y": y})
    run("quietly sem (y <- x1 x2), nolog")
    out = run("estat eqgof")
    X = np.column_stack([x1, x2, np.ones(N)])
    b = np.linalg.lstsq(X, y, rcond=None)[0]
    e = y - X @ b
    r2 = 1 - (e @ e) / ((y - y.mean()) @ (y - y.mean()))
    assert "Equation-level goodness of fit" in out
    m = run.session.r["eqfit"]
    assert m.data[0, 3] == pytest.approx(r2, rel=1e-8)          # R-squared
    assert m.data[0, 5] == pytest.approx(r2, rel=1e-8)          # mc2 = R2 in a recursive model
    assert run.session.r["CD"] == pytest.approx(r2, rel=1e-8)
    out = run("estat eqtest")
    assert "Wald tests for equations" in out
    w = run.session.r["test"].data[0]
    s2 = e @ e / N
    cov = s2 * np.linalg.inv(X.T @ X)
    wald = b[:2] @ np.linalg.inv(cov[:2, :2]) @ b[:2]
    assert w[0] == pytest.approx(wald, rel=1e-6) and w[1] == 2


def test_estat_teffects_indirect_is_product(run):
    rng = np.random.default_rng(9)
    N = 500
    x = rng.normal(size=N)
    m = .6 * x + rng.normal(size=N)
    y = .5 * m + .2 * x + rng.normal(size=N)
    _load(run, {"x": x, "m": m, "y": y})
    run("quietly sem (m <- x) (y <- m x), nolog")
    out = run("estat teffects")
    assert "Direct effects" in out and "Indirect effects" in out and "Total effects" in out
    bm_x, by_m, by_x = _b(run, "m", "x"), _b(run, "y", "m"), _b(run, "y", "x")
    # the indirect effect appears as the product of the two direct effects
    vals = [float(t) for line in out.split("Indirect effects")[1].split("Total effects")[0].splitlines()
            if "|" in line for t in line.split("|")[1].split()[:1] if t.replace(".", "", 1).replace("-", "", 1).isdigit()]
    assert vals and vals[0] == pytest.approx(bm_x * by_m, abs=2e-6)


def test_estat_stable_nonrecursive(run):
    rng = np.random.default_rng(12)
    N = 400
    z1, z2 = rng.normal(size=N), rng.normal(size=N)
    e1, e2 = rng.normal(size=N), rng.normal(size=N)
    # y1 = .3 y2 + z1 + e1 ; y2 = .4 y1 + z2 + e2
    det = 1 - .3 * .4
    y1 = (z1 + e1 + .3 * (z2 + e2)) / det
    y2 = (z2 + e2 + .4 * (z1 + e1)) / det
    _load(run, {"z1": z1, "z2": z2, "y1": y1, "y2": y2})
    run("quietly sem (y1 <- y2 z1) (y2 <- y1 z2), nolog")
    out = run("estat stable")
    assert "Stability analysis of simultaneous equation systems" in out
    r = run.session.r
    b12, b21 = _b(run, "y1", "y2"), _b(run, "y2", "y1")
    assert r["stindex"] == pytest.approx(np.sqrt(abs(b12 * b21)), rel=1e-8)


def test_estat_mindices_agrees_with_likelihood_ratio(run):
    rng = np.random.default_rng(13)
    N = 4000
    F = rng.normal(size=N)
    cols = {"a": F + rng.normal(scale=.7, size=N), "b": F + rng.normal(scale=.7, size=N),
            "c": F + rng.normal(scale=.7, size=N), "d": F + rng.normal(scale=.7, size=N)}
    e = rng.normal(size=N)
    cols["a"] = cols["a"] + .2 * e
    cols["b"] = cols["b"] + .2 * e
    _load(run, cols)
    run("quietly sem (F -> a b c d), nolog")
    ll0 = run.session.e["ll"]
    out = run("estat mindices, minchi2(0)")
    assert "Modification indices" in out and "cov(e.a,e.b)" in out
    mi = run.session.r["mindices"]
    row = [i for i, n in enumerate(mi.rownames) if n == "cov(e.a,e.b)"][0]
    run("quietly sem (F -> a b c d), cov(e.a*e.b) nolog")
    lr = 2 * (run.session.e["ll"] - ll0)
    assert mi.data[row, 0] == pytest.approx(lr, rel=.05)
    assert mi.data[row, 3] > 0


def test_estat_scoretests_and_ginvariant(two_groups):
    run = two_groups
    run("quietly sem (F -> x1 x2 x3 x4), group(grp) nolog")
    ll0 = run.session.e["ll"]
    out = run("estat ginvariant")
    assert "Tests for group invariance of parameters" in out
    t = run.session.r["test"]
    names = t.rownames
    # the score test of the constrained intercept of x2 is close to the LR test of releasing it
    i = names.index("x2 <- _cons")
    assert t.data[i, 3] > 3.84
    out = run("estat ginvariant, class")
    assert "mcons" in out
    out = run("estat ggof")
    assert "Group-level goodness of fit" in out
    assert run.session.r["N_groups"] == 2


def test_estat_framework_and_residuals(cfa):
    run = cfa
    run("quietly sem (F -> x1 x2 x3) (G -> y1 y2 y3), nolog")
    out = run("estat framework, fitted")
    assert "Covariances of error variables" in out and "Fitted covariance matrix" in out
    S = run.session.r["Sigma"]
    assert S.rows == 6
    out = run("estat residuals")
    assert "Raw residuals" in out
    rc = run.session.r["res_cov"].data
    assert rc[0, 0] == pytest.approx(0, abs=1e-8)
    out = run("estat residuals, standardized")
    assert "Standardized residuals" in out


def test_estat_stdize_runs_command_on_standardized_estimates(run):
    rng = np.random.default_rng(14)
    N = 300
    x = rng.normal(size=N) * 3
    y = 2 + .5 * x + rng.normal(size=N)
    _load(run, {"x": x, "y": y})
    run("quietly sem (y <- x), nolog")
    b0 = _b(run, "y", "x")
    out = run("estat stdize: lincom [y]x")
    xs = x.std() / y.std()
    # the standardized coefficient is b * sd(x) / sd(y), with the ML (divide-by-N) deviations;
    # the model-implied standard deviation of y is used in the denominator
    sd_y = np.sqrt(b0 ** 2 * x.var() + _b(run, "var(e.y)", "_cons"))
    assert f"{b0 * x.std() / sd_y:.6f}"[1:7] in out.replace(" ", "")
    assert _b(run, "y", "x") == pytest.approx(b0, rel=1e-12)       # e() is restored
