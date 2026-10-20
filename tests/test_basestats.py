"""Phase 17: base statistics (resampling, ml, model building, postestimation, tests, agreement,
transformations and displays). Reference values come from identities and from independent
computations with numpy and scipy."""

import math

import numpy as np
from scipy import stats as st

DATA = """
clear
set seed 17
set obs 80
gen x = invnormal(runiform())
gen z = invnormal(runiform())
gen g = 1 + mod(_n, 4)
gen b1 = runiform() < invlogit(.5*x)
gen b2 = runiform() < invlogit(.5*x + .3*z)
gen y = 1 + 2*x - z + invnormal(runiform())
gen yp = exp(.8 + .5*x + .3*invnormal(runiform()))
gen byte f = runiform() < .4
"""


def _arr(run, name, cond=""):
    s = run.session
    v = s.data.get(name).data.astype(float)
    return v


def test_bitest_matches_scipy(run):
    run(DATA + "bitest b1 == .5")
    k, n = int(run.session.r["k"]), int(run.session.r["N"])
    p = st.binomtest(k, n, .5).pvalue
    assert abs(run.session.r["p"] - p) < 1e-10
    assert abs(run.session.r["p_u"] - st.binom.sf(k - 1, n, .5)) < 1e-12
    run("bitesti 20 12 .3")
    assert abs(run.session.r["p"] - st.binomtest(12, 20, .3).pvalue) < 1e-10


def test_ksmirnov_two_sample(run):
    out = run(DATA + "ksmirnov x, by(f)")
    x = _arr(run, "x")
    f = _arr(run, "f")
    res = st.ks_2samp(x[f == 0], x[f == 1], method="asymp")
    ne = (f == 0).sum() * (f == 1).sum() / ((f == 0).sum() + (f == 1).sum())
    assert abs(run.session.r["D"] - res.statistic) < 1e-12
    assert abs(run.session.r["p"] - st.kstwobign.sf(math.sqrt(ne) * res.statistic)) < 1e-10
    assert "Combined K-S" in out


def test_ksmirnov_one_sample(run):
    run(DATA + "ksmirnov x = normal(x)")
    x = np.sort(_arr(run, "x"))
    res = st.kstest(x, "norm")
    assert abs(run.session.r["D"] - res.statistic) < 1e-12


def test_median_test(run):
    out = run(DATA + "median x, by(f)")
    x, f = _arr(run, "x"), _arr(run, "f")
    med = np.median(x)
    tab = np.array([[((x <= med) & (f == k)).sum() for k in (0, 1)],
                    [((x > med) & (f == k)).sum() for k in (0, 1)]])
    chi2 = st.chi2_contingency(tab, correction=False)[0]
    assert abs(run.session.r["chi2"] - chi2) < 1e-10
    assert "Median test" in out


def test_nptrend_cuzick_equals_spearman_like_trend(run):
    run(DATA + "nptrend y, by(g)")
    y, g = _arr(run, "y"), _arr(run, "g")
    N = len(y)
    r = st.rankdata(y)
    L = (g * r).sum()
    EL = (N + 1) / 2 * g.sum()
    var = ((r - (N + 1) / 2) ** 2).sum() / (N - 1) * (N * (g ** 2).sum() - g.sum() ** 2) / N
    assert abs(run.session.r["z"] - (L - EL) / math.sqrt(var)) < 1e-10


def test_runtest_counts_runs(run):
    run("clear\ninput v\n1\n2\n8\n9\n3\n4\n7\n9\nend\nruntest v, threshold(5)")
    r = run.session.r
    assert r["N_below"] == 4 and r["N_above"] == 4
    assert r["n_runs"] == 4            # below below above above below below above above
    assert abs(r["mu"] - 5.0) < 1e-12


def test_ktau_matches_scipy(run):
    run(DATA + "ktau x y")
    x, y = _arr(run, "x"), _arr(run, "y")
    res = st.kendalltau(x, y)
    assert abs(run.session.r["tau_b"] - res.statistic) < 1e-10
    n = len(x)
    assert abs(run.session.r["tau_a"] - run.session.r["score"] / (n * (n - 1) / 2)) < 1e-12


def test_pcorr_equals_residual_correlation(run):
    run(DATA + "pcorr y x z")
    y, x, z = _arr(run, "y"), _arr(run, "x"), _arr(run, "z")

    def resid(a, others):
        Z = np.column_stack([np.ones(len(a))] + others)
        return a - Z @ np.linalg.lstsq(Z, a, rcond=None)[0]
    ry, rx = resid(y, [z]), resid(x, [z])
    pc = run.session.r["p_corr"].data[0, 0]
    assert abs(pc - np.corrcoef(ry, rx)[0, 1]) < 1e-10
    rx2 = resid(x, [z])
    sc = run.session.r["sp_corr"].data[0, 0]
    assert abs(sc - np.corrcoef(y, rx2)[0, 1]) < 1e-10


def test_tetrachoric_recovers_rho():
    from opendta.core import basestats as B
    rho, h, k = .6, .2, -.3
    p00 = B.bvn_cdf(h, k, rho)
    p0_, p_0 = st.norm.cdf(h), st.norm.cdf(k)
    P = np.array([[p00, p0_ - p00], [p_0 - p00, 1 - p0_ - p_0 + p00]])
    res = B.tetrachoric(P * 1e6)
    assert abs(res["rho"] - rho) < 1e-3
    # bivariate normal probability against scipy
    ref = st.multivariate_normal.cdf([h, k], mean=[0, 0], cov=[[1, rho], [rho, 1]])
    assert abs(p00 - ref) < 1e-6


def test_tetrachoric_command(run):
    out = run(DATA + "tetrachoric b1 b2")
    assert "Tetrachoric" in out
    assert -1 < run.session.r["rho"] < 1


def test_ameans(run):
    run(DATA + "ameans yp")
    yp = _arr(run, "yp")
    r = run.session.r
    assert abs(r["mean"] - yp.mean()) < 1e-10
    assert abs(r["mean_g"] - math.exp(np.log(yp).mean())) < 1e-10
    assert abs(r["mean_h"] - 1 / (1 / yp).mean()) < 1e-10


SHROUT = """
clear
input r1 r2 r3 r4
9 2 5 8
6 1 3 2
8 4 6 8
7 1 2 6
10 5 6 9
6 2 4 7
end
gen target = _n
reshape long r, i(target) j(judge)
"""


def test_icc_shrout_fleiss(run):
    run(SHROUT + "icc r target judge, absolute")
    assert round(run.session.r["icc_i"], 2) == 0.29 and round(run.session.r["icc_avg"], 2) == 0.62
    run("icc r target judge")
    assert round(run.session.r["icc_i"], 2) == 0.71 and round(run.session.r["icc_avg"], 2) == 0.91
    run("icc r target")
    assert round(run.session.r["icc_i"], 2) == 0.17 and round(run.session.r["icc_avg"], 2) == 0.44
    lo, hi = run.session.r["icc_i_lb"], run.session.r["icc_i_ub"]
    assert lo < run.session.r["icc_i"] < hi


def test_loneway_matches_oneway(run):
    out = run(DATA + "loneway y g\noneway y g")
    F = st.f_oneway(*[_arr(run, "y")[_arr(run, "g") == k] for k in (1, 2, 3, 4)]).statistic
    out2 = run("loneway y g")
    assert abs(run.session.r["F"] - F) < 1e-10
    assert run.session.r["lb"] < run.session.r["rho"] < run.session.r["ub"]
    assert "Intraclass" in out2


def test_esize_cohen_hedges(run):
    run(DATA + "esize twosample y, by(f) all")
    y, f = _arr(run, "y"), _arr(run, "f")
    a, b = y[f == 0], y[f == 1]
    sp = math.sqrt(((len(a) - 1) * a.var(ddof=1) + (len(b) - 1) * b.var(ddof=1)) / (len(a) + len(b) - 2))
    d = (a.mean() - b.mean()) / sp
    r = run.session.r
    assert abs(r["d"] - d) < 1e-10
    df = len(a) + len(b) - 2
    J = math.gamma(df / 2) / (math.sqrt(df / 2) * math.gamma((df - 1) / 2))
    assert abs(r["g"] - J * d) < 1e-10
    assert r["lb_d"] < r["d"] < r["ub_d"]
    run("esizei 30 5 2 30 6 2")
    assert abs(run.session.r["d"] - (-0.5)) < 1e-12


def test_lnskew0_zero_skewness(run):
    run(DATA + "lnskew0 ly = yp + 0")
    from opendta.core import basestats as B
    v = _arr(run, "ly")
    assert abs(B.skewness(v)) < 1e-6
    assert abs(run.session.r["skewness"]) < 1e-6


def test_ladder_identity_matches_sktest(run):
    out = run(DATA + "ladder yp\nsktest yp")
    out = run("ladder yp")
    # the identity row: adjusted chi2 as in sktest
    row = [ln for ln in out.splitlines() if ln.startswith("identity")][0]
    chi2 = float(row.split()[2])
    run("sktest yp")
    assert abs(chi2 - run.session.r["chi2"]) < 0.01
    assert "cubic" in out and "1/cubic" in out


def test_boxcox_matches_scipy(run):
    run(DATA + "boxcox yp, nolog")
    yp = _arr(run, "yp")
    lam = st.boxcox(yp)[1]
    assert abs(run.session.e["theta"] - lam) < 1e-3
    assert run.session.e["cmd"] == "boxcox"


def test_cumul(run):
    run("clear\ninput v\n3\n1\n2\n2\nend\ncumul v, gen(c)\ncumul v, gen(c2) equal")
    assert list(_arr(run, "c")) == [1.0, .25, .75, .75]
    assert list(_arr(run, "c2")) == [1.0, .25, .5, .75]


def test_smooth_median_of_three(run):
    run("clear\ninput v\n1\n5\n2\n6\n3\nend\nsmooth 3 v, gen(s)")
    assert list(_arr(run, "s")) == [1, 2, 5, 3, 3]
    run("smooth 3, twice v, gen(s2)")
    assert run.session.data.has("s2")


def test_stem_and_lv(run):
    out = run("clear\ninput v\n10\n11\n12\n25\n26\n40\nend\nstem v")
    assert "Stem-and-leaf plot for v" in out
    out = run("clear\nset obs 9\ngen v = _n\nlv v")
    assert "M" in out and "inner fence" in out


def test_jackknife_mean_se(run):
    out = run(DATA + "jackknife m=r(mean), double: summarize y")
    n = 80
    y = _arr(run, "y")
    se = run.session.e["V"].data[0, 0] ** .5
    assert abs(se - y.std(ddof=1) / math.sqrt(n)) < 1e-10
    assert "Jackknife" in out


def test_permute_and_simulate(run):
    run(DATA + "set seed 3\npermute f d=(r(mu_2)-r(mu_1)), reps(30) nodots: ttest y, by(f)")
    run("program define sim1, rclass\n drop _all\n set obs 20\n gen u = runiform()\n summarize u\n return scalar m = r(mean)\nend\n"
        "simulate m=r(m), reps(15) nodots: sim1")
    assert run.session.data.nobs == 15 and run.session.data.has("m")


def test_ml_lf_matches_logit(run):
    run(DATA + "logit b1 x z")
    b = run.session.e["b"].data.copy()
    run("program define mylogit\n args lnf xb\n quietly replace `lnf' = ln(invlogit(`xb')) if $ML_y1 == 1\n"
        " quietly replace `lnf' = ln(invlogit(-`xb')) if $ML_y1 == 0\nend\n"
        "ml model lf mylogit (b1 = x z)\nml maximize, nolog")
    b2 = run.session.e["b"].data
    assert np.allclose(b, b2, atol=1e-5)
    run("mlexp (b1*ln(invlogit({xb: x z _cons})) + (1-b1)*ln(invlogit(-{xb:})))")
    assert np.allclose(b, run.session.e["b"].data, atol=1e-5)


def test_ml_d0_and_d1(run):
    run(DATA + "logit b1 x")
    b = run.session.e["b"].data.copy()
    run("program define d0log\n args todo b lnfj g H\n tempvar xb\n mleval `xb' = `b'\n"
        " mlsum `lnfj' = cond($ML_y1==1, ln(invlogit(`xb')), ln(invlogit(-`xb')))\nend\n"
        "ml model d0 d0log (b1 = x)\nml maximize, nolog")
    assert np.allclose(b, run.session.e["b"].data, atol=1e-4)


def test_nestreg_final_model_and_block_test(run):
    out = run(DATA + "nestreg: regress y (x) (z)")
    assert "Block" in out
    r2 = run.session.e["r2"]
    run("regress y x z")
    assert abs(run.session.e["r2"] - r2) < 1e-12


def test_stepwise_keeps_signal(run):
    run(DATA + "gen noise = invnormal(runiform())\nstepwise, pr(.05): regress y x z noise")
    assert run.session.e["cmd"] == "regress"
    names = list(run.session.e["b"].colnames)
    assert "x" in names and "z" in names and "noise" not in names


def test_mkspline_linear_pieces_add_up(run):
    run(DATA + "mkspline s1 -.5 s2 .5 s3 = x")
    x = _arr(run, "x")
    tot = _arr(run, "s1") + _arr(run, "s2") + _arr(run, "s3")
    assert np.allclose(tot, x, atol=1e-6)
    run("mkspline c = x, cubic nknots(4)")
    assert run.session.data.has("c1") and run.session.data.has("c3")


def test_xi_creates_dummies(run):
    run(DATA + "xi: regress y i.g x")
    names = run.session.data.names
    assert sum(1 for n in names if n.startswith("_Ig_")) == 3


def test_fp_and_mfp_run(run):
    out = run(DATA + "gen w = exp(x)\nfp <w>: regress y <w>")
    assert "Deviance" in out or "deviance" in out
    out = run("mfp: regress y w z")
    assert run.session.e["cmd"] == "regress"


def test_testnl_linear_equals_wald(run):
    run(DATA + "regress y x z\ntestnl _b[x] = 2")
    b = run.session.e["b"].data[0, 0]
    se = run.session.e["V"].data[0, 0] ** .5
    assert abs(run.session.r["F"] - ((b - 2) / se) ** 2) < 1e-8


def test_predictnl_delta_method(run):
    run(DATA + "regress y x z\npredictnl h = 2*_b[x], se(seh)")
    se = run.session.e["V"].data[0, 0] ** .5
    assert abs(_arr(run, "seh")[0] - 2 * se) < 1e-6


def test_linktest_adds_hatsq(run):
    out = run(DATA + "regress y x z\nlinktest")
    assert "_hatsq" in out


def test_contrast_and_pwcompare(run):
    out = run(DATA + "regress y i.g\ncontrast g")
    assert "g" in out
    run("regress y i.g\ntest 2.g 3.g 4.g")
    F = run.session.r["F"]
    run("regress y i.g\ncontrast g")
    assert abs(run.session.r["F"].data[0, 0] - F) < 1e-8
    out = run("pwcompare g, mcompare(bonferroni)")
    assert "Bonferroni" in out or "bonferroni" in out


def test_pwmean_and_suest(run):
    out = run(DATA + "pwmean y, over(g) mcompare(noadjust)")
    y, g = _arr(run, "y"), _arr(run, "g")
    assert "Pairwise" in out or "pairwise" in out
    run("logit b1 x\nestimates store a\nlogit b2 x\nestimates store b\nsuest a b")
    assert run.session.e["cmd"] == "suest"


def test_margins_after_tobit_and_ologit(run):
    run(DATA + "gen yc = max(y, 0)\ntobit yc x z, ll(0)\nmargins, dydx(x)")
    assert run.rc == 0
    run("gen o = 1 + (y > 0) + (y > 2)\nologit o x z\nmargins, predict(outcome(2))")
    assert run.rc == 0


def test_graph_commands_build_graphs(run):
    from opendta.graphics.store import store
    run(DATA + "gen se = .3\nquantile x")
    g = store(run.session).graphs["Graph"]
    assert len(g.panels[0].layers[0].x) == 80
    run("symplot x")
    run("serrbar y se x")
    run("dotplot x")
    run("sunflower y x")
    run("cusum y x, gen(cs)")
    assert run.session.data.has("cs")
    run("gladder yp")
    g = store(run.session).graphs["Graph"]
    assert g.family == "combine" and len(g.children) == 9


def test_lpoly_recovers_a_line(run):
    run("clear\nset obs 40\ngen x = _n\ngen y = 2 + 3*x\nlpoly y x, degree(1) bwidth(10) nograph generate(gx gy) n(10)")
    gx, gy = _arr(run, "gx")[:10], _arr(run, "gy")[:10]
    assert np.allclose(gy, 2 + 3 * gx, atol=1e-8)
