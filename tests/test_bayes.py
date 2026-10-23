"""bayesmh, the bayes prefix and the posterior summaries, checked against conjugate posteriors,
maximum likelihood and the exact marginal likelihood."""

import math

import numpy as np
import pytest

from opendta.core import mcmc


def _load(run, cols):
    n = len(next(iter(cols.values())))
    run("clear\ninput " + " ".join("double " + k for k in cols) + "\n"
        + "\n".join(" ".join(repr(float(cols[k][i])) for k in cols) for i in range(n)) + "\nend")


@pytest.fixture
def data(run):
    rng = np.random.default_rng(3)
    N = 60
    x = rng.normal(size=N)
    z = rng.normal(size=N)
    y = 2 + 1.5 * x - .7 * z + rng.normal(scale=1.3, size=N)
    b = (rng.random(N) < 1 / (1 + np.exp(-(.4 + .9 * x)))).astype(float)
    run.cols = dict(x=x, z=z, y=y, b=b)
    _load(run, run.cols)
    return run


def _mat(run, name):
    return run.session.e[name]


def test_ess_of_an_ar1_chain():
    rng = np.random.default_rng(0)
    phi, n = 0.8, 20000
    e = rng.normal(size=n)
    x = np.zeros(n)
    for i in range(1, n):
        x[i] = phi * x[i - 1] + e[i]
    theory = n * (1 - phi) / (1 + phi)
    assert mcmc.ess(x) == pytest.approx(theory, rel=0.3)
    assert mcmc.ess(rng.normal(size=5000)) == pytest.approx(5000, rel=0.15)


def test_intervals_of_a_normal_sample():
    rng = np.random.default_rng(1)
    x = rng.normal(size=100000)
    lo, hi = mcmc.equal_tailed(x, 95)
    assert (lo, hi) == pytest.approx((-1.96, 1.96), abs=0.03)
    lo, hi = mcmc.hpd(x, 95)
    assert (lo, hi) == pytest.approx((-1.96, 1.96), abs=0.05)


def test_conjugate_normal_mean_and_marginal_likelihood(data):
    rng = np.random.default_rng(7)
    y = rng.normal(1.0, 1.0, size=25)
    _load(data, {"w": y})
    data("quietly bayesmh w = ({mu}), likelihood(normal(1)) prior({mu}, normal(0,4)) "
         "rseed(11) mcmcsize(20000) nomodelsummary")
    n = len(y)
    v = 1 / (n + 1 / 4)
    m = v * y.sum()
    mean = float(_mat(data, "mean").data[0, 0])
    sd = float(_mat(data, "sd").data[0, 0])
    assert mean == pytest.approx(m, abs=0.02)
    assert sd == pytest.approx(math.sqrt(v), rel=0.05)
    # y ~ N(0, I + 4 * 11'): exact log marginal likelihood
    S = np.eye(n) + 4 * np.ones((n, n))
    exact = -0.5 * (n * math.log(2 * math.pi) + np.linalg.slogdet(S)[1] + y @ np.linalg.solve(S, y))
    assert float(data.session.e["lml_lm"]) == pytest.approx(exact, abs=0.1)


def test_linear_regression_matches_ols(data):
    data("quietly bayesmh y x z, likelihood(normal({s2})) prior({y:}, normal(0,10000)) "
         "prior({s2}, igamma(.01,.01)) rseed(5) mcmcsize(10000)")
    mean = _mat(data, "mean").data[0].copy()
    assert float(data.session.e["N"]) == 60
    data("quietly regress y x z")
    ols = data.session.e["b"].data[0]
    assert mean[:3] == pytest.approx(ols[[0, 1, 2]], abs=0.1)


def test_reproducible_with_rseed_and_replay(data):
    cmd = ("bayesmh y x, likelihood(normal({s2})) prior({y:}, normal(0,100)) "
           "prior({s2}, igamma(.01,.01)) rseed(9) mcmcsize(2000) nomodelsummary")
    a = data("quietly " + cmd)
    m1 = _mat(data, "mean").data.copy()
    data("quietly " + cmd)
    assert np.allclose(m1, _mat(data, "mean").data)
    out = data("bayesmh")
    assert "Mean   Std. Dev.     MCSE     Median" in out and "Equal-tailed" in out
    out = data("bayesmh, hpd clevel(90)")
    assert "HPD" in out and "90% Cred. Interval" in out


def test_prefix_logit_close_to_ml(data):
    out = data("bayes, rseed(2) mcmcsize(10000) nomodelsummary: logit b x")
    assert "Bayesian logistic regression" in out
    data("quietly logit b x")
    ml = data.session.e["b"].data[0]
    data("quietly bayes, rseed(2) mcmcsize(10000): logit b x")
    mean = _mat(data, "mean").data[0]
    assert mean == pytest.approx(ml, abs=0.2)


def test_bayesstats_functions_of_parameters(data):
    data("quietly bayesmh y x, likelihood(normal({s2})) prior({y:}, normal(0,100)) "
         "prior({s2}, igamma(.01,.01)) rseed(3) mcmcsize(3000)")
    out = data("bayesstats summary {y:x} (diff: {y:x}-{y:_cons})")
    assert "diff" in out and "MCMC sample size =     3,000" in out
    out = data("bayesstats ess")
    assert "ESS" in out and "Corr. time" in out
    out = data("bayestest interval {y:x}, lower(0)")
    assert "Interval tests" in out


def test_model_comparison_with_stored_estimates(data, tmp_path):
    base = ("likelihood(normal({s2})) prior({y:}, normal(0,100)) prior({s2}, igamma(.01,.01)) "
            "rseed(3) mcmcsize(3000)")
    # as in Stata, estimates store needs the draws saved with saving()
    assert "you must save simulation results" in data(f"capture noisily bayesmh y x z, {base}\n"
                                                      "capture noisily estimates store full")
    f = str(tmp_path / "full").replace("\\", "/")
    data(f"quietly bayesmh y x z, {base} saving({f}, replace)")
    data("estimates store full")
    f = str(tmp_path / "small").replace("\\", "/")
    data(f"quietly bayesmh y z, {base} saving({f}, replace)")
    data("estimates store small")
    out = data("bayesstats ic full small")
    assert "full" in out and "small" in out and "DIC" in out
    out = data("bayestest model full small")
    probs = [float(ln.split()[-1]) for ln in out.splitlines() if ln.strip().startswith(("full", "small"))]
    assert sum(probs) == pytest.approx(1.0, abs=1e-3)
    assert probs[0] > 0.9                                # x matters


def test_saving_draws_and_graphs_keep_the_data(data, tmp_path):
    f = str(tmp_path / "draws").replace("\\", "/")
    data("quietly bayesmh y x, likelihood(normal({s2})) prior({y:}, normal(0,100)) "
         f"prior({{s2}}, igamma(.01,.01)) rseed(3) mcmcsize(500) saving({f}, replace)")
    data(f'describe using "{f}.dta"')
    assert "obs:           500" in data("describe using \"" + f + ".dta\"")
    n, k = data.session.data.nobs, len(data.session.data.vars)
    for g in ("trace {s2}", "histogram {y:x}", "ac {s2}", "kdensity {s2}", "cusum {s2}"):
        data(f"bayesgraph {g}")
        assert data.rc == 0, g
    assert data.session.data.nobs == n and len(data.session.data.vars) == k


def test_errors(data):
    data("bayesmh y x, likelihood(normal({s2})) prior({y:}, normal(0,100))")
    assert data.rc == 198                                 # no prior for {s2}
    data("bayesmh y x, prior({y:}, normal(0,100))")
    assert data.rc == 198                                 # no likelihood
    data("bayes: ivregress 2sls y x")
    assert data.rc == 199
    data("bayesstats summary")
    assert data.rc in (301, 198)
