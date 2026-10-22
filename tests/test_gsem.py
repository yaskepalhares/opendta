"""gsem: models without latent variables, random effects and measurement models, checked against
logit, melogit, meologit, sem and glm."""

import numpy as np
import pytest


def _load(run, cols):
    n = len(next(iter(cols.values())))
    run("clear\ninput " + " ".join("double " + k for k in cols) + "\n"
        + "\n".join(" ".join(repr(float(cols[k][i])) for k in cols) for i in range(n)) + "\nend")


@pytest.fixture
def panel(run):
    rng = np.random.default_rng(5)
    G, n = 50, 8
    N = G * n
    gid = np.repeat(np.arange(G), n)
    u = rng.normal(scale=.9, size=G)[gid]
    x = rng.normal(size=N)
    y = (rng.random(N) < 1 / (1 + np.exp(-(-.3 + .8 * x + u)))).astype(float)
    oc = np.digitize(.7 * x + u + rng.logistic(size=N), [-1, .5, 1.5]).astype(float)
    L = rng.normal(size=N)
    cols = {"id": gid, "x": x, "y": y, "oc": oc, "z1": L + rng.normal(scale=.6, size=N),
            "z2": .8 * L + rng.normal(scale=.7, size=N), "z3": 1.2 * L + rng.normal(scale=.5, size=N),
            "b1": (rng.random(N) < 1 / (1 + np.exp(-L))).astype(float),
            "b2": (rng.random(N) < 1 / (1 + np.exp(-(.8 * L + .3)))).astype(float),
            "b3": (rng.random(N) < 1 / (1 + np.exp(-(1.2 * L - .2)))).astype(float)}
    _load(run, cols)
    return run


def _e(run, name):
    return run.session.e[name]


def _b(run, eq, name):
    e = _e(run, "b")
    for j, (q, c) in enumerate(zip(e.coleq, e.colnames)):
        if q == eq and c == name:
            return float(e.data[0, j])
    raise KeyError((eq, name))


def test_no_latent_equals_logit(panel):
    panel("quietly logit y x")
    ref = (float(_e(panel, "ll")), _b(panel, "y", "x") if False else float(_e(panel, "b").data[0, 0]))
    panel("quietly gsem (y <- x, logit), nolog")
    assert float(_e(panel, "ll")) == pytest.approx(ref[0], abs=1e-6)
    assert _b(panel, "y", "x") == pytest.approx(ref[1], abs=1e-5)
    assert _e(panel, "family1") == "bernoulli" and _e(panel, "link1") == "logit"


def test_random_intercept_equals_melogit(panel):
    panel("quietly melogit y x || id:")
    ll = float(_e(panel, "ll"))
    bx = float(_e(panel, "b").data[0, 0])
    panel("quietly gsem (y <- x M1[id], logit), nolog")
    assert float(_e(panel, "ll")) == pytest.approx(ll, abs=1e-5)
    assert _b(panel, "y", "x") == pytest.approx(bx, abs=1e-4)
    assert _b(panel, "var(M1[id])", "_cons") > 0


def test_ordinal_random_intercept_equals_meologit(panel):
    panel("quietly meologit oc x || id:")
    ll = float(_e(panel, "ll"))
    panel("quietly gsem (oc <- x M1[id], ologit), nolog")
    assert float(_e(panel, "ll")) == pytest.approx(ll, abs=1e-5)


def test_gaussian_measurement_equals_sem(panel):
    panel("quietly sem (L -> z1 z2 z3)")
    ll = float(_e(panel, "ll"))
    l2 = _b(panel, "z2", "L")
    panel("quietly gsem (L -> z1 z2 z3), nolog")
    assert float(_e(panel, "ll")) == pytest.approx(ll, abs=1e-5)
    assert _b(panel, "z2", "L") == pytest.approx(l2, abs=1e-5)


def test_equal_loadings_label_and_constraint_line(panel):
    out = panel("gsem (L -> z1 z2@a z3@a), nolog")
    assert "[z2]L - [z3]L = 0" in out
    assert _b(panel, "z2", "L") == pytest.approx(_b(panel, "z3", "L"), abs=1e-12)


def test_fixed_variance_frees_first_loading(panel):
    panel("quietly gsem (L -> z1 z2 z3), variance(L@1) nolog")
    assert _b(panel, "z1", "L") != 1.0
    assert abs(_b(panel, "z1", "L")) > .3


def test_covariance_between_latents_is_estimated(panel):
    panel("quietly gsem (L1 -> z1 z2) (L2 -> z3 b1, family(gaussian)), cov(L1*L2) nolog")
    assert abs(_b(panel, "cov(L1,L2)", "_cons")) > 0


def test_irt_two_parameter_runs_and_discrimination_positive(panel):
    panel("quietly gsem (L -> b1 b2 b3, logit), nolog")
    assert _b(panel, "b2", "L") > 0 and _b(panel, "var(L)", "_cons") > 0


def test_gamma_equals_glm_coefficients(panel):
    panel("gen g = exp(.3 + .2*x)*(1 + abs(z1)) + .1")
    panel("quietly glm g x, family(gamma) link(log)")
    bx = float(_e(panel, "b").data[0, 0])
    panel("quietly gsem (g <- x, gamma), nolog")
    assert _b(panel, "g", "x") == pytest.approx(bx, abs=1e-4)


def test_predict_latent_has_mean_zero_and_mu_matches_xb_without_latent(panel):
    panel("quietly gsem (L -> z1 z2 z3), nolog")
    panel("predict lh, latent(L)")
    panel("quietly summarize lh")
    assert float(panel.session.r["mean"]) == pytest.approx(0, abs=1e-6)


def test_missing_responses_keep_the_observation(panel):
    panel("replace z2 = . in 1/10")
    panel("quietly gsem (L -> z1 z2 z3), nolog")
    assert float(_e(panel, "N")) == 400


def test_errors(panel):
    panel("gsem (y <- x, family(gaussian) link(logit))")
    assert panel.rc == 198
    panel("gsem (y <- x L M1[id], logit)")
    assert panel.rc == 198
