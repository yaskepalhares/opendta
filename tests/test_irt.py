"""irt: models checked against gsem, against the marginal likelihood computed by brute force at
the estimates, and the postestimation commands."""

import numpy as np
import pytest


def _load(run, cols):
    n = len(next(iter(cols.values())))
    run("clear\ninput " + " ".join("double " + k for k in cols) + "\n"
        + "\n".join(" ".join(repr(float(cols[k][i])) for k in cols) for i in range(n)) + "\nend")


@pytest.fixture
def items(run):
    rng = np.random.default_rng(8)
    N = 300
    th = rng.normal(size=N)
    cols = {}
    for i, (a, b) in enumerate([(1.2, -.8), (.8, -.2), (1.5, .3), (1.0, .9)]):
        cols[f"q{i + 1}"] = (rng.random(N) < 1 / (1 + np.exp(-a * (th - b)))).astype(float)
    for i in range(3):
        lat = th + rng.logistic(size=N) * .8
        cols[f"r{i + 1}"] = np.digitize(lat, [-1 + .4 * i, .3, 1.4]).astype(float)
    run.cols = cols
    _load(run, cols)
    return run


def _b(run, eq, name):
    e = run.session.e["b"]
    for j, (q, c) in enumerate(zip(e.coleq, e.colnames)):
        if q == eq and c == name:
            return float(e.data[0, j])
    raise KeyError((eq, name))


def _marginal(prob_fn, y, grid=np.linspace(-6, 6, 241)):
    """sum_i log integral prod_j P_j(y_ij | t) phi(t) dt on a fine grid."""
    w = np.exp(-grid ** 2 / 2) / np.sqrt(2 * np.pi) * (grid[1] - grid[0])
    L = np.ones((y.shape[0], len(grid)))
    for j in range(y.shape[1]):
        P = prob_fn(j, grid)
        L *= P[y[:, j].astype(int)]
    return float(np.sum(np.log(L @ w)))


def test_2pl_equals_gsem_and_brute_force(items):
    items("quietly irt 2pl q1-q4, nolog")
    ll = float(items.session.e["ll"])
    a = [_b(items, f"q{j}", "Theta") for j in range(1, 5)]
    c = [_b(items, f"q{j}", "_cons") for j in range(1, 5)]
    y = np.column_stack([items.cols[f"q{j}"] for j in range(1, 5)])
    f = lambda j, t: np.vstack([1 - 1 / (1 + np.exp(-(a[j] * t + c[j]))), 1 / (1 + np.exp(-(a[j] * t + c[j])))])   # noqa: E731
    assert ll == pytest.approx(_marginal(f, y), abs=.03), (ll, _marginal(f, y))
    items("quietly gsem (Theta -> q1 q2 q3 q4, logit), variance(Theta@1) nolog")
    assert float(items.session.e["ll"]) == pytest.approx(ll, abs=1e-6)
    # the likelihood is maximal: a perturbation lowers it
    for d in (.05, -.05):
        a2 = list(a)
        a2[1] += d
        g = lambda j, t: np.vstack([1 - 1 / (1 + np.exp(-(a2[j] * t + c[j]))), 1 / (1 + np.exp(-(a2[j] * t + c[j])))])   # noqa: E731
        assert _marginal(g, y) < ll


def test_2pl_table_reports_discrimination_and_difficulty(items):
    out = items("irt 2pl q1-q3, nolog")
    assert "Two-parameter logistic model" in out and "Number of obs     =        300" in out
    a, c = _b(items, "q1", "Theta"), _b(items, "q1", "_cons")
    first = [ln for ln in out.splitlines() if ln.strip().startswith("Diff")][0]
    assert float(first.split("|")[1].split()[0]) == pytest.approx(-c / a, abs=1e-5)
    assert items.session.e["cmd2"] == "irt" and items.session.e["model1"] == "2pl"


def test_1pl_has_one_common_discrimination(items):
    out = items("irt 1pl q1-q4, nolog")
    assert "One-parameter logistic model" in out
    assert _b(items, "q1", "Theta") == pytest.approx(_b(items, "q4", "Theta"), abs=1e-12)
    assert sum(1 for ln in out.splitlines() if "Discrim" in ln) == 1


def test_3pl_guessing_in_unit_interval_and_not_worse_than_2pl(items):
    items("quietly irt 2pl q1-q4, nolog")
    ll2 = float(items.session.e["ll"])
    out = items("irt 3pl q1-q4, nolog")
    assert float(items.session.e["ll"]) >= ll2 - 1e-6
    for ln in out.splitlines():
        if ln.strip().startswith("Guess"):
            g = float(ln.split("|")[1].split()[0])
            assert 0 <= g < 1


def test_3pl_common_guessing_unless_sepguessing(items):
    out = items("irt 3pl q1-q4, nolog")
    gs = {ln.split("|")[1].split()[0] for ln in out.splitlines() if ln.strip().startswith("Guess")}
    assert len(gs) == 1
    ll_common = float(items.session.e["ll"])
    out = items("irt 3pl q1-q4, sepguessing nolog")
    assert sum(ln.strip().startswith("Guess") for ln in out.splitlines()) == 4
    assert float(items.session.e["ll"]) >= ll_common - 1e-4


def test_grm_matches_brute_force(items):
    items("quietly irt grm r1-r3, nolog")
    ll = float(items.session.e["ll"])
    a = [_b(items, f"r{j}", "Theta") for j in range(1, 4)]
    k = [_b(items, f"r{j}", f"cut{c}") for j in (1, 2, 3) for c in (1, 2, 3)]
    y = np.column_stack([items.cols[f"r{j}"] for j in range(1, 4)])

    def f(j, t):
        kk = k[3 * j:3 * j + 3]
        cum = [np.ones_like(t)] + [1 / (1 + np.exp(-(a[j] * t - x))) for x in kk] + [np.zeros_like(t)]
        return np.vstack([cum[i] - cum[i + 1] for i in range(4)])
    assert ll == pytest.approx(_marginal(f, y), abs=.03), (ll, _marginal(f, y))


def test_gpcm_pcm_rsm_nrm_ordering_of_fit(items):
    lls = {}
    for m in ("rsm", "pcm", "gpcm", "nrm"):
        items(f"quietly irt {m} r1-r3, nolog")
        lls[m] = float(items.session.e["ll"])
    # nested models: rsm in pcm in gpcm in nrm
    assert lls["rsm"] <= lls["pcm"] + 1e-6 <= lls["gpcm"] + 2e-6 <= lls["nrm"] + 3e-6


def test_gpcm_matches_brute_force(items):
    items("quietly irt gpcm r1-r3, nolog")
    ll = float(items.session.e["ll"])
    e = items.session.e["b"]
    get = lambda q, c: float(e.data[0, [i for i, (a, b) in enumerate(zip(e.coleq, e.colnames)) if a == q and b == c][0]])   # noqa: E731
    y = np.column_stack([items.cols[f"r{j}"] for j in range(1, 4)])

    def f(j, t):
        a = get(f"r{j + 1}", "Theta")
        d = [get(f"r{j + 1}", f"{k}._cons") for k in (1, 2, 3)]
        eta = np.vstack([a * k * t + sum(d[:k]) for k in range(4)])
        p = np.exp(eta - eta.max(axis=0))
        return p / p.sum(axis=0)
    assert ll == pytest.approx(_marginal(f, y), abs=.03), (ll, _marginal(f, y))


def test_hybrid_and_estat_report(items):
    items("quietly irt hybrid (2pl q1-q3) (grm r1), nolog")
    # as in Stata, one items#/model# pair per group of items
    assert items.session.e["model1"] == "2pl" and items.session.e["model2"] == "grm"
    out = items("estat report q3 r1")
    assert "q3" in out and "r1" in out and "q1" not in out
    out = items("estat report, byparm sort(a) descending")
    first = [ln for ln in out.splitlines() if "|" in ln and ln.split("|")[1].strip()][1]
    vals = [float(ln.split("|")[1].split()[0]) for ln in out.splitlines()
            if ln.strip()[:2] in ("q1", "q2", "q3", "r1") and ln.split("|")[1].strip()][:4]
    assert vals == sorted(vals, reverse=True)


def test_irtgraph_keeps_the_data(items):
    items("quietly irt 2pl q1-q3, nolog")
    n, k = items.session.data.nobs, len(items.session.data.vars)
    for sub in ("icc", "iif q1 q2", "tcc", "tif"):
        items(f"irtgraph {sub}")
        assert items.rc == 0
    assert items.session.data.nobs == n and len(items.session.data.vars) == k


def test_errors(items):
    items("irt zzz q1-q3")
    assert items.rc == 198
    items("irt 2pl q1-q3, group(q4)")
    assert items.rc == 198
    items("irtgraph icc")
    assert items.rc in (301, 198)
