"""set numerics stata|precise."""

from fractions import Fraction as Fr

import numpy as np


def _exact_ols(X, y):
    k = X.shape[1]
    Xf = [[Fr(v) for v in r] for r in X.tolist()]
    yf = [Fr(v) for v in y.tolist()]
    n = len(Xf)
    A = [[sum(Xf[r][i] * Xf[r][j] for r in range(n)) for j in range(k)] + [sum(Xf[r][i] * yf[r] for r in range(n))]
         for i in range(k)]
    for i in range(k):
        for j in range(k):
            if j != i:
                f = A[j][i] / A[i][i]
                A[j] = [a - f * b for a, b in zip(A[j], A[i])]
    return np.array([float(A[i][k] / A[i][i]) for i in range(k)])


def test_precise_regress_is_more_accurate(run):
    run("clear\nset obs 30\ngen double x = 1000 + _n\ngen double x2 = x^2\ngen double x3 = x^3\n"
        "gen double y = 1 + x/7 + x2/1e3 + x3/1e7 + mod(_n,3)")
    ds = run.session.data
    X = np.column_stack([ds.get(n).data for n in ("x", "x2", "x3")] + [np.ones(30)])
    exact = _exact_ols(X, ds.get("y").data)
    errs = {}
    for mode in ("stata", "precise"):
        run(f"quietly set numerics {mode}\nquietly regress y x x2 x3")
        b = run.session.e["b"].data.ravel()
        errs[mode] = np.max(np.abs(b - exact) / np.abs(exact))
    assert errs["precise"] < 1e-8 and errs["precise"] < errs["stata"] / 100


def test_precise_differences_and_default_type(run):
    run("clear\nset obs 5\ngen t = _n\ngen x = t/10\nqui tsset t")
    run("set numerics stata\ngen a = 1")
    assert run.session.data.get("a").vtype == "float"
    run("set numerics precise\ngen b = 1")
    assert run.session.data.get("b").vtype == "double"
    from opendta.estimation import fvars
    comps = fvars._component_alts(run.session.data, "D.x", False)
    d_precise = fvars.comp_values(run.session.data, comps[0])
    run("set numerics stata")
    d_stata = fvars.comp_values(run.session.data, comps[0])
    x = run.session.data.get("x").data
    assert d_precise[1] == x[1] - x[0]
    assert d_stata[1] == float(np.float32(x[1] - x[0]))
    out = run("set numerics maybe")
    assert "invalid syntax" in out


def test_lags_vectorized_match_panels(run):
    run("clear\nset obs 9\ngen id = ceil(_n/3)\ngen t = mod(_n-1,3)*2\ngen y = _n\nqui xtset id t\n"
        "gen l = L2.y\ngen f = F2.y")
    ds = run.session.data
    l, f = ds.get("l").data, ds.get("f").data
    SYS = 8.98e307
    assert l[0] > SYS and l[1] == 1 and l[2] == 2 and l[3] > SYS
    assert f[0] == 2 and f[2] > SYS


def test_precise_ml_reaches_exact_maximum(run):
    run("clear\nset seed 611\nset obs 80\ngen x1 = round(runiform()*10, .1)\n"
        "gen y = runiform() < invlogit(-1 + .3*x1)")
    ds = run.session.data
    X = np.column_stack([ds.get("x1").data, np.ones(80)])
    y = ds.get("y").data
    grads = {}
    for mode in ("stata", "precise"):
        out = run(f"quietly set numerics {mode}\nlogit y x1")
        b = run.session.e["b"].data.ravel()
        p = 1 / (1 + np.exp(-X @ b))
        grads[mode] = np.abs(X.T @ (y - p)).max()
        assert "Iteration 0:" in out
    assert grads["precise"] < 1e-10 < grads["stata"]
