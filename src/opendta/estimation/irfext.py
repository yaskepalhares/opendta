"""irf after svar and vec, bootstrap standard errors, and irf ctable, cgraph, ograph and add
([TS] irf create, [TS] irf ctable, [TS] irf cgraph, [TS] irf ograph, [TS] irf add)."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..lang.words import split_words

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


# ---------------------------------------------------------------------------
# pieces for other models
# ---------------------------------------------------------------------------

def vec_inputs(est) -> tuple[list[str], list[np.ndarray], np.ndarray]:
    """(variables, A_1..A_p of the levels VAR, residual covariance) implied by a VECM."""
    from .vec import var_matrices
    v = est.vec
    K = len(v["deps"])
    A = var_matrices(v["coef"], v["alpha"], v["beta"], K, v["lags"])
    return v["deps"], A, v["Om"]


def orthogonalizer(est, Sig: np.ndarray) -> np.ndarray:
    """Matrix P with Sigma = P P' used for orthogonalized responses: A^-1 B after svar, the
    Cholesky factor otherwise."""
    if hasattr(est, "svar"):
        return np.linalg.solve(est.svar["A"], est.svar["B"])
    return np.linalg.cholesky(Sig)


# ---------------------------------------------------------------------------
# bootstrap
# ---------------------------------------------------------------------------

def _simulate_var(V: dict, U: np.ndarray) -> np.ndarray:
    """Dynamic simulation of the fitted VAR from its first observed lags with innovations U."""
    deps, xnames, B = V["deps"], V["xnames"], V["B"]
    K = len(deps)
    p = max(l for _v, l in xnames if l > 0)
    first = V["X"][0]
    hist = [np.zeros(K) for _ in range(p)]                  # hist[l - 1] = y_{t0 - l}
    for i, (vn, l) in enumerate(xnames):
        if l > 0:
            hist[l - 1][deps.index(vn)] = first[i]
    seq = hist[::-1]                                         # oldest first
    T = U.shape[0]
    out = np.zeros((T, K))
    for t in range(T):
        x = np.array([1.0 if vn == "_cons" else seq[-l][deps.index(vn)] for vn, l in xnames])
        y = x @ B + U[t]
        out[t] = y
        seq.append(y)
    return out


def _refit_var(V: dict, Ysim: np.ndarray):
    """OLS re-estimation on the simulated sample (same lags and constant)."""
    deps, xnames = V["deps"], V["xnames"]
    p = max(l for _v, l in xnames if l > 0)
    first = V["X"][0]
    K = len(deps)
    hist = [np.zeros(K) for _ in range(p)]
    for i, (vn, l) in enumerate(xnames):
        if l > 0:
            hist[l - 1][deps.index(vn)] = first[i]
    full = np.vstack(hist[::-1] + [Ysim])
    T = Ysim.shape[0]
    X = np.zeros((T, len(xnames)))
    for i, (vn, l) in enumerate(xnames):
        X[:, i] = 1.0 if vn == "_cons" else full[p - l:p - l + T, deps.index(vn)]
    B = np.linalg.lstsq(X, Ysim, rcond=None)[0]
    U = Ysim - X @ B
    Sig = U.T @ U / (T - X.shape[1]) if V["dfk"] else U.T @ U / T
    return B, Sig


def _simulate_vecm(v: dict, U: np.ndarray) -> np.ndarray:
    """Dynamic simulation of the fitted VECM; returns the whole path (presample included)."""
    Y, p, trend, coef, beta = v["Y"], v["lags"], v["trend"], v["coef"], v["beta"]
    n, K = Y.shape
    path = [Y[i].copy() for i in range(p)]
    r = beta.shape[1]
    for t in range(p, n):
        z1 = path[-1]
        if trend == "rconstant":
            z1 = np.append(z1, 1.0)
        elif trend == "rtrend":
            z1 = np.append(z1, float(t + 1))
        w = list(z1 @ beta)
        for i in range(1, p):
            w += list(path[-i] - path[-i - 1])
        if trend in ("constant", "trend", "rtrend"):
            w.append(1.0)
        if trend == "trend":
            w.append(float(t + 1))
        dy = np.array(w) @ coef + U[t - p]
        path.append(path[-1] + dy)
    return np.array(path)


def bootstrap_se(s: "Session", est, steps: int, perm: list[int], reps: int, parametric: bool,
                 dots: bool, stats_fn, seed_rng=None) -> dict[str, np.ndarray]:
    """Standard errors of the irf statistics by simulation. `stats_fn(A, Sig, P)` returns the
    statistics dictionary of arrays (steps + 1, K, K) for the lag matrices A."""
    from .resample import _dot, _dots
    rng = seed_rng or np.random.default_rng()
    draws = []
    if dots:
        _dots(s, "Bootstrap replications", reps, {})
    kind = est.cmd
    for rep in range(reps):
        ok = True
        try:
            if kind == "vec":
                v = est.vec
                U0 = v["J"]["Z0"] - np.hstack([v["J"]["Z1"] @ v["beta"], v["J"]["Z2"]]) @ v["coef"]
                T = U0.shape[0]
                U = rng.multivariate_normal(np.zeros(U0.shape[1]), v["Om"], size=T) if parametric \
                    else U0[rng.integers(0, T, T)]
                from .vec import fit_vecm, var_matrices
                F = fit_vecm(_simulate_vecm(v, U), v["lags"], v["trend"], v["r"])
                K = len(v["deps"])
                A = var_matrices(F["coef"], F["alpha"], F["beta"], K, v["lags"])
                draws.append(stats_fn(A, F["Om"], None))
            else:
                V = est.var
                U0 = V["U"]
                T = U0.shape[0]
                U = rng.multivariate_normal(np.zeros(U0.shape[1]), V["Sig"], size=T) if parametric \
                    else U0[rng.integers(0, T, T)]
                Ysim = _simulate_var(V, U)
                B, Sig = _refit_var(V, Ysim)
                from .irf import _lag_matrices
                A = _lag_matrices(B, V["xnames"], V["deps"], V["lags"])
                P = None
                if hasattr(est, "svar"):
                    from .svar import fit_svar, unpack
                    sv = est.svar
                    res, _V = fit_svar(Sig, T, sv["Afix"], sv["Bfix"])
                    Am, Bm, *_ = unpack(res.b, sv["Afix"], sv["Bfix"])
                    P = np.linalg.solve(Am, Bm)
                draws.append(stats_fn(A, Sig, P))
        except (np.linalg.LinAlgError, StataError):
            ok = False
        if dots:
            _dot(s, ok, rep, reps)
    if not draws:
        raise StataError(498, "bootstrap failed")      # VERIFY
    return {k: np.std(np.array([d[k] for d in draws]), axis=0, ddof=1) for k in draws[0]}


# ---------------------------------------------------------------------------
# irf add, ctable, cgraph, ograph
# ---------------------------------------------------------------------------

def parse_specs(text: str) -> tuple[list[list[str]], str]:
    """'(a y1 y2 oirf) (b y2 y1 irf, noci) , opts' -> ([['a','y1','y2','oirf'], ...], 'opts')."""
    specs, depth, cur, rest = [], 0, "", ""
    i = 0
    t = text.strip()
    while i < len(t):
        ch = t[i]
        if depth == 0 and ch == "(":
            depth, cur = 1, ""
        elif depth == 1 and ch == ")":
            depth = 0
            body = cur.split(",")[0]
            specs.append(split_words(body))
        elif depth == 1:
            cur += ch
        elif ch == ",":
            rest = t[i + 1:]
            break
        elif not ch.isspace():
            raise StataError(198, "invalid syntax")
        i += 1
    for sp in specs:
        if len(sp) != 4:
            raise StataError(198, "each group needs irfname impulse response statistic")   # VERIFY
    if not specs:
        raise StataError(198, "invalid syntax")
    return specs, rest


def cmd_add(s: "Session", rest: str) -> None:
    from ..lang.syntax import match_options, parse_standard
    from .irf import _active, _irf_path, _load, _need_file, _save
    p = parse_standard(rest)
    o = match_options(p.options, {"using": 1})
    if not isinstance(o.get("using"), str):
        raise StataError(198, "option using() required")
    src = _irf_path(o["using"])
    if not src.exists():
        raise StataError(601, f"file {src.name} not found")
    other = _load(src)
    data = _need_file(s)
    names = split_words(p.varlist) or list(other)
    if names == ["_all"]:
        names = list(other)
    for n in names:
        if n not in other:
            raise StataError(111, f"irfname {n} not found in {src.name}")     # VERIFY
        if n in data:
            raise StataError(110, f"irfname {n} already exists")               # VERIFY
        data[n] = other[n]
    _save(_active(s), data)
    s.output.write(f"(file {_active(s).name} updated)\n", "text")


def cmd_ctable(s: "Session", rest: str) -> None:
    from .irf import STATS, _need_file, _table_render
    specs, opts = parse_specs(rest)
    data = _need_file(s)
    picks = []
    for n, a, r, k in specs:
        if n not in data:
            raise StataError(111, f"irfname {n} not found")       # VERIFY
        deps = data[n]["meta"]["vars"]
        if a not in deps or r not in deps:
            raise StataError(111, f"{a} or {r} not found in {n}")  # VERIFY
        picks.append((n, a, r))
    from ..lang.syntax import match_options
    o = match_options(opts, {"noci": 4, "level": 1, "stderror": 3, "individual": 3, "title": 3, "step": 4}) \
        if opts.strip() else {}
    for _n, _a, _r, k in specs:
        if k.lower() not in STATS:
            raise StataError(198, f"statistic {k} not allowed")       # VERIFY
    _table_render(s, data, picks, [[k.lower()] for _n, _a, _r, k in specs], o)


def cmd_cgraph(s: "Session", rest: str, overlay: bool) -> None:
    from ..commands.graphcmd import finish
    from ..graphics.model import Graph, Panel
    from ..graphics.opts import graph_options, split_layer_graph
    from .irf import STATS, _need_file, _spec_layers
    specs, opts = parse_specs(rest)
    lopts, graw = split_layer_graph(opts, {"noci": 4, "level": 1, "individual": 3})
    data = _need_file(s)
    level = float(lopts["level"]) if isinstance(lopts.get("level"), str) else float(s.settings.get("level", "95"))
    panels, layers = [], []
    for q, (n, a, r, k) in enumerate(specs):
        if n not in data:
            raise StataError(111, f"irfname {n} not found")       # VERIFY
        if k.lower() not in STATS:
            raise StataError(198, f"statistic {k} not allowed")  # VERIFY
        ls = _spec_layers(data, n, a, r, k.lower(), q if overlay else 0, level, bool(lopts.get("noci")))
        if overlay:
            layers += ls
        else:
            panels.append(Panel(f"{n}, {a}, {r}", ls))
    if overlay:
        panels = [Panel("", layers)]
    gopts = graph_options(graw)
    gopts["_xtitle"] = "step"
    finish(s, Graph(family="twoway", opts=gopts, panels=panels), graw)
