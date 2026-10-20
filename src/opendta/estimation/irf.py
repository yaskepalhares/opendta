"""irf ([TS] irf): impulse-response functions and forecast-error variance decompositions after
var.

    irf set filename [, replace]          active results file (.irf)
    irf create irfname [, set(filename) replace step(#) order(varlist) nose level(#)]
    irf table [stats] [, irf(names) impulse(varlist) response(varlist) noci stderror level(#)]
    irf graph [stats] [, irf(names) impulse(varlist) response(varlist) noci level(#)]
    irf describe | irf drop names | irf rename old new

Statistics: irf (simple), oirf (orthogonalized with the Cholesky factor of Sigma in the order
of the variables, or order()), cirf and coirf (cumulative), fevd (Cholesky forecast-error
variance decomposition). From the VAR coefficients, Phi_0 = I and Phi_i = sum_j Phi_{i-j} A_j
(Luetkepohl 2005, sec. 2.3); asymptotic standard errors by the delta method with the variance
of the coefficients and of vech(Sigma), 2 D+ (Sigma x Sigma) D+' / T, computed with numerical
derivatives (Methods and formulas of [TS] irf create cite the same results).

The .irf file is a .dta dataset with one observation per (irfname, impulse, response, step);
its variables are OpenDTA's own (VERIFY: Stata's layout of .irf files is undocumented), as is
the layout of the tables (case compat/do/1303_irf.do)."""

from __future__ import annotations

import json
from pathlib import Path
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.formats import format_value
from ..lang.syntax import match_options, parse_standard
from ..lang.words import split_words
from ..lazy import lazy

st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS
STATS = ("irf", "oirf", "cirf", "coirf", "fevd")
_LABELS = {"irf": "irf", "oirf": "oirf", "cirf": "cirf", "coirf": "coirf", "fevd": "fevd"}


def _g(x):
    if x is None or not np.isfinite(x) or abs(x) >= SYS:
        return "."
    return format_value(float(x), "%9.0g", pad=False).strip()


# ---------------------------------------------------------------------------
# computation
# ---------------------------------------------------------------------------

def _lag_matrices(B: np.ndarray, xnames, deps, lags) -> list[np.ndarray]:
    """A_1..A_p (K x K) from the coefficient matrix B (rows = regressors, cols = equations)."""
    K = len(deps)
    p = max(lags)
    A = [np.zeros((K, K)) for _ in range(p)]
    for r, (vn, l) in enumerate(xnames):
        if l > 0 and vn in deps:
            A[l - 1][:, deps.index(vn)] = B[r, :]
    return A


def _phis(A: list[np.ndarray], steps: int) -> list[np.ndarray]:
    K = A[0].shape[0] if A else 1
    out = [np.eye(K)]
    for i in range(1, steps + 1):
        acc = np.zeros((K, K))
        for j in range(1, min(i, len(A)) + 1):
            acc += out[i - j] @ A[j - 1]
        out.append(acc)
    return out


def _all_stats(A, Sig, steps, P=None):
    """Arrays [step, response, impulse] for each statistic. P: matrix with Sigma = P P' for the
    orthogonalized responses (default: the Cholesky factor; svar: A^-1 B)."""
    K = Sig.shape[0]
    Ph = np.array(_phis(A, steps))
    if P is None:
        try:
            P = np.linalg.cholesky(Sig)
        except np.linalg.LinAlgError:
            P = np.linalg.cholesky(Sig + 1e-12 * np.eye(K))
    Th = Ph @ P
    res = {"irf": Ph, "oirf": Th, "cirf": np.cumsum(Ph, axis=0), "coirf": np.cumsum(Th, axis=0)}
    fevd = np.zeros_like(Ph)
    sq = np.cumsum(Th ** 2, axis=0)          # sum_{i<=h} theta_{jk,i}^2
    mse = sq.sum(axis=2)
    for h in range(1, steps + 1):
        # step h: forecast horizon h, i.e. the first h MA terms (step 0 = 0)
        fevd[h] = sq[h - 1] / mse[h - 1][:, None]
    res["fevd"] = fevd
    return res


def _vech_idx(K):
    return [(i, j) for j in range(K) for i in range(j, K)]


def _stats_of(A, Sig, steps, perm, P=None):
    """Statistics with the variables ordered by `perm` for the Cholesky factorization."""
    Ap = [a[np.ix_(perm, perm)] for a in A]
    out = _all_stats(Ap, Sig[np.ix_(perm, perm)], steps, P)
    inv = np.argsort(perm)
    return {k: x[:, inv][:, :, inv] for k, x in out.items()}


def _compute(est, steps: int, order: list[str], se: bool):
    if est.cmd == "vec":
        from .irfext import vec_inputs
        deps, A, Sig = vec_inputs(est)
        point = _stats_of(A, Sig, steps, [deps.index(n) for n in order])
        return deps, point, {k: np.full_like(x, np.nan) for k, x in point.items()}
    v = est.var
    deps = list(v["deps"])
    lags = list(v["lags"])
    xnames = v["xnames"]
    B = v["B"]
    Sig = v["Sig"]
    T = v["T"]
    perm = [deps.index(n) for n in order]
    # coefficients of the lags, as a vector, with their covariance from e(V)
    lag_rows = [r for r, (vn, l) in enumerate(xnames) if l > 0 and vn in deps]
    m = len(xnames)
    K = len(deps)
    pos = [j * m + r for j in range(K) for r in lag_rows]
    theta_a = np.array([B[r, j] for j in range(K) for r in lag_rows])
    Vall = est.V if est.cmd == "var" else np.kron(Sig, v["XtXi"])
    Va = Vall[np.ix_(pos, pos)]
    vi = _vech_idx(K)
    sig_v = np.array([Sig[i, j] for i, j in vi])
    sv = getattr(est, "svar", None)
    # var(vech Sigma) = 2 D+ (Sigma x Sigma) D+' / T
    Vs = np.zeros((len(vi), len(vi)))
    for a, (i, j) in enumerate(vi):
        for b, (k, l) in enumerate(vi):
            Vs[a, b] = (Sig[i, k] * Sig[j, l] + Sig[i, l] * Sig[j, k]) / T

    def build(ta, sv_):
        Bc = B.copy()
        q = 0
        for j in range(K):
            for r in lag_rows:
                Bc[r, j] = ta[q]
                q += 1
        A = _lag_matrices(Bc, xnames, deps, lags)
        if sv is not None:
            # svar: orthogonalize with A^-1 B (order() does not apply)
            from .svar import unpack
            Am, Bm, *_ = unpack(sv_, sv["Afix"], sv["Bfix"])
            return _stats_of(A, Sig, steps, list(range(K)), np.linalg.solve(Am, Bm))
        S = np.zeros((K, K))
        for val, (i, j) in zip(sv_, vi):
            S[i, j] = S[j, i] = val
        # order(): permute the variables before the Cholesky factorization
        return _stats_of(A, S, steps, perm)

    if sv is not None:
        sig_v, Vs = sv["phi"], sv["V"]
    point = build(theta_a, sig_v)
    ses = {k: np.full_like(x, np.nan) for k, x in point.items()}
    if se:
        full = np.concatenate([theta_a, sig_v])
        na = len(theta_a)
        Vfull = np.zeros((len(full), len(full)))
        Vfull[:na, :na] = Va
        Vfull[na:, na:] = Vs
        flat0 = np.concatenate([point[k].ravel() for k in STATS])
        J = np.zeros((len(flat0), len(full)))
        for i in range(len(full)):
            h = 1e-6 * (abs(full[i]) + 1e-3)
            e = np.zeros(len(full))
            e[i] = h
            up = build((full + e)[:na], (full + e)[na:])
            dn = build((full - e)[:na], (full - e)[na:])
            J[:, i] = (np.concatenate([up[k].ravel() for k in STATS])
                       - np.concatenate([dn[k].ravel() for k in STATS])) / (2 * h)
        var = np.einsum("ij,jk,ik->i", J, Vfull, J)
        sd = np.sqrt(np.maximum(var, 0))
        off = 0
        for k in STATS:
            n = point[k].size
            ses[k] = sd[off:off + n].reshape(point[k].shape)
            off += n
        # fixed values at step 0 (irf = I, fevd = 0): no sampling variance
        ses["irf"][0] = 0.0
        ses["cirf"][0] = 0.0
        ses["fevd"][0] = 0.0
    return deps, point, ses


# ---------------------------------------------------------------------------
# results file
# ---------------------------------------------------------------------------

def _active(s: "Session") -> Path | None:
    p = getattr(s, "_irf_file", None)
    return Path(p) if p else None


def _irf_path(name: str) -> Path:
    from ..lang.words import strip_outer_quotes
    t = strip_outer_quotes(name.strip())
    p = Path(t)
    if not p.suffix:
        p = p.with_suffix(".irf")
    return p


def _load(path: Path) -> dict:
    """{irfname: {"meta": {...}, "stats": {stat: array}, "se": {...}}}"""
    if not path.exists():
        return {}
    from ..io.dta import read_dta
    ds = read_dta(path)
    meta = json.loads(ds.chars.get("_dta", {}).get("opendta_irf", "{}") or "{}")
    out = {}
    names = list(ds.get("irfname").raw) if ds.nobs else []
    for name, m_ in meta.items():
        deps = m_["vars"]
        K, steps = len(deps), m_["step"]
        sel = np.array([n == name for n in names], dtype=bool)
        imp = np.asarray(ds.get("impulse").raw)[sel]
        resp = np.asarray(ds.get("response").raw)[sel]
        stp = ds.get("step").data[sel].astype(int)
        res = {"meta": m_, "stats": {}, "se": {}}
        for k in STATS:
            arr = np.full((steps + 1, K, K), np.nan)
            sarr = np.full((steps + 1, K, K), np.nan)
            vals = ds.get(k).data[sel]
            svals = ds.get("std" + k).data[sel]
            for q in range(len(stp)):
                i, j = deps.index(resp[q]), deps.index(imp[q])
                arr[stp[q], i, j] = vals[q] if vals[q] < SYS else np.nan
                sarr[stp[q], i, j] = svals[q] if svals[q] < SYS else np.nan
            res["stats"][k] = arr
            res["se"][k] = sarr
        out[name] = res
    return out


def _save(path: Path, data: dict) -> None:
    from ..core.dataset import Dataset, Variable
    from ..io.dta import write_dta
    rows = []
    for name, res in data.items():
        deps = res["meta"]["vars"]
        steps = res["meta"]["step"]
        for j, imp in enumerate(deps):
            for i, resp in enumerate(deps):
                for h in range(steps + 1):
                    rows.append((name, imp, resp, h, [(res["stats"][k][h, i, j], res["se"][k][h, i, j])
                                                      for k in STATS]))
    ds = Dataset()
    ds.nobs = len(rows)

    def strvar(nm, vals):
        w = max([len(v) for v in vals] + [1])
        return Variable(nm, f"str{w}", np.array(vals, dtype=object))
    ds.vars.append(strvar("irfname", [r[0] for r in rows]))
    ds.vars.append(strvar("impulse", [r[1] for r in rows]))
    ds.vars.append(strvar("response", [r[2] for r in rows]))
    ds.vars.append(Variable("step", "int", np.array([float(r[3]) for r in rows])))
    for q, k in enumerate(STATS):
        for prefix, idx in (("", 0), ("std", 1)):
            vals = np.array([r[4][q][idx] for r in rows], dtype=float)
            vals = np.where(np.isfinite(vals), vals, SYS)
            ds.vars.append(Variable(prefix + k, "double", vals))
    ds.chars.setdefault("_dta", {})["opendta_irf"] = json.dumps({n: r["meta"] for n, r in data.items()})
    ds.label = "irf results"
    write_dta(ds, path)


# ---------------------------------------------------------------------------
# command
# ---------------------------------------------------------------------------

_SUBS = (("set", 3), ("create", 2), ("table", 1), ("graph", 1), ("describe", 1), ("drop", 4),
         ("rename", 3), ("ctable", 2), ("cgraph", 2), ("ograph", 1), ("add", 3))


def cmd_irf(s: "Session", args: str) -> None:
    t = args.strip()
    if not t:
        raise StataError(198, "subcommand required")   # VERIFY
    head, _, rest = t.partition(" ")
    if "," in head:
        head, _, tail = head.partition(",")
        rest = "," + tail + " " + rest
    sub = next((n for n, m in _SUBS if len(head) >= m and n.startswith(head)), None)
    if sub is None:
        raise StataError(198, f"invalid subcommand {head}")   # VERIFY
    if sub in ("ctable", "cgraph", "ograph", "add"):
        from . import irfext
        if sub == "add":
            irfext.cmd_add(s, rest)
        elif sub == "ctable":
            irfext.cmd_ctable(s, rest)
        else:
            irfext.cmd_cgraph(s, rest, overlay=(sub == "ograph"))
        return
    {"set": _set, "create": _create, "table": _table, "graph": _graph, "describe": _describe,
     "drop": _drop, "rename": _rename}[sub](s, rest)


def _set(s: "Session", rest: str) -> None:
    p = parse_standard(rest)
    o = match_options(p.options, {"replace": 7, "clear": 5})
    if o.get("clear"):
        s._irf_file = None
        return
    name = p.varlist.strip() or (p.using or "")
    if not name:
        cur = _active(s)
        s.output.write(f"\n(active irf file: {cur if cur else 'none'})\n", "text")   # VERIFY
        return
    path = _irf_path(name)
    if o.get("replace") and path.exists():
        path.unlink()
    if path.exists():
        s.output.write(f"(file {path.name} now active)\n", "text")
    else:
        s.output.write(f"(file {path.name} created)\n(file {path.name} now active)\n", "text")
        _save(path, {})
    s._irf_file = str(path)


def _create(s: "Session", rest: str) -> None:
    from .postest import current
    p = parse_standard(rest)
    o = match_options(p.options, {"set": 3, "replace": 7, "step": 4, "order": 3, "nose": 4,
                                  "level": 1, "bs": 2, "bsp": 3, "reps": 4, "nodots": 4})
    name = p.varlist.strip()
    if not name or len(name.split()) != 1:
        raise StataError(198, "irfname required")   # VERIFY
    est = current(s)
    if est.cmd not in ("var", "svar", "vec"):
        raise StataError(321, "irf create works after var, svar and vec")   # VERIFY
    if isinstance(o.get("set"), str):
        _set(s, o["set"])
    path = _active(s)
    if path is None:
        raise StataError(198, "no irf file active; use irf set or the set() option")   # VERIFY
    steps = int(float(o["step"])) if isinstance(o.get("step"), str) else 8
    mv = est.vec if est.cmd == "vec" else est.var
    deps = list(mv["deps"])
    order = s.expand_varlist(o["order"]) if isinstance(o.get("order"), str) else deps
    if sorted(order) != sorted(deps):
        raise StataError(198, "order() must contain the endogenous variables")   # VERIFY
    boot = bool(o.get("bs") or o.get("bsp"))
    deps, point, ses = _compute(est, steps, order, not (o.get("nose") or boot))
    if boot:
        from .irfext import bootstrap_se
        perm = [deps.index(n) for n in order]
        reps = int(float(o["reps"])) if isinstance(o.get("reps"), str) else 200
        ses = bootstrap_se(s, est, steps, perm, reps, bool(o.get("bsp")), not o.get("nodots"),
                           lambda A, Sig, P: _stats_of(A, Sig, steps, perm if P is None else list(range(len(deps))), P))
        for k in ("irf", "cirf", "fevd"):
            ses[k][0] = 0.0
    std = "bootstrap" if boot else ("none" if o.get("nose") or est.cmd == "vec" else "asymptotic")
    data = _load(path)
    if name in data and not o.get("replace"):
        raise StataError(110, f"irfname {name} already exists in {path.name}")   # VERIFY
    meta = {"model": est.cmd, "vars": deps, "order": order, "step": steps, "lags": mv["lags"],
            "exog": mv.get("exog", []), "std": std, "cmdline": s.e.get("cmdline", "")}
    if o.get("replace") and name not in data:
        s.output.write(f"irfname {name} not found in {path.name}\n", "text")
    data[name] = {"meta": meta, "stats": point, "se": ses}
    _save(path, data)
    s.output.write(f"(file {path.name} updated)\n", "text")


def _select(s, data, o):
    names = split_words(str(o["irf"])) if isinstance(o.get("irf"), str) else list(data)
    for n in names:
        if n not in data:
            raise StataError(111, f"irfname {n} not found")   # VERIFY
    picks = []
    for n in names:
        deps = data[n]["meta"]["vars"]
        imp = split_words(str(o["impulse"])) if isinstance(o.get("impulse"), str) else deps
        resp = split_words(str(o["response"])) if isinstance(o.get("response"), str) else deps
        for a in imp:
            if a not in deps:
                raise StataError(111, f"impulse {a} not found in {n}")   # VERIFY
        for r in resp:
            if r not in deps:
                raise StataError(111, f"response {r} not found in {n}")   # VERIFY
        for a in imp:
            for r in resp:
                picks.append((n, a, r))
    return picks


def _need_file(s) -> dict:
    path = _active(s)
    if path is None or not path.exists():
        raise StataError(198, "no irf file active")   # VERIFY
    data = _load(path)
    if not data:
        raise StataError(198, f"no irf results in {path.name}")   # VERIFY
    return data


def _six(x: float) -> str:
    """irf table's numbers: six digits after the leading "0." of numbers below 1, six
    significant digits otherwise, trailing zeros dropped (.479762, .06157, 1.53619;
    compat 1303; VERIFY numbers of 1e5 and more)."""
    if not np.isfinite(x):
        return "."
    ax = abs(x)
    ndig = 0 if ax < 1 else len(str(int(ax)))
    t = f"{x:.{max(6 - ndig, 0)}f}"
    if "." in t:
        t = t.rstrip("0").rstrip(".")
    if t.startswith("0."):
        t = t[1:]
    elif t.startswith("-0."):
        t = "-" + t[2:]
    return "0" if t in ("", "-0", "-") else t


def _table(s: "Session", rest: str) -> None:
    p = parse_standard(rest)
    o = match_options(p.options, {"irf": 3, "impulse": 3, "response": 3, "noci": 4,
                                  "stderror": 3, "level": 1, "individual": 3, "title": 3, "step": 4})
    stats = [w.lower() for w in p.varlist.split()] or list(STATS)
    for k in stats:
        if k not in STATS:
            raise StataError(198, f"statistic {k} not allowed")   # VERIFY: sirf, dm, cdm...
    data = _need_file(s)
    picks = _select(s, data, o)
    _table_render(s, data, picks, [stats] * len(picks), o)


def _table_render(s: "Session", data: dict, picks: list, stat_lists: list, o: dict) -> None:
    level = float(o["level"]) if isinstance(o.get("level"), str) else float(s.settings.get("level", "95"))
    z = float(st.norm.isf((1 - level / 100) / 2))
    out = s.output
    names = list(dict.fromkeys(n for n, _a, _r in picks))
    # one group of columns per (irf result, statistic): the estimate, then its limits or its
    # standard error; groups separated by "|" (compat 1303)
    groups = []
    for q, pick in enumerate(picks, start=1):
        for k in stat_lists[q - 1]:
            cols = [(q, k, "est")]
            if not o.get("noci"):
                cols += [(q, k, "lo"), (q, k, "hi")]
            if o.get("stderror"):
                cols.append((q, k, "se"))
            groups.append(cols)
    steps = min(data[n]["meta"]["step"] for n in names)
    if isinstance(o.get("step"), str):
        steps = min(steps, int(float(o["step"])))
    per = 6                      # VERIFY: columns per table
    tables, cur = [], []
    for g in groups:
        if cur and sum(len(x) for x in cur) + len(g) > per:
            tables.append(cur)
            cur = []
        cur.append(g)
    tables.append(cur)
    heads = {"est": None, "lo": "Lower", "hi": "Upper", "se": "S.E."}

    def cells(texts):
        # 12-character cells, the last of each group one shorter
        return "".join(t.ljust(12) for t in texts)[:12 * len(texts) - 1]
    title = f"Results from {' '.join(names)}"
    for ti, tab in enumerate(tables):
        width = 10 + sum(12 * len(g) for g in tab)
        if ti == 0:
            out.write("\n" + " " * max(0, (width - len(title)) // 2 - 1) + title + "\n", "text")
        out.write("\n+" + "-" * (width - 2) + "+\n", "text")
        out.write("|        |" + "|".join(cells([f"{'(' + str(q) + ')':>7}" for q, _k, _t in g]) for g in tab)
                  + "|\n", "text")
        # names start at (10 - length) // 2 within the cell
        out.write("|  step  |" + "|".join(cells([" " * ((10 - len(nm)) // 2) + nm for nm in
                                                 [k if heads[t] is None else heads[t] for _q, k, t in g]])
                                         for g in tab) + "|\n", "text")
        out.write("|--------+" + "+".join("-" * (12 * len(g) - 1) for g in tab) + "|\n", "text")
        for h in range(steps + 1):
            parts = []
            for g in tab:
                vals = []
                for q, k, t in g:
                    n, a, r = picks[q - 1]
                    res = data[n]
                    deps = res["meta"]["vars"]
                    i, j = deps.index(r), deps.index(a)
                    b = res["stats"][k][h, i, j]
                    se_ = res["se"][k][h, i, j]
                    v = {"est": b, "lo": b - z * se_, "hi": b + z * se_, "se": se_}[t]
                    vals.append(" " + _six(v))
                parts.append(cells(vals))
            out.write(f"|{h:<8}|", "text")
            out.write("|".join(parts), "result")
            out.write("|\n", "text")
        out.write("+" + "-" * (width - 2) + "+\n", "text")
    if not o.get("noci"):
        out.write(f"{level:g}% lower and upper bounds reported\n", "text")
    for q, (n, a, r) in enumerate(picks, start=1):
        out.write(f"({q}) irfname = {n}, impulse = {a}, and response = {r}\n", "text")


def _spec_layers(data: dict, n: str, a: str, r: str, k: str, pal: int, level: float, noci: bool) -> list:
    """Confidence band and line of one statistic of one (irfname, impulse, response)."""
    from ..graphics.model import Layer
    z = float(st.norm.isf((1 - level / 100) / 2))
    res = data[n]
    deps = res["meta"]["vars"]
    i, j = deps.index(r), deps.index(a)
    xs = list(range(res["meta"]["step"] + 1))
    b = res["stats"][k][:, i, j]
    se_ = res["se"][k][:, i, j]
    layers = []
    if not noci and np.isfinite(se_).any():
        layers.append(Layer("rarea", xs, list(b - z * se_), {"y2": list(b + z * se_)},
                            {"_pal": pal, "_ci": True}, f"{level:g}% CI"))
    layers.append(Layer("line", xs, list(b), {}, {"_pal": pal}, _LABELS[k]))
    return layers


def _graph(s: "Session", rest: str) -> None:
    from ..commands.graphcmd import finish
    from ..graphics.model import Graph, Panel
    from ..graphics.opts import graph_options, split_layer_graph
    p = parse_standard(rest)
    lopts, graw = split_layer_graph(p.options, {"irf": 3, "impulse": 3, "response": 3, "noci": 4,
                                                "level": 1, "stderror": 3, "individual": 3})
    stats = [w.lower() for w in p.varlist.split()] or ["irf"]
    for k in stats:
        if k not in STATS:
            raise StataError(198, f"statistic {k} not allowed")   # VERIFY
    data = _need_file(s)
    picks = _select(s, data, lopts)
    level = float(lopts["level"]) if isinstance(lopts.get("level"), str) else float(s.settings.get("level", "95"))
    panels = []
    for n, a, r in picks:
        layers = []
        for pal, k in enumerate(stats):
            layers += _spec_layers(data, n, a, r, k, pal, level, bool(lopts.get("noci")))
        panels.append(Panel(f"{n}, {a}, {r}", layers))
    opts = graph_options(graw)
    opts["_xtitle"] = "step"
    opts.setdefault("note", {"lines": ["Graphs by irfname, impulse variable, and response variable"]})
    g = Graph(family="twoway", opts=opts, panels=panels)
    finish(s, g, graw)


def _describe(s: "Session", rest: str) -> None:
    data = _need_file(s)
    path = _active(s)
    out = s.output
    import datetime
    stamp = datetime.datetime.fromtimestamp(path.stat().st_mtime).strftime("%d %b %Y %H:%M")
    width = int(s.settings.get("linesize", 80))
    # rules as wide as the line (compat 1303)
    out.write(f"Contains irf results from {path.name} (dated {stamp})\n\n", "text")
    out.write(f"{'irfname':>14} | model    endogenous variables and order (*)\n", "text")
    out.write("-" * 15 + "+" + "-" * (width - 16) + "\n", "text")
    for n, res in data.items():
        m_ = res["meta"]
        out.write(f"{n:>14} | ", "text")
        out.write(f"{'sr svar' if m_['model'] == 'svar' else m_['model']:<9}" + " ".join(m_["order"]) + "\n", "result")
    out.write("-" * width + "\n", "text")
    out.write("  (*) order is relevant only when model is var\n", "text")
    s.r = {"irfnames": " ".join(data)}


def _drop(s: "Session", rest: str) -> None:
    data = _need_file(s)
    path = _active(s)
    names = split_words(rest.split(",")[0])
    if not names:
        raise StataError(198, "irfname required")   # VERIFY
    for n in names:
        if n not in data:
            raise StataError(111, f"irfname {n} not found")   # VERIFY
        del data[n]
        s.output.write(f"({n} dropped)\n", "text")
    _save(path, data)
    s.output.write(f"file {path.name} updated\n", "text")


def _rename(s: "Session", rest: str) -> None:
    data = _need_file(s)
    path = _active(s)
    words = split_words(rest.split(",")[0])
    if len(words) != 2:
        raise StataError(198, "invalid syntax")
    old, new = words
    if old not in data:
        raise StataError(111, f"irfname {old} not found")   # VERIFY
    if new in data:
        raise StataError(110, f"irfname {new} already exists")   # VERIFY
    data = {(new if k == old else k): v for k, v in data.items()}
    _save(path, data)
    s.output.write(f"({len(data[new]['meta']['vars']) ** 2 * (data[new]['meta']['step'] + 1)} "
                   f"real changes made)\n", "text")   # VERIFY
    s.output.write(f"{old} renamed to {new}\n", "text")
