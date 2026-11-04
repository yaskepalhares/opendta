"""The ml command and its programming tools: ml model/init/search/maximize/display/score/clear
and mleval, mlsum, mlvecsum, mlmatsum ([R] ml).

Evaluator methods: lf, lf0 (log likelihood per observation), d0, d1 and d2 (total log
likelihood, with the gradient and Hessian the program gives, or numerical ones). The
maximization is estimation/ml.py. Numerical derivatives are taken with respect to the linear
predictors (lf) or the coefficients (d0, d1), with a fourth-order stencil.

Stata's search, rescale and display details are only partly documented: lines marked VERIFY
come from the layout of the manual's examples (compat/do/1903)."""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.stats import chi2_tail
from ..lang.syntax import find_top, match_options, parse_standard
from .ml import maximize, options as ml_options, _g10
from .results import CoefRow, Estimates, coef_table, g9, post
from .mlmodel import Eq, Model, _LFAMILY, _full, _fun, _mats

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS
_METHODS = ("lf", "d0", "d1", "d2", "lf0", "lf1", "lf2", "gf0", "gf1", "gf2")


# ---------------------------------------------------------------------------
# parsing of ml model
# ---------------------------------------------------------------------------

def _balanced(text: str, i: int) -> int:
    depth = 0
    for j in range(i, len(text)):
        if text[j] == "(":
            depth += 1
        elif text[j] == ")":
            depth -= 1
            if depth == 0:
                return j
    raise StataError(198, "unbalanced parentheses")


def _parse_eq(s: "Session", txt: str, index: int) -> Eq:
    if txt.startswith("/"):
        return Eq(txt[1:] or f"eq{index}", [], [], True, scalar=True)
    inner = txt[1:-1].strip()
    name = ""
    m = re.match(r"^([A-Za-z_]\w*)\s*:(?!:)\s*(.*)$", inner, re.S)
    if m:
        name, inner = m.group(1), m.group(2)
    c = find_top(inner, ",")
    eopt = ""
    if c != -1:
        inner, eopt = inner[:c], inner[c + 1:]
    eo = match_options(eopt, {"noconstant": 3, "offset": 6, "exposure": 8}) if eopt.strip() else {}
    dep, eq_sign, ind = inner.partition("=")
    if not eq_sign:
        dep, ind = "", inner
    dv = s.expand_varlist(dep) if dep.strip() else []
    iv = s.expand_varlist(ind) if ind.strip() else []
    return Eq(name or f"eq{index}", dv, iv, not eo.get("noconstant"),
              eo["offset"] if isinstance(eo.get("offset"), str) else "",
              eo["exposure"] if isinstance(eo.get("exposure"), str) else "")


_OPT = {"vce": 3, "robust": 3, "cluster": 2, "title": 5, "maximize": 4, "difficult": 4, "technique": 4,
        "iterate": 2, "tolerance": 3, "ltolerance": 3, "nrtolerance": 3, "nonrtolerance": 5, "trace": 5,
        "gradient": 4, "showstep": 5, "hessian": 4, "nolog": 3, "log": 3, "showtolerance": 8,
        "lf0": 3, "continue": 4, "waldtest": 4, "nopreserve": 5, "collinear": 4, "missing": 4,
        "init": 4, "search": 3, "level": 3, "constraints": 3, "svy": 3, "subpop": 3, "diparm": 3,
        "obs": 3, "nocnsnotes": 5, "from": 4, "nowarning": 5, "noclear": 5, "nooutput": 5, "noheader": 6,
        "title2": 6, "crittype": 5, "noconstant": 6, "copy": 4, "skip": 4}


def _eqs_tail(s: "Session", text: str):
    """Splits `eq eq ... [weight] [if] [in]` into the equations and the rest."""
    eqs, i, n = [], 0, len(text)
    while i < n:
        if text[i].isspace():
            i += 1
        elif text[i] == "(":
            j = _balanced(text, i)
            eqs.append(_parse_eq(s, text[i:j + 1], len(eqs) + 1))
            i = j + 1
        elif text[i] == "/":
            m = re.match(r"/(\w*)", text[i:])
            eqs.append(_parse_eq(s, m.group(0), len(eqs) + 1))
            i += len(m.group(0))
        else:
            break
    return eqs, text[i:]


_SUBS = (("model", "mod"), ("init", "ini"), ("search", "sea"), ("maximize", "max"), ("display", "di"),
         ("clear", "cl"), ("score", "sc"))
_LATER = ("check", "query", "report", "graph", "plot", "trace", "footnote", "count")


def cmd_ml(s: "Session", args: str) -> None:
    mm = re.match(r"\s*(\w+)(.*)$", args, re.S)
    if not mm:
        raise StataError(198, "invalid syntax")
    sub, rest = mm.group(1).lower(), mm.group(2).strip()
    fns = {"model": _model, "init": _init, "search": _search, "maximize": _maximize,
           "display": _display_cmd, "clear": _clear, "score": _score}
    for full, mn in _SUBS:
        if len(sub) >= len(mn) and full.startswith(sub):
            return fns[full](s, rest)
    for full in _LATER:
        if full.startswith(sub) and len(sub) >= 2:
            raise StataError(198, f"ml {full} is not available yet in OpenDTA")
    raise StataError(199, f"unrecognized command: ml {sub}")


def _active(s: "Session") -> Model:
    m = getattr(s, "_ml", None)
    if m is None:
        raise StataError(301, "last estimates not found")   # VERIFY
    return m


def _model(s: "Session", args: str) -> None:
    c = find_top(args, ",")
    head, opt_text = (args, "") if c == -1 else (args[:c], args[c + 1:])
    words = head.split(None, 2)
    if len(words) < 3:
        raise StataError(198, "invalid syntax")
    method, prog, tail = words[0].lower(), words[1], words[2]
    if method not in _METHODS:
        raise StataError(198, f"unknown method {method}")
    if prog not in s.programs:
        raise StataError(111, f"program {prog} not found")   # VERIFY
    eqs, rest = _eqs_tail(s, tail)
    if not eqs:
        raise StataError(198, "equations required")
    p = parse_standard(rest + (", " + opt_text if opt_text.strip() else ""))
    o = match_options(p.options, _OPT) if p.options.strip() else {}
    from ..commands._util import touse
    from ..commands.summarize import weights
    mask = touse(s, p)
    w, wtype, mask = weights(s, p, mask, ("fweight", "aweight", "pweight", "iweight"))
    design = None
    if o.get("svy"):
        # svy: the design weights (pweights) and, after the fit, the linearized design variance
        from .svy import DesignData, require_design
        if p.weight is not None:
            raise StataError(101, "weights not allowed with the svy option")      # VERIFY
        design = DesignData(s, require_design(s), np.ones(s.data.nobs, bool))
        mask = mask & design.mask
        w = np.zeros(s.data.nobs)
        w[design.mask] = design.w
        wtype = "pweight"
    ds = s.data
    # observations with missing values in any variable of the equations are out of the sample
    used = [v for e in eqs for v in e.depvars + e.indep]
    used += [e.offset for e in eqs if e.offset] + [e.exposure for e in eqs if e.exposure]
    cl = None
    clustvar = ""
    vce = o.get("vce")
    if isinstance(o.get("cluster"), str):
        clustvar = o["cluster"]
    elif isinstance(vce, str) and vce.split()[0].lower().startswith("cl"):
        clustvar = vce.split()[1]
    if clustvar:
        used.append(clustvar)
    if not o.get("missing"):
        for v in used:
            var = ds.get(v)
            mask = mask & (np.array([x != "" for x in var.raw]) if var.is_string else var.data < SYS)
    m = Model(method, prog, eqs, mask, w, wtype, o, f"ml model {args}")
    m.design = design
    m.N = int(mask.sum())
    if m.N == 0:
        raise StataError(2000, "no observations")
    if clustvar:
        cv = ds.get(clustvar)
        m.cluster = np.asarray(cv.raw if cv.is_string else cv.data)[mask]
        m.clustvar = clustvar
    pos = 0
    for e in eqs:
        cols = [ds.get(v).data[mask] for v in e.indep]
        if e.cons:
            cols.append(np.ones(m.N))
        e.X = np.column_stack(cols) if cols else np.zeros((m.N, 0))
        off = np.zeros(m.N)
        if e.offset:
            off = off + ds.get(e.offset).data[mask]
        if e.exposure:
            off = off + np.log(ds.get(e.exposure).data[mask])
        e.off = off
        e.start = pos
        pos += e.k
        for v in e.indep:
            m.names.append((e.name, v))
        if e.cons:
            m.names.append((e.name, "_cons"))
    m.ysets = [v for e in eqs for v in e.depvars]
    if m.w is not None:
        m.w = w[mask]
    m.cns = _constraint_rows(s, m, o.get("constraints")) if o.get("constraints") not in (None, True) else None
    s._ml = m
    if isinstance(o.get("from"), str) or isinstance(o.get("init"), str):
        _init(s, o.get("from") or o.get("init"))
    if o.get("maximize"):
        _maximize(s, "")


# ---------------------------------------------------------------------------
# initial values
# ---------------------------------------------------------------------------

def _init(s: "Session", args: str) -> None:
    m = _active(s)
    c = find_top(args, ",")
    body, opt = (args, "") if c == -1 else (args[:c], args[c + 1:])
    io = match_options(opt, {"copy": 3, "skip": 3}) if opt.strip() else {}
    b = np.zeros(m.nb) if m.b0 is None else m.b0.copy()
    toks = body.split()
    if len(toks) == 1 and re.fullmatch(r"[A-Za-z_]\w*", toks[0]) and toks[0] in _mats(s):
        mat = _mats(s)[toks[0]]
        vals = mat.data.ravel()
        if io.get("skip"):
            names = {(q, c2): i for i, (q, c2) in enumerate(m.names)}
            for j, cn in enumerate(mat.colnames):
                ce = (getattr(mat, "coleq", None) or ["_"] * mat.cols)[j]
                if (ce, cn) in names:
                    b[names[(ce, cn)]] = vals[j]
        else:
            if len(vals) != m.nb:
                raise StataError(503, "conformability error")
            b = vals.astype(float).copy()
    elif io.get("copy") or all(re.fullmatch(r"[-+0-9.eE]+", t) for t in toks):
        vals = [float(t) for t in toks]
        if len(vals) != m.nb:
            raise StataError(503, "conformability error")
        b = np.array(vals)
    else:
        for t in toks:
            mm = re.fullmatch(r"(?:/(\w+)|(?:(\w+):)?(\w+))=(.+)", t)
            if not mm:
                raise StataError(198, "invalid syntax")
            val = float(s.eval(mm.group(4)))
            if mm.group(1):
                eqn, col = mm.group(1), "_cons"
            else:
                eqn, col = mm.group(2) or m.eqs[0].name, mm.group(3)
            try:
                b[m.names.index((eqn, col))] = val
            except ValueError:
                raise StataError(111, f"{eqn}:{col} not found") from None   # VERIFY
    m.b0 = b


# ---------------------------------------------------------------------------
# evaluation
# ---------------------------------------------------------------------------


def _constraint_rows(s: "Session", m: Model, numlist) -> np.ndarray:
    """The rows (R, r) of the constraints numbered in constraints(), read against the
    coefficients of the model."""
    from types import SimpleNamespace
    from ..core.constraints import constraint_transform
    from .postest import _constraint
    from .sysmodels import _constraints, _numlist
    defined = _constraints(s)
    fake = SimpleNamespace(names=[c for _, c in m.names], eqnames=[q for q, _ in m.names])
    rows = []
    for n in _numlist(str(numlist)):
        if n not in defined:
            raise StataError(111, f"constraint {n} not found")
        coef, rhs = _constraint(s, fake, defined[n])
        row = np.zeros(m.nb + 1)
        for j, v in coef.items():
            row[j] += v
        row[-1] = rhs
        rows.append(row)
    C = np.array(rows).reshape(-1, m.nb + 1)
    m.cns_T, m.cns_a = constraint_transform(C)
    m.cns_a = m.cns_a.ravel()
    return C


# ---------------------------------------------------------------------------
# search, rescale and maximize
# ---------------------------------------------------------------------------

def _line(out, lead: str, ll: float, note: str = "", label: str = "log likelihood") -> None:
    shown = "-<inf>" if ll == -np.inf or not np.isfinite(ll) else _g10(ll)
    out.write(f"{lead:<14} {label} = {shown:>10}{note}\n", "text")


def _pseudo(m: Model, mo: dict) -> bool:
    """The log is of a pseudolikelihood with robust or cluster variances, pweights and svy
    (compat 1913)."""
    vce = str(mo["vce"]).split()[0].lower() if isinstance(mo.get("vce"), str) else ""
    return (bool(mo.get("robust")) or vce.startswith(("rob", "cl")) or bool(m.clustvar)
            or m.wtype == "pweight" or getattr(m, "design", None) is not None)


def _start(s: "Session", m: Model, fun, out, quiet: bool, label: str = "log likelihood") -> np.ndarray:
    """Initial point: the given one, else zeros; if the evaluator fails there, a search for a
    feasible point (VERIFY: Stata's random search is undocumented), then rescaling."""
    b = np.zeros(m.nb) if m.b0 is None else m.b0.copy()
    ll = fun(b, 0)[0]
    if not quiet:
        _line(out, "initial:", ll, "  (could not be evaluated)" if not np.isfinite(ll) else "", label)
    if not np.isfinite(ll):
        rng = np.random.default_rng(1)
        best = None
        for trial in range(400):
            cand = b.copy()
            scale = 10.0 ** (trial % 4 - 1)
            for e in m.eqs:
                if e.cons:
                    cand[e.start + e.k - 1] = scale * (1 + trial // 4 * 0.1) if trial < 40 else rng.normal() * scale
                if trial >= 40:
                    cand[e.start:e.start + e.k - (1 if e.cons else 0)] = rng.normal(size=e.k - (1 if e.cons else 0)) * scale
            v = fun(cand, 0)[0]
            if np.isfinite(v):
                best = (cand, v)
                break
        if best is None:
            raise StataError(491, "could not find feasible values")
        b, ll = best
        if not quiet:
            _line(out, "feasible:", ll, label=label)
    elif m.b0 is None and any(e.cons for e in m.eqs):
        # alternative: constants set to -.5 or .5, the better of the two, kept even when it is
        # worse than the initial point (compat 1913: with pweights the alternative -279.10 follows
        # the initial -277.95, and rescale starts from it; VERIFY: Stata's rule is not documented)
        best = None
        for c in (-0.5, 0.5):
            cand = b.copy()
            for e in m.eqs:
                if e.cons:
                    cand[e.start + e.k - 1] = c
            v = fun(cand, 0)[0]
            if np.isfinite(v) and (best is None or v > best[1]):
                best = (cand, v)
        if best is not None:
            b, ll = best
            if not quiet:
                _line(out, "alternative:", ll, label=label)
    # rescale: best multiple of the whole vector, then of each equation
    def rescale(idx_sets):
        nonlocal b, ll
        for idx in idx_sets:
            bestc, bestv = 1.0, ll
            for c in (0.0625, 0.125, 0.25, 0.5, 0.75, 1.5, 2.0, 4.0, 8.0, 16.0):
                cand = b.copy()
                cand[idx] = b[idx] * c
                v = fun(cand, 0)[0]
                if np.isfinite(v) and v > bestv + 1e-12 * abs(bestv):
                    bestc, bestv = c, v
            if bestc != 1.0:
                b = b.copy()
                b[idx] = b[idx] * bestc
                ll = bestv
    rescale([np.arange(m.nb)])
    if not quiet:
        _line(out, "rescale:", ll, label=label)
    if m.K > 1:
        rescale([np.arange(e.start, e.start + e.k) for e in m.eqs])
        if not quiet:
            _line(out, "rescale eq:", ll, label=label)
    return b


def _cleanup(s: "Session", m: Model) -> None:
    """Drops the temporary variables of the evaluator."""
    names = [n for n in list(getattr(m, "_tmp", None) or []) + [m.sampvar, m.wvar] if n and s.data.has(n)]
    if names:
        changed = s.data.changed
        s.data.drop_vars(names)
        s.data.changed = changed
    m._tmp = None
    m.sampvar = m.wvar = ""
    s._ml_active = None


def _search(s: "Session", args: str) -> None:
    m = _active(s)
    try:
        fun = _fun(s, m)
        m.b0 = _start(s, m, fun, s.output, False)
    finally:
        _cleanup(s, m)


def _maximize(s: "Session", args: str) -> None:
    m = _active(s)
    out = s.output
    o = match_options(args.lstrip(", "), _OPT) if args.strip() else {}
    mo = dict(m.opts)
    mo.update(o)
    nolog = bool(mo.get("nolog"))
    label = "log pseudolikelihood" if _pseudo(m, mo) else "log likelihood"
    try:
        fun = _fun(s, m)
        b = _start(s, m, fun, out, nolog, label)
        if getattr(m, "cns", None) is not None:
            # constraints(): b = T bc + a, maximized over bc
            T, a = m.cns_T, m.cns_a

            def cfun(bc, todo):
                ll, g, H = fun(T @ bc + a, todo)
                if g is None:
                    return ll, g, H
                return ll, T.T @ g, (T.T @ H @ T if H is not None else None)
            r = maximize(cfun, T.T @ (b - a), **ml_options(mo))
            r.b = T @ r.b + a
        else:
            r = maximize(fun, b, **ml_options(mo))
        if not nolog:
            for it, ll, note in r.log:
                lead = f"Iteration {it}:"
                out.write(f"{lead:<14} {label} = {_g10(ll):>10}{note}  \n", "text")
            if not r.converged:
                out.write("convergence not achieved\n", "error")
        m.result = r
        _variance(s, m, fun, mo)
    finally:
        _cleanup(s, m)
    _post(s, m, mo)
    if not mo.get("nooutput"):
        _display(s, m, mo)


def _variance(s: "Session", m: Model, fun, mo: dict) -> None:
    r = m.result
    ll, g, H = fun(r.b, 2)
    r.g = g
    try:
        if getattr(m, "cns", None) is not None:
            T = m.cns_T
            Vh = T @ np.linalg.inv(-(T.T @ H @ T)) @ T.T
        else:
            Vh = np.linalg.inv(-H)
    except np.linalg.LinAlgError:
        Vh = np.full_like(H, np.nan)
    vce = str(mo["vce"]).split()[0].lower() if isinstance(mo.get("vce"), str) else "oim"
    robust = bool(mo.get("robust")) or vce.startswith("rob") or bool(m.clustvar) or m.wtype == "pweight"
    m.vcetype, m.vce = "OIM", "oim"
    m.Vh = Vh
    if vce.startswith("opg"):
        robust = False
    if getattr(m, "design", None) is not None:
        # svy: linearized design variance of the total scores (VERIFY: replication methods)
        dd = m.design
        S = _scores(m)
        Z = np.zeros((dd.n, m.nb))
        pos = np.full(len(m.touse), -1)
        pos[np.flatnonzero(dd.mask)] = np.arange(dd.n)
        Z[pos[m.touse]] = S
        m.V = Vh @ dd.variance(Z) @ Vh
        # for e(V_srs): the linearized scores and the weights of the sample on the design rows
        wd = np.zeros(dd.n)
        wd[pos[m.touse]] = dd.w[pos[m.touse]]
        from .svy import srs_variances
        m._srs = srs_variances(Z @ Vh, wd, dd, int(m.touse.sum()))
        m.vcetype, m.vce = "Linearized", "linearized"
        return
    if robust or vce.startswith("opg"):
        if m.method not in ("lf", "mlexp") + _LFAMILY:
            raise StataError(198, f"vce({vce}) is not allowed with method {m.method}")   # VERIFY
        S = _scores(m)
        if vce.startswith("opg"):
            m.V = np.linalg.inv(S.T @ S)
            m.vcetype, m.vce = "OPG", "opg"
        else:
            if m.cluster is not None:
                ids = np.unique(m.cluster, return_inverse=True)[1]
                S = np.add.reduceat(S[np.argsort(ids, kind="stable")],
                                    np.r_[0, np.cumsum(np.bincount(ids))[:-1]], axis=0)
                gcount = len(np.unique(ids))
                adj = gcount / (gcount - 1) if gcount > 1 else 1.0
                m.nclust = gcount
            else:
                adj = m.N / (m.N - 1) if m.N > 1 else 1.0
            m.V = adj * (Vh @ (S.T @ S) @ Vh)
            m.vcetype, m.vce = "Robust", "robust"      # VERIFY: small-sample adjustment
    else:
        m.V = Vh


def _scores(m: Model) -> np.ndarray:
    """Weighted per-observation scores with respect to the coefficients (n x nb)."""
    if m.method in _LFAMILY:
        return m._S
    if m.method not in ("lf", "mlexp"):
        raise StataError(198, f"method {m.method} does not give the scores this variance needs")   # VERIFY
    w = m.w if m.w is not None else np.ones(m.N)
    G = getattr(m, "_G", None)
    return np.column_stack([e.X * (w * G[:, j])[:, None] for j, e in enumerate(m.eqs)])


def _post(s: "Session", m: Model, mo: dict) -> None:
    from ..commands.matrix import Matrix
    r = m.result
    b, V = r.b, m.V
    names = [c for _, c in m.names]
    eqs = [q for q, _ in m.names]
    mlx = m.method == "mlexp"
    est = Estimates(cmd="mlexp" if mlx else "ml", depvar=m.ysets[0] if m.ysets else "", names=names, b=b, V=V,
                    N=m.N, stat="z")
    est.eqnames = eqs
    est.vcetype = m.vcetype
    est.title = mo["title"].strip().strip('"') if isinstance(mo.get("title"), str) else \
        ("Maximum likelihood estimation" if mlx else "")
    est.level = float(mo["level"]) if isinstance(mo.get("level"), str) else float(s.settings.get("level", 95))
    k_eq = m.K
    skip = isinstance(mo.get("waldtest"), str) and float(mo["waldtest"]) == -1
    # Wald test of the first equation without its constant (VERIFY: waldtest(#) with # > 0)
    keep = [] if mlx else [i for i, (q, c) in enumerate(m.names) if q == m.eqs[0].name and c != "_cons"]
    dfm = len(keep)
    chi2 = SYS
    if dfm and not skip and getattr(m, "cns", None) is not None:
        # with constraints the test has the rank of the constrained variance (VERIFY)
        Vk = V[np.ix_(keep, keep)]
        dfm = int(np.linalg.matrix_rank(Vk))
        chi2 = float(b[keep] @ np.linalg.pinv(Vk) @ b[keep]) if dfm else SYS
    elif dfm and not skip:
        try:
            chi2 = float(b[keep] @ np.linalg.solve(V[np.ix_(keep, keep)], b[keep]))
        except np.linalg.LinAlgError:
            chi2 = SYS
    m.chi2, m.dfm = chi2, dfm
    scalars = [("N", float(m.N)), ("k", float(m.nb)), ("k_eq", float(k_eq)), ("k_dv", float(len(m.ysets))),
               ("converged", 1.0 if r.converged else 0.0), ("rc", 0.0 if r.converged else 430.0),
               ("ic", float(r.iterations)), ("rank", float(np.linalg.matrix_rank(V))), ("ll", r.ll)]
    if chi2 < SYS:
        scalars += [("df_m", float(dfm)), ("chi2", chi2), ("p", float(chi2_tail(chi2, dfm)))]
    if not mlx:
        # order of the reference log: rank N ic k k_eq k_dv converged rc ll k_eq_model df_m chi2 p
        # (VERIFY: e(k_eq_model) when there is no Wald test)
        sc = dict(scalars)
        sc["k_eq_model"] = float(k_eq)
        scalars = [(k, sc[k]) for k in ("p", "chi2", "df_m", "k_eq_model", "ll", "rc", "converged", "k_dv",
                                        "k_eq", "k", "ic", "N", "rank") if k in sc]
    if m.clustvar:
        scalars.append(("N_clust", float(getattr(m, "nclust", 0))))
    macros = [("properties", "b V"), ("depvar", " ".join(m.ysets)), ("which", "max"), ("technique", "nr"),
              ("ml_method", m.method), ("user", m.prog), ("vce", m.vce),
              ("vcetype", m.vcetype if m.vcetype != "OIM" else ""), ("predict", "ml_p"), ("opt", "moptimize"),
              ("chi2type", "Wald" if chi2 < SYS else ""), ("cmd", "ml")]
    if mlx:
        macros = [("cmd", "mlexp"), ("cmdline", m.cmdline), ("lexp", m.lexp), ("title", est.title),
                  ("technique", "nr"), ("vce", m.vce), ("vcetype", m.vcetype), ("opt", "moptimize"),
                  ("which", "max"), ("crittype", "log likelihood"), ("properties", "b V")]
    if m.clustvar:
        macros.append(("clustvar", m.clustvar))
    if m.wtype:
        macros.append(("wtype", m.wtype))
    ilog = np.zeros((1, 20))
    vals = [ll for it, ll, _ in r.log if it >= 0][:20]
    ilog[0, :len(vals)] = vals
    extra = [("ilog", Matrix(ilog, ["r1"], [str(i) for i in range(1, 21)])),
             ("gradient", Matrix(np.asarray(r.g).reshape(1, -1), ["r1"], names))]
    if m.vcetype in ("Robust", "Linearized"):
        extra.append(("V_modelbased", Matrix(m.Vh, names, names)))
    if getattr(m, "cns", None) is not None:
        cm = Matrix(m.cns.copy(), [f"r{i}" for i in range(1, len(m.cns) + 1)], names + ["_r"])
        # the last column is _Cns:_r, and e(Cns) comes right after e(V) (compat 1913)
        cm.coleq = eqs + ["_Cns" if any(q not in ("", "_") for q in eqs) else "_"]
        extra.insert(0, ("Cns", cm))
    dd = getattr(m, "design", None)
    if dd is not None:
        # ml, svy: the results in the order of Stata's ereturn list, with the design-based F test
        # of the first equation and no e(ll) nor e(rank) (compat 1913)
        from .svyreg import _adjusted_F
        d = dd.d
        F, df1, df2 = _adjusted_F(r.b, V, keep, dd.df) if keep else (SYS, 0, 0)
        m.svyF = (F, df1, df2)
        scalars = [("ic", float(r.iterations)), ("k", float(m.nb)), ("k_eq", float(k_eq)),
                   ("k_dv", float(len(m.ysets))), ("converged", 1.0 if r.converged else 0.0),
                   ("rc", 0.0 if r.converged else 430.0), ("N_psu", float(dd.n_psu)),
                   ("N_strata", float(dd.n_strata)), ("df_r", float(dd.df)), ("N", float(m.N)),
                   ("N_pop", float(m.w.sum())), ("census", 0.0), ("singleton", 0.0), ("N_strata_omit", 0.0),
                   ("stages", 1.0), ("k_eq_model", float(k_eq))]
        if F < SYS:
            scalars += [("df_m", float(df1)), ("F", F), ("p", float(st_f_sf(F, df1, df2)))]
        macros = [("cmd", "ml"), ("svyml", "svyml"), ("opt", "moptimize"), ("predict", "ml_p"),
                  ("singleunit", d.singleunit)]
        if d.strata:
            macros.append(("strata1", d.strata))
        if d.psu:
            macros.append(("su1", d.psu))
        wv = getattr(m, "wvar_posted", "")
        if wv:
            # Stata's is a temporary variable (__000004 in compat 1913)
            macros += [("wvar", wv), ("wexp", f"= {wv}")]
        macros += [("wtype", "pweight"), ("vcetype", m.vcetype), ("vce", m.vce), ("prefix", "svy"),
                   ("user", m.prog), ("ml_method", m.method), ("technique", "nr"), ("which", "max"),
                   ("depvar", " ".join(m.ysets)), ("properties", "b V")]
        Swor, Swr = m._srs
        extra += [("V_srssubwr", Matrix(Swr, names, names)), ("V_srssub", Matrix(Swor, names, names)),
                  ("V_srswr", Matrix(Swr, names, names)), ("V_srs", Matrix(Swor, names, names)),
                  ("_N_strata_certain", Matrix(np.zeros((1, 1)))), ("_N_strata_single", Matrix(np.zeros((1, 1)))),
                  ("_N_strata", Matrix(np.array([[float(dd.n_strata)]])))]
        post(s, est, scalars, macros, m.touse.copy(), extra_mats=extra)
        s._estimates = est
        return
    sample = m.touse.copy()
    post(s, est, list(reversed(scalars)), list(reversed(macros)), sample, extra_mats=extra)
    s._estimates = est


def st_f_sf(F: float, d1: float, d2: float) -> float:
    from scipy.stats import f
    return float(f.sf(F, d1, d2))


def _display(s: "Session", m: Model, mo: dict) -> None:
    out = s.output
    est = s._estimates
    level = float(mo["level"]) if isinstance(mo.get("level"), str) else est.level
    lll = "Log pseudolikelihood" if m.vcetype == "Robust" else "Log likelihood"
    if not mo.get("noheader") and m.method == "mlexp":
        out.write(f"\n{est.title}\n\n", "text")
        out.write(f"{lll + ' = ' + _g10(m.result.ll):<48}{'Number of obs':<18}= ", "text")
        out.write(f"{m.N:>10,}\n\n", "result")
    elif not mo.get("noheader") and getattr(m, "design", None) is not None:
        # svy: the design header and the design-based F test of the first equation (VERIFY)
        from .svy import format_size
        dd = m.design
        F, df1, df2 = m.svyF
        out.write(f"\n{est.title}\n\n" if est.title else "\n", "text")
        # the layout of svy: regress (compat 1913)
        rows = [(f"Number of strata   = {dd.n_strata:>9,}", "Number of obs", f"{m.N:,}"),
                (f"Number of PSUs     = {dd.n_psu:>9,}", "Population size", format_size(float(m.w.sum()))),
                ("", "Design df", f"{dd.df:,}")]
        if F < SYS:
            rows += [("", f"F({df1:>4},{df2:>7})", g9(F, "%10.2f")),
                     ("", "Prob > F", g9(float(st_f_sf(F, df1, df2)), "%10.4f"))]
        for left, lab, val in rows:
            out.write(f"{left:<48}{lab:<18}= ", "text")
            out.write(f"{val:>10}\n", "result")
        out.write("\n", "text")
    elif not mo.get("noheader"):
        out.write("\n", "text")
        rows = [(est.title, "Number of obs", f"{m.N:,}")]
        if m.chi2 < SYS:
            rows.append(("", f"Wald chi2({m.dfm})", g9(m.chi2, "%10.2f")))
            rows.append((f"{lll} = {_g10(m.result.ll):>10}", "Prob > chi2",
                         g9(float(chi2_tail(m.chi2, m.dfm)), "%10.4f")))
        else:
            rows.append((f"{lll} = {_g10(m.result.ll):>10}", "", ""))
        for left, lab, val in rows:
            out.write(f"{left:<48}{lab:<18}{'= ' if lab else ''}", "text")
            out.write(f"{val:>10}\n" if lab else "\n", "result")
        out.write("\n", "text")
    V = m.V
    se = np.sqrt(np.maximum(np.diag(V), 0))
    rows_ = []
    multi = m.K > 1
    for j, e in enumerate(m.eqs):
        if e.scalar and m.method == "mlexp":
            rows_.append(CoefRow("coef", "/" + e.name, float(m.result.b[e.start]), float(se[e.start]), name=e.name))
            continue
        if multi:
            if rows_ and j:
                rows_.append(CoefRow("sep"))
            rows_.append(CoefRow("eq", e.name))
        for i in range(e.k):
            idx = e.start + i
            nm = m.names[idx][1]
            rows_.append(CoefRow("coef", nm, float(m.result.b[idx]), float(se[idx]), name=nm))
    if m.vcetype == "Robust" and m.clustvar:
        out.write(f"{'':>48}(Std. Err. adjusted for {getattr(m, 'nclust', 0)} clusters in {m.clustvar})\n", "text")
    if getattr(m, "cns", None) is not None:
        from .postest import _ctext
        from types import SimpleNamespace
        fake = SimpleNamespace(names=[c for _, c in m.names], eqnames=[q for q, _ in m.names])
        for i, row in enumerate(m.cns, 1):
            coef = {j: float(v) for j, v in enumerate(row[:-1]) if v != 0}
            out.write(f" ({i:>2})  {_ctext(fake, coef, float(row[-1]))}\n", "text")
    dep = m.ysets[0] if m.ysets else ""
    dd = getattr(m, "design", None)
    coef_table(s, "" if m.method == "mlexp" else dep, rows_, stat="t" if dd is not None else "z",
               df=dd.df if dd is not None else None, level=level,
               vcetype="Robust" if m.vcetype == "Robust" else "" if m.vcetype in ("OIM",) else m.vcetype)


def _display_cmd(s: "Session", args: str) -> None:
    m = _active(s)
    if m.result is None:
        raise StataError(301, "last estimates not found")
    o = match_options(args.lstrip(", "), _OPT) if args.strip() else {}
    _display(s, m, o)


def _clear(s: "Session", args: str) -> None:
    m = getattr(s, "_ml", None)
    if m is not None and s.data.has(m.sampvar):
        pass
    s._ml = None
    for k in [k for k in s.macros.globals if k.startswith("ML_")]:
        del s.macros.globals[k]


def _score(s: "Session", args: str) -> None:
    m = _active(s)
    if m.method != "lf":
        raise StataError(198, "ml score requires an lf-type evaluator")   # VERIFY
    from ..core.dataset import Variable
    names = args.replace(",", " ").split()
    if len(names) != m.K:
        raise StataError(198, "one new variable per equation required")
    try:
        fun = _fun(s, m)
        fun(m.result.b, 2)
    finally:
        _cleanup(s, m)
    for j, nm in enumerate(names):
        s.data.add(Variable(nm, "double", _full(m, m._G[:, j], s.data.nobs)))
