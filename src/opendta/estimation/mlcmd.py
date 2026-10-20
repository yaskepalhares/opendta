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
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.stats import chi2_tail
from ..lang.syntax import find_top, match_options, parse_standard
from .ml import MLResult, maximize, options as ml_options, _g10
from .results import CoefRow, Estimates, coef_table, g9, post

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS
_METHODS = ("lf", "d0", "d1", "d2")
_UNSUPPORTED = ("lf0", "lf1", "lf2", "gf0", "gf1", "gf2")


@dataclass
class Eq:
    name: str
    depvars: list[str]
    indep: list[str]
    cons: bool = True
    offset: str = ""
    exposure: str = ""
    scalar: bool = False                # /name: constant only, a parameter of the likelihood
    X: np.ndarray | None = None         # design matrix over the sample (constant last)
    off: np.ndarray | None = None
    start: int = 0

    @property
    def k(self) -> int:
        return len(self.indep) + (1 if self.cons else 0)


@dataclass
class Model:
    method: str
    prog: str
    eqs: list[Eq]
    touse: np.ndarray
    w: np.ndarray | None
    wtype: str
    opts: dict
    cmdline: str
    cluster: np.ndarray | None = None
    clustvar: str = ""
    b0: np.ndarray | None = None
    result: MLResult | None = None
    V: np.ndarray | None = None
    N: int = 0
    names: list[tuple[str, str]] = field(default_factory=list)    # (equation, column)
    ysets: list[str] = field(default_factory=list)
    sampvar: str = ""
    wvar: str = ""

    @property
    def K(self) -> int:
        return len(self.eqs)

    @property
    def nb(self) -> int:
        return sum(e.k for e in self.eqs)


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
    if method in _UNSUPPORTED:
        raise StataError(198, f"method {method} is not available yet in OpenDTA")
    if method not in _METHODS:
        raise StataError(198, f"unknown method {method}")
    from ..lang.programs import Program      # noqa: F401
    if prog not in s.programs:
        raise StataError(111, f"program {prog} not found")   # VERIFY
    eqs, rest = _eqs_tail(s, tail)
    if not eqs:
        raise StataError(198, "equations required")
    p = parse_standard(rest + (", " + opt_text if opt_text.strip() else ""))
    o = match_options(p.options, _OPT) if p.options.strip() else {}
    if o.get("constraints") or o.get("svy"):
        raise StataError(198, "option constraints() and svy are not available yet in OpenDTA")
    from ..commands._util import touse
    from ..commands.summarize import weights
    mask = touse(s, p)
    w, wtype, mask = weights(s, p, mask, ("fweight", "aweight", "pweight", "iweight"))
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


def _mats(s: "Session") -> dict:
    from ..commands.matrix import _store
    return _store(s)


# ---------------------------------------------------------------------------
# evaluation
# ---------------------------------------------------------------------------

def _theta(m: Model, b: np.ndarray) -> np.ndarray:
    th = np.empty((m.N, m.K))
    for j, e in enumerate(m.eqs):
        th[:, j] = e.X @ b[e.start:e.start + e.k] + e.off
    return th


def _setup(s: "Session", m: Model) -> None:
    """Globals and sample variable read by the evaluator programs."""
    from ..commands.syntaxcmd import _temp_name
    from ..core.dataset import Variable
    ds = s.data
    if not m.sampvar or not ds.has(m.sampvar):
        m.sampvar = _temp_name(s)
        ds.add(Variable(m.sampvar, "double", np.zeros(ds.nobs)))
    ds.get(m.sampvar).data = m.touse.astype(float)
    g = s.macros.globals
    g["ML_samp"] = m.sampvar
    for i, y in enumerate(m.ysets, start=1):
        g[f"ML_y{i}"] = y
    if m.w is not None:
        if not m.wvar or not ds.has(m.wvar):
            m.wvar = _temp_name(s)
            ds.add(Variable(m.wvar, "double", np.zeros(ds.nobs)))
        full = np.zeros(ds.nobs)
        full[m.touse] = m.w
        ds.get(m.wvar).data = full
        g["ML_w"] = m.wvar
    else:
        g.pop("ML_w", None)
    s._ml_active = m


def _full(m: Model, vec: np.ndarray, nobs: int) -> np.ndarray:
    out = np.full(nobs, SYS)
    out[m.touse] = vec
    return out


def _call(s: "Session", text: str) -> None:
    with s.output.quietly():
        s.interp._execute_expanded(text)


def _lnf_obs(s: "Session", m: Model, th: np.ndarray) -> np.ndarray:
    """lf: the program fills lnf for the given linear predictors; returns lnf per observation."""
    if m.method == "mlexp":
        return _mlexp_obs(s, m, th)
    from ..commands.syntaxcmd import _temp_name
    from ..core.dataset import Variable
    ds = s.data
    names = getattr(m, "_tmp", None)
    if names is None or any(not ds.has(n) for n in names):
        names = []
        for _ in range(m.K + 1):
            nm = _temp_name(s)
            ds.add(Variable(nm, "double", np.full(ds.nobs, SYS)))
            names.append(nm)
        m._tmp = names
    for j in range(m.K):
        ds.get(names[j + 1]).data = _full(m, th[:, j], ds.nobs)
    ds.get(names[0]).data = np.full(ds.nobs, SYS)
    _call(s, f"{m.prog} " + " ".join(names))
    v = ds.get(names[0]).data[m.touse]
    return np.where(v >= SYS, np.nan, v)


def _total(m: Model, lnf: np.ndarray) -> float:
    if np.any(np.isnan(lnf)):
        return float("nan")
    return float((lnf if m.w is None else lnf * m.w).sum())


def _eval_d(s: "Session", m: Model, b: np.ndarray, todo: int):
    mats = _mats(s)
    from ..commands.matrix import Matrix
    bm = Matrix(b.reshape(1, -1), ["y1"], [c for _, c in m.names])
    bm.coleq = [q for q, _ in m.names]
    mats["__mlb"] = bm
    s.scalars["__mllnf"] = SYS
    extra = ""
    if todo >= 1:
        mats.pop("__mlg", None)
        extra = " __mlg"
    if todo >= 2:
        mats.pop("__mlH", None)
        extra += " __mlH"
    _call(s, f"{m.prog} {todo} __mlb __mllnf{extra}")
    ll = float(s.scalars.get("__mllnf", SYS))
    ll = float("nan") if ll >= SYS else ll
    g = H = None
    if todo >= 1 and "__mlg" in mats:
        g = np.asarray(mats["__mlg"].data, dtype=float).ravel()
    if todo >= 2 and "__mlH" in mats:
        H = np.asarray(mats["__mlH"].data, dtype=float)
    return ll, g, H


def _num_fun(f0, nb: int):
    """fun(b, todo) with fourth-order numerical gradient and Hessian of the scalar f0(b)."""
    def num_grad(b):
        g = np.zeros(nb)
        for j in range(nb):
            hj = 5e-3 * max(1.0, abs(b[j]))
            e = np.zeros(nb)
            e[j] = hj
            g[j] = (-f0(b + 2 * e) + 8 * f0(b + e) - 8 * f0(b - e) + f0(b - 2 * e)) / (12 * hj)
        return g

    def num_hess(b, f00):
        H = np.zeros((nb, nb))
        hs = [5e-3 * max(1.0, abs(b[j])) for j in range(nb)]
        for j in range(nb):
            e = np.zeros(nb)
            e[j] = hs[j]
            H[j, j] = (-f0(b + 2 * e) + 16 * f0(b + e) - 30 * f00 + 16 * f0(b - e) - f0(b - 2 * e)) / (12 * hs[j] ** 2)
            for k in range(j):
                def D(a):
                    ej = np.zeros(nb)
                    ek = np.zeros(nb)
                    ej[j], ek[k] = a * hs[j], a * hs[k]
                    return (f0(b + ej + ek) - f0(b + ej - ek) - f0(b - ej + ek) + f0(b - ej - ek)) / (4 * a * a * hs[j] * hs[k])
                H[j, k] = H[k, j] = (4 * D(1) - D(2)) / 3
        return H

    def fun(b, todo):
        ll = f0(b)
        if todo == 0:
            return (ll if np.isfinite(ll) else -np.inf), None, None
        if not np.isfinite(ll):
            return -np.inf, np.zeros(nb), np.full((nb, nb), np.nan)
        return ll, num_grad(b), num_hess(b, ll)
    return fun


def _d_fun(s: "Session", m: Model):
    """fun(b, todo) -> (ll, g, H) for maximize()."""
    nb = m.nb
    use_g = m.method in ("d1", "d2")
    use_H = m.method == "d2"
    if not use_g:
        return _num_fun(lambda b: _eval_d(s, m, b, 0)[0], nb)

    def fun(b, todo):
        if todo == 0:
            ll = _eval_d(s, m, b, 0)[0]
            return (ll if np.isfinite(ll) else -np.inf), None, None
        ll, g, H = _eval_d(s, m, b, 2 if use_H else 1)
        if g is None:
            raise StataError(198, f"{m.prog} did not return the gradient")   # VERIFY
        if not np.isfinite(ll):
            return -np.inf, g, np.full((nb, nb), np.nan)
        if H is None:
            H = np.zeros((nb, nb))
            for j in range(nb):
                hj = 5e-3 * max(1.0, abs(b[j]))
                e = np.zeros(nb)
                e[j] = hj
                gp, gm = _eval_d(s, m, b + e, 1)[1], _eval_d(s, m, b - e, 1)[1]
                gp2, gm2 = _eval_d(s, m, b + 2 * e, 1)[1], _eval_d(s, m, b - 2 * e, 1)[1]
                H[:, j] = (-gp2 + 8 * gp - 8 * gm + gm2) / (12 * hj)
            H = (H + H.T) / 2
        return ll, g, H
    return fun


def _lf_derivs(s: "Session", m: Model, th: np.ndarray):
    """Per-observation first and second derivatives of lnf with respect to the linear
    predictors (fourth-order central differences)."""
    n, K = th.shape
    h = 5e-3 * np.maximum(1.0, np.abs(th))

    def f(t):
        return _lnf_obs(s, m, t)

    def shift(cols):
        t = th.copy()
        for j, a in cols:
            t[:, j] = t[:, j] + a * h[:, j]
        return t
    f0 = f(th)
    G = np.zeros((n, K))
    H = np.zeros((n, K, K))
    for j in range(K):
        fp1, fm1 = f(shift([(j, 1)])), f(shift([(j, -1)]))
        fp2, fm2 = f(shift([(j, 2)])), f(shift([(j, -2)]))
        G[:, j] = (-fp2 + 8 * fp1 - 8 * fm1 + fm2) / (12 * h[:, j])
        H[:, j, j] = (-fp2 + 16 * fp1 - 30 * f0 + 16 * fm1 - fm2) / (12 * h[:, j] ** 2)
        for k in range(j):
            def D(a):
                return (f(shift([(j, a), (k, a)])) - f(shift([(j, a), (k, -a)])) - f(shift([(j, -a), (k, a)]))
                        + f(shift([(j, -a), (k, -a)]))) / (4 * a * a * h[:, j] * h[:, k])
            H[:, j, k] = H[:, k, j] = (4 * D(1) - D(2)) / 3
    return f0, G, H


def _lf_fun(s: "Session", m: Model):
    w = m.w if m.w is not None else np.ones(m.N)

    def fun(b, todo):
        th = _theta(m, b)
        if todo == 0:
            ll = _total(m, _lnf_obs(s, m, th))
            return (ll if np.isfinite(ll) else -np.inf), None, None
        lnf, G, Hh = _lf_derivs(s, m, th)
        ll = _total(m, lnf)
        if not np.isfinite(ll) or not np.all(np.isfinite(G)):
            return -np.inf, np.zeros(m.nb), np.full((m.nb, m.nb), np.nan)
        g = np.zeros(m.nb)
        H = np.zeros((m.nb, m.nb))
        for j, ej in enumerate(m.eqs):
            sj = slice(ej.start, ej.start + ej.k)
            g[sj] = ej.X.T @ (w * G[:, j])
            for k, ek in enumerate(m.eqs):
                sk = slice(ek.start, ek.start + ek.k)
                H[sj, sk] = ej.X.T @ (ek.X * (w * Hh[:, j, k])[:, None])
        m._G = G
        return ll, g, H
    return fun


def _fun(s: "Session", m: Model):
    _setup(s, m)
    return _lf_fun(s, m) if m.method in ("lf", "mlexp") else _d_fun(s, m)


# ---------------------------------------------------------------------------
# search, rescale and maximize
# ---------------------------------------------------------------------------

def _line(out, lead: str, ll: float, note: str = "") -> None:
    shown = "-<inf>" if ll == -np.inf or not np.isfinite(ll) else _g10(ll)
    out.write(f"{lead:<14} log likelihood = {shown:>10}{note}\n", "text")


def _start(s: "Session", m: Model, fun, out, quiet: bool) -> np.ndarray:
    """Initial point: the given one, else zeros; if the evaluator fails there, a search for a
    feasible point (VERIFY: Stata's random search is undocumented), then rescaling."""
    b = np.zeros(m.nb) if m.b0 is None else m.b0.copy()
    ll = fun(b, 0)[0]
    if not quiet:
        _line(out, "initial:", ll, "  (could not be evaluated)" if not np.isfinite(ll) else "")
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
            _line(out, "feasible:", ll)
    elif m.b0 is None and any(e.cons for e in m.eqs):
        # alternative: constants set to -.5 or .5, kept if the likelihood improves (VERIFY: the
        # reference log shows constant -.5 for a logit with 43% of ones; Stata's rule is not documented)
        best = None
        for c in (-0.5, 0.5):
            cand = b.copy()
            for e in m.eqs:
                if e.cons:
                    cand[e.start + e.k - 1] = c
            v = fun(cand, 0)[0]
            if np.isfinite(v) and v > ll and (best is None or v > best[1]):
                best = (cand, v)
        if best is not None:
            b, ll = best
            if not quiet:
                _line(out, "alternative:", ll)
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
        _line(out, "rescale:", ll)
    if m.K > 1:
        rescale([np.arange(e.start, e.start + e.k) for e in m.eqs])
        if not quiet:
            _line(out, "rescale eq:", ll)
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
    try:
        fun = _fun(s, m)
        b = _start(s, m, fun, out, nolog)
        r = maximize(fun, b, **ml_options(mo))
        if not nolog:
            for it, ll, note in r.log:
                lead = f"Iteration {it}:"
                out.write(f"{lead:<14} log likelihood = {_g10(ll):>10}{note}  \n", "text")
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
    try:
        Vh = np.linalg.inv(-H)
    except np.linalg.LinAlgError:
        Vh = np.full_like(H, np.nan)
    vce = str(mo["vce"]).split()[0].lower() if isinstance(mo.get("vce"), str) else "oim"
    robust = bool(mo.get("robust")) or vce.startswith("rob") or bool(m.clustvar) or m.wtype == "pweight"
    m.vcetype, m.vce = "OIM", "oim"
    m.Vh = Vh
    if vce.startswith("opg"):
        robust = False
    if robust or vce.startswith("opg"):
        if m.method not in ("lf", "mlexp"):
            raise StataError(198, f"vce({vce}) is not allowed with method {m.method}")   # VERIFY
        w = m.w if m.w is not None else np.ones(m.N)
        G = getattr(m, "_G", None)
        S = np.column_stack([e.X * (w * G[:, j])[:, None] for j, e in enumerate(m.eqs)])
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
    if dfm and not skip:
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
    if m.vcetype == "Robust":
        extra.append(("V_modelbased", Matrix(m.Vh, names, names)))
    sample = m.touse.copy()
    post(s, est, list(reversed(scalars)), list(reversed(macros)), sample, extra_mats=extra)
    s._estimates = est


def _display(s: "Session", m: Model, mo: dict) -> None:
    out = s.output
    est = s._estimates
    level = float(mo["level"]) if isinstance(mo.get("level"), str) else est.level
    lll = "Log pseudolikelihood" if m.vcetype == "Robust" else "Log likelihood"
    if not mo.get("noheader") and m.method == "mlexp":
        out.write(f"\n{est.title}\n\n", "text")
        out.write(f"{lll + ' = ' + _g10(m.result.ll):<48}{'Number of obs':<18}= ", "text")
        out.write(f"{m.N:>10,}\n\n", "result")
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
    dep = m.ysets[0] if m.ysets else ""
    coef_table(s, "" if m.method == "mlexp" else dep, rows_, stat="z", df=None, level=level,
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


# ---------------------------------------------------------------------------
# programming tools
# ---------------------------------------------------------------------------

def _tool_split(args: str, opts: dict):
    c = find_top(args, ",")
    body, opt = (args, "") if c == -1 else (args[:c], args[c + 1:])
    return body, (match_options(opt, opts) if opt.strip() else {})


def _eq_index(m: Model, v) -> int:
    k = int(v) if isinstance(v, str) else 1
    if not 1 <= k <= m.K:
        raise StataError(198, f"equation {k} not found")
    return k - 1


def cmd_mleval(s: "Session", args: str) -> None:
    body, o = _tool_split(args, {"eq": 2, "scalar": 6})
    m = getattr(s, "_ml_active", None)
    if m is None:
        raise StataError(198, "mleval can only be used within an ml evaluator")
    mm = re.fullmatch(r"\s*(\S+?)\s*=\s*(\S+)\s*", body)
    if not mm:
        raise StataError(198, "invalid syntax")
    new, bname = mm.groups()
    e = m.eqs[_eq_index(m, o.get("eq"))]
    bvec = np.asarray(_mats(s)[bname].data, dtype=float).ravel()
    coef = bvec[e.start:e.start + e.k]
    if o.get("scalar"):
        s.scalars[new] = float(coef[-1] if e.cons else 0.0)
        return
    from ..core.dataset import Variable
    xb = e.X @ coef + e.off
    full = _full(m, xb, s.data.nobs)
    if s.data.has(new):
        s.data.get(new).data = full
    else:
        s.data.add(Variable(new, "double", full))


def _expr(s: "Session", m: Model, text: str):
    mm = re.fullmatch(r"\s*(.+?)\s*", text)
    body, _, ifp = mm.group(1).partition(" if ")
    from ..commands._util import eval_vector
    v = np.broadcast_to(np.asarray(eval_vector(s, body), dtype=float), (s.data.nobs,))
    sel = m.touse.copy()
    if ifp.strip():
        from ..lang.vexpr import truth
        t = np.broadcast_to(np.asarray(truth(eval_vector(s, ifp)), dtype=bool), (s.data.nobs,))
        sel &= t
    return v, sel


def cmd_mlsum(s: "Session", args: str) -> None:
    body, o = _tool_split(args, {"noweight": 3})
    m = getattr(s, "_ml_active", None)
    if m is None:
        raise StataError(198, "mlsum can only be used within an ml evaluator")
    name, eq_sign, expr = body.partition("=")
    if not eq_sign:
        raise StataError(198, "invalid syntax")
    v, sel = _expr(s, m, expr)
    x = v[sel]
    if np.any(x >= SYS):
        s.scalars[name.strip()] = SYS
        return
    wfull = np.ones(s.data.nobs)
    if m.w is not None and not o.get("noweight"):
        wfull = np.zeros(s.data.nobs)
        wfull[m.touse] = m.w
    s.scalars[name.strip()] = float((x * wfull[sel]).sum())


def _x_full(m: Model, e: Eq, nobs: int) -> np.ndarray:
    X = np.zeros((nobs, e.k))
    X[m.touse] = e.X
    return X


def cmd_mlvecsum(s: "Session", args: str) -> None:
    body, o = _tool_split(args, {"eq": 2})
    m = getattr(s, "_ml_active", None)
    if m is None:
        raise StataError(198, "mlvecsum can only be used within an ml evaluator")
    left, eq_sign, expr = body.partition("=")
    names = left.split()
    if len(names) != 2 or not eq_sign:
        raise StataError(198, "invalid syntax")
    lnf, vec = names
    e = m.eqs[_eq_index(m, o.get("eq"))]
    v, sel = _expr(s, m, expr)
    x = v[sel]
    if np.any(x >= SYS):
        s.scalars[lnf] = SYS
        return
    wfull = np.ones(s.data.nobs)
    if m.w is not None:
        wfull = np.zeros(s.data.nobs)
        wfull[m.touse] = m.w
    X = _x_full(m, e, s.data.nobs)[sel]
    from ..commands.matrix import Matrix
    row = (x * wfull[sel]) @ X
    _mats(s)[vec] = Matrix(row.reshape(1, -1), ["r1"],
                           [c for q, c in m.names if q == e.name])


def cmd_mlmatsum(s: "Session", args: str) -> None:
    body, o = _tool_split(args, {"eq": 2})
    m = getattr(s, "_ml_active", None)
    if m is None:
        raise StataError(198, "mlmatsum can only be used within an ml evaluator")
    left, eq_sign, expr = body.partition("=")
    names = left.split()
    if len(names) != 2 or not eq_sign:
        raise StataError(198, "invalid syntax")
    lnf, mat = names
    ev = o.get("eq")
    if isinstance(ev, str) and "," in ev:
        i1, i2 = (int(x) for x in ev.split(","))
    else:
        i1 = int(ev) if isinstance(ev, str) else 1
        i2 = i1
    e1, e2 = m.eqs[i1 - 1], m.eqs[i2 - 1]
    v, sel = _expr(s, m, expr)
    x = v[sel]
    if np.any(x >= SYS):
        s.scalars[lnf] = SYS
        return
    wfull = np.ones(s.data.nobs)
    if m.w is not None:
        wfull = np.zeros(s.data.nobs)
        wfull[m.touse] = m.w
    X1 = _x_full(m, e1, s.data.nobs)[sel]
    X2 = _x_full(m, e2, s.data.nobs)[sel]
    from ..commands.matrix import Matrix
    _mats(s)[mat] = Matrix(X1.T @ (X2 * (x * wfull[sel])[:, None]),
                           [c for q, c in m.names if q == e1.name], [c for q, c in m.names if q == e2.name])


# ---------------------------------------------------------------------------
# mlexp
# ---------------------------------------------------------------------------

_BRACE = re.compile(r"\{([^{}]*)\}")


def cmd_mlexp(s: "Session", args: str) -> None:
    t = args.strip()
    if not t.startswith("("):
        raise StataError(198, "invalid syntax")
    j = _balanced(t, 0)
    expr, tail = t[1:j].strip(), t[j + 1:]
    p = parse_standard(tail)
    o = match_options(p.options, dict(_OPT, initial=4, variables=3)) if p.options.strip() else {}
    if o.get("constraints") or o.get("svy"):
        raise StataError(198, "option constraints() and svy are not available yet in OpenDTA")
    ds = s.data
    combos: dict[str, Eq] = {}
    scalars: dict[str, float] = {}
    order: list[str] = []

    def sub(mt):
        body = mt.group(1).strip()
        if ":" in body:
            name, _, rhs = body.partition(":")
            name = name.strip()
            lst, _, opt = rhs.partition(",")
            if name not in combos:
                words = lst.split()
                has_cons = "_cons" in words      # VERIFY: the manual lists _cons explicitly
                vars_ = s.expand_varlist(" ".join(w for w in words if w != "_cons")) if words else []
                combos[name] = Eq(name, [], vars_, has_cons or "noconstant" not in opt)
                order.append(name)
            return f"\x00c{name}\x00"
        name, _, init = body.partition("=")
        name = name.strip()
        if name not in scalars:
            scalars[name] = float(s.eval(init)) if init.strip() else 0.0
            order.append("/" + name)
        return f"\x00s{name}\x00"
    template = _BRACE.sub(sub, expr)
    if not order:
        raise StataError(198, "no parameters in the expression")
    from ..commands._util import touse
    from ..commands.summarize import weights
    mask = touse(s, p)
    w, wtype, mask = weights(s, p, mask, ("fweight", "aweight", "pweight", "iweight"))
    used = {v for e in combos.values() for v in e.indep}
    used |= {x for x in re.findall(r"[A-Za-z_]\w*", _BRACE.sub("", expr)) if ds.has(x)}
    clustvar = ""
    vce = o.get("vce")
    if isinstance(o.get("cluster"), str):
        clustvar = o["cluster"]
    elif isinstance(vce, str) and vce.split()[0].lower().startswith("cl"):
        clustvar = vce.split()[1]
    if clustvar:
        used.add(clustvar)
    for v in used:
        var = ds.get(v)
        mask = mask & (np.array([x != "" for x in var.raw]) if var.is_string else var.data < SYS)
    eqs: list[Eq] = []
    pos = 0
    for key in order:
        if key.startswith("/"):
            e = Eq(key[1:], [], [], True, scalar=True)
        else:
            e = combos[key]
        cols = [ds.get(v).data[mask] for v in e.indep]
        if e.cons:
            cols.append(np.ones(int(mask.sum())))
        e.X = np.column_stack(cols) if cols else np.zeros((int(mask.sum()), 0))
        e.off = np.zeros(int(mask.sum()))
        e.start = pos
        pos += e.k
        eqs.append(e)
    m = Model("mlexp", "", eqs, mask, None if w is None else w[mask], wtype, o, f"mlexp {args}")
    m.N = int(mask.sum())
    if m.N == 0:
        raise StataError(2000, "no observations")
    if clustvar:
        cv = ds.get(clustvar)
        m.cluster = np.asarray(cv.raw if cv.is_string else cv.data)[mask]
        m.clustvar = clustvar
    for e in eqs:
        if e.scalar:
            m.names.append(("/", e.name))
        else:
            m.names += [(e.name, v) for v in e.indep] + ([(e.name, "_cons")] if e.cons else [])
    m.template, m.lexp = template, expr
    m.b0 = np.zeros(m.nb)
    for e in eqs:
        if e.scalar:
            m.b0[e.start] = scalars[e.name]
    if isinstance(o.get("from"), str):
        s._ml = m
        _init(s, o["from"])
    elif isinstance(o.get("initial"), str):
        s._ml = m
        _init(s, " ".join(f"{k}={v}" if "=" not in k else k for k, v in
                          zip(o["initial"].split()[::2], o["initial"].split()[1::2])))
    elif not m.b0.any():
        m.b0 = None          # no starting values given: the search steps apply
    s._ml = m
    _maximize(s, "")


def _mlexp_obs(s: "Session", m: Model, th: np.ndarray) -> np.ndarray:
    """mlexp: the expression evaluated for the given linear predictors (one column per
    parameter group, like lf); returns the log likelihood per observation (NaN when missing)."""
    from ..commands._util import eval_vector
    from ..commands.syntaxcmd import _temp_name
    from ..core.dataset import Variable
    ds = s.data
    names = getattr(m, "_names", None)
    if names is None or any(not ds.has(n) for n in names.values()):
        names = {e.name: _temp_name(s) for e in m.eqs if not e.scalar}
        for nm in names.values():
            ds.add(Variable(nm, "double", np.full(ds.nobs, SYS)))
        m._names = names
        m._tmp = list(names.values())
    text = m.template
    for j, e in enumerate(m.eqs):
        if e.scalar:
            text = text.replace(f"\x00s{e.name}\x00", f"({float(th[0, j])!r})")
        else:
            ds.get(names[e.name]).data = _full(m, th[:, j], ds.nobs)
            text = text.replace(f"\x00c{e.name}\x00", names[e.name])
    v = np.broadcast_to(np.asarray(eval_vector(s, text), dtype=float), (ds.nobs,))[m.touse]
    return np.where(v >= SYS, np.nan, v)
