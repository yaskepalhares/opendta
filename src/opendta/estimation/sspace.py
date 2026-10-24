"""sspace and dfactor: linear Gaussian state-space models ([TS] sspace, [TS] dfactor), estimated by
maximum likelihood with the Kalman filter of ssm.py.

Both commands translate their syntax into a ModelSpec (parameter list plus a function that builds
the system matrices) and share the fitting and display code."""

from __future__ import annotations

import math
import re
from dataclasses import dataclass, field
from typing import Callable, TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.stats import chi2_tail
from ..lang.syntax import find_top, match_options, parse_standard
from .ml import maximize, options as ml_options, print_log
from .results import CoefRow, Estimates, coef_table, post
from .ssm import kalman_loglik
from .tsops import require_ts

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


@dataclass
class Param:
    name: str                 # name in e(b)
    eq: str                   # equation shown in the table ("" for the variances)
    label: str                # row label
    kind: str = "id"          # id or var
    start: float = 0.0
    hdr: str = ""             # variable shown on its own row above a lagged term (L1. under it)


@dataclass
class ModelSpec:
    cmd: str
    title: str
    deps: list[str]
    params: list[Param]
    build: Callable            # build(theta: dict[name, value], X: ndarray) -> D, A, S, H, c, d
    xcols: list[np.ndarray] = field(default_factory=list)
    diffuse: bool = False


def _split_eqs(args: str) -> tuple[list[str], str]:
    """Parenthesized equations at the top level, and the text after them."""
    eqs, i, n = [], 0, len(args)
    while i < n:
        if args[i] == "(":
            depth, j = 0, i
            while j < n:
                depth += args[j] == "("
                depth -= args[j] == ")"
                j += 1
                if depth == 0:
                    break
            eqs.append(args[i + 1:j - 1])
            i = j
        elif args[i].isspace():
            i += 1
        else:
            return eqs, args[i:]       # [if] [in] [, options]
    return eqs, ""


def parse_sspace(s: "Session", args: str):
    """Equations of sspace: list of dicts (dep, terms, state, noerror, noconstant)."""
    eqs, rest = _split_eqs(args)
    if not eqs:
        raise StataError(198, "equations required")
    out = []
    for e in eqs:
        c = find_top(e, ",")
        body, op = (e, "") if c == -1 else (e[:c], e[c + 1:])
        toks = body.split()
        if len(toks) < 1:
            raise StataError(198, "invalid equation")
        w = op.split()
        out.append({"dep": toks[0], "terms": toks[1:], "state": "state" in w,
                    "noerror": "noerror" in w, "noconstant": "noconstant" in w})
    return out, rest


def build_sspace(s: "Session", eqs: list[dict], covstate: str, covobs: str):
    if covstate not in ("diagonal", "identity") or covobs != "diagonal":
        raise StataError(198, "covstate() other than diagonal and identity, and covobserved() other than diagonal, "
                              "are not implemented in OpenDTA")   # OpenDTA
    states = [e["dep"] for e in eqs if e["state"]]
    obs = [e for e in eqs if not e["state"]]
    if not states or not obs:
        raise StataError(198, "at least one state and one observation equation required")
    r, m = len(states), len(obs)
    params: list[Param] = []
    pidx: dict[str, int] = {}        # shared parameter names
    fixed: dict[tuple, float] = {}
    entries = []                     # (kind, row, col or xindex, ref)  ref: ("p", name) | ("f", value)
    xnames: list[str] = []
    errs_state = [None] * r
    errs_obs = [None] * m

    def coef(eq_label, label, ename, spec, start, hdr=""):
        if spec is None:
            nm = f"{eq_label}:{ename}"
            params.append(Param(nm, eq_label, label, "id", start, hdr))
            return ("p", nm)
        try:
            return ("f", float(spec))
        except ValueError:
            if spec not in pidx:
                pidx[spec] = 1
                params.append(Param(spec, eq_label, label, "id", start))
            return ("p", spec)

    for k, e in enumerate(eqs if False else states):
        pass
    for kind, lst in (("S", [e for e in eqs if e["state"]]), ("O", obs)):
        for i, e in enumerate(lst):
            for t in e["terms"]:
                spec = None
                if "@" in t:
                    t, spec = t.split("@", 1)
                if t.lower().startswith("e."):
                    if spec is not None:
                        raise StataError(198, "error terms may not carry coefficients")
                    nm = t[2:]
                    if nm != e["dep"]:
                        raise StataError(198, f"error term {t} must belong to equation {e['dep']}")   # VERIFY
                    (errs_state if kind == "S" else errs_obs)[i] = nm
                    continue
                lagm = re.fullmatch(r"[lL](\d*)\.(\w+)", t)
                base = lagm.group(2) if lagm else t
                lag = (int(lagm.group(1) or 1) if lagm else 0)
                if kind == "S":
                    if lag != 1 or base not in states:
                        raise StataError(198, "state equations may contain only lagged states and error terms in OpenDTA")   # OpenDTA
                    entries.append(("A", i, states.index(base), coef(e["dep"], "L1.", t, spec, 0.1, base)))
                else:
                    if lag == 0 and base in states:
                        entries.append(("D", i, states.index(base), coef(e["dep"], base, t, spec, 1.0)))
                    else:
                        if base in states or lag:
                            raise StataError(198, "observation equations may contain only current states and exogenous variables in OpenDTA")   # OpenDTA
                        if t not in xnames:
                            xnames.append(t)
                        entries.append(("X", i, xnames.index(t), coef(e["dep"], t, t, spec, 0.0)))
            if not e["noconstant"]:
                entries.append(("C" if kind == "S" else "K", i, 0, coef(e["dep"], "_cons", "_cons", None, 0.0)))
    implicit = not any(nm for nm in errs_state + errs_obs)
    if implicit:
        # covariance form: every equation has an error unless it says noerror
        errs_state = [None if e["noerror"] else e["dep"] for e in eqs if e["state"]]
        errs_obs = [None if e["noerror"] else e["dep"] for e in obs]
    unit = set()
    for i, nm in enumerate(errs_state):
        if nm and not [e for e in eqs if e["state"]][i]["noerror"]:
            if covstate == "identity":
                unit.add(i)
                continue
            params.append(Param(f"var(e.{states[i]})", "", f"var({states[i]})" if implicit else f"var(e.{states[i]})", "var", 1.0))
    for i, nm in enumerate(errs_obs):
        if nm and not obs[i]["noerror"]:
            params.append(Param(f"var(e.{obs[i]['dep']})", "", f"var({obs[i]['dep']})" if implicit else f"var(e.{obs[i]['dep']})", "var", 1.0))
    # an equation without the e. term is exact only when noerror is given or the term is omitted
    deps = [e["dep"] for e in obs]
    pnames = {p.name for p in params}

    def val(ref, th):
        return th[ref[1]] if ref[0] == "p" else ref[1]

    def build(th, X):
        T = X.shape[0] if X.size else th["__T"]
        A = np.zeros((r, r)); D = np.zeros((m, r)); c = np.zeros(r)
        k0 = np.zeros(m); Bx = np.zeros((m, max(len(xnames), 1)))
        for kd, i, j, ref in entries:
            v = val(ref, th)
            if kd == "A": A[i, j] = v
            elif kd == "D": D[i, j] = v
            elif kd == "X": Bx[i, j] = v
            elif kd == "C": c[i] = v
            else: k0[i] = v
        S = np.zeros((r, r)); H = np.zeros((m, m))
        for i in range(r):
            if i in unit:
                S[i, i] = 1.0
            elif f"var(e.{states[i]})" in th:
                S[i, i] = th[f"var(e.{states[i]})"]
        for i in range(m):
            if f"var(e.{deps[i]})" in th:
                H[i, i] = th[f"var(e.{deps[i]})"]
        d = k0[None, :] + (X @ Bx[:, :X.shape[1]].T if xnames else 0.0)
        return D, A, S, H, c, np.broadcast_to(d, (T, m)).copy()

    return ModelSpec("sspace", "State-space model", deps, params, build), xnames


def build_dfactor(s: "Session", obs_eq: tuple, fac_eq: tuple):
    """obs_eq = (deps, exog, noconstant); fac_eq = (name, arlags)."""
    deps, exog, nocons = obs_eq
    fname, arl = fac_eq
    p = max(arl) if arl else 0
    r = max(p, 1)
    m = len(deps)
    params: list[Param] = []
    for l in arl:
        params.append(Param(f"{fname}:L{l}.{fname}", fname, f"L{l}.", "id", 0.3 if l == 1 else 0.0, fname))
    for i, dv in enumerate(deps):
        params.append(Param(f"{dv}:{fname}", dv, fname, "id", 0.5 + 0.1 * i))
        for x in exog:
            params.append(Param(f"{dv}:{x}", dv, x, "id", 0.0))
        if not nocons:
            params.append(Param(f"{dv}:_cons", dv, "_cons", "id", 0.0))
    for dv in deps:
        params.append(Param(f"var(e.{dv})", "", f"var(e.{dv})", "var", 1.0))

    def build(th, X):
        T = X.shape[0]
        A = np.zeros((r, r))
        for l in arl:
            A[0, l - 1] = th[f"{fname}:L{l}.{fname}"]
        for i in range(1, r):
            A[i, i - 1] = 1.0
        S = np.zeros((r, r)); S[0, 0] = 1.0           # the factor variance is normalized to one
        D = np.zeros((m, r)); H = np.zeros((m, m))
        d = np.zeros((T, m))
        for i, dv in enumerate(deps):
            D[i, 0] = th[f"{dv}:{fname}"]
            H[i, i] = th[f"var(e.{dv})"]
            for k, x in enumerate(exog):
                d[:, i] += th[f"{dv}:{x}"] * X[:, k]
            if not nocons:
                d[:, i] += th[f"{dv}:_cons"]
        return D, A, S, H, np.zeros(r), d

    return ModelSpec("dfactor", "Dynamic-factor model", deps, params, build), fname


def _to_nat(u, kind):
    return math.exp(min(u, 700)) if kind == "var" else u


def _fit(s: "Session", spec: ModelSpec, xnames: list[str], mask: np.ndarray, o: dict, cmdline: str,
         arg_extra: dict | None = None) -> None:
    from .tsmodels import _check_gaps, _numderiv, _sample_text
    ds = s.data
    info = require_ts(ds)
    Y = np.column_stack([ds.get(v).data.astype(float) for v in spec.deps])
    Y[Y >= SYS] = np.nan
    X = np.column_stack([ds.get(x).data.astype(float) for x in xnames]) if xnames else np.zeros((ds.nobs, 0))
    ok = mask & np.any(np.isfinite(Y), axis=1)
    if xnames:
        ok &= np.all(X < SYS, axis=1)
    t0, t1 = _check_gaps(s, ok)
    order = np.argsort(ds.get(info.tvar).data[ok], kind="stable")
    y, Xo = Y[ok][order], X[ok][order]
    T = y.shape[0]
    names = [p.name for p in spec.params]
    kinds = [p.kind for p in spec.params]
    if T < 2:
        raise StataError(2001, "insufficient observations")

    def theta(u):
        d = {p.name: _to_nat(v, p.kind) for p, v in zip(spec.params, u)}
        d["__T"] = T
        return d

    def ll(u):
        if any(abs(v) > 25 for v, k in zip(u, kinds) if k == "var"):
            return -np.inf            # variances beyond exp(+-25) are numerically meaningless
        try:
            D, A, S, H, c, d = spec.build(theta(u), Xo)
            sm = spec.diffuse
            dm = np.ones(A.shape[0], dtype=bool) if sm else None
            r = kalman_loglik(y, D, A, S, H, c, d, dm)
        except (np.linalg.LinAlgError, OverflowError):
            return -np.inf
        return -np.inf if r is None else float(r.sum())

    def fun(u, todo):
        if todo == 0:
            return ll(u), None, None
        return _numderiv(ll, u)

    vy = float(np.nanvar(y)) or 1.0
    u0 = np.array([math.log(vy * 0.5) if p.kind == "var" else p.start for p in spec.params])
    if spec.cmd == "sspace":
        # constants of the observation equations start at the sample means
        for i, p in enumerate(spec.params):
            if p.label == "_cons" and p.eq in spec.deps:
                u0[i] = float(np.nanmean(y[:, spec.deps.index(p.eq)]))
    if spec.cmd == "dfactor":
        for i, p in enumerate(spec.params):
            if p.label == "_cons" and p.eq in spec.deps:
                u0[i] = float(np.nanmean(y[:, spec.deps.index(p.eq)]))
            if p.name.startswith("var(e."):
                u0[i] = math.log(vy * 0.5)
    if not np.isfinite(ll(u0)):
        # a smaller autoregressive start usually restores feasibility
        for i, p in enumerate(spec.params):
            if p.kind == "id" and ("L" in p.label and p.label != "_cons"):
                u0[i] = 0.0
    Tm = a_c = Cmat = None
    if o.get("constraints") not in (None, True):
        # constraints() on the coefficients (not on the variances, which are estimated on the log
        # scale): u = Tm uc + a (core/constraints.py)
        from types import SimpleNamespace
        from ..core.constraints import constraint_transform
        from .postest import _constraint
        from .sysmodels import _constraints, _numlist
        fake = SimpleNamespace(names=[p.label for p in spec.params], eqnames=[p.eq or "" for p in spec.params])
        defined = _constraints(s)
        crow = []
        for cn_ in _numlist(str(o["constraints"])):
            if cn_ not in defined:
                raise StataError(111, f"constraint {cn_} not found")
            coef, rhs = _constraint(s, fake, defined[cn_])
            if any(spec.params[j].kind == "var" for j in coef):
                raise StataError(198, "constraints on variances are not available in OpenDTA")   # OpenDTA
            row = np.zeros(len(u0) + 1)
            for j, v in coef.items():
                row[j] += v
            row[-1] = rhs
            crow.append(row)
        Cmat = np.array(crow)
        Tm, a_c = constraint_transform(Cmat)
        a_c = a_c.ravel()
        u0 = a_c + Tm @ (Tm.T @ (u0 - a_c))
    if not np.isfinite(ll(u0)):
        raise StataError(430, "initial values not feasible")
    if Tm is not None:
        def cfun(uc, todo):
            l_, g_, H_ = fun(Tm @ uc + a_c, todo)
            return (l_, g_, H_) if g_ is None else (l_, Tm.T @ g_, Tm.T @ H_ @ Tm)
        r = maximize(cfun, Tm.T @ (u0 - a_c), **ml_options(o))
        r.b = Tm @ r.b + a_c
    else:
        r = maximize(fun, u0, **ml_options(o))
    out = s.output
    if not o.get("nolog"):
        out.write("\n", "text")
        print_log(out, r, "log likelihood")
    u = r.b
    _l, _g, Hm = _numderiv(ll, u)
    Tq = Tm if Tm is not None else np.eye(len(u))
    try:
        Vu = Tq @ np.linalg.inv(-(Tq.T @ Hm @ Tq)) @ Tq.T
    except np.linalg.LinAlgError:
        Vu = Tq @ np.linalg.pinv(-(Tq.T @ Hm @ Tq)) @ Tq.T
    nat = np.array([_to_nat(v, k) for v, k in zip(u, kinds)])
    J = np.array([nat[i] if kinds[i] == "var" else 1.0 for i in range(len(u))])
    V = J[:, None] * Vu * J[None, :]
    V = (V + V.T) / 2
    se = np.sqrt(np.maximum(np.diag(V), 0))
    D, A, S, H, c, d = spec.build(theta(u), Xo)
    # Wald test of all the coefficients other than the variances
    idx = [i for i, p in enumerate(spec.params) if p.kind == "id"]
    df_w = len(idx)
    try:
        if Tm is not None and idx:
            # the test has the rank of the constrained variance (VERIFY)
            df_w = int(np.linalg.matrix_rank(V[np.ix_(idx, idx)]))
            chi2 = float(nat[idx] @ np.linalg.pinv(V[np.ix_(idx, idx)]) @ nat[idx]) if df_w else SYS
        else:
            chi2 = float(nat[idx] @ np.linalg.solve(V[np.ix_(idx, idx)], nat[idx])) if idx else SYS
    except np.linalg.LinAlgError:
        chi2 = SYS
    pv = float(chi2_tail(chi2, df_w)) if chi2 < SYS else SYS
    rows: list[CoefRow] = []
    last, lasth = None, ""
    for i, p in enumerate(spec.params):
        if p.eq and p.eq != last:
            if last is not None:
                rows.append(CoefRow("sep"))
            rows.append(CoefRow("eq", p.eq))
            last, lasth = p.eq, ""
        elif not p.eq and last is not None:
            rows.append(CoefRow("sep"))
            last, lasth = None, ""
        if p.hdr and p.hdr != lasth:
            if lasth:
                rows.append(CoefRow("blank"))
            rows.append(CoefRow("header", p.hdr))
        lasth = p.hdr
        kind_ = "cons" if Tm is not None and se[i] < 1e-12 else "coef"     # fixed by a constraint
        rows.append(CoefRow(kind_, p.label, nat[i], se[i], name=p.name, tight=p.kind == "var", zero_lo=p.kind == "var"))
    est = Estimates(spec.cmd, " ".join(spec.deps), names, nat, V, T, "z", None)
    est.eqnames = [p.eq or "_" for p in spec.params]
    est.rows = rows
    est.level = float(o["level"]) if isinstance(o.get("level"), str) else float(s.settings.get("level", "95"))
    est.vcetype = "OIM"
    cns_text = []
    if Cmat is not None:
        from .postest import _ctext
        cns_text = [_ctext(fake, {j: float(v) for j, v in enumerate(row[:-1]) if v}, float(row[-1])) for row in Cmat]
    est.extra = {"N": float(T), "ll": r.ll, "chi2": chi2, "p": pv, "df_m": df_w, "title": spec.title,
                 "cns_text": cns_text,
                 "sample": _sample_text(s, t0, t1)}
    est.ssm = {"spec": spec, "u": u, "y": y, "X": Xo, "mask": ok.copy(), "order": order}
    est.display = display_ssm
    est.predict = _predict_ssm
    nvar = sum(1 for k in kinds if k == "var")
    scalars = [("N", float(T)), ("k", float(len(u))), ("k_aux", float(nvar)), ("k_eq", float(len(spec.deps))),
               ("k_dv", float(len(spec.deps))), ("df_m", float(len(idx))), ("ll", r.ll), ("chi2", chi2), ("p", pv),
               ("tmin", t0), ("tmax", t1), ("rank", float(len(u))), ("ic", float(r.iterations)), ("rc", 0.0),
               ("converged", 1.0 if r.converged else 0.0)]
    macros = [("cmd", spec.cmd), ("cmdline", cmdline), ("depvar", " ".join(spec.deps)), ("title", spec.title),
              ("tvar", info.tvar), ("vce", "oim"), ("vcetype", "OIM"), ("chi2type", "Wald"),
              ("properties", "b V")]
    post(s, est, scalars, macros, ok)
    display_ssm(s, est)


def display_ssm(s: "Session", est: Estimates, *, header: bool = True, table: bool = True,
                level: float | None = None, **_kw) -> None:
    from .tsmodels import _g, _right
    out = s.output
    x = est.extra
    if header:
        out.write(f"\n{x['title']}\n\n", "text")
        _right(out, x["sample"].replace("Sample:  ", "Sample: "), "Number of obs", f"{int(x['N']):,}")
        _right(out, "", f"Wald chi2({int(x['df_m'])})", _g(x["chi2"], "%10.2f"))
        _right(out, f"Log likelihood = {_g(x['ll'], '%10.0g'):>10}", "Prob > chi2", _g(x["p"], "%10.4f"))
        if x.get("cns_text"):
            out.write("\n", "text")
            for i, t_ in enumerate(x["cns_text"], 1):
                out.write(f" ({i:>2})  {t_}\n", "text")              # VERIFY: place of the list
    if table:
        coef_table(s, est.depvar.split()[0] if est.cmd == "sspace" else "", est.rows, stat="z", df=None, level=level or est.level, vcetype=est.vcetype)
        out.write("Note: Tests of variances against zero are one sided, and the two-sided\n"
                  "      confidence intervals are truncated at zero.\n", "text")


_OPTS = {"covstate": 3, "covobserved": 3, "constraints": 3, "from": 3, "vce": 3, "robust": 1, "level": 1,
         "nolog": 3, "iterate": 4, "tolerance": 3, "ltolerance": 4, "nrtolerance": 4, "technique": 4,
         "difficult": 4, "trace": 5, "diffuse": 4, "nocnsreport": 5, "collinear": 4, "covariance": 3,
         "noconstant": 4}


def _prep(s, args, name):
    from .postest import replay
    from ..commands._util import touse
    if not args.strip() or args.strip().startswith(","):
        replay(s, name, args)
        return None
    return True


def cmd_sspace(s: "Session", args: str) -> None:
    from .regress import cmdline
    if not _prep(s, args, "sspace"):
        return
    # "if"/"in" may follow the equations before the comma
    eqs, rest = parse_sspace(s, args)
    c = find_top(rest, ",") if rest else -1
    pre, opts = (rest, "") if c == -1 else (rest[:c], rest[c + 1:])
    p = parse_standard(pre) if pre.strip() else parse_standard("")
    o = match_options(opts, _OPTS) if opts.strip() else {}
    for bad in ("from", "vce"):
        if o.get(bad):
            raise StataError(198, f"option {bad}() not implemented in OpenDTA")   # OpenDTA
    from ..commands._util import touse
    mask = touse(s, p)
    spec, xnames = build_sspace(s, eqs, str(o.get("covstate", "diagonal")), str(o.get("covobserved", "diagonal")))
    spec.diffuse = bool(o.get("diffuse"))
    _fit(s, spec, xnames, mask, o, cmdline("sspace", args))


def cmd_dfactor(s: "Session", args: str) -> None:
    from .regress import cmdline
    from ..commands._util import touse
    if not _prep(s, args, "dfactor"):
        return
    eqs, rest = _split_eqs(args)
    if len(eqs) != 2:
        raise StataError(198, "one observed-variables equation and one factor equation required")   # OpenDTA
    c = find_top(rest, ",") if rest else -1
    pre, opts = (rest, "") if c == -1 else (rest[:c], rest[c + 1:])
    p = parse_standard(pre) if pre.strip() else parse_standard("")
    o = match_options(opts, _OPTS) if opts.strip() else {}
    for bad in ("from", "vce"):
        if o.get(bad):
            raise StataError(198, f"option {bad}() not implemented in OpenDTA")   # OpenDTA
    if o.get("covariance") and str(o["covariance"]) != "diagonal":
        raise StataError(198, "covariance() other than diagonal not implemented in OpenDTA")   # OpenDTA
    parts = []
    for e in eqs:
        cc = find_top(e, ",")
        body, op = (e, "") if cc == -1 else (e[:cc], e[cc + 1:])
        if "=" not in body:
            raise StataError(198, "equations must have the form (varlist = ...)")
        lhs, rhs = body.split("=", 1)
        parts.append((lhs.split(), rhs.split(), op))
    (od, ox, oop), (fd, fx, fop) = parts
    if len(fd) != 1 or fx:
        raise StataError(198, "factor equation must be (factorname = , ar(numlist))")   # OpenDTA
    odeps = s.expand_varlist(" ".join(od))
    exog = s.expand_varlist(" ".join(ox)) if ox else []
    m_ = re.search(r"\bar\(([^)]*)\)", fop)
    from .tsmodels import _lags
    arl = _lags(m_.group(1)) if m_ else []
    if re.search(r"\bma\(|\barstructure|\bcovstructure|\bcovariance", oop + fop):
        raise StataError(198, "ma(), arstructure() and error autoregressions not implemented in OpenDTA")   # OpenDTA
    spec, _f = build_dfactor(s, (odeps, exog, "noconstant" in oop.split()), (fd[0], arl))
    _fit(s, spec, exog, touse(s, p), o, cmdline("dfactor", args))


def _predict_ssm(s: "Session", est: Estimates, stat: str, mask: np.ndarray):
    """One-step predictions (xb, the default), residuals and standardized residuals of one observed
    variable (equation(), the first by default) over the estimation sample. VERIFY: states, the
    smethod() and dynamic() options."""
    from .ssm import kalman_onestep
    M_ = est.ssm
    spec, u = M_["spec"], M_["u"]
    th = {p.name: _to_nat(v, p.kind) for p, v in zip(spec.params, u)}
    th["__T"] = M_["y"].shape[0]
    D, A, S, H, c, d = spec.build(th, M_["X"])
    dm = np.ones(A.shape[0], dtype=bool) if spec.diffuse else None
    yhat, fv = kalman_onestep(M_["y"], D, A, S, H, c, d, dm)
    eq = str((getattr(est, "_pred_opts", {}) or {}).get("equation") or "").strip().lstrip("#")
    j = 0
    if eq:
        j = int(eq) - 1 if eq.isdigit() else (spec.deps.index(eq) if eq in spec.deps else -1)
        if not 0 <= j < len(spec.deps):
            raise StataError(198, f"equation {eq} not found")                     # VERIFY
    rows = np.nonzero(M_["mask"])[0][M_["order"]]
    out = np.full(s.data.nobs, SYS)
    if stat in ("", "xb"):
        out[rows] = yhat[:, j]
        return out, "" if stat else "(option xb assumed; fitted values)"
    e = M_["y"][:, j] - yhat[:, j]
    if stat == "residuals":
        out[rows] = np.where(np.isfinite(e), e, SYS)
    elif stat == "rstandard":
        out[rows] = np.where(np.isfinite(e), e / np.sqrt(fv[:, j]), SYS)
    else:
        raise StataError(198, f"option {stat} not allowed")                      # VERIFY
    return out, ""
