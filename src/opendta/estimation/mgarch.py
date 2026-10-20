"""mgarch: multivariate GARCH models ([TS] mgarch): diagonal vech, constant conditional correlation
and dynamic conditional correlation, Gaussian errors, maximum likelihood with an OIM covariance."""

from __future__ import annotations

import math
import re
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.stats import chi2_tail
from ..lang.syntax import find_top, match_options, parse_standard
from .ml import maximize_switching, options as ml_options, print_log
from .results import CoefRow, Estimates, coef_table, post
from .sspace import _split_eqs
from .tsops import require_ts, shift

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS
TITLES = {"ccc": "Constant conditional correlation MGARCH model",
          "dcc": "Dynamic conditional correlation MGARCH model",
          "dvech": "Diagonal vech MGARCH model"}


@dataclass
class MGSpec:
    kind: str
    deps: list[str]
    xnames: list[list[str]]          # regressors of each equation (without the constant)
    cons: list[bool]
    arch: list[int]
    garch: list[int]


def _term(ds, t: str) -> np.ndarray:
    m = re.fullmatch(r"[lL](\d*)\.(\w+)", t)
    if m:
        return shift(ds, ds.get(m.group(2)).data.astype(float), int(m.group(1) or 1))
    return ds.get(t).data.astype(float)


def layout(spec: MGSpec):
    """Parameter list: (name, eq, label, group) in the order of e(b)."""
    m = len(spec.deps)
    P = []
    for i, dv in enumerate(spec.deps):
        for x in spec.xnames[i]:
            P.append((f"{dv}:{x}", dv, x, "mean"))
        if spec.cons[i]:
            P.append((f"{dv}:_cons", dv, "_cons", "mean"))
    if spec.kind == "dvech":
        pairs = [(i, j) for i in range(m) for j in range(i + 1)]
        for grp, lags, lab in (("omega", [0], "Sigma0"), ("arch", spec.arch, "L.ARCH"), ("garch", spec.garch, "L.GARCH")):
            for l in lags:
                eq = lab if grp == "omega" or l == 1 else lab.replace("L.", f"L{l}.")
                for i, j in pairs:
                    P.append((f"{eq}:{i + 1}_{j + 1}", eq, f"{i + 1}_{j + 1}", f"{grp}{l}"))
        return P
    for i, dv in enumerate(spec.deps):
        eq = f"ARCH_{dv}"
        for l in spec.arch:
            P.append((f"{eq}:L{l}.arch", eq, ("arch", l), f"arch{i}"))
        for l in spec.garch:
            P.append((f"{eq}:L{l}.garch", eq, ("garch", l), f"garch{i}"))
        P.append((f"{eq}:_cons", eq, "_cons", f"omega{i}"))
    for j in range(m):
        for i in range(j + 1, m):
            P.append((f"corr({spec.deps[i]},{spec.deps[j]})", "", f"corr({spec.deps[i]},{spec.deps[j]})", "corr"))
    if spec.kind == "dcc":
        P.append(("Adjustment:lambda1", "Adjustment", "lambda1", "lam"))
        P.append(("Adjustment:lambda2", "Adjustment", "lambda2", "lam"))
    return P


def _unpack(spec: MGSpec, b: np.ndarray, P):
    """Structured parameters from the flat vector."""
    m = len(spec.deps)
    g = {p[0]: v for p, v in zip(P, b)}
    beta = []
    for i, dv in enumerate(spec.deps):
        cols = list(spec.xnames[i]) + (["_cons"] if spec.cons[i] else [])
        beta.append(np.array([g[f"{dv}:{c}"] for c in cols]))
    out = {"beta": beta}
    if spec.kind == "dvech":
        pairs = [(i, j) for i in range(m) for j in range(i + 1)]

        def mat(eq):
            Mx = np.zeros((m, m))
            for i, j in pairs:
                Mx[i, j] = Mx[j, i] = g[f"{eq}:{i + 1}_{j + 1}"]
            return Mx
        out["omega"] = mat("Sigma0")
        out["arch"] = {l: mat("L.ARCH" if l == 1 else f"L{l}.ARCH") for l in spec.arch}
        out["garch"] = {l: mat("L.GARCH" if l == 1 else f"L{l}.GARCH") for l in spec.garch}
        return out
    out["a"] = [{l: g[f"ARCH_{dv}:L{l}.arch"] for l in spec.arch} for dv in spec.deps]
    out["g"] = [{l: g[f"ARCH_{dv}:L{l}.garch"] for l in spec.garch} for dv in spec.deps]
    out["w"] = [g[f"ARCH_{dv}:_cons"] for dv in spec.deps]
    R = np.eye(m)
    for j in range(m):
        for i in range(j + 1, m):
            R[i, j] = R[j, i] = g[f"corr({spec.deps[i]},{spec.deps[j]})"]
    out["R"] = R
    if spec.kind == "dcc":
        out["l1"], out["l2"] = g["Adjustment:lambda1"], g["Adjustment:lambda2"]
    return out


def residuals(spec: MGSpec, par, y, X):
    E = np.zeros_like(y)
    for i in range(len(spec.deps)):
        E[:, i] = y[:, i] - X[i] @ par["beta"][i]
    return E


def mgarch_loglik(spec: MGSpec, par, E: np.ndarray, return_h: bool = False):
    """Log likelihood contributions (length T), or None when the parameters are infeasible."""
    T, m = E.shape
    log2pi = math.log(2 * math.pi)
    out = np.zeros(T)
    maxl = max(list(spec.arch) + list(spec.garch) + [1])
    Hs = np.zeros((T, m, m)) if return_h else None
    if spec.kind == "dvech":
        S0 = E.T @ E / T                       # VERIFY: presample matrix is the average of e e'
        Hh = np.zeros((T, m, m))
        if np.min(np.linalg.eigvalsh(par["omega"])) <= 0:
            return None
        for t in range(T):
            H = par["omega"].copy()
            for l, Al in par["arch"].items():
                e = E[t - l] if t - l >= 0 else None
                H = H + Al * (np.outer(e, e) if e is not None else S0)
            for l, Bl in par["garch"].items():
                H = H + Bl * (Hh[t - l] if t - l >= 0 else S0)
            Hh[t] = H
            sign, ld = np.linalg.slogdet(H)
            if sign <= 0:
                return None
            try:
                q = float(E[t] @ np.linalg.solve(H, E[t]))
            except np.linalg.LinAlgError:
                return None
            out[t] = -0.5 * (m * log2pi + ld + q)
        return (out, Hh) if return_h else out
    # univariate variances
    h = np.zeros((T, m))
    for i in range(m):
        s0 = float(np.mean(E[:, i] ** 2))
        if par["w"][i] <= 0 or any(v < 0 for v in par["a"][i].values()) or any(v < 0 for v in par["g"][i].values()):
            return None
        for t in range(T):
            v = par["w"][i]
            for l, al in par["a"][i].items():
                v += al * (E[t - l, i] ** 2 if t - l >= 0 else s0)
            for l, gl in par["g"][i].items():
                v += gl * (h[t - l, i] if t - l >= 0 else s0)
            h[t, i] = v
    if np.any(h <= 0) or not np.all(np.isfinite(h)):
        return None
    U = E / np.sqrt(h)
    R = par["R"]
    if spec.kind == "ccc":
        if np.min(np.linalg.eigvalsh(R)) <= 1e-8:
            return None
        sign, ldR = np.linalg.slogdet(R)
        Ri = np.linalg.inv(R)
        q = np.einsum("ti,ij,tj->t", U, Ri, U)
        out = -0.5 * (m * log2pi + np.log(h).sum(axis=1) + ldR + q)
        if return_h:
            for t in range(T):
                d = np.sqrt(h[t])
                Hs[t] = R * np.outer(d, d)
            return out, Hs
        return out
    l1, l2 = par["l1"], par["l2"]
    if l1 < 0 or l2 < 0 or l1 + l2 >= 1 or np.min(np.linalg.eigvalsh(R)) <= 1e-8:
        return None
    Q = R.copy()
    uu_prev = R.copy()                         # VERIFY: presample u u' is the correlation matrix
    for t in range(T):
        Q = (1 - l1 - l2) * R + l1 * uu_prev + l2 * Q
        dq = 1 / np.sqrt(np.diag(Q))
        Rt = Q * np.outer(dq, dq)
        sign, ldR = np.linalg.slogdet(Rt)
        if sign <= 0:
            return None
        q = float(U[t] @ np.linalg.solve(Rt, U[t]))
        out[t] = -0.5 * (m * log2pi + np.log(h[t]).sum() + ldR + q)
        uu_prev = np.outer(U[t], U[t])
        if return_h:
            d = np.sqrt(h[t])
            Hs[t] = Rt * np.outer(d, d)
    return (out, Hs) if return_h else out


def start_values(spec: MGSpec, P, y, X):
    m = len(spec.deps)
    b = np.zeros(len(P))
    idx = {p[0]: k for k, p in enumerate(P)}
    E = np.zeros_like(y)
    for i, dv in enumerate(spec.deps):
        cols = list(spec.xnames[i]) + (["_cons"] if spec.cons[i] else [])
        bi = np.linalg.lstsq(X[i], y[:, i], rcond=None)[0] if X[i].shape[1] else np.zeros(0)
        for c, v in zip(cols, bi):
            b[idx[f"{dv}:{c}"]] = v
        E[:, i] = y[:, i] - X[i] @ bi
    S = E.T @ E / len(E)
    na, ng = len(spec.arch), len(spec.garch)
    if spec.kind == "dvech":
        pairs = [(i, j) for i in range(m) for j in range(i + 1)]
        for i, j in pairs:
            b[idx[f"Sigma0:{i + 1}_{j + 1}"]] = S[i, j] * 0.1
            for l in spec.arch:
                eq = "L.ARCH" if l == 1 else f"L{l}.ARCH"
                b[idx[f"{eq}:{i + 1}_{j + 1}"]] = 0.1 / na if i == j else 0.02 / na
            for l in spec.garch:
                eq = "L.GARCH" if l == 1 else f"L{l}.GARCH"
                b[idx[f"{eq}:{i + 1}_{j + 1}"]] = 0.7 / ng if i == j else 0.3 / ng
        return b
    for i, dv in enumerate(spec.deps):
        for l in spec.arch:
            b[idx[f"ARCH_{dv}:L{l}.arch"]] = 0.1 / na
        for l in spec.garch:
            b[idx[f"ARCH_{dv}:L{l}.garch"]] = 0.7 / ng
        b[idx[f"ARCH_{dv}:_cons"]] = S[i, i] * 0.2
    C = np.corrcoef(E.T) if m > 1 else np.eye(1)
    for j in range(m):
        for i in range(j + 1, m):
            b[idx[f"corr({spec.deps[i]},{spec.deps[j]})"]] = float(np.clip(C[i, j], -0.9, 0.9))
    if spec.kind == "dcc":
        b[idx["Adjustment:lambda1"]] = 0.05
        b[idx["Adjustment:lambda2"]] = 0.8
    return b


def _scores(contrib, b: np.ndarray) -> np.ndarray:
    """Observation scores by central differences; one-sided where the parameters are infeasible on
    one side (a variance or adjustment parameter on its boundary)."""
    c0 = contrib(b)
    S = np.zeros((len(c0), b.size))
    for j in range(b.size):
        h = 1e-6 * (abs(b[j]) + 1)
        e = np.zeros(b.size)
        e[j] = h
        cp, cm = contrib(b + e), contrib(b - e)
        if cp is not None and cm is not None:
            S[:, j] = (cp - cm) / (2 * h)
        elif cp is not None:
            S[:, j] = (cp - c0) / h
        elif cm is not None:
            S[:, j] = (c0 - cm) / h
    return S


_OPTS = {"arch": 3, "garch": 3, "distribution": 3, "het": 3, "constraints": 3, "from": 3, "vce": 3, "robust": 1,
         "level": 1, "nolog": 3, "iterate": 4, "tolerance": 3, "ltolerance": 4, "nrtolerance": 4,
         "technique": 4, "difficult": 4, "trace": 5, "nocnsreport": 5, "collinear": 4, "noconstant": 4}


def cmd_mgarch(s: "Session", args: str) -> None:
    from ..commands._util import touse
    from .postest import replay
    from .regress import cmdline
    from .tsmodels import _check_gaps, _lags, _numderiv, _sample_text
    if not args.strip() or args.strip().startswith(","):
        replay(s, "mgarch", args)
        return
    w = args.split(None, 1)
    kind = w[0].lower()
    if kind in ("vcc",):
        raise StataError(198, "mgarch vcc not implemented in OpenDTA")   # OpenDTA
    if kind not in TITLES:
        raise StataError(198, "invalid syntax; mgarch requires dvech, ccc, dcc or vcc")
    eqs, rest = _split_eqs(w[1] if len(w) > 1 else "")
    if not eqs:
        raise StataError(198, "equation required")
    c = find_top(rest, ",") if rest else -1
    pre, opts = (rest, "") if c == -1 else (rest[:c], rest[c + 1:])
    p = parse_standard(pre) if pre.strip() else parse_standard("")
    o = match_options(opts, _OPTS) if opts.strip() else {}
    for bad in ("constraints", "from", "vce", "het"):
        if o.get(bad):
            raise StataError(198, f"option {bad}() not implemented in OpenDTA")   # OpenDTA
    if isinstance(o.get("distribution"), str) and o["distribution"].split()[0].lower() not in ("gaussian", "normal"):
        raise StataError(198, "distribution() other than gaussian not implemented in OpenDTA")   # OpenDTA
    deps, xnames, cons = [], [], []
    for e in eqs:
        cc = find_top(e, ",")
        body, op = (e, "") if cc == -1 else (e[:cc], e[cc + 1:])
        if "=" in body:
            lhs, rhs = body.split("=", 1)
        else:
            lhs, rhs = body, ""
        for dv in s.expand_varlist(lhs):
            deps.append(dv)
            xnames.append(rhs.split())
            cons.append("noconstant" not in op.split())
    if len(deps) < 2:
        raise StataError(198, "at least two dependent variables required")
    arch = _lags(o["arch"]) if isinstance(o.get("arch"), str) else [1]
    garch = _lags(o["garch"]) if isinstance(o.get("garch"), str) else [1]
    spec = MGSpec(kind, deps, xnames, cons, arch, garch)
    ds = s.data
    info = require_ts(ds)
    mask = touse(s, p)
    n = ds.nobs
    Y = np.column_stack([ds.get(v).data.astype(float) for v in deps])
    Xf = []
    for i in range(len(deps)):
        cols = [_term(ds, x) for x in xnames[i]]
        if cons[i]:
            cols.append(np.ones(n))
        Xf.append(np.column_stack(cols) if cols else np.zeros((n, 0)))
    ok = mask & np.all(Y < SYS, axis=1)
    for Xi in Xf:
        if Xi.size:
            ok &= np.all(Xi < SYS, axis=1)
    t0, t1 = _check_gaps(s, ok)
    order = np.argsort(ds.get(info.tvar).data[ok], kind="stable")
    y = Y[ok][order]
    X = [Xi[ok][order] for Xi in Xf]
    T, m = y.shape
    P = layout(spec)
    names = [q[0] for q in P]

    def contrib(b):
        try:
            par = _unpack(spec, b, P)
            r = mgarch_loglik(spec, par, residuals(spec, par, y, X))
        except (np.linalg.LinAlgError, ValueError, OverflowError):
            return None
        return None if r is None or not np.all(np.isfinite(r)) else r

    def ll(b):
        c_ = contrib(b)
        return -np.inf if c_ is None else float(c_.sum())

    b0 = start_values(spec, P, y, X)
    if not np.isfinite(ll(b0)):
        raise StataError(430, "initial values not feasible")
    r = maximize_switching(ll, lambda bb: _scores(contrib, bb), b0, **ml_options(o))
    out = s.output
    if not o.get("nolog"):
        out.write("\n", "text")
        print_log(out, r, "log likelihood")
    b = r.b
    _l, _g, Hm = _numderiv(ll, b)
    try:
        V = np.linalg.inv(-Hm)
    except np.linalg.LinAlgError:
        V = np.linalg.pinv(-Hm)
    V = (V + V.T) / 2
    se = np.sqrt(np.maximum(np.diag(V), 0))
    meanidx = [k for k, q in enumerate(P) if q[3] == "mean" and q[2] != "_cons"]
    chi2 = SYS
    if meanidx:
        try:
            chi2 = float(b[meanidx] @ np.linalg.solve(V[np.ix_(meanidx, meanidx)], b[meanidx]))
        except np.linalg.LinAlgError:
            pass
    pv = float(chi2_tail(chi2, len(meanidx))) if chi2 < SYS else SYS
    rows: list[CoefRow] = []
    last = None
    lastgrp = None
    for k, q in enumerate(P):
        name, eq, lab, grp = q
        if eq and eq != last:
            if last is not None or rows:
                rows.append(CoefRow("sep"))
            rows.append(CoefRow("eq", eq))
            last, lastgrp = eq, None
        elif not eq and last is not None:
            rows.append(CoefRow("sep"))
            last = None
        if isinstance(lab, tuple):
            if lab[0] != lastgrp:
                if lastgrp is not None:
                    pass
                rows.append(CoefRow("header", lab[0]))
                lastgrp = lab[0]
            rows.append(CoefRow("coef", f"L{lab[1]}.", b[k], se[k], name=name))
        else:
            rows.append(CoefRow("coef", lab, b[k], se[k], name=name))
    est = Estimates("mgarch", " ".join(deps), names, b, V, T, "z", None)
    est.eqnames = [q[1] or "_" for q in P]
    est.rows = rows
    est.level = float(o["level"]) if isinstance(o.get("level"), str) else float(s.settings.get("level", "95"))
    est.vcetype = "OIM"
    est.extra = {"N": float(T), "ll": r.ll, "chi2": chi2, "p": pv, "df_m": len(meanidx), "title": TITLES[kind],
                 "sample": _sample_text(s, t0, t1)}
    est.mgarch = {"spec": spec, "y": y, "X": X, "mask": ok.copy(), "order": order, "layout": P}
    est.display = display_mgarch
    scalars = [("N", float(T)), ("k", float(len(b))), ("k_aux", 0.0), ("k_eq", float(m + (kind != "dvech") + (kind == "dcc"))),
               ("k_dv", float(m)), ("df_m", float(len(meanidx))), ("ll", r.ll), ("chi2", chi2), ("p", pv),
               ("tmin", t0), ("tmax", t1), ("rank", float(len(b))), ("ic", float(r.iterations)), ("rc", 0.0),
               ("converged", 1.0 if r.converged else 0.0)]
    macros = [("cmd", "mgarch"), ("cmdline", cmdline("mgarch", args)), ("model", kind), ("depvar", " ".join(deps)),
              ("title", TITLES[kind]), ("tvar", info.tvar), ("distribution", "gaussian"),
              ("arch", " ".join(map(str, arch))), ("garch", " ".join(map(str, garch))),
              ("vce", "oim"), ("vcetype", "OIM"), ("chi2type", "Wald"), ("properties", "b V")]
    post(s, est, scalars, macros, ok)
    display_mgarch(s, est)


def display_mgarch(s: "Session", est: Estimates, *, header: bool = True, table: bool = True,
                   level: float | None = None, **_kw) -> None:
    from .tsmodels import _g, _right
    out = s.output
    x = est.extra
    if header:
        out.write(f"\n{x['title']}\n\n", "text")
        _right(out, x["sample"], "Number of obs", f"{int(x['N']):,}")
        dfm = f"({int(x['df_m'])})" if x["df_m"] else "(.)"
        _right(out, "Distribution: Gaussian", f"Wald chi2{dfm}", _g(x["chi2"], "%10.2f"))
        _right(out, f"Log likelihood = {_g(x['ll'], '%9.0g')}", "Prob > chi2", _g(x["p"], "%10.4f"))
        out.write("\n", "text")
    if table:
        coef_table(s, est.depvar, est.rows, stat="z", df=None, level=level or est.level, vcetype=est.vcetype)
