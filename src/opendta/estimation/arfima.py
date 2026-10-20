"""arfima: autoregressive fractionally integrated moving-average models, exact Gaussian maximum
likelihood ([TS] arfima). The autocovariances are those of an ARMA process driven by fractionally
integrated noise (the convolution of the ARMA autocovariances with those of the fractional noise);
the likelihood uses the Durbin-Levinson recursions."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np
from scipy import signal, special

from ..core import missing as M
from ..core.errors import StataError
from ..core.stats import chi2_tail
from ..lang.syntax import find_top, match_options, parse_standard
from . import fvars
from .ml import maximize, options as ml_options, print_log
from .results import CoefRow, Estimates, coef_table, post, rows_from_columns
from .tsops import require_ts

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def fracnoise_acov(d: float, n: int) -> np.ndarray:
    """Autocovariances at lags 0..n-1 of fractionally integrated unit-variance white noise:
    g(0) = G(1-2d)/G(1-d)^2 and g(k) = g(k-1) (k-1+d)/(k-d)."""
    g = np.empty(n)
    g[0] = math.exp(special.gammaln(1 - 2 * d) - 2 * special.gammaln(1 - d))
    for k in range(1, n):
        g[k] = g[k - 1] * (k - 1 + d) / (k - d)
    return g


def arma_acov(ar: np.ndarray, ma: np.ndarray, nlags: int) -> np.ndarray:
    """Autocovariances at lags 0..nlags of the unit-variance ARMA(p,q) process, from its psi weights."""
    a = np.concatenate([[1.0], -np.asarray(ar, dtype=float)])
    b = np.concatenate([[1.0], np.asarray(ma, dtype=float)])
    L = 4000
    imp = np.zeros(L)
    imp[0] = 1.0
    psi = signal.lfilter(b, a, imp)
    out = np.empty(nlags + 1)
    for m in range(nlags + 1):
        out[m] = float(psi[:L - m] @ psi[m:]) if m < L else 0.0
    return out


def arfima_acov(d: float, ar: np.ndarray, ma: np.ndarray, sigma2: float, n: int) -> np.ndarray:
    """Autocovariances of the ARFIMA(p,d,q) process at lags 0..n-1."""
    if d == 0:
        return sigma2 * arma_acov(ar, ma, n - 1)
    M_ = min(max(60, 20 * (len(ar) + len(ma)) + 40), 400)
    ga = arma_acov(ar, ma, M_)
    gsym = np.concatenate([ga[:0:-1], ga])                # lags -M..M
    gx = fracnoise_acov(d, n + M_)
    gxs = np.concatenate([gx[:0:-1], gx])                 # lags -(n+M-1)..(n+M-1)
    full = np.convolve(gxs, gsym)                          # lag index offset (n+M-1) + M
    center = (n + M_ - 1) + M_
    return sigma2 * full[center:center + n]


def toeplitz_loglik(y: np.ndarray, g: np.ndarray) -> float:
    """Gaussian log likelihood of y with autocovariances g (Durbin-Levinson); -inf if not positive
    definite."""
    n = len(y)
    v = g[0]
    if v <= 0:
        return -np.inf
    ll = -0.5 * (math.log(2 * math.pi * v) + y[0] ** 2 / v)
    phi = np.zeros(0)
    for t in range(1, n):
        if t == 1:
            k = g[1] / v
            phi = np.array([k])
        else:
            k = (g[t] - phi @ g[t - 1:0:-1]) / v
            phi = np.concatenate([phi - k * phi[::-1], [k]])
        v = v * (1 - k * k)
        if v <= 0 or not np.isfinite(v):
            return -np.inf
        e = y[t] - float(phi @ y[t - 1::-1][:t])
        ll += -0.5 * (math.log(2 * math.pi * v) + e * e / v)
    return ll


_OPTS = {"ar": 2, "ma": 2, "noconstant": 4, "mle": 3, "mpl": 3, "smemory": 2, "vce": 3, "robust": 1,
         "level": 1, "nolog": 3, "iterate": 4, "tolerance": 3, "ltolerance": 4, "nrtolerance": 4,
         "technique": 4, "constraints": 3, "from": 3, "collinear": 4}


def cmd_arfima(s: "Session", args: str) -> None:
    from .postest import replay
    from .regress import Sample, cmdline
    from .tsmodels import _check_gaps, _lags, _numderiv, _sample_text
    if not args.strip() or args.strip().startswith(","):
        replay(s, "arfima", args)
        return
    c = find_top(args, ",")
    body, opts = (args, "") if c == -1 else (args[:c], args[c + 1:])
    o = match_options(opts, _OPTS) if opts.strip() else {}
    for bad in ("mpl", "smemory", "constraints", "from"):
        if o.get(bad):
            raise StataError(198, f"option {bad} not implemented in OpenDTA")   # OpenDTA
    ar, ma = _lags(o.get("ar")), _lags(o.get("ma"))
    p = parse_standard(body)
    extra = " ".join(f"{k}({v})" if isinstance(v, str) else k for k, v in o.items() if k in ("noconstant", "level"))
    text = p.varlist + (f" if {p.if_}" if p.if_ else "") + (f" in {p.in_}" if p.in_ else "")
    smp = Sample(s, text + (f", {extra}" if extra else ""), opts_spec={"noconstant": 4}, weights=())
    t0, t1 = _check_gaps(s, smp.mask)
    fvars.drop_collinear(smp.cols, smp.constant)
    est_cols = [col for col in smp.cols if col.values is not None]
    n = smp.n
    k = len(est_cols) + (1 if smp.constant else 0)
    info = require_ts(s.data)
    order = np.argsort(s.data.get(info.tvar).data[smp.mask], kind="stable")
    X = np.column_stack([col.values for col in est_cols] + ([np.ones(n)] if smp.constant else [])) \
        if k else np.zeros((n, 0))
    X, y = X[order], smp.y.astype(float)[order]
    na, nm = len(ar), len(ma)

    def parts(b):
        phi = np.zeros(max(ar + [0]))
        for lag, v in zip(ar, b[k:k + na]):
            phi[lag - 1] = v
        th = np.zeros(max(ma + [0]))
        for lag, v in zip(ma, b[k + na:k + na + nm]):
            th[lag - 1] = v
        return phi, th, b[k + na + nm], b[-1]

    def ll(b):
        phi, th, d, s2 = parts(b)
        if not -0.5 < d < 0.5 or s2 <= 0:
            return -np.inf
        if len(phi):
            comp = np.zeros((len(phi), len(phi)))
            comp[0] = phi
            if len(phi) > 1:
                comp[1:, :-1] = np.eye(len(phi) - 1)
            if np.max(np.abs(np.linalg.eigvals(comp))) >= 0.9999:
                return -np.inf
        try:
            return toeplitz_loglik(y - X @ b[:k], arfima_acov(d, phi, th, s2, n))
        except (FloatingPointError, np.linalg.LinAlgError, ValueError):
            return -np.inf

    def fun(b, todo):
        if todo == 0:
            return ll(b), None, None
        return _numderiv(ll, b)

    beta0 = np.linalg.lstsq(X, y, rcond=None)[0] if k else np.zeros(0)
    var0 = float(np.var(y - X @ beta0)) or 1.0
    b0 = np.concatenate([beta0, np.zeros(na + nm), [0.1], [var0]])
    r = maximize(fun, b0, **ml_options(o))
    out = s.output
    if not o.get("nolog"):
        out.write("\n", "text")
        print_log(out, r, "log likelihood")
    b = r.b
    _l, _g_, H = _numderiv(ll, b)
    try:
        V = np.linalg.inv(-H)
    except np.linalg.LinAlgError:
        V = np.linalg.pinv(-H)
    if o.get("vce") not in (None, "oim") or o.get("robust"):
        raise StataError(198, "only vce(oim) is implemented in OpenDTA")   # OpenDTA
    V = (V + V.T) / 2
    se = np.sqrt(np.maximum(np.diag(V), 0))
    idx = list(range(len(est_cols))) + list(range(k, k + na + nm + 1))     # the constant is not tested
    try:
        chi2 = float(b[idx] @ np.linalg.solve(V[np.ix_(idx, idx)], b[idx])) if idx else SYS
    except np.linalg.LinAlgError:
        chi2 = SYS
    pv = float(chi2_tail(chi2, len(idx))) if idx and chi2 < SYS else SYS
    eqname = smp.depname.split(".")[-1]
    bmap = {col.name: b[i] for i, col in enumerate(est_cols)}
    semap = {col.name: se[i] for i, col in enumerate(est_cols)}
    if smp.constant:
        bmap["_cons"], semap["_cons"] = b[len(est_cols)], se[len(est_cols)]
    rows, names, eqs = [], [], []
    if k:
        rows.append(CoefRow("eq", eqname))
        rows += rows_from_columns([c_ for c_ in smp.cols if c_.values is not None], bmap, semap, s.data,
                                  constant=smp.constant)
        names += [c_.name for c_ in est_cols] + (["_cons"] if smp.constant else [])
        eqs += [eqname] * k
        rows.append(CoefRow("sep"))
    rows.append(CoefRow("eq", "ARFIMA"))
    j = k
    for lab, lags in (("ar", ar), ("ma", ma)):
        if not lags:
            continue
        if rows[-1].kind != "eq":
            rows.append(CoefRow("blank"))
        rows.append(CoefRow("header", lab))
        for lg in lags:
            rows.append(CoefRow("coef", f"L{lg}.", b[j], se[j], name=f"ARFIMA:L{lg}.{lab}"))
            names.append(f"L{lg}.{lab}")
            eqs.append("ARFIMA")
            j += 1
    if rows[-1].kind != "eq":
        rows.append(CoefRow("blank"))
    rows.append(CoefRow("coef", "d", b[j], se[j], name="ARFIMA:d"))
    names.append("d"); eqs.append("ARFIMA")
    rows.append(CoefRow("sep"))
    rows.append(CoefRow("coef", "/sigma2", b[-1], se[-1], name="sigma2:_cons", zero_lo=True))
    names.append("_cons"); eqs.append("sigma2")
    kk = b.size
    est = Estimates("arfima", smp.depname, names, b, V, n, "z", None)
    est.eqnames = eqs
    est.rows = rows
    est.cols, est.constant, est.terms, est.depcomp = smp.cols, smp.constant, smp.terms, smp.depcomp
    est.level = float(o["level"]) if o.get("level") else smp.level
    est.vcetype = "OIM"
    est.extra = {"N": float(n), "ll": r.ll, "chi2": chi2, "p": pv, "df_m": len(idx), "sample": _sample_text(s, t0, t1)}
    phi, th, d, s2 = parts(b)
    g = arfima_acov(d, phi, th, s2, 2)
    est.arfima = {"ar": {lg: float(phi[lg - 1]) for lg in ar}, "ma": {lg: float(th[lg - 1]) for lg in ma},
                  "sigma2": float(s2), "d": float(d), "gamma0": float(g[0])}
    est.display = display_arfima
    scalars = [("N", float(n)), ("k", float(kk)), ("k_eq", float(len(dict.fromkeys(eqs)))), ("k_dv", 1.0),
               ("df_m", float(len(idx))), ("ll", r.ll), ("chi2", chi2), ("p", pv), ("tmin", t0), ("tmax", t1),
               ("rank", float(kk)), ("ic", float(r.iterations)), ("rc", 0.0), ("converged", 1.0 if r.converged else 0.0)]
    macros = [("cmd", "arfima"), ("cmdline", cmdline("arfima", args)), ("depvar", smp.depname),
              ("eqnames", " ".join(dict.fromkeys(eqs))), ("title", "ARFIMA regression"),
              ("ar", " ".join(map(str, ar))), ("ma", " ".join(map(str, ma))), ("vce", "oim"),
              ("vcetype", "OIM"), ("chi2type", "Wald"), ("tvar", info.tvar), ("properties", "b V")]
    post(s, est, scalars, macros, smp.mask)
    display_arfima(s, est)


def display_arfima(s: "Session", est: Estimates, *, header: bool = True, table: bool = True,
                   level: float | None = None, **_kw) -> None:
    from .tsmodels import _g, _right
    out = s.output
    x = est.extra
    if header:
        out.write("\nARFIMA regression\n\n", "text")
        _right(out, x["sample"].replace("Sample:  ", "Sample: "), "Number of obs", f"{int(x['N']):,}")
        _right(out, "", f"Wald chi2({int(x['df_m'])})", _g(x["chi2"], "%10.2f"))
        _right(out, f"Log likelihood = {_g(x['ll'], '%10.0g'):>10}", "Prob > chi2", _g(x["p"], "%10.4f"))
    if table:
        coef_table(s, est.depvar, est.rows, stat="z", df=None, level=level or est.level, vcetype=est.vcetype)
        out.write("Note: The test of the variance against zero is one sided, and the two-sided\n"
                  "      confidence interval is truncated at zero.\n", "text")
