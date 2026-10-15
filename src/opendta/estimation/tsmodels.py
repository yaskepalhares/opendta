"""Time-series models and tests: arima, var, varsoc, vargranger, dfuller,
corrgram and newey ([TS] arima, [TS] var, [TS] varsoc, [TS] vargranger,
[TS] dfuller, [TS] corrgram, [TS] newey).

arima: exact maximum likelihood with the Kalman filter (Harvey's state
space form, stationary initial covariance), ar() and ma() lag lists or
arima(p,d,q); OPG variance by default, as in Stata. Stata optimizes with
BHHH and then BFGS (technique(bhhh 5 bfgs 10)), as OpenDTA does; the step
sizes are undocumented, so the iteration values may differ, but the maximum
is the same (VERIFY: compat/do/1101_ts.do).

var: OLS equation by equation with the ML covariance of the residuals.
dfuller: Fuller (1976) critical values interpolated in the number of
observations, as in Stata, and MacKinnon (1994) approximate p-values.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np
from scipy import linalg, stats as st

from ..core import missing as M
from ..core.errors import StataError
from ..core.formats import format_value
from ..core.stats import chi2_tail
from ..lang.syntax import find_top, match_options
from . import fvars
from .ml import maximize_switching, print_log
from .results import CoefRow, Estimates, coef_table, post, rows_from_columns
from .tsops import _fmt_time, require_ts

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def _g(x: float, fmt: str = "%9.0g") -> str:
    if x is None or not np.isfinite(x) or abs(x) >= SYS:
        return "."
    return format_value(float(x), fmt, pad=False).strip()


def _right(out, left: str, lab: str, val: str) -> None:
    out.write(f"{left:<48}{lab:<18}= ", "text")
    out.write(f"{val:>10}\n", "result")


def _check_gaps(s: "Session", mask: np.ndarray) -> tuple[float, float]:
    """The sample must be a run of consecutive periods; returns (first, last) time."""
    info = require_ts(s.data)
    if info.panel:
        raise StataError(498, "panel data not allowed; use a single time series")   # VERIFY
    t = s.data.get(info.tvar).data[mask].astype(float)
    if t.size == 0:
        raise StataError(2000, "no observations")
    ts = np.sort(t)
    if np.any(np.abs(np.diff(ts) - info.delta) > 1e-9):
        raise StataError(498, "sample may not contain gaps")   # VERIFY
    return float(ts[0]), float(ts[-1])


def _sample_text(s: "Session", t0: float, t1: float) -> str:
    info = require_ts(s.data)
    return f"Sample:  {_fmt_time(s.data, info.tvar, t0)} - {_fmt_time(s.data, info.tvar, t1)}"


def _numderiv(f, x, h=None):
    from .mixed import _numderiv as nd
    return nd(f, x)


# ---------------------------------------------------------------------------
# arima
# ---------------------------------------------------------------------------

def _scores(contrib, b: np.ndarray) -> np.ndarray:
    """Observation scores by central differences of the log-likelihood contributions."""
    kk = b.size
    c0 = contrib(b)
    S = np.zeros((len(c0) if c0 is not None else 0, kk))
    for j in range(kk):
        h = 1e-6 * (abs(b[j]) + 1)
        e = np.zeros(kk)
        e[j] = h
        cp, cm = contrib(b + e), contrib(b - e)
        if cp is None or cm is None:
            cp = contrib(b + e / 100)
            cm = contrib(b - e / 100)
            h = h / 100
        S[:, j] = (cp - cm) / (2 * h)
    return S


def _lags(text) -> list[int]:
    if not isinstance(text, str):
        return []
    from ..lang.words import parse_numlist
    return sorted({int(v) for v in parse_numlist(text)})


class _ARMA:
    def __init__(self, ar: list[int], ma: list[int]):
        self.ar, self.ma = ar, ma
        self.P = max(ar) if ar else 0
        self.Q = max(ma) if ma else 0
        self.r = max(self.P, self.Q + 1)

    def contributions(self, e: np.ndarray, phi: np.ndarray, theta: np.ndarray, sigma: float):
        """Log likelihood of each observation by the Kalman filter (None if not stationary)."""
        r = self.r
        if sigma <= 0:
            return None
        T = np.zeros((r, r))
        for lag, c in zip(self.ar, phi):
            T[lag - 1, 0] = c
        if r > 1:
            T[:-1, 1:] += np.eye(r - 1)
        R = np.zeros(r)
        R[0] = 1.0
        for lag, c in zip(self.ma, theta):
            R[lag] = c
        if self.P and np.max(np.abs(np.linalg.eigvals(T))) >= 1:
            return None
        RR = np.outer(R, R) * sigma ** 2
        try:
            Pm = linalg.solve_discrete_lyapunov(T, RR)
        except (np.linalg.LinAlgError, ValueError):
            return None
        a = np.zeros(r)
        out = np.empty(len(e))
        pred = np.empty(len(e))
        for t, y in enumerate(e):
            F = Pm[0, 0]
            if F <= 0:
                return None
            v = y - a[0]
            pred[t] = a[0]
            out[t] = -0.5 * (math.log(2 * math.pi) + math.log(F) + v * v / F)
            K = T @ Pm[:, 0] / F
            a = T @ a + K * v
            Pm = T @ Pm @ T.T - np.outer(K, K) * F + RR
        self.pred = pred
        return out


def _hr_start(e: np.ndarray, ar: list[int], ma: list[int]) -> np.ndarray:
    """Starting values of the ARMA coefficients by the Hannan-Rissanen regressions: a long
    autoregression gives the innovations, then e is regressed on its lags and the lagged
    innovations. VERIFY: Stata's starting values are undocumented."""
    na, nm = len(ar), len(ma)
    if not na and not nm:
        return np.zeros(0)
    n = len(e)
    m = min(max(8, (max(ar) if ar else 0) + (max(ma) if ma else 0)), max(n // 4, 1))
    a = np.zeros(n)
    if nm:
        Z = np.column_stack([e[m - j:n - j] for j in range(1, m + 1)])
        c, *_ = np.linalg.lstsq(Z, e[m:], rcond=None)
        a[m:] = e[m:] - Z @ c
    start = max([m if nm else 0] + [lg for lg in ar] + [m + lg for lg in ma])
    if start >= n - (na + nm) - 2:
        return np.zeros(na + nm)
    cols = [e[start - lg:n - lg] for lg in ar] + [a[start - lg:n - lg] for lg in ma]
    c, *_ = np.linalg.lstsq(np.column_stack(cols), e[start:], rcond=None)
    # keep the start inside the stationary and invertible region
    return np.clip(c, -0.95, 0.95)


def cmd_arima(s: "Session", args: str) -> None:
    from .postest import replay
    from .regress import Sample, cmdline
    if not args.strip() or args.strip().startswith(","):
        replay(s, "arima", args)
        return
    c = find_top(args, ",")
    body, opts = (args, "") if c == -1 else (args[:c], args[c + 1:])
    spec = {"ar": 2, "ma": 2, "arima": 5, "sarima": 6, "noconstant": 4, "vce": 3, "nolog": 4, "log": 3,
            "iterate": 4, "level": 1, "condition": 4, "diffuse": 4, "technique": 4, "tolerance": 3,
            "ltolerance": 4, "nrtolerance": 4, "robust": 1}
    o = match_options(opts, spec) if opts.strip() else {}
    if o.get("sarima") or o.get("condition") or o.get("diffuse"):
        raise StataError(198, "sarima(), condition and diffuse are not available in OpenDTA")   # OpenDTA
    ar, ma = _lags(o.get("ar")), _lags(o.get("ma"))
    from ..lang.syntax import parse_standard
    p = parse_standard(body)
    toks = p.varlist.split()
    if isinstance(o.get("arima"), str):
        nums = [int(float(v)) for v in o["arima"].replace(",", " ").split()]
        if len(nums) != 3:
            raise StataError(198, "arima() invalid")   # VERIFY
        pp, d, q = nums
        ar, ma = list(range(1, pp + 1)), list(range(1, q + 1))
        if d:
            pre = "D." if d == 1 else f"D{d}."
            toks = [pre + t if not t.startswith(("i.", "c.")) else t for t in toks]
    rest_opts = " ".join(f"{k}({v})" if isinstance(v, str) else k for k, v in o.items()
                         if k in ("noconstant", "level", "vce", "robust"))
    text = " ".join(toks)
    if p.if_:
        text += f" if {p.if_}"
    if p.in_:
        text += f" in {p.in_}"
    smp = Sample(s, text + (f", {rest_opts}" if rest_opts else ""),
                 opts_spec={"noconstant": 4}, weights=())
    t0, t1 = _check_gaps(s, smp.mask)
    fvars.drop_collinear(smp.cols, smp.constant)
    est_cols = [col for col in smp.cols if col.values is not None]
    n = smp.n
    X = np.column_stack([col.values for col in est_cols] + ([np.ones(n)] if smp.constant else [])) \
        if (est_cols or smp.constant) else np.zeros((n, 0))
    info = require_ts(s.data)
    order = np.argsort(s.data.get(info.tvar).data[smp.mask], kind="stable")
    X = X[order]
    y = smp.y.astype(float)[order]
    k = X.shape[1]
    arma = _ARMA(ar, ma)
    na, nm = len(ar), len(ma)

    def unpack(b):
        return b[:k], b[k:k + na], b[k + na:k + na + nm], b[-1]

    def contrib(b):
        beta, phi, theta, sigma = unpack(b)
        return arma.contributions(y - X @ beta, phi, theta, sigma)

    def ll(b):
        c_ = contrib(b)
        return -np.inf if c_ is None else float(c_.sum())

    def fun(b, todo):
        if todo == 0:
            return ll(b), None, None
        return _numderiv(ll, b)

    beta0 = np.linalg.lstsq(X, y, rcond=None)[0] if k else np.zeros(0)
    res0 = y - X @ beta0
    b0 = np.concatenate([beta0, _hr_start(res0, ar, ma), [float(np.std(res0)) or 1.0]])
    if not np.isfinite(ll(b0)):
        b0 = np.concatenate([beta0, np.zeros(na + nm), [float(np.std(res0)) or 1.0]])
    kw = dict(maxiter=int(o["iterate"]) if o.get("iterate") else 16000,
              tol=float(o["tolerance"]) if o.get("tolerance") else 1e-6,
              ltol=float(o["ltolerance"]) if o.get("ltolerance") else 1e-7,
              nrtol=float(o["nrtolerance"]) if o.get("nrtolerance") else 1e-5)
    r = maximize_switching(ll, lambda bb: _scores(contrib, bb), b0, **kw)
    # Newton steps outside the log take the estimates to the maximum (the BFGS approximation
    # stops up to 1e-4 away from it); the variance needs the Hessian there
    _llm, gm, Hm = _numderiv(ll, r.b)
    for _ in range(5):
        try:
            bn = r.b - np.linalg.solve(Hm, gm)
        except np.linalg.LinAlgError:
            break
        lln = ll(bn)
        if not np.isfinite(lln) or lln < r.ll:
            break
        r.b, r.ll = bn, lln
        _llm, gm, Hm = _numderiv(ll, r.b)
        if abs(float(gm @ np.linalg.solve(Hm, gm))) < 1e-14:
            break
    r.H = Hm
    out = s.output
    vce = str(o.get("vce") or ("robust" if o.get("robust") else "opg")).strip().lower()
    vce = "robust" if vce.startswith("r") else ("oim" if vce == "oim" else "opg")
    if not o.get("nolog"):
        out.write("\n", "text")
        print_log(out, r, "log likelihood")
    b = r.b
    kk = b.size
    S = _scores(contrib, b)
    try:
        Hinv = np.linalg.inv(-r.H)
    except np.linalg.LinAlgError:
        Hinv = np.linalg.pinv(-r.H)
    if vce == "opg":
        V = np.linalg.inv(S.T @ S)
    elif vce == "oim":
        V = Hinv
    else:
        V = Hinv @ (S.T @ S) @ Hinv
    V = (V + V.T) / 2
    # Wald test: everything except the constant and sigma
    idx = [i for i in range(kk - 1) if not (smp.constant and i == k - 1)]
    try:
        chi2 = float(b[idx] @ np.linalg.solve(V[np.ix_(idx, idx)], b[idx])) if idx else SYS
    except np.linalg.LinAlgError:
        chi2 = SYS
    pv = float(chi2_tail(chi2, len(idx))) if idx and chi2 < SYS else SYS
    base = smp.depcomp.name if hasattr(smp.depcomp, "name") else smp.depname.split(".")[-1]
    eqname = smp.depname.split(".")[-1]
    names = [col.name for col in est_cols] + (["_cons"] if smp.constant else []) + \
        [f"L{l}.ar" if l > 1 else "L.ar" for l in ar] + [f"L{l}.ma" if l > 1 else "L.ma" for l in ma] + ["_cons"]
    eqs = [eqname] * k + ["ARMA"] * (na + nm) + ["sigma"]
    se = np.sqrt(np.maximum(np.diag(V), 0))
    bmap = {col.name: b[i] for i, col in enumerate(est_cols)}
    semap = {col.name: se[i] for i, col in enumerate(est_cols)}
    if smp.constant:
        bmap["_cons"], semap["_cons"] = b[k - 1], se[k - 1]
    rows = []
    if k:
        rows.append(CoefRow("eq", eqname))
        rows += rows_from_columns([col for col in smp.cols if col.values is not None], bmap, semap, s.data,
                                  constant=smp.constant)
        rows.append(CoefRow("sep"))
    if na + nm:
        rows.append(CoefRow("eq", "ARMA"))
        j = k
        for lab, lags in (("ar", ar), ("ma", ma)):
            if not lags:
                continue
            if rows[-1].kind != "eq":
                rows.append(CoefRow("blank"))
            rows.append(CoefRow("header", lab))
            for l in lags:
                rows.append(CoefRow("coef", f"L{l}.", b[j], se[j], name=f"ARMA:L{l}.{lab}"))
                j += 1
        rows.append(CoefRow("sep"))
    rows.append(CoefRow("coef", "/sigma", b[-1], se[-1], name="sigma:_cons"))
    est = Estimates("arima", smp.depname, names, b, V, n, "z", None)
    est.eqnames = eqs
    est.rows = rows
    est.cols, est.constant, est.terms, est.depcomp = smp.cols, smp.constant, smp.terms, smp.depcomp
    est.level = float(o["level"]) if o.get("level") else smp.level
    est.vcetype = {"opg": "OPG", "robust": "Semirobust", "oim": "OIM"}[vce]
    est.robust = vce == "robust"
    est.extra = {"N": float(n), "ll": r.ll, "chi2": chi2, "p": pv, "df_m": len(idx),
                 "sample": _sample_text(s, t0, t1)}
    est.arima = {"arma": arma, "k": k, "ar": ar, "ma": ma, "order": order, "mask": smp.mask.copy()}
    est.display = display_arima
    est.predict = _predict_arima
    est.predict_opts = {"y": 1, "mse": 3, "structural": 6}
    scalars = [("N", float(n)), ("k", float(kk)), ("k_eq", float(1 + (na + nm > 0) + 1)),
               ("k_dv", 1.0), ("ll", r.ll), ("chi2", chi2), ("p", pv), ("df_m", float(len(idx))),
               ("sigma", float(b[-1])), ("tmin", t0), ("tmax", t1), ("ar_max", float(arma.P)),
               ("ma_max", float(arma.Q)), ("rank", float(kk)), ("ic", float(r.iterations)),
               ("converged", 1.0 if r.converged else 0.0)]
    macros = [("cmd", "arima"), ("cmdline", cmdline("arima", args)), ("depvar", smp.depname),
              ("eqnames", " ".join(dict.fromkeys(eqs))), ("title", "ARIMA regression"),
              ("ar", " ".join(map(str, ar))), ("ma", " ".join(map(str, ma))), ("vce", vce),
              ("vcetype", est.vcetype), ("chi2type", "Wald"), ("tvar", info.tvar),
              ("predict", "arima_p"), ("properties", "b V")]
    post(s, est, scalars, macros, smp.mask)
    del base
    display_arima(s, est)


def display_arima(s: "Session", est: Estimates, *, header: bool = True, table: bool = True,
                  level: float | None = None, **_kw) -> None:
    out = s.output
    x = est.extra
    if header:
        out.write("\nARIMA regression\n\n", "text")
        _right(out, x["sample"], "Number of obs", f"{int(x['N']):,}")
        _right(out, "", f"Wald chi2({int(x['df_m'])})", _g(x["chi2"], "%10.2f"))
        _right(out, f"Log likelihood = {_g(x['ll'], '%9.0g')}", "Prob > chi2", _g(x["p"], "%10.4f"))
        out.write("\n", "text")
    if table:
        coef_table(s, est.depvar, est.rows, stat="z", df=None, level=level or est.level,
                   vcetype=est.vcetype)
        out.write("Note: The test of the variance against zero is one sided, and the two-sided\n"
                  "      confidence interval is truncated at zero.\n", "text")


def _predict_arima(s: "Session", est: Estimates, stat: str, mask: np.ndarray):
    """One-step predictions of the (differenced) dependent variable over the estimation sample."""
    A = est.arima
    ds = s.data
    k = A["k"]
    b = est.b
    rows = np.nonzero(A["mask"])[0][A["order"]]
    from .postest import xb_all
    xb_full, ok = xb_all(s, _fixed_part(est, k))
    yv = fvars.comp_values(ds, est.depcomp).astype(float)
    e = yv[rows] - xb_full[rows]
    na, nm = len(A["ar"]), len(A["ma"])
    arma = A["arma"]
    c = arma.contributions(e, b[k:k + na], b[k + na:k + na + nm], b[-1])
    if c is None:
        raise StataError(498, "prediction failed")   # VERIFY
    pred = np.full(ds.nobs, SYS)
    pred[rows] = xb_full[rows] + arma.pred
    if stat in ("", "xb"):
        note = "" if stat else "(option xb assumed; predicted values)"   # VERIFY
        if stat == "" and False:
            pass
        return pred, note
    if stat == "residuals":
        r = np.full(ds.nobs, SYS)
        r[rows] = yv[rows] - pred[rows]
        return r, ""
    if stat == "structural":
        return np.where(ok, xb_full, SYS), ""
    if stat == "y":
        name = est.depvar
        if not name.startswith("D"):
            return pred, ""
        from .tsops import shift
        base = name.split(".", 1)[1]
        lagged = shift(ds, ds.get(base).data.astype(float), 1)
        return np.where((pred < SYS) & (lagged < SYS), pred + lagged, SYS), ""
    raise StataError(198, f"option {stat} not allowed")


def _fixed_part(est: Estimates, k: int) -> Estimates:
    e2 = Estimates(est.cmd, est.depvar, est.names[:k], est.b[:k], est.V[:k, :k], est.N)
    e2.cols, e2.constant, e2.terms = est.cols, est.constant, est.terms
    return e2


# ---------------------------------------------------------------------------
# var
# ---------------------------------------------------------------------------

def _var_fit(s: "Session", deps: list[str], lags: list[int], exog: list[str], constant: bool,
             mask: np.ndarray | None = None):
    """Returns the design for all equations over the common sample."""
    from .tsops import shift
    ds = s.data
    n = ds.nobs
    Y = np.column_stack([ds.get(v).data.astype(float) for v in deps])
    cols, names = [], []
    for dv in deps:
        x = ds.get(dv).data.astype(float)
        for l in lags:
            cols.append(shift(ds, x, l))
            names.append((dv, l))
    for ev in exog:
        cols.append(ds.get(ev).data.astype(float))
        names.append((ev, 0))
    X = np.column_stack(cols) if cols else np.zeros((n, 0))
    ok = np.all(Y < SYS, axis=1) & (np.all(X < SYS, axis=1) if X.size else True)
    if mask is not None:
        ok &= mask
    if constant:
        X = np.column_stack([X, np.ones(n)])
        names.append(("_cons", 0))
    return Y, X, ok, names


def _var_stats(Y, X):
    T = Y.shape[0]
    B, *_ = np.linalg.lstsq(X, Y, rcond=None)
    U = Y - X @ B
    Sml = U.T @ U / T
    K = Y.shape[1]
    m = X.shape[1]
    _s, ld = np.linalg.slogdet(Sml)
    ll = -0.5 * T * (K * math.log(2 * math.pi) + ld + K)
    tp = K * m
    aic = -2 * ll / T + 2 * tp / T
    hqic = -2 * ll / T + 2 * math.log(math.log(T)) * tp / T
    sbic = -2 * ll / T + math.log(T) * tp / T
    fpe = math.exp(ld) * ((T + m) / (T - m)) ** K
    return B, U, Sml, ll, aic, hqic, sbic, fpe, math.exp(ld)


def cmd_var(s: "Session", args: str) -> None:
    from .postest import replay
    from ..commands._util import touse
    from ..lang.syntax import parse_standard
    from .regress import cmdline
    if not args.strip() or args.strip().startswith(","):
        replay(s, "var", args)
        return
    p = parse_standard(args)
    o = match_options(p.options, {"lags": 3, "exog": 3, "noconstant": 4, "dfk": 3, "small": 2,
                                  "level": 1, "lutstats": 3, "noisure": 5, "constraints": 6}) \
        if p.options.strip() else {}
    deps = s.expand_varlist(p.varlist)
    lags = _lags(o.get("lags")) if o.get("lags") else [1, 2]
    exog = s.expand_varlist(o["exog"]) if isinstance(o.get("exog"), str) else []
    constant = not o.get("noconstant")
    Y, X, ok, xnames = _var_fit(s, deps, lags, exog, constant, touse(s, p))
    t0, t1 = _check_gaps(s, ok)
    Yo, Xo = Y[ok], X[ok]
    T, K = Yo.shape
    mpar = Xo.shape[1]
    B, U, Sml, ll, aic, hqic, sbic, fpe, det = _var_stats(Yo, Xo)
    Sig = U.T @ U / (T - mpar) if o.get("dfk") else Sml
    XtXi = np.linalg.inv(Xo.T @ Xo)
    V = np.kron(Sig, XtXi)
    b = B.T.ravel()
    small = bool(o.get("small"))
    names, eqs = [], []
    for dv in deps:
        for vn, l in xnames:
            names.append(vn if l == 0 else (f"L.{vn}" if l == 1 else f"L{l}.{vn}"))
            eqs.append(dv)
    # per-equation statistics
    eqstats = []
    for j, dv in enumerate(deps):
        ssr = float(U[:, j] @ U[:, j])
        yc = Yo[:, j] - Yo[:, j].mean()
        r2 = 1 - ssr / float(yc @ yc) if constant else 1 - ssr / float(Yo[:, j] @ Yo[:, j])
        rmse = math.sqrt(ssr / (T - mpar))      # VERIFY: Stata's RMSE divisor
        sl = [i for i, (vn, _l) in enumerate(xnames) if vn != "_cons"]
        bj = B[sl, j]
        Vj = Sig[j, j] * XtXi[np.ix_(sl, sl)]
        stat = float(bj @ np.linalg.solve(Vj, bj)) if sl else SYS
        if small:
            F = stat / len(sl)
            pv = float(st.f.sf(F, len(sl), T - mpar))
            eqstats.append((dv, mpar, rmse, r2, F, pv))
        else:
            eqstats.append((dv, mpar, rmse, r2, stat, float(chi2_tail(stat, len(sl)))))
    rows = []
    se = np.sqrt(np.maximum(np.diag(V), 0))
    i = 0
    for j, dv in enumerate(deps):
        rows.append(CoefRow("eq", dv))
        last = None
        for vn, l in xnames:
            if l > 0:
                if vn != last:
                    if last is not None:
                        rows.append(CoefRow("blank"))
                    rows.append(CoefRow("header", vn))
                    last = vn
                rows.append(CoefRow("coef", f"L{l}.", b[i], se[i], name=f"{dv}:{names[i]}"))
            else:
                if last is not None and last != "__exog":
                    rows.append(CoefRow("blank"))
                last = "__exog"
                rows.append(CoefRow("coef", vn, b[i], se[i], name=f"{dv}:{names[i]}"))
            i += 1
        if j < K - 1:
            rows.append(CoefRow("sep"))
    est = Estimates("var", " ".join(deps), names, b, V, T, "t" if small else "z", (T - mpar) if small else None)
    est.eqnames = eqs
    est.rows = rows
    est.level = float(o["level"]) if o.get("level") else float(s.settings.get("level", "95"))
    est.vcetype = ""
    est.extra = {"N": float(T), "ll": ll, "aic": aic, "hqic": hqic, "sbic": sbic, "fpe": fpe, "det": det,
                 "eqstats": eqstats, "sample": _sample_text(s, t0, t1), "small": small}
    est.var = {"deps": deps, "lags": lags, "exog": exog, "xnames": xnames, "B": B, "Sig": Sig,
               "XtXi": XtXi, "T": T}
    est.display = display_var
    est.predict = _predict_var
    est.predict_opts = {"equation": 2}
    scalars = [("N", float(T)), ("k", float(b.size)), ("k_eq", float(K)), ("k_dv", float(K)),
               ("df_eq", float(mpar)), ("ll", ll), ("aic", aic), ("hqic", hqic), ("sbic", sbic),
               ("fpe", fpe), ("detsig_ml", det), ("tmin", t0), ("tmax", t1), ("mlag", float(max(lags))),
               ("rank", float(b.size))]
    for j, (dv, mp, rmse, r2, stat, pv) in enumerate(eqstats):
        scalars += [(f"rmse_{j + 1}", rmse), (f"r2_{j + 1}", r2), (f"{'F' if small else 'chi2'}_{j + 1}", stat)]
    from ..commands.matrix import Matrix
    mats = [("Sigma", Matrix(Sig, list(deps), list(deps)))]
    macros = [("cmd", "var"), ("cmdline", cmdline("var", args)), ("depvar", " ".join(deps)),
              ("endog", " ".join(deps)), ("exog", " ".join(exog)), ("lags", " ".join(map(str, lags))),
              ("title", "Vector autoregression"), ("predict", "var_p"), ("properties", "b V"),
              ("tvar", require_ts(s.data).tvar)]
    post(s, est, scalars, macros, ok, extra_mats=mats)
    display_var(s, est)


def display_var(s: "Session", est: Estimates, *, header: bool = True, table: bool = True,
                level: float | None = None, **_kw) -> None:
    out = s.output
    x = est.extra
    if header:
        out.write("\nVector autoregression\n\n", "text")
        _right(out, x["sample"], "Number of obs", f"{int(x['N']):,}")
        for lab, v, rl, rv in (("Log likelihood", x["ll"], "AIC", x["aic"]), ("FPE", x["fpe"], "HQIC", x["hqic"]),
                               ("Det(Sigma_ml)", x["det"], "SBIC", x["sbic"])):
            left = f"{lab:<15}= {_g(v, '%9.0g'):>9}"
            _right(out, left, rl, _g(rv, "%9.0g"))
        out.write("\n", "text")
        stat_h = "F" if x["small"] else "chi2"
        out.write(f"{'Equation':<15}{'Parms':>9}{'RMSE':>10}{'R-sq':>9}{stat_h:>10}{'P>' + stat_h:>11}\n", "text")
        out.write("-" * 64 + "\n", "text")
        for dv, mp, rmse, r2, stat, pv in x["eqstats"]:
            out.write(f"{dv[:15]:<15}", "text")
            out.write(f"{mp:>9}{_g(rmse, '%9.0g'):>12}{r2:>9.4f}{_g(stat, '%9.0g'):>11}{pv:>9.4f}\n", "result")
        out.write("-" * 64 + "\n\n", "text")
    if table:
        coef_table(s, "", est.rows, stat=est.stat, df=est.df_r, level=level or est.level)


def _predict_var(s: "Session", est: Estimates, stat: str, mask: np.ndarray):
    raise StataError(198, "predict after var: use the equation() option with xb or residuals; "
                     "not available in OpenDTA yet")   # OpenDTA


def cmd_varsoc(s: "Session", args: str) -> None:
    from ..commands._util import touse
    from ..lang.syntax import parse_standard
    p = parse_standard(args)
    o = match_options(p.options, {"maxlag": 4, "exog": 3, "noconstant": 4, "lutstats": 3}) \
        if p.options.strip() else {}
    deps = s.expand_varlist(p.varlist)
    maxlag = int(o["maxlag"]) if o.get("maxlag") else 4
    exog = s.expand_varlist(o["exog"]) if isinstance(o.get("exog"), str) else []
    constant = not o.get("noconstant")
    # common sample: the one of the largest model
    _Y, _X, ok, _n = _var_fit(s, deps, list(range(1, maxlag + 1)), exog, constant, touse(s, p))
    t0, t1 = _check_gaps(s, ok)
    K = len(deps)
    res = []
    for lag in range(maxlag + 1):
        Y, X, _ok, _nm = _var_fit(s, deps, list(range(1, lag + 1)), exog, constant)
        stats_ = _var_stats(Y[ok], X[ok] if X.shape[1] else np.zeros((int(ok.sum()), 0)))
        res.append(stats_)
    T = int(ok.sum())
    lls = [r_[3] for r_ in res]
    LR = [None] + [2 * (lls[i] - lls[i - 1]) for i in range(1, maxlag + 1)]
    pv = [None] + [float(chi2_tail(LR[i], K * K)) for i in range(1, maxlag + 1)]
    best = {k: int(np.argmin([r_[{"fpe": 7, "aic": 4, "hqic": 5, "sbic": 6}[k]] for r_ in res]))
            for k in ("fpe", "aic", "hqic", "sbic")}
    lr_best = next((i for i in range(maxlag, 0, -1) if pv[i] < 0.05), None)
    out = s.output
    out.write("\n   Selection-order criteria\n", "text")
    out.write(f"   {_sample_text(s, t0, t1):<45}Number of obs      = ", "text")
    out.write(f"{T:>9,}\n", "result")
    out.write("  +" + "-" * 75 + "+\n", "text")
    out.write("  |lag |    LL      LR      df    p      FPE       AIC      HQIC      SBIC    |\n", "text")
    out.write("  |----+" + "-" * 70 + "|\n", "text")

    def star(key, i):
        return "*" if best.get(key) == i else " "
    for i, r_ in enumerate(res):
        line = f"  |{i:>3} | {r_[3]:>8.3f}"
        if i:
            lr = f"{LR[i]:.5g}"
            lr = lr[1:] if lr.startswith("0.") else lr
            line += f"  {lr:>6}{'*' if lr_best == i else ' '}{K * K:>4}  {pv[i]:5.3f}"
        else:
            line += " " * 22
        line += f"  {r_[7]:>7.2g}{star('fpe', i)} {r_[4]:>8.4f}{star('aic', i)} " \
                f"{r_[5]:>8.4f}{star('hqic', i)} {r_[6]:>8.4f}{star('sbic', i)}"   # VERIFY widths
        out.write(f"{line:<78}|\n", "result")
    out.write("  +" + "-" * 75 + "+\n", "text")
    out.write(f"{'Endogenous:':>15}  {' '.join(deps)}\n", "text")
    out.write(f"{'Exogenous:':>15}  {' '.join(exog + (['_cons'] if constant else []))}\n", "text")
    s.r = {"N": float(T), "mlag": float(maxlag)}


def cmd_vargranger(s: "Session", args: str) -> None:
    from .results import current
    est = current(s)
    if est.cmd != "var":
        raise StataError(301, "vargranger only works with estimates from var")   # VERIFY
    V = est.var
    deps, xnames, B, Sig, XtXi = V["deps"], V["xnames"], V["B"], V["Sig"], V["XtXi"]
    out = s.output
    out.write("\n   Granger causality Wald tests\n", "text")
    out.write("  +" + "-" * 66 + "+\n", "text")
    out.write("  |          Equation           Excluded |   chi2     df Prob > chi2 |\n", "text")
    out.write("  |--------------------------------------+---------------------------|\n", "text")
    for j, dv in enumerate(deps):
        others = [d for d in deps if d != dv]
        groups = [(o_, [i for i, (vn, l) in enumerate(xnames) if vn == o_ and l > 0]) for o_ in others]
        groups.append(("ALL", [i for g in groups for i in g[1]]))
        for lab, idx in groups:
            bj = B[idx, j]
            Vj = Sig[j, j] * XtXi[np.ix_(idx, idx)]
            chi2 = float(bj @ np.linalg.solve(Vj, bj))
            pv = float(chi2_tail(chi2, len(idx)))
            out.write(f"  |{dv[:18]:>18} {lab[:18]:>18} |", "text")
            out.write(f" {_g(chi2, '%7.0g'):>7}{len(idx):>6}{pv:>9.3f}    ", "result")
            out.write("|\n", "text")
        out.write("  |--------------------------------------+---------------------------|\n" if j < len(deps) - 1
                  else "  +" + "-" * 66 + "+\n", "text")


# ---------------------------------------------------------------------------
# dfuller
# ---------------------------------------------------------------------------

# Fuller (1976), Table 8.5.2, by sample size 25, 50, 100, 250, 500, infinity
_FULLER_N = [25, 50, 100, 250, 500, 1e12]
_FULLER = {
    "n": {1: [-2.66, -2.62, -2.60, -2.58, -2.58, -2.58], 5: [-1.95] * 6,
          10: [-1.60, -1.61, -1.61, -1.62, -1.62, -1.62]},
    "c": {1: [-3.75, -3.58, -3.51, -3.46, -3.44, -3.43], 5: [-3.00, -2.93, -2.89, -2.88, -2.87, -2.86],
          10: [-2.63, -2.60, -2.58, -2.57, -2.57, -2.57]},
    "ct": {1: [-4.38, -4.15, -4.04, -3.99, -3.98, -3.96], 5: [-3.60, -3.50, -3.45, -3.43, -3.42, -3.41],
           10: [-3.24, -3.18, -3.15, -3.13, -3.13, -3.12]},
}
# MacKinnon (1994), Tables 3 and 4: approximate asymptotic p-values of the tau statistic
# (one variable), p = Phi(polynomial in tau)
_MK = {"n": (-1.04, -19.04, np.inf, [0.6344, 1.2378, 3.2496e-2], [0.4797, 9.3557e-1, -0.6999e-1, 3.3066e-2]),
       "c": (-1.61, -18.83, 2.74, [2.1659, 1.4412, 3.8269e-2], [1.7339, 9.3202e-1, -1.2745e-1, -1.0368e-2]),
       "ct": (-2.89, -16.18, 0.7, [3.2512, 1.6047, 4.9588e-2], [2.5261, 6.1654e-1, -3.7956e-1, -6.0285e-2])}


def mackinnon_p(tau: float, kind: str) -> float:
    star, mn, mx, small, large = _MK[kind]
    if tau > mx:
        return 1.0
    if tau < mn:
        return 0.0
    coef = small if tau <= star else large
    return float(st.norm.cdf(sum(c * tau ** i for i, c in enumerate(coef))))


def _fuller_crit(kind: str, n: int) -> list[float]:
    out = []
    for lev in (1, 5, 10):
        vals = _FULLER[kind][lev]
        if n <= _FULLER_N[0]:
            out.append(vals[0])
            continue
        for i in range(len(_FULLER_N) - 1):
            a, b_ = _FULLER_N[i], _FULLER_N[i + 1]
            if n <= b_:
                if b_ >= 1e12:
                    out.append(vals[i])      # VERIFY: beyond 500 Stata uses the 500 row?
                else:
                    out.append(vals[i] + (n - a) / (b_ - a) * (vals[i + 1] - vals[i]))
                break
    return out


def cmd_dfuller(s: "Session", args: str) -> None:
    from ..commands._util import touse
    from ..lang.syntax import parse_standard
    from .tsops import shift
    p = parse_standard(args)
    o = match_options(p.options, {"noconstant": 4, "trend": 3, "drift": 3, "regress": 3, "lags": 3}) \
        if p.options.strip() else {}
    names = s.expand_varlist(p.varlist)
    if len(names) != 1:
        raise StataError(103 if names else 102, "too many variables specified" if names else "too few variables specified")
    ds = s.data
    yv = ds.get(names[0]).data.astype(float)
    L1 = shift(ds, yv, 1)
    dy = np.where((yv < SYS) & (L1 < SYS), yv - L1, SYS)
    lags = int(o["lags"]) if o.get("lags") else 0
    cols, cn = [L1], [f"L1.{names[0]}"]
    for k in range(1, lags + 1):
        cols.append(shift(ds, dy, k))
        cn.append(f"L{k}D.{names[0]}")
    info = require_ts(ds)
    if o.get("trend"):
        cols.append(ds.get(info.tvar).data.astype(float))
        cn.append("_trend")
    X = np.column_stack(cols)
    ok = touse(s, p) & (dy < SYS) & np.all(X < SYS, axis=1)
    constant = not o.get("noconstant")
    Xo = X[ok]
    if constant:
        Xo = np.column_stack([Xo, np.ones(int(ok.sum()))])
        cn.append("_cons")
    yo = dy[ok]
    n, k = Xo.shape
    b, *_ = np.linalg.lstsq(Xo, yo, rcond=None)
    e = yo - Xo @ b
    s2 = float(e @ e) / (n - k)
    V = s2 * np.linalg.inv(Xo.T @ Xo)
    tau = float(b[0] / math.sqrt(V[0, 0]))
    kind = "ct" if o.get("trend") else ("n" if not constant else "c")
    out = s.output
    title = "Augmented Dickey-Fuller test for unit root" if lags else "Dickey-Fuller test for unit root"
    out.write(f"\n{title:<51}Number of obs   = ", "text")
    out.write(f"{n:>9,}\n\n", "result")
    if o.get("drift"):
        crit = [float(st.t.ppf(a, n - k)) for a in (0.01, 0.05, 0.10)]
        head = "Z(t) has t-distribution"
    else:
        crit = _fuller_crit(kind, n)
        head = "Interpolated Dickey-Fuller"
    nd = 47 - len(head) - 2          # the title is centred over the three critical values
    out.write(f"{'':31}{'-' * ((nd + 1) // 2)} {head} {'-' * (nd // 2)}\n", "text")
    out.write("                  Test         1% Critical       5% Critical      10% Critical\n", "text")
    out.write("               Statistic           Value             Value             Value\n", "text")
    out.write("-" * 78 + "\n", "text")
    out.write(" Z(t)", "text")
    out.write(f"{tau:>19.3f}{crit[0]:>18.3f}{crit[1]:>18.3f}{crit[2]:>18.3f}\n", "result")
    pval = None
    if kind == "n" and not o.get("drift"):
        # noconstant: Stata stops after the statistic, without the closing line and the
        # MacKinnon p-value (compat 1101)
        pass
    elif not o.get("drift"):
        out.write("-" * 78 + "\n", "text")
        pval = mackinnon_p(tau, kind)
        out.write("MacKinnon approximate p-value for Z(t) = ", "text")
        out.write(f"{pval:.4f}\n", "result")
    else:
        out.write("-" * 78 + "\n", "text")
        pval = float(st.t.cdf(tau, n - k))
        out.write("p-value for Z(t) = ", "text")
        out.write(f"{pval:.4f}\n", "result")
    if o.get("regress"):
        se = np.sqrt(np.diag(V))
        # the lags are grouped under the variable name (L1., LD., L2D.)
        rows = [CoefRow("header", names[0])]
        for i, nm in enumerate(cn):
            lab = nm
            if nm.endswith("." + names[0]):
                lab = nm[:-len(names[0])]
                lab = {"L1D.": "LD."}.get(lab, lab)
            rows.append(CoefRow("coef", lab, b[i], se[i]))
        out.write("\n", "text")
        coef_table(s, f"D.{names[0]}", rows, stat="t", df=n - k, level=float(s.settings.get("level", "95")),
                   depvar_left=True)
    s.r = {"N": float(n), "lags": float(lags), "Zt": tau}
    if pval is not None:
        # r(p) only with a constant or trend in the model ([TS] dfuller, Stored results)
        s.r["p"] = pval


# ---------------------------------------------------------------------------
# corrgram
# ---------------------------------------------------------------------------

def _bar(v: float) -> str:
    k = int(math.floor(abs(v) * 8 + 1e-12))
    k = min(k, 8)
    if v >= 0:
        return " " * 8 + "|" + "-" * k
    return " " * (8 - k) + "-" * k + "|"


def cmd_corrgram(s: "Session", args: str) -> None:
    from ..commands._util import touse
    from ..lang.syntax import parse_standard
    p = parse_standard(args)
    o = match_options(p.options, {"lags": 3, "noplot": 3, "yw": 2}) if p.options.strip() else {}
    names = s.expand_varlist(p.varlist)
    if len(names) != 1:
        raise StataError(103, "too many variables specified")
    ds = s.data
    info = require_ts(ds)
    mask = touse(s, p) & (ds.get(names[0]).data < SYS)
    t0, t1 = _check_gaps(s, mask)
    order = np.argsort(ds.get(info.tvar).data[mask], kind="stable")
    y = ds.get(names[0]).data[mask].astype(float)[order]
    n = len(y)
    maxlag = int(o["lags"]) if o.get("lags") else min(int(n / 2) - 2, 40)
    yc = y - y.mean()
    c0 = float(yc @ yc)
    ac = np.array([float(yc[k:] @ yc[:-k]) / c0 for k in range(1, maxlag + 1)])
    pac = []
    for k in range(1, maxlag + 1):
        if o.get("yw"):
            # Yule-Walker (Durbin-Levinson) partial autocorrelations
            r = np.concatenate([[1.0], ac[:k]])
            R = linalg.toeplitz(r[:k])
            phi = np.linalg.solve(R, r[1:k + 1])
            pac.append(float(phi[-1]))
        else:
            Xk = np.column_stack([y[k - j - 1:n - j - 1] for j in range(k)] + [np.ones(n - k)])
            bk, *_ = np.linalg.lstsq(Xk, y[k:], rcond=None)
            pac.append(float(bk[k - 1]))
    Q = n * (n + 2) * np.cumsum(ac ** 2 / (n - np.arange(1, maxlag + 1)))
    out = s.output
    noplot = bool(o.get("noplot"))
    out.write("\n", "text")
    if not noplot:
        out.write(f"{'':42}-1       0       1 -1       0       1\n", "text")
        out.write(" LAG       AC       PAC      Q     Prob>Q  [Autocorrelation]  [Partial Autocor]\n", "text")
    else:
        out.write(" LAG       AC       PAC      Q     Prob>Q\n", "text")
    out.write("-" * (79 if not noplot else 42) + "\n", "text")
    for k in range(maxlag):
        line = f"{k + 1:<5}{ac[k]:>10.4f}{pac[k]:>9.4f}{Q[k]:>9.5g}{float(chi2_tail(Q[k], k + 1)):>8.4f}"
        if not noplot:
            line += f"  {_bar(ac[k]):<19}{_bar(pac[k])}"
        out.write(line.rstrip() + "\n", "result")
    from ..commands.matrix import Matrix
    s.r = {"lags": float(maxlag), "ac": Matrix(ac.reshape(1, -1)), "pac": Matrix(np.array(pac).reshape(1, -1)),
           "Q": Matrix(Q.reshape(1, -1))}
    for k in range(maxlag):
        s.r[f"ac{k + 1}"] = float(ac[k])
        s.r[f"pac{k + 1}"] = float(pac[k])
        s.r[f"q{k + 1}"] = float(Q[k])


# ---------------------------------------------------------------------------
# newey
# ---------------------------------------------------------------------------

def cmd_newey(s: "Session", args: str) -> None:
    from .postest import replay
    from .regress import Sample, cmdline
    if not args.strip() or args.strip().startswith(","):
        replay(s, "newey", args)
        return
    smp = Sample(s, args, opts_spec={"lag": 3, "force": 5}, weights=("aweight",))
    o = smp.o
    if o.get("lag") is None or o.get("lag") is True:
        raise StataError(198, "option lag() required")
    L = int(o["lag"])
    info = require_ts(s.data)
    fvars.drop_collinear(smp.cols, smp.constant)
    est_cols = [col for col in smp.cols if col.values is not None]
    n = smp.n
    X = np.column_stack([col.values for col in est_cols] + ([np.ones(n)] if smp.constant else []))
    y = smp.y.astype(float)
    t = s.data.get(info.tvar).data[smp.mask].astype(float)
    order = np.argsort(t, kind="stable")
    X, y, t = X[order], y[order], t[order]
    if L > 0 and not o.get("force"):
        if np.any(np.abs(np.diff(t) - info.delta) > 1e-9):
            raise StataError(498, f"{info.tvar} is not regularly spaced -- use the force option to override")   # VERIFY
    k = X.shape[1]
    XtXi = np.linalg.inv(X.T @ X)
    b = XtXi @ X.T @ y
    e = y - X @ b
    xe = X * e[:, None]
    Om = xe.T @ xe
    for l in range(1, L + 1):
        w = 1 - l / (L + 1)
        G = xe[l:].T @ xe[:-l]
        Om += w * (G + G.T)
    V = n / (n - k) * XtXi @ Om @ XtXi
    V = (V + V.T) / 2
    se = np.sqrt(np.diag(V))
    sl = list(range(len(est_cols)))
    df_r = n - k
    F = float(b[sl] @ np.linalg.solve(V[np.ix_(sl, sl)], b[sl])) / len(sl) if sl else SYS
    pF = float(st.f.sf(F, len(sl), df_r)) if sl else SYS
    names = [col.name for col in smp.cols] + (["_cons"] if smp.constant else [])
    pos = [i for i, col in enumerate(smp.cols) if col.values is not None] + ([len(smp.cols)] if smp.constant else [])
    kk = len(names)
    bf, Vf = np.zeros(kk), np.zeros((kk, kk))
    bf[pos] = b
    Vf[np.ix_(pos, pos)] = V
    bmap = {col.name: bf[i] for i, col in enumerate(smp.cols)}
    semap = {col.name: float(np.sqrt(max(Vf[i, i], 0))) for i, col in enumerate(smp.cols)}
    if smp.constant:
        bmap["_cons"], semap["_cons"] = bf[-1], float(np.sqrt(max(Vf[-1, -1], 0)))
    est = Estimates("newey", smp.depname, names, bf, Vf, n, "t", float(df_r))
    est.rows = rows_from_columns(smp.cols, bmap, semap, s.data, constant=smp.constant)
    est.cols, est.constant, est.terms, est.depcomp = smp.cols, smp.constant, smp.terms, smp.depcomp
    est.level = smp.level
    est.vcetype = "Newey-West"
    est.extra = {"N": n, "lag": L, "F": F, "pF": pF, "df_m": len(sl), "df_r": df_r}
    est.display = display_newey
    # e() in the order Stata lists it (compat 1101)
    scalars = [("df_m", float(len(sl))), ("df_r", float(df_r)), ("F", F), ("N", float(n)), ("lag", float(L)),
               ("rank", float(k))]
    macros = [("title", "Regression with Newey-West standard errors"), ("cmd", "newey"),
              ("cmdline", cmdline("newey", args)), ("estat_cmd", "newey_estat"), ("predict", "newey_p"),
              ("vcetype", "Newey-West"), ("depvar", smp.depname), ("properties", "b V")]
    post(s, est, scalars, macros, smp.mask)
    display_newey(s, est)


def display_newey(s: "Session", est: Estimates, *, header: bool = True, table: bool = True,
                  level: float | None = None, **_kw) -> None:
    out = s.output
    x = est.extra
    if header:
        out.write("\n", "text")
        _right(out, "Regression with Newey-West standard errors", "Number of obs", f"{int(x['N']):,}")
        _right(out, f"maximum lag: {x['lag']}", f"F({int(x['df_m']):>3},{int(x['df_r']):>10})",
               _g(x["F"], "%10.2f"))
        _right(out, "", "Prob > F", _g(x["pF"], "%10.4f"))
        out.write("\n", "text")
    if table:
        coef_table(s, est.depvar, est.rows, stat="t", df=est.df_r, level=level or est.level,
                   vcetype="Newey-West")
