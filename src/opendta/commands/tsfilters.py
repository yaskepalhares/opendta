"""tsfilter and tssmooth ([TS] tsfilter, [TS] tssmooth).

The computations are in opendta.core.tsstats. tsfilter prints nothing; the cyclical component
(and the trend) are new variables. tssmooth ma prints the smoother applied; the recursive
smoothers print their coefficients and the sum of squared residuals. The layouts the manuals
do not show are marked VERIFY (compat/do/1904_timeseries.do)."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core import tsstats as T
from ..core.errors import StataError
from ..core.formats import format_value
from ..estimation.ml import maximize, print_log
from ..estimation.tsops import require_ts, shift
from ..lang.syntax import match_options, parse_standard
from ._util import eval_vector, touse
from .registry import command
from .tshelp import (PERIODS_PER_YEAR, append_periods, from_nan, kv, new_variable, panels,
                     split_type, to_nan, unit_code, unit_name)

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def _num(s: "Session", text: str) -> float:
    return float(s.eval(text.strip()))


def _series(s: "Session", tok: str) -> np.ndarray:
    """Values of a variable or of an expression with time-series operators (NaN for missing)."""
    if s.data.has(tok):
        return to_nan(s.data.get(tok).data)
    v = eval_vector(s, tok)
    if not isinstance(v, np.ndarray):
        v = np.full(s.data.nobs, float(v))
    return to_nan(v)


def _sources(s: "Session", exp: str) -> list[tuple[str, str]]:
    """[(name used for stubs, token)] of the right side of `newvars = varlist`."""
    out = []
    for tok in exp.split():
        if "*" in tok or "?" in tok:
            out += [(n, n) for n in s.expand_varlist(tok)]
        else:
            out.append((tok.split(".")[-1], tok))
    if not out:
        raise StataError(100, "varlist required")
    return out


def _targets(spec: str, srcs: list[tuple[str, str]], what: str) -> list[str]:
    """Names of the new variables: a list of the same length as the sources or `stub*`."""
    words = spec.split()
    if len(words) == 1 and words[0].endswith("*"):
        return [words[0][:-1] + n for n, _t in srcs]
    if len(words) != len(srcs):
        raise StataError(198, f"number of {what} variables must match the number of variables")
    return words


# ---------------------------------------------------------------------------
# tsfilter
# ---------------------------------------------------------------------------

_FILTERS = ("bk", "bw", "cf", "hp")
_SPEC = {
    "hp": {"smooth": 1, "trend": 1, "gain": 1},
    "bk": {"minperiod": 3, "maxperiod": 4, "smaorder": 2, "stationary": 3, "trend": 1, "gain": 1},
    "cf": {"minperiod": 3, "maxperiod": 4, "smaorder": 2, "stationary": 3, "drift": 1, "trend": 1,
           "gain": 1},
    "bw": {"maxperiod": 4, "order": 1, "trend": 1, "gain": 1},
}
_METHOD = {"hp": "Hodrick-Prescott", "bk": "Baxter-King", "cf": "Christiano-Fitzgerald",
           "bw": "Butterworth"}


def _years(s: "Session", years: float, other: float) -> float:
    c = unit_code(s)
    return round(years * PERIODS_PER_YEAR[c]) if c else other


@command("tsfilter", "tsfilter")
def cmd_tsfilter(s: "Session", args: str) -> None:
    head, _, rest = args.strip().partition(" ")
    method = next((m for m in _FILTERS if m == head), None)
    if method is None:
        raise StataError(198, "invalid syntax")
    info = require_ts(s.data)
    p = parse_standard(rest)
    o = match_options(p.options, _SPEC[method]) if p.options.strip() else {}
    if p.exp is None:
        raise StataError(198, "invalid syntax")
    vtype, lhs = split_type(p.varlist.split())
    srcs = _sources(s, p.exp)
    news = _targets(" ".join(lhs), srcs, "new")
    trends = _targets(o["trend"], srcs, "trend") if isinstance(o.get("trend"), str) else []
    for nm in news + trends:
        if s.data.has(nm):
            raise StataError(110, f"{nm} already defined")
    code = unit_code(s)
    pairs: list[tuple[str, object]] = []
    if method == "hp":
        if o.get("smooth") is not None:
            lam = _num(s, str(o["smooth"]))
            if lam <= 0:
                raise StataError(198, "smooth() must be greater than 0")
        else:
            lam = 1600.0 * (PERIODS_PER_YEAR[code] / 4) ** 4 if code else 1600.0
        pairs.append(("smooth", lam))
    elif method == "bw":
        mx = _num(s, str(o["maxperiod"])) if o.get("maxperiod") is not None else _years(s, 8, 32)
        order = int(_num(s, str(o["order"]))) if o.get("order") is not None else 2
        if mx <= 2:
            raise StataError(198, "maxperiod() must be greater than 2")
        pairs += [("order", float(order)), ("maxperiod", float(mx))]
    else:
        mn = _num(s, str(o["minperiod"])) if o.get("minperiod") is not None else _years(s, 1.5, 6)
        mx = _num(s, str(o["maxperiod"])) if o.get("maxperiod") is not None else _years(s, 8, 32)
        if mn < 2 or mn >= mx:
            raise StataError(198, "minperiod() must be at least 2 and less than maxperiod()")
        q = None
        if o.get("smaorder") is not None:
            q = int(_num(s, str(o["smaorder"])))
        elif method == "bk":
            q = int(_years(s, 3, 12))
        if method == "cf" and q is None:
            pass
        pairs += ([("smaorder", float(q))] if q is not None else []) + \
            [("minperiod", float(mn)), ("maxperiod", float(mx))]
    stationary = bool(o.get("stationary"))
    outcyc = {nm: np.full(s.data.nobs, np.nan) for nm in news}
    outtr = {nm: np.full(s.data.nobs, np.nan) for nm in trends}
    sample = touse(s, p)
    gain_data = None
    for (_n, tok), nm, trn in zip(srcs, news, trends or [None] * len(news)):
        x = _series(s, tok)
        for _pv, rows, t in panels(s, sample & ~np.isnan(x)):
            if not len(rows):
                continue
            if np.any(np.abs(np.diff(t) - info.delta) > 1e-9):
                raise StataError(498, "sample may not contain gaps")      # VERIFY: error text
            y = x[rows]
            if method == "hp":
                c = T.hp_cycle(y, lam)
            elif method == "bk":
                w = T.bk_weights(mn, mx, q, stationary)
                if len(y) <= 2 * q:
                    raise StataError(498, "smaorder() too large for the sample")      # VERIFY
                c = T.bk_cycle(y, w)
            elif method == "cf":
                c = T.cf_cycle(y, mn, mx, stationary=stationary, drift=bool(o.get("drift")),
                               smaorder=q)
            else:
                c = T.bw_cycle(y, mx, order)
            outcyc[nm][rows] = c
            if trn:
                outtr[trn][rows] = y - c
            if gain_data is None and isinstance(o.get("gain"), str):
                N = len(y)
                om = math.pi * np.arange(1, N + 1) / N
                if method == "hp":
                    g = T.hp_gain(lam, om)
                elif method == "bk":
                    g = T.bk_gain(w, om)
                elif method == "bw":
                    g = T.bw_gain(mx, order, om)
                else:
                    raise StataError(198, "option gain() is not available for cf in OpenDTA")
                gain_data = (om, g)
    for nm in news:
        new_variable(s, nm, outcyc[nm], vtype=vtype)
    for nm in trends:
        new_variable(s, nm, outtr[nm], vtype=vtype)
    if gain_data is not None:
        gv = str(o["gain"]).split()
        if len(gv) != 2:
            raise StataError(198, "gain() requires two variable names")
        n = len(gain_data[0])
        for nm, vals in zip(gv, gain_data[::-1]):
            col = np.full(s.data.nobs, np.nan)
            col[:min(n, s.data.nobs)] = vals[:s.data.nobs]
            new_variable(s, nm, col)
    kv(s, pairs + [("varlist", " ".join(t for _n, t in srcs)), ("filterlist", " ".join(news))]
       + ([("trendlist", " ".join(trends))] if trends else [])
       + [("method", _METHOD[method])]
       + ([("symmetric", "symmetric")] if method == "cf" and o.get("smaorder") is not None else [])
       + ([("drift", "drift")] if method == "cf" and o.get("drift") else [])
       + ([("stationary", "stationary")] if method in ("bk", "cf") else [])
       + [("unit", unit_name(s))])


# ---------------------------------------------------------------------------
# tssmooth
# ---------------------------------------------------------------------------

_SMOOTHERS = ("ma", "exponential", "dexponential", "hwinters", "shwinters", "nl")


def _method(word: str) -> str:
    # minimal abbreviations: ma, e, d, h, s, nl
    for name, mn in (("ma", 2), ("exponential", 1), ("dexponential", 1), ("hwinters", 1),
                     ("shwinters", 1), ("nl", 2)):
        if len(word) >= mn and name.startswith(word):
            return name
    raise StataError(198, "invalid syntax")


def _fmt_w(x: float) -> str:
    return f"{x:g}"


def _ma_text(nl: int, nc: int, nf: int, wl: list, w0: float, wf: list, weighted: bool, exp: str) -> str:
    terms = []
    for k in range(nl, 0, -1):
        terms.append((f"{_fmt_w(wl[nl - k])}*" if weighted else "") + f"x(t-{k})")
    if nc or weighted:
        terms.append(f"{_fmt_w(w0) if weighted else 1}*x(t)")
    for k in range(1, nf + 1):
        terms.append((f"{_fmt_w(wf[k - 1])}*" if weighted else "") + f"x(t+{k})")
    tot = (sum(wl) + (w0 if (nc or weighted) else 0) + sum(wf)) if weighted else len(terms)
    den = f"1/{_fmt_w(tot)}"
    return f"({den})*[{' + '.join(terms)}]; x(t)= {exp}"


def _g(x: float, fmt: str = "%9.0g") -> str:
    return format_value(float(x), fmt, pad=False).strip()


def _numlist(s: "Session", text: str) -> list[float]:
    from ..lang.words import parse_numlist
    return [float(v) for v in parse_numlist(text)]


def _ma(s: "Session", vtype, lhs, p, o, exp: str) -> None:
    info = require_ts(s.data)
    if isinstance(o.get("window"), str) == isinstance(o.get("weights"), str):
        raise StataError(198, "option window() or weights() required")
    weighted = isinstance(o.get("weights"), str)
    if weighted:
        txt = o["weights"]
        if "<" not in txt or ">" not in txt:
            raise StataError(198, "weights() requires the weight of the current term in <>")
        a, rest = txt.split("<", 1)
        c, b = rest.split(">", 1)
        wl, wf = _numlist(s, a) if a.strip() else [], _numlist(s, b) if b.strip() else []
        w0 = float(c)
        nl, nf, nc = len(wl), len(wf), 1
    else:
        v = [int(float(x)) for x in str(o["window"]).split()]
        if not 1 <= len(v) <= 3:
            raise StataError(198, "window() invalid")
        nl = v[0]
        nc = v[1] if len(v) > 1 else 0
        nf = v[2] if len(v) > 2 else 0
        if nc not in (0, 1):
            raise StataError(198, "window(): the second number must be 0 or 1")
        wl, w0, wf = [1.0] * nl, float(nc), [1.0] * nf
    x = _series(s, exp)
    sample = touse(s, p)
    x = np.where(sample, x, np.nan)
    ds = s.data
    tvals = ds.get(info.tvar).data.astype(float)
    out = np.full(ds.nobs, np.nan)
    for _pv, rows, t in panels(s, sample):
        if not len(rows):
            continue
        lo, hi = t.min(), t.max()
        num = np.zeros(ds.nobs)
        den = np.zeros(ds.nobs)
        cnt = np.zeros(ds.nobs)
        # offsets: -k is a lag of k periods, +k a lead
        terms = [(-k, wl[nl - k]) for k in range(nl, 0, -1)] + [(0, w0)] + \
                [(k, wf[k - 1]) for k in range(1, nf + 1)]
        for off, wt in terms:
            if wt == 0:
                continue
            sh = shift(ds, np.where(np.isnan(x), SYS, x), -off)
            ok = sh < SYS
            num += np.where(ok, wt * np.where(ok, sh, 0.0), 0.0)
            den += np.where(ok, wt, 0.0)
            cnt += ok
        edge = (tvals - nl * info.delta >= lo - 1e-9) & (tvals + nf * info.delta <= hi + 1e-9)
        # VERIFY: Stata 14 keeps all observations; terms outside the sample count as zero and
        # the denominator stays the full sum of the weights
        full = sum(wt for _o, wt in terms)
        res = np.where(~edge & (full != 0), num / (full if full != 0 else 1.0), num / np.where(den == 0, 1, den))
        res = np.where(den == 0, np.nan, res)
        out[rows] = res[rows]
    new_variable(s, lhs[0], out, vtype=vtype, replace=bool(o.get("replace")))
    s.output.write("The smoother applied was\n", "text")
    s.output.write("     " + _ma_text(nl, nc, nf, wl, w0, wf, weighted, exp) + "\n", "text")
    N = int((~np.isnan(out)).sum())
    pairs = [("N", float(N)), ("w0", w0)]
    pairs += [(f"wlead{k}", wf[k - 1]) for k in range(1, nf + 1)]
    pairs += [(f"wlag{k}", wl[nl - k]) for k in range(1, nl + 1)]
    pairs += [("method", "ma"), ("exp", exp), ("timevar", info.tvar)]
    if info.panel:
        pairs.append(("panelvar", info.panel))
    kv(s, pairs)


def _grid(info, rows: np.ndarray, t: np.ndarray, x: np.ndarray):
    """Dense series on the time grid of one panel, trimmed to the first and last nonmissing."""
    ok = ~np.isnan(x[rows])
    if not ok.any():
        return None
    first, last = t[ok].min(), t[ok].max()
    n = int(round((last - first) / info.delta)) + 1
    g = np.full(n, np.nan)
    sel = ok
    g[np.round((t[sel] - first) / info.delta).astype(int)] = x[rows][sel]
    return g, first


def _recursive(s: "Session", method: str, vtype, lhs, p, o, exp: str) -> None:
    info = require_ts(s.data)
    ds = s.data
    x = _series(s, exp)
    sample = touse(s, p)
    x = np.where(sample, x, np.nan)
    forecast = int(_num(s, str(o["forecast"]))) if o.get("forecast") is not None else 0
    if not 0 <= forecast <= 500:
        raise StataError(198, "forecast() must be between 0 and 500")
    if o.get("samp0") is not None and o.get("s0") is not None:
        raise StataError(198, "options samp0() and s0() may not be combined")
    samp0 = int(_num(s, str(o["samp0"]))) if o.get("samp0") is not None else None
    parms = [float(v) for v in str(o["parms"]).split()] if isinstance(o.get("parms"), str) else None
    s0 = [float(v) for v in str(o["s0"]).split()] if isinstance(o.get("s0"), str) else None
    npar = {"exponential": 1, "dexponential": 1, "hwinters": 2, "shwinters": 3}[method]
    if parms is not None and len(parms) != npar:
        raise StataError(198, f"parms() requires {npar} number{'s' if npar > 1 else ''}")
    if parms is not None and any(not 0 <= v <= 1 for v in parms):
        raise StataError(198, "parms(): parameters must be between 0 and 1")
    out_vals = np.full(ds.nobs, np.nan)
    extra_entries: list[tuple[float | None, float]] = []
    extra_vals: list[float] = []
    name = lhs[0]
    results = []
    out = s.output
    for pv, rows, t in panels(s, sample):
        g = _grid(info, rows, t, x)
        if g is None:
            continue
        series, first = g
        if method == "exponential":
            if s0 is not None and len(s0) != 1:
                raise StataError(198, "s0() requires one number")
            res = T.exp_smooth(series, alpha=parms[0] if parms else None,
                               s0=s0[0] if s0 else None, samp0=samp0, forecast=forecast)
        elif method == "dexponential":
            if s0 is not None:
                raise StataError(198, "option s0() is not available for dexponential in OpenDTA")
            res = T.dexp_smooth(series, alpha=parms[0] if parms else None, samp0=samp0,
                                forecast=forecast)
        elif method == "hwinters":
            def run(fun, b0):
                return maximize(_numfun(fun), b0)
            if parms is None:
                out.write("computing optimal weights\n\n", "text")
            res = T.hw_smooth(series, parms=tuple(parms) if parms else None,
                              s0=tuple(s0) if s0 else None, samp0=samp0, forecast=forecast,
                              diff=bool(o.get("diff")),
                              start=tuple(float(v) for v in str(o["from"]).split())
                              if isinstance(o.get("from"), str) else (0.5, 0.5), maximize=run)
        else:
            L = int(_num(s, str(o["period"]))) if o.get("period") is not None else \
                int(PERIODS_PER_YEAR.get(unit_code(s), 0))
            if L < 2:
                raise StataError(198, "option period() required")
            sn0 = None
            if isinstance(o.get("sn0_0"), str):
                sv = _series(s, str(o["sn0_0"]))
                sn0 = sv[rows][:L]
            def run(fun, b0):
                return maximize(_numfun(fun), b0)
            if parms is None:
                out.write("computing optimal weights\n\n", "text")
            res = T.shw_smooth(series, L, additive=bool(o.get("additive")),
                               parms=tuple(parms) if parms else None,
                               s0=tuple(s0) if s0 else None, sn0=sn0, samp0=samp0,
                               forecast=forecast, normalize=bool(o.get("normalize")),
                               altstarts=bool(o.get("altstarts")),
                               start=tuple(float(v) for v in str(o["from"]).split())
                               if isinstance(o.get("from"), str) else (0.5, 0.5, 0.5), maximize=run)
        results.append((pv, res))
        fit = res["fit"]
        n = len(series)
        # periods of the grid that exist in the data get their value; the periods after the
        # last one (forecast) are added as observations
        tt = first + info.delta * np.arange(len(fit))
        have = {float(a): int(r) for a, r in zip(t, rows)}
        for a, v in zip(tt, fit):
            if float(a) in have:
                out_vals[have[float(a)]] = v
            elif a > t.max() or forecast > 0 or True:
                extra_entries.append((pv, float(a)))
                extra_vals.append(v)
        _report(s, method, res, parms is None, pv if info.panel else None, o)
    flag = append_periods(s, extra_entries) if extra_entries else None
    if extra_entries:
        # the new observations were appended and the data re-sorted: rebuild the column by
        # matching on (panel, time)
        out_vals = _match_back(s, info, out_vals, extra_entries, extra_vals, flag)
    new_variable(s, name, out_vals, vtype=vtype, replace=bool(o.get("replace")))
    if results:
        pv, res = results[-1]
        pairs = [("N", float(res["N"])), ("alpha", res["alpha"])]
        if method in ("hwinters", "shwinters"):
            pairs.append(("beta", res["beta"]))
        if method == "shwinters":
            pairs.append(("gamma", res["gamma"]))
        if method in ("hwinters", "shwinters"):
            if parms is None:
                pairs.append(("prss", res["prss"]))
        pairs += [("rss", res["rss"]), ("rmse", res["rmse"]), ("N_pre", float(res["N_pre"]))]
        if method != "exponential":
            pairs += [("s2_0", res["s2_0"])]
        pairs += [("s1_0", res["s1_0"])]
        if method != "exponential":
            pairs += [("linear", res["linear"]), ("constant", res["constant"])]
        if method == "shwinters":
            pairs.append(("period", float(res["period"])))
        pairs += [("method", method), ("exp", exp), ("timevar", info.tvar)]
        if info.panel:
            pairs.append(("panelvar", info.panel))
        kv(s, pairs)


def _numfun(fun):
    """Adds numerical gradient and Hessian to a function that returns only the value."""
    from ..estimation.mixed import _numderiv

    def f(b, todo):
        if todo == 0:
            return fun(b, 0)[0], None, None
        return _numderiv(lambda z: fun(z, 0)[0], b)
    return f


def _match_back(s, info, out_vals, entries, vals, flag):
    """Values for the appended observations (flag marks them after the sort); the old
    observations keep theirs, which moved with the sort: recompute by (panel, time)."""
    ds = s.data
    t = ds.get(info.tvar).data.astype(float)
    pv = ds.get(info.panel).data.astype(float) if info.panel else np.zeros(ds.nobs)
    new = np.full(ds.nobs, np.nan)
    key_new = {((e[0] if info.panel else 0.0), float(e[1])): v for e, v in zip(entries, vals)}
    n_old = len(out_vals)
    # the old values: they were indexed by the old row order; rows kept their relative order
    # among themselves (stable sort), so the unflagged rows, in order, are the old rows
    new[~flag] = out_vals
    for i in np.flatnonzero(flag):
        new[i] = key_new.get((float(pv[i]) if info.panel else 0.0, float(t[i])), np.nan)
    del n_old
    return new


def _report(s: "Session", method: str, res: dict, optimal: bool, panel, o: dict) -> None:
    out = s.output
    if panel is not None:
        out.write(f"\nPanel: {_g(panel)}\n", "text")      # VERIFY: panel header
    if method in ("exponential", "dexponential"):
        word = "exponential" if method == "exponential" else "double-exponential"
        # VERIFY: Stata 14 layout (column of "=" and the width of the values)
        if optimal:
            out.write(f"\ncomputing optimal {word}  coefficient (0,1)\n\n", "text")
            lab0, w, val = f"optimal {word} coefficient", 14, lambda x: _g(x, "%10.0g")
        else:
            out.write("\n", "text")
            lab0, w, val = f"{word} coefficient ", 13, lambda x: _g(x, "%12.5g")
        pad = max(len(lab0), 24) + 1
        out.write(f"{lab0:<{pad}}={res['alpha']:>{w}.4f}\n", "text")
        out.write(f"{'sum-of-squared residuals':<{pad}}={val(res['rss']):>{w}}\n", "text")
        out.write(f"{'root mean squared error':<{pad}}={val(res['rmse']):>{w}}\n", "text")
        return
    if res.get("opt") is not None:
        print_log(out, res["opt"], "penalized RSS")
        out.write("\nOptimal weights:\n", "text")
    else:
        out.write("\nSpecified weights:\n", "text")
    out.write(f"{'alpha':>24} = {res['alpha']:.4f}\n", "text")
    out.write(f"{'beta':>24} = {res['beta']:.4f}\n", "text")
    if method == "shwinters":
        out.write(f"{'gamma':>24} = {res['gamma']:.4f}\n", "text")
    if res.get("opt") is not None:
        out.write(f"{'penalized sum-of-squared residuals':>24} = {_g(res['prss'], '%9.0g'):<9}\n", "text")
    out.write(f"{'sum-of-squared residuals':>24} = {_g(res['rss'], '%9.0g'):<9}\n", "text")
    out.write(f"{'root mean squared error':>24} = {_g(res['rmse'], '%9.0g'):<9}\n", "text")


def _nl(s: "Session", vtype, lhs, p, o, exp: str) -> None:
    info = require_ts(s.data)
    if not isinstance(o.get("smoother"), str):
        raise StataError(198, "option smoother() required")
    from ..core import basestats as B
    x = _series(s, exp)
    sample = touse(s, p)
    x = np.where(sample, x, np.nan)
    out = np.full(s.data.nobs, np.nan)
    for _pv, rows, t in panels(s, sample & ~np.isnan(x)):
        if len(rows):
            spec, _, tw = str(o["smoother"]).replace(" ", "").partition(",")
            try:
                out[rows] = B.tukey_smooth(x[rows], spec, tw.startswith("twice"))
            except ValueError as exc:
                raise StataError(198, f"invalid smoother: {exc}")
    new_variable(s, lhs[0], out, vtype=vtype, replace=bool(o.get("replace")))
    kv(s, [("N", float((~np.isnan(out)).sum())), ("method", "nl"), ("exp", exp),
           ("timevar", info.tvar)] + ([("panelvar", info.panel)] if info.panel else []))


@command("tssmooth", "tssmooth")
def cmd_tssmooth(s: "Session", args: str) -> None:
    head, _, rest = args.strip().partition(" ")
    method = _method(head)
    require_ts(s.data)
    p = parse_standard(rest)
    if p.exp is None:
        raise StataError(198, "invalid syntax")
    spec = {"replace": 3, "forecast": 3}
    if method == "ma":
        spec.update({"window": 1, "weights": 1})
    elif method == "nl":
        spec.update({"smoother": 1})
    else:
        spec.update({"parms": 2, "samp0": 3, "s0": 2, "iterate": 4, "from": 4, "tolerance": 3,
                     "ltolerance": 4, "nrtolerance": 4, "nolog": 4, "log": 3})
        if method == "hwinters":
            spec["diff"] = 4
        if method == "shwinters":
            spec.update({"period": 3, "additive": 3, "sn0_0": 6, "sn0_v": 6, "snt_v": 5,
                         "normalize": 4, "altstarts": 3})
    o = match_options(p.options, spec) if p.options.strip() else {}
    vtype, lhs = split_type(p.varlist.split())
    if len(lhs) != 1:
        raise StataError(198, "invalid syntax")
    if s.data.has(lhs[0]) and not o.get("replace"):
        raise StataError(110, f"{lhs[0]} already defined")
    exp = p.exp.strip()
    if method == "ma":
        _ma(s, vtype, lhs, p, o, exp)
    elif method == "nl":
        _nl(s, vtype, lhs, p, o, exp)
    else:
        for bad in ("sn0_v", "snt_v"):
            if o.get(bad):
                raise StataError(198, f"option {bad}() is not available yet in OpenDTA")
        _recursive(s, method, vtype, lhs, p, o, exp)
