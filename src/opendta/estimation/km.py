"""Kaplan-Meier and Nelson-Aalen estimates: sts list, sts generate,
sts test and sts graph ([ST] sts, [ST] sts list, [ST] sts generate,
[ST] sts test, [ST] sts graph).

Formulas from [ST] sts, Methods and formulas: product-limit estimate with
delayed entry, Greenwood standard error and the confidence interval on
the ln(-ln S) scale; the log-rank and Wilcoxon (Breslow) tests use the
hypergeometric variance. Layout details carry the VERIFY mark
(compat/do/0901_st.do).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..lang.syntax import find_top, match_options
from ..lazy import lazy
from .st import STData, _g, _level_label, _qualifiers, header, require_st, st_data

# imported on first use (start-up time; see opendta/lazy.py)
st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


@dataclass
class KM:
    times: np.ndarray       # times where something happens (exit or entry)
    n: np.ndarray           # number at risk just before the time
    d: np.ndarray           # failures at the time
    c: np.ndarray           # censored at the time
    e: np.ndarray           # entries at the time (they join the risk set after it)
    S: np.ndarray           # survivor function just after the time
    gw: np.ndarray          # Greenwood sum of d/(n(n-d))
    H: np.ndarray           # Nelson-Aalen cumulative hazard
    hv: np.ndarray          # sum of d/n^2 (variance of H)


def km(t: np.ndarray, t0: np.ndarray, d: np.ndarray, w: np.ndarray) -> KM:
    ent = t0[t0 > 0]
    times = np.unique(np.concatenate([t, ent]))
    ts, t0s = np.sort(t), np.sort(t0)
    # at risk at time x: t0 < x <= t  ->  #(t >= x) - #(t0 >= x)
    o_t, o_t0 = np.argsort(t, kind="stable"), np.argsort(t0, kind="stable")
    cw_t = np.concatenate([np.cumsum(w[o_t][::-1])[::-1], [0.0]])
    cw_t0 = np.concatenate([np.cumsum(w[o_t0][::-1])[::-1], [0.0]])
    n = cw_t[np.searchsorted(ts, times, "left")] - cw_t0[np.searchsorted(t0s, times, "left")]
    idx = np.searchsorted(times, t)
    dd = np.bincount(idx, weights=w * d, minlength=len(times))
    cc = np.bincount(idx, weights=w * (1 - d), minlength=len(times))
    ee = np.zeros(len(times))
    if ent.size:
        np.add.at(ee, np.searchsorted(times, t0[t0 > 0]), w[t0 > 0])
    with np.errstate(divide="ignore", invalid="ignore"):
        q = np.where(n > 0, 1 - dd / np.where(n > 0, n, 1), 1.0)
        S = np.cumprod(q)
        term = np.where((dd > 0) & (n > dd), dd / (n * (n - dd)), 0.0)
        gw = np.cumsum(term)
        gw = np.where(np.cumsum((dd > 0) & (n <= dd)) > 0, np.inf, gw)
        H = np.cumsum(np.where(n > 0, dd / np.where(n > 0, n, 1), 0.0))
        hv = np.cumsum(np.where(n > 0, dd / np.where(n > 0, n, 1) ** 2, 0.0))
    return KM(times, n, dd, cc, ee, S, gw, H, hv)


def at_times(k: KM, x: np.ndarray, what: str) -> np.ndarray:
    """Step function evaluated at x (value just after x)."""
    j = np.searchsorted(k.times, x, "right") - 1
    arr = getattr(k, what)
    base = {"S": 1.0, "gw": 0.0, "H": 0.0, "hv": 0.0}[what]
    return np.where(j >= 0, arr[np.maximum(j, 0)], base)


def km_ci(S: np.ndarray, gw: np.ndarray, level: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Greenwood standard error and interval on the ln(-ln S) scale ([ST] sts)."""
    z = st.norm.ppf(1 - (1 - level / 100) / 2)
    with np.errstate(divide="ignore", invalid="ignore"):
        se = S * np.sqrt(gw)
        lns = np.log(S)
        sig = np.sqrt(gw) / np.abs(lns)
        lo = S ** np.exp(z * sig)
        hi = S ** np.exp(-z * sig)
    ok = (S > 0) & (S < 1) & np.isfinite(gw)
    # before the first failure (S = 1) the standard error is shown as missing (compat 0901)
    return (np.where(np.isfinite(se) & (S > 0) & (S < 1), se, SYS), np.where(ok, lo, SYS), np.where(ok, hi, SYS))


def quantile_time(k: KM, p: float) -> float:
    """Smallest failure time with S(t) <= 1 - p (VERIFY: <= or <)."""
    target = 1 - p
    hit = np.nonzero((k.d > 0) & (k.S <= target + 1e-12))[0]
    return float(k.times[hit[0]]) if hit.size else SYS


# ---------------------------------------------------------------------------
# sts dispatcher
# ---------------------------------------------------------------------------

def cmd_sts(s: "Session", args: str) -> None:
    text = args.strip()
    word = text.split(None, 1)[0] if text and not text.startswith(",") else ""
    w = word.split(",", 1)[0]
    rest = text[len(w):] if w else text
    if w and len(w) >= 1 and "list".startswith(w):
        sts_list(s, rest)
    elif w and len(w) >= 1 and "generate".startswith(w):
        sts_generate(s, rest)
    elif w and len(w) >= 1 and "test".startswith(w):
        sts_test(s, rest)
    elif w and len(w) >= 1 and "graph".startswith(w):
        sts_graph(s, rest)
    elif not w:
        sts_graph(s, rest)
    else:
        raise StataError(198, f"{word} invalid subcommand")   # VERIFY


def _by_groups(s: "Session", mask: np.ndarray, byvars: list[str]):
    """[(label, STData)] for each combination of by() (a single group without by())."""
    if not byvars:
        return [("", st_data(s, mask))]
    from ..core.varlist import resolve_name
    ds = s.data
    names = [resolve_name(ds, v) for v in byvars]
    base = st_data(s, mask)
    cols = [ds.get(v).data for v in names]
    valid = base.mask.copy()
    for v, c in zip(names, cols):
        if not ds.get(v).is_string:
            valid &= c < SYS
    keys = sorted({tuple(c[i] for c in cols) for i in np.nonzero(valid)[0]})
    out = []
    for key in keys:
        m = valid.copy()
        for c, kv in zip(cols, key):
            m &= c == kv
        lab = ", ".join(f"{v} = {_level_label(s, v, kv)}" for v, kv in zip(names, key))
        out.append((lab, st_data(s, m), names, key))
    return [(g[0], g[1]) for g in out]


# ---------------------------------------------------------------------------
# sts list
# ---------------------------------------------------------------------------

def sts_list(s: "Session", args: str) -> None:
    _vl, mask, o = _qualifiers(s, args, {"by": 2, "failure": 1, "at": 2, "level": 1, "noshow": 6,
                                         "enter": 5, "compare": 3, "survival": 3})
    info = require_st(s)
    if not o.get("noshow"):
        header(s, info)
    level = float(o["level"]) if o.get("level") else float(s.settings.get("level", "95"))
    byv = str(o["by"]).split() if isinstance(o.get("by"), str) else []
    fail = bool(o.get("failure"))
    out = s.output
    lv = f"{level:g}"
    func = "Failure" if fail else "Survivor"
    at = None
    if o.get("at"):
        from ..lang.words import parse_numlist
        at = np.array(parse_numlist(str(o["at"])), dtype=float)
    if at is None:
        out.write(f"{'':11}Beg.          Net{'':12}{func:>8}      Std.\n", "text")
        out.write(f"  Time    Total   Fail   Lost{'':11}Function     Error     [{lv}% Conf. Int.]\n", "text")
    else:
        out.write(f"{'':14}Beg.{'':22}{func:<8}      Std.\n", "text")
        out.write(f"    Time     Total     Fail{'':13}Function     Error     [{lv}% Conf. Int.]\n", "text")
    out.write("-" * 79 + "\n", "text")
    groups = _by_groups(s, mask, byv)
    for lab, sd in groups:
        if lab:
            out.write(f"{lab.replace(' = ', '=')} \n", "text")    # VERIFY
        k = km(sd.t, sd.t0, sd.d, sd.w)
        se, lo, hi = km_ci(k.S, k.gw, level)
        if at is None:
            for i in range(len(k.times)):
                if sd.ids is not None and k.d[i] == 0 and k.c[i] == k.e[i]:
                    continue        # a subject's record boundary: nothing happens (VERIFY)
                S, a, b = k.S[i], lo[i], hi[i]
                if fail:
                    S, a, b = 1 - S, (1 - b if b < SYS else SYS), (1 - a if a < SYS else SYS)
                _list_row(out, k.times[i], k.n[i], k.d[i], k.c[i] - k.e[i], S, se[i], a, b)
        else:
            prev = -np.inf
            for x in at:
                sel = (k.times > prev) & (k.times <= x)
                # Beg. Total is the number at risk at the last listed time on or before x
                # (compat 0901)
                j = np.searchsorted(k.times, x, "right") - 1
                nb = float(k.n[j]) if j >= 0 else at_risk(sd, x)
                S = float(k.S[j]) if j >= 0 else 1.0
                sev = float(se[j]) if j >= 0 else SYS
                a, b = (float(lo[j]), float(hi[j])) if j >= 0 else (SYS, SYS)
                if fail:
                    S, a, b = 1 - S, (1 - b if b < SYS else SYS), (1 - a if a < SYS else SYS)
                _at_row(out, x, nb, float(k.d[sel].sum()), S, sev, a, b)
                prev = x
    out.write("-" * 79 + "\n", "text")
    if at is not None:
        out.write(f"Note: {func} function is calculated over full data and evaluated at indicated times; "
                  "it is not calculated from aggregates shown at left.\n", "text")


def _fmt4(v: float) -> str:
    return f"{v:.4f}" if v < SYS else "."


def _at_row(out, t, n, d, S, se, lo, hi) -> None:
    out.write(f"{_g(t, '%8.0g'):>8}{_g(n, '%9.0g'):>10}{_g(d, '%8.0g'):>9}", "result")
    out.write(f"{_fmt4(S):>20}{_fmt4(se):>10}{_fmt4(lo):>11}{_fmt4(hi):>10}\n", "result")


def at_risk(sd: STData, x: float) -> float:
    return float(np.sum(sd.w[(sd.t0 < x) & (sd.t >= x)]))


def _list_row(out, t, n, d, lost, S, se, lo, hi) -> None:
    out.write(f"{_g(t, '%6.0g'):>6} {_g(n, '%8.0g'):>8} {_g(d, '%6.0g'):>6} {_g(lost, '%6.0g'):>6}", "result")
    out.write(f"{_fmt4(S):>19} {_fmt4(se):>9} {_fmt4(lo):>10} {_fmt4(hi):>9}\n", "result")


# ---------------------------------------------------------------------------
# sts generate
# ---------------------------------------------------------------------------

_GEN_FUNCS = ("s", "se(s)", "lb(s)", "ub(s)", "na", "se(na)", "lb(na)", "ub(na)", "n", "d", "h")


def sts_generate(s: "Session", args: str) -> None:
    from ..core.dataset import Variable
    text = args.strip()
    c = find_top(text, ",")
    body, opts = (text, "") if c == -1 else (text[:c], text[c + 1:])
    o = match_options(opts, {"by": 2, "level": 1, "adjustfor": 3}) if opts.strip() else {}
    level = float(o["level"]) if o.get("level") else float(s.settings.get("level", "95"))
    pairs = []
    toks = body.replace("=", " = ").split()
    i = 0
    while i < len(toks):
        if i + 2 >= len(toks) or toks[i + 1] != "=":
            raise StataError(198, "invalid syntax")
        name, fn = toks[i], toks[i + 2].replace(" ", "").lower()
        if fn not in _GEN_FUNCS:
            raise StataError(198, f"{toks[i + 2]} invalid")   # VERIFY
        pairs.append((name, fn))
        i += 3
    if not pairs:
        raise StataError(198, "invalid syntax")
    ds = s.data
    for name, _ in pairs:
        if ds.has(name):
            raise StataError(110, f"{name} already defined")
    require_st(s)
    byv = str(o["by"]).split() if isinstance(o.get("by"), str) else []
    vals = {name: np.full(ds.nobs, SYS) for name, _ in pairs}
    z = st.norm.ppf(1 - (1 - level / 100) / 2)
    for _lab, sd in _by_groups(s, np.ones(ds.nobs, dtype=bool), byv):
        k = km(sd.t, sd.t0, sd.d, sd.w)
        x = sd.t
        S = at_times(k, x, "S")
        gw = at_times(k, x, "gw")
        # the Greenwood sum goes through a float temporary, as Stata's ado-file appears to do:
        # with a double sum lb(s) differs in the last displayed digit (VERIFY, compat 0901)
        gw = np.where(np.isfinite(gw), gw.astype(np.float32).astype(np.float64), gw)
        H = at_times(k, x, "H")
        hv = at_times(k, x, "hv")
        se, lo, hi = km_ci(S, gw, level)
        j = np.searchsorted(k.times, x)
        for name, fn in pairs:
            if fn == "s":
                v = S
            elif fn == "se(s)":
                v = se
            elif fn == "lb(s)":
                v = lo
            elif fn == "ub(s)":
                v = hi
            elif fn == "na":
                v = H
            elif fn == "se(na)":
                v = np.sqrt(hv)
            elif fn in ("lb(na)", "ub(na)"):
                with np.errstate(divide="ignore", invalid="ignore"):
                    f = np.exp((-1 if fn == "lb(na)" else 1) * z * np.sqrt(hv) / H)
                v = np.where(H > 0, H * f, SYS)   # VERIFY: transformation used by Stata
            elif fn == "n":
                v = k.n[j]
            elif fn == "d":
                v = k.d[j]
            else:   # h: hazard contribution d/n at the time
                v = np.where(k.n[j] > 0, k.d[j] / np.where(k.n[j] > 0, k.n[j], 1), SYS)
            arr = vals[name]
            arr[sd.mask] = np.where(np.isfinite(v), v, SYS)
    for name, _ in pairs:
        v = vals[name]
        ok = v < SYS
        out_v = np.full(ds.nobs, SYS)
        out_v[ok] = v[ok].astype(np.float32)
        ds.add(Variable(name, "float", out_v))
    s.notify_state()


# ---------------------------------------------------------------------------
# sts test
# ---------------------------------------------------------------------------

def sts_test(s: "Session", args: str) -> None:
    from ..core.stats import chi2_tail
    from ..core.varlist import resolve_name
    vl, mask, o = _qualifiers(s, args, {"logrank": 1, "wilcoxon": 1, "tware": 2, "peto": 2,
                                        "strata": 3, "detail": 3, "noshow": 6, "trend": 3})
    info = require_st(s)
    names = [resolve_name(s.data, v) for v in s.expand_varlist(vl)] if vl else []
    if not names:
        raise StataError(100, "varlist required")
    test = "wilcoxon" if o.get("wilcoxon") else "tware" if o.get("tware") else \
        "peto" if o.get("peto") else "logrank"
    if not o.get("noshow"):
        header(s, info)
    ds = s.data
    sd = st_data(s, mask)
    cols = [ds.get(v).data[sd.mask] for v in names]
    ok = np.ones(len(sd.t), dtype=bool)
    for v, c in zip(names, cols):
        if not ds.get(v).is_string:
            ok &= c < SYS
    keys = [tuple(c[i] for c in cols) for i in range(len(sd.t))]
    levels = sorted({kk for kk, good in zip(keys, ok) if good})
    G = len(levels)
    grp = np.array([levels.index(kk) if good else -1 for kk, good in zip(keys, ok)])
    strat = np.zeros(len(sd.t), dtype=int)
    svars: list[str] = []
    if isinstance(o.get("strata"), str):
        svars = [resolve_name(ds, v) for v in s.expand_varlist(o["strata"])]
        sk = [tuple(ds.get(v).data[sd.mask][i] for v in svars) for i in range(len(sd.t))]
        su = sorted(set(sk))
        strat = np.array([su.index(x) for x in sk])
    O = np.zeros(G)
    E = np.zeros(G)
    V = np.zeros((G, G))
    R = np.zeros(G)
    for h in np.unique(strat):
        m = (strat == h) & (grp >= 0)
        t, t0, d, w, g = sd.t[m], sd.t0[m], sd.d[m], sd.w[m], grp[m]
        ft = np.unique(t[d > 0])
        km_all = km(t, t0, d, w) if test in ("peto",) else None
        for x in ft:
            risk = (t0 < x) & (t >= x)
            nt = float(w[risk].sum())
            dt = float(w[(t == x) & (d > 0)].sum())
            ng = np.bincount(g[risk], weights=w[risk], minlength=G)
            dg = np.bincount(g[(t == x) & (d > 0)], weights=w[(t == x) & (d > 0)], minlength=G)
            if test == "logrank":
                wt = 1.0
            elif test == "wilcoxon":
                wt = nt
            elif test == "tware":
                wt = np.sqrt(nt)
            else:
                j = np.searchsorted(km_all.times, x) - 1
                wt = float(km_all.S[j]) if j >= 0 else 1.0
                wt = wt * nt / (nt + 1)   # VERIFY: Peto-Peto-Prentice weight
            eg = dt * ng / nt
            O += wt * dg if test != "logrank" else dg
            E += wt * eg if test != "logrank" else eg
            R += wt * (dg - eg)
            if nt > 1:
                fac = dt * (nt - dt) / (nt - 1) / nt
                V += wt * wt * fac * (np.diag(ng) - np.outer(ng, ng) / nt)
    if test == "logrank":
        O_disp, E_disp = O, E
    else:
        # the table shows the unweighted counts; the statistic uses the weighted sums
        O_disp, E_disp = _unweighted_oe(sd, grp, strat, G)
    diff = R if test != "logrank" else O - E
    k = G - 1
    try:
        chi2 = float(diff[:k] @ np.linalg.solve(V[:k, :k], diff[:k])) if k > 0 else SYS
    except np.linalg.LinAlgError:
        chi2 = float(diff[:k] @ np.linalg.pinv(V[:k, :k]) @ diff[:k])
    p = float(chi2_tail(chi2, k)) if k > 0 else SYS
    out = s.output
    title = {"logrank": "Log-rank test for equality of survivor functions",
             "wilcoxon": "Wilcoxon (Breslow) test for equality of survivor functions",
             "tware": "Tarone-Ware test for equality of survivor functions",
             "peto": "Peto-Peto test for equality of survivor functions"}[test]
    nstrata = len(np.unique(strat))
    if nstrata > 1 or isinstance(o.get("strata"), str):
        title = "Stratified " + title[0].lower() + title[1:]
    out.write(f"\n{title}\n", "text")
    if isinstance(o.get("strata"), str):
        nofail = sum(1 for h in np.unique(strat) if not np.any(sd.d[(strat == h) & (grp >= 0)] > 0))
        if nofail:
            out.write(f"Note: {nofail} strat{'um' if nofail == 1 else 'a'} omitted because of no failures\n",
                      "text")   # VERIFY: singular form
    out.write("\n", "text")
    stratified = isinstance(o.get("strata"), str)
    labs = [", ".join(_level_label(s, v, kv) for v, kv in zip(names, kk)) for kk in levels]
    head = " ".join(names)
    L = max(len(head), len("Total"), *(len(x) for x in labs)) + 1
    extra = test != "logrank"
    out.write(f"{'':<{L}}|   Events         Events" + ("        Sum of" if extra else "") + "\n", "text")
    out.write(f"{head:<{L}}|  observed    expected(*)" if stratified else f"{head:<{L}}|  observed       expected",
              "text")
    out.write(("        ranks" if extra else "") + "\n", "text")
    width = 25 + (13 if extra else 0)
    out.write("-" * L + "+" + "-" * width + "\n", "text")
    for i, lab in enumerate(labs):
        out.write(f"{lab:<{L}}|", "text")
        out.write(f"{_g(O_disp[i], '%9.0g'):>10}{E_disp[i]:>15.2f}" +
                  (f"{_g(R[i], '%12.0g'):>13}" if extra else "") + "\n", "result")
    out.write("-" * L + "+" + "-" * width + "\n", "text")
    out.write(f"{'Total':<{L}}|", "text")
    out.write(f"{_g(O_disp.sum(), '%9.0g'):>10}{E_disp.sum():>15.2f}" +
              (f"{_g(R.sum() if abs(R.sum()) > 1e-9 else 0.0, '%12.0g'):>13}" if extra else "") + "\n", "result")
    out.write("\n", "text")
    if stratified:
        svl = " ".join(svars)
        out.write(f"(*) sum over calculations within {svl}\n\n", "text")
    out.write(f"{'chi2(' + str(k) + ') =':>21}", "text")
    out.write(f"{chi2:>11.2f}\n", "result")
    out.write(f"{'Pr>chi2 =':>21}", "text")
    out.write(f"{p:>11.4f}\n", "result")
    s.r = {"df": float(k), "chi2": chi2}


def _unweighted_oe(sd: STData, grp, strat, G):
    O = np.zeros(G)
    E = np.zeros(G)
    for h in np.unique(strat):
        m = (strat == h) & (grp >= 0)
        t, t0, d, w, g = sd.t[m], sd.t0[m], sd.d[m], sd.w[m], grp[m]
        for x in np.unique(t[d > 0]):
            risk = (t0 < x) & (t >= x)
            nt = float(w[risk].sum())
            dsel = (t == x) & (d > 0)
            ng = np.bincount(g[risk], weights=w[risk], minlength=G)
            O += np.bincount(g[dsel], weights=w[dsel], minlength=G)
            E += float(w[dsel].sum()) * ng / nt
    return O, E


# ---------------------------------------------------------------------------
# sts graph
# ---------------------------------------------------------------------------

def _steps(x0: float, times: np.ndarray, y: np.ndarray, y0: float, tmax: float):
    xs, ys = [x0], [y0]
    cur = y0
    for t, v in zip(times, y):
        if t < x0:
            cur = v
            continue
        xs += [float(t), float(t)]
        ys += [cur, float(v)]
        cur = float(v)
    if xs[-1] < tmax:
        xs.append(float(tmax))
        ys.append(cur)
    return xs, ys


def sts_graph(s: "Session", args: str) -> None:
    from ..commands._util import touse
    from ..commands.graphcmd import finish
    from ..graphics.model import Graph, Layer, Panel
    from ..graphics.opts import axis_label, graph_options, split_layer_graph
    from ..graphics.twoway import layer_spec
    from ..lang.syntax import parse_standard
    p = parse_standard(args)
    lopts, graw = split_layer_graph(p.options, {**layer_spec(), "by": 2, "failure": 1, "na": 2,
                                                "cumhaz": 4, "ci": 2, "noshow": 6, "level": 1,
                                                "separate": 3, "censored": 4})
    info = require_st(s)
    if not lopts.get("noshow"):
        header(s, info)
    mask = touse(s, p)
    byv = str(lopts["by"]).split() if isinstance(lopts.get("by"), str) else []
    na = bool(lopts.get("na") or lopts.get("cumhaz"))
    fail = bool(lopts.get("failure"))
    level = float(lopts["level"]) if lopts.get("level") else float(s.settings.get("level", "95"))
    layers = []
    for lab, sd in _by_groups(s, mask, byv):
        k = km(sd.t, sd.t0, sd.d, sd.w)
        if na:
            y, y0 = k.H, 0.0
        elif fail:
            y, y0 = 1 - k.S, 0.0
        else:
            y, y0 = k.S, 1.0
        x0 = float(min(0.0, sd.t0.min())) if len(sd.t0) else 0.0
        tmax = float(sd.t.max()) if len(sd.t) else 0.0
        xs, ys = _steps(x0, k.times, y, y0, tmax)
        if lopts.get("ci") and not na:
            _se, lo, hi = km_ci(k.S, k.gw, level)
            # before the first failure the band is undefined; it starts at the first value
            lo = np.where(lo < SYS, lo, np.nan)
            hi = np.where(hi < SYS, hi, np.nan)
            if fail:
                lo, hi = 1 - hi, 1 - lo
            bx, blo = _steps(x0, k.times, lo, np.nan, tmax)
            _, bhi = _steps(x0, k.times, hi, np.nan, tmax)
            layers.append(Layer("rarea", bx, blo, {"y2": bhi}, {"_pal": len(layers) // 2},
                                f"{level:g}% CI"))
        layers.append(Layer("line", xs, ys, {}, {"_pal": len(layers) // 2} if lopts.get("ci") else {},
                            lab or ("Nelson-Aalen" if na else "Kaplan-Meier")))
    opts = graph_options(graw)
    if na:
        title = "Nelson-Aalen cumulative hazard estimate" + ("s" if byv else "")
    else:
        title = f"Kaplan-Meier {'failure' if fail else 'survival'} estimate" + ("s" if byv else "")
    opts.setdefault("title", {"lines": [title]})
    opts["_xtitle"] = "analysis time"
    opts["_ytitle"] = ""
    if not na and "ylabel" not in graw:
        try:
            opts["ylabel"] = axis_label("0(.25)1, format(%4.2f)")
        except Exception:   # noqa: BLE001
            pass
    if not byv:
        opts.setdefault("legend", {"off": True})
    g = Graph(family="twoway", opts=opts, panels=[Panel("", layers)])
    finish(s, g, graw)
