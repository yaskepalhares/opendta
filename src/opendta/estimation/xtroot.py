"""xtunitroot: panel unit-root tests ([XT] xtunitroot): ht, breitung, hadri and fisher.

The statistics are in core/panelroot.py. The layouts follow the manual's summary and carry a
VERIFY mark; the case is in compat 1905. llc and ips need the tabulated moments of the original
papers and are not available yet."""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core import panelroot as pr
from ..core.errors import StataError
from ..core.formats import format_value
from ..lang.syntax import match_options, parse_standard
from .tsmodels import mackinnon_p
from .tsops import panel_info, ts_info

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS
METHODS = ("llc", "ips", "fisher", "hadri", "breitung", "ht")


def _g(x, fmt="%9.4f") -> str:
    return format_value(float(x), fmt, pad=False).strip()


def _matrix(s: "Session", name: str, mask: np.ndarray, require_balanced: bool):
    """Panels as rows of consecutive periods. Returns (Y, list of series for unbalanced use)."""
    ds = s.data
    pvar = panel_info(ds)
    info = ts_info(ds)
    if not pvar or info is None or not info.tvar:
        raise StataError(459, "panel variable not set; use xtset panelvar timevar ...")
    v = ds.get(name).data.astype(float)
    ok = mask & (v < SYS) & (ds.get(pvar).data < SYS) & (ds.get(info.tvar).data < SYS)
    ids = ds.get(pvar).data[ok]
    t = np.rint((ds.get(info.tvar).data[ok] - ds.get(info.tvar).data[ok].min()) / (info.delta or 1.0)).astype(int)
    y = v[ok]
    if len(y) == 0:
        raise StataError(2000, "no observations")
    uniq, inv = np.unique(ids, return_inverse=True)
    series = []
    for g in range(len(uniq)):
        m = inv == g
        o = np.argsort(t[m])
        tg, yg = t[m][o], y[m][o]
        if np.any(np.diff(tg) != 1):
            raise StataError(498, "panels contain gaps")           # VERIFY
        series.append(yg)
    lens = {len(x) for x in series}
    if require_balanced:
        if len(lens) != 1 or len({}) != 0:
            raise StataError(498, "xtunitroot requires strongly balanced data")   # VERIFY
        starts = set()
        for g in range(len(uniq)):
            starts.add(int(t[inv == g].min()))
        if len(starts) != 1:
            raise StataError(498, "xtunitroot requires strongly balanced data")   # VERIFY
        return np.vstack(series), series
    return None, series


def _header(out, title: str, name: str, left: list[str], right: list[str], extra_rule: int = 0) -> None:
    out.write(f"\n{title} for {name}\n", "text")
    out.write("-" * (len(f"{title} for {name}") + extra_rule) + "\n", "text")
    for lft, rgt in zip(left, right):
        out.write(f"{lft:<44}{rgt}".rstrip() + "\n", "text")


def cmd_xtunitroot(s: "Session", args: str) -> None:
    from ..commands._util import touse
    m = re.match(r"^\s*(\w+)\s*(.*)$", args, re.S)
    if not m:
        raise StataError(198, "subcommand required")        # VERIFY
    sub = m.group(1).lower()
    rest = m.group(2)
    full = [x for x in METHODS if x.startswith(sub)]
    if sub in METHODS:
        full = [sub]
    if len(full) != 1:
        raise StataError(198, f"invalid subcommand {m.group(1)}")           # VERIFY
    sub = full[0]
    if sub in ("llc", "ips"):
        raise StataError(198, f"xtunitroot {sub} is not available yet in OpenDTA")
    spec = {"trend": 3, "noconstant": 4, "demean": 3, "lags": 3, "drift": 3, "robust": 3, "dfuller": 2,
            "pperron": 2, "noconstant": 4, "kernel": 3, "altt": 3}
    p = parse_standard(rest)
    o = match_options(p.options, spec) if p.options.strip() else {}
    names = s.expand_varlist(p.varlist)
    if len(names) != 1:
        raise StataError(103 if names else 102, "too many variables specified" if names else "varlist required")
    name = names[0]
    if o.get("pperron") or o.get("kernel") or o.get("altt"):
        raise StataError(198, "option not available yet in OpenDTA")
    if o.get("drift"):
        raise StataError(198, "option drift is not available yet in OpenDTA")
    if sub == "fisher" and not o.get("dfuller") and not o.get("pperron"):
        raise StataError(198, "must specify one of dfuller or pperron")
    mask = touse(s, p)
    out = s.output
    trend, const = bool(o.get("trend")), not o.get("noconstant")
    if trend and not const:
        raise StataError(198, "options trend and noconstant may not be combined")   # VERIFY
    Y, series = _matrix(s, name, mask, sub != "fisher")
    if Y is not None and o.get("demean"):
        Y = Y - Y.mean(axis=0, keepdims=True)
    N = len(series)
    T = Y.shape[1] if Y is not None else None
    means = "Included" if const else "Not included"
    tr = "Included" if trend else "Not included"
    if sub == "ht":
        if trend:
            raise StataError(198, "option trend is not available yet for ht in OpenDTA")
        r = pr.harris_tzavalis(Y, const)
        _header(out, "Harris-Tzavalis unit-root test", name,
                ["Ho: Panels contain unit roots", "Ha: Panels are stationary", "",
                 "AR parameter: Common", f"Panel means:  {means}", f"Time trend:   {tr}"],
                [f"Number of panels  = {N:>6}", f"Number of periods = {T:>6}", "",
                 "Asymptotics: N -> Infinity", "             T Fixed", ""])
        out.write("-" * 78 + "\n", "text")
        out.write("                    Statistic         z         p-value\n", "text")
        out.write("-" * 78 + "\n", "text")
        out.write(f"{' rho':<20}{r['rho']:>8.4f}{r['z']:>14.4f}{r['p']:>13.4f}\n", "result")
        out.write("-" * 78 + "\n", "text")
        s.r = {"N_g": float(N), "N_t": float(T), "rho": r["rho"], "z": r["z"], "p": r["p"]}
    elif sub == "breitung":
        if o.get("lags"):
            raise StataError(198, "option lags() is not available yet for breitung in OpenDTA")
        if not const:
            raise StataError(198, "option noconstant is not available for breitung")    # VERIFY
        r = pr.breitung(Y, trend)
        _header(out, "Breitung unit-root test", name,
                ["Ho: Panels contain unit roots", "Ha: Panels are stationary", "",
                 "AR parameter: Common", f"Panel means:  {means}", f"Time trend:   {tr}"],
                [f"Number of panels  = {N:>6}", f"Number of periods = {T:>6}", "",
                 "Asymptotics:  T,N -> Infinity", " " * 18 + "sequentially", "Prewhitening: Not performed"],
                extra_rule=2)
        out.write("-" * 78 + "\n", "text")
        out.write("                    Statistic      p-value\n", "text")
        out.write("-" * 78 + "\n", "text")
        out.write(f"{' lambda':<20}{r['lambda']:>8.4f}{r['p']:>14.4f}\n", "result")
        out.write("-" * 78 + "\n", "text")
        s.r = {"N_g": float(N), "N_t": float(T), "lambda": r["lambda"], "p": r["p"]}
    elif sub == "hadri":
        r = pr.hadri(Y, trend, bool(o.get("robust")))
        _header(out, "Hadri LM test", name,
                ["Ho: All panels are stationary", "Ha: Some panels contain unit roots", "",
                 f"{'Time trend:':<20}{tr}", f"{'Heteroskedasticity:':<20}" + ("Robust" if o.get("robust") else "Not robust"),
                 f"{'LR variance:':<20}(not used)"],
                [f"Number of panels  = {N:>6}", f"Number of periods = {T:>6}", "",
                 "Asymptotics: T, N -> Infinity", " " * 18 + "sequentially", ""], extra_rule=2)
        out.write("-" * 78 + "\n", "text")
        out.write("                    Statistic      p-value\n", "text")
        out.write("-" * 78 + "\n", "text")
        out.write(f"{' z':<20}{r['z']:>8.4f}{r['p']:>14.4f}\n", "result")
        out.write("-" * 78 + "\n", "text")
        s.r = {"N_g": float(N), "N_t": float(T), "z": r["z"], "p": r["p"]}
    else:
        lags = int(o["lags"]) if o.get("lags") not in (None, False, True) else 0
        ps = []
        for y in series:
            yy = y - 0.0
            if o.get("demean"):
                raise StataError(198, "option demean is not available yet for fisher in OpenDTA")
            try:
                tau, _ = pr.adf_tau(yy, lags, const, trend)
            except ValueError:
                raise StataError(2001, "insufficient observations")
            ps.append(mackinnon_p(tau, "ct" if trend else ("c" if const else "n")))
        r = pr.fisher_combine(np.array(ps))
        tbar = float(np.mean([len(y) for y in series]))
        _header(out, "Fisher-type unit-root test", name,
                ["Based on augmented Dickey-Fuller tests", "-" * 36, "Ho: All panels contain unit roots",
                 "Ha: At least one panel is stationary", "",
                 "AR parameter: Panel-specific", f"Panel means:  {means}", f"Time trend:   {tr}",
                 "Drift term:   Not included"],
                ["", "", f"{'Number of panels':<25}= {N:>6}", f"{'Avg. number of periods':<25}= {tbar:>6.2f}", "",
                 "Asymptotics: T -> Infinity", "", "", f"ADF regressions: {lags} lags"])
        out.write("-" * 78 + "\n", "text")
        out.write("                                  Statistic      p-value\n", "text")
        out.write("-" * 78 + "\n", "text")
        out.write(f" Inverse chi-squared({r['df_P']})".ljust(26) + " P ", "text")
        out.write(f"{r['P']:>12.4f}   {r['p_P']:>9.4f}\n", "result")
        out.write(" Inverse normal".ljust(26) + " Z ", "text")
        out.write(f"{r['Z']:>12.4f}   {r['p_Z']:>9.4f}\n", "result")
        out.write(f" Inverse logit t({r['df_L']})".ljust(26) + " L*", "text")
        out.write(f"{r['L']:>12.4f}   {r['p_L']:>9.4f}\n", "result")
        out.write(" Modified inv. chi-squared".ljust(26) + " Pm", "text")
        out.write(f"{r['Pm']:>12.4f}   {r['p_Pm']:>9.4f}\n", "result")
        out.write("-" * 78 + "\n", "text")
        out.write("P statistic requires number of panels to be finite.\n", "text")
        out.write("Other statistics are suitable for finite or infinite number of panels.\n", "text")
        s.r = {"N_g": float(N), "lags": float(lags), "P": r["P"], "Z": r["Z"], "L": r["L"], "Pm": r["Pm"],
               "p_P": r["p_P"], "p_Z": r["p_Z"], "p_L": r["p_L"], "p_Pm": r["p_Pm"]}
