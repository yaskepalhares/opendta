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
from ..lazy import lazy
from .tsmodels import mackinnon_p
from .tsops import panel_info, ts_info

st = lazy("scipy.stats")

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
    if o.get("altt"):
        raise StataError(198, "option not available yet in OpenDTA")
    if o.get("kernel") and sub != "hadri":
        raise StataError(198, "option kernel() not allowed")                      # VERIFY
    if sub == "fisher" and o.get("lags") in (None, True):
        raise StataError(198, "must specify lags()")                                # compat 1914
    if o.get("drift") and (sub != "fisher" or o.get("pperron") or o.get("trend") or o.get("noconstant")):
        raise StataError(198, "option drift may only be specified with dfuller, without trend or noconstant")  # VERIFY
    if sub == "fisher" and not o.get("dfuller") and not o.get("pperron"):
        raise StataError(198, "must specify one of dfuller or pperron")
    if sub == "hadri" and o.get("kernel") and o.get("robust"):
        raise StataError(198, "cannot specify both kernel() and robust")            # compat 1914
    mask = touse(s, p)
    out = s.output
    trend, const = bool(o.get("trend")), not o.get("noconstant")
    if trend and not const:
        raise StataError(198, "options trend and noconstant may not be combined")   # VERIFY
    Y, series = _matrix(s, name, mask, sub != "fisher" or bool(o.get("demean")))
    if o.get("demean"):
        if Y is None:
            raise StataError(498, "demean requires strongly balanced data")         # VERIFY
        Y = Y - Y.mean(axis=0, keepdims=True)
        series = list(Y)
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
        kname, klags = "", 0
        if o.get("kernel") not in (None, True):
            kw = str(o["kernel"]).split()
            kname = {"ba": "bartlett", "bartlett": "bartlett", "pa": "parzen", "parzen": "parzen",
                     "qs": "quadraticspectral", "quadraticspectral": "quadraticspectral"}.get(kw[0].lower(), "")
            if not kname or len(kw) > 2 or (len(kw) == 2 and not kw[1].isdigit()):
                # VERIFY: the automatic lag choices (nwest, ...) are not available
                raise StataError(198, "kernel() must be bartlett, parzen or quadraticspectral with a number of lags")
            # without a number of lags Stata takes 1 with T = 20 (compat 1914; VERIFY: its rule)
            klags = int(kw[1]) if len(kw) == 2 else 1
        r = pr.hadri(Y, trend, bool(o.get("robust")), kname, klags)
        # with a kernel the heteroskedasticity line says Robust (compat 1914)
        lrv = {"bartlett": "Bartlett kernel", "parzen": "Parzen kernel",
               "quadraticspectral": "Quad. Spectral kernel"}.get(kname, "")
        if kname == "quadraticspectral":
            lrv += f", \n{'':<20}{klags:.2f} lags"
        elif kname:
            lrv += f", {klags} lag{'s' if klags != 1 else ''}"
        _header(out, "Hadri LM test", name,
                ["Ho: All panels are stationary", "Ha: Some panels contain unit roots", "",
                 f"{'Time trend:':<20}{tr}",
                 f"{'Heteroskedasticity:':<20}" + ("Robust" if o.get("robust") or kname else "Not robust"),
                 f"{'LR variance:':<20}" + (lrv if kname else "(not used)")],
                [f"Number of panels  = {N:>6}", f"Number of periods = {T:>6}", "",
                 "Asymptotics: T, N -> Infinity", " " * 18 + "sequentially", ""], extra_rule=2)
        out.write("-" * 78 + "\n", "text")
        out.write("                    Statistic      p-value\n", "text")
        out.write("-" * 78 + "\n", "text")
        out.write(f"{' z':<20}{r['z']:>8.4f}{r['p']:>14.4f}\n", "result")
        out.write("-" * 78 + "\n", "text")
        # in the order of return list (compat 1914)
        xi, zeta2 = (1 / 15, 11 / 6300) if trend else (1 / 6, 1 / 45)
        res = ([("lags", float(klags))] if kname else []) + [
            ("var", zeta2), ("mu", xi), ("p", r["p"]), ("z", r["z"]), ("N_t", float(T)), ("N_g", float(N)),
            ("N", float(N * T)), ("test", "hadri")] + ([("kernel", kname)] if kname else []) + [
            ("deterministics", "trend" if trend else "constant")]
        s.r = dict(reversed(res))
    else:
        pperron = bool(o.get("pperron"))
        drift = bool(o.get("drift"))
        if pperron and o.get("lags") in (None, False, True):
            raise StataError(198, "option lags() required with pperron")            # VERIFY
        lags = int(o["lags"]) if o.get("lags") not in (None, False, True) else 0
        ps = []
        for y in series:
            yy = y - 0.0
            try:
                if pperron:
                    from .tsunit import pp_stats
                    tau = pp_stats(yy, lags, trend, const)[1]
                elif drift:
                    tau, nn = pr.adf_tau(yy, lags, True, False)
                    # dfuller, drift: the t distribution with n - k degrees of freedom (VERIFY)
                    ps.append(float(st.t.cdf(tau, nn - lags - 2)))
                    continue
                else:
                    tau, _ = pr.adf_tau(yy, lags, const, trend)
            except (ValueError, np.linalg.LinAlgError):
                raise StataError(2001, "insufficient observations")
            ps.append(mackinnon_p(tau, "ct" if trend else ("c" if const else "n")))
        r = pr.fisher_combine(np.array(ps))
        # the layout of compat 1914
        based = "Based on Phillips-Perron tests" if pperron else "Based on augmented Dickey-Fuller tests"
        out.write(f"\nFisher-type unit-root test for {name}\n{based}\n" + "-" * (len(based) + (2 if pperron else 0))
                  + "\n", "text")
        lens = [len(y) for y in series]
        per = (f"{'Number of periods':<18}= {lens[0]:>6}" if min(lens) == max(lens) else
               f"{'Avg. number of periods':<23}= {np.mean(lens):>6.2f}")                 # VERIFY: unbalanced
        lab = 17 if pperron else 14
        rows = [("Ho: All panels contain unit roots", f"{'Number of panels':<18}= {N:>6}"),
                ("Ha: At least one panel is stationary", per), ("", ""),
                (f"{'AR parameter:':<{lab}}Panel-specific", "Asymptotics: T -> Infinity"),
                (f"{'Panel means:':<{lab}}{means}", ""),
                (f"{'Time trend:':<{lab}}{tr}", "Cross-sectional means removed" if o.get("demean") else "")]
        ls = "" if lags == 1 else "s"
        rows.append((f"Newey-West lags: {lags} lag{ls}", "") if pperron else
                    (f"{'Drift term:':<{lab}}{'Included' if drift else 'Not included'}", f"ADF regressions: {lags} lag{ls}"))
        for lft, rgt in rows:
            out.write(f"{lft:<44}{rgt}".rstrip() + "\n", "text")
        out.write("-" * 78 + "\n" + " " * 34 + "Statistic      p-value\n" + "-" * 78 + "\n", "text")
        for lab_, sym, key in ((f" Inverse chi-squared({r['df_P']})", "P", "P"), (" Inverse normal", "Z", "Z"),
                               (f" Inverse logit t({r['df_L']})", "L*", "L"),
                               (" Modified inv. chi-squared", "Pm", "Pm")):
            out.write(f"{lab_:<27}{sym:<3}", "text")
            out.write(f"{r[key]:>13.4f}{r['p_' + key]:>13.4f}\n", "result")
        out.write("-" * 78 + "\n", "text")
        out.write(" P statistic requires number of panels to be finite.\n"
                  " Other statistics are suitable for finite or infinite number of panels.\n" + "-" * 78 + "\n", "text")
        opts_ = [w for w in ("noconstant", "trend", "drift", "demean") if o.get(w)] + [f"lags({lags})"]   # VERIFY: order
        # return list shows them from the last stored to the first (compat 1914)
        res = [("N", float(sum(lens))), ("N_t", float(max(lens))), ("N_g", float(N)), ("p_Pm", r["p_Pm"]),
               ("Pm", r["Pm"]), ("df_L", float(r["df_L"])), ("p_L", r["p_L"]), ("L", r["L"]), ("p_Z", r["p_Z"]),
               ("Z", r["Z"]), ("df_P", float(r["df_P"])), ("p_P", r["p_P"]), ("P", r["P"]),
               ("options", " ".join(opts_)), ("urtest", "pperron" if pperron else "dfuller"), ("test", "fisher")]
        s.r = dict(reversed(res))
