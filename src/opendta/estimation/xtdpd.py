"""Dynamic panel-data estimation: xtabond, xtdpdsys and xtdpd ([XT] xtabond, [XT] xtdpdsys,
[XT] xtdpd), with estat abond and estat sargan.

The three commands share one engine: an instrument specification (GMM-type and standard
instruments for the differenced equation and, in the system estimators, for the level equation)
built per panel on the grid of periods. The algebra is in core/panelgmm.py. The layouts follow the
manual's examples where they are shown; the rest carries a VERIFY mark (compat 1905)."""

from __future__ import annotations

import math
import re
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core import panelgmm as pg
from ..core.errors import StataError
from ..core.formats import format_value
from ..core.stats import chi2_tail
from ..lang.syntax import match_options
from . import fvars
from .regress import Sample, cmdline, comma
from .results import CoefRow, Estimates, coef_table, post
from .tsops import panel_info, ts_info

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def _g(x, fmt="%10.0g") -> str:
    if x is None or not np.isfinite(x) or abs(x) >= SYS:
        return "."
    return format_value(float(x), fmt, pad=False).strip()


def _varlist(s: "Session", txt) -> list[str]:
    if txt in (None, False, True) or not str(txt).strip():
        return []
    body = str(txt).split(",")[0]
    return s.expand_varlist(body)


def _lagrange(txt, default: tuple[int, int | None]) -> tuple[int, int | None]:
    m = re.search(r"lagrange\(\s*(\d+)\s*(\d+|\.)?\s*\)", str(txt))
    if not m:
        return default
    a = int(m.group(1))
    b = m.group(2)
    return a, (a if b is None else (None if b == "." else int(b)))


def _lab(a: int, b: int | None, name: str, diff: bool = False) -> str:
    d = "D." if diff else ""
    return f"L({a}/{'.' if b is None else b}).{d}{name}" if (b is None or b != a) else (f"L{a}.{d}{name}" if a > 1 else f"L.{d}{name}")


class _Spec:
    """Instruments and regressors in terms of series on the grid of periods."""

    def __init__(self):
        self.dgmm: list[tuple[str, int, int | None]] = []     # (series key, first lag, last lag)
        self.div: list[str] = []
        self.lgmm: list[tuple[str, int, int | None]] = []
        self.liv: list[str] = []
        self.level = False
        self.cons = True


def cmd_xtabond(s: "Session", args: str) -> None:
    _run(s, "xtabond", args)


def cmd_xtdpdsys(s: "Session", args: str) -> None:
    _run(s, "xtdpdsys", args)


def cmd_xtdpd(s: "Session", args: str) -> None:
    _run(s, "xtdpd", args)


def _run(s: "Session", cmd: str, args: str) -> None:
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, cmd, args)
        return
    spec_o = {"twostep": 3, "noconstant": 4, "inst": 3, "pre": 3, "endogenous": 4, "lags": 3, "maxldep": 4,
              "maxlags": 4, "diffvars": 5, "teffects": 2, "small": 2, "hascons": 4, "nolevel": 4,
              "dgmmiv": 4, "lgmmiv": 4, "div": 3, "liv": 3, "h": 1, "artests": 3, "coldiag": 4, "vce": 3}
    if cmd == "xtdpd":
        spec_o = {k: v for k, v in spec_o.items() if k in ("twostep", "noconstant", "dgmmiv", "lgmmiv", "div", "liv",
                                                           "h", "artests", "hascons", "vce")}
    else:
        spec_o = {k: v for k, v in spec_o.items() if k not in ("dgmmiv", "lgmmiv", "div", "liv", "h")}
    # first pass: the option names without the model, to read the variable lists
    m = re.match(r"^(.*?)(,.*)?$", args, re.S)
    body, optxt = m.group(1), (m.group(2) or "")
    o0 = match_options(optxt[1:], dict(spec_o, robust=1, level=1)) if optxt.strip(" ,") else {}
    for bad in ("diffvars", "teffects", "small", "hascons", "coldiag", "artests", "h"):
        if o0.get(bad):
            raise StataError(198, f"option {bad} is not available yet in OpenDTA")
    pre = _varlist(s, o0.get("pre"))
    endo = _varlist(s, o0.get("endogenous"))
    inst = _varlist(s, o0.get("inst"))
    extra = " ".join(pre + endo + inst + _varlist(s, o0.get("dgmmiv")) + _varlist(s, o0.get("lgmmiv"))
                     + _varlist(s, o0.get("div")) + _varlist(s, o0.get("liv")))
    opts = {k: v for k, v in spec_o.items() if k not in ("noconstant", "vce")}
    smp = Sample(s, args, opts_spec=opts, weights=(), extra_vars=extra)
    o = smp.o
    if smp.vce == "cluster":
        raise StataError(198, "option cluster() is not available yet in OpenDTA")
    robust = smp.vce == "robust"
    twostep = bool(o.get("twostep"))
    ds = s.data
    pvar, info = panel_info(ds), ts_info(ds)
    if not pvar or info is None or not info.tvar:
        raise StataError(459, "panel variable not set; use xtset panelvar timevar ...")
    ids = ds.get(pvar).data[smp.mask]
    uniq, inv = np.unique(ids, return_inverse=True)
    G = len(uniq)
    tv = ds.get(info.tvar).data[smp.mask].astype(float)
    # the time axis starts at the first period of the data, not at the first period of the sample
    # (the sample loses the first period when a lag is a regressor)
    t_full = ds.get(info.tvar).data.astype(float)
    t_ok = (t_full < SYS) & (ds.get(pvar).data < SYS)
    t0 = float(t_full[t_ok].min())
    tix = np.rint((tv - t0) / (info.delta or 1.0)).astype(int)
    Tg = int(np.rint((t_full[t_ok].max() - t0) / (info.delta or 1.0))) + 1

    def grid(values):
        g = np.full((G, Tg), np.nan)
        g[inv, tix] = values
        return g
    y = grid(smp.y.astype(float))
    xcols = [c for c in smp.cols if c.values is not None]
    series: dict[str, np.ndarray] = {}
    reg_names: list[str] = []
    reg_series: list[np.ndarray] = []

    # instruments and lagged series come from every observation of the panel, not only from the
    # estimation sample (L.y is missing in the first period, which must still be a valid instrument)
    ids_all, t_all = ds.get(pvar).data, ds.get(info.tvar).data.astype(float)
    pos = np.searchsorted(uniq, ids_all)
    pos_ok = (pos < G) & (uniq[np.minimum(pos, G - 1)] == ids_all) & (ids_all < SYS) & (t_all < SYS)
    tix_all = np.rint((t_all - t0) / (info.delta or 1.0))
    pos_ok &= (tix_all >= 0) & (tix_all < Tg)

    def raw(name):
        if name not in series:
            g_ = np.full((G, Tg), np.nan)
            v_ = ds.get(name).data.astype(float)
            ok_ = pos_ok & (v_ < SYS)
            g_[pos[ok_], tix_all[ok_].astype(int)] = v_[ok_]
            series[name] = g_
        return series[name]

    def shift(a, j):
        out = np.full_like(a, np.nan)
        out[:, j:] = a[:, :Tg - j] if j else a
        return out
    sp = _Spec()
    sp.cons = not o.get("noconstant")
    lab_dgmm: list[str] = []
    lab_div: list[str] = []
    lab_lgmm: list[str] = []
    lab_liv: list[str] = []
    dname = smp.depname
    series[dname] = raw(dname)
    p = 1
    if cmd == "xtdpd":
        for c in xcols:
            reg_names.append(c.name)
            reg_series.append(grid(c.values))
        dg, lg = o0.get("dgmmiv"), o0.get("lgmmiv")
        for v in _varlist(s, dg):
            a, b = _lagrange(dg, (2, None))
            sp.dgmm.append((v, a, b))
            lab_dgmm.append(_lab(a, b, v))
        for v in _varlist(s, lg):
            a, b = _lagrange(lg, (1, 1))
            sp.lgmm.append((v, a, b))
            lab_lgmm.append(_lab(a, b, v, True))
        for v in _varlist(s, o0.get("div")):
            sp.div.append(v)
            lab_div.append("D." + v)
        for v in _varlist(s, o0.get("liv")):
            sp.liv.append(v)
            lab_liv.append(v)
        sp.level = bool(sp.lgmm or sp.liv)
        for v in sp.dgmm:
            raw(v[0])
        title = "Dynamic panel-data estimation"
    else:
        p = int(o["lags"]) if o.get("lags") not in (None, False, True) else 1
        if p < 1:
            raise StataError(198, "lags() must be positive")                # VERIFY
        maxldep = int(o["maxldep"]) if o.get("maxldep") not in (None, False, True) else None
        maxlags = int(o["maxlags"]) if o.get("maxlags") not in (None, False, True) else None
        for j in range(1, p + 1):
            reg_names.append(f"L{j}.{dname}" if j > 0 else dname)
            reg_series.append(shift(y, j))
        sp.dgmm.append((dname, 2, None if maxldep is None else maxldep + 1))
        lab_dgmm.append(_lab(2, None if maxldep is None else maxldep + 1, dname))
        for v in pre:
            reg_names.append(v)
            reg_series.append(raw(v))
            sp.dgmm.append((v, 1, maxlags))
            lab_dgmm.append(_lab(1, maxlags, v))
        for v in endo:
            reg_names.append(v)
            reg_series.append(raw(v))
            sp.dgmm.append((v, 2, None if maxlags is None else maxlags + 1))
            lab_dgmm.append(_lab(2, None if maxlags is None else maxlags + 1, v))
        # VERIFY: predetermined and endogenous variables come before the strictly exogenous ones
        for c in xcols:
            reg_names.append(c.name)
            reg_series.append(grid(c.values))
            sp.div.append(c.name)
            series[c.name] = reg_series[-1]
            lab_div.append("D." + c.name)
        for v in inst:
            sp.div.append(v)
            raw(v)
            lab_div.append("D." + v)
        if cmd == "xtdpdsys":
            sp.level = True
            sp.lgmm.append((dname, 1, 1))
            lab_lgmm.append(f"LD.{dname}")
            title = "System dynamic panel-data estimation"
        else:
            title = "Arellano-Bond dynamic panel-data estimation"
    k = len(reg_names) + (1 if sp.cons else 0)
    # ---- equations and instruments per panel
    gmm_keys: dict[tuple, int] = {}
    std_cols: list[str] = []
    nstd = len(sp.div) + len(sp.liv) + (1 if sp.cons else 0)
    rows_X, rows_y, rows_Z, rows_H, rows_t = [], [], [], [], []
    sparse: list[list[tuple[int, int, float]]] = []     # (row, col, value) per panel
    samp_counts: list[int] = []                         # observations per panel as the header reports them
    for g in range(G):
        X_i, y_i, t_i, ent = [], [], [], []
        lvl_rows = []
        for t in range(1, Tg):
            vals = [r[g, t] - r[g, t - 1] for r in reg_series]
            dy = y[g, t] - y[g, t - 1]
            if np.isnan(dy) or any(np.isnan(v) for v in vals):
                continue
            row = vals + ([0.0] if sp.cons else [])
            X_i.append(row)
            y_i.append(dy)
            t_i.append(t)
            e = []
            for key, (v, a, b) in enumerate(sp.dgmm):
                arr = series[v]
                hi = t if b is None else min(b, t)
                for l in range(a, hi + 1):
                    val = arr[g, t - l]
                    if not np.isnan(val):
                        col = gmm_keys.setdefault(("d", key, t, l), len(gmm_keys))
                        e.append((col, val))
            ent.append(e)
        nd = len(X_i)
        lvl_t = []
        n_samp = nd
        if sp.level or sp.cons:
            n_samp = 0
            # VERIFY: without a level equation the constant gets its own level equation, instrumented by
            # the constant alone (it is listed under "Instruments for level equation" in the 1905 log)
            for t in range(1, Tg):
                vals = [r[g, t] for r in reg_series]
                if np.isnan(y[g, t]) or any(np.isnan(v) for v in vals):
                    continue
                n_samp += 1
                inst_ok = True
                e = []
                for key, (v, a, b) in enumerate(sp.lgmm):
                    arr = series[v] if v in series else raw(v)
                    hi = t - 1 if b is None else min(b, t - 1)
                    got = False
                    for l in range(a, hi + 1):
                        if t - l - 1 >= 0:
                            val = arr[g, t - l] - arr[g, t - l - 1]
                            if not np.isnan(val):
                                col = gmm_keys.setdefault(("l", key, t, l), len(gmm_keys))
                                e.append((col, val))
                                got = True
                    inst_ok = inst_ok and got
                X_i.append(vals + ([1.0] if sp.cons else []))
                y_i.append(y[g, t])
                lvl_t.append(t)
                ent.append(e)
        if not X_i:
            continue
        rows_X.append(np.array(X_i, dtype=float))
        rows_y.append(np.array(y_i, dtype=float))
        rows_t.append(np.array(t_i + [-1] * len(lvl_t)))
        # standard instruments
        std = np.zeros((len(X_i), nstd))
        for r, t in enumerate(t_i):
            j = 0
            for v in sp.div:
                a = series[v] if v in series else raw(v)
                d = a[g, t] - a[g, t - 1]
                std[r, j] = 0.0 if np.isnan(d) else d
                j += 1
            j += len(sp.liv)
        for r2, t in enumerate(lvl_t):
            r = nd + r2
            j = len(sp.div)
            for v in sp.liv:
                a = series[v] if v in series else raw(v)
                std[r, j] = 0.0 if np.isnan(a[g, t]) else a[g, t]
                j += 1
            if sp.cons:
                std[r, j] = 1.0
        sparse.append((ent, std, nd, len(X_i)))
        samp_counts.append(nd if cmd == "xtabond" else n_samp)
    if not rows_X:
        raise StataError(2001, "insufficient observations")
    ncols = len(gmm_keys)
    Zs, Hs = [], []
    for (ent, std, nd, n_i), t_i in zip(sparse, rows_t):
        Zg = np.zeros((n_i, ncols))
        for r, e in enumerate(ent):
            for col, val in e:
                Zg[r, col] = val
        Zs.append(np.hstack([Zg, std]))
        Hd = pg.tridiagonal(t_i[:nd])
        H = np.zeros((n_i, n_i))
        H[:nd, :nd] = Hd
        if n_i > nd:
            H[nd:, nd:] = np.eye(n_i - nd)
        Hs.append(H)
    # drop instruments that are zero in every panel
    keep = np.flatnonzero(np.any([np.any(z != 0, axis=0) for z in Zs], axis=0))
    Zs = [z[:, keep] for z in Zs]
    X, yv = rows_X, rows_y
    if len(keep) < k:
        raise StataError(481, "equation not identified; fewer instruments than parameters")   # VERIFY
    fit = pg.gmm_fit(X, yv, Zs, Hs, twostep, robust, nd=[sp_[2] for sp_ in sparse])
    b, V = fit["b"], fit["V"]
    Gu = len(X)
    # VERIFY: xtabond reports the rows of the differenced equation, xtdpd and xtdpdsys the observations of
    # the sample (1905 log)
    counts = np.array(samp_counts)
    nobs = int(counts.sum())
    names = list(reg_names) + (["_cons"] if sp.cons else [])
    slopes = list(range(len(reg_names)))
    try:
        chi2 = float(b[slopes] @ np.linalg.solve(V[np.ix_(slopes, slopes)], b[slopes]))
    except np.linalg.LinAlgError:
        chi2 = SYS
    pchi = chi2_tail(chi2, len(slopes)) if chi2 < SYS else SYS
    # table rows: the lags of the dependent variable under its name, as Stata does
    se = np.sqrt(np.maximum(np.diag(V), 0))
    rows: list[CoefRow] = []
    last_var = None
    for j, nm in enumerate(names):
        mm = re.match(r"^L(\d*)\.(.+)$", nm)
        var, lab = (mm.group(2), f"L{mm.group(1) or 1}.") if mm else (None, nm)
        if var is not None and var != last_var:
            rows.append(CoefRow("header", var))
        elif var is None and last_var is not None:
            rows.append(CoefRow("blank"))
        rows.append(CoefRow("coef", lab, b[j], se[j]))
        last_var = var
    est = Estimates(cmd, dname, names, b, V, nobs, "z", None)
    est.rows = rows
    est.eqnames = [dname] * len(names)
    est.level, est.vcetype, est.constant = smp.level, "", sp.cons
    est.cols, est.terms, est.depcomp = [], [], smp.depcomp
    L = Zs[0].shape[1]
    x = dict(N=nobs, G=Gu, tmin=int(counts.min()), tavg=float(counts.mean()), tmax=int(counts.max()), L=L,
             chi2=chi2, p=pchi, df=len(slopes), twostep=twostep, robust=robust, title=title, pvar=pvar, cmd=cmd,
             tvar=info.tvar, dgmm=lab_dgmm, div=lab_div, lgmm=lab_lgmm, liv=lab_liv, cons=sp.cons, level=sp.level)
    est.extra = x
    est.display = _display
    est.estat = {"abond": (2, lambda ss, ee, rest: _estat_abond(ss, ee, rest, fit, X, Zs, rows_t)),
                 "sargan": (3, lambda ss, ee, rest: _estat_sargan(ss, ee, fit, Zs, twostep, robust))}
    sc = [("N", float(nobs)), ("N_g", float(Gu)), ("df_m", float(len(slopes))), ("g_min", float(counts.min())),
          ("g_avg", float(counts.mean())), ("g_max", float(counts.max())), ("t_min", float(Tg)), ("chi2", chi2),
          ("sig2", fit["sigma2"]), ("zrank", float(L)), ("rank", float(k))]
    mc = [("cmdline", cmdline(cmd, args)), ("cmd", cmd), ("depvar", dname), ("ivar", pvar), ("tvar", info.tvar),
          ("vce", "robust" if robust else "gmm"), ("twostep", "twostep" if twostep else "onestep"),
          ("title", title), ("properties", "b V"), ("estat_cmd", "xtdpd_estat"), ("predict", "xtdpd_p")]
    post(s, est, list(reversed(sc)), list(reversed(mc)), smp.mask)
    _display(s, est)


def _display(s, est, *, header=True, table=True, level=None, **_kw):
    out, x = s.output, est.extra
    if header:
        out.write("\n", "text")
        out.write(f"{x['title']:<48}{'Number of obs':<18}= ", "text")
        out.write(f"{comma(x['N']):>10}\n", "result")
        out.write(f"{'Group variable: ' + x['pvar']:<48}{'Number of groups':<18}= ", "text")
        out.write(f"{comma(x['G']):>10}\n", "result")
        out.write(f"Time variable: {x['tvar']}\n\n", "text")
        out.write(f"{'':<48}Obs per group:\n", "text")
        for lab, v in (("min", _g(x["tmin"])), ("avg", _g(x["tavg"])), ("max", _g(x["tmax"]))):
            out.write(f"{'':<62}{lab} = ", "text")
            out.write(f"{v:>10}\n", "result")
        out.write("\n", "text")
        out.write(f"{'Number of instruments = ' + format(x['L'], '6d'):<48}"
                  f"{'Wald chi2(' + str(x['df']) + ')':<18}= ", "text")
        out.write(f"{_g(x['chi2'], '%10.2f'):>10}\n", "result")
        out.write(f"{'':<48}{'Prob > chi2':<18}= ", "text")
        out.write(f"{_g(x['p'], '%10.4f'):>10}\n", "result")
        out.write(("Two-step results" if x["twostep"] else "One-step results") + "\n", "text")
        if x["robust"] and x["cmd"] == "xtabond":
            out.write(f"{'(Std. Err. adjusted for clustering on ' + x['pvar'] + ')':>78}\n", "text")   # VERIFY
    if not table:
        return
    if x["twostep"] and not x["robust"]:
        out.write("Warning: gmm two-step standard errors are biased; robust standard errors are recommended.\n",
                  "text")        # VERIFY
    coef_table(s, est.depvar, est.rows, stat="z", df=None, level=level or est.level,
               vcetype="WC-Robust" if (x["twostep"] and x["robust"]) else ("Robust" if x["robust"] else ""),
               footer=False)
    out.write("-" * 78 + "\n", "text")
    out.write("Instruments for differenced equation\n", "text")
    if x["dgmm"]:
        out.write("        GMM-type: " + " ".join(x["dgmm"]) + "\n", "text")
    if x["div"]:
        out.write("        Standard: " + " ".join(x["div"]) + "\n", "text")
    if x["level"] or x["cons"]:
        out.write("Instruments for level equation\n", "text")
        if x["lgmm"]:
            out.write("        GMM-type: " + " ".join(x["lgmm"]) + "\n", "text")
        std = x["liv"] + (["_cons"] if x["cons"] else [])
        if std:
            out.write("        Standard: " + " ".join(std) + "\n", "text")


def _estat_abond(s, est, rest, fit, X, Z, times) -> None:
    m = re.search(r"artests\(\s*(\d+)\s*\)", rest)
    top = int(m.group(1)) if m else 2
    x_ = est.extra
    if (x_["level"] or x_["cons"]) and not x_["twostep"] and not x_["robust"]:
        raise StataError(321, "artests not computed for one-step system estimator with vce(gmm)")   # VERIFY
    out = s.output
    out.write("\nArellano-Bond test for zero autocorrelation in first-differenced errors\n", "text")
    out.write("  +-----------------------+\n  |Order |  z     Prob > z|\n  |------+----------------|\n", "text")
    res = {}
    for j in range(1, top + 1):
        z, _ = pg.ar_test(j, fit, X, Z, times)
        pv = float(2 * (0.5 * math.erfc(abs(z) / math.sqrt(2)))) if np.isfinite(z) else SYS
        out.write(f"  |{j:>4}  |{_g(z, '%7.0g'):>7}  {pv:>6.4f} |\n", "result")
        res[j] = (z, pv)
    out.write("  +-----------------------+\n   H0: no autocorrelation \n", "text")
    s.r = {}
    for j, (z, pv) in res.items():
        s.r[f"ar{j}"] = z
        s.r[f"ar{j}_p"] = pv


def _estat_sargan(s, est, fit, Z, twostep, robust) -> None:
    if robust:
        raise StataError(321, "estat sargan is not valid after vce(robust)")       # VERIFY
    stat, df = pg.overid_test(fit, Z, twostep)
    pv = chi2_tail(stat, df) if df > 0 else SYS
    out = s.output
    out.write("\nSargan test of overidentifying restrictions\n", "text")
    out.write("        H0: overidentifying restrictions are valid\n\n", "text")
    out.write(f"        {'chi2(' + str(df) + ')':<13}= ", "text")
    out.write(f"{format(stat, '.6g'):>9}\n", "result")
    out.write(f"        Prob > chi2  = ", "text")
    out.write(f"{_g(pv, '%9.4f'):>9}\n", "result")
    s.r = {"chi2": stat, "df": float(df), "p": pv}
