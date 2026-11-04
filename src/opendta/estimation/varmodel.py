"""var, varsoc and vargranger ([TS] var, varsoc, vargranger): vector autoregressions fitted by
equation-wise least squares, lag-order selection statistics and Granger causality tests."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.stats import chi2_tail
from ..lang.syntax import match_options
from ..lazy import lazy
from .results import CoefRow, Estimates, coef_table, g9, post
from .tsmodels import _check_gaps, _lags, _right, _sample_text
from .tsops import require_ts

if TYPE_CHECKING:
    from ..session import Session

# imported on first use (start-up time; see opendta/lazy.py)
st = lazy("scipy.stats")

SYS = M.SYSMISS


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
    from ..commands._util import touse
    from ..lang.syntax import parse_standard
    from .postest import replay
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
               "XtXi": XtXi, "T": T, "U": U, "X": Xo, "Y": Yo, "Sml": Sml, "t0": t0, "t1": t1,
               "constant": constant, "small": small, "dfk": bool(o.get("dfk"))}
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
            left = f"{lab:<15}= {g9(v, '%9.0g'):>10}"
            _right(out, left, rl, g9(rv, "%9.0g"))
        out.write("\n", "text")
        stat_h = "F" if x["small"] else "chi2"
        out.write(f"{'Equation':<15}{'Parms':>9}{'RMSE':>10}{'R-sq':>9}{stat_h:>10}{'P>' + stat_h:>11}\n", "text")
        out.write("-" * 64 + "\n", "text")
        for dv, mp, rmse, r2, stat, pv in x["eqstats"]:
            out.write(f"{dv[:15]:<15}", "text")
            out.write(f"{mp:>8}{g9(rmse, '%8.0g'):>12}{r2:>9.4f}{g9(stat, '%9.0g'):>11}{pv:>9.4f}\n", "result")
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
        # widths and digits from compat 1101 (six significant digits for FPE and the criteria),
        # each followed by its column for the asterisk
        line = f"  |{i:>3} | {g9(r_[3], '%8.0g'):>8}"
        if i:
            line += f"{g9(LR[i], '%7.0g'):>8}{'*' if lr_best == i else ' '}{K * K:>4}{pv[i]:>7.3f}"
        else:
            line += " " * 20
        line += f"{g9(r_[7], '%8.0g'):>9}{star('fpe', i)}{g9(r_[4], '%8.0g'):>9}{star('aic', i)}" \
                f"{g9(r_[5], '%8.0g'):>9}{star('hqic', i)}{g9(r_[6], '%8.0g'):>9}{star('sbic', i)}"
        out.write(f"{line} |\n", "result")
    out.write("  +" + "-" * 75 + "+\n", "text")
    out.write(f"{'Endogenous:':>14}  {' '.join(deps)}\n", "text")
    out.write(f"{'Exogenous:':>14}  {' '.join(exog + (['_cons'] if constant else []))}\n", "text")
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
            out.write(f" {g9(chi2, '%7.0g'):>7}{len(idx):>6}{pv:>9.3f}    ", "result")
            out.write("|\n", "text")
        out.write("  |--------------------------------------+---------------------------|\n" if j < len(deps) - 1
                  else "  +" + "-" * 66 + "+\n", "text")
