"""Unit-root tests: dfgls and pperron ([TS] dfgls, [TS] pperron).

dfgls: GLS detrending (Elliott, Rothenberg and Stock 1996) and the augmented Dickey-Fuller
regression without deterministic terms on the detrended series, for every lag up to the maximum;
the lag is chosen by the Ng-Perron sequential t, the Schwarz criterion and the modified Akaike
criterion. Critical values are interpolated from the ERS tables (VERIFY: Stata's default 5% and
10% values come from a response surface that the manual does not give).

pperron: Phillips-Perron Z(rho) and Z(t) with the Newey-West long-run variance; critical values
interpolated from Fuller's tables, p-value from MacKinnon (1994).
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..lang.syntax import match_options, parse_standard
from ..lazy import lazy
from .results import CoefRow, coef_table
from .tsmodels import mackinnon_p

st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS

# ERS (1996), Table 1, trend case (c = -13.5): N = 50, 100, 200
_ERS_TREND = {1: [-3.77, -3.58, -3.46], 5: [-3.19, -3.03, -2.93], 10: [-2.89, -2.74, -2.64]}
# mean-only case (c = -7); VERIFY: the manual does not print this table, the asymptotic
# values are used for every N
_ERS_MEAN = {1: [-2.58, -2.58, -2.58], 5: [-1.95, -1.95, -1.95], 10: [-1.62, -1.62, -1.62]}
_ERS_N = [50, 100, 200]


def ers_critical(trend: bool, n: int) -> list[float]:
    tab = _ERS_TREND if trend else _ERS_MEAN
    out = []
    for lev in (1, 5, 10):
        v = tab[lev]
        if n <= _ERS_N[0]:
            out.append(v[0])
        elif n >= _ERS_N[-1]:
            out.append(v[-1])
        else:
            i = 0 if n < _ERS_N[1] else 1
            n1 = _ERS_N[i]
            out.append(v[i] + (n - n1) / n1 * (v[i + 1] - v[i]))
    return out


def gls_detrend(y: np.ndarray, trend: bool = True) -> np.ndarray:
    T = len(y)
    c = -13.5 if trend else -7.0
    a = 1 + c / T
    ye = np.concatenate([[y[0]], y[1:] - a * y[:-1]])
    x = np.concatenate([[1.0], np.full(T - 1, 1 - a)])
    cols = [x]
    if trend:
        t = np.arange(1, T + 1, dtype=float)
        cols.append(np.concatenate([[1.0], t[1:] - a * t[:-1]]))
    Z = np.column_stack(cols)
    d, *_ = np.linalg.lstsq(Z, ye, rcond=None)
    fit = d[0] + (d[1] * np.arange(1, T + 1) if trend else 0.0)
    return y - fit


def adf_gls(ystar: np.ndarray, k: int, kmax: int):
    """Regression of dy*_t on y*_{t-1} and k lagged differences, no constant, on the sample
    t = kmax+1..T. Returns (coefficients, standard errors, rmse, N, sum y*_{t-1}^2)."""
    T = len(ystar)
    dy = np.diff(ystar)                       # dy[i] = y*_{i+1} - y*_i  (period i+2 in 1..T)
    rows = np.arange(kmax, T - 1)             # dy index of the dependent variable (VERIFY: N = T - 1 - kmax)
    y = dy[rows]
    cols = [ystar[rows]]                      # y*_{t-1}
    for j in range(1, k + 1):
        cols.append(dy[rows - j])
    X = np.column_stack(cols)
    n, p = X.shape
    b, *_ = np.linalg.lstsq(X, y, rcond=None)
    e = y - X @ b
    s2 = float(e @ e) / (n - p)
    V = s2 * np.linalg.inv(X.T @ X)
    return b, np.sqrt(np.diag(V)), math.sqrt(s2), n, float(np.sum(X[:, 0] ** 2))


def cmd_dfgls(s: "Session", args: str) -> None:
    from ..commands.tsspectral import _sample_series
    p = parse_standard(args)
    o = match_options(p.options, {"maxlag": 4, "notrend": 3, "ers": 3}) if p.options.strip() else {}
    names = s.expand_varlist(p.varlist) if p.varlist.strip() else []
    if len(names) != 1:
        raise StataError(103 if names else 102, "too many variables specified" if names else "too few variables specified")
    y, _rows = _sample_series(s, p, p.varlist.strip())
    T0 = len(y)
    trend = not o.get("notrend")
    schwert = o.get("maxlag") is None
    kmax = int(math.floor(12 * (T0 / 100) ** 0.25)) if schwert else int(float(s.eval(str(o["maxlag"]))))
    if kmax < 0 or kmax > T0 - 4:
        raise StataError(198, "maxlag() out of range")      # VERIFY
    ystar = gls_detrend(y, trend)
    res = []
    for k in range(0, kmax + 1):
        b, se, rmse, n, sy2 = adf_gls(ystar, k, kmax)
        tau = float(b[0] / se[0])
        tlast = float(b[k] / se[k]) if k >= 1 else math.nan
        res.append({"k": k, "tau": tau, "rmse": rmse, "n": n, "b0": float(b[0]), "tlast": tlast, "sy2": sy2})
    N = res[0]["n"]
    for r in res:
        r["sc"] = math.log(r["rmse"] ** 2) + (r["k"] + 1) * math.log(N) / N
        tk = r["b0"] ** 2 * r["sy2"] / r["rmse"] ** 2
        r["maic"] = math.log(r["rmse"] ** 2) + 2 * (tk + r["k"]) / N
    zc = float(st.norm.ppf(0.95))
    opt = 0
    for r in reversed(res):
        if r["k"] >= 1 and abs(r["tlast"]) > zc:
            opt = r["k"]
            break
    sc_best = min(res, key=lambda r: r["sc"])
    maic_best = min(res, key=lambda r: r["maic"])
    crit = ers_critical(trend, T0)      # VERIFY: interpolation uses the number of observations
    out = s.output
    from ..core.formats import format_value

    def g7(x: float) -> str:
        return format_value(x, "%9.7g", pad=False).strip()
    name = p.varlist.strip()
    out.write(f"\nDF-GLS for {name:<46}Number of obs = ", "text")
    out.write(f"{N:>5}\n\n", "result")
    lab = "tau" if trend else "mu"
    out.write(f"{'':15}{'DF-GLS ' + lab:<16}1% Critical       5% Critical      10% Critical\n", "text")
    out.write(f"  [lags]     Test Statistic        Value             Value             Value\n", "text")
    out.write("-" * 78 + "\n", "text")
    for r in reversed(res):
        if r["k"] == 0 and kmax > 0:
            continue
        out.write(f"{r['k']:>5}{r['tau']:>18.3f}{crit[0]:>17.3f}{crit[1]:>18.3f}{crit[2]:>18.3f}\n", "result")
    out.write("\n", "text")
    if kmax > 0:
        rm = next(r["rmse"] for r in res if r["k"] == opt)
        out.write(f"Opt Lag (Ng-Perron seq t) = {opt:>2} with RMSE {g7(rm):>9}\n", "text")
        out.write(f"Min SC   = {g7(sc_best['sc']):>9} at lag {sc_best['k']:>2} with RMSE {g7(sc_best['rmse']):>9}\n", "text")
        out.write(f"Min MAIC = {g7(maic_best['maic']):>9} at lag {maic_best['k']:>2} with RMSE {g7(maic_best['rmse']):>9}\n", "text")
    from ..commands.matrix import Matrix
    if kmax == 0:
        s.r = {"rmse0": res[0]["rmse"], "dft0": res[0]["tau"]}
    else:
        tab = np.array([[r["k"], r["maic"], r["sc"], r["rmse"], r["tau"]] for r in res])
        s.r = {"maxlag": float(kmax), "N": float(N), "sclag": float(sc_best["k"]),
               "maiclag": float(maic_best["k"]), "optlag": float(opt),
               "results": Matrix(tab, [f"r{i + 1}" for i in range(len(res))],
                                 ["k", "MAIC", "SIC", "RMSE", "DFGLS"])}


# ---------------------------------------------------------------------------
# pperron
# ---------------------------------------------------------------------------

_FULLER_N = [25, 50, 100, 250, 500]
# Fuller (1976), Table 8.5.1: n(rho - 1), by 1%, 5%, 10%
_FULLER_RHO = {
    "n": {1: [-11.9, -12.9, -13.3, -13.6, -13.7], 5: [-7.3, -7.7, -7.9, -8.0, -8.0],
          10: [-5.3, -5.5, -5.6, -5.7, -5.7]},
    "c": {1: [-17.2, -18.9, -19.8, -20.3, -20.5], 5: [-12.5, -13.3, -13.7, -14.0, -14.0],
          10: [-10.2, -10.7, -11.0, -11.2, -11.2]},
    "ct": {1: [-22.5, -25.8, -27.4, -28.4, -28.9], 5: [-17.9, -19.8, -20.7, -21.3, -21.5],
           10: [-15.6, -16.8, -17.5, -18.0, -18.1]},
}
# Fuller (1976), Table 8.5.2: tau
_FULLER_TAU = {
    "n": {1: [-2.66, -2.62, -2.60, -2.58, -2.58], 5: [-1.95] * 5, 10: [-1.60, -1.61, -1.61, -1.62, -1.62]},
    "c": {1: [-3.75, -3.58, -3.51, -3.46, -3.44], 5: [-3.00, -2.93, -2.89, -2.88, -2.87],
          10: [-2.63, -2.60, -2.58, -2.57, -2.57]},
    "ct": {1: [-4.38, -4.15, -4.04, -3.99, -3.98], 5: [-3.60, -3.50, -3.45, -3.43, -3.42],
           10: [-3.24, -3.18, -3.15, -3.13, -3.13]},
}


def _interp(tab: dict, kind: str, n: int) -> list[float]:
    out = []
    for lev in (1, 5, 10):
        v = tab[kind][lev]
        if n <= _FULLER_N[0]:
            out.append(v[0])
        elif n >= _FULLER_N[-1]:
            out.append(v[-1])
        else:
            i = max(j for j in range(len(_FULLER_N) - 1) if _FULLER_N[j] < n)
            a, b = _FULLER_N[i], _FULLER_N[i + 1]
            out.append(v[i] + (n - a) / (b - a) * (v[i + 1] - v[i]))
    return out


def pp_stats(y: np.ndarray, q: int, trend: bool, constant: bool, t0: np.ndarray | None = None):
    """Phillips-Perron statistics of the regression y_t = a + rho y_{t-1} (+ d t) + u_t."""
    n = len(y) - 1
    cols = [y[:-1]]
    if trend:
        cols.append(t0[1:] if t0 is not None else np.arange(1, n + 1, dtype=float))
    if constant:
        cols.append(np.ones(n))
    X = np.column_stack(cols)
    yy = y[1:]
    b, *_ = np.linalg.lstsq(X, yy, rcond=None)
    u = yy - X @ b
    k = X.shape[1]
    s2 = float(u @ u) / (n - k)
    XtXi = np.linalg.inv(X.T @ X)
    sig2 = s2 * XtXi[0, 0]
    sig = math.sqrt(sig2)
    g = [float(u[j:] @ u[:n - j]) / n for j in range(q + 1)]
    lam2 = g[0] + 2 * sum((1 - j / (q + 1)) * g[j] for j in range(1, q + 1))
    rho = float(b[0])
    sn = math.sqrt(s2)
    zrho = n * (rho - 1) - (n * n * sig2 / (2 * s2)) * (lam2 - g[0])
    zt = math.sqrt(g[0] / lam2) * (rho - 1) / sig - 0.5 * (lam2 - g[0]) * n * sig / (math.sqrt(lam2) * sn)
    return zrho, zt, n, b, np.sqrt(np.diag(s2 * XtXi)), X


def cmd_pperron(s: "Session", args: str) -> None:
    from ..commands.tsspectral import _sample_series
    p = parse_standard(args)
    o = match_options(p.options, {"noconstant": 4, "trend": 3, "regress": 3, "lags": 3}) \
        if p.options.strip() else {}
    names = s.expand_varlist(p.varlist) if p.varlist.strip() else []
    if len(names) != 1:
        raise StataError(103 if names else 102, "too many variables specified" if names else "too few variables specified")
    if o.get("noconstant") and o.get("trend"):
        raise StataError(198, "options noconstant and trend may not be combined")
    y, rows = _sample_series(s, p, p.varlist.strip())
    constant = not o.get("noconstant")
    trend = bool(o.get("trend"))
    nobs = len(y) - 1
    q = int(float(s.eval(str(o["lags"])))) if isinstance(o.get("lags"), str) else int(4 * (nobs / 100) ** (2 / 9))
    zrho, zt, n, b, se, X = pp_stats(y, q, trend, constant)
    kind = "ct" if trend else ("c" if constant else "n")
    crho, ctau = _interp(_FULLER_RHO, kind, n), _interp(_FULLER_TAU, kind, n)
    out = s.output
    out.write(f"\n{'Phillips-Perron test for unit root':<51}Number of obs   = ", "text")
    out.write(f"{n:>9,}\n", "result")
    out.write(f"{'':51}Newey-West lags = ", "text")
    out.write(f"{q:>9}\n\n", "result")
    head = "Interpolated Dickey-Fuller"
    nd = 47 - len(head) - 2
    out.write(f"{'':31}{'-' * ((nd + 1) // 2)} {head} {'-' * (nd // 2)}\n", "text")
    out.write("                  Test         1% Critical       5% Critical      10% Critical\n", "text")
    out.write("               Statistic           Value             Value             Value\n", "text")
    out.write("-" * 78 + "\n", "text")
    out.write(" Z(rho)", "text")
    out.write(f"{zrho:>17.3f}{crho[0]:>18.3f}{crho[1]:>18.3f}{crho[2]:>18.3f}\n", "result")
    out.write(" Z(t)", "text")
    out.write(f"{zt:>19.3f}{ctau[0]:>18.3f}{ctau[1]:>18.3f}{ctau[2]:>18.3f}\n", "result")
    pval = None
    if constant:
        out.write("-" * 78 + "\n", "text")
        pval = mackinnon_p(zt, kind)
        out.write("MacKinnon approximate p-value for Z(t) = ", "text")
        out.write(f"{pval:.4f}\n", "result")
    if o.get("regress"):
        nm = p.varlist.strip()
        labels = [f"L1."] + (["_trend"] if trend else []) + (["_cons"] if constant else [])
        rws = [CoefRow("header", nm)]
        for i, lab in enumerate(labels):
            rws.append(CoefRow("coef", lab if lab != "L1." else "L1.", b[i], se[i]))
        out.write("\n", "text")
        coef_table(s, nm, rws, stat="t", df=n - X.shape[1], level=float(s.settings.get("level", "95")),
                   depvar_left=True)
    s.r = {"N": float(n), "lags": float(q)}
    if pval is not None:
        s.r["pval"] = pval
    s.r["Zt"] = zt
    s.r["Zrho"] = zrho
