"""prais: linear regression with AR(1) errors (Prais-Winsten and Cochrane-Orcutt)."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy import stats as st

from ..core import missing as M
from ..core.errors import StataError
from . import fvars
from .results import Estimates, post, rows_from_columns, coef_table, CoefRow
from .tsops import require_ts

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS
RHOTYPES = ("regress", "freg", "tscorr", "dw", "theil", "nagar")


def transform(X: np.ndarray, y: np.ndarray, lagged: np.ndarray, rho: float, pw: bool):
    """AR(1) transformation. `lagged[i]` says that row i has a valid previous period (the previous
    row); rows without it keep the first-observation treatment: scaled by sqrt(1-rho^2) in
    Prais-Winsten, dropped in Cochrane-Orcutt. Returns the transformed data and the kept rows."""
    Z = np.column_stack([X, y])
    T = Z.copy()
    prev = np.vstack([np.zeros((1, Z.shape[1])), Z[:-1]])
    T[lagged] = Z[lagged] - rho * prev[lagged]
    keep = np.ones(len(y), dtype=bool)
    if pw:
        T[~lagged] = Z[~lagged] * np.sqrt(1 - rho * rho)
    else:
        keep = lagged
    return T[keep, :-1], T[keep, -1], keep


def ols(X: np.ndarray, y: np.ndarray):
    XtXi = np.linalg.pinv(X.T @ X)
    b = XtXi @ X.T @ y
    return b, y - X @ b, XtXi


def estimate_rho(e: np.ndarray, lagged: np.ndarray, rhotype: str, k: int) -> float:
    """Autocorrelation of the residuals `e` (in time order) by the chosen estimator."""
    cur = e[lagged]
    prev = np.vstack([[0.0], e[:-1].reshape(-1, 1)]).ravel()[lagged]
    n = len(e)
    if rhotype == "regress":
        return float(cur @ prev / (prev @ prev))
    if rhotype == "freg":
        nxt = np.append(lagged[1:], False)          # rows whose successor has a valid lag
        c2 = e[np.flatnonzero(nxt) + 1]
        p2 = e[nxt]
        return float(p2 @ c2 / (p2 @ p2))
    tscorr = float(cur @ prev / (e @ e))
    if rhotype == "tscorr":
        return tscorr
    if rhotype == "theil":
        return tscorr * (n - k) / n
    d = e[lagged] - prev
    rdw = 1 - float(d @ d) / (e @ e) / 2
    if rhotype == "dw":
        return rdw
    return (rdw * n * n + k * k) / (n * n - k * k)      # nagar


def durbin_watson(e: np.ndarray, lagged: np.ndarray) -> float:
    d = e[lagged] - np.concatenate([[0.0], e[:-1]])[lagged]
    return float(d @ d / (e @ e))


def fit_prais(X: np.ndarray, y: np.ndarray, lagged: np.ndarray, *, pw: bool, rhotype: str,
              twostep: bool, ssesearch: bool, maxiter: int, tol: float, k_slopes: int, log):
    """Returns rho, the OLS estimates of the transformed model and the iteration count."""
    b0, e0, _ = ols(X, y)
    k = X.shape[1]
    sse = lambda r: float(np.sum(ols(*transform(X, y, lagged, r, pw)[:2])[1] ** 2))
    if ssesearch:
        grid = np.linspace(-0.999, 0.999, 400)
        r = grid[int(np.argmin([sse(g) for g in grid]))]
        lo, hi = max(r - 0.005, -0.9999), min(r + 0.005, 0.9999)
        for _ in range(60):
            m1, m2 = lo + (hi - lo) * 0.382, lo + (hi - lo) * 0.618
            lo, hi = (lo, m2) if sse(m1) < sse(m2) else (m1, hi)
        rho, it = (lo + hi) / 2, 0
    else:
        rho, it = 0.0, 0
        log(0, rho)
        e = e0
        while it < maxiter:
            new = estimate_rho(e, lagged, rhotype, k)
            it += 1
            log(it, new)
            done = abs(new - rho) < tol
            rho = new
            Xt, yt, _ = transform(X, y, lagged, rho, pw)
            _, et, _ = ols(Xt, yt)
            if twostep or done:
                break
            # residuals of the original model at the transformed estimates
            bt = ols(Xt, yt)[0]
            e = y - X @ bt
    Xt, yt, keep = transform(X, y, lagged, rho, pw)
    b, et, XtXi = ols(Xt, yt)
    return rho, b, et, XtXi, Xt, yt, keep, e0, it


def cmd_prais(s: "Session", args: str) -> None:
    from ..commands.tshelp import panels
    from .postest import replay
    from .regress import Sample, cmdline
    if not args.strip() or args.strip().startswith(","):
        replay(s, "prais", args)
        return
    smp = Sample(s, args, opts_spec={"rhotype": 3, "corc": 4, "ssesearch": 3, "twostep": 3,
                                     "hascons": 4, "nodw": 4, "iterate": 4, "tolerance": 3,
                                     "nolog": 3, "noheader": 6, "notable": 5},
                 weights=())
    o = smp.o
    if smp.vce != "ols":
        raise StataError(198, "option vce() not allowed")   # VERIFY: robust variance is not implemented
    rhotype = str(o.get("rhotype") or "regress").lower()
    if rhotype not in RHOTYPES:
        raise StataError(198, "option rhotype() incorrectly specified")
    pw = not o.get("corc")
    maxiter = int(o["iterate"]) if o.get("iterate") else 100
    tol = float(o["tolerance"]) if o.get("tolerance") else 1e-6
    info = require_ts(s.data)
    fvars.drop_collinear(smp.cols, smp.constant)
    est_cols = [c for c in smp.cols if c.values is not None]
    n = smp.n
    X = np.column_stack([c.values for c in est_cols] + ([np.ones(n)] if smp.constant else []))
    y = smp.y.astype(float)
    # time order within panel
    rows_all = np.flatnonzero(smp.mask)
    pos_in = {r: i for i, r in enumerate(rows_all)}
    order, lagged_l = [], []
    gaps = 0
    for _k, rows, tt in panels(s, smp.mask):
        for j, r in enumerate(rows):
            order.append(pos_in[r])
            ok = j > 0 and abs(tt[j] - tt[j - 1] - info.delta) < 1e-9
            if j > 0 and not ok:
                gaps += 1
            lagged_l.append(ok)
    order = np.array(order)
    if len(order) != n:
        raise StataError(2000, "no observations")
    X, y, lagged = X[order], y[order], np.array(lagged_l)
    k = X.shape[1]
    k_slopes = len(est_cols)
    out = s.output

    def log(it, r):
        if not o.get("nolog"):
            out.write(f"Iteration {it}:  rho = {r:.4f}\n", "text")
    rho, b, et, XtXi, Xt, yt, keep, e0, it = fit_prais(
        X, y, lagged, pw=pw, rhotype=rhotype, twostep=bool(o.get("twostep")),
        ssesearch=bool(o.get("ssesearch")), maxiter=maxiter, tol=tol, k_slopes=k_slopes, log=log)
    nn = len(yt)
    rss = float(et @ et)
    tss = float(((yt - yt.mean()) ** 2).sum()) if smp.constant else float(yt @ yt)
    mss = tss - rss
    df_r = nn - k
    df_m = k_slopes
    s2 = rss / df_r
    r2 = mss / tss if tss > 0 else SYS
    r2_a = 1 - (1 - r2) * (nn - (1 if smp.constant else 0)) / df_r if tss > 0 else SYS
    F = (mss / df_m) / s2 if df_m > 0 else 0.0
    rmse = float(np.sqrt(s2))
    ll = -0.5 * nn * (np.log(2 * np.pi) + np.log(rss / nn) + 1) + (0.5 * np.log(1 - rho * rho) if pw else 0.0)
    Vb = s2 * XtXi
    dw0 = durbin_watson(y - X @ ols(X, y)[0], lagged)
    dw = durbin_watson(et, _lag_after_keep(lagged, keep))
    names = [c.name for c in smp.cols] + (["_cons"] if smp.constant else [])
    pos = [i for i, c in enumerate(smp.cols) if c.values is not None] + ([len(smp.cols)] if smp.constant else [])
    kk = len(names)
    bf, Vf = np.zeros(kk), np.zeros((kk, kk))
    bf[pos] = b
    Vf[np.ix_(pos, pos)] = Vb
    bmap = {c.name: bf[i] for i, c in enumerate(smp.cols)}
    semap = {c.name: float(np.sqrt(max(Vf[i, i], 0))) for i, c in enumerate(smp.cols)}
    if smp.constant:
        bmap["_cons"], semap["_cons"] = bf[-1], float(np.sqrt(max(Vf[-1, -1], 0)))
    est = Estimates("prais", smp.depname, names, bf, Vf, nn, "t", float(df_r))
    est.rows = rows_from_columns(smp.cols, bmap, semap, s.data, constant=smp.constant)
    est.rows.append(CoefRow("value", "rho", rho))
    est.cols, est.constant, est.terms, est.depcomp = smp.cols, smp.constant, smp.terms, smp.depcomp
    est.level = smp.level
    est.vcetype = ""
    method = "twostep" if o.get("twostep") else ("ssesearch" if o.get("ssesearch") else "iterated")
    method_t = {"iterated": "iterated", "twostep": "twostep", "ssesearch": "SSE search"}[method]
    title = ("Prais-Winsten AR(1) regression" if pw else "Cochrane-Orcutt AR(1) regression") + \
        f" -- {method_t} estimates"
    est.extra = {"N": nn, "df_m": df_m, "anova_df_m": df_m, "df_r": df_r, "F": F, "r2": r2,
                 "r2_a": r2_a, "rmse": rmse, "mss": mss, "rss": rss, "tss": tss, "rho": rho,
                 "dw": dw, "dw0": dw0, "title": title, "nodw": bool(o.get("nodw"))}
    est.display = display_prais
    scalars = [("N", nn), ("N_gaps", gaps), ("mss", mss), ("df_m", df_m), ("rss", rss), ("df_r", df_r),
               ("r2", r2), ("r2_a", r2_a), ("F", F), ("rmse", rmse), ("ll", ll), ("rho", rho),
               ("dw", dw), ("dw_0", dw0), ("rank", k), ("tol", tol), ("max_ic", maxiter), ("ic", it)]
    macros = [("cmd", "prais"), ("cmdline", cmdline("prais", args)), ("depvar", smp.depname),
              ("title", title), ("method", method), ("rhotype", rhotype),
              ("tranmeth", "prais" if pw else "corc"), ("properties", "b V"), ("estat_cmd", "regress_estat")]
    post(s, est, scalars, macros, smp.mask)
    display_prais(s, est, header=not o.get("noheader"), table=not o.get("notable"))


def _lag_after_keep(lagged: np.ndarray, keep: np.ndarray) -> np.ndarray:
    """Lag indicator for the rows of the transformed sample (the first row never has one)."""
    lag = lagged[keep].copy()
    if len(lag):
        lag[0] = False
    return lag


def display_prais(s: "Session", est: Estimates, *, header: bool = True, table: bool = True,
                  level: float | None = None, **_kw) -> None:
    from .regress import _anova
    out = s.output
    x = est.extra
    if header:
        out.write("\n", "text")
        out.write(x["title"] + "\n", "text")
        out.write("\n", "text")
        _anova(s, est)
        out.write("\n", "text")
    if table:
        coef_table(s, est.depvar, est.rows, stat="t", df=est.df_r, level=level or est.level)
        if not x["nodw"]:
            out.write(f"Durbin-Watson statistic (original)    {x['dw0']:.6f}\n", "text")
            out.write(f"Durbin-Watson statistic (transformed) {x['dw']:.6f}\n", "text")
