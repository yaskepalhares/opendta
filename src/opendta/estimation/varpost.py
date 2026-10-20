"""Postestimation of var: varlmar, varnorm, varstable, varwle and fcast ([TS] varlmar, [TS] varnorm,
[TS] varstable, [TS] varwle, [TS] fcast). The statistics are plain functions of arrays; the
commands read the stored VAR, print and store r()."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.stats import chi2_tail
from ..lang.syntax import match_options, parse_standard
from .results import current

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


# ---------------------------------------------------------------------------
# statistics
# ---------------------------------------------------------------------------

def lm_autocorr(U: np.ndarray, X: np.ndarray, lag: int) -> tuple[float, int]:
    """Johansen LM statistic for autocorrelation of order `lag`: the VAR is augmented with the K
    residual series lagged `lag` times (zeros for the missing initial values) and the statistic
    is (T - d - 0.5) ln(|Sigma| / |Sigma_s|). Returns (statistic, df = K^2)."""
    T, K = U.shape
    Ul = np.zeros_like(U)
    Ul[lag:] = U[:-lag]
    Xa = np.column_stack([X, Ul])
    B, *_ = np.linalg.lstsq(Xa, U, rcond=None)
    Ua = U - Xa @ B
    d = X.shape[1] + K          # VERIFY: Stata's d in the statistic
    s0 = np.linalg.slogdet(U.T @ U / T)[1]
    s1 = np.linalg.slogdet(Ua.T @ Ua / T)[1]
    return (T - d - 0.5) * (s0 - s1), K * K


def normality(U: np.ndarray, Sigma: np.ndarray, P: np.ndarray | None = None) -> dict:
    """Skewness, kurtosis and Jarque-Bera tests of the residuals orthogonalised with P (default:
    Cholesky factor of Sigma). Returns the coefficients and the statistics by equation."""
    T, K = U.shape
    if P is None:
        P = np.linalg.cholesky(Sigma)
    W = np.linalg.solve(P, U.T).T
    b1 = (W ** 3).mean(axis=0)
    b2 = (W ** 4).mean(axis=0)
    sk = T * b1 ** 2 / 6
    ku = T * (b2 - 3) ** 2 / 24
    return {"skew": b1, "kurt": b2, "chi_skew": sk, "chi_kurt": ku, "chi_jb": sk + ku,
            "all_skew": float(T / 6 * (b1 @ b1)), "all_kurt": float(T / 24 * ((b2 - 3) @ (b2 - 3)))}


def companion(Alist: list[np.ndarray]) -> np.ndarray:
    K = Alist[0].shape[0]
    p = len(Alist)
    A = np.zeros((K * p, K * p))
    A[:K] = np.hstack(Alist)
    if p > 1:
        A[K:, :-K] = np.eye(K * (p - 1))
    return A


def lag_matrices(B: np.ndarray, xnames: list, deps: list[str], p: int) -> list[np.ndarray]:
    """A_1..A_p (K x K) from the coefficient matrix B (rows: regressors, columns: equations)."""
    K = len(deps)
    out = [np.zeros((K, K)) for _ in range(p)]
    for i, (vn, l) in enumerate(xnames):
        if l > 0:
            out[l - 1][:, deps.index(vn)] = B[i]
    return out


def wald_exclusion(B, Sig, XtXi, idx: list[int], j: int | None):
    """Wald chi2 and df that the coefficients in rows `idx` are zero in equation j (None: in all)."""
    if j is not None:
        b, V = B[idx, j], Sig[j, j] * XtXi[np.ix_(idx, idx)]
    else:
        b, V = B[idx, :].T.ravel(), np.kron(Sig, XtXi[np.ix_(idx, idx)])
    return float(b @ np.linalg.solve(V, b)), len(b)


def forecast_path(y_hist: np.ndarray, B: np.ndarray, xnames: list, deps: list[str], constant: bool,
                  steps: int) -> np.ndarray:
    """Dynamic forecasts: y_hist holds the observed values before the first forecast (rows in time
    order, at least max lag rows); each step uses the previous forecasts."""
    K = len(deps)
    hist = [r for r in y_hist]
    out = []
    for _ in range(steps):
        x = []
        for vn, l in xnames:
            if vn == "_cons":
                x.append(1.0)
            else:
                x.append(hist[-l][deps.index(vn)])
        yhat = np.array(x) @ B
        out.append(yhat)
        hist.append(yhat)
    return np.array(out).reshape(steps, K)


def forecast_mse(B: np.ndarray, Sig: np.ndarray, XtXi: np.ndarray, X: np.ndarray, T: int,
                 A: list[np.ndarray], constant: bool, steps: int) -> np.ndarray:
    """Asymptotic MSE matrices of the h-step forecasts (Lutkepohl 2005, section 3.5.1):
    Sigma_y(h) + Omega(h)/T. Returns an array (steps, K, K)."""
    K = Sig.shape[0]
    p = len(A)
    Phi = [np.eye(K)]
    for i in range(1, steps):
        Phi.append(sum(Phi[i - j] @ A[j - 1] for j in range(1, min(i, p) + 1)))
    # companion matrix of (1, y_t, ..., y_{t-p+1}); the regressors of the sample are Z_{t-1}
    n = K * p + (1 if constant else 0)
    Bc = np.zeros((n, n))
    o = 1 if constant else 0
    if constant:
        Bc[0, 0] = 1.0
        Bc[o:o + K, 0] = B[-1]
    Bc[o:o + K, o:] = np.hstack(A)
    if p > 1:
        Bc[o + K:, o:o + K * (p - 1)] = np.eye(K * (p - 1))
    # the regressor order of X may differ from Z: reorder to (1, lag1 block, lag2 block ...)
    Gam = X.T @ X / T
    Gi = np.linalg.inv(Gam)
    out = np.zeros((steps, K, K))
    acc = np.zeros((K, K))
    for h in range(1, steps + 1):
        acc = acc + Phi[h - 1] @ Sig @ Phi[h - 1].T
        Om = np.zeros((K, K))
        for i in range(h):
            for j in range(h):
                tr = np.trace(np.linalg.matrix_power(Bc.T, h - 1 - i) @ Gi
                              @ np.linalg.matrix_power(Bc, h - 1 - j) @ Gam)
                Om += tr * Phi[i] @ Sig @ Phi[j].T
        out[h - 1] = acc + Om / T
    return out


# ---------------------------------------------------------------------------
# commands
# ---------------------------------------------------------------------------

def _var(s: "Session", name: str):
    est = current(s)
    if est.cmd not in ("var", "svar"):
        raise StataError(301, "last estimates not found")
    return est, est.var


def _opts(args: str, spec: dict) -> dict:
    p = parse_standard(args)
    if p.varlist.strip():
        raise StataError(198, "invalid syntax")
    return match_options(p.options, spec) if p.options.strip() else {}


def cmd_varlmar(s: "Session", args: str) -> None:
    from ..commands.matrix import Matrix
    est, V = _var(s, "varlmar")
    o = _opts(args, {"mlag": 2, "estimates": 3, "separator": 3})
    m = int(o["mlag"]) if o.get("mlag") and o["mlag"] is not True else 2
    if m < 1:
        raise StataError(198, "mlag() must be greater than 0")   # VERIFY
    out = s.output
    rows = []
    for lag in range(1, m + 1):
        chi, df = lm_autocorr(V["U"], V["X"], lag)
        rows.append((lag, chi, df, float(chi2_tail(chi, df))))
    sep = int(o["separator"]) if o.get("separator") and o["separator"] is not True else 0
    out.write("\nLagrange-multiplier test\n\n", "text")     # VERIFY layout
    out.write(" lag       chi2   df       Prob > chi2\n", "text")
    out.write("-" * 38 + "\n", "text")
    for i, (lag, chi, df, pv) in enumerate(rows):
        out.write(f"{lag:>4}{chi:>11.4f}{df:>5}{pv:>16.5f}\n", "result")
        if sep and (i + 1) % sep == 0 and i + 1 < len(rows):
            out.write("-" * 38 + "\n", "text")
    out.write("-" * 38 + "\n", "text")
    out.write("H0: no autocorrelation at lag order\n", "text")
    s.r = {"lm": Matrix(np.array([[r[1], r[2], r[3]] for r in rows]),
                        [str(r[0]) for r in rows], ["chi2", "df", "p"])}


def cmd_varnorm(s: "Session", args: str) -> None:
    from ..commands.matrix import Matrix
    est, V = _var(s, "varnorm")
    o = _opts(args, {"jbera": 1, "skewness": 1, "kurtosis": 1, "estimates": 3, "separator": 3, "cholesky": 4})
    if o.get("cholesky"):
        raise StataError(198, "option cholesky not allowed")    # only for svar
    which = [k for k in ("jbera", "skewness", "kurtosis") if o.get(k)] or ["jbera", "skewness", "kurtosis"]
    sep = int(o["separator"]) if o.get("separator") and o["separator"] is not True else 0
    deps = V["deps"]
    K = len(deps)
    r = normality(V["U"], V["Sig"])
    out = s.output
    tabs = {
        "jbera": ("Jarque-Bera test", None, r["chi_jb"], 2, r["all_skew"] + r["all_kurt"], 2 * K),
        "skewness": ("Skewness test", r["skew"], r["chi_skew"], 1, r["all_skew"], K),
        "kurtosis": ("Kurtosis test", r["kurt"], r["chi_kurt"], 1, r["all_kurt"], K),
    }
    res = {}
    for key in ("jbera", "skewness", "kurtosis"):
        if key not in which:
            continue
        title, coef, chi, df, chi_all, df_all = tabs[key]
        out.write(f"\n  {title}\n\n", "text")                      # VERIFY layout
        head = f"{'Equation':>18}" + (f"{'Skewness' if key == 'skewness' else 'Kurtosis':>11}" if coef is not None else "")
        out.write(head + f"{'chi2':>11}{'df':>5}{'Prob > chi2':>14}\n", "text")
        width = len(head) + 30
        out.write("-" * width + "\n", "text")
        mat = []
        for i, dv in enumerate(deps):
            pv = float(chi2_tail(chi[i], df))
            mat.append([chi[i], df, pv])
            out.write(f"{dv[:18]:>18}" + (f"{coef[i]:>11.5f}" if coef is not None else "") +
                      f"{chi[i]:>11.3f}{df:>5}{pv:>14.5f}\n", "result")
            if sep and (i + 1) % sep == 0:
                out.write("-" * width + "\n", "text")
        pv = float(chi2_tail(chi_all, df_all))
        mat.append([chi_all, df_all, pv])
        out.write(f"{'ALL':>18}" + (" " * 11 if coef is not None else "") +
                  f"{chi_all:>11.3f}{df_all:>5}{pv:>14.5f}\n", "result")
        out.write("-" * width + "\n", "text")
        res[key] = Matrix(np.array(mat), deps + ["ALL"], ["chi2", "df", "p"])
    if V["dfk"]:
        out.write("dfk estimator used in computations\n", "text")
    s.r = {k: res[k] for k in ("jbera", "kurtosis", "skewness") if k in res}
    if "jbera" in s.r:
        s.r["jb"] = s.r.pop("jbera")
    if V["dfk"]:
        s.r = {"dfk": "dfk", **s.r}
    s.r = dict(reversed(list(s.r.items())))


def cmd_varstable(s: "Session", args: str) -> None:
    from ..commands.matrix import Matrix
    est, V = _var(s, "varstable")
    o = _opts(args, {"estimates": 3, "amat": 2, "graph": 1, "dlabel": 2, "modlabel": 3, "nogrid": 3,
                     "pgrid": 2, "rlopts": 2, "addplot": 3})
    p = max(V["lags"])
    A = companion(lag_matrices(V["B"], V["xnames"], V["deps"], p))
    ev = np.linalg.eigvals(A)
    order = np.argsort(-np.abs(ev), kind="stable")
    ev = ev[order]
    mod = np.abs(ev)
    out = s.output
    from .tsmodels import _g
    out.write("\nEigenvalue stability condition\n\n", "text")
    out.write(f"{'Eigenvalue':<34}Modulus\n", "text")
    for z, m in zip(ev, mod):
        if abs(z.imag) < 1e-12:
            lab = _g(z.real)
        else:
            lab = f"{_g(z.real)} {'+' if z.imag > 0 else '-'} {_g(abs(z.imag))}i"
        out.write(f"{lab:<34}{f'{m:.6f}'.replace('0.', '.', 1) if m < 1 else f'{m:.6f}'}\n", "result")
    if np.all(mod < 1):
        out.write("\nAll the eigenvalues lie inside the unit circle.\nVAR satisfies stability condition.\n", "text")
    else:
        out.write("\nAt least one eigenvalue is outside the unit circle.\nVAR does not satisfy stability condition.\n", "text")
    s.r = {"Modulus": Matrix(mod.reshape(-1, 1)), "Im": Matrix(ev.imag.reshape(-1, 1)),
           "Re": Matrix(ev.real.reshape(-1, 1))}
    if o.get("amat") and o["amat"] is not True:
        from ..commands.matrix import _store
        _store(s)[str(o["amat"])] = Matrix(A)
    if o.get("graph"):
        _stability_graph(s, ev, mod, o, args)


def _stability_graph(s, ev, mod, o, args) -> None:
    from ..commands.graphcmd import finish
    from ..graphics.model import Graph, Layer, Panel
    from ..graphics.opts import graph_options, split_layer_graph
    p = parse_standard(args)
    _l, graw = split_layer_graph(p.options, {"estimates": 3, "amat": 2, "graph": 1, "dlabel": 2,
                                             "modlabel": 3, "nogrid": 3, "pgrid": 2, "rlopts": 2,
                                             "addplot": 3})
    th = np.linspace(0, 2 * math.pi, 181)
    layers = []
    if not o.get("nogrid"):
        for r in np.arange(0.1, 1.0, 0.1):
            layers.append(Layer("line", list(r * np.cos(th)), list(r * np.sin(th)), {},
                                {"lcolor": "gs14"}, ""))
    layers.append(Layer("line", list(np.cos(th)), list(np.sin(th)), {}, {"lcolor": "navy"}, ""))
    layers.append(Layer("scatter", list(ev.real), list(ev.imag), {}, {}, ""))
    opts = graph_options(graw)
    opts["_xtitle"] = "Real"
    opts["_ytitle"] = "Imaginary"
    opts.setdefault("title", {"lines": ["Roots of the companion matrix"]})
    finish(s, Graph(family="twoway", opts=opts, panels=[Panel("", layers)]), graw)


def cmd_varwle(s: "Session", args: str) -> None:
    from ..commands.matrix import Matrix
    from scipy.stats import f as _f
    est, V = _var(s, "varwle")
    o = _opts(args, {"estimates": 3, "separator": 3})
    sep = int(o["separator"]) if o.get("separator") and o["separator"] is not True else 0
    deps, xnames, B, Sig, XtXi = V["deps"], V["xnames"], V["B"], V["Sig"], V["XtXi"]
    K, small = len(deps), V["small"]
    dfr = V["T"] - len(xnames)
    out = s.output
    lags = sorted({l for _v, l in xnames if l > 0})
    groups = [(dv, j) for j, dv in enumerate(deps)] + [("All", None)]
    last = None
    for name, j in groups:
        out.write(f"\n   Equation: {name}\n\n", "text")                   # VERIFY layout
        hdr = f"{'lag':>5}{'F':>12}{'df':>5}{'df_r':>6}{'Prob > F':>12}" if small else \
            f"{'lag':>5}{'chi2':>12}{'df':>5}{'Prob > chi2':>14}"
        out.write(hdr + "\n", "text")
        out.write("-" * len(hdr) + "\n", "text")
        rows = []
        for i, l in enumerate(lags):
            idx = [k for k, (_v, ll) in enumerate(xnames) if ll == l]
            chi, df = wald_exclusion(B, Sig, XtXi, idx, j)
            if small:
                F = chi / df
                pv = float(_f.sf(F, df, dfr))
                rows.append((l, F, df, dfr, pv))
                out.write(f"{l:>5}{F:>12.5f}{df:>5}{dfr:>6}{pv:>12.4f}\n", "result")
            else:
                pv = float(chi2_tail(chi, df))
                rows.append((l, chi, df, pv))
                out.write(f"{l:>5}{chi:>12.5f}{df:>5}{pv:>14.3f}\n", "result")
            if sep and (i + 1) % sep == 0 and i + 1 < len(lags):
                out.write("-" * len(hdr) + "\n", "text")
        out.write("-" * len(hdr) + "\n", "text")
        last = rows
    # r(): the last table ("All"), as a matrix lags x columns
    names = [str(r[0]) for r in last]
    if small:
        s.r = {"p": Matrix(np.array([[r[4]] for r in last]), names, ["p"]),
               "df": Matrix(np.array([[r[3]] for r in last]), names, ["df"]),
               "df_r": Matrix(np.array([[r[2]] for r in last]), names, ["df_r"]),
               "F": Matrix(np.array([[r[1]] for r in last]), names, ["F"])}
    else:
        s.r = {"p": Matrix(np.array([[r[3]] for r in last]), names, ["p"]),
               "df": Matrix(np.array([[r[2]] for r in last]), names, ["df"]),
               "chi2": Matrix(np.array([[r[1]] for r in last]), names, ["chi2"])}


# ---------------------------------------------------------------------------
# fcast
# ---------------------------------------------------------------------------

def cmd_fcast(s: "Session", args: str) -> None:
    sub, _, rest = args.strip().partition(" ")
    sub = sub.lower().rstrip(",")
    if sub == "compute" or (len(sub) >= 1 and "compute".startswith(sub) and sub[0] == "c"):
        _fcast_compute(s, rest)
    elif sub and "graph".startswith(sub):
        _fcast_graph(s, rest)
    else:
        raise StataError(198, "invalid syntax")


def _fcast_compute(s: "Session", rest: str) -> None:
    from ..commands.tshelp import append_periods, new_variable
    from .tsops import require_ts
    from scipy import stats as st
    p = parse_standard(rest)
    prefix = p.varlist.strip()
    if not prefix or " " in prefix:
        raise StataError(198, "invalid syntax")
    o = match_options(p.options, {"step": 2, "dynamic": 3, "estimates": 3, "replace": 3, "nose": 2,
                                  "bs": 2, "bsp": 3, "bscentile": 3, "reps": 3, "nodots": 3,
                                  "saving": 3, "level": 1}) if p.options.strip() else {}
    if o.get("bs") or o.get("bsp") or o.get("bscentile"):
        raise StataError(198, "options bs, bsp and bscentile are not implemented in OpenDTA")   # OpenDTA
    est, V = _var(s, "fcast")
    if V["exog"]:
        raise StataError(198, "fcast after var with exog() is not implemented in OpenDTA")   # OpenDTA
    deps, xnames, B = V["deps"], V["xnames"], V["B"]
    K = len(deps)
    lags = V["lags"]
    pmax = max(lags)
    if lags != list(range(1, pmax + 1)):
        raise StataError(198, "fcast requires consecutive lags")   # OpenDTA
    steps = int(o["step"]) if o.get("step") and o["step"] is not True else 1
    level = float(o["level"]) if o.get("level") and o["level"] is not True else float(s.settings.get("level", "95"))
    info = require_ts(s.data)
    if info.panel:
        raise StataError(198, "panel data not allowed")   # VERIFY
    t0, t1, dl = V["t0"], V["t1"], info.delta
    start = t1 + dl
    if o.get("dynamic") and o["dynamic"] is not True:
        start = float(s.eval(str(o["dynamic"])))
        if start < t0 + pmax * dl or start > t1 + dl:
            raise StataError(198, "dynamic() must be within the estimation sample or one period after")   # VERIFY
    names = {}
    for dv in deps:
        for suf in ("", "_LB", "_UB", "_SE"):
            if suf and o.get("nose"):
                continue
            nm = f"{prefix}{dv}{suf}"
            if s.data.has(nm) and not o.get("replace"):
                raise StataError(110, f"{nm} already defined")
            names[(dv, suf)] = nm
    times = [start + h * dl for h in range(steps)]
    ds = s.data
    tv = ds.get(info.tvar).data.astype(float)
    have = {float(v) for v in tv[tv < SYS]}
    new = [(None, tt) for tt in times if tt not in have]
    if new:
        append_periods(s, new)
        ds = s.data
        tv = ds.get(info.tvar).data.astype(float)
    pos = {float(v): i for i, v in enumerate(tv) if v < SYS}
    hist = []
    for k in range(pmax, 0, -1):
        row = pos.get(start - k * dl)
        if row is None:
            raise StataError(198, "observations needed before the forecast period are missing")   # VERIFY
        hist.append([float(ds.get(dv).data[row]) for dv in deps])
    hist = np.array(hist)
    if np.any(hist >= SYS):
        raise StataError(198, "observations needed before the forecast period are missing")   # VERIFY
    fc = forecast_path(hist, B, xnames, deps, V["constant"], steps)
    cols = {}
    if not o.get("nose"):
        A = lag_matrices(B, xnames, deps, pmax)
        perm = ([xnames.index(("_cons", 0))] if V["constant"] else []) + \
            [xnames.index((dv, l)) for l in range(1, pmax + 1) for dv in deps]
        mse = forecast_mse(B[perm], V["Sig"], V["XtXi"], V["X"][:, perm], V["T"], A,
                           V["constant"], steps)
        z = float(st.norm.isf((1 - level / 100) / 2))
        se = np.sqrt(np.array([np.diag(m) for m in mse]))
    n = ds.nobs
    for j, dv in enumerate(deps):
        def put(suf, vals):
            col = np.full(n, np.nan)
            for tt, v in zip(times, vals):
                col[pos[tt]] = v
            new_variable(s, names[(dv, suf)], col, vtype="float", replace=True)
        put("", fc[:, j])
        if not o.get("nose"):
            put("_LB", fc[:, j] - z * se[:, j])
            put("_UB", fc[:, j] + z * se[:, j])
            put("_SE", se[:, j])
    s.notify_state()


def _fcast_graph(s: "Session", rest: str) -> None:
    from ..commands.graphcmd import finish
    from ..commands.tshelp import to_nan
    from ..graphics.model import Graph, Layer, Panel
    from ..graphics.opts import graph_options, split_layer_graph
    from .tsops import require_ts
    p = parse_standard(rest)
    lopts, graw = split_layer_graph(p.options, {"noci": 4, "observed": 3, "level": 1})
    est, V = _var(s, "fcast graph")
    info = require_ts(s.data)
    ds = s.data
    tv = to_nan(ds.get(info.tvar).data)
    names = s.expand_varlist(p.varlist)
    panels_ = []
    for nm in names:
        dep = next((dv for dv in V["deps"] if nm.endswith(dv)), None)
        if dep is None:
            raise StataError(198, f"{nm} is not a forecast variable")   # VERIFY
        f = to_nan(ds.get(nm).data)
        ok = np.isfinite(f) & np.isfinite(tv)
        order = np.argsort(tv[ok])
        xs = tv[ok][order]
        layers = []
        lb, ub = f"{nm}_LB", f"{nm}_UB"
        if not lopts.get("noci") and ds.has(lb) and ds.has(ub):
            lo = to_nan(ds.get(lb).data)[ok][order]
            hi = to_nan(ds.get(ub).data)[ok][order]
            layers.append(Layer("rarea", xs.tolist(), lo.tolist(), {"y2": hi.tolist()},
                                {"_pal": 0, "_ci": True}, "95% CI"))
        layers.append(Layer("line", xs.tolist(), f[ok][order].tolist(), {}, {"_pal": 0}, "forecast"))
        if lopts.get("observed"):
            y = to_nan(ds.get(dep).data)
            ob = np.isfinite(y) & np.isfinite(tv) & (tv >= xs.min() if xs.size else False)
            oo = np.argsort(tv[ob])
            layers.append(Layer("line", tv[ob][oo].tolist(), y[ob][oo].tolist(), {}, {"_pal": 1}, "observed"))
        panels_.append(Panel(nm, layers))
    opts = graph_options(graw)
    opts["_xtitle"] = info.tvar
    finish(s, Graph(family="twoway", opts=opts, panels=panels_), graw)
