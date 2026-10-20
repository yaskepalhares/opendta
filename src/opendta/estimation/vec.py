"""vecrank and vec: Johansen's cointegration rank test and the vector error-correction model
([TS] vecrank, [TS] vec). The statistics are plain functions of arrays; the commands read the
data, print and store e()."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.stats import chi2_tail
from ..lang.syntax import match_options, parse_standard
from .results import CoefRow, Estimates, coef_table, post
from .tsops import _fmt_time, require_ts

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS
TRENDS = ("constant", "rconstant", "trend", "rtrend", "none")

# Osterwald-Lenum (1992) critical values by number of common trends K - r = 1..10. VERIFY: only the
# unrestricted-constant case (Stata's default) has been checked against published output.
_CV = {
    ("constant", "trace", 5): [3.76, 15.41, 29.68, 47.21, 68.52, 94.15, 124.24, 156.00, 192.89, 233.13],
    ("constant", "trace", 1): [6.65, 20.04, 35.65, 54.46, 76.07, 103.18, 133.57, 168.36, 204.95, 247.18],
    ("constant", "max", 5): [3.76, 14.07, 20.97, 27.07, 33.46, 39.37, 45.28, 51.42, 57.12, 62.81],
    ("constant", "max", 1): [6.65, 18.63, 25.52, 32.24, 38.77, 45.10, 51.57, 57.69, 62.80, 69.09],
    ("rconstant", "trace", 5): [9.24, 19.96, 34.91, 53.12, 76.07, 102.14, 131.70, 165.58, 202.92, 244.15],
    ("rconstant", "trace", 1): [12.97, 24.60, 41.07, 60.16, 84.45, 111.01, 143.09, 177.20, 215.74, 257.68],
    ("rconstant", "max", 5): [9.24, 15.67, 22.00, 28.14, 34.40, 40.30, 46.45, 52.00, 57.42, 63.57],
    ("rconstant", "max", 1): [12.97, 20.20, 26.81, 33.24, 39.79, 46.82, 51.91, 57.95, 63.71, 69.94],
}


def critical_value(trend: str, kind: str, level: int, ntrends: int) -> float:
    """Critical value with `ntrends` = K - r common trends; missing when not tabulated."""
    tab = _CV.get((trend, kind, level))
    if tab is None or not 1 <= ntrends <= len(tab):
        return SYS
    return tab[ntrends - 1]


# ---------------------------------------------------------------------------
# estimation
# ---------------------------------------------------------------------------

def _moments(Y: np.ndarray, p: int, trend: str):
    """Residuals R0 (differences) and R1 (lagged levels plus restricted terms) after partialling out
    the lagged differences and the unrestricted deterministic terms; T = n - p observations."""
    n, K = Y.shape
    dY = np.diff(Y, axis=0)                          # row i = Y[i+1] - Y[i]
    T = n - p
    Z0 = dY[p - 1:]                                  # Delta y_t, t = p..n-1
    Z1 = Y[p - 1:n - 1]
    tt = np.arange(p + 1, n + 1, dtype=float).reshape(-1, 1)    # VERIFY: trend counted from the first observation
    one = np.ones((T, 1))
    if trend in ("rconstant",):
        Z1 = np.hstack([Z1, one])
    elif trend == "rtrend":
        Z1 = np.hstack([Z1, tt])
    lagged = [dY[p - 1 - i:n - 1 - i] for i in range(1, p)]
    unr = {"constant": [one], "rconstant": [], "trend": [one, tt], "rtrend": [one], "none": []}[trend]
    Z2 = np.hstack(lagged + unr) if (lagged or unr) else np.zeros((T, 0))
    if Z2.shape[1]:
        P = Z2 @ np.linalg.pinv(Z2)
        R0, R1 = Z0 - P @ Z0, Z1 - P @ Z1
    else:
        R0, R1 = Z0, Z1
    return Z0, Z1, Z2, R0, R1, T, K


def johansen(Y: np.ndarray, p: int, trend: str) -> dict:
    """Eigenvalues (descending) and eigenvectors v (v' S11 v = I) of the Johansen problem."""
    Z0, Z1, Z2, R0, R1, T, K = _moments(Y, p, trend)
    S00, S11, S01 = R0.T @ R0 / T, R1.T @ R1 / T, R0.T @ R1 / T
    L = np.linalg.cholesky(S11)
    Li = np.linalg.inv(L)
    Mx = Li @ S01.T @ np.linalg.solve(S00, S01) @ Li.T
    lam, U = np.linalg.eigh((Mx + Mx.T) / 2)
    o = np.argsort(-lam)[:K]
    lam, V = lam[o], (Li.T @ U)[:, o]
    return {"lam": lam, "V": V, "S00": S00, "S11": S11, "S01": S01, "T": T, "K": K, "Z0": Z0,
            "Z1": Z1, "Z2": Z2, "R0": R0, "R1": R1}


def loglik(S00: np.ndarray, lam: np.ndarray, r: int, T: int, K: int) -> float:
    return -0.5 * T * (K * (math.log(2 * math.pi) + 1) + np.linalg.slogdet(S00)[1]
                       + float(np.sum(np.log1p(-lam[:r]))))


def fit_vecm(Y: np.ndarray, p: int, trend: str, r: int) -> dict:
    """Johansen estimates for rank r: beta (normalized), alpha, the coefficients of the short-run
    equations given beta (rows: error-correction terms, lagged differences, unrestricted terms) and
    the ML covariance of the residuals."""
    J = johansen(Y, p, trend)
    T = J["T"]
    beta = normalize_beta(J["V"], r)
    alpha = J["S01"] @ beta @ np.linalg.inv(beta.T @ J["S11"] @ beta)
    W = np.hstack([J["Z1"] @ beta, J["Z2"]])
    coef = np.linalg.lstsq(W, J["Z0"], rcond=None)[0]
    U = J["Z0"] - W @ coef
    return {"J": J, "T": T, "beta": beta, "alpha": alpha, "W": W, "coef": coef, "U": U,
            "Om": U.T @ U / T}


def var_matrices(coef: np.ndarray, alpha: np.ndarray, beta: np.ndarray, K: int, p: int) -> list[np.ndarray]:
    """A_1..A_p of the levels VAR implied by the VECM: A_1 = I + Pi + Gamma_1,
    A_i = Gamma_i - Gamma_{i-1}, A_p = -Gamma_{p-1}."""
    r = alpha.shape[1]
    Pi = alpha @ beta[:K].T
    G = [coef[r + i * K:r + (i + 1) * K].T for i in range(p - 1)]       # Gamma_i[j, k]
    A = [np.eye(K) + Pi + (G[0] if p > 1 else 0)]
    for i in range(1, p - 1):
        A.append(G[i] - G[i - 1])
    if p > 1:
        A.append(-G[p - 2])
    return A


def n_parms(K: int, p: int, trend: str, r: int) -> int:
    n_unr = {"constant": 1, "rconstant": 0, "trend": 2, "rtrend": 1, "none": 0}[trend]
    m1 = {"rconstant": 1, "rtrend": 1}.get(trend, 0)
    return K * K * (p - 1) + K * n_unr + (K + m1) * r + K * r - r * r


def info_criteria(ll: float, T: int, k: int) -> tuple[float, float, float]:
    return (-2 * ll / T + 2 * k / T, -2 * ll / T + math.log(T) * k / T,
            -2 * ll / T + 2 * math.log(math.log(T)) * k / T)         # AIC, SBIC, HQIC


def normalize_beta(V: np.ndarray, r: int) -> np.ndarray:
    """Johansen normalization beta = (I_r, beta_breve')'."""
    Vr = V[:, :r]
    return Vr @ np.linalg.inv(Vr[:r])


# ---------------------------------------------------------------------------
# commands
# ---------------------------------------------------------------------------

def _data(s: "Session", p, deps_text: str, lags: int, trend: str):
    from ..commands._util import touse
    from .tsmodels import _check_gaps
    deps = s.expand_varlist(deps_text)
    if len(deps) < 2:
        raise StataError(102, "too few variables specified")
    ds = s.data
    mask = touse(s, p)
    cols = [ds.get(v).data.astype(float) for v in deps]
    for c in cols:
        mask &= c < SYS
    t0, t1 = _check_gaps(s, mask)
    info = require_ts(ds)
    tv = ds.get(info.tvar).data.astype(float)
    rows = np.flatnonzero(mask)
    rows = rows[np.argsort(tv[rows], kind="stable")]
    Y = np.column_stack([c[rows] for c in cols])
    if Y.shape[0] <= lags + 1 + len(deps) * lags:
        raise StataError(2001, "insufficient observations")
    return deps, Y, rows, tv[rows], mask


def _trend_opt(o: dict) -> str:
    t = o.get("trend")
    if not t or t is True:
        return "constant"
    t = str(t).lower()
    for full in TRENDS:
        if full.startswith(t) and len(t) >= 1:
            if t == "r":
                break
            return full
    raise StataError(198, "option trend() incorrectly specified")


def cmd_vecrank(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"lags": 2, "trend": 2, "sindicators": 2, "noreduce": 3,
                                  "notrace": 4, "max": 3, "ic": 2, "level99": 7, "levela": 6}) \
        if p.options.strip() else {}
    if o.get("sindicators"):
        raise StataError(198, "option sindicators() not implemented in OpenDTA")   # OpenDTA
    lags = int(o["lags"]) if o.get("lags") and o["lags"] is not True else 2
    if lags < 1:
        raise StataError(198, "lags() must be greater than zero")
    trend = _trend_opt(o)
    deps, Y, rows, times, _m = _data(s, p, p.varlist, lags, trend)
    J = johansen(Y, lags, trend)
    K, T, lam = J["K"], J["T"], J["lam"]
    ll = [loglik(J["S00"], lam, r, T, K) for r in range(K + 1)]
    parms = [n_parms(K, lags, trend, r) for r in range(K + 1)]
    trace = [-T * float(np.sum(np.log1p(-lam[r:]))) for r in range(K)]
    mx = [-T * math.log1p(-lam[r]) for r in range(K)]
    ic = [info_criteria(ll[r], T, parms[r]) for r in range(K + 1)]
    ds = s.data
    info = require_ts(ds)
    first, last = times[lags], times[-1]
    out = s.output
    show_trace = not o.get("notrace")
    two = bool(o.get("levela"))
    lev = 1 if o.get("level99") else 5
    out.write("\n                       Johansen tests for cointegration\n", "text")
    out.write(f"Trend: {trend:<48}Number of obs = {T:>7}\n", "text")
    out.write(f"Sample: {_fmt_time(ds, info.tvar, first)} - {_fmt_time(ds, info.tvar, last):<{max(1, 40 - len(_fmt_time(ds, info.tvar, first)))}}"
              f"Lags          = {lags:>7}\n", "text")
    stats = [("trace", trace)] if show_trace else []
    if o.get("max"):
        stats.append(("max", mx))
    sel = {}
    if stats:
        for kind, vals in stats:
            levels = [5, 1] if two else [lev]
            out.write("\n", "text")
            lab = "trace" if kind == "trace" else "max"
            cols = "".join(f"{(str(l) + '% critical'):>13}" for l in levels)
            out.write(f"  maximum{'':>56}{lab:>10}{cols}\n", "text")
            out.write(f"  rank       parms       LL         eigenvalue{'':>14}statistic" +
                      "".join(f"{'value':>13}" for _ in levels) + "\n", "text")
            chosen = None
            for r in range(K + 1):
                eig = f"{lam[r - 1]:>10.5f}" if r > 0 else f"{'.':>10}"
                line = f"  {r:>3}   {parms[r]:>8}{ll[r]:>14.4f}  {eig}"
                if r < K:
                    cvs = [critical_value(trend, kind, l, K - r) for l in levels]
                    mark = ""
                    if kind == "trace" and chosen is None and cvs[0] < SYS and vals[r] < cvs[0]:
                        chosen, mark = r, "*"
                    line += f"{vals[r]:>22.4f}{mark:<1}" + "".join(
                        f"{c:>13.2f}" if c < SYS else f"{'.':>13}" for c in cvs)
                out.write(line + "\n", "result")
            sel[kind] = chosen
    if o.get("ic"):
        out.write("\n  maximum\n  rank       parms       LL         eigenvalue         SBIC          HQIC          AIC\n", "text")
        best = [int(np.argmin([c[i] for c in ic])) for i in (1, 2, 0)]
        for r in range(K + 1):
            eig = f"{lam[r - 1]:>10.5f}" if r > 0 else f"{'.':>10}"
            cells = "".join(f"{ic[r][i]:>14.5f}{'*' if best[j] == r else ' '}" for j, i in enumerate((1, 2, 0)))
            out.write(f"  {r:>3}   {parms[r]:>8}{ll[r]:>14.4f}  {eig}{cells}\n", "result")
    from ..commands.matrix import Matrix
    k95 = sel.get("trace")
    if k95 is None and show_trace:
        k95 = K
    k_sbic = int(np.argmin([c[1] for c in ic]))
    k_hqic = int(np.argmin([c[2] for c in ic]))
    est = Estimates("vecrank", " ".join(deps), [], np.zeros(0), np.zeros((0, 0)), T, "z", None)
    scalars = [("N", T), ("k_eq", K), ("k_dv", K), ("tmin", first), ("tmax", last), ("n_lags", lags),
               ("k_ce95", k95 if k95 is not None else SYS), ("k_ce99", SYS), ("k_cesbic", k_sbic),
               ("k_cehqic", k_hqic)]
    mats = [("trace", Matrix(np.array(trace + [np.nan]).reshape(1, -1))),
            ("max", Matrix(np.array(mx + [np.nan]).reshape(1, -1))),
            ("lambda", Matrix(lam.reshape(1, -1))),
            ("k_rank", Matrix(np.array(parms, dtype=float).reshape(1, -1))),
            ("hqic", Matrix(np.array([c[2] for c in ic]).reshape(1, -1))),
            ("sbic", Matrix(np.array([c[1] for c in ic]).reshape(1, -1))),
            ("aic", Matrix(np.array([c[0] for c in ic]).reshape(1, -1)))]
    from .regress import cmdline
    macros = [("cmd", "vecrank"), ("cmdline", cmdline("vecrank", args)), ("trend", trend),
              ("tsfmt", ds.get(info.tvar).fmt)]
    est.names = []
    post(s, est, scalars, macros, _mask_of(s, rows), extra_mats=mats)
    s.e.pop("b", None)
    s.e.pop("V", None)


def _mask_of(s: "Session", rows: np.ndarray) -> np.ndarray:
    m = np.zeros(s.data.nobs, dtype=bool)
    m[rows] = True
    return m


def cmd_vec(s: "Session", args: str) -> None:
    from .postest import replay
    from .regress import cmdline
    if not args.strip() or args.strip().startswith(","):
        replay(s, "vec", args)
        return
    p = parse_standard(args)
    o = match_options(p.options, {"rank": 2, "lags": 2, "trend": 2, "bconstraints": 2,
                                  "aconstraints": 2, "sindicators": 2, "noreduce": 3, "level": 1,
                                  "nobtable": 3, "noidtest": 4, "alpha": 1, "pi": 2, "noptable": 4,
                                  "mai": 3, "noetable": 4, "dforce": 2, "nocnsreport": 4,
                                  "iterate": 4, "nolog": 3, "trace": 5, "tolerance": 3,
                                  "ltolerance": 4, "afrom": 2, "bfrom": 2}) if p.options.strip() else {}
    for bad in ("bconstraints", "aconstraints", "sindicators", "afrom", "bfrom"):
        if o.get(bad):
            raise StataError(198, f"option {bad}() not implemented in OpenDTA")   # OpenDTA
    r = int(o["rank"]) if o.get("rank") and o["rank"] is not True else 1
    lags = int(o["lags"]) if o.get("lags") and o["lags"] is not True else 2
    trend = _trend_opt(o)
    deps, Y, rows, times, _m = _data(s, p, p.varlist, lags, trend)
    K = Y.shape[1]
    if not 0 < r < K:
        raise StataError(198, "rank() must be between 1 and the number of variables minus one")   # VERIFY
    F_ = fit_vecm(Y, lags, trend, r)
    J, T, beta, alpha, W, coef, U, Om = (F_[k] for k in ("J", "T", "beta", "alpha", "W", "coef", "U", "Om"))
    S00, S11 = J["S00"], J["S11"]
    mpar = W.shape[1]
    nparms = n_parms(K, lags, trend, r)
    ll = loglik(S00, J["lam"], r, T, K)
    d = nparms // K
    fac = T / (T - d)
    XtXi = np.linalg.inv(W.T @ W)
    level = float(o["level"]) if o.get("level") and o["level"] is not True else float(s.settings.get("level", "95"))
    ds = s.data
    info = require_ts(ds)
    first, last = times[lags], times[-1]
    aic, sbic, hqic = info_criteria(ll, T, nparms)
    # coefficient names (the columns of W: error-correction terms, lagged differences, terms)
    ce = [f"_ce{i + 1}" for i in range(r)]
    unr = {"constant": ["_cons"], "rconstant": [], "trend": ["_cons", "_trend"], "rtrend": ["_cons"],
           "none": []}[trend]
    cn = [f"L._ce{i + 1}" for i in range(r)]
    for i in range(1, lags):
        cn += [f"{'L' if i == 1 else 'L' + str(i)}D.{dv}" for dv in deps]
    cn += unr
    names = list(cn) * K
    eqs = [f"D_{dv}" for dv in deps for _ in cn]
    b = coef.T.ravel()
    V = np.kron(Om, XtXi) * fac
    se = np.sqrt(np.maximum(np.diag(V), 0))
    eqstats = []
    for j, dv in enumerate(deps):
        cols = list(range(mpar - len(unr)))
        bj = coef[cols, j]
        Vj = Om[j, j] * fac * XtXi[np.ix_(cols, cols)]
        chi = float(bj @ np.linalg.solve(Vj, bj))
        rss = float(U[:, j] @ U[:, j])
        z0 = J["Z0"][:, j]
        tss = float(((z0 - z0.mean()) ** 2).sum())
        eqstats.append((f"D_{dv}", mpar, math.sqrt(rss / (T - mpar)) if T > mpar else SYS,
                        1 - rss / tss if tss > 0 else SYS, chi, float(chi2_tail(chi, len(cols)))))
    rows_t = _arrange_rows(deps, unr, b, se, K, mpar, ce, lags)
    # variance of beta (free elements) and the table of the cointegrating equations
    nb = beta.shape[0]
    Ominv = np.linalg.inv(Om)
    info_mat = np.kron(alpha.T @ Ominv @ alpha, S11)       # information per observation, vec(beta)
    free = [i + e * nb for e in range(r) for i in range(nb) if i >= r]
    Vb_free = np.linalg.inv(info_mat[np.ix_(free, free)]) / (T - d)
    bvec = beta.T.ravel()
    sevec = np.zeros(nb * r)
    sevec[free] = np.sqrt(np.maximum(np.diag(Vb_free), 0))
    brn = deps + (["_cons"] if trend == "rconstant" else ["_trend"] if trend == "rtrend" else [])
    est = Estimates("vec", " ".join(f"D_{d_}" for d_ in deps), names, b, V, T, "z", None)
    est.eqnames = eqs
    est.rows = rows_t
    est.level = level
    est.vcetype = ""
    est.extra = {"N": T, "ll": ll, "aic": aic, "hqic": hqic, "sbic": sbic, "det": float(np.linalg.det(Om)),
                 "sample": f"Sample:  {_fmt_time(ds, info.tvar, first)} - {_fmt_time(ds, info.tvar, last)}",
                 "eqstats": eqstats, "beta": beta, "sebeta": sevec, "free": free, "Vb": Vb_free,
                 "r": r, "brn": brn, "ce": ce, "opts": o, "nb": nb}
    est.display = display_vec
    est.vec = {"deps": deps, "lags": lags, "trend": trend, "r": r, "alpha": alpha, "beta": beta, "Om": Om,
               "coef": coef, "T": T, "Y": Y, "J": J, "d": d, "XtXi": XtXi, "times": times, "t0": times[0],
               "t1": times[-1]}
    from ..commands.matrix import Matrix
    scalars = [("N", T), ("k_rank", r), ("k_eq", K), ("k_dv", K), ("k_ce", r), ("n_lags", lags),
               ("df_m", nparms - K * 0), ("ll", ll), ("aic", aic), ("hqic", hqic), ("sbic", sbic),
               ("tmin", first), ("tmax", last), ("detsig_ml", float(np.linalg.det(Om))), ("rank", nparms)]
    for j, (dv, mp, rm, r2, chi, pv) in enumerate(eqstats):
        scalars += [(f"rmse_{j + 1}", rm), (f"r2_{j + 1}", r2), (f"chi2_{j + 1}", chi), (f"k_{j + 1}", mp)]
    mats = [("beta", Matrix(beta.T, ce, brn)), ("alpha", Matrix(alpha, deps, ce)),
            ("pi", Matrix(alpha @ beta.T, deps, brn)), ("omega", Matrix(Om, deps, deps))]
    macros = [("cmd", "vec"), ("cmdline", cmdline("vec", args)), ("title", "Vector error-correction model"),
              ("trend", trend), ("endog", " ".join(deps)), ("tvar", info.tvar),
              ("cenames", " ".join(ce)), ("eqnames", " ".join(f"D_{d_}" for d_ in deps)),
              ("properties", "b V")]
    post(s, est, scalars, macros, _mask_of(s, rows), extra_mats=mats)
    display_vec(s, est)


def _arrange_rows(deps, unr, b, se, K, mpar, ce, lags):
    """Rows of the short-run table: for each equation the error-correction terms, the lagged
    differences of each variable and the unrestricted terms."""
    out = []
    nce, nlag = len(ce), lags - 1
    for j, dv in enumerate(deps):
        out.append(CoefRow("eq", f"D_{dv}"))
        base = j * mpar
        for k, c in enumerate(ce):
            out.append(CoefRow("header", c))
            out.append(CoefRow("coef", "L1.", b[base + k], se[base + k]))
        if nlag:
            for vi, vn in enumerate(deps):
                out.append(CoefRow("blank"))
                out.append(CoefRow("header", vn))
                for l in range(1, nlag + 1):
                    pos = base + nce + (l - 1) * K + vi
                    out.append(CoefRow("coef", "LD." if l == 1 else f"L{l}D.", b[pos], se[pos]))
        pos0 = base + nce + nlag * K
        for q, nm in enumerate(unr):
            if q == 0:
                out.append(CoefRow("blank"))
            out.append(CoefRow("coef", nm, b[pos0 + q], se[pos0 + q]))
        if j < K - 1:
            out.append(CoefRow("sep"))
    return out


def display_vec(s: "Session", est: Estimates, *, header: bool = True, table: bool = True,
                level: float | None = None, **_kw) -> None:
    from .tsmodels import _g, _right
    out = s.output
    x = est.extra
    o = x["opts"]
    lev = level or est.level
    if header:
        out.write("\nVector error-correction model\n\n", "text")
        out.write(f"{x['sample']:<48}{'No. of obs':<18}= {int(x['N']):>9}\n", "text")
        out.write(f"{'':<48}{'AIC':<18}= {_g(x['aic'], '%9.0g'):>10}\n", "text")
        out.write(f"{'Log likelihood = ' + _g(x['ll'], '%9.0g'):<48}{'HQIC':<18}= {_g(x['hqic'], '%9.0g'):>10}\n", "text")
        out.write(f"{'Det(Sigma_ml)  = ' + _g(x['det'], '%9.0g'):<48}{'SBIC':<18}= {_g(x['sbic'], '%9.0g'):>10}\n\n", "text")
        out.write(f"{'Equation':<15}{'Parms':>9}{'RMSE':>10}{'R-sq':>9}{'chi2':>10}{'P>chi2':>11}\n", "text")
        out.write("-" * 64 + "\n", "text")
        for dv, mp, rm, r2, chi, pv in x["eqstats"]:
            out.write(f"{dv[:15]:<15}", "text")
            out.write(f"{mp:>8}{_g(rm, '%8.0g'):>12}{r2:>9.4f}{_g(chi, '%9.0g'):>11}{pv:>9.4f}\n", "result")
        out.write("-" * 64 + "\n\n", "text")
    if table and not o.get("noetable"):
        coef_table(s, "", est.rows, stat="z", df=None, level=lev)
    if table and not o.get("nobtable"):
        r, nb, beta = x["r"], x["nb"], x["beta"]
        out.write("\nCointegrating equations\n\n", "text")
        out.write(f"{'Equation':<15}{'Parms':>9}{'chi2':>10}{'P>chi2':>11}\n", "text")
        out.write("-" * 45 + "\n", "text")
        free, Vb = x["free"], x["Vb"]
        bvec = beta.T.ravel()
        for e_, c in enumerate(x["ce"]):
            ids = [free.index(i + e_ * nb) for i in range(r, nb)]
            bj = bvec[[i + e_ * nb for i in range(r, nb)]]
            chi = float(bj @ np.linalg.solve(Vb[np.ix_(ids, ids)], bj))
            out.write(f"{c:<15}{len(ids):>8}{_g(chi, '%11.0g'):>12}{float(chi2_tail(chi, len(ids))):>9.4f}\n", "result")
        out.write("-" * 45 + "\n\n", "text")
        if not o.get("noidtest"):
            out.write("Identification:  beta is exactly identified\n\n", "text")       # VERIFY layout
            out.write(f"{'Johansen normalization restriction imposed':>60}\n", "text")
        rows = []
        for e_, c in enumerate(x["ce"]):
            rows.append(CoefRow("eq", c))
            for i, nm in enumerate(x["brn"]):
                k = i + e_ * nb
                if i < r:
                    rows.append(CoefRow("coef", nm, float(bvec[k]), float("nan")))
                else:
                    rows.append(CoefRow("coef", nm, float(bvec[k]),
                                        float(est.extra["sebeta"][k])))
        coef_table(s, "beta", rows, stat="z", df=None, level=lev)
