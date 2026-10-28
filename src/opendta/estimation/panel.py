"""Panel and absorbed effects: xtreg (fe, re, be), areg ([XT] xtreg, [R] areg),
and instrumental variables: ivregress 2sls ([R] ivregress)."""

from __future__ import annotations

import math

import re
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.formats import format_value
from ..core.stats import chi2_tail as _chi2_tail
from ..lang.syntax import match_options
from ..lazy import lazy
from . import fvars
from .regress import Sample, cmdline, comma, sig5, wexp_macros
from .results import CoefRow, Estimates, coef_table, post, rows_from_columns
from .tsops import panel_info

# imported on first use (start-up time; see opendta/lazy.py)
st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def g(x: float, fmt: str) -> str:
    if x is None or not np.isfinite(x) or abs(x) >= SYS:
        return "."
    return format_value(float(x), fmt, pad=False).strip()


def _demean(v: np.ndarray, ids: np.ndarray, inv: np.ndarray, counts: np.ndarray) -> np.ndarray:
    sums = np.bincount(inv, weights=v, minlength=len(counts))
    return v - (sums / counts)[inv]


def _group_means(v, inv, counts):
    return np.bincount(inv, weights=v, minlength=len(counts)) / counts


def _maps(cols, b, V, constant: bool):
    bmap, semap = {}, {}
    j = 0
    for c in cols:
        if c.values is not None:
            bmap[c.name] = float(b[j])
            semap[c.name] = float(np.sqrt(max(V[j, j], 0)))
            j += 1
        else:
            bmap[c.name], semap[c.name] = 0.0, 0.0
    if constant:
        bmap["_cons"] = float(b[-1])
        semap["_cons"] = float(np.sqrt(max(V[-1, -1], 0)))
    return bmap, semap


def _names(cols, constant):
    out = [c.name if not c.omitted else fvars.omitted_name(c.name) for c in cols]
    return out + (["_cons"] if constant else [])


def _full(cols, constant, b, V):
    names = _names(cols, constant)
    kk = len(names)
    pos = [i for i, c in enumerate(cols) if c.values is not None] + ([kk - 1] if constant else [])
    bf = np.zeros(kk)
    Vf = np.zeros((kk, kk))
    bf[pos] = b
    Vf[np.ix_(pos, pos)] = V
    return names, bf, Vf


def _cluster_vce(Xs, e, ids, A, N, k, g_count_adjust=True):
    uniq = np.unique(ids)
    G = len(uniq)
    sc = Xs * e[:, None]
    meat = np.zeros((Xs.shape[1], Xs.shape[1]))
    for u in uniq:
        sg = sc[ids == u].sum(axis=0)
        meat += np.outer(sg, sg)
    q = (G / (G - 1)) * ((N - 1) / (N - k)) if g_count_adjust else G / (G - 1)
    return q * A @ meat @ A, G


# ---------------------------------------------------------------------------
# xtreg
# ---------------------------------------------------------------------------

def cmd_xtreg(s: "Session", args: str) -> None:
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, "xtreg", args)
        return
    from .xtgee import run_pa
    if run_pa(s, "xtreg", args):
        return
    from .xtgee import has_option
    if has_option(args, "mle"):
        from .panel2 import cmd_xtreg_mle
        cmd_xtreg_mle(s, args)
        return
    smp = Sample(s, args, opts_spec={"fe": 2, "re": 2, "be": 2, "mle": 3, "i": 1, "sa": 2, "wls": 3},
                 weights=("aweight", "fweight", "pweight"))
    o = smp.o
    ds = s.data
    pvar = str(o.get("i")).strip() if o.get("i") else panel_info(ds)
    if not pvar:
        raise StataError(459, "panel variable not set; use xtset varname ...")   # VERIFY
    ids_all = ds.get(pvar).data
    mask = smp.mask & (ids_all < SYS)
    if not np.array_equal(mask, smp.mask):
        raise StataError(459, f"{pvar} missing in sample")   # VERIFY
    ids = ids_all[smp.mask]
    uniq, inv, counts = np.unique(ids, return_inverse=True, return_counts=True)
    G = len(uniq)
    y = smp.y.astype(float)
    notes = fvars.drop_collinear(smp.cols, True)
    out = s.output
    for nm in notes:
        out.write(f"note: {nm} omitted because of collinearity\n", "text")
    est_cols = [c for c in smp.cols if c.values is not None]
    Xs = np.column_stack([c.values for c in est_cols]) if est_cols else np.zeros((len(y), 0))
    N = len(y)
    k = Xs.shape[1]
    model = "re" if not (o.get("fe") or o.get("be")) else ("fe" if o.get("fe") else "be")
    if model == "fe":
        _xtreg_fe(s, smp, args, y, Xs, est_cols, ids, inv, counts, G, N, k, pvar)
    elif model == "re":
        _xtreg_re(s, smp, args, y, Xs, est_cols, ids, inv, counts, G, N, k, pvar)
    else:
        _xtreg_be(s, smp, args, y, Xs, est_cols, ids, inv, counts, G, N, k, pvar)


def _xtreg_be(s, smp, args, y, Xs, est_cols, ids, inv, counts, G, N, k, pvar):
    """Between estimator: OLS of the panel means of y on the panel means of x ([XT] xtreg;
    wls weights each panel by its number of observations)."""
    if smp.vce in ("robust", "cluster"):
        raise StataError(198, "option vce() not allowed")   # VERIFY
    yb = _group_means(y, inv, counts)
    Xb = np.column_stack([_group_means(Xs[:, j], inv, counts) for j in range(k)] + [np.ones(G)])
    wg = counts.astype(float) if smp.o.get("wls") else np.ones(G)
    sw = np.sqrt(wg)
    A = np.linalg.pinv((Xb * wg[:, None]).T @ Xb)
    b = A @ ((Xb * wg[:, None]).T @ yb)
    e = (yb - Xb @ b) * sw
    df_r = G - (k + 1)
    if df_r <= 0:
        raise StataError(2001, "insufficient observations")
    rss = float(e @ e)
    s2 = rss / df_r
    if smp.o.get("wls"):
        s2 = s2 * G / wg.sum()        # VERIFY: scale of the wls variance
    V = s2 * A
    V = (V + V.T) / 2
    F = float(b[:k] @ np.linalg.solve(V[:k, :k], b[:k])) / k if k else SYS
    pF = float(st.f.sf(F, k, df_r)) if F < SYS else SYS
    xb_all = Xs @ b[:k]
    yw = _demean(y, ids, inv, counts)
    r2_w = float(np.corrcoef(_demean(xb_all, ids, inv, counts), yw)[0, 1] ** 2) if k else 0.0
    r2_b = float(np.corrcoef(_group_means(xb_all, inv, counts), yb)[0, 1] ** 2) if k else 0.0
    r2_o = float(np.corrcoef(xb_all, y)[0, 1] ** 2) if k else 0.0
    names, bf, Vf = _full(smp.cols, True, b, V)
    bmap, semap = _maps(smp.cols, b, V, True)
    rows = rows_from_columns(smp.cols, bmap, semap, s.data, constant=True)
    est = Estimates("xtreg", smp.depname, names, bf, Vf, N, "t", df_r)
    est.rows = rows
    est.cols, est.constant, est.terms, est.depcomp = smp.cols, True, smp.terms, smp.depcomp
    est.level = smp.level
    est.vcetype, est.clustvar, est.n_clust = "", "", None
    est.model = "be"
    sd = math.sqrt(s2)
    est.extra = dict(N=N, G=G, tmin=int(counts.min()), tmax=int(counts.max()), tbar=N / G, r2_w=r2_w,
                     r2_b=r2_b, r2_o=r2_o, F=F, pF=pF, df_m=k, df_r=df_r, pvar=pvar, sd_ue=sd)
    est.display = display_xtreg
    est.predict = _predict_xt
    est.predict_opts = {"u": 1, "e": 1, "ue": 2, "xbu": 3}
    est.panel = (pvar, ids)
    scalars = [("N", float(N)), ("N_g", float(G)), ("mss", float(((Xb @ b - yb.mean()) ** 2 * wg).sum())),
               ("df_m", float(k)), ("rss", rss), ("df_r", float(df_r)), ("r2", r2_b), ("r2_a", SYS),
               ("F", F), ("rmse", sd), ("r2_w", r2_w), ("r2_o", r2_o), ("r2_b", r2_b),
               ("g_min", float(counts.min())), ("g_avg", N / G), ("g_max", float(counts.max())),
               ("Tbar", N / G), ("rank", float(k + 1))]
    macros = [("cmdline", cmdline("xtreg", args)), ("cmd", "xtreg"), ("depvar", smp.depname),
              ("ivar", pvar), ("model", "be"), ("title", "Between regression"),
              ("predict", "xtrefe_p"), ("properties", "b V")]
    post(s, est, scalars, macros, smp.mask)
    display_xtreg(s, est)


def _drop_within_collinear(smp, Xw, est_cols, out):
    """Variables constant within groups drop out in fe (collinearity note)."""
    keep = []
    for j, c in enumerate(est_cols):
        if np.allclose(Xw[:, j], 0, atol=1e-12 * (np.abs(c.values).max() + 1)):
            c.omitted = True
            c.values = None
            out.write(f"note: {c.name} omitted because of collinearity\n", "text")
        else:
            keep.append(j)
    return keep


def _xtreg_fe(s, smp, args, y, Xs, est_cols, ids, inv, counts, G, N, k, pvar):
    out = s.output
    yw = _demean(y, ids, inv, counts)
    Xw = np.column_stack([_demean(Xs[:, j], ids, inv, counts) for j in range(k)]) if k else Xs
    keep = _drop_within_collinear(smp, Xw, est_cols, out)
    Xw, Xs = Xw[:, keep], Xs[:, keep]
    est_cols = [est_cols[j] for j in keep]
    k = Xs.shape[1]
    A = np.linalg.inv(Xw.T @ Xw) if k else np.zeros((0, 0))
    b = A @ (Xw.T @ yw) if k else np.zeros(0)
    e = yw - Xw @ b
    rss = float(e @ e)
    df_r = N - G - k
    s2 = rss / df_r
    xbar, ybar = Xs.mean(axis=0), y.mean()
    cons = ybar - float(xbar @ b)
    # V with the constant: regression of (y - ȳi + ȳ) on (x - x̄i + x̄) and 1
    Xa = np.column_stack([Xw + xbar, np.ones(N)])
    Aa = np.linalg.inv(Xa.T @ Xa)
    V = s2 * Aa
    bb = np.append(b, cons)
    vcetype = ""
    n_clust = None
    if smp.vce in ("robust", "cluster"):
        # xtreg, vce(robust) = cluster on the panel (Stata)
        cl = smp.cluster if smp.vce == "cluster" else ids
        ea = yw - Xw @ b
        uniq = np.unique(cl)
        n_clust = len(uniq)
        sc = Xa * ea[:, None]
        meat = np.zeros((k + 1, k + 1))
        for u in uniq:
            sg = sc[cl == u].sum(axis=0)
            meat += np.outer(sg, sg)
        q = (n_clust / (n_clust - 1)) * ((N - 1) / (N - k - 1))
        V = q * Aa @ meat @ Aa
        vcetype = "Robust"
        df_r = n_clust - 1
    V = (V + V.T) / 2
    # medidas
    xb_all = Xs @ b
    u_i = _group_means(y, inv, counts) - _group_means(xb_all, inv, counts) - cons
    sigma_u = float(np.std(u_i, ddof=1)) if G > 1 else 0.0
    sigma_e = float(np.sqrt(s2 if smp.vce == "ols" else rss / (N - G - k)))
    rho = sigma_u ** 2 / (sigma_u ** 2 + sigma_e ** 2)
    r2_w = float(np.corrcoef(Xw @ b, yw)[0, 1] ** 2) if k else 0.0
    r2_b = float(np.corrcoef(_group_means(xb_all, inv, counts), _group_means(y, inv, counts))[0, 1] ** 2) if k else 0.0
    r2_o = float(np.corrcoef(xb_all, y)[0, 1] ** 2) if k else 0.0
    corr = float(np.corrcoef(u_i[inv], xb_all)[0, 1]) if k else 0.0
    if smp.vce == "ols":
        F = float(b @ np.linalg.solve(V[:k, :k], b)) / k if k else SYS
    else:
        F = float(b @ np.linalg.solve(V[:k, :k], b)) / k if k else SYS
    pF = float(st.f.sf(F, k, df_r)) if F < SYS else SYS
    # F test that all u_i = 0: pooled OLS vs fe
    Xp = np.column_stack([Xs, np.ones(N)])
    bp = np.linalg.lstsq(Xp, y, rcond=None)[0]
    rss_p = float(np.sum((y - Xp @ bp) ** 2))
    F_f = ((rss_p - rss) / (G - 1)) / (rss / (N - G - k))
    p_f = float(st.f.sf(F_f, G - 1, N - G - k))
    names, bf, Vf = _full(smp.cols, True, bb, V)
    bmap, semap = _maps(smp.cols, bb, V, True)
    rows = rows_from_columns(smp.cols, bmap, semap, s.data, constant=True)
    est = Estimates("xtreg", smp.depname, names, bf, Vf, N, "t", df_r)
    est.rows = rows
    est.cols, est.constant, est.terms, est.depcomp = smp.cols, True, smp.terms, smp.depcomp
    est.level = smp.level
    est.vcetype, est.clustvar, est.n_clust = vcetype, (smp.clustvar or pvar) if vcetype else "", n_clust
    est.model = "fe"
    est.extra = dict(N=N, G=G, tmin=int(counts.min()), tmax=int(counts.max()), tbar=N / G, r2_w=r2_w,
                     r2_b=r2_b, r2_o=r2_o, F=F, pF=pF, df_m=k, df_r=df_r, corr=corr, sigma_u=sigma_u,
                     sigma_e=sigma_e, rho=rho, F_f=F_f, p_f=p_f, pvar=pvar)
    est.display = display_xtreg
    est.predict = _predict_xt
    est.predict_opts = {"u": 1, "e": 1, "ue": 2, "xbu": 3}
    # e(): names and order observed in Stata 14 (compat 0510)
    tss_w = float(yw @ yw)
    ll = -0.5 * N * (np.log(2 * np.pi) + np.log(rss / N) + 1)
    ll_0 = -0.5 * N * (np.log(2 * np.pi) + np.log(tss_w / N) + 1)
    tss = float(np.sum((y - ybar) ** 2))
    r2_a = 1 - (1 - r2_w) * (N - 1) / (N - G - k)
    scalars = [("r2_w", r2_w), ("F_f", F_f), ("rank", float(k + 1)), ("Tbar", N / G),
               ("Tcon", 1.0 if counts.min() == counts.max() else 0.0), ("g_min", float(counts.min())),
               ("g_avg", N / G), ("g_max", float(counts.max())), ("N_g", float(G)), ("sigma_u", sigma_u),
               ("corr", corr), ("r2_o", r2_o), ("r2_b", r2_b), ("sigma_e", sigma_e),
               ("sigma", float(np.sqrt(sigma_u ** 2 + sigma_e ** 2))), ("rho", rho), ("df_m", float(k + G - 1)),
               ("df_b", float(k)), ("ll_0", ll_0), ("ll", ll), ("tss", tss), ("df_a", float(G - 1)),
               ("r2_a", r2_a), ("rss", rss), ("mss", tss_w - rss), ("rmse", sigma_e), ("r2", r2_w), ("F", F),
               ("df_r", float(df_r)), ("N", float(N))]
    if n_clust is not None:
        scalars.insert(scalars.index(("N", float(N))), ("N_clust", float(n_clust)))   # VERIFY position
    macros = [("cmdline", cmdline("xtreg", args)), ("cmd", "xtreg"), ("marginsnotok", "E U UE SCore STDP XBU"),
              ("predict", "xtrefe_p"), ("model", "fe"), ("vce", smp.vce if vcetype else "conventional")]
    if vcetype:
        macros += [("vcetype", vcetype), ("clustvar", est.clustvar)]
    macros += [("ivar", pvar), ("depvar", smp.depname), ("properties", "b V"), ("absvar", pvar)]
    post(s, est, scalars, macros, smp.mask)
    display_xtreg(s, est)


def _xtreg_re(s, smp, args, y, Xs, est_cols, ids, inv, counts, G, N, k, pvar):
    # Swamy-Arora: σe² from within, σu² from between (harmonic T̄)
    yw = _demean(y, ids, inv, counts)
    Xw = np.column_stack([_demean(Xs[:, j], ids, inv, counts) for j in range(k)]) if k else Xs
    keepw = [j for j in range(k) if not np.allclose(Xw[:, j], 0)]
    Aw = np.linalg.inv(Xw[:, keepw].T @ Xw[:, keepw]) if keepw else np.zeros((0, 0))
    bw = Aw @ (Xw[:, keepw].T @ yw) if keepw else np.zeros(0)
    ew = yw - Xw[:, keepw] @ bw
    sig_e2 = float(ew @ ew) / (N - G - len(keepw))
    yb = _group_means(y, inv, counts)
    Xb = np.column_stack([_group_means(Xs[:, j], inv, counts) for j in range(k)] + [np.ones(G)])
    bb_, *_ = np.linalg.lstsq(Xb, yb, rcond=None)
    eb = yb - Xb @ bb_
    Tbar = G / np.sum(1.0 / counts)
    sig_u2 = max(float(eb @ eb) / (G - (k + 1)) - sig_e2 / Tbar, 0.0)
    theta_i = 1 - np.sqrt(sig_e2 / (counts * sig_u2 + sig_e2))
    th = theta_i[inv]
    ys = y - th * yb[inv]
    Xstar = np.column_stack([Xs - th[:, None] * Xb[inv, :k], 1 - th])
    A = np.linalg.inv(Xstar.T @ Xstar)
    b = A @ (Xstar.T @ ys)
    e = ys - Xstar @ b
    V = (float(e @ e) / (N - (k + 1))) * A   # VERIFY: s² of the transformed regression
    vcetype = ""
    n_clust = None
    if smp.vce in ("robust", "cluster"):
        cl = smp.cluster if smp.vce == "cluster" else ids
        V, n_clust = _cluster_vce(Xstar, e, cl, A, N, k + 1, g_count_adjust=False)
        vcetype = "Robust"
    V = (V + V.T) / 2
    chi2 = float(b[:k] @ np.linalg.solve(V[:k, :k], b[:k])) if k else 0.0
    p = float(_chi2_tail(chi2, k)) if k else SYS
    xb_all = Xs @ b[:k]
    r2_w = float(np.corrcoef(_demean(xb_all, ids, inv, counts), yw)[0, 1] ** 2) if k else 0.0
    r2_b = float(np.corrcoef(_group_means(xb_all, inv, counts), yb)[0, 1] ** 2) if k else 0.0
    r2_o = float(np.corrcoef(xb_all, y)[0, 1] ** 2) if k else 0.0
    sigma_u, sigma_e = np.sqrt(sig_u2), np.sqrt(sig_e2)
    rho = sig_u2 / (sig_u2 + sig_e2)
    names, bf, Vf = _full(smp.cols, True, b, V)
    bmap, semap = _maps(smp.cols, b, V, True)
    rows = rows_from_columns(smp.cols, bmap, semap, s.data, constant=True)
    est = Estimates("xtreg", smp.depname, names, bf, Vf, N, "z", None)
    est.rows = rows
    est.cols, est.constant, est.terms, est.depcomp = smp.cols, True, smp.terms, smp.depcomp
    est.level = smp.level
    est.vcetype, est.clustvar, est.n_clust = vcetype, (smp.clustvar or pvar) if vcetype else "", n_clust
    est.model = "re"
    est.extra = dict(N=N, G=G, tmin=int(counts.min()), tmax=int(counts.max()), tbar=N / G, r2_w=r2_w,
                     r2_b=r2_b, r2_o=r2_o, chi2=chi2, p=p, df_m=k, sigma_u=sigma_u, sigma_e=sigma_e, rho=rho,
                     pvar=pvar, theta=float(theta_i[0]) if counts.min() == counts.max() else None)
    est.display = display_xtreg
    est.predict = _predict_xt
    est.predict_opts = {"u": 1, "e": 1, "ue": 2, "xbu": 3}
    scalars = [("N", float(N)), ("N_g", float(G)), ("df_m", float(k)), ("sigma_u", sigma_u), ("sigma_e", sigma_e),
               ("r2_w", r2_w), ("r2_o", r2_o), ("r2_b", r2_b), ("chi2", chi2), ("rho", rho),
               ("Tbar", Tbar), ("g_min", float(counts.min())), ("g_avg", N / G), ("g_max", float(counts.max())),
               ("rank", float(k + 1))]
    macros = [("cmdline", cmdline("xtreg", args)), ("cmd", "xtreg"), ("depvar", smp.depname), ("ivar", pvar),
              ("model", "re"), ("title", "Random-effects GLS regression"), ("chi2type", "Wald"),
              ("vce", smp.vce if vcetype else "conventional"), ("properties", "b V"), ("predict", "xtrefe_p")]
    post(s, est, scalars, macros, smp.mask)
    display_xtreg(s, est)


def display_xtreg(s: "Session", est: Estimates, *, header: bool = True, table: bool = True,
                  level: float | None = None, **_kw) -> None:
    out = s.output
    x = est.extra
    if header:
        out.write("\n", "text")
        title = {"fe": "Fixed-effects (within) regression", "re": "Random-effects GLS regression",
                 "be": "Between regression (regression on group means)"}[est.model]
        right_top = [("Number of obs", comma(x["N"])), ("Number of groups", comma(x["G"]))]
        left = [title, f"Group variable: {x['pvar']}"]
        for lft, (rl, rv) in zip(left, right_top):
            out.write(f"{lft:<48}{rl:<18}= ", "text")
            out.write(f"{rv:>10}\n", "result")
        out.write("\n", "text")
        out.write(f"{'R-sq:':<48}Obs per group:\n", "text")
        stats = [("within", x["r2_w"], "min", f"{x['tmin']}"), ("between", x["r2_b"], "avg", f"{x['tbar']:.1f}"),
                 ("overall", x["r2_o"], "max", f"{x['tmax']}")]
        for nm, val, gl, gv in stats:
            out.write(f"     {nm:<7} = ", "text")
            out.write(f"{val:.4f}", "result")
            out.write(" " * (62 - 21) + f"{gl} = ", "text")      # compat 0510
            out.write(f"{gv:>10}\n", "result")
        out.write("\n", "text")
        if est.model in ("fe", "be"):
            fl = f"F({int(x['df_m'])},{int(x['df_r'])})"
            out.write(f"{'':<48}{fl:<18}= ", "text")
            out.write(f"{g(x['F'], '%10.2f'):>10}\n", "result")
            left = ('corr(u_i, Xb)  = ' + format(x['corr'], '.4f')) if est.model == "fe" else \
                ("sd(u_i + avg(e_i.))= " + format_value(x["sd_ue"], "%9.0g"))
            out.write(f"{left:<48}{'Prob > F':<18}= ", "text")
            out.write(f"{g(x['pF'], '%10.4f'):>10}\n", "result")
        else:
            fl = f"Wald chi2({int(x['df_m'])})"
            out.write(f"{'':<48}{fl:<18}= ", "text")
            out.write(f"{g(x['chi2'], '%10.2f'):>10}\n", "result")
            out.write(f"{'corr(u_i, X)   = 0 (assumed)':<48}{'Prob > chi2':<18}= ", "text")
            out.write(f"{g(x['p'], '%10.4f'):>10}\n", "result")
        out.write("\n", "text")
    if not table:
        return
    if est.vcetype and est.clustvar:
        out.write(f"{'(Std. Err. adjusted for ' + format(int(est.n_clust), ',') + ' clusters in ' + est.clustvar + ')':>78}\n",
                  "text")
    if est.model == "be":
        coef_table(s, est.depvar, est.rows, stat=est.stat, df=est.df_r, level=level or est.level)
        return
    coef_table(s, est.depvar, est.rows, stat=est.stat, df=est.df_r, level=level or est.level,
               vcetype=est.vcetype, footer=False)
    out.write("-" * 13 + "+" + "-" * 64 + "\n", "text")
    for nm, val in (("sigma_u", x["sigma_u"]), ("sigma_e", x["sigma_e"]), ("rho", x["rho"])):
        out.write(f"{nm:>12} | ", "text")
        txt = format_value(val, "%10.0g", pad=False).strip() if val == val else "."
        out.write(f"{txt:>10}", "result")
        if nm == "rho":
            out.write("   (fraction of variance due to u_i)", "text")
        out.write("\n", "text")
    out.write("-" * 78 + "\n", "text")
    if est.model == "fe" and not est.vcetype:
        out.write(f"F test that all u_i=0: F({int(x['G'] - 1)}, {int(x['df_r'])}) = {x['F_f']:.2f}".ljust(61)
                  + f"Prob > F = {x['p_f']:.4f}\n", "text")   # VERIFY alinhamento


def _predict_xt(s, est, stat, mask):
    """predict after xtreg: xb, stdp, ue (u_i + e_it), u, e and xbu ([XT] xtreg postestimation).
    u_i: mean residual of the panel (fe and be) or its shrunken estimate T s_u^2 / (T s_u^2 + s_e^2)
    times the mean residual (re; VERIFY)."""
    from .postest import design_all, xb_all
    xb, ok = xb_all(s, est)
    if stat in ("", "xb"):
        note = "(option xb assumed; fitted values)" if not stat else ""
        return np.where(ok, xb, SYS), note
    if stat == "stdp":
        X, _ = design_all(s, est)
        return np.where(ok, np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", X, est.V, X), 0)), SYS), ""
    if stat not in ("ue", "u", "e", "xbu"):
        raise StataError(198, f"option {stat} not allowed")
    ds = s.data
    yv = fvars.comp_values(ds, est.depcomp)
    pvar = est.extra["pvar"]
    ids = ds.get(pvar).data
    smp = ds.esample if getattr(ds, "esample", None) is not None and len(ds.esample) == ds.nobs \
        else np.ones(ds.nobs, bool)
    good = ok & (yv < SYS) & (ids < SYS)
    res = np.where(good, yv - xb, 0.0)
    use = good & smp
    uniq, inv = np.unique(np.where(ids < SYS, ids, -1e300), return_inverse=True)
    tot = np.bincount(inv, weights=np.where(use, res, 0.0), minlength=len(uniq))
    cnt = np.bincount(inv, weights=use.astype(float), minlength=len(uniq))
    with np.errstate(invalid="ignore", divide="ignore"):
        u = np.where(cnt > 0, tot / cnt, 0.0)
    if est.model == "re":
        su2, se2 = est.extra["sigma_u"] ** 2, est.extra["sigma_e"] ** 2
        with np.errstate(invalid="ignore", divide="ignore"):
            u = u * (cnt * su2 / (cnt * su2 + se2))
    ui = u[inv]
    if stat == "ue":
        return np.where(good, yv - xb, SYS), ""
    if stat == "u":
        return np.where(ok & (ids < SYS) & (cnt[inv] > 0), ui, SYS), ""
    if stat == "e":
        return np.where(good & (cnt[inv] > 0), yv - xb - ui, SYS), ""
    return np.where(ok & (ids < SYS) & (cnt[inv] > 0), xb + ui, SYS), ""


# ---------------------------------------------------------------------------
# areg
# ---------------------------------------------------------------------------

def cmd_areg(s: "Session", args: str) -> None:
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, "areg", args)
        return
    smp = Sample(s, args, opts_spec={"absorb": 1}, weights=("aweight", "fweight", "pweight"))
    o = smp.o
    if not o.get("absorb"):
        raise StataError(198, "option absorb() required")
    from ..core.varlist import resolve_name
    ds = s.data
    avar = resolve_name(ds, str(o["absorb"]).strip())
    av = ds.get(avar).data
    if np.any(av[smp.mask] >= SYS):
        smp.mask &= av < SYS   # VERIFY
    ids = av[smp.mask]
    uniq, inv, counts = np.unique(ids, return_inverse=True, return_counts=True)
    G = len(uniq)
    y = smp.y.astype(float)
    out = s.output
    notes = fvars.drop_collinear(smp.cols, True)
    for nm in notes:
        out.write(f"note: {nm} omitted because of collinearity\n", "text")
    est_cols = [c for c in smp.cols if c.values is not None]
    Xs = np.column_stack([c.values for c in est_cols]) if est_cols else np.zeros((len(y), 0))
    N, k = len(y), Xs.shape[1]
    yw = _demean(y, ids, inv, counts)
    Xw = np.column_stack([_demean(Xs[:, j], ids, inv, counts) for j in range(k)]) if k else Xs
    A = np.linalg.inv(Xw.T @ Xw) if k else np.zeros((0, 0))
    b = A @ (Xw.T @ yw) if k else np.zeros(0)
    e = yw - Xw @ b
    rss = float(e @ e)
    df_r = N - G - k
    s2 = rss / df_r
    xbar, ybar = Xs.mean(axis=0), y.mean()
    cons = ybar - float(xbar @ b)
    Xa = np.column_stack([Xw + xbar, np.ones(N)])
    Aa = np.linalg.inv(Xa.T @ Xa)
    V = s2 * Aa
    vcetype, n_clust = "", None
    if smp.vce in ("robust", "cluster"):
        if smp.vce == "cluster":
            V, n_clust = _cluster_vce(Xa, e, smp.cluster, Aa, N, k + G)
            df_r = n_clust - 1
        else:
            sc = Xa * e[:, None]
            V = (N / df_r) * Aa @ (sc.T @ sc) @ Aa
        vcetype = "Robust"
    V = (V + V.T) / 2
    bb = np.append(b, cons)
    tss = float(np.sum((y - ybar) ** 2))
    r2 = 1 - rss / tss
    r2_a = 1 - (1 - r2) * (N - 1) / df_r
    F = float(b @ np.linalg.solve(V[:k, :k], b)) / k if k else SYS
    pF = float(st.f.sf(F, k, df_r)) if k else SYS
    Xp = np.column_stack([Xs, np.ones(N)])
    bp = np.linalg.lstsq(Xp, y, rcond=None)[0]
    rss_p = float(np.sum((y - Xp @ bp) ** 2))
    F_abs = ((rss_p - rss) / (G - 1)) / (rss / (N - G - k))
    p_abs = float(st.f.sf(F_abs, G - 1, N - G - k))
    names, bf, Vf = _full(smp.cols, True, bb, V)
    bmap, semap = _maps(smp.cols, bb, V, True)
    rows = rows_from_columns(smp.cols, bmap, semap, ds, constant=True)
    est = Estimates("areg", smp.depname, names, bf, Vf, N, "t", df_r)
    est.rows = rows
    est.cols, est.constant, est.terms, est.depcomp = smp.cols, True, smp.terms, smp.depcomp
    est.level = smp.level
    est.vcetype, est.clustvar, est.n_clust = vcetype, smp.clustvar, n_clust
    est.extra = dict(N=N, G=G, F=F, pF=pF, df_m=k, df_r=df_r, r2=r2, r2_a=r2_a, rmse=np.sqrt(s2),
                     F_abs=F_abs, p_abs=p_abs, avar=avar)
    est.display = display_areg
    est.predict = _predict_xt
    est.predict_opts = {"u": 1, "e": 1, "ue": 2, "xbu": 3}
    scalars = [("N", float(N)), ("tss", tss), ("df_m", float(k)), ("rss", rss), ("df_r", float(df_r)),
               ("r2", r2), ("r2_a", r2_a), ("F", F), ("rmse", np.sqrt(s2)), ("df_a", float(G - 1)),
               ("F_absorb", F_abs), ("rank", float(k + 1))]
    macros = [("cmdline", cmdline("areg", args)), ("cmd", "areg"), ("depvar", smp.depname), ("absvar", avar),
              ("title", "Linear regression, absorbing indicators"),
              ("vce", smp.vce), ("properties", "b V"), ("predict", "areg_p")]
    post(s, est, scalars, macros, smp.mask)
    display_areg(s, est)


def display_areg(s, est, *, header=True, table=True, level=None, **_kw):
    out = s.output
    x = est.extra
    if header:
        out.write("\n", "text")
        right = [("Number of obs", comma(x["N"])), (f"F({int(x['df_m']):>4},{int(x['df_r']):>7})", g(x["F"], "%10.2f")),
                 ("Prob > F", g(x["pF"], "%10.4f")), ("R-squared", g(x["r2"], "%10.4f")),
                 ("Adj R-squared", g(x["r2_a"], "%10.4f")), ("Root MSE", g(x["rmse"], "%10.4f"))]
        for i, (rl, rv) in enumerate(right):
            lft = "Linear regression, absorbing indicators" if i == 0 else ""
            out.write(f"{lft:<48}{rl:<18}= ", "text")
            out.write(f"{rv:>10}\n", "result")
        out.write("\n", "text")
    if not table:
        return
    coef_table(s, est.depvar, est.rows, stat="t", df=est.df_r, level=level or est.level, vcetype=est.vcetype,
               footer=False)
    out.write("-" * 13 + "+" + "-" * 64 + "\n", "text")
    out.write(f"{x['avar']:>12} |", "text")
    # compat 0510: no closing line after the test of the indicators
    out.write(f"{'':10}F({int(x['G'] - 1)}, {int(x['N'] - x['G'] - x['df_m'])}) = {x['F_abs']:10.3f}"
              f"{x['p_abs']:8.3f}{'':10}({int(x['G'])} categories)\n", "result")


# ---------------------------------------------------------------------------
# ivregress 2sls
# ---------------------------------------------------------------------------

def cmd_ivregress(s: "Session", args: str) -> None:
    from .postest import replay
    t = args.strip()
    if not t or t.startswith(","):
        replay(s, "ivregress", args)
        return
    m = re.match(r"^(\w+)\s+(.*)$", t, re.S)
    if not m or m.group(1) not in ("2sls", "liml", "gmm"):
        raise StataError(198, "estimator must be 2sls, liml, or gmm")   # VERIFY
    estimator = m.group(1)
    rest = m.group(2)
    mp = re.search(r"\(([^=()]+)=([^()]*)\)", rest)
    if not mp:
        raise StataError(198, "equation not specified")   # VERIFY
    endog_txt, inst_txt = mp.group(1), mp.group(2)
    exog_txt = rest[:mp.start()] + " " + rest[mp.end():]
    # sample: dependent + exogenous + endogenous + instruments without missing
    smp = Sample(s, exog_txt.strip(), opts_spec={"small": 2, "first": 3, "wmatrix": 2, "center": 6, "igmm": 4},
                 weights=("aweight", "fweight", "pweight"), extra_vars=endog_txt + " " + inst_txt)
    if estimator != "gmm" and smp.o.get("wmatrix"):
        raise StataError(198, "option wmatrix() is allowed only with gmm")   # VERIFY
    for opt_ in ("center", "igmm"):
        if estimator != "gmm" and smp.o.get(opt_):
            raise StataError(198, f"option {opt_} is allowed only with gmm")   # VERIFY
    ds = s.data
    endog_terms = fvars.expand_fv(ds, endog_txt)
    inst_terms = fvars.expand_fv(ds, inst_txt)
    mask = smp.mask
    y = smp.y.astype(float)
    N = len(y)
    endog_cols = fvars.build_design(ds, endog_terms, mask)
    inst_cols = fvars.build_design(ds, inst_terms, mask)
    ex_cols = [c for c in smp.cols if c.values is not None and not c.base]
    Xn = [c for c in endog_cols if c.values is not None and not c.base]
    Zc = [c for c in inst_cols if c.values is not None and not c.base]
    X = np.column_stack([c.values for c in Xn] + [c.values for c in ex_cols] + [np.ones(N)])
    Z = np.column_stack([c.values for c in Zc] + [c.values for c in ex_cols] + [np.ones(N)])
    if Z.shape[1] < X.shape[1]:
        raise StataError(481, "equation not identified; must have at least as many instruments not in\n"
                         "the regression as there are instrumented variables")
    k = X.shape[1]
    wfit = None
    if smp.wtype:
        # weighted estimation: every row scaled by the square root of its (normalized) weight,
        # which gives weighted 2SLS, LIML and GMM (case 1909: svy: ivregress)
        wfit, N = smp.weights_for_fit()
        sq = np.sqrt(wfit)
        y, X, Z = y * sq, X * sq[:, None], Z * sq[:, None]
    PzX = Z @ np.linalg.lstsq(Z, X, rcond=None)[0]
    kappa = 1.0
    QX = PzX
    if estimator == "liml":
        kappa = _liml_kappa(y, Xn, ex_cols, X, Z)
        QX = X - kappa * (X - PzX)         # (I - kappa M_Z) X
    A = np.linalg.inv(QX.T @ X)
    b = A @ (QX.T @ y)
    e = y - X @ b
    rss = float(e @ e)
    small = bool(smp.o.get("small"))
    s2 = rss / (N - k) if small else rss / N
    V = s2 * A
    vcetype = ""
    n_clust = None
    wm_label = ""
    if estimator == "gmm":
        b, e, V, wm_label, n_clust_w, ginfo = _gmm_two_step(smp, y, X, Z, A @ (QX.T @ y), small)
        if smp.o.get("igmm") and not smp.o.get("nolog"):
            for i_, (db, dw) in enumerate(ginfo["trace"], 1):
                s.output.write(f"Iteration {i_}: change in beta = {db:10.3e}   change in W = {dw:10.3e}\n", "text")
        rss = float(e @ e)
        s2 = rss / (N - k) if small else rss / N
        if smp.vce in ("robust", "cluster") or wm_label != "Unadjusted":
            vcetype = "HAC" if wm_label.startswith("HAC") else "Robust"   # VERIFY the HAC label
            n_clust = n_clust_w if (smp.vce == "cluster" or wm_label.startswith("Cluster")) else None
    elif smp.vce in ("robust", "cluster") and getattr(s, "_svy_ctx", None) is not None:
        # svy: the linearized variance from the scores (QX e, already weighted) carried to the
        # design observations, as the other estimators under svy_generic (case 1909)
        ctx = s._svy_ctx
        dd = ctx["dd"]
        Zl = (QX * e[:, None]) @ A.T
        pos = ctx["pos"][np.flatnonzero(mask)]
        Zf = np.zeros((dd.n, Zl.shape[1]))
        Zf[pos] = Zl
        ctx["Z"] = Zf
        V = dd.variance(Zf)
        vcetype = "Robust"
    elif smp.vce in ("robust", "cluster"):
        sc = QX * e[:, None]
        if smp.vce == "cluster":
            uniq = np.unique(smp.cluster)
            n_clust = len(uniq)
            meat = sum(np.outer(sc[smp.cluster == u].sum(0), sc[smp.cluster == u].sum(0)) for u in uniq)
            q = n_clust / (n_clust - 1) * ((N - 1) / (N - k) if small else 1)
        else:
            meat = sc.T @ sc
            q = N / (N - k) if small else 1.0
        V = q * A @ meat @ A
        vcetype = "Robust"
    V = (V + V.T) / 2
    slopes = list(range(k - 1))
    W = float(b[slopes] @ np.linalg.solve(V[np.ix_(slopes, slopes)], b[slopes]))
    if wfit is None:
        ybar = y.mean()
        tss = float(np.sum((y - ybar) ** 2))
    else:
        yo = y / np.sqrt(wfit)
        ybar = float(np.sum(wfit * yo) / np.sum(wfit))
        tss = float(np.sum(wfit * (yo - ybar) ** 2))
    r2 = 1 - rss / tss
    rmse = np.sqrt(rss / (N - k) if small else rss / N)
    names = [c.name for c in Xn] + [c.name for c in ex_cols] + ["_cons"]
    bmap = dict(zip(names, b))
    semap = {nm: float(np.sqrt(max(V[i, i], 0))) for i, nm in enumerate(names)}
    rows = rows_from_columns(Xn + ex_cols, bmap, semap, ds, constant=True)
    stat = "t" if small else "z"
    df_r = N - k if small else None
    est = Estimates("ivregress", smp.depname, names, b, V, N, stat, df_r)
    est.rows = rows
    est.cols = Xn + ex_cols
    est.terms = endog_terms + smp.terms
    _retarget(est, Xn, ex_cols, endog_terms, smp.terms)
    est.constant, est.depcomp = True, smp.depcomp
    est.level = smp.level
    est.vcetype, est.clustvar, est.n_clust = vcetype, smp.clustvar, n_clust
    est.extra = dict(N=N, W=W, df_m=len(slopes), r2=r2, rmse=rmse, small=small, estimator=estimator,
                     wm_label=wm_label,
                     F=W / len(slopes), instd=" ".join(c.name for c in Xn),
                     insts=" ".join([c.name for c in ex_cols] + [c.name for c in Zc]))
    est.display = display_iv
    est.predict = _predict_xt
    est.predict_opts = {"u": 1, "e": 1, "ue": 2, "xbu": 3}
    r2_a = 1 - (1 - r2) * (N - 1) / (N - k)
    est.extra.update(mss=tss - rss, rss=rss, tss=tss, r2_a=r2_a, k=k)
    # e(): order observed in Stata 14 with small (compat 0510); without small, VERIFY
    gmm = estimator == "gmm"
    # gmm: J in place of kappa, the iterations, the HAC lags and S and W (compat 1918; the other
    # weight matrices follow the HAC case, VERIFY)
    scalars = [("J", ginfo["J"]) if gmm else ("kappa", kappa), ("rss", rss), ("N", float(N)),
               ("df_m", float(len(slopes)))]
    scalars += [("df_r", float(N - k))] if small else []
    scalars += [("rmse", rmse), ("mss", tss - rss), ("r2", r2), ("r2_a", r2_a)]
    scalars += [("F", W / len(slopes))] if small else [("chi2", W)]
    if n_clust is not None:
        scalars.append(("N_clust", float(n_clust)))
    scalars += [("iterations", float(ginfo["iterations"]) if gmm else 0.0)]
    if gmm and ginfo["hac"] is not None:
        scalars.append(("hac_lag", float(ginfo["hac"][1])))
    scalars.append(("rank", float(k)))
    exogr = " ".join(c.name for c in ex_cols)
    macros = [("cmdline", cmdline("ivregress", args)), ("cmd", "ivregress"), ("estat_cmd", "ivregress_estat")]
    if gmm:
        wtext = " ".join(str(smp.o.get("wmatrix") or "robust").split())
        if ginfo["hac"] is not None:
            macros.append(("hac_kernel", {"quadraticspectral": "Quadratic spectral"}.get(ginfo["hac"][0], ginfo["hac"][0].capitalize())))
            wtext = f"hac {ginfo['hac'][0]} {ginfo['hac'][1]}"
        macros += [("wmatrix", wtext), ("vce", wtext if ginfo["hac"] is not None else (smp.vce if vcetype else "unadjusted"))]
    else:
        macros.append(("vce", smp.vce if vcetype else "unadjusted"))
    if vcetype:
        macros.append(("vcetype", vcetype))
    if smp.clustvar:
        macros.append(("clustvar", smp.clustvar))
    macros += [("estimator", estimator), ("footnote", "ivreg_footnote"), ("marginsnotok", "Residuals SCores"),
               ("marginsok", "XB default"), ("predict", "ivregress_p"), ("depvar", smp.depname)]
    if exogr:
        macros.append(("exogr", exogr))
    macros += [("insts", est.extra["insts"]), ("instd", est.extra["instd"])]
    if small:
        macros.append(("small", "small"))
    title = f"Instrumental variables ({estimator.upper()}) regression"
    est.extra["title"] = title
    macros += [("title", title), ("properties", "b V")]
    mats = []
    if gmm:
        from ..commands.matrix import Matrix
        zn = [c.name for c in Zc] + [c.name for c in ex_cols] + ["_cons"]
        mats = [("S", Matrix(ginfo["S"], zn, zn)), ("W", Matrix(ginfo["W"], zn, zn))]
        est.extra["hac"] = ginfo["hac"]
    post(s, est, scalars, macros, mask, extra_mats=mats)
    display_iv(s, est)


def _retarget(est, Xn, ex_cols, endog_terms, ex_terms):
    """The endogenous and exogenous columns use term indices of the joined list."""
    off = len(endog_terms)
    for c in ex_cols:
        c.term = c.term + off if not getattr(c, "_retargeted", False) else c.term
        c._retargeted = True


def display_iv(s, est, *, header=True, table=True, level=None, **_kw):
    out = s.output
    x = est.extra
    if header and x["small"] and not est.vcetype and not x["wm_label"]:
        # with small, ANOVA table in the old regress layout (compat 0510)
        out.write(f"\n\n{x['title']}\n\n", "text")
        df_m, df_r, N = int(x["df_m"]), int(est.df_r), x["N"]
        right = [("Number of obs", comma(N)), (f"F({df_m:>3},{df_r:>6})", g(x["F"], "%10.2f")),
                 ("Prob > F", g(st.f.sf(x["F"], df_m, df_r), "%10.4f")), ("R-squared", g(x["r2"], "%10.4f")),
                 ("Adj R-squared", g(x["r2_a"], "%10.4f")), ("Root MSE", g(x["rmse"], "%10.4f"))]

        def row(nm, ss, df, ms):
            return f"{nm:>12} | {g(ss, '%11.0g'):>11}{int(df):>6}{g(ms, '%11.0g'):>12}"
        left = [f"{'Source':>12} |       SS       df       MS", "-" * 13 + "+" + "-" * 30,
                row("Model", x["mss"], df_m, x["mss"] / df_m), row("Residual", x["rss"], df_r, x["rss"] / df_r),
                "-" * 13 + "+" + "-" * 30, row("Total", x["tss"], N - 1, x["tss"] / (N - 1))]
        for lft, (rl, rv) in zip(left, right):
            out.write(f"{lft:<50}{rl:<16}= ", "text")
            out.write(f"{rv:>10}\n", "result")
        out.write("\n", "text")
    elif header:
        out.write("\n", "text")
        if x["small"]:
            right = [("Number of obs", comma(x["N"])), (f"F({int(x['df_m'])}, {int(est.df_r)})", g(x["F"], "%10.2f")),
                     ("Prob > F", g(st.f.sf(x["F"], x["df_m"], est.df_r), "%10.4f"))]
        else:
            right = [("Number of obs", comma(x["N"])), (f"Wald chi2({int(x['df_m'])})", g(x["W"], "%10.2f")),
                     ("Prob > chi2", g(_chi2_tail(x["W"], x["df_m"]), "%10.4f"))]
        right += [("R-squared", g(x["r2"], "%10.4f")), ("Root MSE", sig5(x["rmse"]))]
        for i, (rl, rv) in enumerate(right):
            lft = x["title"] if i == 0 else ""
            if x["wm_label"] and i == len(right) - 1:
                lft = f"GMM weight matrix: {x['wm_label']}"
            out.write(f"{lft:<50}{rl:<16}= ", "text")
            out.write(f"{rv:>10}\n", "result")
        out.write("\n", "text")
    if not table:
        return
    coef_table(s, est.depvar, est.rows, stat=est.stat, df=est.df_r, level=level or est.level,
               vcetype=est.vcetype)
    out.write(f"Instrumented:  {x['instd']}\n", "text")
    out.write(f"Instruments:   {x['insts']}\n", "text")
    if x.get("hac") is not None:
        kern, lags = x["hac"][0], x["hac"][1]
        kname = {"quadraticspectral": "Quadratic spectral"}.get(kern, kern.capitalize())
        out.write(f"HAC VCE:       {kname} kernel with {lags} lags\n", "text")


def _liml_kappa(y, Xn, ex_cols, X, Z) -> float:
    """Smallest eigenvalue of (Y' M_Z Y)^-1 (Y' M_X1 Y), Y = [y, endogenous regressors], X1 =
    included exogenous variables with the constant (kappa >= 1)."""
    from scipy.linalg import eigh
    g = len(Xn)
    Y = np.column_stack([y, X[:, :g]])
    X1 = X[:, g:]

    def resid(M, B):
        return B - M @ np.linalg.lstsq(M, B, rcond=None)[0]
    A_ = resid(X1, Y)
    B_ = resid(Z, Y)
    return float(eigh(A_.T @ A_, B_.T @ B_, eigvals_only=True)[0])


def _hac_weights(kernel: str, lags: int, N: int) -> list[float]:
    """Weights of the autocovariances 1..N-1 for wmatrix(hac kernel #) ([R] ivregress)."""
    out = []
    for j in range(1, N):
        if kernel == "quadraticspectral":
            x = j / (lags + 1)                                # bandwidth # + 1, as the other kernels (compat 1918)
            if not np.isfinite(x):
                out.append(0.0)
                continue
            u = 6 * np.pi * x / 5
            out.append(float(25 / (12 * np.pi ** 2 * x ** 2) * (np.sin(u) / u - np.cos(u))))
            continue
        x = j / (lags + 1)
        if x > 1:
            out.append(0.0)
        elif kernel == "bartlett":
            out.append(1 - x)
        else:                                             # parzen
            out.append(1 - 6 * x * x + 6 * x ** 3 if x <= .5 else 2 * (1 - x) ** 3)
    return out


def _gmm_two_step(smp, y, X, Z, b1, small):
    """Two-step (or, with igmm, iterated) efficient GMM: the weight matrix comes from the 2SLS
    residuals (robust, cluster, HAC or unadjusted; center centers the moments); returns b,
    residuals, V (unadjusted: the efficient variance; otherwise the sandwich with the same kind of
    meat), the weight matrix label and the number of clusters."""
    N = len(y)
    o = smp.o
    wm = str(o.get("wmatrix") or ("cluster " + smp.clustvar if smp.vce == "cluster" else "robust")).split()
    kind = wm[0].lower()
    center = bool(o.get("center"))
    e1 = y - X @ b1
    n_cl = None
    hac = None
    if kind == "hac":
        names = {"bartlett": "bartlett", "nwest": "bartlett", "parzen": "parzen", "gallant": "parzen",
                 "quadraticspectral": "quadraticspectral", "andrews": "quadraticspectral"}
        if len(wm) < 2 or wm[1].lower() not in names:
            raise StataError(198, "wmatrix(hac) requires a kernel")   # VERIFY
        kern = names[wm[1].lower()]
        if len(wm) > 2 and wm[2].lower() == "opt":
            raise StataError(198, "wmatrix(hac kernel opt) is not available yet in OpenDTA")   # OpenDTA
        lags = int(float(wm[2])) if len(wm) > 2 else N - 2          # VERIFY the default
        hac = (kern, lags, _hac_weights(kern, lags, N))

    def meat(e, kind_):
        nonlocal n_cl
        G = Z * e[:, None]
        if center and not kind_.startswith("una"):
            G = G - G.mean(axis=0)
        if kind_ == "hac":
            S = G.T @ G
            for j, wj in enumerate(hac[2], start=1):
                if wj != 0:
                    Gj = G[j:].T @ G[:-j]
                    S += wj * (Gj + Gj.T)
            return S
        if kind_.startswith("rob"):
            return G.T @ G
        if kind_.startswith("una"):
            return float(e @ e / N) * (Z.T @ Z)
        if kind_.startswith("clu"):
            from ..core.varlist import resolve_name
            cv = smp.cluster if smp.vce == "cluster" and len(wm) < 2 else \
                smp.s.data.get(resolve_name(smp.s.data, wm[1])).data[smp.mask]
            u = np.unique(cv)
            n_cl = len(u)
            T_ = np.array([G[cv == v].sum(axis=0) for v in u])
            if center:
                T_ = T_ - T_.mean(axis=0)
            return T_.T @ T_
        raise StataError(198, f"wmatrix({kind_}) is not available yet in OpenDTA")
    ZX = Z.T @ X
    W = np.linalg.inv(meat(e1, kind))
    b = None
    trace = []

    def mreldif(a_, b_):
        return float(np.max(np.abs(a_ - b_) / (np.abs(b_) + 1)))
    for _it in range(16000 if o.get("igmm") else 1):          # VERIFY igmm's limits
        Ainv = np.linalg.inv(ZX.T @ W @ ZX)
        b_new = Ainv @ (ZX.T @ W @ (Z.T @ y))
        if b is not None:
            W_new = W
            trace.append((mreldif(b_new, b), mreldif(W_new * N, W_old * N)))
        b_prev, b = b, b_new
        if not o.get("igmm"):
            break
        if trace and trace[-1][0] < 1e-6 and trace[-1][1] < 1e-6:      # VERIFY Stata's rule
            break
        W_old = W
        W = np.linalg.inv(meat(y - X @ b, kind))
    e = y - X @ b
    V = Ainv
    if smp.vce in ("robust", "cluster") or kind[:3] in ("rob", "clu", "hac"):
        # a robust or cluster weight matrix gives the sandwich variance (compat 1907)
        k2 = "hac" if kind == "hac" else ("cluster" if (smp.vce == "cluster" or kind.startswith("clu")) else "robust")
        S2 = meat(e, k2)
        V = Ainv @ (ZX.T @ W @ S2 @ W @ ZX) @ Ainv
    if hac is not None:
        label = f"HAC {hac[0].capitalize() if hac[0] != 'quadraticspectral' else 'quadratic spectral'} {hac[1]}"   # compat 1918
    else:
        label = {"rob": "Robust", "una": "Unadjusted", "clu": f"Cluster ({wm[1] if len(wm) > 1 else smp.clustvar})"}[kind[:3]]
    Ze = Z.T @ e
    J = float(Ze @ W @ Ze)              # with the weight matrix of the estimates (compat 1918)
    info = {"J": J, "S": meat(e, kind) / N, "W": W * N, "trace": trace, "iterations": max(len(trace), 1),
            "hac": hac}
    return b, e, (V + V.T) / 2, label, n_cl, info
