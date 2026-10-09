"""Painel e efeitos absorvidos: xtreg (fe, re, be), areg ([XT] xtreg, [R] areg),
e variáveis instrumentais: ivregress 2sls ([R] ivregress)."""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

import numpy as np
from scipy import stats as st

from ..core import missing as M
from ..core.errors import StataError
from ..core.formats import format_value
from ..lang.syntax import match_options
from . import fvars
from .regress import Sample, cmdline, comma, sig5, wexp_macros
from .results import CoefRow, Estimates, coef_table, post, rows_from_columns
from .tsops import panel_info

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
    smp = Sample(s, args, opts_spec={"fe": 2, "re": 2, "be": 2, "mle": 3, "i": 1, "sa": 2},
                 weights=("aweight", "fweight", "pweight"))
    o = smp.o
    ds = s.data
    pvar = str(o.get("i")).strip() if o.get("i") else panel_info(ds)
    if not pvar:
        raise StataError(459, "panel variable not set; use xtset varname ...")   # VERIFICAR
    ids_all = ds.get(pvar).data
    mask = smp.mask & (ids_all < SYS)
    if not np.array_equal(mask, smp.mask):
        raise StataError(459, f"{pvar} missing in sample")   # VERIFICAR
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
        raise StataError(198, "xtreg, be not yet supported by OpenDTA")   # VERIFICAR


def _drop_within_collinear(smp, Xw, est_cols, out):
    """Variáveis constantes dentro dos grupos somem no fe (nota de colinearidade)."""
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
    # V com a constante: regressão de (y - ȳi + ȳ) em (x - x̄i + x̄) e 1
    Xa = np.column_stack([Xw + xbar, np.ones(N)])
    Aa = np.linalg.inv(Xa.T @ Xa)
    V = s2 * Aa
    bb = np.append(b, cons)
    vcetype = ""
    n_clust = None
    if smp.vce in ("robust", "cluster"):
        # xtreg, vce(robust) = cluster no painel (Stata)
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
    # F de que todos os u_i = 0: MQO agrupado vs fe
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
    scalars = [("N", float(N)), ("N_g", float(G)), ("df_m", float(k)), ("sigma", sigma_e), ("sigma_u", sigma_u),
               ("sigma_e", sigma_e), ("r2_w", r2_w), ("r2_o", r2_o), ("r2_b", r2_b), ("corr", corr),
               ("rho", rho), ("F", F), ("F_f", F_f), ("df_a", float(G - 1)), ("df_b", float(k)),
               ("df_r", float(df_r)), ("Tbar", N / G), ("Tcon", 1.0 if counts.min() == counts.max() else 0.0),
               ("rss", rss), ("g_min", float(counts.min())), ("g_avg", N / G), ("g_max", float(counts.max())),
               ("rank", float(k + 1))]
    if n_clust is not None:
        scalars.append(("N_clust", float(n_clust)))
    macros = [("cmdline", cmdline("xtreg", args)), ("cmd", "xtreg"), ("depvar", smp.depname),
              ("ivar", pvar), ("model", "fe"), ("title", "Fixed-effects (within) regression"),
              ("vce", smp.vce if vcetype else "conventional"), ("properties", "b V"),
              ("predict", "xtrefe_p"), ("marginsok", "XB XBU U")]
    if vcetype:
        macros += [("vcetype", vcetype), ("clustvar", est.clustvar)]
    post(s, est, scalars, macros, smp.mask)
    display_xtreg(s, est)


def _xtreg_re(s, smp, args, y, Xs, est_cols, ids, inv, counts, G, N, k, pvar):
    # Swamy-Arora: σe² do within, σu² do between (T̄ harmônico)
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
    V = sig_e2 * A
    vcetype = ""
    n_clust = None
    if smp.vce in ("robust", "cluster"):
        cl = smp.cluster if smp.vce == "cluster" else ids
        V, n_clust = _cluster_vce(Xstar, e, cl, A, N, k + 1, g_count_adjust=False)
        vcetype = "Robust"
    V = (V + V.T) / 2
    chi2 = float(b[:k] @ np.linalg.solve(V[:k, :k], b[:k])) if k else 0.0
    p = float(st.chi2.sf(chi2, k)) if k else SYS
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
        title = "Fixed-effects (within) regression" if est.model == "fe" else "Random-effects GLS regression"
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
            out.write(" " * (48 - 5 - 7 - 3 - 6) + f"{gl:>18}{'':<0} = ", "text")
            out.write(f"{gv:>10}\n", "result")
        out.write("\n", "text")
        if est.model == "fe":
            fl = f"F({int(x['df_m'])},{int(x['df_r'])})"
            out.write(f"{'':<48}{fl:<18}= ", "text")
            out.write(f"{g(x['F'], '%10.2f'):>10}\n", "result")
            out.write(f"{'corr(u_i, Xb)  = ' + format(x['corr'], '.4f'):<48}{'Prob > F':<18}= ", "text")
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
    coef_table(s, est.depvar, est.rows, stat=est.stat, df=est.df_r, level=level or est.level,
               vcetype=est.vcetype, footer=False)
    out.write("-" * 13 + "+" + "-" * 64 + "\n", "text")
    for nm, val in (("sigma_u", x["sigma_u"]), ("sigma_e", x["sigma_e"]), ("rho", x["rho"])):
        out.write(f"{nm:>12} |  ", "text")
        txt = format_value(val, "%9.0g", pad=False).strip() if val == val else "."
        out.write(f"{txt:>9}", "result")
        if nm == "rho":
            out.write("   (fraction of variance due to u_i)", "text")
        out.write("\n", "text")
    out.write("-" * 78 + "\n", "text")
    if est.model == "fe" and not est.vcetype:
        out.write(f"F test that all u_i=0: F({int(x['G'] - 1)}, {int(x['df_r'])}) = {x['F_f']:.2f}".ljust(61)
                  + f"Prob > F = {x['p_f']:.4f}\n", "text")   # VERIFICAR alinhamento


def _predict_xt(s, est, stat, mask):
    from .postest import xb_all
    xb, ok = xb_all(s, est)
    if stat in ("", "xb"):
        note = "(option xb assumed; fitted values)" if not stat else ""
        return np.where(ok, xb, SYS), note
    raise StataError(198, f"option {stat} not yet supported by OpenDTA")


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
        smp.mask &= av < SYS   # VERIFICAR
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
                 ("Adj R-squared", g(x["r2_a"], "%10.4f")), ("Root MSE", sig5(x["rmse"]))]
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
    out.write(f"   F({int(x['G'] - 1)}, {int(x['N'] - x['G'] - x['df_m'])}) = {x['F_abs']:9.3f}   {x['p_abs']:.3f}"
              f"          ({int(x['G'])} categories)\n", "result")   # VERIFICAR layout
    out.write("-" * 78 + "\n", "text")


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
        raise StataError(198, "estimator must be 2sls, liml, or gmm")   # VERIFICAR
    estimator = m.group(1)
    if estimator != "2sls":
        raise StataError(198, f"ivregress {estimator} not yet supported by OpenDTA")
    rest = m.group(2)
    mp = re.search(r"\(([^=()]+)=([^()]*)\)", rest)
    if not mp:
        raise StataError(198, "equation not specified")   # VERIFICAR
    endog_txt, inst_txt = mp.group(1), mp.group(2)
    exog_txt = rest[:mp.start()] + " " + rest[mp.end():]
    # amostra: dependente + exógenas + endógenas + instrumentos sem missing
    smp = Sample(s, exog_txt.strip(), opts_spec={"small": 2, "first": 3},
                 weights=("aweight", "fweight", "pweight"), extra_vars=endog_txt + " " + inst_txt)
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
    PzX = Z @ np.linalg.lstsq(Z, X, rcond=None)[0]
    A = np.linalg.inv(PzX.T @ X)
    b = A @ (PzX.T @ y)
    e = y - X @ b
    rss = float(e @ e)
    small = bool(smp.o.get("small"))
    s2 = rss / (N - k) if small else rss / N
    V = s2 * A
    vcetype = ""
    n_clust = None
    if smp.vce in ("robust", "cluster"):
        sc = PzX * e[:, None]
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
    ybar = y.mean()
    tss = float(np.sum((y - ybar) ** 2))
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
    est.extra = dict(N=N, W=W, df_m=len(slopes), r2=r2, rmse=rmse, small=small,
                     F=W / len(slopes), instd=" ".join(c.name for c in Xn),
                     insts=" ".join([c.name for c in ex_cols] + [c.name for c in Zc]))
    est.display = display_iv
    est.predict = _predict_xt
    scalars = [("N", float(N)), ("mss", tss - rss), ("df_m", float(len(slopes))), ("rss", rss), ("r2", r2),
               ("rmse", rmse), ("rank", float(k))]
    scalars += [("F", W / len(slopes)), ("df_r", float(N - k))] if small else [("chi2", W)]
    macros = [("cmdline", cmdline("ivregress", args)), ("cmd", "ivregress"), ("estimator", "2sls"),
              ("depvar", smp.depname), ("instd", est.extra["instd"]), ("insts", est.extra["insts"]),
              ("title", "Instrumental variables (2SLS) regression"), ("vce", smp.vce if vcetype else "unadjusted"),
              ("properties", "b V"), ("predict", "ivreg_p")]
    post(s, est, scalars, macros, mask)
    display_iv(s, est)


def _retarget(est, Xn, ex_cols, endog_terms, ex_terms):
    """As colunas das endógenas e exógenas usam índices de termo da lista unida."""
    off = len(endog_terms)
    for c in ex_cols:
        c.term = c.term + off if not getattr(c, "_retargeted", False) else c.term
        c._retargeted = True


def display_iv(s, est, *, header=True, table=True, level=None, **_kw):
    out = s.output
    x = est.extra
    if header:
        out.write("\n", "text")
        if x["small"]:
            right = [("Number of obs", comma(x["N"])), (f"F({int(x['df_m'])}, {int(est.df_r)})", g(x["F"], "%10.2f")),
                     ("Prob > F", g(st.f.sf(x["F"], x["df_m"], est.df_r), "%10.4f"))]
        else:
            right = [("Number of obs", comma(x["N"])), (f"Wald chi2({int(x['df_m'])})", g(x["W"], "%10.2f")),
                     ("Prob > chi2", g(st.chi2.sf(x["W"], x["df_m"]), "%10.4f"))]
        right += [("R-squared", g(x["r2"], "%10.4f")), ("Root MSE", sig5(x["rmse"]))]
        for i, (rl, rv) in enumerate(right):
            lft = "Instrumental variables (2SLS) regression" if i == 0 else ""
            out.write(f"{lft:<50}{rl:<16}= ", "text")
            out.write(f"{rv:>10}\n", "result")
        out.write("\n", "text")
    if not table:
        return
    coef_table(s, est.depvar, est.rows, stat=est.stat, df=est.df_r, level=level or est.level,
               vcetype=est.vcetype)
    out.write(f"Instrumented:  {x['instd']}\n", "text")
    out.write(f"Instruments:   {x['insts']}\n", "text")
