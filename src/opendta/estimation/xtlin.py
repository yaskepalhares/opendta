"""Linear panel commands: xtgls, xtpcse, xtregar, xtrc, xthtaylor and xtivreg ([XT] xtgls, [XT] xtpcse,
[XT] xtregar, [XT] xtrc, [XT] xthtaylor, [XT] xtivreg).

The estimators are in core/panellin.py. Layouts follow the manual's examples; where the manual
does not settle them the code carries a VERIFY mark and the case is in compat 1905."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core import panellin as pl
from ..core.errors import StataError
from ..core.formats import format_value
from ..core.stats import chi2_tail
from ..lazy import lazy
from . import fvars
from .regress import Sample, cmdline, comma
from .results import CoefRow, Estimates, coef_table, crit, post, rows_from_columns
from .tsops import panel_info, ts_info

st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def _g(x, fmt="%10.0g") -> str:
    if x is None or not np.isfinite(x) or abs(x) >= SYS:
        return "."
    return format_value(float(x), fmt, pad=False).strip()


class _Lin:
    """Sample of a linear panel command: y, the columns (with the constant last), and the panel
    index with the time position of every observation."""

    def __init__(self, s: "Session", cmd: str, args: str, opts: dict, need_time: bool = False,
                 need_dep: bool = True, extra_vars: str = "", constant: bool | None = None):
        self.s = s
        opts = dict(opts)
        opts.update({"i": 1, "t": 1})
        self.smp = smp = Sample(s, args, opts_spec=opts, weights=(), need_dep=need_dep, extra_vars=extra_vars)
        self.o = smp.o
        ds = s.data
        pvar = str(self.o["i"]).strip() if isinstance(self.o.get("i"), str) else panel_info(ds)
        if not pvar:
            raise StataError(459, "panel variable not set; use xtset varname ...")
        self.pvar = pvar
        ids = ds.get(pvar).data[smp.mask]
        uniq, inv = np.unique(ids, return_inverse=True)
        self.tvar = ""
        tinfo = ts_info(ds)
        t = None
        if isinstance(self.o.get("t"), str):
            self.tvar, delta = str(self.o["t"]).strip(), 1.0
            t = ds.get(self.tvar).data[smp.mask].astype(float)
        elif tinfo is not None and tinfo.tvar:
            self.tvar, delta = tinfo.tvar, tinfo.delta
            t = ds.get(self.tvar).data[smp.mask].astype(float)
        if t is None:
            if need_time:
                raise StataError(459, "time variable not set; use xtset panelvar timevar ...")
            self.tidx = np.zeros(len(inv), dtype=int)
            order = np.lexsort((np.arange(len(inv)), inv))
            pos = np.zeros(len(inv), dtype=int)
            cnt = np.bincount(inv)
            starts = np.concatenate([[0], np.cumsum(cnt)[:-1]])
            pos[order] = np.arange(len(inv)) - starts[inv[order]]
            self.tidx = pos
            self.has_time = False
        else:
            self.tidx = np.rint((t - t.min()) / delta).astype(int)
            self.has_time = True
        self.pi = pl.PanelIndex(inv, self.tidx)
        self.y = smp.y.astype(float) if smp.y is not None else None
        w = np.ones(smp.n)
        self.notes = fvars.drop_collinear(smp.cols, smp.constant if constant is None else constant, w)
        self.est_cols = [c for c in smp.cols if c.values is not None]
        self.constant = smp.constant if constant is None else constant
        Xs = np.column_stack([c.values for c in self.est_cols]) if self.est_cols else np.zeros((smp.n, 0))
        self.Xs = Xs
        self.X = np.column_stack([Xs, np.ones(smp.n)]) if self.constant else Xs
        self.k = self.X.shape[1]
        self.ntimes = int(len(np.unique(self.tidx)))

    def write_notes(self, out) -> None:
        for nm in self.notes:
            out.write(f"note: {nm} omitted because of collinearity\n", "text")

    def rows(self, b, V, extra=()):
        smp = self.smp
        bmap, semap = {}, {}
        j = 0
        for c in smp.cols:
            if c.values is not None:
                bmap[c.name] = float(b[j])
                semap[c.name] = float(np.sqrt(max(V[j, j], 0)))
                j += 1
            else:
                bmap[c.name], semap[c.name] = 0.0, 0.0
        if self.constant:
            bmap["_cons"], semap["_cons"] = float(b[-1]), float(np.sqrt(max(V[-1, -1], 0)))
        rows = rows_from_columns(smp.cols, bmap, semap, self.s.data, constant=self.constant)
        return rows + list(extra)

    def layout(self):
        smp = self.smp
        names = [c.name if not c.omitted else fvars.omitted_name(c.name) for c in smp.cols]
        free = [c.values is not None for c in smp.cols]
        if self.constant:
            names.append("_cons")
            free.append(True)
        return names, free

    def expand(self, b, V, extra_names=()):
        names, free = self.layout()
        kk = len(names) + len(extra_names)
        pos = [i for i, f in enumerate(free) if f] + list(range(len(names), kk))
        bf, Vf = np.zeros(kk), np.zeros((kk, kk))
        bf[pos] = b
        Vf[np.ix_(pos, pos)] = V
        return list(names) + ["_cons"] * len(extra_names), bf, Vf

    def wald(self, b, V):
        nslope = self.k - (1 if self.constant else 0)
        if nslope == 0:
            return SYS, 0, SYS
        idx = list(range(nslope))
        try:
            chi2 = float(b[idx] @ np.linalg.solve(V[np.ix_(idx, idx)], b[idx]))
        except np.linalg.LinAlgError:
            return SYS, nslope, SYS
        return chi2, nslope, float(chi2_tail(chi2, nslope))


def _opt_word(o, key, choices, default):
    v = o.get(key)
    if v in (None, False):
        return default
    if v is True:
        raise StataError(198, f"option {key}() requires a value")
    w = str(v).strip().lower()
    for c in choices:
        if c.startswith(w) and w:
            return c
    raise StataError(198, f"option {key}() incorrectly specified")


def _two_cols(out, rows_left, rows_right, width=48):
    for lft, (rl, rv) in zip(rows_left, rows_right):
        out.write(f"{lft:<{width}}{rl:<18}= " if rl else f"{lft:<{width}}", "text")
        if rl:
            out.write(f"{rv:>10}\n", "result")
        else:
            out.write("\n", "text")


# ---------------------------------------------------------------------------
# xtgls
# ---------------------------------------------------------------------------

def cmd_xtgls(s: "Session", args: str) -> None:
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, "xtgls", args)
        return
    L = _Lin(s, "xtgls", args, {"panels": 1, "corr": 4, "rhotype": 4, "igls": 4, "nmk": 3, "force": 5,
                                "iterate": 4, "tolerance": 3, "nolog": 4, "log": 3})
    o, out, pi = L.o, s.output, L.pi
    L.write_notes(out)
    panels = _opt_word(o, "panels", ("iid", "heteroskedastic", "correlated"), "iid")
    corr = _opt_word(o, "corr", ("independent", "ar1", "psar1"), "independent")
    rhot = _opt_word(o, "rhotype", pl.RHOTYPES, "regress")
    pkey = {"iid": "iid", "heteroskedastic": "hetero", "correlated": "correlated"}[panels]
    if corr != "independent" and not L.has_time:
        raise StataError(459, "time variable not set; use xtset panelvar timevar ...")
    y, X, k = L.y, L.X, L.k
    igls, nmk = bool(o.get("igls")), bool(o.get("nmk"))
    if panels == "correlated" and not L.has_time:
        raise StataError(459, "time variable not set; use xtset panelvar timevar ...")
    tol = float(o["tolerance"]) if o.get("tolerance") else 1e-7
    maxit = int(o["iterate"]) if o.get("iterate") else 100
    b_ols = np.linalg.lstsq(X, y, rcond=None)[0]
    e = y - X @ b_ols
    rho_g = None
    log: list[tuple[int, float]] = []

    def step(e_orig):
        nonlocal rho_g
        if corr == "independent":
            return y, X, e_orig, pl.gls_fit(y, X, pi, pkey, e_orig, nmk=nmk)
        rho_i = pl.rho_from_residuals(e_orig, pi, rhot, k)
        rho_g = rho_i if corr == "psar1" else np.full(pi.G, pl.common_rho(rho_i, pi))
        ro = rho_g[pi.inv]
        ys, Xs = pl.prais_winsten(y, pi, ro), pl.prais_winsten(X, pi, ro)
        bo = np.linalg.lstsq(Xs, ys, rcond=None)[0]
        es = ys - Xs @ bo
        return ys, Xs, es, pl.gls_fit(ys, Xs, pi, pkey, es, nmk=nmk)
    ys, Xs, es, fit = step(e)
    b = fit["b"]
    if igls:
        for it in range(1, maxit + 1):
            e_new = y - X @ b
            ys, Xs, es, fit = step(e_new)
            bn = fit["b"]
            tolv = float(np.max(np.abs(bn - b) / (np.abs(b) + 1e-10)))
            log.append((it, tolv))
            b = bn
            if tolv < tol:
                break
    ysf = ys - Xs @ b
    ll = pl.gls_loglik(ysf, pi, pkey, fit["sigma"], rho_g if corr != "independent" else None) if igls else None
    V = fit["V"]
    chi2, df, p = L.wald(b, V)
    if igls and not o.get("nolog"):
        for it, tv in log:
            out.write(f"Iteration {it}: tolerance = {_g(tv)}\n", "text")
        out.write("\n", "text") if log else None
    ncorr = 0 if corr == "independent" else (1 if corr == "ar1" else pi.G)
    pan_txt = {"iid": "homoskedastic", "heteroskedastic": "heteroskedastic",
               "correlated": "heteroskedastic with cross-sectional correlation"}[panels]
    if corr == "independent":
        cor_txt = "no autocorrelation"
    elif corr == "ar1":
        cor_txt = f"common AR(1) coefficient for all panels  ({_g(float(rho_g[0]), '%6.4f')})"
    else:
        cor_txt = "panel-specific AR(1)"
    N = len(y)
    x = dict(N=N, G=pi.G, T=L.ntimes, ncov=fit["ncov"], ncorr=ncorr, ncf=k, chi2=chi2, df=df, p=p, ll=ll,
             pan_txt=pan_txt, cor_txt=cor_txt)
    b_full_names, bf, Vf = L.expand(b, V)
    est = Estimates("xtgls", L.smp.depname, b_full_names, bf, Vf, N, "z", None)
    est.rows = L.rows(b, V)
    est.eqnames = [L.smp.depname] * len(b_full_names)
    est.cols, est.constant, est.terms, est.depcomp = L.smp.cols, L.constant, L.smp.terms, L.smp.depcomp
    est.level, est.vcetype = L.smp.level, ""
    est.extra = x
    est.display = _display_xtgls
    sc = [("N", float(N)), ("N_g", float(pi.G)), ("N_t", float(L.ntimes)), ("n_cf", float(k)), ("n_cv", float(fit["ncov"])),
          ("n_cr", float(ncorr)), ("df_m", float(df)), ("chi2", chi2), ("p", p), ("rank", float(np.linalg.matrix_rank(Vf)))]
    if ll is not None:
        sc.append(("ll", ll))
    if corr == "ar1":
        sc.append(("rho", float(rho_g[0])))
    mc = [("cmdline", cmdline("xtgls", args)), ("cmd", "xtgls"), ("depvar", L.smp.depname), ("ivar", L.pvar),
          ("tvar", L.tvar), ("rhotype", rhot if corr != "independent" else ""), ("corr", corr), ("panels", pkey),
          ("chi2type", "Wald"), ("vce", "oim"), ("title", "Cross-sectional time-series FGLS regression"),
          ("properties", "b V"), ("predict", "xtgls_p")]
    post(s, est, list(reversed(sc)), list(reversed(mc)), L.smp.mask)
    _display_xtgls(s, est)


def _display_xtgls(s, est, *, header=True, table=True, level=None, **_kw):
    out, x = s.output, est.extra
    if header:
        out.write("\nCross-sectional time-series FGLS regression\n\n", "text")
        out.write("Coefficients:  generalized least squares\n", "text")
        out.write(f"Panels:        {x['pan_txt']}\n", "text")
        out.write(f"Correlation:   {x['cor_txt']}\n\n", "text")
        left = [("Estimated covariances", comma(x["ncov"])), ("Estimated autocorrelations", comma(x["ncorr"])),
                ("Estimated coefficients", comma(x["ncf"])), ("", ""), ("", "")]
        right = [("Number of obs", comma(x["N"])), ("Number of groups", comma(x["G"])),
                 ("Time periods", comma(x["T"])), (f"Wald chi2({x['df']})", _g(x["chi2"], "%10.2f")),
                 ("Prob > chi2", _g(x["p"], "%10.4f"))]
        for i, ((ll_, lv), (rl, rv)) in enumerate(zip(left, right)):
            if x["ll"] is not None and rl == "Prob > chi2":
                out.write(f"{'Log likelihood':<27}= {_g(x['ll']):>9}{'':>11}", "text")
            elif ll_:
                out.write(f"{ll_:<27}= {lv:>9}{'':>11}", "text")
            else:
                out.write(f"{'':<48}", "text")
            out.write(f"{rl:<18}= ", "text")
            out.write(f"{rv:>10}\n", "result")
        out.write("\n", "text")
    if table:
        coef_table(s, est.depvar, est.rows, stat="z", df=None, level=level or est.level)


# ---------------------------------------------------------------------------
# xtpcse
# ---------------------------------------------------------------------------

def cmd_xtpcse(s: "Session", args: str) -> None:
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, "xtpcse", args)
        return
    L = _Lin(s, "xtpcse", args, {"correlation": 4, "rhotype": 4, "hetonly": 3, "independent": 3,
                                 "pairwise": 4, "np1": 3, "force": 5}, need_time=True)
    o, out, pi = L.o, s.output, L.pi
    L.write_notes(out)
    corr = _opt_word(o, "correlation", ("independent", "ar1", "psar1"), "independent")
    rhot = _opt_word(o, "rhotype", ("regress", "freg", "tscorr", "dw"), "regress")
    mode = "hetonly" if o.get("hetonly") else ("independent" if o.get("independent") else "correlated")
    y, X, k = L.y, L.X, L.k
    bo = np.linalg.lstsq(X, y, rcond=None)[0]
    e = y - X @ bo
    rho_g = None
    ys, Xs = y, X
    if corr != "independent":
        rho_i = pl.rho_from_residuals(e, pi, rhot, k)
        rho_g = rho_i if corr == "psar1" else np.full(pi.G, pl.common_rho(rho_i, pi, bool(o.get("np1"))))
        ro = rho_g[pi.inv]
        ys, Xs = pl.prais_winsten(y, pi, ro), pl.prais_winsten(X, pi, ro)
        bo = np.linalg.lstsq(Xs, ys, rcond=None)[0]
        e = ys - Xs @ bo
    V, ncov = pl.pcse_covariance(Xs, e, pi, mode, bool(o.get("pairwise")))
    V = (V + V.T) / 2
    b = bo
    chi2, df, p = L.wald(b, V)
    N = len(y)
    rss = float(e @ e)
    tss = float(((ys - ys.mean()) ** 2).sum()) if L.constant else float(ys @ ys)
    r2 = 1 - rss / tss if tss else SYS
    ncorr = 0 if corr == "independent" else (1 if corr == "ar1" else pi.G)
    bal = "balanced" if pi.balanced else "unbalanced"
    pan_txt = {"correlated": f"correlated ({bal})", "hetonly": f"heteroskedastic ({bal})",
               "independent": f"uncorrelated ({bal})"}[mode]
    cor_txt = {"independent": "no autocorrelation", "ar1": "common AR(1)", "psar1": "panel-specific AR(1)"}[corr]
    title = {"correlated": "Linear regression, correlated panels corrected standard errors (PCSEs)",
             "hetonly": "Linear regression, heteroskedastic panels corrected standard errors",
             "independent": "Linear regression, homoskedastic panels corrected standard errors"}[mode]
    x = dict(N=N, G=pi.G, tmin=int(pi.counts.min()), tavg=N / pi.G, tmax=int(pi.counts.max()), ncov=ncov, ncorr=ncorr,
             ncf=k, chi2=chi2, df=df, p=p, r2=r2, pan_txt=pan_txt, cor_txt=cor_txt, title=title, pvar=L.pvar, tvar=L.tvar)
    names, bf, Vf = L.expand(b, V)
    est = Estimates("xtpcse", L.smp.depname, names, bf, Vf, N, "z", None)
    est.rows = L.rows(b, V)
    est.eqnames = [L.smp.depname] * len(names)
    est.cols, est.constant, est.terms, est.depcomp = L.smp.cols, L.constant, L.smp.terms, L.smp.depcomp
    est.level, est.vcetype = L.smp.level, "Panel-corrected"
    est.extra = x
    est.display = _display_xtpcse
    sc = [("N", float(N)), ("N_g", float(pi.G)), ("N_t", float(L.ntimes)), ("n_cf", float(k)), ("n_cv", float(ncov)),
          ("n_cr", float(ncorr)), ("g_min", float(x["tmin"])), ("g_avg", x["tavg"]), ("g_max", float(x["tmax"])),
          ("df_m", float(df)), ("chi2", chi2), ("r2", r2), ("rss", rss), ("mss", tss - rss), ("rmse", math.sqrt(rss / N)),
          ("rank", float(np.linalg.matrix_rank(Vf)))]
    if rho_g is not None and corr == "ar1":
        sc.append(("rho", float(rho_g[0])))
    mc = [("cmdline", cmdline("xtpcse", args)), ("cmd", "xtpcse"), ("depvar", L.smp.depname), ("ivar", L.pvar),
          ("tvar", L.tvar), ("corr", corr), ("chi2type", "Wald"), ("vce", "independent" if mode == "independent" else mode),
          ("vcetype", "Panel-corrected"), ("title", title), ("properties", "b V"), ("predict", "xtpcse_p")]
    post(s, est, list(reversed(sc)), list(reversed(mc)), L.smp.mask)
    _display_xtpcse(s, est)


def _display_xtpcse(s, est, *, header=True, table=True, level=None, **_kw):
    out, x = s.output, est.extra
    if header:
        out.write(f"\n{x['title']}\n\n", "text")
        left = [f"{'Group variable:':<18}{x['pvar']}", f"{'Time variable:':<18}{x['tvar']}",
                f"{'Panels:':<18}{x['pan_txt']}", f"{'Autocorrelation:':<18}{x['cor_txt']}", "", ""]
        right = [("Number of obs", comma(x["N"])), ("Number of groups", comma(x["G"])), ("Obs per group:", ""),
                 ("            min", _g(x["tmin"])), ("            avg", _g(x["tavg"])),
                 ("            max", _g(x["tmax"]))]
        for lft, (rl, rv) in zip(left, right):
            out.write(f"{lft:<48}", "text")
            if rl.strip() in ("min", "avg", "max"):
                out.write(f"{'':>14}{rl.strip()} = ", "text")
                out.write(f"{rv:>10}\n", "result")
            elif rv == "":
                out.write(f"{rl}\n", "text")
            else:
                out.write(f"{rl:<18}= ", "text")
                out.write(f"{rv:>10}\n", "result")
        tail = [("Estimated covariances", comma(x["ncov"])), ("Estimated autocorrelations", comma(x["ncorr"])),
                ("Estimated coefficients", comma(x["ncf"]))]
        rtail = [("R-squared", _g(x["r2"], "%10.4f")), (f"Wald chi2({x['df']})", _g(x["chi2"], "%10.2f")),
                 ("Prob > chi2", _g(x["p"], "%10.4f"))]
        for (ll_, lv), (rl, rv) in zip(tail, rtail):
            out.write(f"{ll_:<27}= {lv:>9}{'':>11}{rl:<18}= ", "text")
            out.write(f"{rv:>10}\n", "result")
        out.write("\n", "text")
    if table:
        coef_table(s, est.depvar, est.rows, stat="z", df=None, level=level or est.level, vcetype=est.vcetype)


# ---------------------------------------------------------------------------
# xtregar
# ---------------------------------------------------------------------------

def _co(v, pi, rho):
    lagv = pi.lagged(v) if v.ndim == 1 else np.column_stack([pi.lagged(v[:, j]) for j in range(v.shape[1])])
    return (v - rho * lagv)


def cmd_xtregar(s: "Session", args: str) -> None:
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, "xtregar", args)
        return
    L = _Lin(s, "xtregar", args, {"fe": 2, "re": 2, "rhotype": 4, "rhof": 4, "twostep": 3, "lbi": 3},
             need_time=True, constant=True)
    o, out, pi = L.o, s.output, L.pi
    if o.get("lbi"):
        raise StataError(198, "option lbi is not available yet in OpenDTA")
    L.write_notes(out)
    fe = bool(o.get("fe"))
    rhot = _opt_word(o, "rhotype", ("dw", "regress", "freg", "tscorr", "theil", "nagar"), "dw")
    y, Xs = L.y, L.Xs
    n, k = len(y), Xs.shape[1]
    twostep = bool(o.get("twostep"))
    rho_fixed = float(o["rhof"]) if o.get("rhof") not in (None, False, True) else None
    keep = pi.has_lag
    G_all = pi.G

    def sub_index():
        inv = pi.inv[keep]
        _, inv2 = np.unique(inv, return_inverse=True)
        return pl.PanelIndex(inv2, pi.tidx[keep])
    pk = sub_index()
    Nk = int(keep.sum())

    def fit_given_rho(rho):
        ys = _co(y, pi, rho)[keep]
        Xk = _co(Xs, pi, rho)[keep]
        c = np.full(Nk, 1 - rho)
        if fe:
            yw, Xw = pk.within(ys), pk.within(Xk)
            A = np.linalg.inv(Xw.T @ Xw)
            b = A @ (Xw.T @ yw)
            e = yw - Xw @ b
            rss = float(e @ e)
            df = Nk - pk.G - k
            s2 = rss / df
            cons_t = float(ys.mean() - Xk.mean(axis=0) @ b)
            return dict(b=b, A=A, s2=s2, rss=rss, df=df, cons_t=cons_t, ys=ys, Xk=Xk, e=e, Xw=Xw, yw=yw)
        # random effects (Swamy-Arora variance components on the transformed data)
        yw, Xw = pk.within(ys), pk.within(Xk)
        bw = np.linalg.solve(Xw.T @ Xw, Xw.T @ yw)
        ew = yw - Xw @ bw
        s_e2 = float(ew @ ew) / (Nk - pk.G - k)
        yb = pk.group_means(ys)
        Xb = np.column_stack([pk.group_means(Xk), np.full(pk.G, 1 - rho)])
        bb = np.linalg.lstsq(Xb, yb, rcond=None)[0]
        eb = yb - Xb @ bb
        Tbar = pk.G / float(np.sum(1.0 / pk.counts))
        s_u2 = max(float(eb @ eb) / (pk.G - (k + 1)) - s_e2 / Tbar, 0.0)
        th = 1 - np.sqrt(s_e2 / (s_e2 + pk.counts * s_u2))
        thi = th[pk.inv]
        R = np.column_stack([Xk, c])
        Rs = R - thi[:, None] * np.column_stack([pk.group_means(Xk), np.full(pk.G, 1 - rho)])[pk.inv]
        ysr = ys - thi * yb[pk.inv]
        A = np.linalg.inv(Rs.T @ Rs)
        b = A @ (Rs.T @ ysr)
        e = ysr - Rs @ b
        return dict(b=b, A=A, s_e2=s_e2, s_u2=s_u2, theta=th, ys=ys, Xk=Xk, e=e, Rs=Rs, ysr=ysr)

    def rho_update(f, rho):
        # residuals of the original-scale model, deviated from their panel means
        if fe:
            b = f["b"]
            cons = f["cons_t"] / (1 - rho)
            r = y - Xs @ b - cons
        else:
            b = f["b"]
            r = y - Xs @ b[:k] - b[k] / (1 - rho) * 1.0
        r = pi.within(r)
        return _pooled_rho(r, pi, rhot, k)
    if rho_fixed is not None:
        rho = rho_fixed
        f = fit_given_rho(rho)
    else:
        rho = 0.0
        f = fit_given_rho(rho)
        for it in range(100):
            new = rho_update(f, rho)
            f = fit_given_rho(new)
            done = abs(new - rho) < 1e-9
            rho = new
            if twostep or done:
                break
    N = Nk
    G = pk.G
    z = crit("z", L.smp.level, None)
    if fe:
        b, A, s2 = f["b"], f["A"], f["s2"]
        V = s2 * A
        xbar = f["Xk"].mean(axis=0)
        cons_t = f["cons_t"]
        var_c = s2 / N + float(xbar @ V @ xbar)
        cons = cons_t / (1 - rho)
        se_c = math.sqrt(var_c) / (1 - rho)
        bb = np.append(b, cons)
        Vf = np.zeros((k + 1, k + 1))
        Vf[:k, :k] = V
        Vf[k, k] = se_c ** 2
        xb_t = f["Xk"] @ b
        u_i = pk.group_means(f["ys"]) - pk.group_means(xb_t) - cons_t
        sigma_u = float(np.std(u_i, ddof=1)) / (1 - rho) if G > 1 else 0.0
        sigma_e = math.sqrt(s2)
        corr = float(np.corrcoef(u_i[pk.inv], xb_t)[0, 1]) if k else 0.0
        F = float(b @ np.linalg.solve(V, b)) / k if k else SYS
        pF = float(st.f.sf(F, k, f["df"])) if F < SYS else SYS
        Xp = np.column_stack([f["Xk"], np.full(N, 1 - rho)])
        bp = np.linalg.lstsq(Xp, f["ys"], rcond=None)[0]
        rss_p = float(((f["ys"] - Xp @ bp) ** 2).sum())
        F_f = ((rss_p - f["rss"]) / (G - 1)) / (f["rss"] / f["df"])
        p_f = float(st.f.sf(F_f, G - 1, f["df"]))
        r2w = float(np.corrcoef(f["Xw"] @ b, f["yw"])[0, 1] ** 2) if k else 0.0
        r2b = float(np.corrcoef(pk.group_means(xb_t), pk.group_means(f["ys"]))[0, 1] ** 2) if k else 0.0
        r2o = float(np.corrcoef(xb_t, f["ys"])[0, 1] ** 2) if k else 0.0
        chi2 = pchi = None
        theta = None
        stat, df_r = "t", float(f["df"])
        title = "FE (within) regression with AR(1) disturbances"
    else:
        b_all, A = f["b"], f["A"]
        b = b_all[:k]
        cons = b_all[k] / (1 - rho)
        V_all = f["s_e2"] * A
        Vf = np.zeros((k + 1, k + 1))
        Vf[:k, :k] = V_all[:k, :k]
        Vf[k, k] = V_all[k, k] / (1 - rho) ** 2
        Vf[:k, k] = Vf[k, :k] = V_all[:k, k] / (1 - rho)
        V = Vf[:k, :k]
        bb = np.append(b, cons)
        sigma_u = math.sqrt(f["s_u2"]) / (1 - rho)
        sigma_e = math.sqrt(f["s_e2"])
        chi2 = float(b @ np.linalg.solve(V, b)) if k else SYS
        pchi = float(chi2_tail(chi2, k)) if k else SYS
        xb_t = f["Xk"] @ b
        r2w = float(np.corrcoef(pk.within(xb_t), pk.within(f["ys"]))[0, 1] ** 2) if k else 0.0
        r2b = float(np.corrcoef(pk.group_means(xb_t), pk.group_means(f["ys"]))[0, 1] ** 2) if k else 0.0
        r2o = float(np.corrcoef(xb_t, f["ys"])[0, 1] ** 2) if k else 0.0
        corr = F = pF = F_f = p_f = None
        theta = float(f["theta"][0]) if pk.balanced else None
        stat, df_r = "z", None
        title = "RE GLS regression with AR(1) disturbances"
    rho_fov = sigma_u ** 2 / (sigma_u ** 2 + sigma_e ** 2)
    names, bf, Vff = L.expand(bb, Vf)
    est = Estimates("xtregar", L.smp.depname, names, bf, Vff, N, stat, df_r)
    est.rows = L.rows(bb, Vf)
    est.eqnames = [L.smp.depname] * len(names)
    est.cols, est.constant, est.terms, est.depcomp = L.smp.cols, True, L.smp.terms, L.smp.depcomp
    est.level, est.vcetype = L.smp.level, ""
    est.model = "fe" if fe else "re"
    est.extra = dict(N=N, G=G, tmin=int(pk.counts.min()), tavg=N / G, tmax=int(pk.counts.max()), r2w=r2w, r2b=r2b, r2o=r2o,
                     k=k, df_r=f.get("df"), F=F, pF=pF, corr=corr, chi2=chi2, p=pchi, F_f=F_f, p_f=p_f, title=title,
                     pvar=L.pvar, rho=rho, sigma_u=sigma_u, sigma_e=sigma_e, rho_fov=rho_fov, theta=theta)
    est.display = _display_xtregar
    sc = [("N", float(N)), ("N_g", float(G)), ("df_m", float(k)), ("rho_ar", rho), ("sigma_u", sigma_u),
          ("sigma_e", sigma_e), ("rho_fov", rho_fov), ("r2_w", r2w), ("r2_b", r2b), ("r2_o", r2o),
          ("g_min", float(pk.counts.min())), ("g_avg", N / G), ("g_max", float(pk.counts.max())),
          ("rank", float(k + 1))]
    if fe:
        sc += [("df_r", float(f["df"])), ("F", F), ("F_f", F_f), ("corr", corr)]
    else:
        sc += [("chi2", chi2), ("theta", theta if theta is not None else SYS)]
    mc = [("cmdline", cmdline("xtregar", args)), ("cmd", "xtregar"), ("depvar", L.smp.depname), ("ivar", L.pvar),
          ("tvar", L.tvar), ("model", "fe" if fe else "re"), ("rhotype", rhot), ("title", title),
          ("properties", "b V"), ("predict", "xtregar_p")]
    post(s, est, list(reversed(sc)), list(reversed(mc)), L.smp.mask & True)
    _display_xtregar(s, est)


def _pooled_rho(r: np.ndarray, pi: pl.PanelIndex, kind: str, k: int) -> float:
    """One AR(1) coefficient from all panels together (pairs of consecutive observations)."""
    m = pi.has_lag
    a = np.zeros_like(r)
    a[m] = r[pi.lag[m]]
    b = np.where(m, r, 0.0)
    ab, a2, b2, e2 = float((a * b).sum()), float((a * a).sum()), float((b * b).sum()), float((r * r).sum())
    d2 = float((((b - a) ** 2) * m).sum())
    n = len(r)
    if kind == "regress":
        return ab / a2
    if kind == "freg":
        return ab / b2
    if kind == "tscorr":
        return ab / e2
    if kind == "dw":
        return 1 - d2 / e2 / 2
    if kind == "theil":
        return ab / e2 * (n - k) / n
    rdw = 1 - d2 / e2 / 2
    return (rdw * n ** 2 + k ** 2) / (n ** 2 - k ** 2)


def _display_xtregar(s, est, *, header=True, table=True, level=None, **_kw):
    out, x = s.output, est.extra
    if header:
        out.write("\n", "text")
        out.write(f"{x['title']:<48}{'Number of obs':<19}= ", "text")
        out.write(f"{comma(x['N']):>10}\n", "result")
        out.write(f"{'Group variable: ' + x['pvar']:<48}{'Number of groups':<19}= ", "text")
        out.write(f"{comma(x['G']):>10}\n\n", "result")
        rs = [("within ", x["r2w"], "min", _g(x["tmin"])), ("between", x["r2b"], "avg", _g(x["tavg"])),
              ("overall", x["r2o"], "max", _g(x["tmax"]))]
        for i, (nm, v, gl, gv) in enumerate(rs):
            lead = "R-sq:  " if i == 0 else "       "
            out.write(f"{lead}{nm} = ", "text")
            out.write(f"{v:.4f}", "result")
            tail = "Obs per group: " if i == 0 else "               "
            out.write(" " * (48 - 7 - 10 - 6) + f"{tail}{gl} = ", "text")
            out.write(f"{gv:>10}\n", "result")
        out.write("\n", "text")
        if est.model == "fe":
            out.write(f"{'':<48}{'F(' + str(x['k']) + ',' + str(int(x['df_r'])) + ')':<19}= ", "text")
            out.write(f"{_g(x['F'], '%10.2f'):>10}\n", "result")
            out.write(f"{'corr(u_i, Xb)  = ' + format(x['corr'], '.4f'):<48}{'Prob > F':<19}= ", "text")
            out.write(f"{_g(x['pF'], '%10.4f'):>10}\n\n", "result")
        else:
            out.write(f"{'':<48}{'Wald chi2(' + str(x['k']) + ')':<19}= ", "text")
            out.write(f"{_g(x['chi2'], '%10.2f'):>10}\n", "result")
            out.write(f"{'corr(u_i, Xb)  = 0 (assumed)':<48}{'Prob > chi2':<19}= ", "text")
            out.write(f"{_g(x['p'], '%10.4f'):>10}\n\n", "result")
    if not table:
        return
    coef_table(s, est.depvar, est.rows, stat=est.stat, df=est.df_r, level=level or est.level, footer=False)
    out.write("-" * 13 + "+" + "-" * 64 + "\n", "text")
    notes = {"rho_ar": "(estimated autocorrelation coefficient)",
             "rho_fov": "(fraction of variance because of u_i)" if est.model == "fe" else "(fraction of variance due to u_i)"}
    rows = [("rho_ar", x["rho"]), ("sigma_u", x["sigma_u"]), ("sigma_e", x["sigma_e"]), ("rho_fov", x["rho_fov"])]
    if x.get("theta") is not None:
        rows.append(("theta", x["theta"]))
    for nm, val in rows:
        out.write(f"{nm:>12} | ", "text")
        out.write(f"{_g(val):>10}", "result")
        if nm in notes and (nm == "rho_fov" or est.model == "re"):
            out.write("   " + notes[nm], "text")
        out.write("\n", "text")
    out.write("-" * 78 + "\n", "text")
    if est.model == "fe":
        left = f"F test that all u_i=0: F({int(x['G'] - 1)},{int(x['df_r'])}) = {x['F_f']:.2f}"
        out.write(left.ljust(61) + f"Prob > F = {x['p_f']:.4f}\n", "text")


# ---------------------------------------------------------------------------
# xtrc
# ---------------------------------------------------------------------------

def cmd_xtrc(s: "Session", args: str) -> None:
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, "xtrc", args)
        return
    L = _Lin(s, "xtrc", args, {"betas": 4, "offset": 3})
    o, out, pi = L.o, s.output, L.pi
    if o.get("betas"):
        raise StataError(198, "option betas is not available yet in OpenDTA")
    L.write_notes(out)
    y, X, k = L.y, L.X, L.k
    small = pi.counts <= k
    if small.any():
        keep = ~small[pi.inv]
        if not keep.any():
            raise StataError(2001, "insufficient observations")
        _, inv = np.unique(pi.inv[keep], return_inverse=True)
        pi2 = pl.PanelIndex(inv, pi.tidx[keep])
        yk, Xk = y[keep], X[keep]
    else:
        keep, pi2, yk, Xk = np.ones(len(y), bool), pi, y, X
    if pi2.G < 2:
        raise StataError(2001, "insufficient observations")
    r = pl.swamy_rc(yk, Xk, pi2)
    b, V = r["b"], r["V"]
    chi2, df, p = L.wald(b, V)
    N = len(yk)
    pc = float(chi2_tail(r["chi2"], r["df"]))
    x = dict(N=N, G=pi2.G, tmin=int(pi2.counts.min()), tavg=N / pi2.G, tmax=int(pi2.counts.max()), chi2=chi2, df=df, p=p,
             c_chi2=r["chi2"], c_df=r["df"], c_p=pc, pvar=L.pvar)
    names, bf, Vf = L.expand(b, V)
    est = Estimates("xtrc", L.smp.depname, names, bf, Vf, N, "z", None)
    est.rows = L.rows(b, V)
    est.eqnames = [L.smp.depname] * len(names)
    est.cols, est.constant, est.terms, est.depcomp = L.smp.cols, L.constant, L.smp.terms, L.smp.depcomp
    est.level, est.vcetype = L.smp.level, ""
    est.extra = x
    est.display = _display_xtrc
    sc = [("N", float(N)), ("N_g", float(pi2.G)), ("df_m", float(df)), ("chi2", chi2), ("chi2_c", r["chi2"]),
          ("g_min", float(x["tmin"])), ("g_avg", x["tavg"]), ("g_max", float(x["tmax"])), ("rank", float(np.linalg.matrix_rank(Vf)))]
    mc = [("cmdline", cmdline("xtrc", args)), ("cmd", "xtrc"), ("depvar", L.smp.depname), ("ivar", L.pvar),
          ("chi2type", "Wald"), ("title", "Random-coefficients regression"), ("properties", "b V"),
          ("predict", "xtrc_p")]
    post(s, est, list(reversed(sc)), list(reversed(mc)), L.smp.mask)
    _display_xtrc(s, est)


def _display_xtrc(s, est, *, header=True, table=True, level=None, **_kw):
    out, x = s.output, est.extra
    if header:
        out.write("\n", "text")
        out.write(f"{'Random-coefficients regression':<48}{'Number of obs':<19}= ", "text")
        out.write(f"{comma(x['N']):>10}\n", "result")
        out.write(f"{'Group variable: ' + x['pvar']:<48}{'Number of groups':<19}= ", "text")
        out.write(f"{comma(x['G']):>10}\n\n", "result")
        out.write(f"{'':<48}Obs per group:\n", "text")
        for lab, v in (("min", _g(x["tmin"])), ("avg", _g(x["tavg"], "%10.1f")), ("max", _g(x["tmax"]))):
            out.write(f"{'':<62}{lab} = ", "text")
            out.write(f"{v:>10}\n", "result")
        out.write("\n", "text")
        out.write(f"{'':<48}{'Wald chi2(' + str(x['df']) + ')':<19}= ", "text")
        out.write(f"{_g(x['chi2'], '%10.2f'):>10}\n", "result")
        out.write(f"{'':<48}{'Prob > chi2':<19}= ", "text")
        out.write(f"{_g(x['p'], '%10.4f'):>10}\n\n", "result")
    if table:
        coef_table(s, est.depvar, est.rows, stat="z", df=None, level=level or est.level)
        out.write(f"Test of parameter constancy:    chi2({x['c_df']}) = {x['c_chi2']:>8.2f}".ljust(58)
                  + f"Prob > chi2 = {x['c_p']:.4f}\n", "text")


# ---------------------------------------------------------------------------
# xtivreg
# ---------------------------------------------------------------------------

def _pairs_r2(xb, y, pi):
    def c2(a, b):
        if len(a) < 2 or np.std(a) == 0 or np.std(b) == 0:
            return 0.0
        return float(np.corrcoef(a, b)[0, 1] ** 2)
    return c2(pi.within(xb), pi.within(y)), c2(pi.group_means(xb), pi.group_means(y)), c2(xb, y)


def cmd_xtivreg(s: "Session", args: str) -> None:
    import re
    from .postest import replay
    t = args.strip()
    if not t or t.startswith(","):
        replay(s, "xtivreg", args)
        return
    mp = re.search(r"\(([^=()]+)=([^()]*)\)", t)
    if not mp:
        raise StataError(198, "equation not specified")   # VERIFY
    endog_txt, inst_txt = mp.group(1), mp.group(2)
    rest = t[:mp.start()] + " " + t[mp.end():]
    L = _Lin(s, "xtivreg", rest.strip(), {"re": 2, "be": 2, "fe": 2, "fd": 2, "ec2sls": 3, "nosa": 3, "small": 2,
                                          "regress": 3, "first": 3}, extra_vars=endog_txt + " " + inst_txt,
             constant=True)
    o, out, pi, ds = L.o, s.output, L.pi, s.data
    if o.get("first") or o.get("regress"):
        raise StataError(198, f"option {'first' if o.get('first') else 'regress'} is not available yet in OpenDTA")
    sel = [m for m in ("be", "fe", "fd") if o.get(m)]
    if len(sel) > 1 or (sel and o.get("re")):
        raise StataError(198, "only one of re, be, fe or fd may be specified")   # VERIFY
    if o.get("ec2sls"):
        method = "ec2sls"
    else:
        method = sel[0] if sel else "g2sls"
    if method == "fd" and not L.has_time:
        raise StataError(459, "time variable not set; use xtset panelvar timevar ...")
    L.write_notes(out)
    mask = L.smp.mask
    endog_terms = fvars.expand_fv(ds, endog_txt)
    inst_terms = fvars.expand_fv(ds, inst_txt)
    endog_cols = [c for c in fvars.build_design(ds, endog_terms, mask, constant=False) if c.values is not None]
    inst_cols = [c for c in fvars.build_design(ds, inst_terms, mask, constant=False) if c.values is not None]
    ex_cols = L.est_cols
    Xn = endog_cols
    X = np.column_stack([c.values for c in Xn + ex_cols])
    Z = np.column_stack([c.values for c in inst_cols + ex_cols])
    k = X.shape[1]
    if Z.shape[1] < k:
        raise StataError(481, "equation not identified; must have at least as many instruments not in\n"
                         "the regression as there are instrumented variables")
    y = L.y
    r = pl.panel_iv(y, X, Z, pi, method, bool(o.get("nosa")))
    small = bool(o.get("small"))
    has_cons = method != "fd"
    b = r["b"]
    s2 = r["s2"]
    V = s2 * r["A"]
    if method == "fe":
        V_full = np.zeros((k + 1, k + 1))
        V_full[:k, :k] = V
        xbar = r["xbar"]
        V_full[k, k] = s2 / r["n"] + float(xbar @ V @ xbar)
        V_full[:k, k] = V_full[k, :k] = -V @ xbar
        V = V_full
    names_all = [c.name for c in Xn + ex_cols] + (["_cons"] if has_cons else [])
    nobs = r["n"]
    G = pi.G
    Gn = G
    if method == "fd":
        rows_i = r["rows"]
        Gn = int(len(np.unique(pi.inv[rows_i])))
    sl = list(range(k))
    W = float(b[sl] @ np.linalg.solve(V[np.ix_(sl, sl)], b[sl])) if k else SYS
    stat = "t" if small else "z"
    df_r = float(r["df"]) if small else None
    # VERIFY: Wald chi2 and z without small, F and t with small, for every estimator
    use_F = small
    Fv = W / k
    if use_F:
        pv = float(st.f.sf(Fv, k, r["df"]))
    else:
        pv = float(chi2_tail(W, k))
    xb_all = X @ b[:k]
    r2w, r2b, r2o = _pairs_r2(xb_all, y, pi)
    x = dict(N=nobs, G=Gn, tmin=int(pi.counts.min()), tavg=len(y) / G, tmax=int(pi.counts.max()), k=k, W=W, F=Fv, p=pv,
             use_F=use_F, df_r=r["df"], r2w=r2w, r2b=r2b, r2o=r2o, method=method, pvar=L.pvar,
             instd=" ".join(c.name for c in Xn), insts=" ".join([c.name for c in ex_cols] + [c.name for c in inst_cols]))
    if method in ("g2sls", "ec2sls"):
        su, se_ = math.sqrt(r["sigma_u2"]), math.sqrt(r["sigma_e2"])
        x.update(sigma_u=su, sigma_e=se_, rho=su ** 2 / (su ** 2 + se_ ** 2), theta=float(r["theta"][0]) if pi.balanced else None)
    elif method == "fe":
        u_i = pi.group_means(y) - pi.group_means(xb_all) - r["cons"]
        su = float(np.std(u_i, ddof=1)) if G > 1 else 0.0
        se_ = math.sqrt(s2)
        corr = float(np.corrcoef(u_i[pi.inv], xb_all)[0, 1]) if k else 0.0
        x.update(sigma_u=su, sigma_e=se_, rho=su ** 2 / (su ** 2 + se_ ** 2), corr=corr)
        rss_w = float(r["e"] @ r["e"])
        # F test that all u_i = 0: pooled 2SLS against the within fit (VERIFY)
        Xp = np.column_stack([X, np.ones(len(y))])
        Zp = np.column_stack([Z, np.ones(len(y))])
        fp = pl.tsls(y, Xp, Zp)
        rss_p = float(fp["e"] @ fp["e"])
        Ff = ((rss_p - rss_w) / (G - 1)) / (rss_w / r["df"])
        x.update(F_f=Ff, p_f=float(st.f.sf(Ff, G - 1, r["df"])))
    # coefficient order: endogenous first, then exogenous (as in ivregress)
    bmap = {nm: float(b[i]) for i, nm in enumerate(names_all)}
    semap = {nm: float(math.sqrt(max(V[i, i], 0))) for i, nm in enumerate(names_all)}
    rows = rows_from_columns(Xn + ex_cols, bmap, semap, ds, constant=has_cons)
    est = Estimates("xtivreg", L.smp.depname, names_all, b, V, nobs, stat, df_r)
    est.rows = rows
    est.eqnames = [L.smp.depname] * len(names_all)
    est.cols, est.constant, est.terms, est.depcomp = Xn + ex_cols, has_cons, endog_terms + L.smp.terms, L.smp.depcomp
    est.level, est.vcetype = L.smp.level, ""
    est.extra = x
    est.display = _display_xtivreg
    ttl = {"g2sls": "G2SLS random-effects IV regression", "ec2sls": "EC2SLS random-effects IV regression",
           "fe": "Fixed-effects (within) IV regression", "be": "Between regression (regressing means)",
           "fd": "First-differenced IV regression"}[method]
    x["title"] = ttl
    sc = [("N", float(nobs)), ("N_g", float(Gn)), ("df_m", float(k)), ("g_min", float(x["tmin"])),
          ("g_avg", x["tavg"]), ("g_max", float(x["tmax"])), ("r2_w", r2w), ("r2_b", r2b), ("r2_o", r2o),
          ("rank", float(k + (1 if has_cons else 0)))]
    sc += [("F", Fv)] if use_F else [("chi2", W)]
    if "sigma_u" in x:
        sc += [("sigma_u", x["sigma_u"]), ("sigma_e", x["sigma_e"]), ("rho", x["rho"])]
    mc = [("cmdline", cmdline("xtivreg", args)), ("cmd", "xtivreg"), ("depvar", L.smp.depname), ("ivar", L.pvar),
          ("model", method), ("instd", x["instd"]), ("insts", x["insts"]), ("title", ttl), ("properties", "b V"),
          ("predict", "xtivreg_p")]
    post(s, est, list(reversed(sc)), list(reversed(mc)), mask)
    _display_xtivreg(s, est)


def _display_xtivreg(s, est, *, header=True, table=True, level=None, **_kw):
    out, x = s.output, est.extra
    if header:
        out.write("\n", "text")
        out.write(f"{x['title']:<48}{'Number of obs':<19}= ", "text")
        out.write(f"{comma(x['N']):>10}\n", "result")
        out.write(f"{'Group variable: ' + x['pvar']:<48}{'Number of groups':<19}= ", "text")
        out.write(f"{comma(x['G']):>10}\n\n", "result")
        rs = [("within ", x["r2w"], "min", _g(x["tmin"])), ("between", x["r2b"], "avg", _g(x["tavg"], "%10.1f")),
              ("overall", x["r2o"], "max", _g(x["tmax"]))]
        for i, (nm, v, gl, gv) in enumerate(rs):
            lead = "R-sq:  " if i == 0 else "       "
            out.write(f"{lead}{nm} = ", "text")
            out.write(f"{v:.4f}", "result")
            tail = "Obs per group: " if i == 0 else "               "
            out.write(" " * (48 - 7 - 10 - 6) + f"{tail}{gl} = ", "text")
            out.write(f"{gv:>10}\n", "result")
        out.write("\n", "text")
        if x["use_F"]:
            lab = f"F({x['k']},{int(x['df_r'])})"
            out.write(f"{'':<48}{lab:<19}= ", "text")
            out.write(f"{_g(x['F'], '%10.2f'):>10}\n", "result")
            cl = f"corr(u_i, Xb)  = {x['corr']:.4f}" if "corr" in x else ""
            out.write(f"{cl:<48}{'Prob > F':<19}= ", "text")
            out.write(f"{_g(x['p'], '%10.4f'):>10}\n\n", "result")
        else:
            out.write(f"{'':<48}{'Wald chi2(' + str(x['k']) + ')':<19}= ", "text")
            out.write(f"{_g(x['W'], '%10.2f'):>10}\n", "result")
            cl = f"corr(u_i, Xb)  = {x['corr']:.4f}" if "corr" in x else "corr(u_i, X)       = 0 (assumed)"
            out.write(f"{cl if x['method'] in ('fe', 'g2sls', 'ec2sls') else '':<48}{'Prob > chi2':<19}= ", "text")
            out.write(f"{_g(x['p'], '%10.4f'):>10}\n\n", "result")
    if not table:
        return
    coef_table(s, est.depvar, est.rows, stat=est.stat, df=est.df_r, level=level or est.level, footer=False)
    if "sigma_u" in x:
        out.write("-" * 13 + "+" + "-" * 64 + "\n", "text")
        for nm in ("sigma_u", "sigma_e", "rho"):
            out.write(f"{nm:>12} | ", "text")
            out.write(f"{_g(x[nm]):>10}", "result")
            if nm == "rho":
                out.write("   (fraction of variance due to u_i)", "text")
            out.write("\n", "text")
    out.write("-" * 78 + "\n", "text")
    if x["method"] == "fe":
        left = f"F test that all u_i=0:     F({int(x['G'] - 1)}, {int(x['df_r'])}) = {x['F_f']:.2f}"
        out.write(left.ljust(58) + f"Prob > F = {x['p_f']:.4f}\n", "text")
    out.write(f"Instrumented:   {x['instd']}\n", "text")
    out.write(f"Instruments:    {x['insts']}\n", "text")
    out.write("-" * 78 + "\n", "text")


# ---------------------------------------------------------------------------
# xthtaylor
# ---------------------------------------------------------------------------

def cmd_xthtaylor(s: "Session", args: str) -> None:
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, "xthtaylor", args)
        return
    L = _Lin(s, "xthtaylor", args, {"endog": 3, "constant": 5, "varying": 4, "amacurdy": 2, "small": 2},
             constant=True)
    o, out, pi, ds = L.o, s.output, L.pi, s.data
    if not o.get("endog"):
        raise StataError(198, "option endog() required")   # VERIFY
    L.write_notes(out)
    cols = L.est_cols
    names = [c.name for c in cols]

    def named(key):
        txt = str(o.get(key) or "").strip()
        if not txt:
            return None
        got = [c.name for c in fvars.build_design(ds, fvars.expand_fv(ds, txt), L.smp.mask, constant=False)
               if c.values is not None]
        bad = [g for g in got if g not in names]
        if bad:
            raise StataError(198, f"{bad[0]} not found among the independent variables")   # VERIFY
        return set(got)
    endog = named("endog")
    const = named("constant")
    vary = named("varying")
    allnames = set(names)
    if const is None and vary is None:
        const = {c.name for c in cols if np.allclose(pi.within(c.values), 0)}
        vary = allnames - const
    elif const is None:
        const = allnames - vary
    elif vary is None:
        vary = allnames - const
    if const & vary:
        raise StataError(198, "a variable cannot be both constant and varying")   # VERIFY
    cl = {c.name: c for c in cols}
    grp = {"TVexogenous": [n for n in names if n in vary and n not in endog],
           "TVendogenous": [n for n in names if n in vary and n in endog],
           "TIexogenous": [n for n in names if n in const and n not in endog],
           "TIendogenous": [n for n in names if n in const and n in endog]}
    order = grp["TVexogenous"] + grp["TVendogenous"] + grp["TIexogenous"] + grp["TIendogenous"]

    def mat(ns):
        return np.column_stack([cl[n].values for n in ns]) if ns else np.zeros((len(L.y), 0))
    X1, X2, Z1, Z2 = (mat(grp[k]) for k in ("TVexogenous", "TVendogenous", "TIexogenous", "TIendogenous"))
    if X1.shape[1] < Z2.shape[1]:
        raise StataError(481, "equation not identified; at least as many time-varying exogenous variables "
                         "as time-invariant endogenous variables are required")   # VERIFY
    r = pl.hausman_taylor(L.y, X1, X2, Z1, Z2, pi, bool(o.get("amacurdy")))
    b, V = r["b"], r["V"] * 1.0
    k = len(order)
    N, G = len(L.y), pi.G
    W = float(b[:k] @ np.linalg.solve(V[:k, :k], b[:k]))
    p = float(chi2_tail(W, k))
    su, se_ = r["sigma_u"], r["sigma_e"]
    rho = su ** 2 / (su ** 2 + se_ ** 2)
    bmap = {nm: float(b[i]) for i, nm in enumerate(order + ["_cons"])}
    semap = {nm: float(math.sqrt(max(V[i, i], 0))) for i, nm in enumerate(order + ["_cons"])}
    ordered_cols = [cl[n] for n in order]
    rows = rows_from_columns(ordered_cols, bmap, semap, ds, constant=True)
    # block headers
    pos = 0
    out_rows = []
    for gname in ("TVexogenous", "TVendogenous", "TIexogenous", "TIendogenous"):
        n_g = len(grp[gname])
        if n_g:
            out_rows.append(CoefRow("header", gname))
            out_rows.extend(rows[pos:pos + n_g])
            pos += n_g
    out_rows.extend(rows[pos:])
    est = Estimates("xthtaylor", L.smp.depname, order + ["_cons"], b, V, N, "z", None)
    est.rows = out_rows
    est.eqnames = [L.smp.depname] * (k + 1)
    est.cols, est.constant, est.terms, est.depcomp = ordered_cols, True, L.smp.terms, L.smp.depcomp
    est.level, est.vcetype = L.smp.level, ""
    est.extra = dict(N=N, G=G, tmin=int(pi.counts.min()), tavg=N / G, tmax=int(pi.counts.max()), W=W, k=k, p=p,
                     sigma_u=su, sigma_e=se_, rho=rho, pvar=L.pvar)
    est.display = _display_xthtaylor
    sc = [("N", float(N)), ("N_g", float(G)), ("df_m", float(k)), ("chi2", W), ("sigma_u", su), ("sigma_e", se_),
          ("rho", rho), ("g_min", float(est.extra["tmin"])), ("g_avg", est.extra["tavg"]),
          ("g_max", float(est.extra["tmax"])), ("rank", float(k + 1))]
    mc = [("cmdline", cmdline("xthtaylor", args)), ("cmd", "xthtaylor"), ("depvar", L.smp.depname), ("ivar", L.pvar),
          ("model", "amacurdy" if o.get("amacurdy") else "ht"), ("title", "Hausman-Taylor estimation"),
          ("properties", "b V"), ("predict", "xthtaylor_p")]
    post(s, est, list(reversed(sc)), list(reversed(mc)), L.smp.mask)
    _display_xthtaylor(s, est)


def _display_xthtaylor(s, est, *, header=True, table=True, level=None, **_kw):
    out, x = s.output, est.extra
    if header:
        out.write("\n", "text")
        out.write(f"{'Hausman-Taylor estimation':<48}{'Number of obs':<19}= ", "text")
        out.write(f"{comma(x['N']):>10}\n", "result")
        out.write(f"{'Group variable: ' + x['pvar']:<48}{'Number of groups':<19}= ", "text")
        out.write(f"{comma(x['G']):>10}\n\n", "result")
        out.write(f"{'':<48}Obs per group: min = ", "text")
        out.write(f"{_g(x['tmin']):>10}\n", "result")
        out.write(f"{'':<63}avg = ", "text")
        out.write(f"{_g(x['tavg'], '%10.1f'):>10}\n", "result")
        out.write(f"{'':<63}max = ", "text")
        out.write(f"{_g(x['tmax']):>10}\n\n", "result")
        out.write(f"{'Random effects u_i ~ i.i.d.':<48}{'Wald chi2(' + str(x['k']) + ')':<19}= ", "text")
        out.write(f"{_g(x['W'], '%10.2f'):>10}\n", "result")
        out.write(f"{'':<48}{'Prob > chi2':<19}= ", "text")
        out.write(f"{_g(x['p'], '%10.4f'):>10}\n\n", "result")
    if not table:
        return
    coef_table(s, est.depvar, est.rows, stat="z", df=None, level=level or est.level, footer=False)
    out.write("-" * 13 + "+" + "-" * 64 + "\n", "text")
    for nm in ("sigma_u", "sigma_e", "rho"):
        out.write(f"{nm:>12} | ", "text")
        out.write(f"{_g(x[nm]):>10}", "result")
        if nm == "rho":
            out.write("   (fraction of variance due to u_i)", "text")
        out.write("\n", "text")
    out.write("-" * 78 + "\n", "text")
    out.write("note: TV refers to time varying; TI refers to time invariant.\n", "text")
