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
    # Stata prints the log likelihood for the OLS fit (iid, no autocorrelation) and for igls without autocorrelation
    want_ll = corr == "independent" and (panels == "iid" or igls)
    ll = pl.gls_loglik(ysf, pi, pkey, fit["sigma"], None) if want_ll else None
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
             pan_txt=pan_txt, cor_txt=cor_txt, singular=bool(fit.get("singular")))
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
        left = [("Estimated covariances", str(x["ncov"])), ("Estimated autocorrelations", str(x["ncorr"])),
                ("Estimated coefficients", str(x["ncf"])), ("", ""), ("", "")]
        right = [("Number of obs", comma(x["N"])), ("Number of groups", comma(x["G"])),
                 ("Time periods", comma(x["T"])), (f"Wald chi2({x['df']})", _g(x["chi2"], "%10.2f")),
                 ("Prob > chi2", _g(x["p"], "%10.4f"))]
        for i, ((ll_, lv), (rl, rv)) in enumerate(zip(left, right)):
            if x["ll"] is not None and rl == "Prob > chi2":
                out.write(f"{'Log likelihood':<27}= {_g(x['ll'], '%9.0g'):>9}{'':>10}", "text")
            elif ll_:
                out.write(f"{ll_:<27}= {lv:>9}{'':>10}", "text")
            else:
                out.write(f"{'':<48}", "text")
            out.write(f"{rl:<18}= ", "text")
            out.write(f"{rv:>10}\n", "result")
        out.write("\n", "text")
    if table:
        coef_table(s, est.depvar, est.rows, stat="z", df=None, level=level or est.level)
        if x.get("singular"):
            out.write("Note: When the number of periods is greater than or equal to the number\n"
                      "      of panels, results are based on a generalized inverse of a\n"
                      "      singular matrix.\n", "text")
        if x.get("ncov", 0) >= x["N"]:
            out.write("Note: You estimated at least as many quantities as you have observations.\n", "text")


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
        if np.any(np.abs(rho_i) > 1):
            # the panel estimates are bounded before they are averaged (compat 1905)
            out.write("(note: estimates of rho outside [-1,1] bounded to be in the range [-1,1])\n", "text")
            rho_i = np.clip(rho_i, -1.0, 1.0)
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
    if corr != "independent":
        title = title.replace("Linear regression", "Prais-Winsten regression")     # VERIFY for hetonly/independent
    x = dict(N=N, G=pi.G, tmin=int(pi.counts.min()), tavg=N / pi.G, tmax=int(pi.counts.max()), ncov=ncov, ncorr=ncorr,
             ncf=k, chi2=chi2, df=df, p=p, r2=r2, pan_txt=pan_txt, cor_txt=cor_txt, title=title, pvar=L.pvar, tvar=L.tvar)
    names, bf, Vf = L.expand(b, V)
    est = Estimates("xtpcse", L.smp.depname, names, bf, Vf, N, "z", None)
    est.rows = L.rows(b, V)
    if corr == "ar1":
        est.rows += [CoefRow("sep"), CoefRow("value", "rho", float(rho_g[0]))]
    est.eqnames = [L.smp.depname] * len(names)
    est.cols, est.constant, est.terms, est.depcomp = L.smp.cols, L.constant, L.smp.terms, L.smp.depcomp
    est.level, est.vcetype = L.smp.level, "Het-corrected" if mode == "hetonly" else "Panel-corrected"
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
        tail = [("Estimated covariances", str(x["ncov"])), ("Estimated autocorrelations", str(x["ncorr"])),
                ("Estimated coefficients", str(x["ncf"]))]
        rtail = [("R-squared", _g(x["r2"], "%10.4f")), (f"Wald chi2({x['df']})", _g(x["chi2"], "%10.2f")),
                 ("Prob > chi2", _g(x["p"], "%10.4f"))]
        for (ll_, lv), (rl, rv) in zip(tail, rtail):
            out.write(f"{ll_:<27}= {lv:>9}{'':>10}{rl:<18}= ", "text")
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
    z = crit("z", L.smp.level, None)

    # rho: Stata runs prais on the equation with the panel means removed (no constant), with rhotype()
    if rho_fixed is not None:
        rho = rho_fixed
    else:
        from .prais import fit_prais
        od = pi.order
        yd, Xd = pi.within(y), pi.within(Xs)
        rho = float(fit_prais(Xd[od], yd[od], pi.has_lag[od], pw=True, rhotype=rhot, twostep=twostep,
                              ssesearch=False, maxiter=100, tol=1e-6, k_slopes=k, log=lambda *_a: None)[0])
    # VERIFY: the data are transformed with -rho and the constant is divided by (1 - rho), as found in
    # compat 1905 (coefficients, sigma_e and F agree with Stata to 8 digits with this sign)
    r = -rho
    G = pi.G

    def pw_tr(v):
        lagv = pi.lagged(v) if v.ndim == 1 else np.column_stack([pi.lagged(v[:, j]) for j in range(v.shape[1])])
        mm = keep if v.ndim == 1 else keep[:, None]
        return np.where(mm, v - r * lagv, np.sqrt(1 - r ** 2) * v)
    if fe:
        inv_k = pi.inv[keep]
        _, inv2 = np.unique(inv_k, return_inverse=True)
        pk = pl.PanelIndex(inv2, pi.tidx[keep])
        N = int(keep.sum())
        ys = pw_tr(y)[keep]
        Xk = pw_tr(Xs)[keep]
        yw, Xw = pk.within(ys), pk.within(Xk)
        A = np.linalg.inv(Xw.T @ Xw)
        b = A @ (Xw.T @ yw)
        e = yw - Xw @ b
        rss = float(e @ e)
        df = N - pk.G - k
        s2 = rss / df
        V = s2 * A
        xbar = Xk.mean(axis=0)
        cons_t = float(ys.mean() - xbar @ b)
        cons = cons_t / (1 - rho)
        se_c = math.sqrt(s2 / N + float(xbar @ V @ xbar))          # VERIFY: not divided by (1 - rho), as in Stata
        bb = np.append(b, cons)
        Vf = np.zeros((k + 1, k + 1))
        Vf[:k, :k] = V
        Vf[k, k] = se_c ** 2
        yo, Xo = y[keep], Xs[keep]
        xb_o = Xo @ b
        u_i = pk.group_means(yo) - pk.group_means(xb_o)
        sigma_u = float(np.std(u_i, ddof=1)) if pk.G > 1 else 0.0
        sigma_e = math.sqrt(s2)
        corr = float(np.corrcoef(u_i[pk.inv], xb_o)[0, 1]) if k else 0.0
        F = float(b @ np.linalg.solve(V, b)) / k if k else SYS
        pF = float(st.f.sf(F, k, df)) if F < SYS else SYS
        Xp = np.column_stack([Xk, np.ones(N)])
        bp = np.linalg.lstsq(Xp, ys, rcond=None)[0]
        rss_p = float(((ys - Xp @ bp) ** 2).sum())
        F_f = ((rss_p - rss) / (pk.G - 1)) / s2
        p_f = float(st.f.sf(F_f, pk.G - 1, df))
        r2w = 1 - rss / float(yw @ yw) if k else 0.0
        r2b = float(np.corrcoef(pk.group_means(xb_o), pk.group_means(yo))[0, 1] ** 2) if k else 0.0
        r2o = float(np.corrcoef(xb_o, yo)[0, 1] ** 2) if k else 0.0
        chi2 = pchi = None
        theta = None
        stat, df_r = "t", float(df)
        title = "FE (within) regression with AR(1) disturbances"
        tmin, tmax, tavg, Ng = int(pk.counts.min()), int(pk.counts.max()), N / pk.G, pk.G
        f = {"df": df}
    else:
        # Baltagi-Wu GLS: the transformed constant is c = (sqrt(1 - r^2), 1 - r, ..., 1 - r)
        N = len(y)
        ys = pw_tr(y)
        Xt = pw_tr(Xs)
        c = np.where(keep, 1 - r, np.sqrt(1 - r ** 2))
        R = np.column_stack([Xt, c])
        cc = pi.group_sums(c * c)
        b0 = np.linalg.lstsq(R, ys, rcond=None)[0]
        u = ys - R @ b0
        cm = pi.group_sums(c * u) / cc
        proj = u - c * cm[pi.inv]
        s_e2 = float(proj @ proj) / (N - G)
        s_u2 = max(float(np.mean(cm ** 2)) - s_e2 * float(np.mean(1 / cc)), 0.0)
        th = 1 - np.sqrt(s_e2 / (s_e2 + cc * s_u2))
        thi = th[pi.inv]
        Rm = pi.group_sums(R * c[:, None]) / cc[:, None]
        Rs = R - thi[:, None] * c[:, None] * Rm[pi.inv]
        ym = pi.group_sums(c * ys) / cc
        ysr = ys - thi * c * ym[pi.inv]
        Ai = np.linalg.inv(Rs.T @ Rs)
        bb = Ai @ (Rs.T @ ysr)
        Vf = s_e2 * Ai
        b, V = bb[:k], Vf[:k, :k]
        sigma_u, sigma_e = math.sqrt(s_u2), math.sqrt(s_e2)
        chi2 = float(bb @ np.linalg.solve(Vf, bb))
        pchi = float(chi2_tail(chi2, k + 1))
        xb_o = Xs @ b
        r2w = float(np.corrcoef(pi.within(xb_o), pi.within(y))[0, 1] ** 2) if k else 0.0
        r2b = float(np.corrcoef(pi.group_means(xb_o), pi.group_means(y))[0, 1] ** 2) if k else 0.0
        r2o = float(np.corrcoef(xb_o, y)[0, 1] ** 2) if k else 0.0
        corr = F = pF = F_f = p_f = None
        theta = float(th[0]) if pi.balanced else None
        stat, df_r = "z", None
        title = "RE GLS regression with AR(1) disturbances"
        tmin, tmax, tavg, Ng = int(pi.counts.min()), int(pi.counts.max()), N / G, G
        f = {}
    rho_fov = sigma_u ** 2 / (sigma_u ** 2 + sigma_e ** 2)
    names, bf, Vff = L.expand(bb, Vf)
    est = Estimates("xtregar", L.smp.depname, names, bf, Vff, N, stat, df_r)
    est.rows = L.rows(bb, Vf)
    est.eqnames = [L.smp.depname] * len(names)
    est.cols, est.constant, est.terms, est.depcomp = L.smp.cols, True, L.smp.terms, L.smp.depcomp
    est.level, est.vcetype = L.smp.level, ""
    est.model = "fe" if fe else "re"
    est.extra = dict(N=N, G=Ng, tmin=tmin, tavg=tavg, tmax=tmax, r2w=r2w, r2b=r2b, r2o=r2o,
                     k=k, df_r=f.get("df"), F=F, pF=pF, corr=corr, chi2=chi2, p=pchi, F_f=F_f, p_f=p_f, title=title,
                     pvar=L.pvar, rho=rho, sigma_u=sigma_u, sigma_e=sigma_e, rho_fov=rho_fov, theta=theta)
    est.display = _display_xtregar
    sc = [("N", float(N)), ("N_g", float(Ng)), ("df_m", float(k)), ("rho_ar", rho), ("sigma_u", sigma_u),
          ("sigma_e", sigma_e), ("rho_fov", rho_fov), ("r2_w", r2w), ("r2_b", r2b), ("r2_o", r2o),
          ("g_min", float(tmin)), ("g_avg", tavg), ("g_max", float(tmax)),
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
        out.write(f"{x['title']:<48}{'Number of obs':<18}= ", "text")
        out.write(f"{comma(x['N']):>10}\n", "result")
        out.write(f"{'Group variable: ' + x['pvar']:<48}{'Number of groups':<18}= ", "text")
        out.write(f"{comma(x['G']):>10}\n", "result")
        out.write(f"{'R-sq:':<48}Obs per group:\n", "text")
        rs = [("within ", x["r2w"], "min", _g(x["tmin"])), ("between", x["r2b"], "avg", _g(x["tavg"], "%10.1f")),
              ("overall", x["r2o"], "max", _g(x["tmax"]))]
        for nm, v, gl, gv in rs:
            out.write(f"{'     ' + nm + ' = ' + format(v, '.4f'):<62}{gl} = ", "text")
            out.write(f"{gv:>10}\n", "result")
        if est.model == "fe":
            out.write(f"{'':<48}{'F(' + str(x['k']) + ',' + str(int(x['df_r'])) + ')':<18}= ", "text")
            out.write(f"{_g(x['F'], '%10.2f'):>10}\n", "result")
            out.write(f"{'corr(u_i, Xb)  = ' + format(x['corr'], '.4f'):<48}{'Prob > F':<18}= ", "text")
            out.write(f"{_g(x['pF'], '%10.4f'):>10}\n\n", "result")
        else:
            out.write(f"{'':<48}{'Wald chi2(' + str(x['k'] + 1) + ')':<18}= ", "text")
            out.write(f"{_g(x['chi2'], '%10.2f'):>10}\n", "result")
            out.write(f"{'corr(u_i, Xb)      = 0 (assumed)':<48}{'Prob > chi2':<18}= ", "text")
            out.write(f"{_g(x['p'], '%10.4f'):>10}\n\n", "result")
    if not table:
        return
    coef_table(s, est.depvar, est.rows, stat=est.stat, df=est.df_r, level=level or est.level, footer=False)
    out.write("-------------+----------------------------------------------------------------\n", "text")
    notes = {"rho_ar": "(estimated autocorrelation coefficient)",
             "rho_fov": "(fraction of variance because of u_i)" if est.model == "fe" else "(fraction of variance due to u_i)"}
    rows = [("rho_ar", x["rho"]), ("sigma_u", x["sigma_u"]), ("sigma_e", x["sigma_e"]), ("rho_fov", x["rho_fov"])]
    if x.get("theta") is not None:
        rows.append(("theta", x["theta"]))
    for nm, val in rows:
        out.write(f"{nm:>12}" + " | ", "text")
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
    # e() in the order of Stata's ereturn list (compat 1914)
    sc = [("N", float(N)), ("rank", float(np.linalg.matrix_rank(Vf))), ("chi2", chi2), ("df_m", float(df)),
          ("g_min", float(x["tmin"])), ("g_avg", x["tavg"]), ("g_max", float(x["tmax"])), ("N_g", float(pi2.G)),
          ("chi2_c", r["chi2"]), ("df_chi2c", float(r["df"]))]
    mc = [("cmdline", cmdline("xtrc", args)), ("cmd", "xtrc"), ("predict", "xtrc_p"),
          ("marginsnotok", "GRoup(passthru) STDP"), ("tvar", L.tvar), ("ivar", L.pvar), ("depvar", L.smp.depname),
          ("chi2type", "Wald"), ("title", "Random-coefficients regression"), ("vce", "conventional"),
          ("properties", "b V")]
    from ..commands.matrix import Matrix
    cn = [nm for nm, f in zip(*L.layout()) if f]
    extra_m = [("Sigma", Matrix(r["delta"], cn, cn))]
    if o.get("betas"):
        bb, VV = pl.swamy_blup(r)
        pv = s.data.get(L.pvar).data[L.smp.mask][keep]
        labels = [format_value(float(v), "%9.0g", pad=False).strip() for v in np.unique(pv)]
        x["betas"] = (labels, bb, VV, cn)
        # e(V_ps): one row per panel with its covariance matrix row after row; e(beta_ps): rows
        # r1, r2, ... (compat 1914)
        rn = [f"r{g}" for g in range(1, len(labels) + 1)]
        extra_m = [("V_ps", Matrix(np.vstack([v.ravel() for v in VV]), rn,
                                   [f"{a}:{c}" for a in cn for c in cn])),
                   ("beta_ps", Matrix(bb, rn, cn))] + extra_m
    for _, m_ in extra_m:
        m_.keep_eq = True                     # no equation names (compat 1914)
    post(s, est, sc, mc, L.smp.mask, extra_mats=extra_m)
    _display_xtrc(s, est)


def _display_xtrc(s, est, *, header=True, table=True, level=None, **_kw):
    out, x = s.output, est.extra
    if header:
        out.write("\n", "text")
        out.write(f"{'Random-coefficients regression':<48}{'Number of obs':<18}= ", "text")
        out.write(f"{comma(x['N']):>10}\n", "result")
        out.write(f"{'Group variable: ' + x['pvar']:<48}{'Number of groups':<18}= ", "text")
        out.write(f"{comma(x['G']):>10}\n\n", "result")
        out.write(f"{'':<48}Obs per group:\n", "text")
        for lab, v in (("min", _g(x["tmin"])), ("avg", _g(x["tavg"], "%10.1f")), ("max", _g(x["tmax"]))):
            out.write(f"{'':<62}{lab} = ", "text")
            out.write(f"{v:>10}\n", "result")
        out.write("\n", "text")
        out.write(f"{'':<48}{'Wald chi2(' + str(x['df']) + ')':<18}= ", "text")
        out.write(f"{_g(x['chi2'], '%10.2f'):>10}\n", "result")
        out.write(f"{'':<48}{'Prob > chi2':<18}= ", "text")
        out.write(f"{_g(x['p'], '%10.4f'):>10}\n\n", "result")
    if table:
        coef_table(s, est.depvar, est.rows, stat="z", df=None, level=level or est.level)
        out.write(f"Test of parameter constancy:    chi2({x['c_df']}) = {x['c_chi2']:>8.2f}".ljust(58)
                  + f"Prob > chi2 = {x['c_p']:.4f}\n", "text")
        if "betas" in x:
            # the layout of compat 1914
            labels, bb, VV, cn = x["betas"]
            lv = level or est.level
            c = float(st.norm.ppf(0.5 + lv / 200))
            out.write("\n" + " " * 24 + "Group-specific coefficients\n" + "-" * 78 + "\n", "text")
            out.write(" " * 20 + f"Coef.   Std. Err.      z    P>|z|     [{lv:g}% Conf. Interval]\n" + "-" * 78 + "\n",
                      "text")
            for g, lab in enumerate(labels):
                if g:
                    out.write("-" * 13 + "+" + "-" * 64 + "\n", "text")
                out.write(f"{'Group ' + lab:>12} | \n" + "-" * 13 + "+" + "-" * 64 + "\n", "text")
                for j, nm in enumerate(cn):
                    b_, se_ = float(bb[g, j]), float(np.sqrt(max(VV[g][j, j], 0)))
                    zz = b_ / se_
                    out.write(f"{nm:>12} |", "text")
                    g9_ = lambda v: _g(v, "%9.0g")             # noqa: E731
                    out.write(f"  {g9_(b_):>9}  {g9_(se_):>9}  {zz:>7.2f}  {2 * st.norm.sf(abs(zz)):>6.3f}    "
                              f"{g9_(b_ - c * se_):>9}   {g9_(b_ + c * se_):>9}\n", "result")
            out.write("-" * 78 + "\n", "text")


# ---------------------------------------------------------------------------
# xtivreg
# ---------------------------------------------------------------------------

def _pairs_r2_inv(xb, y, inv):
    """Within, between and overall R-sq (squared correlations) for rows labelled 0..G-1 by `inv`."""
    def c2(a, b):
        if len(a) < 2 or np.std(a) == 0 or np.std(b) == 0:
            return 0.0
        return float(np.corrcoef(a, b)[0, 1] ** 2)
    cnt = np.bincount(inv)

    def gm(v):
        return (np.bincount(inv, weights=v) / cnt)
    return c2(xb - gm(xb)[inv], y - gm(y)[inv]), c2(gm(xb), gm(y)), c2(xb, y)


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
    if mp:
        endog_txt, inst_txt = mp.group(1), mp.group(2)
        rest = t[:mp.start()] + " " + t[mp.end():]
    else:
        raise StataError(198, "equation not specified")   # VERIFY
    L = _Lin(s, "xtivreg", rest.strip(), {"re": 2, "be": 2, "fe": 2, "fd": 2, "ec2sls": 3, "nosa": 3, "small": 2,
                                          "regress": 3, "first": 3}, extra_vars=endog_txt + " " + inst_txt,
             constant=True)
    o, out, pi, ds = L.o, s.output, L.pi, s.data
    if o.get("first") and o.get("regress"):
        raise StataError(198, "options first and regress may not be specified together")
    regress = bool(o.get("regress"))
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
    Z = X.copy() if regress else np.column_stack([c.values for c in inst_cols + ex_cols])
    k = X.shape[1]
    if Z.shape[1] < k:
        raise StataError(481, "equation not identified; must have at least as many instruments not in\n"
                         "the regression as there are instrumented variables")
    y = L.y
    r = pl.panel_iv(y, X, Z, pi, method, bool(o.get("nosa")))
    if o.get("first"):
        _xtiv_first_stage(s, method, endog_cols, ex_cols, inst_cols, pi, r, L.smp.level)
    small = bool(o.get("small"))
    has_cons = True
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
    # VERIFY: the within Wald statistic includes the constant (k degrees of freedom)
    sl = list(range(k + 1)) if method == "fe" else list(range(k))
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
    if method == "fd":
        # levels of the rows that have a lag, as the differenced panel is reported (VERIFY)
        rws = r["rows"]
        inv_r = np.unique(pi.inv[rws], return_inverse=True)[1]
        r2w, r2b, r2o = _pairs_r2_inv(xb_all[rws], y[rws], inv_r)
        tcounts = np.bincount(inv_r)
    else:
        r2w, r2b, r2o = _pairs_r2(xb_all, y, pi)
        tcounts = pi.counts
        if method == "be":
            # 1 - rss/tss of the regression on the group means (compat 1914)
            yb_ = pi.group_means(y)
            rb_ = yb_ - pi.group_means(xb_all) - float(b[k])
            r2b = 1.0 - float(rb_ @ rb_) / float(((yb_ - yb_.mean()) ** 2).sum())
            r2b = r2b if r2b >= 0 else None          # a negative value is shown as . (compat 1905)
        if method == "fe":
            yw_ = pi.within(y)
            r2w = 1.0 - float(r["e"] @ r["e"]) / float(yw_ @ yw_)   # VERIFY: 1 - rss/tss
    if method == "be":
        nobs = len(y)
    x = dict(N=nobs, G=Gn, tmin=int(tcounts.min()), tavg=float(tcounts.mean()), tmax=int(tcounts.max()), k=k, W=W, F=Fv,
             p=pv, tvar=L.tvar,
             use_F=use_F, df_r=r["df"], r2w=r2w, r2b=r2b, r2o=r2o, method=method, pvar=L.pvar,
             instd=" ".join(c.name for c in Xn),
             # with regress the endogenous variables are their own instruments (compat 1914)
             insts=" ".join([c.name for c in ex_cols] + [c.name for c in (Xn if regress else inst_cols)]))
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
        # F test that all u_i = 0: the Wald test of the group dummies in the dummy-variable 2SLS
        # fit, with s2 = rss/df (compat 1914: 11.14, and 11.77 with regress)
        D = (pi.inv[:, None] == np.arange(1, G)[None, :]).astype(float)
        Xd = np.column_stack([X, D, np.ones(len(y))])
        Zd = np.column_stack([Z, D, np.ones(len(y))])
        Pz = Zd @ np.linalg.pinv(Zd.T @ Zd)
        Ad = np.linalg.pinv((Xd.T @ Pz) @ (Zd.T @ Xd))
        bd = Ad @ ((Xd.T @ Pz) @ (Zd.T @ y))
        dd_ = bd[k:k + G - 1]
        Ff = float(dd_ @ np.linalg.solve(s2 * Ad[k:k + G - 1, k:k + G - 1], dd_)) / (G - 1)
        x.update(F_f=Ff, p_f=float(st.f.sf(Ff, G - 1, r["df"])), rss=rss_w)
    elif method == "fd":
        # VERIFY: u_i, sigma_u and corr(u_i, Xb) come from the levels of the rows that have a lag
        cnt = np.bincount(inv_r)
        u_i = np.bincount(inv_r, weights=y[rws] - xb_all[rws]) / cnt - float(b[k])
        su = float(np.std(u_i, ddof=1)) if len(u_i) > 1 else 0.0
        se_ = math.sqrt(s2)
        x.update(sigma_u=su, sigma_e=se_, rho=su ** 2 / (su ** 2 + se_ ** 2),
                 corr=float(np.corrcoef(u_i[inv_r], xb_all[rws])[0, 1]) if k else 0.0, depdisp="D." + L.smp.depname)
    elif method == "be":
        x.update(sd_be=math.sqrt(s2))
    # coefficient order: endogenous first, then exogenous (as in ivregress)
    bmap = {nm: float(b[i]) for i, nm in enumerate(names_all)}
    semap = {nm: float(math.sqrt(max(V[i, i], 0))) for i, nm in enumerate(names_all)}
    rows = rows_from_columns(Xn + ex_cols, bmap, semap, ds, constant=has_cons)
    if method == "fd":
        rows = []
        for c in Xn + ex_cols:
            rows += [CoefRow("header", c.name), CoefRow("coef", "D1.", bmap[c.name], semap[c.name], name=c.name),
                     CoefRow("blank")]
        rows.append(CoefRow("coef", "_cons", bmap["_cons"], semap["_cons"], name="_cons"))
    est = Estimates("xtivreg", L.smp.depname, names_all, b, V, nobs, stat, df_r)
    est.rows = rows
    est.eqnames = [L.smp.depname] * len(names_all)
    est.cols, est.constant, est.terms, est.depcomp = Xn + ex_cols, has_cons, endog_terms + L.smp.terms, L.smp.depcomp
    est.level, est.vcetype = L.smp.level, ""
    est.extra = x
    est.display = _display_xtivreg
    ttl = {"g2sls": "G2SLS random-effects IV regression", "ec2sls": "EC2SLS random-effects IV regression",
           "fe": "Fixed-effects (within) IV regression", "be": "Between-effects IV regression:",
           "fd": "First-differenced IV regression"}[method]
    x["title"] = ttl
    x["regress"] = regress
    sc = [("N", float(nobs)), ("N_g", float(Gn)), ("df_m", float(k)), ("g_min", float(x["tmin"])),
          ("g_avg", x["tavg"]), ("g_max", float(x["tmax"])), ("r2_w", r2w), ("r2_b", r2b), ("r2_o", r2o),
          ("rank", float(k + (1 if has_cons else 0)))]
    sc += [("F", Fv)] if use_F else [("chi2", W)]
    if "sigma_u" in x:
        sc += [("sigma_u", x["sigma_u"]), ("sigma_e", x["sigma_e"]), ("rho", x["rho"])]
    mc = [("cmdline", cmdline("xtivreg", args)), ("cmd", "xtivreg"), ("depvar", L.smp.depname), ("ivar", L.pvar),
          ("model", method), ("instd", x["instd"]), ("insts", x["insts"]), ("title", ttl), ("properties", "b V"),
          ("predict", "xtivreg_p")]
    if method == "fe":
        # Stata's e() of the within estimator, in the order of ereturn list (compat 1914, with
        # regress; VERIFY without it)
        sc = [("rank", float(k + 1)), ("df_rz", float(r["df"])), ("N", float(nobs)), ("rss", x["rss"]),
              ("N_g", float(Gn)), ("N_clust", float(Gn)), ("chi2", W), ("chi2_p", pv), ("sigma_u", x["sigma_u"]),
              ("corr", x["corr"]), ("r2_o", r2o), ("r2_b", r2b), ("sigma_e", x["sigma_e"]),
              ("sigma", math.sqrt(x["sigma_u"] ** 2 + x["sigma_e"] ** 2)), ("rho", x["rho"]), ("r2_w", r2w),
              ("df_b", float(k)), ("df_m", float(k + G)), ("df_a", float(G - 1)), ("g_min", float(x["tmin"])),
              ("F_fp", x["p_f"]), ("F_f", x["F_f"]), ("g_avg", x["tavg"]), ("g_max", float(x["tmax"]))]
        mc = [("cmdline", cmdline("xtivreg", args)), ("cmd", "xtivreg"), ("insts", x["insts"]),
              ("instd", x["instd"]), ("marginsok", "XB default"), ("marginsnotok", "ue xbu u e"),
              ("predict", "xtivp_1"), ("depvar", L.smp.depname), ("model", "fe"), ("tvar", L.tvar),
              ("ivar", L.pvar), ("vce", "conventional"), ("properties", "b V")]
        post(s, est, sc, mc, mask)
    else:
        post(s, est, list(reversed(sc)), list(reversed(mc)), mask)
    _display_xtivreg(s, est)


def _xtiv_first_stage(s, method, endog_cols, ex_cols, inst_cols, pi, r, level) -> None:
    """xtivreg, first: each endogenous variable on the exogenous regressors and the instruments
    (compat 1914): xtreg, fe (in its older layout) or be, regress of the differences for fd, and
    for the random-effects models the regression of the GLS-transformed variables with the theta
    of the IV fit."""
    out = s.output
    regs = " ".join(c.name for c in ex_cols + inst_cols)
    for e in endog_cols:
        if method in ("fe", "be"):
            out.write("\nFirst-stage " + ("within" if method == "fe" else "between") + " regression\n", "text")
            s._xtreg_old = method == "fe"
            try:
                s.interp._execute_expanded(f"xtreg {e.name} {regs}, {method}")
            finally:
                s._xtreg_old = False
            continue
        if method == "fd":
            out.write("\nFirst-stage first-differenced regression\n", "text")
            s.interp._execute_expanded(f"regress D.{e.name} " + " ".join(f"D.{c.name}" for c in ex_cols + inst_cols))
            continue
        # G2SLS (EC2SLS: VERIFY) first stage: OLS of x - theta mean(x) on the transformed
        # regressors and 1 - theta, with s2 = rss/(n - k)
        th = np.asarray(r["theta"], dtype=float)
        th = th[pi.inv] if th.size == pi.G else np.broadcast_to(th, (len(pi.inv),))
        tr = lambda v: v - th * pi.group_means(v)[pi.inv]          # noqa: E731
        Xf = np.column_stack([tr(c.values) for c in ex_cols + inst_cols] + [1 - th])
        yf = tr(e.values)
        A = np.linalg.inv(Xf.T @ Xf)
        b = A @ (Xf.T @ yf)
        res = yf - Xf @ b
        V = float(res @ res) / (len(yf) - Xf.shape[1]) * A
        kk = Xf.shape[1] - 1
        W = float(b[:kk] @ np.linalg.solve(V[:kk, :kk], b[:kk]))
        out.write("\nFirst-stage " + ("G2SLS" if method == "g2sls" else "EC2SLS") + " regression\n", "text")
        for lab, val in (("Number of obs", f"{len(yf):,}"), (f"Wald chi({kk})", f"{W:.0f}"),
                         ("Prob > chi2", _g(float(chi2_tail(W, kk)), "%10.4f"))):
            out.write(f"{'':<49}{lab:<17}= ", "text")
            out.write(f"{val:>10}\n", "result")
        out.write("\n", "text")
        se = np.sqrt(np.maximum(np.diag(V), 0))
        names = [c.name for c in ex_cols + inst_cols] + ["_cons"]
        coef_table(s, e.name, [CoefRow("coef", nm, float(b[i]), float(se[i]), name=nm) for i, nm in enumerate(names)],
                   stat="z", df=None, level=level)


def _display_xtivreg(s, est, *, header=True, table=True, level=None, **_kw):
    out, x = s.output, est.extra
    if header:
        out.write("\n", "text")
        if x["method"] == "fd":
            out.write(f"{x['title']}\n", "text")
            l1, l2 = f"{'Group variable:':<19}{x['pvar']}", f"{'Time variable:':<19}{x['tvar']}"
        else:
            l1, l2 = x["title"], f"Group variable: {x['pvar']}"
        out.write(f"{l1:<48}{'Number of obs':<18}= ", "text")
        out.write(f"{comma(x['N']):>10}\n", "result")
        out.write(f"{l2:<48}{'Number of groups':<18}= ", "text")
        out.write(f"{comma(x['G']):>10}\n", "result")
        out.write(f"{'R-sq:':<48}Obs per group:\n", "text")
        rs = [("within ", x["r2w"], "min", _g(x["tmin"])), ("between", x["r2b"], "avg", _g(x["tavg"], "%10.1f")),
              ("overall", x["r2o"], "max", _g(x["tmax"]))]
        for nm, v, gl, gv in rs:
            vs = "     ." if v is None else f"{v:.4f}"
            out.write(f"     {nm} = {vs}".ljust(62), "text")
            out.write(f"{gl} = ", "text")
            out.write(f"{gv:>10}\n", "result")
        stat_lab, stat_val, pl_, pv_ = (f"F({x['k']},{int(x['df_r'])})", x["F"], "Prob > F", x["p"]) if x["use_F"] else (
            f"Wald chi2({x['k']})", x["W"], "Prob > chi2", x["p"])
        out.write(f"{'':<48}{stat_lab:<18}= ", "text")
        out.write(f"{_g(stat_val, '%10.2f'):>10}\n", "result")
        if "corr" in x:
            cl = f"corr(u_i, Xb)  = {x['corr']:.4f}"
        elif x["method"] == "be":
            cl = "sd(u_i + avg(e_i.))= " + format_value(x["sd_be"], "%9.0g")
        elif x["method"] in ("g2sls", "ec2sls"):
            cl = "corr(u_i, X)       = 0 (assumed)"
        else:
            cl = ""
        out.write(f"{cl:<48}{pl_:<18}= ", "text")
        out.write(f"{_g(pv_, '%10.4f'):>10}\n", "result")
    if not table:
        return
    coef_table(s, x.get("depdisp", est.depvar), est.rows, stat=est.stat, df=est.df_r, level=level or est.level,
               footer=False)
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
        out.write(f"F  test that all u_i=0:     F({int(x['G'] - 1)},{int(x['df_r'])}) = {x['F_f']:8.2f}          "
                  f"Prob > F    = {x['p_f']:6.4f}\n", "text")
        out.write("-" * 78 + "\n", "text")
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
            out_rows.append(CoefRow("header", gname, left=True))
            out_rows.extend(rows[pos:pos + n_g])
            pos += n_g
    if out_rows and out_rows[-1].kind != "blank" and rows[pos:] and rows[pos].label == "_cons" and grp["TIexogenous"]:
        out_rows.append(CoefRow("blank"))      # VERIFY: blank line before _cons after a TI block
    out_rows.extend(rows[pos:])
    est = Estimates("xthtaylor", L.smp.depname, order + ["_cons"], b, V, N, "z", None)
    est.rows = out_rows
    est.eqnames = [L.smp.depname] * (k + 1)
    est.cols, est.constant, est.terms, est.depcomp = ordered_cols, True, L.smp.terms, L.smp.depcomp
    est.level, est.vcetype = L.smp.level, ""
    est.extra = dict(N=N, G=G, tmin=int(pi.counts.min()), tavg=N / G, tmax=int(pi.counts.max()), W=W, k=k, p=p,
                     sigma_u=su, sigma_e=se_, rho=rho, pvar=L.pvar, tvar=L.tvar,
                     title="Amemiya-MaCurdy estimation" if o.get("amacurdy") else "Hausman-Taylor estimation")
    est.display = _display_xthtaylor
    sc = [("N", float(N)), ("N_g", float(G)), ("df_m", float(k)), ("chi2", W), ("sigma_u", su), ("sigma_e", se_),
          ("rho", rho), ("g_min", float(est.extra["tmin"])), ("g_avg", est.extra["tavg"]),
          ("g_max", float(est.extra["tmax"])), ("rank", float(k + 1))]
    mc = [("cmdline", cmdline("xthtaylor", args)), ("cmd", "xthtaylor"), ("depvar", L.smp.depname), ("ivar", L.pvar),
          ("model", "amacurdy" if o.get("amacurdy") else "ht"), ("title", est.extra["title"]),
          ("properties", "b V"), ("predict", "xthtaylor_p")]
    post(s, est, list(reversed(sc)), list(reversed(mc)), L.smp.mask)
    _display_xthtaylor(s, est)


def _display_xthtaylor(s, est, *, header=True, table=True, level=None, **_kw):
    out, x = s.output, est.extra
    if header:
        out.write("\n", "text")
        out.write(f"{x['title']:<48}{'Number of obs':<18}= ", "text")
        out.write(f"{comma(x['N']):>10}\n", "result")
        out.write(f"{'Group variable: ' + x['pvar']:<48}{'Number of groups':<18}= ", "text")
        out.write(f"{comma(x['G']):>10}\n\n", "result")
        out.write(f"{('Time variable: ' + x['tvar']) if x['title'].startswith('Amemiya') else '':<48}Obs per group:\n",
                  "text")
        for lab, v in (("min", x["tmin"]), ("avg", x["tavg"]), ("max", x["tmax"])):
            out.write(f"{'':<62}{lab} = ", "text")
            out.write(f"{_g(v):>10}\n", "result")
        out.write(f"{'Random effects u_i ~ i.i.d.':<48}{'Wald chi2(' + str(x['k']) + ')':<18}= ", "text")
        out.write(f"{_g(x['W'], '%10.2f'):>10}\n", "result")
        out.write(f"{'':<48}{'Prob > chi2':<18}= ", "text")
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
    out.write("Note: TV refers to time varying; TI refers to time invariant.\n", "text")
