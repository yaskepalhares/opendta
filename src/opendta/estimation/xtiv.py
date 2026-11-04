"""xtivreg and xthtaylor ([XT] xtivreg, xthtaylor): instrumental-variables estimators for panel data
(between, fixed, random effects G2SLS and EC2SLS, first differences) and the Hausman-Taylor
estimator."""

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
from .regress import cmdline, comma
from .results import CoefRow, Estimates, coef_table, post, rows_from_columns
from .xtlin import _Lin, _g

if TYPE_CHECKING:
    from ..session import Session

# imported on first use (start-up time; see opendta/lazy.py)
st = lazy("scipy.stats")

SYS = M.SYSMISS


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
