"""Modelos de máxima verossimilhança: logit, logistic, probit, poisson
([R] logit, [R] logistic, [R] probit, [R] poisson)."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy import special, stats as st

from ..core import missing as M
from ..core.errors import StataError
from ..core.formats import format_value
from . import fvars
from .ml import maximize, print_log
from .regress import Sample, cmdline, comma, wexp_macros
from .results import CoefRow, Estimates, coef_table, post, rows_from_columns

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def g(x: float, fmt: str) -> str:
    if x is None or not np.isfinite(x) or abs(x) >= SYS:
        return "."
    return format_value(float(x), fmt, pad=False).strip()


# ---------------------------------------------------------------------------
# verossimilhanças: (ll, gradiente, Hessiana) e escores por observação
# ---------------------------------------------------------------------------

def _logit(y, X, w):
    def fun(b, todo):
        xb = X @ b
        # ln F(xb) e ln(1 - F(xb)) estáveis
        lp = -np.logaddexp(0, -xb)
        lq = -np.logaddexp(0, xb)
        ll = float(np.sum(w * np.where(y > 0, lp, lq)))
        if todo == 0:
            return ll, None, None
        p = special.expit(xb)
        r = (y > 0) - p
        gr = X.T @ (w * r)
        H = -(X * (w * p * (1 - p))[:, None]).T @ X
        return ll, gr, H

    def scores(b):
        p = special.expit(X @ b)
        return (y > 0) - p
    return fun, scores


def _probit(y, X, w):
    q = np.where(y > 0, 1.0, -1.0)

    def fun(b, todo):
        xb = X @ b
        ll = float(np.sum(w * st.norm.logcdf(q * xb)))
        if todo == 0:
            return ll, None, None
        lam = q * np.exp(st.norm.logpdf(q * xb) - st.norm.logcdf(q * xb))
        gr = X.T @ (w * lam)
        H = -(X * (w * lam * (lam + xb))[:, None]).T @ X
        return ll, gr, H

    def scores(b):
        xb = X @ b
        return q * np.exp(st.norm.logpdf(q * xb) - st.norm.logcdf(q * xb))
    return fun, scores


def _poisson(y, X, w, offset):
    lgy = special.gammaln(y + 1)

    def fun(b, todo):
        xb = X @ b + offset
        mu = np.exp(np.clip(xb, -700, 700))
        ll = float(np.sum(w * (-mu + y * xb - lgy)))
        if todo == 0:
            return ll, None, None
        gr = X.T @ (w * (y - mu))
        H = -(X * (w * mu)[:, None]).T @ X
        return ll, gr, H

    def scores(b):
        return y - np.exp(X @ b + offset)
    return fun, scores


# ---------------------------------------------------------------------------
# comando comum
# ---------------------------------------------------------------------------

_TITLES = {"logit": "Logistic regression", "logistic": "Logistic regression",
           "probit": "Probit regression", "poisson": "Poisson regression"}


def _fit(s: "Session", cmd: str, args: str) -> None:
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, "logit" if cmd == "logistic" else cmd, args)
        return
    spec = {"or": 2, "irr": 3, "nolog": 4, "log": 3, "iterate": 4, "tolerance": 3, "ltolerance": 4,
            "nrtolerance": 4, "asis": 4, "offset": 3, "exposure": 3, "noheader": 6, "coef": 4,
            "difficult": 4, "technique": 4}
    smp = Sample(s, args, opts_spec=spec, weights=("fweight", "iweight", "pweight"))
    o = smp.o
    out = s.output
    ds = s.data
    w = smp.w[smp.mask].astype(np.float64) if smp.w is not None else np.ones(smp.n)
    if smp.wtype == "pweight":
        w = w.copy()
    N = float(w.sum()) if smp.wtype in ("fweight",) else float(smp.n)
    if smp.wtype == "iweight":
        N = float(w.sum())
    y = smp.y.astype(np.float64)
    if cmd in ("logit", "logistic", "probit"):
        yb = (y != 0).astype(float)
        if np.all(yb == yb[0]):
            raise StataError(2000, "outcome does not vary; remember:\n"
                             "                                  0 = negative outcome,\n"
                             "        all other nonmissing values = positive outcome")
        y = yb
    elif cmd == "poisson":
        if np.any(y < 0):
            raise StataError(459, f"{smp.depname} must be greater than or equal to zero")   # VERIFICAR
    notes = fvars.drop_collinear(smp.cols, smp.constant, w)
    for nm in notes:
        out.write(f"note: {nm} omitted because of collinearity\n", "text")
    est_cols = [c for c in smp.cols if c.values is not None]
    X = np.column_stack([c.values for c in est_cols] + ([np.ones(smp.n)] if smp.constant else [])) \
        if (est_cols or smp.constant) else np.zeros((smp.n, 0))
    k = X.shape[1]
    offset = np.zeros(smp.n)
    off_name = ""
    if cmd == "poisson":
        if o.get("offset"):
            off_name = str(o["offset"]).strip()
            offset = ds.get(off_name).data[smp.mask].astype(float)
        elif o.get("exposure"):
            off_name = str(o["exposure"]).strip()
            ev = ds.get(off_name).data[smp.mask].astype(float)
            if np.any(ev <= 0):
                raise StataError(459, f"{off_name} <= 0 in sample")   # VERIFICAR
            offset = np.log(ev)
            off_name = f"ln({off_name})"
    # modelo só com a constante (ll_0) e valores iniciais
    ybar = float(np.sum(w * y) / np.sum(w))
    if cmd in ("logit", "logistic"):
        fun, scores = _logit(y, X, w)
        c0 = np.log(ybar / (1 - ybar))
        ll0 = float(np.sum(w * (y * np.log(ybar) + (1 - y) * np.log(1 - ybar))))
    elif cmd == "probit":
        fun, scores = _probit(y, X, w)
        c0 = float(st.norm.ppf(ybar))
        ll0 = float(np.sum(w * (y * np.log(ybar) + (1 - y) * np.log(1 - ybar))))
    else:
        fun, scores = _poisson(y, X, w, offset)
        if np.any(offset != 0):
            c0 = float(np.log(np.sum(w * y) / np.sum(w * np.exp(offset))))
            mu0 = np.exp(c0 + offset)
            ll0 = float(np.sum(w * (-mu0 + y * np.log(mu0) - special.gammaln(y + 1))))
        else:
            c0 = np.log(ybar) if ybar > 0 else -20.0
            ll0 = float(np.sum(w * (-ybar + y * np.log(ybar if ybar > 0 else 1) - special.gammaln(y + 1))))
    b0 = np.zeros(k)
    if smp.constant and k:
        b0[-1] = c0
    if not smp.constant:
        ll0 = SYS
    maxiter = int(o["iterate"]) if o.get("iterate") else 300
    # VERIFICAR: valores iniciais do poisson no Stata (o log de iterações pode diferir)
    res = maximize(fun, b0, maxiter=maxiter,
                   tol=float(o["tolerance"]) if o.get("tolerance") else 1e-6,
                   ltol=float(o["ltolerance"]) if o.get("ltolerance") else 1e-7,
                   nrtol=float(o["nrtolerance"]) if o.get("nrtolerance") else 1e-5)
    robust = smp.vce in ("robust", "cluster")
    lab = "log pseudolikelihood" if robust else "log likelihood"
    if cmd != "logistic" and not o.get("nolog"):
        out.write("\n", "text")
        print_log(out, res, lab)
    try:
        Vb = np.linalg.inv(-res.H) if k else np.zeros((0, 0))
    except np.linalg.LinAlgError:
        Vb = np.linalg.pinv(-res.H)
    Vb = (Vb + Vb.T) / 2
    n_clust = None
    if robust:
        sc = X * (w * scores(res.b))[:, None]
        if smp.vce == "cluster":
            ids = smp.cluster
            uniq = np.unique(ids)
            n_clust = len(uniq)
            meat = np.zeros((k, k))
            for u in uniq:
                sg = sc[ids == u].sum(axis=0)
                meat += np.outer(sg, sg)
            q = n_clust / (n_clust - 1)
        else:
            if smp.wtype == "fweight":
                meat = (X * (w * scores(res.b) ** 2)[:, None]).T @ X
            else:
                meat = sc.T @ sc
            q = N / (N - 1)
        Vb = q * Vb @ meat @ Vb
        Vb = (Vb + Vb.T) / 2
    slopes = list(range(len(est_cols)))
    df_m = len(slopes)
    if robust:
        chi2type = "Wald"
        try:
            chi2 = float(res.b[slopes] @ np.linalg.solve(Vb[np.ix_(slopes, slopes)], res.b[slopes])) \
                if slopes else 0.0
        except np.linalg.LinAlgError:
            chi2 = SYS
    else:
        chi2type = "LR"
        chi2 = 2 * (res.ll - ll0) if ll0 < SYS else SYS
    pval = float(st.chi2.sf(chi2, df_m)) if chi2 < SYS and df_m > 0 else SYS
    r2_p = 1 - res.ll / ll0 if ll0 < SYS and ll0 != 0 else SYS
    names = [c.name if not c.omitted else fvars.omitted_name(c.name) for c in smp.cols]
    if smp.constant:
        names.append("_cons")
    kk = len(names)
    b_full = np.zeros(kk)
    V_full = np.zeros((kk, kk))
    pos = [i for i, c in enumerate(smp.cols) if c.values is not None] + ([kk - 1] if smp.constant else [])
    b_full[pos] = res.b
    V_full[np.ix_(pos, pos)] = Vb
    bmap, semap = {}, {}
    for i, c in enumerate(smp.cols):
        bmap[c.name] = b_full[i]
        semap[c.name] = np.sqrt(V_full[i, i]) if V_full[i, i] > 0 else 0.0
    if smp.constant:
        bmap["_cons"], semap["_cons"] = b_full[-1], np.sqrt(max(V_full[-1, -1], 0))
    base_cmd = "logit" if cmd == "logistic" else cmd
    est = Estimates(base_cmd, smp.depname, names, b_full, V_full, smp.n, "z", None)
    est.rows = rows_from_columns(smp.cols, bmap, semap, ds, constant=smp.constant)
    est.cols, est.constant, est.terms, est.depcomp = smp.cols, smp.constant, smp.terms, smp.depcomp
    est.level = smp.level
    est.vcetype = "Robust" if robust else ""
    est.clustvar = smp.clustvar
    est.n_clust = n_clust
    est.title = _TITLES[cmd]
    est.offset_name = off_name
    est.offset_kind = "exposure" if o.get("exposure") else "offset" if o.get("offset") else ""
    est.extra = {"N": N, "ll": res.ll, "ll_0": ll0, "chi2": chi2, "p": pval, "df_m": df_m,
                 "r2_p": r2_p, "chi2type": chi2type}
    est.robust = robust
    est.display = display_ml
    est.predict = predict_ml
    est.eform_default = {"logistic": "Odds Ratio"}.get(cmd, "")
    est.ml_cmd = cmd
    scalars = [("rank", float(k)), ("N", N), ("ic", float(res.iterations)), ("k", float(kk)),
               ("k_eq", 1.0), ("k_dv", 1.0), ("converged", 1.0 if res.converged else 0.0), ("rc", 0.0),
               ("ll", res.ll)]
    if n_clust is not None:
        scalars.append(("N_clust", float(n_clust)))
    scalars += [("k_eq_model", 1.0), ("ll_0", ll0), ("df_m", float(df_m)), ("chi2", chi2), ("p", pval)]
    if cmd != "poisson":
        scalars.append(("r2_p", r2_p))
    else:
        scalars.append(("r2_p", r2_p))
    macros = [("cmdline", cmdline(cmd, args)), ("cmd", base_cmd if cmd != "logistic" else "logistic"),
              ("title", _TITLES[cmd]), ("chi2type", chi2type), ("vce", smp.vce if robust else "oim")]
    if robust:
        macros.append(("vcetype", "Robust"))
    if smp.clustvar:
        macros.append(("clustvar", smp.clustvar))
    macros += wexp_macros(smp)
    macros += [("depvar", smp.depname), ("opt", "moptimize"), ("which", "max"),
               ("ml_method", "d2"), ("technique", "nr"), ("properties", "b V"),
               ("predict", f"{base_cmd}_p"), ("estat_cmd", f"{base_cmd}_estat")]
    if off_name:
        macros.append((est.offset_kind, off_name))
    post(s, est, scalars, macros, smp.mask)
    eform = ""
    if cmd == "logistic" and not o.get("coef"):
        eform = "Odds Ratio"
    elif o.get("or") and cmd == "logit":
        eform = "Odds Ratio"
    elif o.get("irr") and cmd == "poisson":
        eform = "IRR"
    display_ml(s, est, eform=eform, header=not o.get("noheader"))
    if not res.converged:
        out.write("convergence not achieved\n", "error")   # VERIFICAR


def display_ml(s: "Session", est: Estimates, *, eform: str | None = None, header: bool = True,
               table: bool = True, level: float | None = None, **_kw) -> None:
    out = s.output
    x = est.extra
    lev = level if level is not None else est.level
    if eform is None:
        eform = est.eform_default
    if header:
        out.write("\n", "text")
        chi_lab = f"{x['chi2type']} chi2({int(x['df_m'])})"
        right = [("Number of obs", comma(x["N"])), (chi_lab, g(x["chi2"], "%10.2f")),
                 ("Prob > chi2", g(x["p"], "%10.4f")), ("Pseudo R2", g(x["r2_p"], "%10.4f"))]
        lab = "Log pseudolikelihood" if est.robust else "Log likelihood"
        left = [est.title, "", "", f"{lab} = {g(x['ll'], '%10.0g')}"]
        for lft, (rl, rv) in zip(left, right):
            out.write(f"{lft:<48}{rl:<18}= ", "text")
            out.write(f"{rv:>10}\n", "result")
        out.write("\n", "text")
    if not table:
        return
    if est.vcetype and est.clustvar:
        out.write(f"{'(Std. Err. adjusted for ' + format(int(est.n_clust), ',') + ' clusters in ' + est.clustvar + ')':>78}\n",
                  "text")
    rows = est.rows
    title = "Coef."
    if eform:
        title = eform
        rows = [_eform_row(r) for r in rows if not (r.kind == "coef" and r.label == "_cons" and False)]
    coef_table(s, est.depvar, rows, stat="z", df=None, level=lev, vcetype=est.vcetype,
               coef_title=title, eform=bool(eform))
    if est.offset_name:
        pass


def _eform_row(r: CoefRow) -> CoefRow:
    return r


# ---------------------------------------------------------------------------
# predict
# ---------------------------------------------------------------------------

def predict_ml(s: "Session", est: Estimates, stat: str, mask: np.ndarray):
    from .postest import design_all
    X, ok = design_all(s, est)
    xb = X @ est.b
    if est.offset_name and stat not in ("xb_nooffset",):
        ds = s.data
        name = est.offset_name
        if name.startswith("ln("):
            v = ds.get(name[3:-1]).data
            off = np.where(v < SYS, np.log(np.where(v > 0, v, 1)), 0)
        else:
            v = ds.get(name).data
            off = np.where(v < SYS, v, 0)
        xb = xb + off
    note = ""
    cmd = est.ml_cmd
    if stat == "xb":
        return np.where(ok, xb, SYS), ""
    if stat == "stdp":
        sp = np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", X, est.V, X), 0))
        return np.where(ok, sp, SYS), ""
    if cmd in ("logit", "logistic", "probit"):
        if stat in ("", "pr", "p"):
            if not stat:
                note = f"(option pr assumed; Pr({est.depvar}))"
            p = special.expit(xb) if cmd != "probit" else st.norm.cdf(xb)
            return np.where(ok, p, SYS), note
    if cmd == "poisson":
        if stat in ("", "n"):
            if not stat:
                note = "(option n assumed; predicted number of events)"
            return np.where(ok, np.exp(xb), SYS), note
        if stat == "ir":
            return np.where(ok, np.exp(X @ est.b), SYS), ""
    raise StataError(198, f"option {stat} not allowed")


def cmd_logit(s, args):
    _fit(s, "logit", args)


def cmd_logistic(s, args):
    _fit(s, "logistic", args)


def cmd_probit(s, args):
    _fit(s, "probit", args)


def cmd_poisson(s, args):
    _fit(s, "poisson", args)
