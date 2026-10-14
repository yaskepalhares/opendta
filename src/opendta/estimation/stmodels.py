"""Survival regression: stcox and streg ([ST] stcox, [ST] streg,
[ST] stcox postestimation, [ST] stcox PH-assumption tests).

stcox maximizes the Cox partial likelihood with the Breslow (default) or
Efron approximation for ties; delayed entry and strata enter through the
risk sets. The robust variance uses the efficient score residuals (Lin and
Wei 1989). streg fits the exponential, Weibull, Gompertz, lognormal and
loglogistic models with the full likelihood of [ST] streg, Methods and
formulas, including delayed entry. Layout details carry the VERIFY mark
(compat/do/0902_stcox_streg.do).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np
from scipy import special, stats as st

from ..core import missing as M
from ..core.errors import StataError
from ..core.formats import format_value
from ..core.stats import chi2_tail
from . import fvars
from .ml import maximize, print_log
from .results import CoefRow, Estimates, coef_table, post, rows_from_columns
from .st import header, n_subjects, require_st, st_data

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def _g(x: float, fmt: str = "%10.0g") -> str:
    if x is None or not np.isfinite(x) or abs(x) >= SYS:
        return "."
    return format_value(float(x), fmt, pad=False).strip()


# ---------------------------------------------------------------------------
# sample shared by stcox and streg
# ---------------------------------------------------------------------------

class STSample:
    def __init__(self, s: "Session", args: str, spec: dict):
        from .regress import Sample
        info = require_st(s)
        ds = s.data
        self.info = info
        st_ok = ds.get("_st").data == 1
        # the regressors come after a dummy depvar (_t) so that Sample builds the design matrix
        wt = ""
        if info.wtype:
            wt = f" [{info.wtype}={info.wvar}]"
        body, comma, opts = _split_opts(args)
        if "[" in body:
            raise StataError(101, "weights not allowed; specify weights with stset")   # VERIFY
        text = f"_t {body}{wt}" + (f", {opts}" if comma else "")
        extra = " ".join(v for v in (info.idvar,) if v and not ds.get(v).is_string)
        smp = Sample(s, text, opts_spec=spec, weights=("fweight", "pweight", "iweight"),
                     extra_vars=extra)
        m = smp.mask & st_ok
        if not np.array_equal(m, smp.mask):
            smp = Sample(s, _with_if(text, "_st == 1"), opts_spec=spec,
                         weights=("fweight", "pweight", "iweight"), extra_vars=extra)
        self.smp = smp
        self.o = smp.o
        self.mask = smp.mask
        sd = st_data(s, smp.mask)
        self.sd = sd
        self.t, self.t0, self.d = sd.t, sd.t0, sd.d
        self.w = smp.w[smp.mask].astype(np.float64) if smp.w is not None else np.ones(smp.n)
        self.wtype = smp.wtype or ""
        self.n = smp.n
        notes = fvars.drop_collinear(smp.cols, False, self.w)
        self.notes = notes
        self.est_cols = [c for c in smp.cols if c.values is not None]
        self.X = np.column_stack([c.values for c in self.est_cols]) if self.est_cols else np.zeros((smp.n, 0))


def _split_opts(args: str):
    from ..lang.syntax import find_top
    c = find_top(args, ",")
    if c == -1:
        return args.strip(), False, ""
    return args[:c].strip(), True, args[c + 1:].strip()


def _with_if(text: str, cond: str) -> str:
    body, comma, opts = _split_opts(text)
    from ..lang.syntax import _find_keyword
    k = _find_keyword(body, "if")
    if k == -1:
        body = f"{body} if {cond}"
    else:
        body = f"{body[:k]} if ({body[k + 2:].strip()}) & {cond}"
    return body + (f", {opts}" if comma else "")


def _st_header_block(s: "Session", title: str, est: Estimates) -> None:
    x = est.extra
    out = s.output
    out.write(f"\n{title}\n\n", "text")
    robust = est.robust
    llab = "Log pseudolikelihood" if robust else "Log likelihood"
    lw = len(llab) if robust else 15
    chi = f"{x['chi2type']} chi2({int(x['df_m'])})"
    left = [("No. of subjects", _g(x["N_sub"], "%12.0gc")), ("No. of failures", _g(x["N_fail"], "%12.0gc")),
            ("Time at risk", _g(x["risk"], "%12.0g")), None, (llab, _g(x["ll"], "%10.0g"))]
    right = [("Number of obs", _g(x["N"], "%12.0gc")), None, None, (chi, f"{x['chi2']:.2f}" if x["chi2"] < SYS else "."),
             ("Prob > chi2", f"{x['p']:.4f}" if x["p"] < SYS else ".")]
    for lft, rgt in zip(left, right):
        line = ""
        if lft is not None:
            line = f"{lft[0]:<{lw}} ={lft[1]:>13}"
        if rgt is not None:
            out.write(f"{line:<48}{rgt[0]:<17}=", "text")
            out.write(f"{rgt[1]:>12}\n", "result")
        else:
            out.write(line + "\n", "text")
    out.write("\n", "text")


def _st_scalars(stx: STSample) -> dict:
    sd = stx.sd
    fw = stx.w if stx.wtype == "fweight" else np.ones(stx.n)
    return {"N_sub": n_subjects(sd), "N_fail": float(np.sum(stx.d * fw)),
            "risk": float(np.sum((stx.t - stx.t0) * fw)),
            "N": float(fw.sum()) if stx.wtype == "fweight" else float(stx.n)}


def _cluster_ids(stx: STSample):
    """Clusters of the robust variance: vce(cluster) or, with robust and stset id(), the subjects."""
    smp = stx.smp
    if smp.vce == "cluster":
        return smp.cluster, smp.clustvar
    if smp.vce == "robust" and stx.sd.ids is not None:
        return stx.sd.ids, stx.info.idvar
    return None, ""


# ---------------------------------------------------------------------------
# Cox partial likelihood
# ---------------------------------------------------------------------------

class _Stratum:
    def __init__(self, t, t0, d, X, w):
        self.t, self.t0, self.d, self.X, self.w = t, t0, d, X, w
        self.ot = np.argsort(t, kind="stable")
        self.ts = t[self.ot]
        self.ot0 = np.argsort(t0, kind="stable")
        self.t0s = t0[self.ot0]
        fail = d > 0
        self.ft = np.unique(t[fail])                      # distinct failure times
        self.fidx = np.searchsorted(self.ft, t[fail])
        self.fail = fail
        m = len(self.ft)
        self.D = np.bincount(self.fidx, weights=w[fail], minlength=m)         # weighted failures
        self.sx = np.zeros((m, X.shape[1]))
        np.add.at(self.sx, self.fidx, X[fail] * w[fail][:, None])
        self.ia = np.searchsorted(self.ts, self.ft, "left")
        self.ib = np.searchsorted(self.t0s, self.ft, "left")

    def sums(self, r: np.ndarray, order: int):
        """S0, S1, S2 over the risk sets of the failure times (weights w·r)."""
        wr = self.w * r
        X = self.X

        def rev(a, o):
            c = np.cumsum(a[o][::-1], axis=0)[::-1]
            return np.concatenate([c, np.zeros((1,) + c.shape[1:])], axis=0)

        S0 = rev(wr, self.ot)[self.ia] - rev(wr, self.ot0)[self.ib]
        S1 = S2 = None
        if order >= 1:
            a1 = X * wr[:, None]
            S1 = rev(a1, self.ot)[self.ia] - rev(a1, self.ot0)[self.ib]
        if order >= 2:
            a2 = X[:, :, None] * X[:, None, :] * wr[:, None, None]
            S2 = rev(a2, self.ot)[self.ia] - rev(a2, self.ot0)[self.ib]
        return S0, S1, S2

    def fsums(self, r: np.ndarray):
        """Sums over the failures tied at each time (for Efron)."""
        f = self.fail
        wr = (self.w * r)[f]
        m = len(self.ft)
        X = self.X[f]
        F0 = np.bincount(self.fidx, weights=wr, minlength=m)
        F1 = np.zeros((m, X.shape[1]))
        np.add.at(F1, self.fidx, X * wr[:, None])
        F2 = np.zeros((m, X.shape[1], X.shape[1]))
        np.add.at(F2, self.fidx, X[:, :, None] * X[:, None, :] * wr[:, None, None])
        return F0, F1, F2


def _cox_fun(strata: list[_Stratum], ties: str, k: int):
    def fun(b, todo):
        ll = 0.0
        g = np.zeros(k)
        H = np.zeros((k, k))
        for sr in strata:
            xb = sr.X @ b
            r = np.exp(xb)
            S0, S1, S2 = sr.sums(r, 0 if todo == 0 else 2)
            ll += float(np.sum((sr.w * xb)[sr.fail]))
            if ties == "efron":
                F0, F1, F2 = sr.fsums(r)
                cnt = np.bincount(sr.fidx, weights=sr.w[sr.fail], minlength=len(sr.ft))
                for j in range(len(sr.ft)):
                    m = int(round(cnt[j]))
                    for l_ in range(m):
                        c = l_ / m
                        s0 = S0[j] - c * F0[j]
                        ll -= np.log(s0)
                        if todo:
                            s1 = S1[j] - c * F1[j]
                            s2 = S2[j] - c * F2[j]
                            a = s1 / s0
                            g -= a
                            H -= s2 / s0 - np.outer(a, a)
                if todo:
                    g += sr.sx.sum(axis=0)
            else:
                ll -= float(np.sum(sr.D * np.log(S0)))
                if todo:
                    a = S1 / S0[:, None]
                    g += sr.sx.sum(axis=0) - (sr.D[:, None] * a).sum(axis=0)
                    H -= np.einsum("j,jkl->kl", sr.D / S0, S2) - np.einsum("j,jk,jl->kl", sr.D, a, a)
        if todo == 0:
            return ll, None, None
        return ll, g, H
    return fun


def _cox_score_resid(strata_obs, b, ties):
    """Efficient score residuals per observation (Breslow form; with Efron ties the
    risk-set means of the Breslow form are used: VERIFY)."""
    out = []
    for sr, idx in strata_obs:
        r = np.exp(sr.X @ b)
        S0, S1, _ = sr.sums(r, 1)
        xbar = S1 / S0[:, None]
        haz = sr.D / S0                       # baseline hazard increments
        # cumulative sums over failure times of haz and haz*xbar, evaluated at t and t0
        ch = np.concatenate([[0.0], np.cumsum(haz)])
        chx = np.concatenate([np.zeros((1, sr.X.shape[1])), np.cumsum(haz[:, None] * xbar, axis=0)])
        jt = np.searchsorted(sr.ft, sr.t, "right")
        jt0 = np.searchsorted(sr.ft, sr.t0, "right")
        A = ch[jt] - ch[jt0]
        AX = chx[jt] - chx[jt0]
        U = -r[:, None] * (sr.X * A[:, None] - AX)
        f = sr.fail
        jf = np.searchsorted(sr.ft, sr.t[f])
        U[f] += sr.X[f] - xbar[jf]
        out.append((idx, U))
    return out


def cmd_stcox(s: "Session", args: str) -> None:
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, "cox", args)
        return
    spec = {"nohr": 4, "breslow": 6, "efron": 5, "exactm": 6, "exactp": 6, "strata": 3, "nolog": 4,
            "log": 3, "noshow": 6, "iterate": 4, "tolerance": 3, "ltolerance": 4, "nrtolerance": 4,
            "estimate": 3, "basesurv": 6, "basechazard": 6, "basehc": 5, "noheader": 6}
    stx = STSample(s, args, spec)
    o = stx.o
    out = s.output
    if o.get("exactm") or o.get("exactp"):
        raise StataError(198, "exactm and exactp ties are not available in OpenDTA")   # OpenDTA
    ties = "efron" if o.get("efron") else "breslow"
    if not o.get("noshow"):
        header(s, stx.info)
    for nm in stx.notes:
        out.write(f"note: {nm} omitted because of collinearity\n", "text")
    X = stx.X
    k = X.shape[1]
    if k == 0:
        raise StataError(102, "too few variables specified")   # VERIFY: stcox, estimate
    strata_var = []
    if isinstance(o.get("strata"), str):
        from ..core.varlist import resolve_name
        strata_var = [resolve_name(s.data, v) for v in s.expand_varlist(o["strata"])]
    sgrp = np.zeros(stx.n, dtype=int)
    if strata_var:
        keys = [tuple(s.data.get(v).data[stx.mask][i] for v in strata_var) for i in range(stx.n)]
        u = sorted(set(keys))
        sgrp = np.array([u.index(kk) for kk in keys])
    strata = []
    strata_obs = []
    for h in np.unique(sgrp):
        idx = np.nonzero(sgrp == h)[0]
        sr = _Stratum(stx.t[idx], stx.t0[idx], stx.d[idx], X[idx], stx.w[idx])
        strata.append(sr)
        strata_obs.append((sr, idx))
    tied = any(np.any(np.bincount(sr.fidx, weights=sr.w[sr.fail], minlength=len(sr.ft)) > 1)
               for sr in strata)
    fun = _cox_fun(strata, ties, k)
    maxiter = int(o["iterate"]) if o.get("iterate") else 16000
    res = maximize(fun, np.zeros(k), maxiter=maxiter,
                   tol=float(o["tolerance"]) if o.get("tolerance") else 1e-6,
                   ltol=float(o["ltolerance"]) if o.get("ltolerance") else 1e-7,
                   nrtol=float(o["nrtolerance"]) if o.get("nrtolerance") else 1e-5)
    ll0 = res.log[0][1]
    clusters, clustvar = _cluster_ids(stx)
    robust = stx.smp.vce in ("robust", "cluster") or stx.wtype == "pweight"
    lab = "log pseudolikelihood" if robust else "log likelihood"
    if not o.get("nolog"):
        print_log(out, res, lab)
        out.write("Refining estimates:\n", "text")
        out.write(f"{'Iteration 0:':<14} {lab} = {_g(res.ll):>10}  \n", "text")
    try:
        V = np.linalg.inv(-res.H)
    except np.linalg.LinAlgError:
        V = np.linalg.pinv(-res.H)
    V = (V + V.T) / 2
    n_clust = None
    if robust:
        U = np.zeros((stx.n, k))
        for idx, Ui in _cox_score_resid(strata_obs, res.b, ties):
            U[idx] = Ui
        wU = U * stx.w[:, None] if stx.wtype != "fweight" else U
        if clusters is not None:
            uq, inv = np.unique(clusters, return_inverse=True)
            n_clust = len(uq)
            G = np.zeros((n_clust, k))
            np.add.at(G, inv, wU)
            meat = G.T @ G
        elif stx.wtype == "fweight":
            meat = (U * stx.w[:, None]).T @ U
        else:
            meat = wU.T @ wU
        V = V @ meat @ V
        V = (V + V.T) / 2
    b = res.b
    df_m = k
    if robust:
        chi2type = "Wald"
        try:
            chi2 = float(b @ np.linalg.solve(V, b))
        except np.linalg.LinAlgError:
            chi2 = SYS
    else:
        chi2type = "LR"
        chi2 = 2 * (res.ll - ll0)
    pval = float(chi2_tail(chi2, df_m)) if chi2 < SYS else SYS
    smp = stx.smp
    names = [c.name if not c.omitted else fvars.omitted_name(c.name) for c in smp.cols]
    kk = len(names)
    pos = [i for i, c in enumerate(smp.cols) if c.values is not None]
    bf = np.zeros(kk)
    Vf = np.zeros((kk, kk))
    bf[pos] = b
    Vf[np.ix_(pos, pos)] = V
    bmap = {c.name: bf[i] for i, c in enumerate(smp.cols)}
    semap = {c.name: float(np.sqrt(max(Vf[i, i], 0))) for i, c in enumerate(smp.cols)}
    est = Estimates("cox", "_t", names, bf, Vf, stx.n, "z", None)
    est.rows = rows_from_columns(smp.cols, bmap, semap, s.data, constant=False)
    est.cols, est.constant, est.terms, est.depcomp = smp.cols, False, smp.terms, smp.depcomp
    est.level = smp.level
    est.vcetype = "Robust" if robust else ""
    est.clustvar = clustvar
    est.n_clust = n_clust
    est.robust = robust
    x = _st_scalars(stx)
    x.update({"ll": res.ll, "ll_0": ll0, "chi2": chi2, "p": pval, "df_m": df_m, "chi2type": chi2type})
    est.extra = x
    est.ties = ties
    est.tied = tied
    est.strata_vars = strata_var
    est.hr = not o.get("nohr")
    est.display = _display_cox
    est.predict = _predict_cox
    est.predict_opts = {"hr": 2, "basesurv": 5, "basechazard": 5, "basehc": 5, "csnell": 3,
                        "mgale": 2, "schoenfeld": 3, "scaledsch": 3, "esr": 3, "partial": 4}
    est.estat = {"phtest": (3, _phtest)}
    est.cox = {"strata": strata_obs, "sgrp": sgrp, "mask": stx.mask.copy()}
    from ..commands.matrix import Matrix
    scalars = [("N", x["N"]), ("N_sub", x["N_sub"]), ("N_fail", x["N_fail"]), ("risk", x["risk"]),
               ("ll", res.ll), ("ll_0", ll0), ("chi2", chi2), ("df_m", float(df_m)), ("p", pval),
               ("rank", float(k)),
               ("converged", 1.0 if res.converged else 0.0), ("ic", float(res.iterations))]
    if n_clust is not None:
        scalars.append(("N_clust", float(n_clust)))
    from .regress import cmdline
    macros = [("cmdline", cmdline("stcox", args)), ("cmd", "cox"), ("cmd2", "stcox"), ("depvar", "_t"),
              ("t0", "_t0"), ("ties", "breslow" if ties == "breslow" else "efron"),
              ("chi2type", chi2type), ("vce", "robust" if clusters is None and robust else
                                       "cluster" if robust else "oim"),
              ("vcetype", "Robust" if robust else ""), ("clustvar", clustvar),
              ("strata", " ".join(strata_var)), ("predict", "stcox_p"), ("estat_cmd", "stcox_estat"),
              ("properties", "b V"), ("crittype", "log pseudolikelihood" if robust else "log likelihood")]
    post(s, est, scalars, macros, stx.mask,
         extra_mats=[("gradient", Matrix(res.g.reshape(1, -1), ["r1"], [names[i] for i in pos]))])
    _stcox_basevars(s, est, o)
    _display_cox(s, est, header=not o.get("noheader"))


def _stcox_basevars(s, est, o):
    from ..core.dataset import Variable
    for opt, stat in (("basesurv", "basesurv"), ("basechazard", "basechazard"), ("basehc", "basehc")):
        if isinstance(o.get(opt), str):
            name = o[opt].strip()
            if s.data.has(name):
                raise StataError(110, f"{name} already defined")
            vals, _ = _predict_cox(s, est, stat, np.ones(s.data.nobs, dtype=bool))
            ok = vals < SYS
            v = np.full(s.data.nobs, SYS)
            v[ok] = vals[ok].astype(np.float32)
            s.data.add(Variable(name, "float", v))


def _display_cox(s: "Session", est: Estimates, *, header: bool = True, table: bool = True,
                 level: float | None = None, eform=None, **_kw) -> None:
    out = s.output
    if header:
        if not est.tied:
            title = "Cox regression -- no ties in data"
        else:
            title = f"Cox regression -- {'Breslow' if est.ties == 'breslow' else 'Efron'} method for ties"
        if est.strata_vars:
            title = title   # VERIFY: Stata adds "Stratified" to the title?
        _st_header_block(s, title, est)
    if not table:
        return
    if est.clustvar:
        msg = f"(Std. Err. adjusted for {int(est.n_clust):,} clusters in {est.clustvar})"
        out.write(f"{msg:>78}\n", "text")
    hr = est.hr if eform is None else bool(eform)
    coef_table(s, "_t", est.rows, stat="z", df=None, level=level or est.level, vcetype=est.vcetype,
               coef_title="Haz. Ratio" if hr else "Coef.", eform=hr, footer=not est.strata_vars)
    if est.strata_vars:
        out.write("-" * 78 + "\n", "text")
        out.write(f"{'Stratified by ' + ' '.join(est.strata_vars):>78}\n", "text")


def _baseline(s, est):
    """Baseline cumulative hazard (covariates at zero) at each observation's _t and _t0,
    and the hazard contribution at failure times, by stratum."""
    ds = s.data
    info = require_st(s)
    n = ds.nobs
    H_t = np.full(n, SYS)
    H_t0 = np.full(n, SYS)
    h_c = np.full(n, SYS)
    sgrp_all = np.full(n, -1)
    cox = est.cox
    st_ok = (ds.get("_st").data == 1)
    t_all = ds.get("_t").data.astype(float)
    t0_all = ds.get("_t0").data.astype(float)
    d_all = ds.get("_d").data
    mask = cox["mask"]
    rows = np.nonzero(mask)[0]
    # strata of all st observations: those of the estimation sample, by value of the strata vars
    if est.strata_vars:
        keys_all = [tuple(ds.get(v).data[i] for v in est.strata_vars) for i in range(n)]
        est_keys = sorted({keys_all[i] for i in rows})
        for i in range(n):
            if keys_all[i] in est_keys:
                sgrp_all[i] = est_keys.index(keys_all[i])
    else:
        sgrp_all[:] = 0
    b = est.b[[i for i, c in enumerate(est.cols) if c.values is not None]]
    for h, (sr, idx) in enumerate(cox["strata"]):
        r = np.exp(sr.X @ b)
        S0, _, _ = sr.sums(r, 0)
        if est.ties == "efron":
            F0, _, _ = sr.fsums(r)
            cnt = np.bincount(sr.fidx, weights=sr.w[sr.fail], minlength=len(sr.ft))
            haz = np.array([sum(1.0 / (S0[j] - (l_ / int(round(cnt[j]))) * F0[j])
                                for l_ in range(int(round(cnt[j])))) for j in range(len(sr.ft))])
        else:
            haz = sr.D / S0
        ch = np.concatenate([[0.0], np.cumsum(haz)])
        m = st_ok & (sgrp_all == h) & (t_all < SYS)
        H_t[m] = ch[np.searchsorted(sr.ft, t_all[m], "right")]
        H_t0[m] = ch[np.searchsorted(sr.ft, t0_all[m], "right")]
        mf = m & (d_all == 1)
        j = np.searchsorted(sr.ft, t_all[mf])
        okj = (j < len(sr.ft)) & (sr.ft[np.minimum(j, len(sr.ft) - 1)] == t_all[mf])
        vals = np.where(okj, haz[np.minimum(j, len(haz) - 1)], SYS)
        h_c[mf] = vals
    return H_t, H_t0, h_c, info


def _predict_cox(s: "Session", est: Estimates, stat: str, mask: np.ndarray):
    from .postest import design_all, xb_all
    xb, ok = xb_all(s, est)
    ds = s.data
    if stat in ("", "hr"):
        return np.where(ok, np.exp(xb), SYS), ""
    if stat == "xb":
        return np.where(ok, xb, SYS), ""
    if stat == "stdp":
        X, _ = design_all(s, est)
        sp = np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", X, est.V, X), 0))
        return np.where(ok, sp, SYS), ""
    H_t, H_t0, h_c, info = _baseline(s, est)
    if stat == "basechazard":
        return H_t, ""
    if stat == "basesurv":
        return np.where(H_t < SYS, np.exp(-H_t), SYS), ""
    if stat == "basehc":
        return h_c, ""
    if stat in ("csnell", "mgale"):
        cs = np.where(ok & (H_t < SYS), np.exp(xb) * (H_t - H_t0), SYS)
        d = ds.get("_d").data
        mg = np.where(cs < SYS, d - cs, SYS)
        val = cs if stat == "csnell" else mg
        if info.idvar and not getattr(est, "_partial", False):
            # one value per subject, on its last record (VERIFY: Stata's default without partial)
            ids = ds.get(info.idvar).data
            tt = ds.get("_t").data
            res = np.full(ds.nobs, SYS)
            acc: dict = {}
            last: dict = {}
            for i in np.nonzero(val < SYS)[0]:
                acc[ids[i]] = acc.get(ids[i], 0.0) + val[i]
                if ids[i] not in last or tt[i] >= tt[last[ids[i]]]:
                    last[ids[i]] = i
            for kk, i in last.items():
                res[i] = acc[kk]
            val = res
        return val, ""
    if stat in ("schoenfeld", "scaledsch", "esr"):
        raise StataError(198, f"use estat phtest; predict {stat} with several new variables is not "
                         "available in OpenDTA")   # OpenDTA
    raise StataError(198, f"option {stat} not allowed")


# ---------------------------------------------------------------------------
# estat phtest
# ---------------------------------------------------------------------------

def _schoenfeld(est):
    """Unscaled Schoenfeld residuals of each failure (rows) and their times."""
    b = est.b[[i for i, c in enumerate(est.cols) if c.values is not None]]
    res, times = [], []
    for sr, _idx in est.cox["strata"]:
        r = np.exp(sr.X @ b)
        S0, S1, _ = sr.sums(r, 1)
        xbar = S1 / S0[:, None]
        f = np.nonzero(sr.fail)[0]
        j = np.searchsorted(sr.ft, sr.t[f])
        for i, jj in zip(f, j):
            reps = int(round(sr.w[i]))
            for _ in range(max(reps, 1)):
                res.append(sr.X[i] - xbar[jj])
                times.append(sr.t[i])
    return np.array(res), np.array(times)


def _phtest(s: "Session", est: Estimates, rest: str) -> None:
    from ..lang.syntax import match_options
    from .km import km
    o = match_options(rest.lstrip(", "), {"detail": 1, "log": 3, "km": 2, "rank": 4, "time": 4,
                                          "plot": 4}) if rest.strip(" ,") else {}
    if est.robust:
        pass
    R, tt = _schoenfeld(est)
    if o.get("log"):
        g, tlab = np.log(tt), "Log(t)"
    elif o.get("km"):
        sd = st_data(s, est.cox["mask"])
        k = km(sd.t, sd.t0, sd.d, sd.w)
        j = np.searchsorted(k.times, tt, "right") - 1
        g, tlab = 1 - k.S[j], "1-KM"     # VERIFY: Stata's km scale
    elif o.get("rank"):
        g, tlab = st.rankdata(tt), "Rank(t)"
    else:
        g, tlab = tt.astype(float), "Time"
    xx = g - g.mean()
    pos = [i for i, c in enumerate(est.cols) if c.values is not None]
    V = est.V[np.ix_(pos, pos)]
    nd = R.shape[0]
    names = [est.names[i] for i in pos]
    out = s.output
    out.write("\n      Test of proportional-hazards assumption\n\n", "text")
    out.write(f"      Time:  {tlab}\n", "text")
    out.write("      " + "-" * 64 + "\n", "text")
    if o.get("detail"):
        out.write("                  |       rho            chi2       df       Prob>chi2\n", "text")
        out.write("      " + "-" * 12 + "+" + "-" * 51 + "\n", "text")
        r2 = R @ V * nd
        for j, nm in enumerate(names):
            test = float(xx @ r2[:, j])
            z = test ** 2 / (V[j, j] * nd * float(xx @ xx))
            rho = float(np.corrcoef(xx, r2[:, j])[0, 1])
            p = float(chi2_tail(z, 1))
            out.write(f"      {nm[:11]:<11} |", "text")
            out.write(f"  {rho:>9.5f}      {z:>10.2f}{1:>9}{p:>16.4f}\n", "result")
        out.write("      " + "-" * 12 + "+" + "-" * 51 + "\n", "text")
    else:
        out.write("                  |                   chi2       df       Prob>chi2\n", "text")
        out.write("      " + "-" * 12 + "+" + "-" * 51 + "\n", "text")
    test = xx @ R
    zg = float(test @ V @ test) * nd / float(xx @ xx)
    df = len(names)
    pg = float(chi2_tail(zg, df))
    out.write(f"      {'global test':<11} |", "text")
    out.write(f"{'':>17}{zg:>10.2f}{df:>9}{pg:>16.4f}\n", "result")
    out.write("      " + "-" * 64 + "\n", "text")
    s.r = {"chi2": zg, "df": float(df), "p": pg}


# ---------------------------------------------------------------------------
# streg
# ---------------------------------------------------------------------------

_DISTS = {"exponential": 1, "weibull": 1, "gompertz": 1, "lognormal": 4, "lnormal": 2,
          "loglogistic": 4, "llogistic": 2}


def _dist_name(text: str) -> str:
    t = text.strip().lower()
    for full, k in (("exponential", 1), ("weibull", 1), ("gompertz", 1), ("lognormal", 4),
                    ("lnormal", 2), ("loglogistic", 4), ("llogistic", 2), ("ggamma", 2)):
        if len(t) >= k and full.startswith(t):
            return {"lnormal": "lognormal", "llogistic": "loglogistic"}.get(full, full)
    raise StataError(198, f"distribution({text}) invalid")   # VERIFY


def _lf(dist: str, aft: bool, t, t0, d):
    """Per-observation ll and derivatives with respect to the two linear indexes
    (eta = xb, th = ancillary parameter)."""
    lt = np.log(np.where(t > 0, t, 1.0))
    has0 = t0 > 0
    lt0 = np.log(np.where(has0, t0, 1.0))

    def exp_ph(eta):
        A = t - t0
        e = np.exp(eta)
        l = d * eta - e * A
        return l, d - e * A, -e * A

    def weib_ph(eta, th):
        p = np.exp(th)
        tp = np.where(t > 0, t ** p, 0.0)
        t0p = np.where(has0, t0 ** p, 0.0)
        A = tp - t0p
        B = p * (tp * lt - t0p * lt0)
        C = tp * lt ** 2 - t0p * lt0 ** 2
        e = np.exp(eta)
        l = d * (eta + th + (p - 1) * lt) - e * A
        return (l, d - e * A, d * (1 + p * lt) - e * B, -e * A, -e * B, d * p * lt - e * (B + p * p * C))

    def gomp(eta, gam):
        def G(x):
            gx = gam * x
            small = np.abs(gx) < 1e-4
            E = np.exp(gx)
            with np.errstate(divide="ignore", invalid="ignore"):
                g0 = np.where(small, x + gam * x ** 2 / 2 + gam ** 2 * x ** 3 / 6, (E - 1) / gam)
                g1 = np.where(small, x ** 2 / 2 + gam * x ** 3 / 3 + gam ** 2 * x ** 4 / 8,
                              x * E / gam - (E - 1) / gam ** 2)
                g2 = np.where(small, x ** 3 / 3 + gam * x ** 4 / 4,
                              x * x * E / gam - 2 * x * E / gam ** 2 + 2 * (E - 1) / gam ** 3)
            return g0, g1, g2
        a0, a1, a2 = G(t)
        b0, b1, b2 = G(t0)
        A, A1, A2 = a0 - b0, a1 - b1, a2 - b2
        e = np.exp(eta)
        l = d * (eta + gam * t) - e * A
        return (l, d - e * A, d * t - e * A1, -e * A, -e * A1, -e * A2)

    def locscale(eta, th, kind):
        sig = np.exp(th)

        def parts(z, which):
            if kind == "lognormal":
                if which == "f":
                    return -0.5 * z * z - 0.5 * np.log(2 * np.pi), -z, -np.ones_like(z)
                lS = special.log_ndtr(-z)
                lam = np.exp(st.norm.logpdf(z) - lS)
                return lS, -lam, -lam * (lam - z)
            sp = np.logaddexp(0, z)
            sg = special.expit(z)
            if which == "f":
                return z - 2 * sp, 1 - 2 * sg, -2 * sg * (1 - sg)
            return -sp, -sg, -sg * (1 - sg)

        z = (lt - eta) / sig
        Tf, Tf1, Tf2 = parts(z, "f")
        TS, TS1, TS2 = parts(z, "S")
        T = np.where(d > 0, Tf - th - lt, TS)
        T1 = np.where(d > 0, Tf1, TS1)
        T2 = np.where(d > 0, Tf2, TS2)
        l = T
        l_e = -T1 / sig
        l_t = -z * T1 - (d > 0)
        l_ee = T2 / sig ** 2
        l_et = (T2 * z + T1) / sig
        l_tt = z * T1 + z * z * T2
        if np.any(has0):
            z0 = (lt0 - eta) / sig
            S0, S1, S2 = parts(z0, "S")
            S0, S1, S2 = [np.where(has0, a, 0.0) for a in (S0, S1, S2)]
            l = l - S0
            l_e = l_e + S1 / sig
            l_t = l_t + z0 * S1
            l_ee = l_ee - S2 / sig ** 2
            l_et = l_et - (S2 * z0 + S1) / sig
            l_tt = l_tt - (z0 * S1 + z0 * z0 * S2)
            l_t = np.where(has0, l_t, l_t)
        return l, l_e, l_t, l_ee, l_et, l_tt

    if dist == "exponential":
        def f1(u):
            eta = -u if aft else u
            l, a, b = exp_ph(eta)
            return (l, -a, b) if aft else (l, a, b)
        return f1, 1
    if dist == "weibull":
        def f2(u, th):
            if not aft:
                return weib_ph(u, th)
            p = np.exp(th)
            eta = -p * u
            l, le, lt_, lee, let, ltt = weib_ph(eta, th)
            # chain rule: eta = -p u, d eta/du = -p, d eta/dth = eta
            return (l, -p * le, lt_ + le * eta, p * p * lee,
                    -p * (lee * eta + let + le), ltt + 2 * let * eta + lee * eta * eta + le * eta)
        return f2, 2
    if dist == "gompertz":
        return gomp, 2
    return (lambda u, th: locscale(u, th, dist)), 2


_AUX = {"weibull": "ln_p", "gompertz": "gamma", "lognormal": "ln_sig", "loglogistic": "ln_gam"}
_TITLE = {"exponential": "Exponential", "weibull": "Weibull", "gompertz": "Gompertz",
          "lognormal": "Lognormal", "loglogistic": "Loglogistic"}


def cmd_streg(s: "Session", args: str) -> None:
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, "streg", args)
        return
    spec = {"distribution": 1, "time": 4, "nohr": 4, "tr": 2, "noshow": 6, "nolog": 4, "log": 3,
            "iterate": 4, "tolerance": 3, "ltolerance": 4, "nrtolerance": 4, "noheader": 6,
            "nolrtest": 5}
    stx = STSample(s, args, spec)
    o = stx.o
    if not isinstance(o.get("distribution"), str):
        raise StataError(198, "option distribution() required")   # VERIFY
    dist = _dist_name(o["distribution"])
    if dist == "ggamma":
        raise StataError(198, "distribution(ggamma) is not available in OpenDTA")   # OpenDTA
    aft = bool(o.get("time")) or dist in ("lognormal", "loglogistic")
    if dist == "gompertz" and o.get("time"):
        raise StataError(198, "Gompertz has no accelerated failure-time form")   # VERIFY
    out = s.output
    if not o.get("noshow"):
        header(s, stx.info)
    for nm in stx.notes:
        out.write(f"note: {nm} omitted because of collinearity\n", "text")
    t, t0, d, w = stx.t, stx.t0, stx.d, stx.w
    if np.any((t <= 0) & (d > 0)):
        raise StataError(498, "failures at time 0")   # VERIFY
    X = np.column_stack([stx.X, np.ones(stx.n)]) if stx.X.size else np.ones((stx.n, 1))
    kx = X.shape[1]
    lf, npar = _lf(dist, aft, t, t0, d)
    k = kx + (npar - 1)

    def fun_of(Xm):
        kk = Xm.shape[1] + (npar - 1)

        def fun(b, todo):
            u = Xm @ b[:Xm.shape[1]]
            if npar == 1:
                l, l1, l11 = lf(u)
                ll = float(np.sum(w * l))
                if not np.isfinite(ll):
                    return -np.inf, None, None
                if todo == 0:
                    return ll, None, None
                g = Xm.T @ (w * l1)
                H = (Xm * (w * l11)[:, None]).T @ Xm
                return ll, g, H
            th = b[-1]
            with np.errstate(all="ignore"):
                l, le, lt_, lee, let, ltt = lf(u, th)
            ll = float(np.sum(w * l))
            if not np.isfinite(ll):
                return -np.inf, None, None
            if todo == 0:
                return ll, None, None
            g = np.concatenate([Xm.T @ (w * le), [np.sum(w * lt_)]])
            H = np.zeros((kk, kk))
            H[:-1, :-1] = (Xm * (w * lee)[:, None]).T @ Xm
            H[:-1, -1] = H[-1, :-1] = Xm.T @ (w * let)
            H[-1, -1] = np.sum(w * ltt)
            return ll, g, H
        return fun

    # starting values: constant-only model from moments of ln t of the failures
    rate = float(np.sum(w * d) / np.sum(w * (t - t0)))
    if dist == "exponential":
        c0 = np.array([np.log(rate) * (-1 if aft else 1)])
    elif dist == "weibull":
        c0 = np.array([np.log(rate) * (-1 if aft else 1), 0.0])
    elif dist == "gompertz":
        c0 = np.array([np.log(rate), 0.0])
    else:
        # VERIFY: Stata's starting values are undocumented; ln of the mean time at risk per failure
        c0 = np.array([-np.log(rate), 0.0])
    maxiter = int(o["iterate"]) if o.get("iterate") else 16000
    kw = dict(maxiter=maxiter, tol=float(o["tolerance"]) if o.get("tolerance") else 1e-6,
              ltol=float(o["ltolerance"]) if o.get("ltolerance") else 1e-7,
              nrtol=float(o["nrtolerance"]) if o.get("nrtolerance") else 1e-5)
    robust = stx.smp.vce in ("robust", "cluster") or stx.wtype == "pweight"
    lab = "log pseudolikelihood" if robust else "log likelihood"
    res0 = maximize(fun_of(np.ones((stx.n, 1))), c0, **kw)
    b0 = np.zeros(k)
    b0[kx - 1] = res0.b[0]
    if npar == 2:
        b0[-1] = res0.b[1]
    res = maximize(fun_of(X), b0, **kw)
    if not o.get("nolog"):
        out.write("Fitting constant-only model:\n\n", "text")
        print_log(out, res0, lab)
        out.write("\nFitting full model:\n\n", "text")
        print_log(out, res, lab)
    try:
        V = np.linalg.inv(-res.H)
    except np.linalg.LinAlgError:
        V = np.linalg.pinv(-res.H)
    V = (V + V.T) / 2
    clusters, clustvar = _cluster_ids(stx)
    n_clust = None
    if robust:
        # scores per observation: derivative of ll_i with respect to the parameters
        u = X @ res.b[:kx]
        if npar == 1:
            _l, l1, _ = lf(u)
            S = X * l1[:, None]
        else:
            _l, le, lt_, *_ = lf(u, res.b[-1])
            S = np.column_stack([X * le[:, None], lt_])
        S = S * w[:, None] if stx.wtype != "fweight" else S
        if clusters is not None:
            uq, inv = np.unique(clusters, return_inverse=True)
            n_clust = len(uq)
            G = np.zeros((n_clust, k))
            np.add.at(G, inv, S)
            meat = G.T @ G * n_clust / (n_clust - 1)
        else:
            meat = (S * (w[:, None] if stx.wtype == "fweight" else 1)).T @ S * stx.n / (stx.n - 1)
        V = V @ meat @ V
        V = (V + V.T) / 2
    nslope = kx - 1
    if robust:
        chi2type = "Wald"
        sl = list(range(nslope))
        chi2 = float(res.b[sl] @ np.linalg.solve(V[np.ix_(sl, sl)], res.b[sl])) if sl else 0.0
    else:
        chi2type = "LR"
        chi2 = 2 * (res.ll - res0.ll)
    pval = float(chi2_tail(chi2, nslope)) if nslope and chi2 < SYS else SYS
    smp = stx.smp
    xnames = [c.name if not c.omitted else fvars.omitted_name(c.name) for c in smp.cols] + ["_cons"]
    pos = [i for i, c in enumerate(smp.cols) if c.values is not None] + [len(smp.cols)]
    kk = len(xnames) + (npar - 1)
    bf, Vf = np.zeros(kk), np.zeros((kk, kk))
    full_pos = pos + ([kk - 1] if npar == 2 else [])
    bf[full_pos] = res.b
    Vf[np.ix_(full_pos, full_pos)] = V
    names = xnames + (["_cons"] if npar == 2 else [])
    eqs = ["_t"] * len(xnames) + ([_AUX[dist]] if npar == 2 else [])
    bmap = {c.name: bf[i] for i, c in enumerate(smp.cols)}
    semap = {c.name: float(np.sqrt(max(Vf[i, i], 0))) for i, c in enumerate(smp.cols)}
    bmap["_cons"], semap["_cons"] = bf[len(xnames) - 1], float(np.sqrt(max(Vf[len(xnames) - 1, len(xnames) - 1], 0)))
    est = Estimates("streg", "_t", names, bf, Vf, stx.n, "z", None)
    est.eqnames = eqs
    est.rows = rows_from_columns(smp.cols, bmap, semap, s.data, constant=True)
    est.cols, est.constant, est.terms, est.depcomp = smp.cols, True, smp.terms, smp.depcomp
    est.level = smp.level
    est.vcetype = "Robust" if robust else ""
    est.clustvar, est.n_clust, est.robust = clustvar, n_clust, robust
    x = _st_scalars(stx)
    x.update({"ll": res.ll, "ll_0": res0.ll, "chi2": chi2, "p": pval, "df_m": nslope, "chi2type": chi2type})
    est.extra = x
    est.dist, est.aft = dist, aft
    est.aux = (float(bf[-1]), float(np.sqrt(max(Vf[-1, -1], 0)))) if npar == 2 else None
    est.eform_kind = "tr" if (aft and o.get("tr")) else ("hr" if (not aft and not o.get("nohr")) else "")
    est.display = _display_streg
    est.predict = _predict_streg
    est.predict_opts = {"median": 3, "mean": 4, "time": 4, "lntime": 3, "hazard": 3, "hr": 2,
                        "surv": 4, "csurv": 5, "csnell": 3, "mgale": 2}
    from .regress import cmdline
    scalars = [("N", x["N"]), ("N_sub", x["N_sub"]), ("N_fail", x["N_fail"]), ("risk", x["risk"]),
               ("ll", res.ll), ("ll_0", res0.ll), ("chi2", chi2), ("df_m", float(nslope)), ("p", pval),
               ("k_eq", float(npar)), ("k_aux", float(npar - 1)), ("rank", float(k)),
               ("converged", 1.0 if res.converged else 0.0), ("ic", float(res.iterations))]
    if dist == "weibull":
        scalars.append(("aux_p", float(np.exp(bf[-1]))))
    elif dist == "gompertz":
        scalars.append(("gamma", float(bf[-1])))
    elif dist == "lognormal":
        scalars.append(("sigma", float(np.exp(bf[-1]))))
    elif dist == "loglogistic":
        scalars.append(("gamma", float(np.exp(bf[-1]))))
    if n_clust is not None:
        scalars.append(("N_clust", float(n_clust)))
    cmd_map = {"exponential": "ereg", "weibull": "weibull", "gompertz": "gompertz",
               "lognormal": "lnormal", "loglogistic": "llogistic"}
    macros = [("cmdline", cmdline("streg", args)), ("cmd", cmd_map[dist]), ("cmd2", "streg"),
              ("depvar", "_t"), ("t0", "_t0"), ("distribution", dist),
              ("frm2", "time" if aft else "hazard"), ("title", _TITLE[dist] + " regression"),
              ("chi2type", chi2type), ("vce", "robust" if robust else "oim"),
              ("vcetype", "Robust" if robust else ""), ("clustvar", clustvar),
              ("predict", "streg_p"), ("properties", "b V")]
    est.cmd = "streg"
    post(s, est, scalars, macros, stx.mask)
    _display_streg(s, est, header=not o.get("noheader"))


def _display_streg(s: "Session", est: Estimates, *, header: bool = True, table: bool = True,
                   level: float | None = None, **_kw) -> None:
    out = s.output
    form = "accelerated failure-time form" if est.aft else "log relative-hazard form"
    if header:
        _st_header_block(s, f"{_TITLE[est.dist]} regression -- {form}", est)
    if not table:
        return
    if est.clustvar:
        msg = f"(Std. Err. adjusted for {int(est.n_clust):,} clusters in {est.clustvar})"
        out.write(f"{msg:>78}\n", "text")
    lev = level or est.level
    eform = est.eform_kind
    title = {"hr": "Haz. Ratio", "tr": "Tm. Ratio"}.get(eform, "Coef.")   # VERIFY: Time Ratio label
    rows = list(est.rows)
    if est.aux is not None:
        b, se = est.aux
        z = st.norm.ppf(1 - (1 - lev / 100) / 2)
        aux = _AUX[est.dist]
        rows.append(CoefRow("sep"))
        rows.append(CoefRow("coef", f"/{aux}", b, se, name=f"{aux}:_cons", noeform=True))
        rows.append(CoefRow("sep"))
        lo, hi = b - z * se, b + z * se
        if est.dist == "weibull":
            p = np.exp(b)
            rows.append(CoefRow("aux", "p", p, p * se, eq=f"{float(np.exp(lo))!r} {float(np.exp(hi))!r}"))
            rows.append(CoefRow("aux", "1/p", 1 / p, se / p, eq=f"{float(np.exp(-hi))!r} {float(np.exp(-lo))!r}"))
        elif est.dist == "lognormal":
            sg = np.exp(b)
            rows.append(CoefRow("aux", "sigma", sg, sg * se, eq=f"{float(np.exp(lo))!r} {float(np.exp(hi))!r}"))
        elif est.dist == "loglogistic":
            gm = np.exp(b)
            rows.append(CoefRow("aux", "gamma", gm, gm * se, eq=f"{float(np.exp(lo))!r} {float(np.exp(hi))!r}"))
        else:
            rows.pop()      # gompertz: /gamma only
    coef_table(s, "_t", rows, stat="z", df=None, level=lev, vcetype=est.vcetype,
               coef_title=title, eform=bool(eform))


def _predict_streg(s: "Session", est: Estimates, stat: str, mask: np.ndarray):
    from .postest import xb_all
    xb, ok = xb_all(s, est)
    ds = s.data
    dist, aft = est.dist, est.aft
    th = est.aux[0] if est.aux else 0.0
    t = ds.get("_t").data.astype(float)
    t0 = ds.get("_t0").data.astype(float)
    okt = ok & (t < SYS)

    def surv(x):
        x = np.where(x > 0, x, 0.0)
        if dist == "exponential":
            lam = np.exp(-xb if aft else xb)
            return np.exp(-lam * x)
        if dist == "weibull":
            p = np.exp(th)
            lam = np.exp(-p * xb) if aft else np.exp(xb)
            return np.exp(-lam * x ** p)
        if dist == "gompertz":
            g = th
            G = (np.exp(g * x) - 1) / g if abs(g) > 1e-12 else x
            return np.exp(-np.exp(xb) * G)
        sig = np.exp(th)
        with np.errstate(divide="ignore"):
            z = (np.log(np.where(x > 0, x, 1e-300)) - xb) / sig
        if dist == "lognormal":
            return np.where(x > 0, st.norm.sf(z), 1.0)
        return np.where(x > 0, special.expit(-z), 1.0)

    if stat == "xb":
        return np.where(ok, xb, SYS), ""
    if stat in ("", "median"):
        note = "" if stat else "(option median time assumed; predicted median time)"   # VERIFY
        if dist == "exponential":
            lam = np.exp(-xb if aft else xb)
            v = np.log(2) / lam
        elif dist == "weibull":
            p = np.exp(th)
            lam = np.exp(-p * xb) if aft else np.exp(xb)
            v = (np.log(2) / lam) ** (1 / p)
        elif dist == "gompertz":
            g = th
            with np.errstate(invalid="ignore", divide="ignore"):
                v = np.log1p(g * np.log(2) / np.exp(xb)) / g if abs(g) > 1e-12 else np.log(2) / np.exp(xb)
        else:
            v = np.exp(xb)
        if not stat and not aft:
            pass
        return np.where(ok & np.isfinite(v), v, SYS), note
    if stat in ("hr",):
        if aft:
            raise StataError(198, "option hr not allowed after accelerated failure-time models")   # VERIFY
        return np.where(ok, np.exp(xb), SYS), ""
    if stat == "surv":
        return np.where(okt, surv(t) / surv(t0), SYS), ""
    if stat == "csurv":
        return np.where(okt, surv(t), SYS), ""
    if stat in ("csnell", "mgale"):
        cs = -np.log(surv(t)) + np.log(surv(t0))
        if stat == "csnell":
            return np.where(okt, cs, SYS), ""
        d = ds.get("_d").data
        return np.where(okt, d - cs, SYS), ""
    if stat == "hazard":
        if dist == "exponential":
            v = np.exp(-xb if aft else xb) * np.ones_like(t)
        elif dist == "weibull":
            p = np.exp(th)
            lam = np.exp(-p * xb) if aft else np.exp(xb)
            v = lam * p * np.where(t > 0, t, np.nan) ** (p - 1)
        elif dist == "gompertz":
            v = np.exp(xb + th * t)
        else:
            eps = 1e-6 * np.maximum(t, 1e-8)
            S1, S2 = surv(t), surv(t + eps)
            v = (np.log(S1) - np.log(S2)) / eps
        return np.where(okt & np.isfinite(v), v, SYS), ""
    if stat in ("mean", "time", "lntime"):
        if stat == "lntime" and dist in ("lognormal", "loglogistic"):
            return np.where(ok, xb, SYS), ""
        raise StataError(198, f"option {stat} is not available in OpenDTA for this distribution")   # OpenDTA
    raise StataError(198, f"option {stat} not allowed")
