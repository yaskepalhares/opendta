"""Cox proportional-hazards model: stcox, its predictions and estat phtest ([ST] stcox,
[ST] stcox postestimation, [ST] stcox PH-assumption tests).

The Cox partial likelihood with the Breslow (default) or Efron approximation for ties; delayed
entry and strata enter through the risk sets. The robust variance uses the efficient score
residuals (Lin and Wei 1989). Layout details carry the VERIFY mark
(compat/do/0902_stcox_streg.do).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.grouping import group_ids
from ..core.stats import chi2_tail
from ..lazy import lazy
from . import fvars
from .ml import maximize, print_log, refine
from .ml import options as ml_options
from .results import Estimates, coef_table, post, rows_from_columns
from .st import require_st, st_data
from .stfit import STSample, cluster_note, fmt, header_block, place
from .vce import bread, sandwich

# imported on first use (start-up time; see opendta/lazy.py)
st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


# ---------------------------------------------------------------------------
# partial likelihood
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


def _coefs(est: Estimates) -> np.ndarray:
    """The estimated coefficients (without the omitted columns)."""
    return est.b[[i for i, c in enumerate(est.cols) if c.values is not None]]


def _strata_of_all(ds, strata_vars: list[str], mask: np.ndarray) -> np.ndarray:
    """Stratum index (as in the fit) of every observation whose strata values occur in the
    estimation sample; -1 for the others. Without strata, 0 everywhere."""
    if not strata_vars:
        return np.zeros(ds.nobs, dtype=np.int64)
    gid = group_ids(ds, strata_vars)[0]
    present = np.unique(gid[mask])
    h = np.searchsorted(present, gid)
    ok = (h < len(present)) & (present[np.minimum(h, len(present) - 1)] == gid)
    return np.where(ok, h, -1)


def _efron_hazard(sr: "_Stratum", r: np.ndarray, S0: np.ndarray) -> np.ndarray:
    """Baseline hazard increments with Efron's correction for the failures tied at each time."""
    F0, _, _ = sr.fsums(r)
    cnt = np.rint(np.bincount(sr.fidx, weights=sr.w[sr.fail], minlength=len(sr.ft))).astype(int)
    return np.array([sum(1.0 / (S0[j] - (l_ / cnt[j]) * F0[j]) for l_ in range(cnt[j]))
                     for j in range(len(sr.ft))])


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
    if o.get("exactm") or o.get("exactp"):
        raise StataError(198, "exactm and exactp ties are not available in OpenDTA")   # OpenDTA
    ties = "efron" if o.get("efron") else "breslow"
    stx.notes_and_header(s)
    X = stx.X
    k = X.shape[1]
    if k == 0:
        raise StataError(102, "too few variables specified")   # VERIFY: stcox, estimate
    strata_var = s.expand_varlist(o["strata"]) if isinstance(o.get("strata"), str) else []
    sgrp = group_ids(s.data, strata_var, stx.mask)[0][stx.mask] if strata_var else np.zeros(stx.n, dtype=int)
    strata_obs = []
    for h in np.unique(sgrp):
        idx = np.flatnonzero(sgrp == h)
        strata_obs.append((_Stratum(stx.t[idx], stx.t0[idx], stx.d[idx], X[idx], stx.w[idx]), idx))
    strata = [sr for sr, _ in strata_obs]
    tied = any(np.any(np.bincount(sr.fidx, weights=sr.w[sr.fail], minlength=len(sr.ft)) > 1)
               for sr in strata)
    fun = _cox_fun(strata, ties, k)
    # VERIFY: stcox stops one iteration earlier than ml when g H^-1 g' is already tiny (stratified
    # model in compat 0902: 1.9e-8 stops, 7.4e-7 does not)
    res = maximize(fun, np.zeros(k), nr_stop=1e-7, **ml_options(o))
    refine(fun, res)              # "Refining estimates", outside the iteration log
    ll0 = res.log[0][1]
    robust = stx.robust
    lab = "log pseudolikelihood" if robust else "log likelihood"
    if not o.get("nolog"):
        print_log(s.output, res, lab)
        s.output.write("Refining estimates:\n", "text")
        s.output.write(f"{'Iteration 0:':<14} {lab} = {fmt(res.ll):>10}  \n", "text")
    V = bread(res.H)
    clusters, clustvar = stx.clusters()
    n_clust = None
    if robust:
        U = np.zeros((stx.n, k))
        for idx, Ui in _cox_score_resid(strata_obs, res.b, ties):
            U[idx] = Ui
        # compat 0902: the standard errors include the factor n/(n-1)
        V, n_clust = sandwich(V, U, stx.w, fweight=stx.fweight,
                              n=float(stx.w.sum()) if stx.fweight else float(stx.n), clusters=clusters)
    b = res.b
    if robust:
        chi2type = "Wald"
        try:
            chi2 = float(b @ np.linalg.solve(V, b))
        except np.linalg.LinAlgError:
            chi2 = SYS
    else:
        chi2type = "LR"
        chi2 = 2 * (res.ll - ll0)
    pval = float(chi2_tail(chi2, k)) if chi2 < SYS else SYS
    smp = stx.smp
    names, bf, Vf, _pos, bmap, semap = place(smp, b, V, constant=False)
    est = Estimates("cox", "_t", names, bf, Vf, stx.n, "z", None)
    est.rows = rows_from_columns(smp.cols, bmap, semap, s.data, constant=False)
    est.cols, est.constant, est.terms, est.depcomp = smp.cols, False, smp.terms, smp.depcomp
    est.level = smp.level
    est.vcetype = "Robust" if robust else ""
    est.clustvar, est.n_clust, est.robust = clustvar, n_clust, robust
    x = stx.scalars()
    x.update({"ll": res.ll, "ll_0": ll0, "chi2": chi2, "p": pval, "df_m": k, "chi2type": chi2type})
    est.extra = x
    est.ties, est.tied, est.strata_vars = ties, tied, strata_var
    est.hr = not o.get("nohr")
    est.display = _display_cox
    est.predict = _predict_cox
    est.predict_opts = {"hr": 2, "basesurv": 5, "basechazard": 5, "basehc": 5, "csnell": 3,
                        "mgale": 2, "schoenfeld": 3, "scaledsch": 3, "esr": 3, "partial": 4}
    est.estat = {"phtest": (3, _phtest)}
    est.cox = {"strata": strata_obs, "sgrp": sgrp, "mask": stx.mask.copy()}
    _post_cox(s, est, args, stx, res, clusters)
    _stcox_basevars(s, est, o)
    _display_cox(s, est, header=not o.get("noheader"))


def _post_cox(s: "Session", est: Estimates, args: str, stx: STSample, res, clusters) -> None:
    """e() in the order Stata lists it (compat 0902)."""
    from .regress import cmdline
    from .results import datasignature
    x = est.extra
    k = len(res.b)
    ll0 = x["ll_0"]
    robust = est.robust
    scalars = [("N", x["N"]), ("ll_0", ll0), ("ll", res.ll), ("chi2", x["chi2"]),
               ("r2_p", 1 - res.ll / ll0 if ll0 else SYS),
               ("converged", 1.0 if res.converged else 0.0), ("rank", float(k)), ("df_m", float(k))]
    if est.n_clust is not None:
        scalars.append(("N_clust", float(est.n_clust)))   # VERIFY: position
    scalars += [("N_sub", x["N_sub"]), ("N_fail", x["N_fail"]), ("risk", x["risk"])]
    sigvars = ["_t", "_t0", "_d"] + fvars.base_vars(stx.smp.terms)
    macros = [("cmdline", cmdline("stcox", args)), ("datasignaturevars", " ".join(sigvars)),
              ("datasignature", datasignature(s, sigvars, stx.mask)), ("cmd2", "stcox"),
              ("marginsnotok", "CSNell DEViance DFBeta ESR LDisplace LMax MGale SCAledsch SCHoenfeld "
                               "SCores"),
              ("footnote", "stcox_footnote"), ("predict", "stcox_p"), ("estat_cmd", "stcox_estat"),
              # VERIFY: position of clustvar and vcetype
              ("clustvar", est.clustvar), ("vcetype", "Robust" if robust else ""),
              ("vce", "robust" if clusters is None and robust else "cluster" if robust else "oim"),
              ("chi2type", x["chi2type"]), ("method", est.ties), ("ties", est.ties),
              ("strata", " ".join(est.strata_vars)),   # VERIFY: position
              ("t0", "_t0"), ("depvar", "_t"), ("cmd", "cox"), ("properties", "b V")]
    post(s, est, scalars, macros, stx.mask)


def _stcox_basevars(s, est, o):
    from ..core.dataset import Variable
    for opt in ("basesurv", "basechazard", "basehc"):
        if isinstance(o.get(opt), str):
            name = o[opt].strip()
            if s.data.has(name):
                raise StataError(110, f"{name} already defined")
            vals, _ = _predict_cox(s, est, opt, np.ones(s.data.nobs, dtype=bool))
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
            title = "Stratified Cox regr." + title[len("Cox regression"):]
        header_block(s, title, est)
    if not table:
        return
    cluster_note(s, est)
    hr = est.hr if eform is None else bool(eform)
    coef_table(s, "_t", est.rows, stat="z", df=None, level=level or est.level, vcetype=est.vcetype,
               coef_title="Haz. Ratio" if hr else "Coef.", eform=hr, footer=not est.strata_vars)
    if est.strata_vars:
        out.write("-" * 78 + "\n", "text")
        out.write(f"{'Stratified by ' + ' '.join(est.strata_vars):>78}\n", "text")


def _baseline(s, est):
    """Baseline cumulative hazard (covariates at zero) at each observation's _t and _t0,
    the hazard contribution at failure times and the baseline survivor at _t, by stratum."""
    ds = s.data
    info = require_st(s)
    n = ds.nobs
    H_t = np.full(n, SYS)
    H_t0 = np.full(n, SYS)
    S_t = np.full(n, SYS)
    h_c = np.full(n, SYS)
    cox = est.cox
    st_ok = (ds.get("_st").data == 1)
    t_all = ds.get("_t").data.astype(float)
    t0_all = ds.get("_t0").data.astype(float)
    d_all = ds.get("_d").data
    # strata of all st observations: those of the estimation sample, by value of the strata vars
    sgrp_all = _strata_of_all(ds, est.strata_vars, cox["mask"])
    b = _coefs(est)
    for h, (sr, idx) in enumerate(cox["strata"]):
        r = np.exp(sr.X @ b)
        S0, _, _ = sr.sums(r, 0)
        haz = _efron_hazard(sr, r, S0) if est.ties == "efron" else sr.D / S0
        ch = np.concatenate([[0.0], np.cumsum(haz)])
        # baseline survivor: Kalbfleisch-Prentice product of the conditional survival
        # probabilities alpha_j, which solve sum_D r/(1 - alpha^r) = sum_R r ([ST] stcox
        # postestimation, Methods and formulas)
        alpha = np.ones(len(sr.ft))
        fr = r[sr.fail]
        fw = sr.w[sr.fail]
        for j in range(len(sr.ft)):
            rj, wj = fr[sr.fidx == j], fw[sr.fidx == j]
            if len(rj) == 1 and wj[0] == 1:
                alpha[j] = max(1 - rj[0] / S0[j], 0.0) ** (1 / rj[0])
                continue
            lo_, hi_ = 0.0, 1.0
            for _ in range(200):
                mid = (lo_ + hi_) / 2
                f_ = np.sum(wj * rj / (1 - mid ** rj)) - S0[j]
                if f_ > 0:
                    hi_ = mid
                else:
                    lo_ = mid
            alpha[j] = (lo_ + hi_) / 2
        cs_ = np.concatenate([[1.0], np.cumprod(alpha)])
        m = st_ok & (sgrp_all == h) & (t_all < SYS)
        S_t[m] = cs_[np.searchsorted(sr.ft, t_all[m], "right")]
        H_t[m] = ch[np.searchsorted(sr.ft, t_all[m], "right")]
        H_t0[m] = ch[np.searchsorted(sr.ft, t0_all[m], "right")]
        mf = m & (d_all == 1)
        j = np.searchsorted(sr.ft, t_all[mf])
        okj = (j < len(sr.ft)) & (sr.ft[np.minimum(j, len(sr.ft) - 1)] == t_all[mf])
        vals = np.where(okj, haz[np.minimum(j, len(haz) - 1)], SYS)
        h_c[mf] = vals
    return H_t, H_t0, h_c, S_t, info


def _predict_cox(s: "Session", est: Estimates, stat: str, mask: np.ndarray):
    from .postest import design_all, xb_all
    xb, ok = xb_all(s, est)
    ds = s.data
    if stat in ("", "hr"):
        return np.where(ok, np.exp(xb), SYS), ("" if stat else "(option hr assumed; relative hazard)")
    if stat == "xb":
        return np.where(ok, xb, SYS), ""
    if stat == "stdp":
        X, _ = design_all(s, est)
        sp = np.sqrt(np.maximum(np.einsum("ij,jk,ik->i", X, est.V, X), 0))
        return np.where(ok, sp, SYS), ""
    H_t, H_t0, h_c, S_t, info = _baseline(s, est)
    if stat == "basechazard":
        return H_t, ""
    if stat == "basesurv":
        return S_t, ""
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
    b = _coefs(est)
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
            out.write(f"{rho:>13.5f}{z:>13.2f}{1:>9}{p:>15.4f}\n", "result")
        out.write("      " + "-" * 12 + "+" + "-" * 51 + "\n", "text")
    else:
        out.write("                  |                      chi2       df       Prob>chi2\n", "text")
        out.write("      " + "-" * 12 + "+" + "-" * 51 + "\n", "text")
    test = xx @ R
    zg = float(test @ V @ test) * nd / float(xx @ xx)
    df = len(names)
    pg = float(chi2_tail(zg, df))
    out.write(f"      {'global test':<11} |", "text")
    out.write(f"{zg:>26.2f}{df:>9}{pg:>15.4f}\n", "result")
    out.write("      " + "-" * 64 + "\n", "text")
    s.r = {"chi2": zg, "df": float(df), "p": pg}


# ---------------------------------------------------------------------------
# streg
# ---------------------------------------------------------------------------
