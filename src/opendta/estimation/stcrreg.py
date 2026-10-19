"""Competing-risks regression: stcrreg ([ST] stcrreg).

The Fine-Gray subhazard model: the Breslow partial likelihood in which subjects that had a
competing event stay in the risk sets with weight S_c(t)/S_c(t_j), S_c being the Kaplan-Meier
survivor function of the censoring times. The variance is the sandwich estimator with the scores
of the weighted likelihood plus the term for the estimation of S_c (Fine and Gray 1999). Layout
details carry the VERIFY mark (compat/do/1901_survival.do)."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.stats import chi2_tail
from ..lang.words import parse_numlist
from .km import km
from .ml import maximize, print_log
from .ml import options as ml_options
from .results import Estimates, coef_table, post, rows_from_columns
from .st import _fail_text
from .stcox import _Stratum
from .stfit import STSample, cluster_note, fmt, place
from .vce import bread, sandwich

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


class FineGray:
    """Data of the weighted partial likelihood: t, t0, d (failure of interest), comp (competing
    event at the exit of the record), cens (censoring at the exit of the record), weights and
    the design matrix X."""

    def __init__(self, t, t0, d, comp, cens, w, X):
        self.t, self.t0, self.d, self.comp, self.cens, self.w, self.X = t, t0, d, comp, cens, w, X
        self.km = km(t, t0, cens.astype(np.float64), w)
        self.risk = _Stratum(t, t0, d, X, w)
        self.ft = self.risk.ft
        self.sc_ft = self.sc_left(self.ft)
        c = np.flatnonzero(comp)
        c = c[np.argsort(t[c], kind="stable")]
        self.ic = c
        self.tc = t[c]
        self.invsc = 1.0 / self.sc_left(self.tc) if len(c) else np.zeros(0)
        self.cnt = np.searchsorted(self.tc, self.ft, "left")      # competing before each failure time

    def sc_left(self, x: np.ndarray) -> np.ndarray:
        """Kaplan-Meier survivor function of the censoring times just before x."""
        k = self.km
        j = np.searchsorted(k.times, x, "left") - 1
        return np.where(j >= 0, k.S[np.maximum(j, 0)], 1.0)

    def _csum(self, a: np.ndarray) -> np.ndarray:
        return np.concatenate([np.zeros((1,) + a.shape[1:]), np.cumsum(a, axis=0)], axis=0)

    def sums(self, r: np.ndarray, order: int):
        """S0, S1, S2 over the weighted risk sets of the failure times."""
        S0, S1, S2 = self.risk.sums(r, order)
        if len(self.ic):
            wr = (self.w * r)[self.ic] * self.invsc
            Xc = self.X[self.ic]
            S0 = S0 + self.sc_ft * self._csum(wr)[self.cnt]
            if order >= 1:
                S1 = S1 + self.sc_ft[:, None] * self._csum(wr[:, None] * Xc)[self.cnt]
            if order >= 2:
                S2 = S2 + self.sc_ft[:, None, None] * self._csum(
                    wr[:, None, None] * Xc[:, :, None] * Xc[:, None, :])[self.cnt]
        return S0, S1, S2

    def fun(self, b, todo):
        rk = self.risk
        xb = self.X @ b
        S0, S1, S2 = self.sums(np.exp(xb), 0 if todo == 0 else 2)
        ll = float(np.sum((self.w * xb)[rk.fail])) - float(np.sum(rk.D * np.log(S0)))
        if todo == 0:
            return ll, None, None
        a = S1 / S0[:, None]
        g = rk.sx.sum(axis=0) - (rk.D[:, None] * a).sum(axis=0)
        H = -(np.einsum("j,jkl->kl", rk.D / S0, S2) - np.einsum("j,jk,jl->kl", rk.D, a, a))
        return ll, g, H

    def baseline(self, b) -> np.ndarray:
        """Breslow cumulative baseline subhazard at the failure times."""
        S0, _, _ = self.sums(np.exp(self.X @ b), 0)
        return np.cumsum(self.risk.D / S0)

    def scores(self, b) -> np.ndarray:
        """Score of each record: the partial-likelihood term plus the one for the estimation of
        the censoring distribution (VERIFY: derived from Fine and Gray 1999)."""
        t, t0, d, w, X = self.t, self.t0, self.d, self.w, self.X
        rk = self.risk
        n, k = X.shape
        r = np.exp(X @ b)
        S0, S1, _ = self.sums(r, 1)
        a = S1 / S0[:, None]
        c = rk.D / S0
        Cc, Ca = self._csum(c), self._csum(c[:, None] * a)
        lo = np.searchsorted(self.ft, t0, "right")
        hi = np.searchsorted(self.ft, t, "right")
        eta = -r[:, None] * (X * (Cc[hi] - Cc[lo])[:, None] - (Ca[hi] - Ca[lo]))
        f = np.flatnonzero(d > 0)
        eta[f] += X[f] - a[np.searchsorted(self.ft, t[f])]
        cs = c * self.sc_ft
        Ccs, Cas = self._csum(cs), self._csum(cs[:, None] * a)
        ic = self.ic
        if len(ic):
            T0 = Ccs[-1] - Ccs[hi[ic]]
            T1 = Cas[-1] - Cas[hi[ic]]
            eta[ic] -= (r[ic] * self.invsc_of(ic))[:, None] * (X[ic] * T0[:, None] - T1)
        psi = self._psi(r, a, c, cs, Ccs, Cas)
        return eta + psi

    def invsc_of(self, rows: np.ndarray) -> np.ndarray:
        return 1.0 / self.sc_left(self.t[rows])

    def _psi(self, r, a, c, cs, Ccs, Cas) -> np.ndarray:
        t, t0, X, w = self.t, self.t0, self.X, self.w
        n, k = X.shape
        kk = self.km
        grid = np.flatnonzero(kk.d > 0)
        psi = np.zeros((n, k))
        if not len(grid):
            return psi
        s = kk.times[grid]
        ybar, dN = kk.n[grid], kk.d[grid]
        wr = (w * r)[self.ic] * self.invsc
        Xc = self.X[self.ic]
        P0, P1 = self._csum(wr), self._csum(wr[:, None] * Xc)
        up = np.searchsorted(self.tc, s, "right")
        lw = np.searchsorted(self.tc, s, "left")
        nf = np.searchsorted(self.ft, s, "right")
        q = P1[up] * (Ccs[-1] - Ccs[nf])[:, None] - P0[up][:, None] * (Cas[-1] - Cas[nf])
        at = np.searchsorted(self.ft, s, "left")
        edge = (at < len(self.ft)) & (self.ft[np.minimum(at, len(self.ft) - 1)] == s)
        m = at[edge]
        q[edge] += cs[m][:, None] * (P1[lw[edge]] - P0[lw[edge]][:, None] * a[m])
        G = self._csum(q * (dN / ybar ** 2)[:, None])
        hi = np.searchsorted(s, t, "right")
        lo = np.searchsorted(s, t0, "right")
        psi -= G[hi] - G[lo]
        ce = np.flatnonzero(self.cens)
        psi[ce] += q[np.searchsorted(s, t[ce])] / ybar[np.searchsorted(s, t[ce])][:, None]
        return psi


def _cox_start(t, t0, d, w, X) -> np.ndarray:
    """Starting values: the Cox (Breslow) estimates with the competing events censored, which is
    where Stata's log starts (compat 1901)."""
    k = X.shape[1]
    times = np.unique(t[d > 0])
    risks = [(t0 < tf) & (t >= tf) for tf in times]
    fails = [(t == tf) & (d > 0) for tf in times]

    def fun(b, todo):
        eta = X @ b
        ll, g, H = 0.0, np.zeros(k), np.zeros((k, k))
        for risk, fail in zip(risks, fails):
            e = w[risk] * np.exp(eta[risk])
            S0 = e.sum()
            Xr = X[risk]
            S1 = Xr.T @ e
            nd = float(w[fail].sum())
            ll += float((w[fail] * eta[fail]).sum()) - nd * np.log(S0)
            if todo:
                g += X[fail].T @ w[fail] - nd * S1 / S0
                H -= nd * ((Xr * e[:, None]).T @ Xr / S0 - np.outer(S1, S1) / S0 ** 2)
        return ll, g, H
    try:
        return maximize(fun, np.zeros(k)).b
    except (np.linalg.LinAlgError, StataError):
        return np.zeros(k)


def _competing(s: "Session", spec: str, mask: np.ndarray):
    """(indicator over the sample rows of mask, variable, text of the event) for compete()."""
    ds = s.data
    var, sep, rest = spec.partition("==")
    name = s.expand_varlist(var.strip())[0]
    v = ds.get(name)
    if v.is_string:
        raise StataError(109, "type mismatch")
    x = v.data[mask]
    if sep:
        vals = [float(z) for z in parse_numlist(rest.strip())]
        return np.isin(x, vals), name, f"{name} == {' '.join(f'{z:g}' for z in vals)}"
    return (x != 0) & (x < SYS), name, f"{name} != 0 & {name} < ."


_SPEC = {"compete": 4, "noshr": 4, "noshow": 6, "noheader": 6, "nolog": 4, "iterate": 4,
         "tolerance": 3, "ltolerance": 4, "nrtolerance": 4, "noadjust": 5, "tvc": 3, "texp": 4,
         "offset": 3, "constraints": 1, "collinear": 4, "nodisplay": 6, "notable": 6}


def cmd_stcrreg(s: "Session", args: str) -> None:
    if not args.strip() or args.strip().startswith(","):
        from .postest import replay
        replay(s, "stcrreg", args)
        return
    stx = STSample(s, args, _SPEC)
    o = stx.o
    for k in ("tvc", "texp", "offset", "constraints", "collinear", "noadjust"):
        if o.get(k):
            raise StataError(198, f"option {k} is not available yet in OpenDTA")
    if not o.get("compete"):
        raise StataError(198, "option compete() required")
    if str(o.get("vce", "")).strip().startswith(("oim", "ols")):
        raise StataError(198, "vcetype 'oim' not allowed")   # VERIFY
    comp_rec, cname, ctext = _competing(s, str(o["compete"]), stx.mask)
    if np.any(comp_rec & (stx.d > 0)):
        raise StataError(198, "competing events must be censored in stset")   # VERIFY
    stx.notes_and_header(s)
    ids = stx.sd.ids
    t, n = stx.t, stx.n
    last = np.ones(n, dtype=bool)
    comp = comp_rec.copy()
    if ids is not None:
        # only the last record of a subject carries the event (VERIFY: ties in t)
        last[:] = False
        for sid in np.unique(ids):
            j = np.flatnonzero(ids == sid)
            last[j[np.lexsort((j, t[j]))[-1]]] = True
        comp = comp_rec & last
    cens = last & (stx.d == 0) & ~comp
    X = stx.X
    k = X.shape[1]
    if k == 0:
        raise StataError(102, "too few variables specified")
    model = FineGray(t, stx.t0, stx.d, comp, cens, stx.w, X)
    res = maximize(model.fun, _cox_start(t, stx.t0, stx.d, stx.w, X), **ml_options(o))
    ll0 = res.log[0][1]
    if not o.get("nolog"):
        print_log(s.output, res, "log pseudolikelihood")
    b = res.b
    U = model.scores(b)
    cl, clustvar = stx.clusters()
    if stx.smp.vce == "cluster":
        clusters, clustvar = cl, stx.smp.clustvar
    elif ids is not None:
        clusters, clustvar = ids, stx.info.idvar
    else:
        clusters, clustvar = None, ""
    V, n_clust = sandwich(bread(res.H), U, stx.w, fweight=stx.fweight,
                          n=float(stx.w.sum()) if stx.fweight else float(n), clusters=clusters)
    try:
        chi2 = float(b @ np.linalg.solve(V, b))
    except np.linalg.LinAlgError:
        chi2 = SYS
    pval = float(chi2_tail(chi2, k)) if chi2 < SYS else SYS
    smp = stx.smp
    names, bf, Vf, _pos, bmap, semap = place(smp, b, V, constant=False)
    est = Estimates("stcrreg", "_t", names, bf, Vf, n, "z", None)
    est.rows = rows_from_columns(smp.cols, bmap, semap, s.data, constant=False)
    est.cols, est.constant, est.terms, est.depcomp = smp.cols, False, smp.terms, smp.depcomp
    est.level = smp.level
    est.vcetype, est.robust = "Robust", True
    est.clustvar, est.n_clust = clustvar, n_clust
    wv = stx.w if stx.fweight else np.ones(n)
    sc = stx.scalars()
    sc.update({"N_fail": float(np.sum(stx.d * wv)), "N_compete": float(np.sum(comp * wv)),
               "ll": res.ll, "chi2": chi2, "p": pval, "df_m": k})
    sc["N_censor"] = sc["N_sub"] - sc["N_fail"] - sc["N_compete"]
    est.extra = sc
    est.shr = not o.get("noshr")
    est.fevent, est.crevent, est.cname = _fail_text(stx.info), ctext, cname
    est.display = _display
    est.predict = _predict
    est.predict_opts = {"xb": 2, "shr": 3}
    est.fg = {"ft": model.ft, "H0": model.baseline(b), "b": b,
              "cols": [c for c in smp.cols if c.values is not None]}
    _post(s, est, args, stx, res, clusters, o, ll0)
    _display(s, est, header=not o.get("noheader"), table=not o.get("notable"))


def _post(s: "Session", est: Estimates, args: str, stx: STSample, res, clusters, o, ll0) -> None:
    from .regress import cmdline
    x = est.extra
    k = len(res.b)
    scalars = [("rank", float(k)), ("N", x["N"]), ("ic", float(res.iterations)), ("k", float(k)),
               ("k_eq", 1.0), ("k_dv", 1.0), ("converged", 1.0 if res.converged else 0.0),
               ("rc", 0.0), ("ll", res.ll), ("k_eq_model", 1.0)]
    if est.n_clust is not None:
        scalars.append(("N_clust", float(est.n_clust)))
    scalars += [("chi2", x["chi2"]), ("p", x["p"]), ("df_m", float(k)), ("fmult", 0.0),
                ("crmult", 0.0), ("fnz", 0.0), ("crnz", 0.0), ("N_sub", x["N_sub"]),
                ("N_censor", x["N_censor"]), ("N_compete", x["N_compete"]), ("N_fail", x["N_fail"])]
    mainvars = " ".join(c.name for c in est.fg["cols"])
    macros = [("cmdline", cmdline("stcrreg", args)), ("title", "Competing-risks regression"),
              ("mainvars", mainvars), ("fevent", est.fevent), ("crevent", est.crevent),
              ("compete", str(o["compete"])), ("predict", "stcrreg_p"),
              ("marginsnotok", "SCOres ESR DFBeta SCHoenfeld"), ("cmd", "stcrreg"),
              ("clustvar", est.clustvar), ("vce", "cluster" if est.clustvar else "robust"),
              ("vcetype", "Robust"), ("chi2type", "Wald"), ("opt", "moptimize"),
              ("user", "stcrr_lf"), ("ml_method", "d2"), ("technique", "nr"), ("which", "max"),
              ("depvar", "_t"), ("properties", "b V")]   # VERIFY: the e() lists
    post(s, est, scalars, macros, stx.mask)


def _display(s: "Session", est: Estimates, *, header: bool = True, table: bool = True,
             level: float | None = None, eform=None, **_kw) -> None:
    out = s.output
    x = est.extra
    if header:
        left = ["Competing-risks regression", "", f"{'Failure event  :'} {est.fevent}",
                f"{'Competing event:'} {est.crevent}", "", "", f"Log pseudolikelihood = {fmt(x['ll'], '%10.0g'):>10}"]
        right = [("No. of obs", fmt(x["N"], "%12.0gc")), ("No. of subjects", fmt(x["N_sub"], "%12.0gc")),
                 ("No. failed", fmt(x["N_fail"], "%12.0gc")), ("No. competing", fmt(x["N_compete"], "%12.0gc")),
                 ("No. censored", fmt(x["N_censor"], "%12.0gc")), None,
                 (f"Wald chi2({int(x['df_m'])})", f"{x['chi2']:.2f}" if x["chi2"] < SYS else "."),
                 ("Prob > chi2", f"{x['p']:.4f}" if x["p"] < SYS else ".")]
        out.write("\n", "text")
        lines = [(left[0], right[0]), ("", right[1]), (left[2], right[2]), (left[3], right[3]),
                 ("", right[4]), ("", None), ("", right[6]), (left[6], right[7])]
        for lft, rgt in lines:
            if rgt is None:
                out.write(lft + "\n", "text")
                continue
            out.write(f"{lft:<49}{rgt[0]:<17}=", "text")
            out.write(f"{rgt[1]:>11}\n", "result")
        out.write("\n", "text")
    if not table:
        return
    cluster_note(s, est)
    shr = est.shr if eform is None else bool(eform)
    coef_table(s, "_t", est.rows, stat="z", df=None, level=level or est.level, vcetype="Robust",
               coef_title="SHR" if shr else "Coef.", eform=shr, footer=True)


def _predict(s: "Session", est: Estimates, stat: str, mask: np.ndarray):
    from .postest import xb_all
    xb, ok = xb_all(s, est)
    if stat in ("", "xb"):
        note = "" if stat else "(option xb assumed; fitted values)"   # VERIFY
        return np.where(ok, xb, SYS), note
    if stat == "shr":
        return np.where(ok, np.exp(xb), SYS), ""
    raise StataError(198, f"option {stat} not allowed")
