"""Parametric survival models: streg with the exponential, Weibull, Gompertz, lognormal and
loglogistic distributions, and its predictions ([ST] streg, [ST] streg postestimation).

Each distribution is a small class that knows its log likelihood (with delayed entry, [ST]
streg, Methods and formulas), its survivor function, median and hazard, the labels of its
ancillary parameter and the names Stata stores in e(); cmd_streg itself has no per-distribution
branches. Layout details carry the VERIFY mark (compat/do/0902_stcox_streg.do).
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.stats import chi2_tail
from ..lazy import lazy
from .ml import ilog, maximize, print_log
from .ml import options as ml_options
from .results import CoefRow, Estimates, coef_table, post, rows_from_columns
from .stfit import STSample, cluster_note, header_block, place
from .vce import bread, sandwich

# imported on first use (start-up time; see opendta/lazy.py)
special = lazy("scipy.special")
st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS
LN2 = float(np.log(2))


# ---------------------------------------------------------------------------
# distributions
# ---------------------------------------------------------------------------

class Dist:
    """A streg distribution. The log likelihood returns, per observation, l and its first and
    second derivatives with respect to the linear index u and the ancillary parameter th:
    (l, l_u, l_th, l_uu, l_uth, l_thth); without an ancillary parameter the th terms are 0."""

    name = title = cmd = predict_sub = user = ""
    aux = ""                   # ancillary equation (ln_p, gamma...), "" for none
    metrics = ("hazard", "time")
    fits_time_metric = True    # False: fitted in the hazard metric, transformed afterwards

    def loglik(self, t, t0, d, aft: bool):
        raise NotImplementedError

    def start(self, rate: float, aft: bool) -> list[float]:
        """Constant-only starting values from the failure rate D/T."""
        return [float(np.log(rate)), 0.0]

    def aux_scalar(self, th: float) -> tuple[str, float] | None:
        return None

    def aux_rows(self, b: float, se: float, lo: float, hi: float) -> list[CoefRow]:
        """Transformed ancillary parameters shown under /aux (p and 1/p for the Weibull...)."""
        return []

    def surv(self, x, xb, th: float, aft: bool):
        raise NotImplementedError

    def median(self, xb, th: float, aft: bool):
        raise NotImplementedError

    def hazard(self, t, xb, th: float, aft: bool):
        eps = 1e-6 * np.maximum(t, 1e-8)
        return (np.log(self.surv(t, xb, th, aft)) - np.log(self.surv(t + eps, xb, th, aft))) / eps


def _exp_aux(label: str, b, se, lo, hi, sign: int = 1) -> CoefRow:
    v = np.exp(sign * b)
    a, z = (lo, hi) if sign > 0 else (-hi, -lo)
    return CoefRow("aux", label, v, v * se, eq=f"{float(np.exp(a))!r} {float(np.exp(z))!r}")


def _pos_log(x):
    return np.log(np.where(x > 0, x, 1.0))


class _PH(Dist):
    """Proportional-hazards distributions: S(x) = exp(-exp(eta) H0(x)); in the time metric the
    linear index enters as eta(xb)."""

    def eta(self, xb, th: float, aft: bool):
        return xb

    def surv(self, x, xb, th, aft):
        x = np.where(x > 0, x, 0.0)
        return np.exp(-np.exp(self.eta(xb, th, aft)) * self.cumhaz0(x, th))


class Exponential(_PH):
    name, title, cmd, predict_sub, user = "exponential", "Exponential", "ereg", "ereg_p", "ereg_lf0"

    def loglik(self, t, t0, d, aft):
        A = t - t0
        z = np.zeros_like(t)
        sgn = -1.0 if aft else 1.0          # time metric: eta = -u

        def lf(u, th):
            e = np.exp(sgn * u)
            return d * sgn * u - e * A, sgn * (d - e * A), z, -e * A, z, z
        return lf

    def start(self, rate, aft):
        return [float(np.log(rate)) * (-1 if aft else 1)]

    def eta(self, xb, th, aft):
        return -xb if aft else xb

    def cumhaz0(self, x, th):
        return x

    def median(self, xb, th, aft):
        return LN2 / np.exp(self.eta(xb, th, aft))

    def hazard(self, t, xb, th, aft):
        return np.exp(self.eta(xb, th, aft)) * np.ones_like(t)


class Weibull(_PH):
    name, title, cmd, predict_sub, user = "weibull", "Weibull", "weibull", "weibul_p", "weib_lf0"
    aux = "ln_p"
    # fitted in the hazard metric and transformed afterwards: Stata's iteration log is the same in
    # both metrics (compat 0902)
    fits_time_metric = False

    def loglik(self, t, t0, d, aft):
        lt, lt0, has0 = _pos_log(t), _pos_log(t0), t0 > 0

        def lf(eta, th):
            p = np.exp(th)
            tp = np.where(t > 0, t ** p, 0.0)
            t0p = np.where(has0, t0 ** p, 0.0)
            A = tp - t0p
            B = p * (tp * lt - t0p * lt0)
            C = tp * lt ** 2 - t0p * lt0 ** 2
            e = np.exp(eta)
            return (d * (eta + th + (p - 1) * lt) - e * A, d - e * A, d * (1 + p * lt) - e * B,
                    -e * A, -e * B, d * p * lt - e * (B + p * p * C))
        return lf

    def aux_scalar(self, th):
        return ("aux_p", float(np.exp(th)))

    def aux_rows(self, b, se, lo, hi):
        return [_exp_aux("p", b, se, lo, hi), _exp_aux("1/p", b, se, lo, hi, sign=-1)]

    def eta(self, xb, th, aft):
        return -np.exp(th) * xb if aft else xb

    def cumhaz0(self, x, th):
        return x ** np.exp(th)

    def median(self, xb, th, aft):
        return (LN2 / np.exp(self.eta(xb, th, aft))) ** (1 / np.exp(th))

    def hazard(self, t, xb, th, aft):
        p = np.exp(th)
        return np.exp(self.eta(xb, th, aft)) * p * np.where(t > 0, t, np.nan) ** (p - 1)


class Gompertz(_PH):
    name, title, cmd, predict_sub, user = "gompertz", "Gompertz", "gompertz", "gomper_p", "gomp_lf0"
    aux = "gamma"
    metrics = ("hazard",)

    @staticmethod
    def _G(x, gam):
        """Integral of exp(gam s) over (0, x) and its first two derivatives in gam (series near
        gam x = 0)."""
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

    def loglik(self, t, t0, d, aft):
        def lf(eta, gam):
            (a0, a1, a2), (b0, b1, b2) = self._G(t, gam), self._G(t0, gam)
            A, A1, A2 = a0 - b0, a1 - b1, a2 - b2
            e = np.exp(eta)
            return d * (eta + gam * t) - e * A, d - e * A, d * t - e * A1, -e * A, -e * A1, -e * A2
        return lf

    def aux_scalar(self, th):
        return ("gamma", float(th))

    def cumhaz0(self, x, th):
        return (np.exp(th * x) - 1) / th if abs(th) > 1e-12 else x

    def median(self, xb, th, aft):
        with np.errstate(invalid="ignore", divide="ignore"):
            if abs(th) > 1e-12:
                return np.log1p(th * LN2 / np.exp(xb)) / th
            return LN2 / np.exp(xb)

    def hazard(self, t, xb, th, aft):
        return np.exp(xb + th * t)


class _LocationScale(Dist):
    """ln t = u + sigma e, with sigma = exp(th): always in the time metric."""

    metrics = ("time",)

    def start(self, rate, aft):
        # VERIFY: Stata's starting values are undocumented; ln of the mean time at risk per failure
        return [-float(np.log(rate)), 0.0]

    def parts(self, z, density: bool):
        """ln f(z) or ln S(z) of the standard error distribution, and two derivatives in z."""
        raise NotImplementedError

    def loglik(self, t, t0, d, aft):
        lt, lt0, has0, fail = _pos_log(t), _pos_log(t0), t0 > 0, d > 0

        def lf(u, th):
            sig = np.exp(th)
            z = (lt - u) / sig
            Tf, Tf1, Tf2 = self.parts(z, True)
            TS, TS1, TS2 = self.parts(z, False)
            T1 = np.where(fail, Tf1, TS1)
            T2 = np.where(fail, Tf2, TS2)
            l = np.where(fail, Tf - th - lt, TS)
            l_u, l_t = -T1 / sig, -z * T1 - fail
            l_uu, l_ut, l_tt = T2 / sig ** 2, (T2 * z + T1) / sig, z * T1 + z * z * T2
            if np.any(has0):
                # delayed entry: minus ln S at the entry time
                z0 = (lt0 - u) / sig
                S0, S1, S2 = (np.where(has0, a, 0.0) for a in self.parts(z0, False))
                l, l_u, l_t = l - S0, l_u + S1 / sig, l_t + z0 * S1
                l_uu, l_ut = l_uu - S2 / sig ** 2, l_ut - (S2 * z0 + S1) / sig
                l_tt = l_tt - (z0 * S1 + z0 * z0 * S2)
            return l, l_u, l_t, l_uu, l_ut, l_tt
        return lf

    def aux_scalar(self, th):
        return (self.aux_scalar_name, float(np.exp(th)))

    def aux_rows(self, b, se, lo, hi):
        return [_exp_aux(self.aux_scalar_name, b, se, lo, hi)]

    def surv(self, x, xb, th, aft):
        with np.errstate(divide="ignore"):
            z = (np.log(np.where(x > 0, x, 1e-300)) - xb) / np.exp(th)
        return np.where(x > 0, self.sf(z), 1.0)

    def median(self, xb, th, aft):
        return np.exp(xb)


class Lognormal(_LocationScale):
    name, title, cmd, predict_sub, user = "lognormal", "Lognormal", "lnormal", "lnorma_p", "lnorm_lf0"
    aux, aux_scalar_name = "ln_sig", "sigma"

    def parts(self, z, density):
        if density:
            return -0.5 * z * z - 0.5 * np.log(2 * np.pi), -z, -np.ones_like(z)
        lS = special.log_ndtr(-z)
        lam = np.exp(st.norm.logpdf(z) - lS)
        return lS, -lam, -lam * (lam - z)

    def sf(self, z):
        return st.norm.sf(z)


class Loglogistic(_LocationScale):
    name, title, cmd, predict_sub, user = "loglogistic", "Loglogistic", "llogistic", "llogis_p", "llog_lf0"
    aux, aux_scalar_name = "ln_gam", "gamma"

    def parts(self, z, density):
        sp, sg = np.logaddexp(0, z), special.expit(z)
        if density:
            return z - 2 * sp, 1 - 2 * sg, -2 * sg * (1 - sg)
        return -sp, -sg, -sg * (1 - sg)

    def sf(self, z):
        return special.expit(-z)


DISTS: dict[str, Dist] = {d.name: d for d in (Exponential(), Weibull(), Gompertz(), Lognormal(),
                                              Loglogistic())}


def _dist(text: str) -> Dist:
    t = text.strip().lower()
    for full, k in (("exponential", 1), ("weibull", 1), ("gompertz", 1), ("lognormal", 4),
                    ("lnormal", 2), ("loglogistic", 4), ("llogistic", 2), ("ggamma", 2)):
        if len(t) >= k and full.startswith(t):
            if full == "ggamma":
                from .streg_ext import GGAMMA
                return GGAMMA
            return DISTS[{"lnormal": "lognormal", "llogistic": "loglogistic"}.get(full, full)]
    raise StataError(198, f"distribution({text}) invalid")   # VERIFY


# ---------------------------------------------------------------------------
# streg
# ---------------------------------------------------------------------------

def _model(lf, X, w, lnt_c: float, k_aux: int):
    """The maximize() function of a design X: log likelihood, gradient and Hessian in (beta, th),
    without the th row and column when the distribution has no ancillary parameter."""
    kx = X.shape[1]
    k = kx + k_aux

    def fun(b, todo):
        th = b[kx] if k_aux else 0.0
        with np.errstate(all="ignore"):
            l, lu, lt, luu, lut, ltt = lf(X @ b[:kx], th)
        ll = float(np.sum(w * l)) + lnt_c
        if not np.isfinite(ll):
            return -np.inf, None, None
        if todo == 0:
            return ll, None, None
        g = np.concatenate([X.T @ (w * lu), [np.sum(w * lt)]])
        H = np.zeros((kx + 1, kx + 1))
        H[:-1, :-1] = (X * (w * luu)[:, None]).T @ X
        H[:-1, -1] = H[-1, :-1] = X.T @ (w * lut)
        H[-1, -1] = np.sum(w * ltt)
        return ll, g[:k], H[:k, :k]
    return fun


def _scores(lf, X, b, k_aux: int) -> np.ndarray:
    kx = X.shape[1]
    _l, lu, lt, *_ = lf(X @ b[:kx], b[kx] if k_aux else 0.0)
    return np.column_stack([X * lu[:, None], lt])[:, :kx + k_aux]


def _to_time_metric(b: np.ndarray, V: np.ndarray, kx: int):
    """Weibull: gamma = -beta/p, with ln_p right after the coefficients and the parameters that
    follow it unchanged (delta method for V)."""
    p = float(np.exp(b[kx]))
    beta = b[:kx]
    J = np.eye(len(b))
    J[:kx, :kx] = -np.eye(kx) / p
    J[:kx, kx] = beta / p
    V = J @ V @ J.T
    return np.concatenate([-beta / p, b[kx:]]), (V + V.T) / 2


def cmd_streg(s: "Session", args: str) -> None:
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, "streg", args)
        return
    spec = {"distribution": 1, "time": 4, "nohr": 4, "tr": 2, "noshow": 6, "nolog": 4, "log": 3,
            "iterate": 4, "tolerance": 3, "ltolerance": 4, "nrtolerance": 4, "noheader": 6,
            "nolrtest": 5, "frailty": 3, "shared": 2}
    stx = STSample(s, args, spec)
    o = stx.o
    if not isinstance(o.get("distribution"), str):
        raise StataError(198, "option distribution() required")   # VERIFY
    dist = _dist(o["distribution"])
    aft = bool(o.get("time")) or dist.metrics == ("time",)
    if aft and "time" not in dist.metrics:
        raise StataError(198, f"{dist.title} has no accelerated failure-time form")   # VERIFY
    stx.notes_and_header(s)
    t, t0, d, w = stx.t, stx.t0, stx.d, stx.w
    if np.any((t <= 0) & (d > 0)):
        raise StataError(498, "failures at time 0")   # VERIFY
    if dist.name == "ggamma" or o.get("frailty") or o.get("shared"):
        from .streg_ext import fit_extended
        fit_extended(s, stx, dist, o, args, aft)
        return
    X = np.column_stack([stx.X, np.ones(stx.n)])
    kx = X.shape[1]
    k_aux = 1 if dist.aux else 0
    k = kx + k_aux
    fit_aft = aft and dist.fits_time_metric
    lf = dist.loglik(t, t0, d, fit_aft)
    # Stata reports the log likelihood of ln t rather than of t: the densities are multiplied by
    # t, which adds sum(d ln t) (compat 0902: the same constant in every distribution)
    lnt_c = float(np.sum(w * d * _pos_log(t)))
    rate = float(np.sum(w * d) / np.sum(w * (t - t0)))
    kw = ml_options(o)
    res0 = maximize(_model(lf, np.ones((stx.n, 1)), w, lnt_c, k_aux), np.array(dist.start(rate, fit_aft)), **kw)
    b0 = np.zeros(k)
    b0[kx - 1:] = res0.b
    res = maximize(_model(lf, X, w, lnt_c, k_aux), b0, **kw)
    robust = stx.robust
    lab = "log pseudolikelihood" if robust else "log likelihood"
    if not o.get("nolog"):
        # the exponential constant-only model has a closed form: Stata does not show its fit
        if k_aux:
            s.output.write("Fitting constant-only model:\n\n", "text")
            print_log(s.output, res0, lab)
            s.output.write("\nFitting full model:\n\n", "text")
        print_log(s.output, res, lab)
    V = bread(res.H)
    clusters, clustvar = stx.clusters()
    n_clust = None
    if robust:
        V, n_clust = sandwich(V, _scores(lf, X, res.b, k_aux), w, fweight=stx.fweight,
                              n=float(w.sum()) if stx.fweight else float(stx.n), clusters=clusters)
    if aft and not fit_aft:
        res.b, V = _to_time_metric(res.b, V, kx)
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
    names, bf, Vf, full_pos, bmap, semap = place(smp, res.b, V, constant=True, aux=k_aux)
    est = Estimates("streg", "_t", names, bf, Vf, stx.n, "z", None)
    est.eqnames = ["_t"] * (len(names) - k_aux) + [dist.aux] * k_aux
    est.rows = rows_from_columns(smp.cols, bmap, semap, s.data, constant=True)
    est.cols, est.constant, est.terms, est.depcomp = smp.cols, True, smp.terms, smp.depcomp
    est.level = smp.level
    est.vcetype = "Robust" if robust else ""
    est.clustvar, est.n_clust, est.robust = clustvar, n_clust, robust
    x = stx.scalars()
    x.update({"ll": res.ll, "ll_0": res0.ll, "chi2": chi2, "p": pval, "df_m": nslope, "chi2type": chi2type})
    est.extra = x
    est.dist, est.aft = dist.name, aft
    est.aux = (float(bf[-1]), float(np.sqrt(max(Vf[-1, -1], 0)))) if k_aux else None
    est.eform_kind = "tr" if (aft and o.get("tr")) else ("hr" if (not aft and not o.get("nohr")) else "")
    est.display = _display_streg
    est.predict = _predict_streg
    est.predict_opts = {"median": 3, "mean": 4, "time": 4, "lntime": 3, "hazard": 3, "hr": 2,
                        "surv": 4, "csurv": 5, "csnell": 3, "mgale": 2}
    est.cmd = "streg"
    _post_streg(s, est, args, dist, res, res0, full_pos, stx.mask)
    _display_streg(s, est, header=not o.get("noheader"))


def _post_streg(s: "Session", est: Estimates, args: str, dist: Dist, res, res0,
                full_pos: list[int], mask: np.ndarray) -> None:
    """e() in the order Stata lists it after streg, distribution(weibull) (compat 0902); the other
    distributions follow the same pattern (VERIFY, as their predict_sub and user)."""
    from ..commands.matrix import Matrix
    from .regress import cmdline
    x = est.extra
    aft = est.aft
    k = len(full_pos)
    scalars = [("rank", float(k)), ("N", x["N"])]
    if est.n_clust is not None:
        scalars.append(("N_clust", float(est.n_clust)))   # VERIFY: position
    scalars += [("ic", float(res.iterations)), ("k", float(k)), ("k_eq", 2.0 if dist.aux else 1.0),
                ("k_dv", 1.0), ("converged", 1.0 if res.converged else 0.0), ("rc", 0.0),
                ("ll", res.ll), ("k_eq_model", 1.0), ("ll_0", res0.ll), ("df_m", float(x["df_m"])),
                ("chi2", x["chi2"]), ("p", x["p"])]
    if dist.aux:
        scalars += [dist.aux_scalar(float(est.b[-1])), ("k_aux", 1.0)]
    scalars += [("N_sub", x["N_sub"]), ("N_fail", x["N_fail"]), ("risk", x["risk"])]
    macros = [("cmdline", cmdline("streg", args)), ("predict", "streg_p"),
              ("predict_sub", dist.predict_sub), ("cmd2", "streg"), ("title", dist.title + " regression"),
              ("title2", "accelerated failure-time form" if aft else "log relative-hazard form"),
              ("cmd", dist.cmd), ("dead", "_d"), ("frm2", "time" if aft else "hazard"),
              ("stcurve", "stcurve"), ("t0", "_t0"), ("chi2type", x["chi2type"]), ("opt", "moptimize"),
              ("clustvar", est.clustvar), ("vcetype", est.vcetype),
              ("vce", "cluster" if est.n_clust is not None else "robust" if est.robust else "oim"),
              ("user", dist.user), ("ml_method", "e2"), ("technique", "nr"), ("which", "max"),
              ("depvar", "_t"), ("properties", "b V")]
    grad = np.zeros((1, len(est.names)))
    if res.g is not None:
        grad[0, full_pos] = res.g
    mats = [("ilog", Matrix(ilog(res), ["r1"], [f"c{i + 1}" for i in range(20)])),
            ("gradient", Matrix(grad, ["r1"], list(est.names)))]
    post(s, est, scalars, macros, mask, extra_mats=mats)


def _display_streg(s: "Session", est: Estimates, *, header: bool = True, table: bool = True,
                   level: float | None = None, **_kw) -> None:
    dist = DISTS[est.dist]
    form = "accelerated failure-time form" if est.aft else "log relative-hazard form"
    if header:
        header_block(s, f"{dist.title} regression -- {form}", est)
    if not table:
        return
    cluster_note(s, est)
    lev = level or est.level
    rows = list(est.rows)
    if est.aux is not None:
        b, se = est.aux
        z = st.norm.ppf(1 - (1 - lev / 100) / 2)
        rows += [CoefRow("sep"), CoefRow("coef", f"/{dist.aux}", b, se, name=f"{dist.aux}:_cons", noeform=True)]
        extra = dist.aux_rows(b, se, b - z * se, b + z * se)
        if extra:
            rows += [CoefRow("sep")] + extra
    title = {"hr": "Haz. Ratio", "tr": "Tm. Ratio"}.get(est.eform_kind, "Coef.")   # VERIFY: Time Ratio label
    coef_table(s, "_t", rows, stat="z", df=None, level=lev, vcetype=est.vcetype,
               coef_title=title, eform=bool(est.eform_kind))


def _predict_streg(s: "Session", est: Estimates, stat: str, mask: np.ndarray):
    from .postest import xb_all
    xb, ok = xb_all(s, est)
    ds = s.data
    dist, aft = DISTS[est.dist], est.aft
    th = est.aux[0] if est.aux else 0.0
    t = ds.get("_t").data.astype(float)
    t0 = ds.get("_t0").data.astype(float)
    okt = ok & (t < SYS)
    if stat == "xb":
        return np.where(ok, xb, SYS), ""
    if stat in ("", "median"):
        note = "" if stat else "(option median time assumed; predicted median time)"   # VERIFY
        v = dist.median(xb, th, aft)
        return np.where(ok & np.isfinite(v), v, SYS), note
    if stat == "hr":
        if aft:
            raise StataError(198, "option hr not allowed after accelerated failure-time models")   # VERIFY
        return np.where(ok, np.exp(xb), SYS), ""
    if stat in ("surv", "csurv", "csnell", "mgale"):
        S_t = dist.surv(t, xb, th, aft)
        if stat == "csurv":
            return np.where(okt, S_t, SYS), ""
        S_t0 = dist.surv(t0, xb, th, aft)
        if stat == "surv":
            return np.where(okt, S_t / S_t0, SYS), ""
        cs = -np.log(S_t) + np.log(S_t0)
        if stat == "csnell":
            return np.where(okt, cs, SYS), ""
        return np.where(okt, ds.get("_d").data - cs, SYS), ""
    if stat == "hazard":
        v = dist.hazard(t, xb, th, aft)
        return np.where(okt & np.isfinite(v), v, SYS), ""
    if stat in ("mean", "time", "lntime"):
        if stat == "lntime" and isinstance(dist, _LocationScale):
            return np.where(ok, xb, SYS), ""
        raise StataError(198, f"option {stat} is not available in OpenDTA for this distribution")   # OpenDTA
    raise StataError(198, f"option {stat} not allowed")
