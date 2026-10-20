"""Multilevel parametric survival models: mestreg ([ME] mestreg).

The conditional distribution of the failure time given the random effects is one of the streg
distributions (exponential, Weibull, lognormal and loglogistic, from streg.py) or the gamma
distribution; the linear predictor is x b + z u. The exponential and Weibull models are fitted
in the hazard metric (log relative-hazard form) unless `time` asks for the accelerated
failure-time form; the others are always in the time metric. The likelihood of each group is
integrated by the adaptive quadrature of meglm.py. Each observation contributes the density
(failures) or the survivor function (censored) at its time, divided by the survivor function at
its entry time.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..lang.syntax import find_top
from ..lazy import lazy
from . import fvars
from .meglm import Fit, _fit_me, _posterior_eta
from .mixed import build_levels, split_levels
from .results import CoefRow, Estimates
from .stfit import STSample
from .streg import DISTS, Dist

special = lazy("scipy.special")
st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS

_OPTS = {"distribution": 1, "time": 4, "nohr": 4, "tratio": 2, "noshow": 6, "nolog": 4, "iterate": 4,
         "tolerance": 3, "ltolerance": 4, "nrtolerance": 4, "noheader": 6, "nogroup": 5,
         "nolrtest": 5, "intpoints": 4, "intmethod": 7, "covariance": 3, "startgrid": 6,
         "nofetable": 5, "noretable": 5, "level": 1, "log": 3}

# the title of each distribution in the header and in the LR test line
_NAMES = {"exponential": "exponential", "weibull": "Weibull", "lognormal": "lognormal",
          "loglogistic": "loglogistic", "gamma": "gamma"}


def _gamma_family(t, t0, d):
    """Gamma distribution of the failure time with mean exp(eta) and scale s = exp(a): shape
    1/s; density of the failures, survivor function of the censored, divided by the survivor
    function at the entry time. With x = t exp(-eta) / s, S = Q(1/s, x) and d ln S / d eta =
    g / Q, g = x^(1/s) exp(-x) / Gamma(1/s)."""
    has0 = t0 > 0
    fail = d > 0
    lt = np.log(np.where(t > 0, t, 1.0))

    def f(eta, a):
        lsc = 2.0 * a[0]            # /logs is the log of the square root of the scale (compat 1902)
        sc = math.exp(lsc)
        k = 1.0 / sc
        e = np.exp(np.clip(-eta, -700, 700))
        x = t * e / sc
        lg = special.gammaln(k)
        # failures: ln of the gamma density of t
        lf = (k - 1) * lt - x - k * (lsc + eta) - lg
        sf1 = x - k                                          # d ln f / d eta = x - k
        sf2 = -x
        # censored: ln Q(k, x)
        Q = np.maximum(special.gammaincc(k, x), 1e-300)
        lx = np.log(np.maximum(x, 1e-300))
        lgx = k * lx - x - lg
        g = np.exp(lgx)
        sc1 = g / Q
        sc2 = -sc1 * (k - x) - sc1 * sc1
        l = np.where(fail, lf, np.log(Q))
        s1 = np.where(fail, sf1, sc1)
        s2 = np.where(fail, sf2, sc2)
        if np.any(has0):
            x0 = np.where(has0, t0 * e / sc, 0.0)
            Q0 = np.maximum(special.gammaincc(k, x0), 1e-300)
            g0 = np.exp(k * np.log(np.maximum(x0, 1e-300)) - x0 - lg)
            c1 = np.where(has0, g0 / Q0, 0.0)
            c2 = np.where(has0, -c1 * (k - x0) - c1 * c1, 0.0)
            l = l - np.where(has0, np.log(Q0), 0.0)
            s1 = s1 - c1
            s2 = s2 - c2
        return l, s1, s2
    return f, 1


def _dist(text: str):
    t = text.strip().lower()
    for full, k in (("exponential", 1), ("weibull", 1), ("lognormal", 4), ("lnormal", 2),
                    ("loglogistic", 4), ("llogistic", 2), ("gamma", 2)):
        if len(t) >= k and full.startswith(t):
            return {"lnormal": "lognormal", "llogistic": "loglogistic"}.get(full, full)
    raise StataError(198, f"distribution({text}) invalid")   # VERIFY


def cmd_mestreg(s: "Session", args: str) -> None:
    from .meglm import _p_chibar  # noqa: F401  (kept next to the other me helpers)
    from .postest import replay
    from .regress import cmdline  # noqa: F401
    if not args.strip() or args.strip().startswith(","):
        replay(s, "mestreg", args)
        return
    fixed, levels_spec, glob = split_levels(args) if "||" in args else (args, [], "")
    if len(levels_spec) > 1:
        raise StataError(198, "OpenDTA's me commands support one level of random effects")   # OpenDTA
    c = find_top(fixed, ",")
    body, fe_opts = (fixed, "") if c == -1 else (fixed[:c], fixed[c + 1:])
    opts_text = " ".join(x for x in (fe_opts.strip(), glob.strip()) if x)
    if not levels_spec:         # without || the model has no random effects (compat 1902)
        opts_text, fe_opts = fe_opts.strip(), fe_opts
    extra = " ".join([lv for lv, _r, _o in levels_spec] + [r for _l, r, _o in levels_spec if r])
    stx = STSample(s, body.strip() + (f", {opts_text}" if opts_text else ""), _OPTS, extra_vars=extra)
    o = stx.o
    if not isinstance(o.get("distribution"), str):
        raise StataError(198, "option distribution() required")   # VERIFY
    name = _dist(o["distribution"])
    if stx.wtype:
        raise StataError(198, "weights are not available yet with mestreg in OpenDTA")   # OpenDTA
    ph = name in ("exponential", "weibull")
    if o.get("time") and not ph:
        raise StataError(198, "option time not allowed")   # VERIFY
    time_metric = (not ph) or bool(o.get("time"))
    late = "option nohr not allowed" if (o.get("nohr") and time_metric) else ""   # after the fit (compat 1902)
    if o.get("tratio") and not time_metric:
        raise StataError(198, "option tratio not allowed")   # VERIFY
    stx.notes_and_header(s)
    t, t0, d = stx.t, stx.t0, stx.d
    if np.any((t <= 0) & (d > 0)):
        raise StataError(498, "failures at time 0")   # VERIFY
    smp = stx.smp
    n = stx.n
    X = np.column_stack([stx.X, np.ones(n)])
    lv = None
    if levels_spec:
        lvspec = levels_spec[0]
        if lvspec[2].get("covariance") is None and o.get("covariance"):
            lvspec[2]["covariance"] = o["covariance"]
        lv = build_levels(s, smp.mask, [lvspec], n)[0]
        if lv.cov == "identity" and lv.q > 1 and not lvspec[2].get("covariance"):
            lv.cov = "independent"
    npts = int(o["intpoints"]) if o.get("intpoints") else 7
    method = {"mv": "mvaghermite", "mc": "mcaghermite", "gh": "ghermite"}.get(
        str(o.get("intmethod") or "mvaghermite").strip().lower(),
        str(o.get("intmethod") or "mvaghermite").strip().lower())
    # conditional distribution of the failure time
    t1, t01, d1 = t[:, None], t0[:, None], d[:, None]
    rate = float(np.sum(d) / np.sum(t - t0))
    if name == "gamma":
        fam_f, naux = _gamma_family(t1, t01, d1)
        start, aux_name, aux_lab = [-math.log(rate), 0.0], "logs", "/logs"
    elif name in ("lognormal", "loglogistic"):
        dist = DISTS[name]
        lf = dist.loglik(t1, t01, d1, True)

        def fam_f(u, a):
            r = lf(u, a[0])
            return r[0], r[1], r[3]
        naux = 1
        start = dist.start(rate, True)
        aux_name, aux_lab = "logs", "/logs"                  # compat 1902
    elif name == "exponential":
        lf = DISTS[name].loglik(t1, t01, d1, time_metric)

        def fam_f(u, a):
            r = lf(u, 0.0)
            return r[0], r[1], r[3]
        naux = 0
        start = [-math.log(rate) if time_metric else math.log(rate)]
        aux_name = aux_lab = ""
    else:                                          # weibull
        lf = DISTS[name].loglik(t1, t01, d1, False)
        if time_metric:
            def fam_f(u, a):
                p = math.exp(a[0])
                r = lf(-p * u, a[0])
                return r[0], -p * r[1], p * p * r[3]
        else:
            def fam_f(u, a):
                r = lf(u, a[0])
                return r[0], r[1], r[3]
        naux = 1
        start = [-math.log(rate) if time_metric else math.log(rate), 0.0]
        aux_name, aux_lab = "ln_p", "/ln_p"
    if name == "weibull":
        # like streg, the likelihood includes the sum of d ln t; the exponential one does not (compat 1902)
        lnt = np.where(d1 > 0, np.log(np.where(t1 > 0, t1, 1.0)), 0.0)
        fam_raw = fam_f

        def fam_f(u, a):
            r = fam_raw(u, a)
            return (r[0] + lnt,) + tuple(r[1:])
    p = X.shape[1]
    b0 = np.zeros(p + naux)
    b0[p - 1] = start[0]
    if naux:
        b0[p] = start[1]
    dist_ = DISTS.get(name)

    if lv is None:
        _plain_mestreg(s, smp, X, fam_f, naux, b0, o, name, shown_=_NAMES[name], time_metric=time_metric,
                       eform_=("Time Ratio" if (time_metric and o.get("tratio")) else
                               ("Haz. Ratio" if (not time_metric and not o.get("nohr")) else "")),
                       aux_lab=aux_lab, late=late)
        return

    def aux_rows(aux, Va, z):
        b = float(aux[0])
        se = float(np.sqrt(max(Va[0, 0], 0)))
        return [CoefRow("coef", aux_lab, b, se, name="_cons", noeform=True)]   # compat 1902: no p, 1/p, sigma rows

    def edit(macros):
        # the order of e() after mestreg (compat 1902); depvar is _t and the covariates follow
        d = dict(macros)
        order = ["cmd", "estat_cmd", "predict", "cmdline", "cmd2", "frm2", "model", "title", "distribution",
                 "n_quad", "ivars", "chi2type", "datasignaturevars", "datasignature", "intmethod", "depvar",
                 "marginsnotok", "covariates", "opt", "vce", "user", "ml_method", "technique", "which",
                 "properties"]
        d.update({"cmd": "gsem", "cmd2": "mestreg", "frm2": "time" if time_metric else "hazard",
                  "model": _NAMES[name], "distribution": name, "depvar": "_t",
                  "cmdline": "mestreg _t " + args.strip(), "estat_cmd": "mestreg_estat",
                  "predict": "mestreg_p", "datasignaturevars": " ".join(
                      ["_d", "_t", "_t0"] + [lv.var] + [c for c in d.get("covariates", "").split()])})
        from .results import datasignature
        sv = d["datasignaturevars"].split()
        d["datasignature"] = datasignature(s, sv, smp.mask)
        return [(k, d[k]) for k in order if k in d]
    shown = _NAMES[name]
    aft = " -- AFT" if (ph and time_metric) else ""
    sp = {"cmd": "mestreg", "fam": name, "link": "", "title": f"Mixed-effects {shown} regression{aft}",
          "model": "", "vs": f"{shown} model", "disp": "", "nbin": ""}
    eform = "Time Ratio" if (time_metric and o.get("tratio")) else \
        ("Haz. Ratio" if (not time_metric and not o.get("nohr")) else "")
    F = Fit(smp=smp, o=o, sp=sp, X=X, off=np.zeros(n), lv=lv, fam_f=fam_f, naux=naux, b0=b0,
            npts=npts, method=method, family=name, link="", aux_eq=[aux_lab] if naux else [],
            aux_in_table=bool(naux), aux_rows=aux_rows if naux else None, eform=eform,
            predict=_predict_mestreg, macro_edit=edit, late_error=late,
            predict_opts={"xb": 2, "eta": 3, "reffects": 4, "hazard": 3, "median": 3, "mean": 2})
    _fit_me(s, args, "mestreg", F)


def _predict_mestreg(s: "Session", est: Estimates, stat: str, mask: np.ndarray):
    if stat in ("xb", "eta"):
        eta, ok = _posterior_eta(s, est, fixedonly=stat != "eta")
        return np.where(ok, eta, SYS), ""
    sp = est.extra["spec"]
    if stat == "hazard" and sp["fam"] in ("exponential", "weibull"):
        # hazard at the observed _t given the fixed effects and the posterior means of the effects
        eta, ok = _posterior_eta(s, est, fixedonly=False)
        tt = s.data.get("_t").data.astype(float)
        aft = est.extra["spec"]["title"].endswith("-- AFT")
        p_ = math.exp(float(est.glmm["aux"][0])) if sp["fam"] == "weibull" else 1.0
        e_ = np.exp(-p_ * eta) if aft else np.exp(eta)
        h = p_ * np.power(np.maximum(tt, 1e-300), p_ - 1) * e_
        return np.where(ok & (tt < SYS), h, SYS), "(predictions based on fixed effects and posterior means of random effects)"
    raise StataError(198, f"option {stat or 'xb'} is not available yet in OpenDTA after mestreg")   # OpenDTA


def _plain_mestreg(s, smp, X, fam_f, naux, b0, o, name, *, shown_, time_metric, eform_, aux_lab, late):
    """mestreg without random effects: the same model fitted without the || equation, shown as
    Stata does (title without "Mixed-effects", no group lines, no LR test; compat 1902)."""
    from ..core.stats import chi2_tail  # noqa: F401
    from .meglm import _fit_fixed, _g
    from .mixed import _right
    from .ml import options as ml_options
    from .ml import print_log
    from .results import coef_table, rows_from_columns
    out = s.output
    n = smp.n
    F = Fit(smp=smp, o=o, sp={}, X=X, off=np.zeros(n), lv=None, fam_f=fam_f, naux=naux, b0=b0, npts=0,
            method="", family=name, link="", aux_eq=[])
    r = _fit_fixed(F, ml_options(o))
    if not o.get("nolog"):
        out.write("\n", "text")
        print_log(out, r, "log likelihood")
    if late:
        raise StataError(198, late)
    p = X.shape[1]
    b = r.b
    V = np.linalg.inv(-r.H)
    est_cols = [c for c in smp.cols if c.values is not None]
    k = len(est_cols)
    wchi2 = float(b[:k] @ np.linalg.solve(V[:k, :k], b[:k])) if k else SYS
    wp = float(chi2_tail(wchi2, k)) if k else SYS
    names = [c.name for c in est_cols] + ["_cons"]
    bmap = {nm: b[i] for i, nm in enumerate(names)}
    semap = {nm: float(np.sqrt(max(V[i, i], 0))) for i, nm in enumerate(names)}
    rows = rows_from_columns(smp.cols, bmap, semap, s.data, constant=True)
    if naux:
        rows.append(CoefRow("sep"))
        rows.append(CoefRow("coef", aux_lab, float(b[p]), float(np.sqrt(max(V[p, p], 0))), name="_cons",
                            noeform=True))
    title = f"{shown_} regression" + (" -- AFT" if (name in ("exponential", "weibull") and time_metric) else "")
    out.write("\n", "text")
    _right(out, title, "Number of obs", f"{int(n):,}")
    out.write("\n", "text")
    _right(out, "", f"Wald chi2({k})", _g(wchi2, "%10.2f"))
    _right(out, f"Log likelihood = {_g(r.ll, '%10.0g'):>10}", "Prob > chi2", _g(wp, "%10.4f"))
    coef_table(s, "_t", rows, stat="z", df=None, level=smp.level, coef_title=eform_ or "Coef.",
               eform=bool(eform_))
