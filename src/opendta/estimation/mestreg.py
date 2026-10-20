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
        sc = math.exp(a[0])
        k = 1.0 / sc
        e = np.exp(np.clip(-eta, -700, 700))
        x = t * e / sc
        lg = special.gammaln(k)
        # failures: ln of the gamma density of t
        lf = (k - 1) * lt - x - k * (a[0] + eta) - lg
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
    fixed, levels_spec, glob = split_levels(args)
    if len(levels_spec) != 1:
        raise StataError(198, "OpenDTA's me commands support one level of random effects")   # OpenDTA
    c = find_top(fixed, ",")
    body, fe_opts = (fixed, "") if c == -1 else (fixed[:c], fixed[c + 1:])
    opts_text = " ".join(x for x in (fe_opts.strip(), glob.strip()) if x)
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
    if o.get("nohr") and time_metric:
        raise StataError(198, "option nohr not allowed")   # VERIFY
    if o.get("tratio") and not time_metric:
        raise StataError(198, "option tratio not allowed")   # VERIFY
    stx.notes_and_header(s)
    t, t0, d = stx.t, stx.t0, stx.d
    if np.any((t <= 0) & (d > 0)):
        raise StataError(498, "failures at time 0")   # VERIFY
    smp = stx.smp
    n = stx.n
    X = np.column_stack([stx.X, np.ones(n)])
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
        aux_name, aux_lab = dist.aux, f"/{dist.aux}"
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
    p = X.shape[1]
    b0 = np.zeros(p + naux)
    b0[p - 1] = start[0]
    if naux:
        b0[p] = start[1]
    dist_ = DISTS.get(name)

    def aux_rows(aux, Va, z):
        b = float(aux[0])
        se = float(np.sqrt(max(Va[0, 0], 0)))
        rows = [CoefRow("coef", aux_lab, b, se, name="_cons", noeform=True)]
        extra_rows = dist_.aux_rows(b, se, b - z * se, b + z * se) if dist_ is not None else []
        if extra_rows:
            rows += [CoefRow("sep")] + extra_rows
        return rows

    def edit(macros):
        keep = [(k, v) for k, v in macros if k not in ("model", "link", "family", "dispersion", "cmd")]
        keep.insert(0, ("cmd", "mestreg"))
        keep += [("distribution", name), ("frm2", "time" if time_metric else "hazard")]
        return keep
    shown = _NAMES[name]
    sp = {"cmd": "mestreg", "fam": name, "link": "", "title": f"Mixed-effects {shown} regression",
          "model": "", "vs": f"{shown} model", "disp": "", "nbin": ""}
    eform = "Tm. Ratio" if (time_metric and o.get("tratio")) else \
        ("Haz. Ratio" if (not time_metric and not o.get("nohr")) else "")
    F = Fit(smp=smp, o=o, sp=sp, X=X, off=np.zeros(n), lv=lv, fam_f=fam_f, naux=naux, b0=b0,
            npts=npts, method=method, family=name, link="", aux_eq=[aux_lab] if naux else [],
            aux_in_table=bool(naux), aux_rows=aux_rows if naux else None, eform=eform,
            predict=_predict_mestreg, macro_edit=edit)
    _fit_me(s, args, "mestreg", F)


def _predict_mestreg(s: "Session", est: Estimates, stat: str, mask: np.ndarray):
    if stat in ("xb", "eta", "fixedonly"):
        eta, ok = _posterior_eta(s, est, fixedonly=stat != "eta")
        return np.where(ok, eta, SYS), ""
    raise StataError(198, f"option {stat or 'xb'} is not available yet in OpenDTA after mestreg")   # OpenDTA
