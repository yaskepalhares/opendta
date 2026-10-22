"""Treatment-effects commands: teffects (ra, ipw, ipwra, aipw, nnmatch, psmatch) and eteffects
([TE] teffects, [TE] eteffects). The estimators are in core/treatfx.py. Binary treatments only.

Layouts come from the manuals; the lines marked VERIFY wait for the Stata 14 log of case 1908."""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core import treatfx as T
from .mleq import block_from_pr, block_from_text, names_in
from .models2 import _Prep
from .regress import cmdline, comma
from .results import CoefRow, Estimates, coef_table, post

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS

_OMODELS = ("linear", "logit", "probit", "poisson")
_TMODELS = ("logit", "probit")


def _groups(text: str) -> tuple[list[str], str]:
    """The parenthesised groups at the start of `text` and what follows them."""
    out, i, n = [], 0, len(text)
    while True:
        while i < n and text[i].isspace():
            i += 1
        if i >= n or text[i] != "(":
            return out, text[i:]
        depth, j = 0, i
        while j < n:
            depth += text[j] == "("
            depth -= text[j] == ")"
            j += 1
            if depth == 0:
                break
        out.append(text[i + 1:j - 1].strip())
        i = j


def _group(text: str, models: tuple[str, ...], default: str) -> tuple[str, str, bool]:
    """'y x1 x2, logit noconstant' -> ('y x1 x2', model, noconstant)."""
    head, _, tail = text.partition(",")
    model, noc = default, False
    for w in tail.split():
        if "noconstant".startswith(w) and len(w) >= 3:
            noc = True
        elif w in models:
            model = w
        elif w:
            raise StataError(198, f"option {w} not allowed")   # VERIFY: other models are not available
    return head.strip(), model, noc


def _level_label(s, tname, value) -> str:
    return f"{value:g}"


def _estat(o, default="ate") -> str:
    for k in ("ate", "atet", "pomeans"):
        if o.get(k):
            return k
    return default


def _effect_rows(tvar, levels, stat, mu, Vm, lv0, lv1):
    """(b, V, names, equation names, table rows) of the ATE/ATET and POmean estimates."""
    if stat == "pomeans":
        L = np.eye(2)
        b = mu.copy()
        names = [f"0.{tvar}", f"1.{tvar}"] if levels[0] == 0 and levels[1] == 1 else [f"{lv0}.{tvar}", f"{lv1}.{tvar}"]
        eq = ["POmean", "POmean"]
    else:
        L = np.array([[-1.0, 1.0], [1.0, 0.0]])
        b = L @ mu
        names = [f"r{lv1}vs{lv0}.{tvar}", f"{lv0}.{tvar}"]
        eq = ["ATE" if stat == "ate" else "ATET", "POmean"]
    V = L @ Vm @ L.T
    se = np.sqrt(np.diag(V))
    if stat == "pomeans":
        rows = [CoefRow("eq", "POmean"), CoefRow("header", tvar),
                CoefRow("coef", lv0, b[0], se[0], name=names[0]),
                CoefRow("coef", lv1, b[1], se[1], name=names[1])]
    else:
        rows = [CoefRow("eq", eq[0]), CoefRow("header", tvar),
                CoefRow("coef", f"({lv1} vs {lv0})", b[0], se[0], name=names[0], tight=True),
                CoefRow("sep"), CoefRow("eq", "POmean"), CoefRow("header", tvar),
                CoefRow("coef", lv0, b[1], se[1], name=names[1])]
    return b, V, names, eq, rows


def _post(s, cmd, sub, full, ovar, tvar, levels, stat, names, b, V, eq, rows, n, level, mask, o, *, title,
          est_title, omodel_txt, tmodel_txt, match=False, extra_hdr=(), dist="", cens=""):
    est = Estimates(cmd, ovar, names, b, V, n, "z", None)
    est.level, est.vcetype = level, "Robust"
    est.rows = rows
    est.eqnames = eq
    est.extra = dict(N=n)

    def display(s_, e, *, header=True, table=True, level=None, **_kw):
        out = s_.output
        if header:
            out.write("\n", "text")
            out.write(f"{title:<47}{'Number of obs':<19}= {comma(n):>10}\n", "text")
            if match:
                out.write(f"Estimator      : {est_title:<30}{'Matches: requested':<18} = {extra_hdr[0][1]:>10}\n", "text")
                out.write(f"Outcome model  : {'matching':<30}{'min':>18} = {extra_hdr[1][1]:>10}\n", "text")
                out.write(f"{('Distance metric: ' + dist) if dist else ('Treatment model: ' + tmodel_txt):<47}"
                          f"{'max':>18} = {extra_hdr[2][1]:>10}\n", "text")
            else:
                out.write(f"Estimator      : {est_title}\n", "text")
                out.write(f"Outcome model  : {omodel_txt}\n", "text")
                out.write(f"Treatment model: {tmodel_txt}\n", "text")
                if cens:
                    out.write(f"Censoring model: {cens}\n", "text")   # VERIFY
        if table:
            coef_table(s_, ovar, e.rows, stat="z", df=None, level=level or e.level, vcetype=e.vcetype)
    est.display = display
    scalars = [("N", n), ("n_levels", 2), ("k", len(names)), ("treated", float(levels[1])), ("control", float(levels[0]))]
    macros = [("cmdline", cmdline(cmd, full)), ("cmd", cmd), ("subcmd", sub), ("title", title),
              ("stat", stat), ("tvar", tvar), ("depvar", ovar), ("omodel", omodel_txt),
              ("tmodel", tmodel_txt), ("vce", "robust"), ("vcetype", "Robust"), ("properties", "b V")]
    post(s, est, scalars, macros, mask)
    display(s, est, header=not o.get("noheader"), table=not o.get("notable"))


def cmd_teffects(s: "Session", args: str) -> None:
    sub, _, rest = args.strip().partition(" ")
    sub = sub.strip().lower()
    if sub not in ("ra", "ipw", "ipwra", "aipw", "nnmatch", "psmatch"):
        raise StataError(199, f"teffects {sub}: subcommand not found" if sub else "subcommand required")
    _run(s, "teffects", sub, rest.strip(), args)


def cmd_eteffects(s: "Session", args: str) -> None:
    _run(s, "eteffects", "", args.strip(), args)


def _run(s, cmd, sub, rest, full):
    groups, tail = _groups(rest)
    if len(groups) != 2:
        raise StataError(198, "outcome and treatment models required")
    match = sub in ("nnmatch", "psmatch")
    nn_vce = 2
    mm = re.search(r"vce\(\s*robust\s*,\s*nn\(\s*(\d+)\s*\)\s*\)", tail)
    if mm:
        nn_vce = int(mm.group(1))
        tail = tail.replace(mm.group(0), "vce(robust)")
    oden = "linear"
    if sub in ("ipw", "psmatch"):
        ovar, omodel, onoc = groups[0].split()[0] if groups[0].split() else "", "linear", False
        omv = ""
    else:
        head, omodel, onoc = _group(groups[0], _OMODELS if not match else ("linear",), oden)
        ovar, _, omv = head.partition(" ")
        omv = omv.strip()
    thead, tmodel, tnoc = _group(groups[1], _TMODELS, "probit" if cmd == "eteffects" else "logit")
    tvar, _, tmv = thead.partition(" ")
    tmv = tmv.strip()
    if cmd == "eteffects" and (omodel != "linear" or tmodel != "probit"):
        raise StataError(198, "only linear outcome and probit treatment models are available")   # VERIFY
    spec = {"ate": 3, "atet": 4, "pomeans": 3, "control": 3, "tlevel": 3, "nneighbor": 3,
            "metric": 3, "aequations": 4, "osample": 3, "logit": 4, "probit": 4, "vce": 3,
            "level": 1, "noconstant": 4}
    spec.pop("noconstant")
    args2 = f"{ovar} {omv} {tail}"
    pr = _Prep(s, args2, spec, extra_vars=names_in(tvar, tmv), constant=not onoc and not match and sub != "ipw")
    pr.write_notes(s.output)
    o = pr.o
    stat = _estat(o)
    if match and stat == "pomeans":
        raise StataError(198, "option pomeans not allowed")
    if sub == "aipw" and stat == "atet":
        raise StataError(198, "option atet not allowed")
    mask = pr.smp.mask
    tv = s.data.get(tvar).data[mask].astype(float)
    levels = np.unique(tv)
    if len(levels) != 2:
        raise StataError(420, f"treatment variable {tvar} must have 2 levels") if len(levels) < 2 else \
            StataError(198, "multivalued treatments are not available")   # VERIFY
    t = (tv == levels[1]).astype(float)
    y = pr.y
    n = len(y)
    Xo = block_from_pr(pr, ovar).X if sub != "ipw" else np.ones((n, 1))
    if match and sub == "nnmatch":
        Xo = pr.X
    Z = block_from_text(pr, tvar, tmv, constant=not tnoc).X if tmv or not tnoc else np.ones((n, 1))
    clusters = pr.smp.cluster if pr.smp.vce == "cluster" else None
    level = float(o["level"]) if o.get("level") else float(s.settings.get("level", "95"))
    lv0, lv1 = _level_label(s, tvar, levels[0]), _level_label(s, tvar, levels[1])

    est_title = {"ra": "regression adjustment", "ipw": "inverse-probability weights",
                 "ipwra": "IPW regression adjustment", "aipw": "augmented IPW",
                 "nnmatch": "nearest-neighbor matching", "psmatch": "propensity-score matching",
                 "": "extended regression"}[sub]
    extra_hdr: list[str] = []
    if match:
        if sub == "psmatch":
            zb = T.glm_fit(tmodel, Z, t)
            X = T._mean(tmodel, Z @ zb)[0][:, None]
            metric = "euclidean"
        else:
            X = Xo
            m_ = str(o.get("metric", "maha")).lower()
            metric = "maha" if "mahalanobis".startswith(m_) else ("euclidean" if "euclidean".startswith(m_) else "ivariance")
        nn = int(o["nneighbor"]) if o.get("nneighbor") else 1
        tau, se, mu, (mn, mx) = T.nnmatch(y, t, X, stat, nn, metric, nn_vce)
        names = [f"r{lv1}vs{lv0}.{tvar}"]
        eq = ["ATE" if stat == "ate" else "ATET"]
        b = np.array([tau])
        V = np.array([[se * se]])
        mlabel = {"maha": "Mahalanobis", "euclidean": "Euclidean", "ivariance": "Inverse variance"}[metric]
        extra_hdr = [("Matches: requested", str(nn)), ("min", str(mn)), ("max", str(mx))]
        omodel_txt, tmodel_txt = "matching", ("none" if sub == "nnmatch" else tmodel)
        dist = mlabel if sub == "nnmatch" else ""
        rows = [CoefRow("eq", eq[0]), CoefRow("header", tvar),
                CoefRow("coef", f"({lv1} vs {lv0})", tau, se, name=names[0], tight=True)]
    else:
        if cmd == "eteffects":
            th, Vf, _k = T.eteffects(y, t, Xo, Z, stat, clusters)
        else:
            th, Vf, _lay = T.teffects(sub, y, t, Xo, Z, omodel, tmodel, stat, clusters)
        mu, Vm = th[:2], Vf[:2, :2]
        b, V, names, eq, rows = _effect_rows(tvar, levels, stat, mu, Vm, lv0, lv1)
        omodel_txt = omodel if sub != "ipw" else "none"
        tmodel_txt = tmodel if sub not in ("ra",) else "none"
        dist = ""
    title = "Endogenous treatment-effects estimation" if cmd == "eteffects" else "Treatment-effects estimation"
    _post(s, cmd, sub, full, ovar, tvar, levels, stat, names, b, V, eq, rows, n, level, mask, o, title=title,
          est_title=est_title, omodel_txt=omodel_txt, tmodel_txt=tmodel_txt, match=match, extra_hdr=extra_hdr,
          dist=dist)


# ---------------------------------------------------------------------------
# stteffects
# ---------------------------------------------------------------------------

_SMODELS = ("exponential", "weibull", "lognormal")


def _surv_outcome(name, t, t0, d):
    """LikOutcome of a streg distribution fitted to the potential-outcome mean survival time."""
    from scipy.special import gamma as gfun
    from .streg import DISTS
    dist = DISTS[name]
    lf = dist.loglik(t, t0, d, False)
    rate = float(np.sum(d) / np.sum(t - t0))
    start = dist.start(rate, False)
    if name == "exponential":
        return T.LikOutcome(lf, lambda eta, th: np.exp(-eta), 0, start)
    if name == "weibull":
        return T.LikOutcome(lf, lambda eta, th: np.exp(-eta / np.exp(th)) * gfun(1 + 1 / np.exp(th)), 1, start)
    return T.LikOutcome(lf, lambda u, th: np.exp(u + np.exp(2 * th) / 2), 1, start)


def cmd_stteffects(s: "Session", args: str) -> None:
    from .stfit import STSample
    sub, _, rest = args.strip().partition(" ")
    sub = sub.lower()
    if sub not in ("ra", "ipw", "ipwra"):
        raise StataError(199, f"stteffects {sub}: subcommand not found" if sub else "subcommand required")   # VERIFY: wra
    groups, tail = _groups(rest.strip())
    if len(groups) != (1 if sub == "ipw" else 2):
        raise StataError(198, "outcome and treatment models required")
    if sub == "ipw":
        ohead, smodel, onoc = "", "weibull", False
        thead, tmodel, tnoc = _group(groups[0], _TMODELS, "logit")
    else:
        ohead, smodel, onoc = _group(groups[0], _SMODELS, "weibull")
        thead, tmodel, tnoc = _group(groups[1], _TMODELS, "logit")
    tvar, _, tmv = thead.partition(" ")
    tmv = tmv.strip()
    spec = {"ate": 3, "atet": 4, "pomeans": 3, "control": 3, "tlevel": 3, "level": 1, "vce": 3, "noshow": 6,
            "exponential": 3, "weibull": 3, "lognormal": 4, "logit": 4, "probit": 4}
    stx = STSample(s, f"{ohead} {tail}", spec, extra_vars=names_in(tvar, tmv))
    o = stx.o
    for nm in _SMODELS:
        if o.get(nm):
            smodel = nm
    if stx.wtype or stx.sd.info.wtype:
        raise StataError(101, "weights not allowed")   # VERIFY
    stat = _estat(o)
    if sub == "aipw" and stat == "atet":
        raise StataError(198, "option atet not allowed")
    stx.notes_and_header(s)
    mask = stx.mask
    n = stx.n
    tv = s.data.get(tvar).data[mask].astype(float)
    levels = np.unique(tv)
    if len(levels) != 2:
        raise StataError(198, "treatment variable must have 2 levels")   # VERIFY
    t = (tv == levels[1]).astype(float)
    if np.any((stx.t <= 0) & (stx.d > 0)):
        raise StataError(498, "failures at time 0")
    Z = _design(s, stx, tvar, tmv, tnoc)
    Xo = np.column_stack([stx.X, np.ones(n)]) if (sub != "ipw" and not onoc) else np.ones((n, 1))
    om = _surv_outcome(smodel, stx.t, stx.t0, stx.d)
    est = "ipwra" if sub == "ipw" else sub
    clusters, _cv = stx.clusters()
    th, Vf, _lay = T.teffects(est, None, t, Xo, Z, om, tmodel, stat, clusters)
    level = float(o["level"]) if o.get("level") else float(s.settings.get("level", "95"))
    lv0, lv1 = _level_label(s, tvar, levels[0]), _level_label(s, tvar, levels[1])
    b, V, names, eq, rows = _effect_rows(tvar, levels, stat, th[:2], Vf[:2, :2], lv0, lv1)
    est_title = {"ra": "regression adjustment", "ipw": "inverse-probability weights",
                 "ipwra": "IPW regression adjustment"}[sub]
    _post(s, "stteffects", sub, args, "_t", tvar, levels, stat, names, b, V, eq, rows, n, level, mask, o,
          title="Treatment-effects estimation", est_title=est_title, omodel_txt=smodel.capitalize(),
          tmodel_txt=tmodel if sub != "ra" else "none", cens="independent")


def _design(s, stx, tvar, tmv, noc):
    """Treatment-model design (with the constant) on the sample of an STSample."""
    from . import fvars
    ds = s.data
    with fvars.caching():
        terms = fvars.expand_fv(ds, tmv) if tmv.strip() else []
        cols = fvars.build_design(ds, terms, stx.mask, constant=not noc)
    fvars.drop_collinear(cols, not noc, stx.w)
    free = [c.values for c in cols if c.values is not None]
    parts = free + ([np.ones(stx.n)] if not noc else [])
    return np.column_stack(parts)
