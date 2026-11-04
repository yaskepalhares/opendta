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
_ETMODELS = ("linear", "probit", "fractional", "exponential")


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
        eq = ["POmeans", "POmeans"]
    else:
        L = np.array([[-1.0, 1.0], [1.0, 0.0]])
        b = L @ mu
        names = [f"r{lv1}vs{lv0}.{tvar}", f"{lv0}.{tvar}"]
        eq = ["ATE" if stat == "ate" else "ATET", "POmean"]
    V = L @ Vm @ L.T
    se = np.sqrt(np.diag(V))
    if stat == "pomeans":
        rows = [CoefRow("eq", "POmeans"), CoefRow("header", tvar),
                CoefRow("coef", f"{lv0} ", b[0], se[0], name=names[0]),
                CoefRow("coef", f"{lv1} ", b[1], se[1], name=names[1])]
    else:
        rows = [CoefRow("eq", eq[0]), CoefRow("header", tvar),
                CoefRow("coef", f"({lv1} vs {lv0}) ", b[0], se[0], name=names[0]),
                CoefRow("sep"), CoefRow("eq", "POmean"), CoefRow("header", tvar),
                CoefRow("coef", f"{lv0} ", b[1], se[1], name=names[1])]
    return b, V, names, eq, rows


def _effect_rows_multi(tvar, labels, ci, stat, mu, Vm):
    """_effect_rows for L treatment levels: the effects of each level against the control ci (and
    its POmean), or the L POmeans."""
    L = len(mu)
    if stat == "pomeans":
        M_ = np.eye(L)
        names = [f"{lb}.{tvar}" for lb in labels]
        eq = ["POmeans"] * L
    else:
        others = [l for l in range(L) if l != ci]
        M_ = np.zeros((L, L))
        for r, l in enumerate(others):
            M_[r, l], M_[r, ci] = 1.0, -1.0
        M_[L - 1, ci] = 1.0
        names = [f"r{labels[l]}vs{labels[ci]}.{tvar}" for l in others] + [f"{labels[ci]}.{tvar}"]
        eq = ["ATE" if stat == "ate" else "ATET"] * (L - 1) + ["POmean"]
    b = M_ @ mu
    V = M_ @ Vm @ M_.T
    se = np.sqrt(np.diag(V))
    if stat == "pomeans":
        rows = [CoefRow("eq", "POmeans"), CoefRow("header", tvar)]
        rows += [CoefRow("coef", f"{labels[l]} ", b[l], se[l], name=names[l]) for l in range(L)]
    else:
        rows = [CoefRow("eq", eq[0]), CoefRow("header", tvar)]
        rows += [CoefRow("coef", f"({labels[l]} vs {labels[ci]}) ", b[r], se[r], name=names[r])
                 for r, l in enumerate(others)]
        rows += [CoefRow("sep"), CoefRow("eq", "POmean"), CoefRow("header", tvar),
                 CoefRow("coef", f"{labels[ci]} ", b[-1], se[-1], name=names[-1])]
    return b, V, names, eq, rows


def _post(s, cmd, sub, full, ovar, tvar, levels, stat, names, b, V, eq, rows, n, level, mask, o, *, title,
          est_title, omodel_txt, tmodel_txt, match=False, extra_hdr=(), dist="", cens="", control_treated=None,
          et_aux=None, nt=None, multi=None):
    if multi is not None:
        names, b, V, eq = multi[0]
    if et_aux is not None:
        # eteffects posts every parameter (compat 1927: e(b) is 1 x 11 with one regressor)
        names, b, V, eq = et_aux
    est = Estimates(cmd, ovar, names, b, V, n, "z", None)
    est.level, est.vcetype = level, "AI Robust" if match else "Robust"
    est.rows = rows
    est.eqnames = eq
    est.extra = dict(N=n)

    def display(s_, e, *, header=True, table=True, level=None, **_kw):
        out = s_.output
        if header:
            out.write("\n", "text")
            w1, w2 = (47, 19) if match else (48, 18)
            out.write(f"{title:<{w1}}{'Number of obs':<{w2}}= {comma(n):>10}\n", "text")
            if match:
                out.write(f"Estimator      : {est_title:<30}{'Matches: requested':<18} = {extra_hdr[0][1]:>10}\n", "text")
                out.write(f"Outcome model  : {'matching':<30}{'min':>18} = {extra_hdr[1][1]:>10}\n", "text")
                out.write(f"{('Distance metric: ' + dist) if dist else ('Treatment model: ' + tmodel_txt):<47}"
                          f"{'max':>18} = {extra_hdr[2][1]:>10}\n", "text")
            else:
                if cmd != "eteffects":
                    out.write(f"Estimator      : {est_title}\n", "text")
                out.write(f"Outcome model  : {omodel_txt}\n", "text")
                out.write(f"Treatment model: {tmodel_txt}\n", "text")
                if cens:
                    out.write(f"Censoring model: {cens}\n", "text")
        if table:
            coef_table(s_, ovar, e.rows, stat="z", df=None, level=level or e.level, vcetype=e.vcetype)
    est.display = display
    cv, tv_ = control_treated if control_treated else (levels[0], levels[1])
    scalars = [("N", n), ("n_levels", len(levels)), ("k", len(names)), ("treated", float(tv_)), ("control", float(cv))]
    macros = [("cmdline", cmdline(cmd, full)), ("cmd", cmd), ("subcmd", sub), ("title", title),
              ("stat", stat), ("tvar", tvar), ("depvar", ovar), ("omodel", omodel_txt),
              ("tmodel", tmodel_txt), ("vce", "robust"), ("vcetype", "Robust"), ("properties", "b V")]
    if multi is not None:
        # Stata's lists and order with multivalued treatments (compat 1925)
        _, counts, tlevels = multi
        scalars = [("N", n)] + [(f"n{lb}", c) for lb, c in reversed(list(zip(tlevels.split(), counts)))] + \
                  [("k_eq", float(len(dict.fromkeys(eq)))), ("k_levels", float(len(levels))),
                   ("treated", float(tv_)), ("control", float(cv)), ("converged", 1.0)]
        macros = [("cmd", cmd), ("cmdline", cmdline(cmd, full)), ("predict", "teffects_p"),
                  ("estat_cmd", "teffects_estat"), ("marginsnotok", "_ALL"), ("vcetype", "Robust"),
                  ("vce", "robust"), ("depvar", ovar), ("tvar", tvar), ("subcmd", sub), ("omodel", omodel_txt)]
        if sub != "ra":
            macros.append(("tmodel", tmodel_txt))   # VERIFY: place of e(tmodel)
        macros += [("stat", stat), ("title", title), ("tlevels", tlevels), ("properties", "b V")]
    if cmd == "eteffects":
        # Stata's lists and order (compat 1927)
        scalars = [("N", n), ("converged", 1.0), ("rank", float(np.linalg.matrix_rank(V))),
                   ("k_eq", float(len(dict.fromkeys(eq)))), ("k_levels", float(len(levels))),
                   ("n1", float(n - nt)), ("n0", nt)]   # Stata's e(n1) is 900 with 1100 treated (compat 1927; VERIFY)
        macros = [("cmdline", cmdline(cmd, full)), ("marginsnotok", "_ALL"), ("cmd", cmd),
                  ("estat_cmd", "eteffects_estat"), ("predict", "eteffects_p"), ("vcetype", "Robust"),
                  ("vce", "robust"), ("omodel", omodel_txt), ("stat", stat), ("tvar", tvar), ("depvar", ovar),
                  ("properties", "b V")]
    post(s, est, scalars, macros, mask)
    display(s, est, header=not o.get("noheader"), table=not o.get("notable"))


def _et_aux(th, Vf, k1, Xo, Z, t, b, V, names, eq, ovar, tvar, omv, tmv):
    """e(b) and e(V) of eteffects with the auxiliary equations: the outcome equations of each arm
    (OME0, OME1), the treatment probit (TME1) and the coefficients of the residual in each arm
    (TEOM0, TEOM1: VERIFY names and order). The probit block of e(V) is its inverse information,
    with no covariance with the rest (the variance takes it as known)."""
    from scipy import stats
    ko = Xo.shape[1]
    xn = omv.split() + ["_cons"]
    zn = tmv.split() + ["_cons"]
    c0, c1, g = th[2:2 + k1], th[2 + k1:2 + 2 * k1], th[2 + 2 * k1:]
    bb = list(b) + list(c0[:ko]) + list(c1[:ko]) + list(g) + [c0[-1], c1[-1]]
    nm = list(names) + xn + xn + zn + [tvar, tvar]
    eqs = list(eq) + ["OME0"] * ko + ["OME1"] * ko + ["TME1"] * len(zn) + ["TEOM0", "TEOM1"]
    m = len(bb)
    Vb = np.zeros((m, m))
    Vb[:2, :2] = V
    # arms: from the sandwich (mu0, mu1, c0, c1); the first two rows are the reported effects
    sel = list(range(2, 2 + 2 * k1))
    pos = list(range(2, 2 + ko)) + [m - 2] + list(range(2 + ko, 2 + 2 * ko)) + [m - 1]
    Vb[np.ix_(pos, pos)] = Vf[np.ix_(sel, sel)]
    e = Z @ g
    F = np.clip(stats.norm.cdf(e), 1e-15, 1 - 1e-15)
    f = stats.norm.pdf(e)
    info = (Z * (f * f / (F * (1 - F)))[:, None]).T @ Z
    zp = list(range(2 + 2 * ko, 2 + 2 * ko + len(zn)))
    Vb[np.ix_(zp, zp)] = np.linalg.inv(info)
    return nm, np.array(bb), Vb, eqs


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
        omodels = _ETMODELS if cmd == "eteffects" else (_OMODELS if not match else ("linear",))
        head, omodel, onoc = _group(groups[0], omodels, oden)
        ovar, _, omv = head.partition(" ")
        omv = omv.strip()
    thead, tmodel, tnoc = _group(groups[1], _TMODELS, "probit" if cmd == "eteffects" else "logit")
    tvar, _, tmv = thead.partition(" ")
    tmv = tmv.strip()
    if cmd == "eteffects" and tmodel != "probit":
        raise StataError(198, "the treatment model of eteffects is probit")   # VERIFY
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
    if len(levels) < 2:
        raise StataError(420, f"treatment variable {tvar} must have 2 levels")   # VERIFY
    if len(levels) > 2:
        if match or cmd != "teffects":
            raise StataError(198, "multivalued treatments are not allowed with this estimator")   # VERIFY
        return _run_multi(s, sub, full, pr, o, tvar, tmv, tnoc, tmodel, omodel, ovar, levels, tv, mask, stat)
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
    et_aux = None
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
                CoefRow("coef", f"({lv1} vs {lv0}) ", tau, se, name=names[0])]
    else:
        if cmd == "eteffects":
            if omodel in ("probit", "fractional") and np.any((y < 0) | (y > 1)):
                raise StataError(459, f"{ovar} must be between 0 and 1")   # VERIFY
            if omodel == "probit" and np.any((y != 0) & (y != 1)):
                raise StataError(2000, f"{ovar} must be 0 or 1")   # VERIFY
            if omodel == "exponential" and np.any(y < 0):
                raise StataError(459, f"{ovar} must be nonnegative")   # VERIFY
            th, Vf, _k = T.eteffects(y, t, Xo, Z, stat, clusters, omodel=omodel)
        else:
            th, Vf, _lay = T.teffects(sub, y, t, Xo, Z, omodel, tmodel, stat, clusters)
        mu, Vm = th[:2], Vf[:2, :2]
        b, V, names, eq, rows = _effect_rows(tvar, levels, stat, mu, Vm, lv0, lv1)
        et_aux = _et_aux(th, Vf, _k, Xo, Z, t, b, V, names, eq, ovar, tvar, omv, tmv) if cmd == "eteffects" else None
        omodel_txt = {"ipw": "weighted mean", "aipw": "linear by ML"}.get(sub, omodel)

        tmodel_txt = tmodel if sub not in ("ra",) else "none"
        dist = ""
    title = "Endogenous treatment-effects estimation" if cmd == "eteffects" else "Treatment-effects estimation"
    _post(s, cmd, sub, full, ovar, tvar, levels, stat, names, b, V, eq, rows, n, level, mask, o, title=title,
          est_title=est_title, omodel_txt=omodel_txt, tmodel_txt=tmodel_txt, match=match, extra_hdr=extra_hdr,
          dist=dist, et_aux=et_aux, nt=float(t.sum()))


def omv_names(pr, ovar):
    """Names of the outcome-model columns (with _cons when there is one)."""
    blk = block_from_pr(pr, ovar)
    return [c.name for c in blk.cols if c.values is not None] + (["_cons"] if blk.constant else [])


def _level_option(o, key, levels, default):
    if not o.get(key):
        return default
    v = float(o[key])
    if v not in levels:
        raise StataError(198, f"{key}({o[key]}) is not a level of the treatment")   # VERIFY
    return int(np.flatnonzero(levels == v)[0])


def _run_multi(s, sub, full, pr, o, tvar, tmv, tnoc, tmodel, omodel, ovar, levels, tv, mask, stat):
    """teffects ra, ipw, ipwra and aipw with more than two treatment levels: one outcome model per
    level and a multinomial logit treatment model ([TE] teffects multivalued)."""
    if sub != "ra" and tmodel != "logit":
        raise StataError(198, f"treatment-model {tmodel} requires 2 levels in treatment variable {tvar}, "
                              f"but {len(levels)} were found")   # compat 1925 (VERIFY: return code)
    L = len(levels)
    tcode = np.searchsorted(levels, tv)
    ci = _level_option(o, "control", levels, 0)                 # VERIFY: lowest level by default
    tl = _level_option(o, "tlevel", levels, 1 if ci == 0 else 0)   # the first level after the control (compat 1925)
    y = pr.y
    n = len(y)
    Xo = block_from_pr(pr, ovar).X if sub != "ipw" else np.ones((n, 1))
    Z = block_from_text(pr, tvar, tmv, constant=not tnoc).X if tmv or not tnoc else np.ones((n, 1))
    clusters = pr.smp.cluster if pr.smp.vce == "cluster" else None
    level = float(o["level"]) if o.get("level") else float(s.settings.get("level", "95"))
    th, Vf = T.teffects_multi(sub, y, tcode, L, Xo, Z, omodel, stat, tl, clusters)
    labels = [_level_label(s, tvar, v) for v in levels]
    b, V, names, eq, rows = _effect_rows_multi(tvar, labels, ci, stat, th[:L], Vf[:L, :L])
    est_title = {"ra": "regression adjustment", "ipw": "inverse-probability weights",
                 "ipwra": "IPW regression adjustment", "aipw": "augmented IPW"}[sub]
    omodel_txt = {"ipw": "weighted mean", "aipw": "linear by ML"}.get(sub, omodel)
    tmodel_txt = "(multinomial) logit" if sub != "ra" else "none"     # VERIFY
    # e(b) and e(V) with every parameter: the effects, the outcome equations OME0, OME1, ... and the
    # multinomial logit TME1, TME2, ... (compat 1925: 1 x 9 for ra with one regressor; VERIFY the
    # names of the auxiliary equations)
    kb = Xo.shape[1] if sub != "ipw" else 0
    xn = (omv_names(pr, ovar) if sub != "ipw" else [])
    zn = [nm for nm in tmv.split()] + ([] if tnoc else ["_cons"])
    aux_b = list(th[L:])
    aux_names = []
    aux_eq = []
    if kb:
        for l in range(L):
            aux_names += xn
            aux_eq += [f"OME{labels[l]}"] * kb
    if sub != "ra":
        for l in range(1, L):
            aux_names += zn
            aux_eq += [f"TME{labels[l]}"] * len(zn)
    m = len(b) + len(aux_b)
    Vb = np.zeros((m, m))
    Vb[:len(b), :len(b)] = V
    Vb[len(b):, len(b):] = Vf[L:, L:]
    full_post = (list(names) + aux_names, np.concatenate([b, aux_b]), Vb, list(eq) + aux_eq)
    counts = [float(np.sum(tcode == l)) for l in range(L)]
    _post(s, "teffects", sub, full, ovar, tvar, levels, stat, names, b, V, eq, rows, n, level, mask, o,
          title="Treatment-effects estimation", est_title=est_title, omodel_txt=omodel_txt,
          tmodel_txt=tmodel_txt, control_treated=(levels[ci], levels[tl]),
          multi=(full_post, counts, " ".join(labels)))


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
            if sub == "ipw":
                raise StataError(198, f"option {nm} not allowed")
            smodel = nm
    if stx.wtype or stx.sd.info.wtype:
        raise StataError(101, "weights not allowed")   # VERIFY
    stat = _estat(o)
    if sub == "aipw" and stat == "atet":
        raise StataError(198, "option atet not allowed")
    if sub == "ipw" and np.any(stx.d == 0):
        raise StataError(459, "all observations must fail for the IPW model without a censoring distribution; "
                              "this is not the case")   # VERIFY return code
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
          title="Survival treatment-effects estimation", est_title=est_title,
          omodel_txt=("Weibull" if smodel == "weibull" else smodel) if sub != "ipw" else "weighted mean",  # VERIFY lognormal label
          tmodel_txt=tmodel if sub != "ra" else "none", cens="none")


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
