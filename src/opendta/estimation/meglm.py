"""Mixed-effects generalized linear models: meglm and its special cases melogit, meprobit,
mecloglog, mepoisson, menbreg, meologit and meoprobit ([ME] meglm and the others).

Two-level models with random intercepts and slopes. This module parses the syntax, reads the
data, runs the maximization and stores e(); the likelihood is in glmm.py and the output and
predict in meglm_out.py. The likelihood of each group is integrated by mean-variance adaptive
Gauss-Hermite quadrature (intmethod(mvaghermite), intpoints(7) by default; with q random
effects, the product rule with 7^q points); intmethod(ghermite) gives the non-adaptive rule and
intmethod(laplace) a single point. The adaptive parameters (posterior means and Cholesky
factors) are updated at each iteration and held fixed within it; the gradient is analytic and
the Hessian is the numerical derivative of the gradient.

Stata's starting values ("Refining starting values" with a grid search)
are not documented, so the iteration logs may differ; the final
estimates should agree to the precision of the quadrature (VERIFY:
compat/do/1002_melogit.do).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.stats import chi2_tail
from ..lang.syntax import find_top
from ..lazy import lazy
from . import fvars
from .glmm import GLMM, MGLMM, basic_family, family_of, ordinal_family
from .meglm_out import COV_TITLE, display_me, display_qr, pr_multi, predict_me, reffects_me
from .mixed import Level, _p_chibar, build_levels, estat_group, estat_icc, estat_recov, split_levels
from .ml import maximize, print_log
from .ml import options as ml_options
from .results import CoefRow, Estimates, g9, post, rows_from_columns

# imported on first use (start-up time; see opendta/lazy.py)
special = lazy("scipy.special")
st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


# ---------------------------------------------------------------------------
# the commands
# ---------------------------------------------------------------------------

# command -> (family, link, title, model, reference model of the LR test)
_CMDS = {
    "melogit": ("bernoulli", "logit", "Mixed-effects logistic regression", "logistic", "logistic model"),
    "meprobit": ("bernoulli", "probit", "Mixed-effects probit regression", "probit", "probit model"),
    "mecloglog": ("bernoulli", "cloglog", "Mixed-effects cloglog regression", "cloglog", "cloglog model"),
    "mepoisson": ("poisson", "log", "Mixed-effects Poisson regression", "Poisson", "Poisson model"),
    "menbreg": ("nbinomial", "log", "Mixed-effects nbinomial regression", "nbinomial", "nbinomial model"),
    "meologit": ("ordinal", "logit", "Mixed-effects ologit regression", "ologit", "ologit model"),
    "meoprobit": ("ordinal", "probit", "Mixed-effects oprobit regression", "oprobit", "oprobit model"),
    "meglm": (None, None, "Mixed-effects GLM", "", ""),
    "meqrlogit": ("bernoulli", "logit", "Mixed-effects logistic regression", "logistic", "logistic model"),
    "meqrpoisson": ("poisson", "log", "Mixed-effects Poisson regression", "Poisson", "Poisson model"),
}
_FAMILIES = ("gaussian", "bernoulli", "binomial", "gamma", "nbinomial", "ordinal", "poisson")
_DEFAULT_LINK = {"gaussian": "identity", "bernoulli": "logit", "binomial": "logit", "gamma": "log",
                 "nbinomial": "log", "ordinal": "logit", "poisson": "log"}
_LINKS_OK = {"gaussian": ("identity", "log"), "bernoulli": ("logit", "probit", "cloglog"),
             "binomial": ("logit", "probit", "cloglog"), "gamma": ("log",), "nbinomial": ("log",),
             "ordinal": ("logit", "probit"), "poisson": ("log",)}
# reference model of the LR test of meglm (VERIFY the names)
_VS_GLM = {("gaussian", "identity"): "linear model", ("gaussian", "log"): "linear model", ("bernoulli", "logit"): "logistic model",
           ("bernoulli", "probit"): "probit model", ("bernoulli", "cloglog"): "cloglog model",
           ("binomial", "logit"): "logistic model", ("binomial", "probit"): "probit model",
           ("binomial", "cloglog"): "cloglog model", ("gamma", "log"): "gamma model",
           ("nbinomial", "log"): "nbinomial model", ("ordinal", "logit"): "ologit model",
           ("ordinal", "probit"): "oprobit model", ("poisson", "log"): "Poisson model"}
_GLOBAL = {"or": 2, "irr": 3, "intpoints": 4, "intmethod": 7, "nolog": 4, "log": 3, "iterate": 4,
           "tolerance": 3, "ltolerance": 4, "nrtolerance": 4, "level": 1, "noheader": 6, "nogroup": 5,
           "nolrtest": 5, "noconstant": 4, "offset": 3, "exposure": 3, "binomial": 3, "covariance": 3,
           "startgrid": 6, "nofetable": 5, "noretable": 5, "estmetric": 4, "family": 3, "link": 4,
           "dispersion": 4, "asis": 4, "laplace": 3, "stddeviations": 3, "variance": 3, "refineopts": 6}


def _word(text: str, choices, what: str) -> str:
    t = text.strip().lower()
    hit = [c for c in choices if c.startswith(t)] if t else []
    if len(hit) != 1 and t in choices:
        hit = [t]
    if len(hit) != 1:
        raise StataError(198, f"{what}({text.strip()}) invalid")   # VERIFY
    return hit[0]


def _resolve(cmdname: str, o: dict) -> dict:
    """Family, link and titles of the command, from its name and the family()/link() options of
    meglm."""
    fam, link, title, model, vs = _CMDS[cmdname]
    disp = ""
    nb = ""
    if cmdname == "meglm":
        words = str(o.get("family") or "gaussian").split()
        fam = _word(words[0], _FAMILIES, "family")
        rest = words[1:]
        if fam == "nbinomial":
            disp = _word(rest[0], ("mean", "constant"), "family") if rest else "mean"
        elif fam == "binomial":
            nb = rest[0] if rest else "1"
        elif rest:
            raise StataError(198, f"family({' '.join(words)}) invalid")   # VERIFY
        link = _word(str(o["link"]), ("identity", "log", "logit", "probit", "cloglog"), "link") \
            if o.get("link") else _DEFAULT_LINK[fam]
        if link not in _LINKS_OK[fam]:
            raise StataError(198, f"link({link}) with family({fam}) is not available yet in OpenDTA")   # OpenDTA
        vs = _VS_GLM[(fam, link)]
        model = vs[:-6] if vs.endswith(" model") else ""
    elif cmdname == "menbreg":
        disp = _word(str(o["dispersion"]), ("mean", "constant"), "dispersion") if o.get("dispersion") else "mean"
    elif o.get("dispersion"):
        raise StataError(198, "option dispersion() not allowed")
    if o.get("binomial") is not None and o.get("binomial") is not False:
        if cmdname not in ("melogit", "meprobit", "mecloglog"):
            raise StataError(198, "option binomial() not allowed")
        fam, nb = "binomial", str(o["binomial"]).strip()
    if fam == "binomial" and not nb:
        nb = "1"
    return {"cmd": cmdname, "fam": fam, "link": link, "title": title, "model": model, "vs": vs,
            "disp": disp, "nbin": nb}


def _fixed_start(sp, y, off, p, constant, n_nb):
    """Starting values (beta, auxiliary) of the fixed-effects model."""
    fam, link = sp["fam"], sp["link"]
    b0 = np.zeros(p)
    ybar = float(y.mean())
    if fam == "gaussian":
        c = ybar - float(off.mean())
        a0 = [math.log(max(float(y.var()), 1e-8))]
    elif fam == "gamma":
        c = math.log(max(ybar, 1e-10)) - float(off.mean())
        a0 = [0.5 * math.log(max(float(y.var()) / max(ybar * ybar, 1e-10), 1e-2))]
    elif fam == "nbinomial":
        c = math.log(max(ybar, 1e-10)) - float(off.mean())
        over = max(float(y.var()) / max(ybar, 1e-10) - 1.0, 0.05)
        a0 = [math.log(over / max(ybar, 1e-10)) if sp["disp"] == "mean" else math.log(over)]
    elif fam == "poisson":
        c = math.log(max(ybar, 1e-10)) - float(off.mean())
        a0 = []
    else:           # bernoulli, binomial
        pr = min(max(float((y / n_nb).mean()), 1e-4), 1 - 1e-4)
        c = {"logit": math.log(pr / (1 - pr)), "probit": float(st.norm.ppf(pr)),
             "cloglog": math.log(-math.log(1 - pr))}[link]
        a0 = []
    if constant:
        b0[-1] = c
    return b0, np.array(a0, dtype=float)


def _num_gh(gfun, b):
    """Hessian by central differences of the gradient."""
    k = b.size
    H = np.zeros((k, k))
    for i in range(k):
        hi = 1e-5 * (abs(b[i]) + 1)
        e = np.zeros(k)
        e[i] = hi
        H[:, i] = (gfun(b + e) - gfun(b - e)) / (2 * hi)
    return (H + H.T) / 2


@dataclass
class Fit:
    """What a command hands to fit_me once it has read its data: design, offset, level, the
    conditional distribution of the outcome with its auxiliary parameters, and the details that
    differ from one command to the next."""
    smp: object
    o: dict
    sp: dict
    X: np.ndarray
    off: np.ndarray
    lv: Level
    fam_f: object
    naux: int
    b0: np.ndarray                    # starting values of the fixed-effects model (with auxiliary)
    npts: int
    method: str
    family: str
    link: str
    aux_eq: list
    off_name: str = ""
    ordered: bool = False
    analytic: bool = False            # analytic Hessian of the fixed-effects model
    aux_in_table: bool = False
    aux_rows: object = None           # (aux, Va, z) -> rows of the coefficient table
    resid: object = None              # latent residual variance of the outcome (float or f(aux))
    nbin: object = None
    nlev: int = 0
    qr: bool = False
    eform: str = ""
    predict: object = None
    predict_multi: dict = field(default_factory=dict)
    macro_edit: object = None
    cmdline: str = ""                 # e(cmdline) when Stata rebuilds it (meglm, mestreg)
    late_error: str = ""              # an error that Stata only raises after the fit
    predict_opts: dict | None = None
    cats: list = field(default_factory=list)      # outcome values of the ordinal models
    lvs: list = field(default_factory=list)       # all the levels, the highest first (default [lv])
    wobs: object = None               # observation weights ([fweight], [pweight], [iweight])
    wgrp: object = None               # group weights of each level (fweight(), pweight()) or None
    level_opts: list = field(default_factory=list)   # the options of each random-effects equation


def cmd_me(s: "Session", args: str, cmdname: str) -> None:
    from .postest import replay
    from .regress import Sample
    if not args.strip() or args.strip().startswith(","):
        replay(s, cmdname, args)
        return
    fixed, levels_spec, glob = split_levels(args)
    import re
    m_ = re.search(r"\b(exposure|offset)\s*\(\s*\w+\s*\)", glob)
    if m_ and not cmdname.startswith("meqr"):      # compat 1902: these belong to the fixed-effects equation (before the ||)
        raise StataError(198, f"option {re.sub(chr(32), '', m_.group(0))} allowed only in the fixed-effects equation")
    if not levels_spec:
        raise StataError(198, "random-effects equation required")   # VERIFY
    qr_ = cmdname.startswith("meqr")
    c = find_top(fixed, ",")
    fixed_body, fe_opts = (fixed, "") if c == -1 else (fixed[:c], fixed[c + 1:])
    opts_text = " ".join(x for x in (fe_opts.strip(), glob.strip()) if x)
    extra = " ".join([lv for lv, _r, _o in levels_spec] + [r for _l, r, _o in levels_spec if r])
    for _l, _r, _o in levels_spec:
        for _wk in ("fweight", "pweight"):
            if isinstance(_o.get(_wk), str):
                extra += " " + _o[_wk].strip()
    try:
        smp = Sample(s, fixed_body + (f", {opts_text}" if opts_text else ""), opts_spec=_GLOBAL,
                     weights=() if qr_ else ("fweight", "iweight", "pweight"), extra_vars=extra)
    except StataError as e_:
        if e_.message == "aweights not allowed":
            raise StataError(e_.rc, "aweight not allowed") from None     # compat 1916
        raise
    o = smp.o
    sp = _resolve(cmdname, o)
    fam, link = sp["fam"], sp["link"]
    qr = cmdname.startswith("meqr")          # the older commands: same model, layout of mixed
    for _opt in ("stddeviations", "variance"):
        if o.get(_opt) and not qr:      # compat 1902: only the meqr commands accept them
            raise StataError(198, f"option {_opt} not allowed")
    if fam == "ordinal":
        smp.constant = False
    out = s.output
    ds = s.data
    y = smp.y.astype(np.float64)
    n = smp.n
    nbin = np.ones(n)
    if fam == "bernoulli":
        y = (y != 0).astype(float)
    elif fam == "binomial":
        t = sp["nbin"]
        try:
            nbin = np.full(n, float(t))
        except ValueError:
            nbin = ds.get(t).data[smp.mask].astype(float)
        if np.any(y < 0) or np.any(y > nbin):
            raise StataError(459, f"{smp.depname} must be between 0 and the number of trials")   # VERIFY
    elif fam == "gamma":
        if np.any(y <= 0):
            raise StataError(459, f"{smp.depname} must be greater than zero")   # VERIFY
    elif fam in ("poisson", "nbinomial"):
        if np.any(y < 0):
            raise StataError(459, f"{smp.depname} must be greater than or equal to zero")   # VERIFY
    yi = None
    if fam == "ordinal":
        cats = sorted(set(y.tolist()))
        if len(cats) < 2:
            raise StataError(2000, "outcome does not vary")   # VERIFY
        if len(cats) > 50:
            raise StataError(149, f"{smp.depname}:  too many outcomes")   # VERIFY
        yi = np.searchsorted(cats, y)
    fvars.drop_collinear(smp.cols, smp.constant)
    est_cols = [col for col in smp.cols if col.values is not None]
    X = np.column_stack([col.values for col in est_cols] + ([np.ones(n)] if smp.constant else [])) \
        if (est_cols or smp.constant) else np.zeros((n, 0))
    off = np.zeros(n)
    off_name = ""
    if o.get("offset") or o.get("exposure"):
        if o.get("exposure") and cmdname in ("melogit", "meprobit", "mecloglog", "meologit", "meoprobit"):
            raise StataError(198, "option exposure() not allowed")
        nm = str(o.get("offset") or o.get("exposure")).strip()
        v = ds.get(nm).data[smp.mask].astype(float)
        off = np.log(v) if o.get("exposure") else v
        off_name = f"ln({nm})" if o.get("exposure") else nm
    lvspec = levels_spec[-1]
    if lvspec[2].get("covariance") is None and o.get("covariance"):
        lvspec[2]["covariance"] = o["covariance"]
    if any(lv_ == "_all" for lv_, _r, _o in levels_spec):
        raise StataError(198, "crossed random effects (_all:) are not available for this command")   # OpenDTA
    levels = build_levels(s, smp.mask, levels_spec, n)
    for lv_, (_l, _r, ro_) in zip(levels, levels_spec):
        if lv_.cov == "identity" and lv_.q > 1 and not ro_.get("covariance"):
            lv_.cov = "independent"
    lv = levels[0]
    wobs, wgrp = level_weights(s, smp, levels, levels_spec)
    npts = int(o["intpoints"]) if o.get("intpoints") else 7
    method = str(o.get("intmethod") or "mvaghermite").strip().lower()
    method = {"mv": "mvaghermite", "mvaghermite": "mvaghermite", "mc": "mcaghermite",
              "mcaghermite": "mcaghermite", "gh": "ghermite", "ghermite": "ghermite",
              "laplace": "laplace"}.get(method, method)
    if qr and o.get("laplace"):
        method = "laplace"
    # the conditional distribution of the outcome
    if fam == "bernoulli" and link in ("logit", "probit"):
        fam_f, naux = basic_family(link, y), 0
    elif fam == "poisson":
        fam_f, naux = basic_family("poisson", y), 0
    elif fam == "ordinal":
        fam_f, naux = ordinal_family(link, yi)
    else:
        fam_f, naux = family_of("binomial" if fam == "bernoulli" else fam, link, y, nbin=nbin,
                                disp=sp["disp"])
    p = X.shape[1]
    analytic = naux == 0 and fam != "binomial" and not (fam == "bernoulli" and link == "cloglog")
    if analytic:
        b0 = _fixed_start(sp, y, off, p, smp.constant, nbin)[0]
    elif fam == "ordinal":
        cum = np.array([np.mean(yi <= j) for j in range(len(cats) - 1)])
        b0 = np.concatenate([np.zeros(p), (special.logit if link == "logit" else st.norm.ppf)(cum)])
    else:
        b0b, a0 = _fixed_start(sp, y, off, p, smp.constant, nbin)
        if fam == "gaussian" and link == "log":
            b0b = np.zeros(p)
            if smp.constant:
                b0b[p - 1] = math.log(max(float(y.mean()), 0.1))
            a0 = np.array([math.log(max(float(np.var(y)), 1e-8))])
        elif fam == "gaussian" and p:
            b0b = np.linalg.lstsq(X, y - off, rcond=None)[0]
            a0 = np.array([math.log(max(float(np.mean((y - off - X @ b0b) ** 2)), 1e-8))])
        b0 = np.concatenate([b0b, a0])
    if fam == "ordinal":
        aux_eq = [f"/cut{j + 1}" for j in range(naux)]
    elif fam == "nbinomial":
        aux_eq = ["/lnalpha" if sp["disp"] == "mean" else "/lndelta"]
    elif fam == "gamma":
        aux_eq = ["/logs"]                                      # compat 1902
    elif fam == "gaussian":
        aux_eq = [f"/var(e.{smp.depname})"]
    else:
        aux_eq = []
    resid = {("bernoulli", "logit"): math.pi ** 2 / 3, ("bernoulli", "probit"): 1.0,
             ("bernoulli", "cloglog"): math.pi ** 2 / 6, ("binomial", "logit"): math.pi ** 2 / 3,
             ("binomial", "probit"): 1.0, ("binomial", "cloglog"): math.pi ** 2 / 6,
             ("ordinal", "logit"): math.pi ** 2 / 3, ("ordinal", "probit"): 1.0}.get((fam, link))
    if fam == "gaussian":
        resid = lambda a: float(math.exp(a[0]))                 # noqa: E731
    eform = "Odds Ratio" if (link == "logit" and o.get("or")) else \
        ("IRR" if (fam in ("poisson", "nbinomial") and o.get("irr")) else "")
    F = Fit(smp=smp, o=o, sp=sp, X=X, off=off, off_name=off_name, lv=lv, lvs=levels, wobs=wobs, wgrp=wgrp,
            level_opts=[ro_ for _l, _r, ro_ in levels_spec],
            fam_f=fam_f, naux=naux,
            ordered=fam == "ordinal", b0=b0, analytic=analytic, aux_eq=aux_eq,
            aux_in_table=fam in ("ordinal", "nbinomial", "gamma"), resid=resid, nbin=nbin,
            nlev=(len(cats) if fam == "ordinal" else 0), qr=qr, eform=eform, npts=npts, method=method,
            family=fam, link=link)
    if fam == "ordinal":
        F.predict_multi = {"": pr_multi, "pr": pr_multi}
        F.cats = [float(c) for c in cats]
    if cmdname == "meglm":      # Stata rebuilds the command line (compat 1902; VERIFY for other forms)
        F.cmdline = (f"meglm {fixed_body.strip()}   ,  "
                     + "".join(f"|| {lvn}: {rev_.strip()}{' ' if rev_.strip() else ''}" for lvn, rev_, _ro in levels_spec)
                     + ", " + "    ".join(_split_opts(glob)))
    fit_me(s, args, cmdname, F)


def level_weights(s: "Session", smp, levels, levels_spec):
    """Observation weights (from [fweight], [pweight] or [iweight]) and the group weights of
    each level (its fweight() or pweight() option, which must be constant within a group)."""
    wobs = None
    if smp.wtype:
        wobs = np.asarray(smp.w, dtype=float)[smp.mask]
    wgrp = []
    kinds = {smp.wtype} if smp.wtype else set()
    for lv, (_l, _r, ro) in zip(levels, levels_spec):
        nm = next(((k, ro[k]) for k in ("fweight", "pweight") if isinstance(ro.get(k), str)), None)
        if nm is None:
            wgrp.append(None)
            continue
        kinds.add(nm[0])
        from ..core.varlist import resolve_name
        v = s.data.get(resolve_name(s.data, nm[1].strip())).data[smp.mask].astype(float)
        G = int(lv.gid.max()) + 1
        first = np.full(G, np.nan)
        first[lv.gid[::-1]] = v[::-1]
        if np.any(v != first[lv.gid]):
            raise StataError(198, f"weights in variable {nm[1].strip()} not constant within groups defined by: {lv.var}")   # compat 1916
        if np.any(first < 0):
            raise StataError(402, "negative weights encountered")   # VERIFY
        wgrp.append(first)
    if len(kinds - {"iweight"}) > 1:
        raise StataError(198, "frequency weights and sampling weights may not be combined")   # VERIFY
    if all(w is None for w in wgrp):
        wgrp = None
    return wobs, wgrp


def _split_opts(text: str) -> list[str]:
    """The options of a command line as separate words (parentheses kept together)."""
    out, depth, cur = [], 0, ""
    for ch in text:
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth -= 1
        if ch.isspace() and depth == 0:
            if cur:
                out.append(cur)
            cur = ""
        else:
            cur += ch
    if cur:
        out.append(cur)
    return out


def _me_vce(F: "Fit"):
    """vce(oim), vce(robust) or vce(cluster clustvar); sampling weights imply vce(robust)."""
    smp = F.smp
    vce = getattr(smp, "vce", "ols")
    clustvar = getattr(smp, "clustvar", "") or ""
    if vce not in ("ols", "oim", "robust", "cluster"):
        raise StataError(198, f"vcetype '{vce}' not allowed")   # VERIFY
    pw = getattr(smp, "wtype", "") == "pweight" or any(
        isinstance(ro.get("pweight"), str) for ro in (F.level_opts or []))
    if vce in ("ols", "oim"):
        vce = "robust" if pw else "oim"
    if vce != "oim" and F.qr:
        raise StataError(198, f"option vce({vce}) not allowed")   # VERIFY
    if clustvar:
        gtop = (F.lvs or [F.lv])[0].gid
        cl = np.asarray(smp.cluster)
        first = np.zeros(int(gtop.max()) + 1, dtype=cl.dtype)
        first[gtop[::-1]] = cl[::-1]
        if np.any(first[gtop] != cl):
            raise StataError(198, f"highest-level groups are not nested within {clustvar}")   # VERIFY
    return vce, clustvar


def _robust_v(model, b, V, F: "Fit", clustvar: str):
    """Sandwich variance with the scores of the highest-level groups summed within the clusters
    (each highest-level group is a cluster under vce(robust)); multiplier G/(G-1) (VERIFY)."""
    Sc = model.scores(b)
    gtop = model.gtop
    if clustvar:
        cl = np.asarray(F.smp.cluster)
        G = Sc.shape[0]
        first = np.zeros(G, dtype=cl.dtype)
        first[gtop[::-1]] = cl[::-1]
        _, cid = np.unique(first, return_inverse=True)
        Sc = np.array([Sc[cid == c].sum(axis=0) for c in range(cid.max() + 1)])
    g = Sc.shape[0]
    M_ = Sc.T @ Sc * (g / (g - 1) if g > 1 else 1.0)
    Vr = V @ M_ @ V
    return (Vr + Vr.T) / 2, g


def _obs_weight(F: "Fit"):
    """The weight of each observation in the model without random effects: its own weight times
    the weights of its groups."""
    w = None if F.wobs is None else np.asarray(F.wobs, dtype=float).copy()
    for lv, wg in zip(F.lvs or [F.lv], F.wgrp or []):
        if wg is not None:
            w = (np.ones(F.X.shape[0]) if w is None else w) * wg[lv.gid]
    return w


def fit_fixed(F: "Fit", kw: dict):
    """The model without random effects: its log likelihood is the reference of the LR test and
    its estimates start the full model."""
    X, off, naux, b0 = F.X, F.off, F.naux, F.b0
    p = X.shape[1]
    wt = _obs_weight(F)
    if wt is None:
        fam_f = F.fam_f
    else:
        def fam_f(eta, a=None):
            r = F.fam_f(eta, a) if a is not None else F.fam_f(eta)
            return tuple(x * wt[:, None] for x in r)
    if F.analytic:
        def fe_fun(b, todo):
            eta = (X @ b + off)[:, None]
            lf, s1, s2 = fam_f(eta)
            ll = float(lf.sum())
            if todo == 0:
                return ll, None, None
            return ll, X.T @ s1[:, 0], (X * s2[:, 0][:, None]).T @ X
    else:
        def fe_grad(b):
            eta = (X @ b[:p] + off)[:, None]
            lf, s1, s2 = fam_f(eta, b[p:])
            g = np.zeros(b.size)
            g[:p] = X.T @ s1[:, 0]
            for m in range(naux):
                h = 1e-5 * (abs(b[p + m]) + 1)
                e = np.zeros(naux)
                e[m] = h
                g[p + m] = float(np.sum(fam_f(eta, b[p:] + e)[0] - fam_f(eta, b[p:] - e)[0])) / (2 * h)
            return g

        def fe_fun(b, todo):
            if F.ordered and np.any(np.diff(b[p:]) <= 0):
                return -np.inf, None, None
            eta = (X @ b[:p] + off)[:, None]
            ll = float(fam_f(eta, b[p:])[0].sum())
            if todo == 0:
                return ll, None, None
            return ll, fe_grad(b), _num_gh(fe_grad, b)
    return maximize(fe_fun, b0, **kw)


def fit_me(s: "Session", args: str, cmdname: str, F: "Fit") -> None:
    """Fits the mixed-effects model with the conditional distribution F.fam_f and posts and
    shows the results (the part that meglm and mestreg have in common)."""
    from .regress import cmdline
    smp, o, sp, X, off, lv = F.smp, F.o, F.sp, F.X, F.off, F.lv
    fam, link, qr, naux = F.family, F.link, F.qr, F.naux
    off_name, aux_eq, npts, method = F.off_name, F.aux_eq, F.npts, F.method
    out = s.output
    ds = s.data
    n = smp.n
    lvs = F.lvs or [lv]
    if len(lvs) == 1:
        model = GLMM(fam, np.zeros(n), X, off, lv, npts, method, fam=F.fam_f, naux=naux, ordered=F.ordered,
                     wobs=F.wobs, wgrp=(F.wgrp[0] if F.wgrp else None))
    else:
        if method == "laplace":
            raise StataError(198, "intmethod(laplace) with more than one level is not available yet in OpenDTA")   # OpenDTA
        model = MGLMM(X, off, lvs, npts, method, F.fam_f, naux=naux, ordered=F.ordered, wobs=F.wobs,
                      wgrp=F.wgrp)
    vce, clustvar = _me_vce(F)
    robust = vce in ("robust", "cluster")
    kw = ml_options(o)
    p = X.shape[1]
    est_cols = [col for col in smp.cols if col.values is not None]
    rfe = fit_fixed(F, kw)
    # starting values for the variance components: best of a small grid (VERIFY)
    best = None
    for ls in (math.log(0.5), 0.0, math.log(0.25), math.log(1.5)):
        ths = []
        for lv_ in lvs:
            th_ = np.zeros(lv_.npar())
            th_[:min(lv_.q, th_.size) if lv_.cov != "exchangeable" else 1] = ls
            ths.append(th_)
        model.reset()
        bb = np.concatenate([rfe.b] + ths)
        model.adapt(bb)
        ll = model.loglik(bb)
        if best is None or ll > best[0]:
            best = (ll, bb, model.mu.copy(), model.C.copy())
    model.mu, model.C = best[2], best[3]
    lab = "log pseudolikelihood" if robust else "log likelihood"
    if not o.get("nolog"):
        if qr:
            # no fixed-effects fit is shown; the first stage refines the starting values (VERIFY)
            import re as _re
            if not _re.search(r"iterate\(\s*0\s*\)", str(o.get("refineopts") or "")):   # compat 1902
                out.write("\nRefining starting values:\n\n", "text")
                out.write(f"{'Iteration 0:':<14} {lab} = {g9(best[0], '%10.0g'):>10}\n", "text")
            out.write("\nPerforming gradient-based optimization:\n\n", "text")
        else:
            # the fixed-effects fit and the grid show "log likelihood" even with vce(robust) or
            # sampling weights (compat 1916)
            out.write("\nFitting fixed-effects model:\n\n", "text")
            print_log(out, rfe)
            out.write("\nRefining starting values:\n\n", "text")
            out.write(f"{'Grid node 0:':<14} log likelihood = {g9(best[0], '%10.0g'):>10}\n", "text")
            out.write("\nFitting full model:\n\n", "text")
    r = maximize(model.fun, best[1], **kw)
    if not o.get("nolog"):
        print_log(out, r, lab)
    b = r.b
    try:
        V = np.linalg.inv(-r.H)
        V = (V + V.T) / 2
    except np.linalg.LinAlgError:
        V = np.full((b.size, b.size), np.nan)
    n_clust = 0
    V_model = V
    if robust:
        V, n_clust = _robust_v(model, b, V, F, clustvar)
        # the note and e(N_clust) count the highest-level groups also with vce(cluster) (compat
        # 1916: "30 clusters in cl" for 6 clusters of 5 groups; VERIFY)
        n_clust = int(lvs[0].gid.max()) + 1
    beta, th = model.split(b)
    aux = model.aux(b)
    Vb = V[:p, :p]
    Va = V[p:p + naux, p:p + naux]
    Vth = V[p + naux:, p + naux:]
    slopes = list(range(len(est_cols)))
    try:
        wchi2 = float(beta[slopes] @ np.linalg.solve(Vb[np.ix_(slopes, slopes)], beta[slopes])) if slopes else SYS
    except np.linalg.LinAlgError:
        wchi2 = SYS
    wp = float(chi2_tail(wchi2, len(slopes))) if slopes and wchi2 < SYS else SYS
    chi2_c = max(2 * (r.ll - rfe.ll), 0.0)
    names = [col.name if not col.omitted else fvars.omitted_name(col.name) for col in smp.cols]
    if smp.constant:
        names.append("_cons")
    kk = len(names)
    pos = [i for i, col in enumerate(smp.cols) if col.values is not None] + ([kk - 1] if smp.constant else [])
    bf = np.zeros(kk)
    Vf = np.zeros((kk, kk))
    bf[pos] = beta
    Vf[np.ix_(pos, pos)] = Vb
    bmap = {col.name: bf[i] for i, col in enumerate(smp.cols)}
    semap = {col.name: float(np.sqrt(max(Vf[i, i], 0))) for i, col in enumerate(smp.cols)}
    if smp.constant:
        bmap["_cons"], semap["_cons"] = bf[-1], float(np.sqrt(max(Vf[-1, -1], 0)))
    if len(lvs) == 1:
        th_names = [f"{t}:_cons" for t in lv.labels()]
    else:       # VERIFY the names of e(b) with several levels
        th_names = [f"{t.replace('_', str(li + 1) + '_', 1)}:_cons" for li, lv_ in enumerate(lvs) for t in lv_.labels()]
    # auxiliary parameters: names in e(b) and the rows of the table are in F.aux_eq
    ball = np.concatenate([bf, aux, th])
    Vall = np.zeros((ball.size, ball.size))
    Vall[:kk, :kk] = Vf
    Vall[kk:kk + naux, kk:kk + naux] = Va
    Vall[kk:kk + naux, kk + naux:] = V[p:p + naux, p + naux:]
    Vall[kk + naux:, kk:kk + naux] = V[p + naux:, p:p + naux]
    Vall[kk + naux:, kk + naux:] = Vth
    bidx = pos + list(range(kk, ball.size))
    V_model_all = np.zeros_like(Vall)
    V_model_all[np.ix_(bidx, bidx)] = V_model
    est = Estimates(cmdname, smp.depname, names + ["_cons"] * (naux + th.size), ball, Vall, n, "z", None)
    est.eqnames = [smp.depname] * kk + aux_eq + [f"/{t.split(':')[0]}" for t in th_names]
    est.rows = rows_from_columns(smp.cols, bmap, semap, ds, constant=smp.constant)
    if off_name:
        est.rows.append(CoefRow("fixed", off_name, eq="(exposure)" if o.get("exposure") else "(offset)"))
    if F.aux_in_table:
        est.rows.append(CoefRow("sep"))
        if F.aux_rows is not None:
            est.rows += F.aux_rows(aux, Va, st.norm.ppf(1 - (1 - smp.level / 100) / 2))
        else:
            for j, nm_ in enumerate(aux_eq):
                se_a = float(np.sqrt(max(Va[j, j], 0)))
                est.rows.append(CoefRow("coef", nm_, float(aux[j]), se_a, name="_cons", noeform=True))
    est.cols, est.constant, est.terms, est.depcomp = smp.cols, smp.constant, smp.terms, smp.depcomp
    est.level = smp.level
    est.vcetype = "Robust" if robust else ""
    est.robust = robust
    if robust:
        est.n_clust = n_clust
        est.clustvar = clustvar or lvs[0].var
    est.offset_name = off_name
    resid = F.resid(aux) if callable(F.resid) else F.resid
    est.glmm = {"model": model, "theta": th, "Vth": Vth, "beta": beta, "levels": lvs, "aux": aux,
                "Va": Va, "aux_eq": aux_eq, "resid_var": resid, "mask": smp.mask.copy(), "nbin": F.nbin,
                "nlev": F.nlev}
    # frequency weights: the observations of a group are counted with their [fweight] and the
    # groups (and the observations) with the fweight() of their level and of the levels above
    # (compat 1916)
    fwo = np.asarray(F.wobs, dtype=float) if (F.wobs is not None and getattr(smp, "wtype", "") == "fweight") \
        else np.ones(n)
    fwg = np.ones(n)
    groups = []
    for li_, lv_ in enumerate(lvs):
        wg_ = F.wgrp[li_] if (F.wgrp and isinstance(F.level_opts[li_].get("fweight"), str)) else None
        if wg_ is not None:
            fwg = fwg * wg_[lv_.gid]
        cnt = np.bincount(lv_.gid, weights=fwo)
        ng = np.bincount(lv_.gid, weights=fwg) / np.bincount(lv_.gid)    # the group's weight
        groups.append((lv_.var, int(round(ng.sum())), int(round(cnt.min())), float(np.sum(cnt * ng) / ng.sum()),
                       int(round(cnt.max()))))
    nvc = sum(lv_.npar() for lv_ in lvs)
    N_ = float(np.sum(fwo * fwg))
    est.extra = {"N": N_, "ll": r.ll, "chi2": wchi2, "p": wp, "df_m": len(slopes), "ll_c": rfe.ll,
                 "chi2_c": chi2_c, "nvc": nvc, "spec": sp, "method": method, "qr": qr,
                 "sd": bool(o.get("stddeviations")), "robust": robust, "n_clust": n_clust,
                 "clustvar": clustvar or lvs[0].var, "groups": groups}
    est.extra["npts"] = 1 if method == "laplace" else npts
    est.display = display_qr if qr else display_me
    est.predict = F.predict or predict_me
    est.predict_opts = F.predict_opts if F.predict_opts is not None else \
        {"mu": 2, "eta": 3, "reffects": 4, "fixedonly": 5, "pr": 2, "outcome": 3, "conditional": 4}
    est.predict_multi = {"reffects": reffects_me, **F.predict_multi}
    est.estat = {"icc": (3, estat_icc), "group": (3, estat_group)}
    if qr:
        est.estat["recovariance"] = (5, estat_recov)       # compat 1902: unrecognized after the others
    est.icc_ok = cmdname in ("melogit", "meprobit", "meqrlogit", "meqrpoisson") or \
        (cmdname == "meglm" and fam in ("bernoulli", "gaussian"))   # VERIFY for meglm
    est.eform_default = F.eform
    # e() in the order Stata lists it (compat 1002 and 1902)
    nvar_p = sum({"identity": 1, "exchangeable": 1}.get(lv_.cov, lv_.q) for lv_ in lvs)
    ncov_p = nvc - nvar_p
    from ..commands.matrix import Matrix
    from .fvars import base_vars
    from .ml import ilog
    from .results import datasignature
    covs = [c.name for c in smp.cols if c.values is not None]
    family = fam
    names_all = list(est.names)
    grad = np.zeros((1, len(names_all)))
    if r.g is not None and len(r.g) == len(names_all):
        grad[0] = r.g
    ivars = " ".join(lv_.var for lv_ in lvs)
    gmats = [(nm, Matrix(np.array([[g_[j] for g_ in groups]], dtype=float), ["r1"], [lv_.var for lv_ in lvs]))
             for nm, j in (("g_max", 4), ("g_avg", 3), ("g_min", 2), ("N_g", 1))]
    cmdline_ = F.cmdline or cmdline(cmdname, args)
    if qr:
        # meqrlogit and meqrpoisson: the e() of the xtmelogit family. The plain fit lists it in the
        # order below and a fit with refineopts() in the reverse order (compat 1902; VERIFY)
        sigvars = [smp.depname] + [v for v in base_vars(smp.terms) if v != smp.depname] + [lv.var] \
            + [nm for nm in lv.names if nm != "_cons"]
        scalars = [("converged", 1.0 if r.converged else 0.0), ("rc", 0.0), ("p", wp), ("chi2", wchi2),
                   ("df_m", float(len(slopes))), ("N", float(n)), ("p_c", _p_chibar(chi2_c, nvc)),
                   ("chi2_c", chi2_c), ("df_c", float(nvc)), ("ll_c", rfe.ll), ("ll", r.ll),
                   ("k_rc", float(ncov_p)), ("k_rs", float(nvar_p)), ("k_r", float(th.size)),
                   ("k_f", float(len(slopes) + int(smp.constant))), ("k", float(ball.size)),
                   ("reparm_rc", 0.0), ("rank", float(ball.size))]
        macros = [("cmdline", cmdline_), ("datasignaturevars", " ".join(sigvars)),
                  ("datasignature", datasignature(s, sigvars, smp.mask)), ("chi2type", "Wald"),
                  ("n_quad", str(int(est.extra["npts"]))), ("redim", " ".join(str(l_.q) for l_ in lvs)),
                  ("ivars", " ".join(l_.var for l_ in lvs)),   # nested levels: VERIFY
                  ("vartypes", " ".join(COV_TITLE[l_.cov] for l_ in lvs)),
                  ("revars", " ".join(" ".join(l_.names) for l_ in lvs)), ("title", sp["title"]),
                  ("technique", "nr"), ("ml_method", "d0"), ("opt", "moptimize"), ("model", sp["model"]),
                  ("estat_cmd", f"{cmdname}_estat"), ("predict", f"{cmdname}_p"), ("cmd", cmdname),
                  ("marginsdefault", "predict(xb)"),
                  ("marginsnotok", "stdp Pearson Deviance Anscombe REFfects RESEs"), ("method", "ML"),
                  ("depvar", smp.depname), ("properties", "b V")]
        extra_mats = list(gmats)
        if o.get("refineopts"):
            scalars.reverse()
            macros.reverse()
            extra_mats.reverse()
        if method == "laplace":
            macros.append(("laplace", "laplace"))
    else:
        # the level variables from the lowest level up (compat 1916)
        sigvars = [lv_.var for lv_ in reversed(lvs)] + [v for v in base_vars(smp.terms) if v != smp.depname] + [smp.depname]
        ordinal = fam == "ordinal"
        scalars = [("rank", float(ball.size)), ("N", N_), ("ic", float(r.iterations)),
                   ("converged", 1.0 if r.converged else 0.0), ("rc", 0.0), ("ll", r.ll),
                   ("k", float(ball.size)), ("k_eq", float(1 + naux + nvc))]
        if ordinal:
            scalars.append(("k_cat", float(F.nlev)))
        scalars += [("k_eq_model", 1.0), ("ll_c", rfe.ll), ("k_rs", float(nvar_p)), ("k_rc", float(ncov_p)),
                    ("k_dv", 1.0)]
        if len(lvs) > 1:        # compat 1916: with nested levels this part comes in reverse order
            scalars.reverse()
        scalars += [("k_r", float(sum(lv_.q for lv_ in lvs))),
                    ("k_f", float(len(slopes)))]   # VERIFY: k_f without _cons
        scalars += [("p_c", _p_chibar(chi2_c, nvc)), ("df_c", float(nvc)), ("chi2_c", chi2_c)]
        scalars += [("p", wp), ("chi2", wchi2), ("df_m", float(len(slopes)))]
        if robust:
            scalars += [("N_clust", float(n_clust))]
        common = [("n_quad", str(int(est.extra["npts"]))), ("ivars", ivars), ("chi2type", "Wald"),
                  ("datasignaturevars", " ".join(sigvars)),
                  ("datasignature", datasignature(s, sigvars, smp.mask)), ("intmethod", method),
                  ("depvar", smp.depname)] + \
            ([("clustvar", clustvar or lvs[0].var)] if robust else []) + \
            [("marginsnotok", "LATent1(passthru) LATent se(passthru)"), ("covariates", " ".join(covs))]
        tail = [("opt", "moptimize")] + ([("vcetype", "Robust")] if robust else []) + [("vce", vce)] + \
            [("user", "_gsem_eval_chol()"), ("ml_method", "gf2"),
             ("technique", "nr"), ("which", "max"), ("properties", "b V")]   # VERIFY vcetype/clustvar
        if getattr(smp, "wtype", ""):
            from .regress import wexp_macros
            tail[1:1] = wexp_macros(smp)   # VERIFY the place
        if ordinal:
            tail.insert(0, ("marginsdefault", " ".join(f"predict(pr outcome({c:g}))" for c in F.cats)))
        if sp["cmd"] == "meglm":
            macros = [("cmdline", cmdline_), ("cmd2", "meglm"), ("cmd", "gsem"), ("title", sp["title"]),
                      ("family", family), ("model", sp["model"]), ("predict", "meglm_p"),
                      ("estat_cmd", "meglm_estat")] + common + [("link", link)] + tail
        else:
            macros = [("cmd", "meglm"), ("estat_cmd", f"{cmdname}_estat"), ("predict", f"{cmdname}_p"),
                      ("cmdline", cmdline_), ("cmd2", cmdname), ("model", sp["model"]), ("title", sp["title"]),
                      ("dispersion", sp["disp"]), ("link", link), ("family", family)] + common + tail
            if cmdname == "meoprobit":      # compat 1902
                d_ = dict(macros)
                macros = [(k_, d_[k_]) for k_ in ("cmd", "cmdline", "cmd2", "estat_cmd", "predict", "model",
                                                  "title", "link", "family")] + \
                         [m for m in macros if m[0] not in ("cmd", "cmdline", "cmd2", "estat_cmd", "predict",
                                                            "model", "title", "link", "family")]
        extra_mats = list(gmats)
        if ordinal:
            extra_mats.append(("cat", Matrix(np.array([F.cats], dtype=float), ["r1"],
                                             [f"c{i + 1}" for i in range(len(F.cats))])))
        if robust:
            extra_mats.append(("V_modelbased", Matrix(V_model_all, names_all, names_all)))
        extra_mats += [("gradient", Matrix(grad, ["r1"], names_all)),
                       ("ilog", Matrix(ilog(r), ["r1"], [f"c{i + 1}" for i in range(20)]))]
    if F.macro_edit is not None:
        macros = F.macro_edit(macros)
    post(s, est, scalars, macros, smp.mask, extra_mats=extra_mats)
    if F.late_error:
        raise StataError(198, F.late_error)
    if qr:
        display_qr(s, est, header=not o.get("noheader"), group=not o.get("nogroup"),
                    lrtest=not o.get("nolrtest"), table=not o.get("nofetable"),
                    retable=not o.get("noretable"))
    else:
        display_me(s, est, header=not o.get("noheader"), group=not o.get("nogroup"),
                   lrtest=not o.get("nolrtest"))


def cmd_melogit(s, args):
    cmd_me(s, args, "melogit")


def cmd_meprobit(s, args):
    cmd_me(s, args, "meprobit")


def cmd_mepoisson(s, args):
    cmd_me(s, args, "mepoisson")


def cmd_mecloglog(s, args):
    cmd_me(s, args, "mecloglog")


def cmd_menbreg(s, args):
    cmd_me(s, args, "menbreg")


def cmd_meologit(s, args):
    cmd_me(s, args, "meologit")


def cmd_meoprobit(s, args):
    cmd_me(s, args, "meoprobit")


def cmd_meglm(s, args):
    cmd_me(s, args, "meglm")


def cmd_meqrlogit(s, args):
    cmd_me(s, args, "meqrlogit")


def cmd_meqrpoisson(s, args):
    cmd_me(s, args, "meqrpoisson")
