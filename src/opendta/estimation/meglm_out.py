"""Output and prediction after meglm and its special cases: the header, the table of fixed and
random effects (also in the layout of the meqr commands), and predict (mu, eta, the random
effects by empirical Bayes, the probabilities of the ordinal outcomes)."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..lazy import lazy
from .mixed import _right, group_header, lr_line
from .results import CoefRow, Estimates, coef_table, g9

# imported on first use (start-up time; see opendta/lazy.py)
special = lazy("scipy.special")
st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


_FAM_TITLE = {"gaussian": "Gaussian", "bernoulli": "Bernoulli", "binomial": "binomial",
              "gamma": "gamma", "nbinomial": "negative binomial", "ordinal": "ordinal",
              "poisson": "Poisson"}                                  # compat 1902


_LINK_TITLE = {"identity": "identity", "logit": "logit", "probit": "probit", "log": "log",
               "cloglog": "complementary log-log"}                   # compat 1902


def mean_fn(sp: dict):
    """Inverse link of the family: mean of the outcome given the linear predictor."""
    link = sp["link"]
    if link == "logit":
        return special.expit
    if link == "probit":
        return st.norm.cdf
    if link == "cloglog":
        return lambda e: -np.expm1(-np.exp(np.clip(e, -700, 6.5)))
    if link == "log":
        return lambda e: np.exp(np.clip(e, -700, 700))
    return lambda e: e


def display_me(s: "Session", est: Estimates, *, header: bool = True, table: bool = True,
               level: float | None = None, group: bool = True, lrtest: bool = True, eform=None, **_kw) -> None:
    out = s.output
    x = est.extra
    sp = x["spec"]
    lev = level or est.level
    if header:
        ll = f"Log {'pseudolikelihood' if x.get('robust') else 'likelihood'} = {g9(x['ll'], '%10.0g'):>10}"
        meth = x["method"]
        mid = []
        if sp["cmd"] == "meglm":
            mid = [f"{'Family:':<7}{_FAM_TITLE[sp['fam']]:>24}", f"{'Link:':<5}{_LINK_TITLE[sp['link']]:>26}"]
        if sp["fam"] == "nbinomial":
            mid.append(f"{'Overdispersion:':<15}{sp['disp']:>16}")      # VERIFY
        if sp["fam"] == "binomial":
            nb_ = str(sp["nbin"])
            lab_ = "Binomial trials:" if nb_.replace(".", "").isdigit() else "Binomial variable:"
            mid.append(f"{lab_}{nb_:>{31 - len(lab_)}}")
        meth_line = (f"{'Integration method:':<19}{meth:>12}", "", "") if meth == "laplace" else \
            (f"{'Integration method:':<19}{meth:>12}", "Integration pts.", f"{int(x['npts'])}")
        group_header(out, sp["title"], x["groups"] if group else [], x["N"],
                     [("", "", ""),
                      meth_line,
                      ("", "", ""),
                      ("", f"Wald chi2({int(x['df_m'])})", g9(x["chi2"], "%10.2f")),
                      (ll, "Prob > chi2", g9(x["p"], "%10.4f"))], gv_right=True, mid=mid)
    if x.get("robust") and header:
        note = f"(Std. Err. adjusted for {int(x['n_clust']):,} clusters in {x['clustvar']})"   # VERIFY
        out.write(f"{note:>78}\n", "text")
    if not table:
        return
    z = st.norm.ppf(1 - (1 - lev / 100) / 2)
    ef = est.eform_default if eform is None else eform
    rows = list(est.rows)
    G = est.glmm
    from .mixed import _vc_rows
    sd = bool(x.get("sd"))
    k0 = 0
    for li, lv in enumerate(G["levels"]):
        npar = lv.npar()
        t_ = G["theta"][k0:k0 + npar]
        V_ = G["Vth"][k0:k0 + npar, k0:k0 + npar]
        k0 += npar
        # the label of a nested level joins the levels above it: school>class (VERIFY)
        eqlab = ">".join(l_.var for l_ in G["levels"][:li + 1])
        vc = _vc_rows(lv, t_, V_, sd, z)
        rows.append(CoefRow("sep"))
        rows.append(CoefRow("eq", eqlab))
        for lab, v, se, lo, hi in vc:
            if not lab.startswith("cov("):      # stddeviations: sd() and corr() rows, VERIFY the layout
                rows.append(CoefRow("aux", lab, v, se, eq=f"{float(lo)!r} {float(hi)!r}", noeform=True,
                                    tight=True))
        covs = [r_ for r_ in vc if r_[0].startswith("cov(")]
        if covs:
            # covariances: a section of their own, with z tests, and the pair written in reverse
            # order of e(revars) (compat 1002: cov(_cons,x) after var(x) and var(_cons))
            rows.append(CoefRow("sep"))
            rows.append(CoefRow("eq", eqlab))
            for lab, v, se, _lo, _hi in covs:
                a, b_ = lab[4:-1].split(",", 1) if "," in lab else (lab[4:-1], "")
                lab2 = f"cov({b_},{a})" if b_ else lab
                rows.append(CoefRow("coef", lab2, v, se, noeform=True, tight=True))
    if sp["fam"] == "gaussian":
        # residual variance, after the variances of the levels (VERIFY the layout)
        a_, se_a = float(G["aux"][0]), float(np.sqrt(max(G["Va"][0, 0], 0)))
        h_ = 0.5 if sd else 1.0
        v_ = math.exp(h_ * a_)
        rows.append(CoefRow("sep"))
        rows.append(CoefRow("aux", f"{'sd' if sd else 'var'}(e.{est.depvar})", v_, v_ * h_ * se_a,
                            eq=f"{math.exp(h_ * (a_ - z * se_a))!r} {math.exp(h_ * (a_ + z * se_a))!r}",
                            noeform=True, tight=True))
    coef_table(s, est.depvar, rows, stat="z", df=None, level=lev,
               coef_title=ef or "Coef.", eform=bool(ef), vcetype="Robust" if x.get("robust") else "")
    if lrtest and not x.get("robust"):
        lr_line(out, sp["vs"], x["chi2_c"], x["nvc"])


COV_TITLE = {"identity": "Identity", "independent": "Independent", "exchangeable": "Exchangeable",
       "unstructured": "Unstructured"}


def display_qr(s: "Session", est: Estimates, *, header: bool = True, table: bool = True,
                level: float | None = None, group: bool = True, lrtest: bool = True, eform=None,
                retable: bool = True, **_kw) -> None:
    """meqrlogit and meqrpoisson: the layout of mixed (header with the integration points, table
    of random-effects parameters); VERIFY against the Stata run."""
    from .mixed import _re_row
    out = s.output
    x = est.extra
    sp = x["spec"]
    lev = level or est.level
    if header:
        ll = f"Log likelihood = {g9(x['ll'], '%10.0g'):>10}"
        if len(x["groups"]) > 1 and group:
            # nested levels: the group table has a column of integration points (compat 1916)
            out.write("\n", "text")
            _right(out, sp["title"], "Number of obs", f"{int(x['N']):,}")
            out.write("\n" + "-" * 76 + "\n", "text")
            out.write("                |     No. of       Observations per Group       Integration\n", "text")
            out.write(" Group Variable |     Groups    Minimum    Average    Maximum      Points\n", "text")
            out.write("----------------+" + "-" * 59 + "\n", "text")
            for var, ng, mn, av, mx in x["groups"]:
                out.write(f"{var[:14]:>15} |", "text")
                out.write(f" {ng:>10,} {mn:>10,} {av:>10.1f} {mx:>10,} {int(x['npts']):>11}\n", "result")
            out.write("-" * 76 + "\n", "text")
            for left, lab_, val in [("", "", ""), ("", f"Wald chi2({int(x['df_m'])})", g9(x["chi2"], "%10.2f")),
                                    (ll, "Prob > chi2", g9(x["p"], "%10.4f"))]:
                if lab_:
                    _right(out, left, lab_, val)
                else:
                    out.write(left + "\n", "text")
        else:
            group_header(out, sp["title"], x["groups"] if group else [], x["N"],
                         [("", "", ""),
                          (f"Integration points = {int(x['npts']):>3}", f"Wald chi2({int(x['df_m'])})",
                           g9(x["chi2"], "%10.2f")),
                          (ll, "Prob > chi2", g9(x["p"], "%10.4f"))])
    if table:
        ef = est.eform_default if eform is None else eform
        coef_table(s, est.depvar, est.rows, stat="z", df=None, level=lev,
                   coef_title=ef or "Coef.", eform=bool(ef))
    z = st.norm.ppf(1 - (1 - lev / 100) / 2)
    G = est.glmm
    if retable:
        from .mixed import _vc_rows
        out.write("\n" + "-" * 78 + "\n", "text")
        out.write(f"  Random-effects Parameters  |   Estimate   Std. Err.     [{lev:g}% Conf. Interval]\n", "text")
        k0 = 0
        for lv in G["levels"]:
            k1 = k0 + lv.npar()
            out.write("-" * 29 + "+" + "-" * 48 + "\n", "text")
            kind = {"identity": "Identity", "independent": "Independent", "exchangeable": "Exchangeable",
                    "unstructured": "Unstructured"}[lv.cov]
            out.write(f"{lv.var + ': ' + kind:<29}|\n", "text")
            for row in _vc_rows(lv, G["theta"][k0:k1], G["Vth"][k0:k1, k0:k1], x["sd"], z):
                _re_row(out, *row)
            k0 = k1
        out.write("-" * 78 + "\n", "text")
    if lrtest:
        lr_line(out, sp["vs"], x["chi2_c"], x["nvc"])
    if x["method"] == "laplace":      # compat 1902
        out.write("\nNote: Log-likelihood calculations are based on the Laplacian approximation.\n", "text")


def posterior_eta(s: "Session", est: Estimates, fixedonly: bool = False):
    """Linear predictor for every observation of the data: fixed part (with the offset) plus the
    posterior means of the random effects of the group (zero when the group is not in the
    estimation sample)."""
    from .postest import xb_all
    xb, ok = xb_all(s, est)
    ds = s.data
    G = est.glmm
    model = G["model"]
    if est.offset_name:
        nm = est.offset_name
        v = ds.get(nm[3:-1] if nm.startswith("ln(") else nm).data.astype(float)
        xb = xb + np.where(v < SYS, np.log(np.where(v > 0, v, 1)) if nm.startswith("ln(") else v, 0)
    if fixedonly:
        return xb, ok
    b = np.concatenate([G["beta"], G["aux"], G["theta"]])
    model.adapt(b, steps=50)
    eta = xb.copy()
    for li, lv in enumerate(G["levels"]):
        gidx = _group_index(s, G, li)
        re_all = np.where(gidx[:, None] >= 0, model.mus[li][np.maximum(gidx, 0)], np.nan)
        Zcols = []
        for nm in lv.names:
            Zcols.append(np.ones(ds.nobs) if nm == "_cons" else ds.get(nm).data.astype(float))
        Z = np.column_stack(Zcols)
        eta = eta + np.where(np.isfinite(re_all).all(axis=1), np.einsum("nq,nq->n", Z, np.nan_to_num(re_all)), 0.0)
    return eta, ok


def _group_index(s: "Session", G: dict, li: int) -> np.ndarray:
    """The group of level li of every observation of the data (-1 when the group is not in the
    estimation sample); a nested group is the combination of the values of the levels above."""
    ds = s.data
    rows = np.nonzero(G["mask"])[0]
    lvs = G["levels"][:li + 1]
    cols = [ds.get(lv.var).data for lv in lvs]
    keys = list(zip(*cols))
    key_of = {keys[i]: g for i, g in zip(rows, G["levels"][li].gid)}
    return np.array([key_of.get(k, -1) for k in keys])


def predict_me(s: "Session", est: Estimates, stat: str, mask: np.ndarray):
    sp = est.extra["spec"]
    G = est.glmm
    cond = str((getattr(est, "_pred_opts", None) or {}).get("conditional") or "").strip().lower()
    if cond and not ("fixedonly".startswith(cond) or "ebmeans".startswith(cond)):
        raise StataError(198, f"conditional({cond}) invalid")   # ebmodes: VERIFY, not available
    if cond and "fixedonly".startswith(cond):
        stat = "fixedonly" if stat in ("", "mu") else stat
        if stat == "eta":
            stat = "xb"
    if stat == "xb":
        eta, ok = posterior_eta(s, est, fixedonly=True)
        return np.where(ok, eta, SYS), ""
    eta, ok = posterior_eta(s, est, fixedonly=stat == "fixedonly")
    if stat == "eta":
        if est.extra.get("qr") or (cond and "fixedonly".startswith(cond)):
            return np.where(ok, eta, SYS), ""
        return np.where(ok, eta, SYS), f"(predictions based on fixed effects and posterior means of random effects)\n(using {int(est.extra['npts'])} quadrature points)"
    mu = mean_fn(sp)(eta)
    if sp["fam"] == "binomial":
        t = sp["nbin"]
        try:
            mu = mu * float(t)
        except ValueError:
            mu = mu * s.data.get(t).data.astype(float)
    # notes (compat 1902): fixedonly alone says only that mu was assumed; the meqr commands say
    # "predicted means"; conditional(fixedonly) is silent
    base = "(predictions based on fixed effects and posterior means of random effects)"
    quad = f"(using {int(est.extra['npts'])} quadrature points)"
    lines = []
    if est.extra.get("qr"):
        if stat in ("", "mu"):
            lines = ["(option mu assumed; predicted means)"] if not stat else []
    elif cond and "fixedonly".startswith(cond):
        lines = []
    elif stat == "fixedonly":
        lines = ["(option mu assumed)"]
    else:
        lines = [base] + (["(option mu assumed)"] if stat in ("", "outcome") else []) + [quad]
    note = "\n".join(lines)
    return np.where(ok, mu, SYS), note


def pr_multi(s: "Session", est: Estimates):
    """predict, pr after meologit/meoprobit: the probability of each outcome (or of outcome(#))
    given the fixed effects and the posterior means of the random effects (VERIFY)."""
    G = est.glmm
    sp = est.extra["spec"]
    J = G["nlev"]
    oc = getattr(est, "_pred_opts", {}).get("outcome")
    nreq = 1        # compat 1902: pr takes one new variable ("too many variables specified" otherwise)
    s.output.write("(predictions based on fixed effects and posterior means of random effects)\n", "text")
    s.output.write(f"(using {int(est.extra['npts'])} quadrature points)\n", "text")

    def fn(mask):
        eta, ok = posterior_eta(s, est)
        cuts = np.concatenate([[-np.inf], G["aux"], [np.inf]])
        F = special.expit if sp["link"] == "logit" else special.ndtr
        ps = [np.where(ok, F(cuts[j + 1] - eta) - F(cuts[j] - eta), SYS) for j in range(J)]
        if oc:
            t = str(oc).strip().lstrip("#")
            k = int(float(t)) - 1 if str(oc).strip().startswith("#") else None
            if k is None:
                cats = [c for c in range(J)]       # outcome(value): VERIFY with non-consecutive values
                k = int(float(t)) - 1 if 0 <= int(float(t)) - 1 < len(cats) else 0
            return [ps[k]]
        return ps[:1]       # pr with one variable: the first outcome (VERIFY)
    return nreq, fn, "too few variables specified"


def reffects_me(s: "Session", est: Estimates):
    """predict newvarlist, reffects: posterior means, one variable per random effect in the
    order of the level's effects (VERIFY: the order for random slopes)."""
    G = est.glmm
    nreq = sum(lv.q for lv in G["levels"])
    if not est.extra.get("qr"):
        s.output.write("(calculating posterior means of random effects)\n", "text")

    def fn(mask):
        model = G["model"]
        model.adapt(np.concatenate([G["beta"], G["aux"], G["theta"]]), steps=50)
        out = []
        # the highest level first, within a level the order of its effects (VERIFY)
        for li, lv in enumerate(G["levels"]):
            gidx = _group_index(s, G, li)
            for j in range(lv.q):
                vals = np.where(gidx >= 0, model.mus[li][np.maximum(gidx, 0), j], SYS)
                out.append(np.where(np.isfinite(vals), vals, SYS))
        if not est.extra.get("qr"):
            s.output.write(f"(using {int(est.extra['npts'])} quadrature points)\n", "text")
        return out
    # Stata's message here has a single line (compat 1002)
    return nreq, fn, "too few variables specified"
