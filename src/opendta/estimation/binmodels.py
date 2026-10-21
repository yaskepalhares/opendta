"""Binary and count models with an ancillary parameter or a second equation:
scobit, hetprobit, fracreg, zip, zinb, tpoisson, tnbreg, gnbreg and binreg
([R] scobit, [R] hetprobit, [R] fracreg, [R] zip, [R] zinb, [R] tpoisson, [R] tnbreg,
[R] gnbreg, [R] binreg). The log likelihoods are in likelihoods.py; the layout of the
equations in mleq.py."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.stats import chi2_tail
from ..lazy import lazy
from . import likelihoods as L
from .ml import maximize
from .mleq import (Engine, aux_row, block_from_pr, block_from_text, constant_only, display_eq,
                   finish_model, fit_log, lab, make_predict, names_in, offset_values, peek, std_right,
                   sub_options, transformed_row)
from .models import _logit, _poisson, _probit
from .models2 import _header, _Prep, _log, g
from .regress import comma
from .results import CoefRow

special = lazy("scipy.special")
st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def _replay(s, cmd, args) -> bool:
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, cmd, args)
        return True
    return False


def _wsum(w, mask) -> float:
    return float(np.sum(w[mask]))


def _binary_y(pr) -> np.ndarray:
    y = (pr.y != 0).astype(float)
    if np.all(y == y[0]):
        raise StataError(2000, "outcome does not vary; remember:\n"
                         "                                  0 = negative outcome,\n"
                         "        all other nonmissing values = positive outcome")
    return y


def _count_y(pr) -> np.ndarray:
    y = pr.y
    if np.any(y < 0) or np.any(y != np.trunc(y)):
        raise StataError(459, f"{pr.smp.depname} must be a nonnegative integer")   # VERIFY
    return y


def _set_offset(blk, pr, off_opt=None, exp_opt=None) -> None:
    blk.offset, blk.off_name = offset_values(pr, off_opt, exp_opt)


def _chi2p(x: float, df: int) -> float:
    return float(chi2_tail(x, df)) if df > 0 else SYS


def _null_fit(pr, lik, blocks, w, n_theta, b_full):
    """Maximum of the model with the constant only in the first equation (the other equations
    and ancillary parameters free, started at the full model's values)."""
    first = blocks[0]
    if not first.constant:
        return None
    nb = [constant_only(first)] + list(blocks[1:])
    eng0 = Engine(lik, nb, w, n_theta)
    b0 = np.concatenate([[b_full[first.k - 1]], b_full[first.k:]])
    return maximize(eng0.fun, b0, maxiter=300)


def _hdr_logit_style(extra_right=(), left_extra=()):
    def hdr(est):
        x = est.extra
        left = [est.title] + list(left_extra) + [""] * max(0, 2 - len(left_extra)) + \
            [f"{lab(est)} = {g(x['ll'], '%10.0g'):>10}"]
        return left, std_right(x, extra_right)
    return hdr


def _poisson_start(y, X, w, off):
    pf, _ = _poisson(y, X, w, off if off is not None else np.zeros(len(y)))
    ybar = float(np.sum(w * y) / np.sum(w))
    b0 = np.zeros(X.shape[1])
    return maximize(pf, b0), ybar


def _alpha_start(y, mu, w) -> float:
    """Method-of-moments overdispersion for starting values."""
    m = float(np.sum(w * mu))
    a = (float(np.sum(w * ((y - mu) ** 2 - mu))) / float(np.sum(w * mu ** 2))) if m > 0 else 0.5
    return float(np.clip(a, 0.05, 5.0))


# ---------------------------------------------------------------------------
# scobit
# ---------------------------------------------------------------------------

def cmd_scobit(s: "Session", args: str) -> None:
    if _replay(s, "scobit", args):
        return
    o = peek(args, {"offset": 3, "or": 2, "asis": 4})
    pr = _Prep(s, args, {"offset": 3, "or": 2, "asis": 4}, extra_vars=names_in(o.get("offset")))
    pr.write_notes(s.output)
    y = _binary_y(pr)
    main = block_from_pr(pr, pr.smp.depname)
    _set_offset(main, pr, o.get("offset"))
    w = pr.w
    off = main.offset if main.offset is not None else np.zeros(len(y))
    # starting values: the logit fit (alpha = 1), whose log likelihood is the model of the LR test
    rl = maximize(_offset_logit(y, main.X, w, off), np.zeros(main.k))
    rl_res = rl

    def lik(etas, th):
        ll, d_eta, d_la = L.scobit(y, etas[0], th[0])
        return ll, [d_eta], d_la.reshape(-1, 1)

    eng = Engine(lik, [main], w, 1)
    res = pr.maximize(eng.fun, np.concatenate([rl.b, [0.0]]))
    fit_log(s, pr, "Fitting logistic model:", rl_res)
    fit_log(s, pr, "Fitting full model:", res)
    V, n_clust = pr.vce(res.H, eng.scores(res.b) if pr.robust else None)
    lna, se_lna = float(res.b[-1]), float(np.sqrt(V[-1, -1]))
    lr = max(2 * (res.ll - rl.ll), 0.0)
    p_lr = _chi2p(lr, 1)
    n1, n0 = _wsum(w, y > 0), _wsum(w, y == 0)
    null = _null_fit(pr, lik, [main], w, 1, res.b)

    def aux_rows(level):
        return [CoefRow("sep"), aux_row("/lnalpha", lna, se_lna, level), CoefRow("sep"),
                transformed_row("alpha", np.exp, np.exp, lna, se_lna, level)]

    def hdr(est):
        x = est.extra
        ll4 = g(x["ll"], "%9.0g")
        left = [est.title, "", f"{lab(est)} = {ll4}"]
        right = [("Number of obs", comma(x["N"])), ("Zero outcomes", comma(n0)),
                 ("Nonzero outcomes", comma(n1))]
        if est.robust:     # VERIFY: with vce(robust) the Wald test is shown
            left = [est.title, "", "", "", f"{lab(est)} = {ll4}"]
            right += [(f"Wald chi2({int(x['df_m'])})", g(x["chi2"], "%10.2f")), ("Prob > chi2", g(x["p"], "%10.4f"))]
        return left, right

    def footer(s_, est):
        s_.output.write(f"Likelihood-ratio test of alpha=1:   chi2(1) = {lr:9.2f}    Prob > chi2 = {p_lr:6.4f}\n",
                        "text")
        s_.output.write("\nNote: Likelihood-ratio tests are recommended for inference with scobit models.\n", "text")

    est = finish_model(s, pr, "scobit", args, title="Skewed logistic regression", blocks=[main], eng=eng,
                       res=res, V=V, n_clust=n_clust, ll0=null.ll if null else SYS,
                       aux_names=[("lnalpha", "_cons")], aux_rows_fn=aux_rows, hdr_fn=hdr, wald=True,
                       footer=footer, k_eq=2, extra_scalars=[("alpha", float(np.exp(lna)))])
    est.predict = make_predict(lambda what, etas, th, e: (1 - (1 + np.exp(etas[0])) ** -np.exp(th[0]))
                               if what in ("pr", "p") else None, "pr")
    est.pred_label = f"Pr({est.depvar})"
    est.coef_title, est.eform_on = ("Odds Ratio", True) if o.get("or") else ("Coef.", False)
    display_eq(s, est)
    est.coef_title, est.eform_on = "Coef.", False


def _offset_logit(y, X, w, off):
    def fun(b, todo):
        xb = X @ b + off
        ll = float(np.sum(w * np.where(y > 0, -np.logaddexp(0, -xb), -np.logaddexp(0, xb))))
        if todo == 0:
            return ll, None, None
        p = special.expit(xb)
        return ll, X.T @ (w * (y - p)), -(X * (w * p * (1 - p))[:, None]).T @ X
    return fun


# ---------------------------------------------------------------------------
# hetprobit
# ---------------------------------------------------------------------------

def cmd_hetprobit(s: "Session", args: str) -> None:
    if _replay(s, "hetprobit", args):
        return
    spec = {"het": 3, "offset": 3, "asis": 4}
    o = peek(args, spec)
    if o.get("het") in (None, True):
        raise StataError(198, "option het() required")
    het_vars, het_opts = sub_options(str(o["het"]), {"offset": 3})
    pr = _Prep(s, args, spec, extra_vars=names_in(het_vars, o.get("offset"), het_opts.get("offset")))
    pr.write_notes(s.output)
    y = _binary_y(pr)
    main = block_from_pr(pr, pr.smp.depname)
    _set_offset(main, pr, o.get("offset"))
    het = block_from_text(pr, "lnsigma2", het_vars, constant=False)
    _set_offset(het, pr, het_opts.get("offset"))
    for nm in het.notes:
        s.output.write(f"note: {nm} omitted because of collinearity\n", "text")
    w = pr.w
    off = main.offset if main.offset is not None else np.zeros(len(y))
    # the homoskedastic probit is both the starting point and the model of the LR test
    pf, _ = _probit(y, main.X, w)
    if main.offset is not None:
        pf = _offset_probit(y, main.X, w, off)
    rp = maximize(pf, np.zeros(main.k))

    def lik(etas, th):
        ll, d1, d2 = L.hetprobit(y, etas[0], etas[1])
        return ll, [d1, d2], None

    eng = Engine(lik, [main, het], w, 0)
    res = pr.maximize(eng.fun, np.concatenate([rp.b, np.zeros(het.k)]))
    fit_log(s, pr, "Fitting probit model:", rp)
    fit_log(s, pr, "Fitting full model:", res)
    V, n_clust = pr.vce(res.H, eng.scores(res.b) if pr.robust else None)
    lr = max(2 * (res.ll - rp.ll), 0.0)
    p_lr = _chi2p(lr, het.k)
    n1, n0 = _wsum(w, y > 0), _wsum(w, y == 0)

    def hdr(est):
        x = est.extra
        ll4 = g(x["ll"], "%9.0g")
        left = [est.title, "", "", "", "", f"{lab(est)} = {ll4}"]
        right = [("Number of obs", comma(x["N"])), ("Zero outcomes", comma(n0)),
                 ("Nonzero outcomes", comma(n1)), ("", ""), (f"Wald chi2({int(x['df_m'])})", g(x["chi2"], "%10.2f")),
                 ("Prob > chi2", g(x["p"], "%10.4f"))]
        return left, right

    def footer(s_, est):
        s_.output.write(f"Likelihood-ratio test of lnsigma2=0: chi2({het.k}) = {lr:8.2f}   Prob > chi2 = {p_lr:6.4f}\n",
                        "text")

    est = finish_model(s, pr, "hetprobit", args, title="Heteroskedastic probit model", blocks=[main, het],
                       eng=eng, res=res, V=V, n_clust=n_clust, ll0=SYS, aux_names=[],
                       aux_rows_fn=lambda lv: [], hdr_fn=hdr, wald=True, footer=footer, k_eq=2,
                       extra_scalars=[("chi2_c", lr), ("p_c", p_lr)])
    est.predict = make_predict(lambda what, etas, th, e: st.norm.cdf(etas[0] * np.exp(-etas[1]))
                               if what in ("pr", "p") else None, "pr")
    est.pred_label = f"Pr({est.depvar})"
    display_eq(s, est)


def _offset_probit(y, X, w, off):
    q = np.where(y > 0, 1.0, -1.0)

    def fun(b, todo):
        xb = X @ b + off
        ll = float(np.sum(w * st.norm.logcdf(q * xb)))
        if todo == 0:
            return ll, None, None
        lam = q * np.exp(st.norm.logpdf(q * xb) - st.norm.logcdf(q * xb))
        return ll, X.T @ (w * lam), -(X * (w * lam * (lam + xb))[:, None]).T @ X
    return fun


# ---------------------------------------------------------------------------
# fracreg
# ---------------------------------------------------------------------------

_FRAC_LINKS = {"probit": "Fractional probit regression", "logit": "Fractional logistic regression",
               "cloglog": "Fractional complementary log-log regression"}


def cmd_fracreg(s: "Session", args: str) -> None:
    words = args.strip().split(None, 1)
    if not words or words[0].startswith(","):
        _replay(s, "fracreg", args)
        return
    link = words[0].lower()
    if link not in _FRAC_LINKS:
        raise StataError(198, f"{words[0]} invalid link")   # VERIFY
    rest = words[1] if len(words) > 1 else ""
    if _replay(s, "fracreg", rest) and not rest.strip():
        return
    spec = {"offset": 3}
    o = peek(rest, spec)
    pr = _Prep(s, rest, spec, extra_vars=names_in(o.get("offset")))
    pr.write_notes(s.output)
    if not pr.robust:      # fracreg reports robust standard errors by default
        pr.robust, pr.smp.vce = True, "robust"
    y = pr.y
    if np.any(y < 0) or np.any(y > 1):
        raise StataError(459, f"{pr.smp.depname} must be between 0 and 1 inclusive")   # VERIFY
    main = block_from_pr(pr, pr.smp.depname)
    _set_offset(main, pr, o.get("offset"))
    w = pr.w

    def lik(etas, th):
        ll, d = L.fractional(y, etas[0], link)
        return ll, [d], None

    eng = Engine(lik, [main], w, 0)
    ybar = float(np.clip(np.sum(w * y) / np.sum(w), 1e-6, 1 - 1e-6))
    c0 = {"probit": float(st.norm.ppf(ybar)), "logit": float(np.log(ybar / (1 - ybar))),
          "cloglog": float(np.log(-np.log(1 - ybar)))}[link]
    b0 = np.zeros(main.k)
    if main.constant:
        b0[-1] = c0
    null = None
    if main.constant:
        eng0 = Engine(lik, [constant_only(main)], w, 0)
        null = maximize(eng0.fun, np.array([c0]), maxiter=300)
    res = pr.maximize(eng.fun, b0)
    _log(s, pr, res, "log pseudolikelihood")
    V, n_clust = pr.vce(res.H, eng.scores(res.b) if pr.robust else None)
    ll0 = null.ll if null else SYS
    r2p = 1 - res.ll / ll0 if null and ll0 else None

    def hdr(est):
        x = est.extra
        left = [est.title, "", "", f"Log pseudolikelihood = {g(x['ll'], '%10.0g'):>10}"]
        return left, std_right(x)

    est = finish_model(s, pr, "fracreg", args, title=_FRAC_LINKS[link], blocks=[main], eng=eng, res=res,
                       V=V, n_clust=n_clust, ll0=ll0, aux_names=[], aux_rows_fn=lambda lv: [], hdr_fn=hdr,
                       wald=True, r2p=r2p, k_eq=1,
                       extra_macros=[("link", link)])
    F = {"probit": st.norm.cdf, "logit": special.expit, "cloglog": lambda e: 1 - np.exp(-np.exp(e))}[link]
    est.predict = make_predict(lambda what, etas, th, e: F(etas[0]) if what in ("cm", "pr", "p") else None, "cm")
    est.pred_label = f"conditional mean of {est.depvar}"
    est.robust_label = True
    display_eq(s, est)


# ---------------------------------------------------------------------------
# zip and zinb
# ---------------------------------------------------------------------------

def _zero_inflated(s: "Session", args: str, cmd: str) -> None:
    nb = cmd == "zinb"
    if _replay(s, cmd, args):
        return
    spec = {"inflate": 3, "probit": 3, "vuong": 3, "irr": 3, "offset": 3, "exposure": 3}
    o = peek(args, spec)
    if o.get("inflate") in (None, True):
        raise StataError(198, "option inflate() required")
    infl_text = str(o["inflate"])
    infl_vars, infl_opts = sub_options(infl_text, {"offset": 3})
    if infl_vars.strip() == "_cons":
        infl_vars = ""
    pr = _Prep(s, args, spec, extra_vars=names_in(infl_vars, o.get("offset"), o.get("exposure"),
                                                  infl_opts.get("offset")))
    pr.write_notes(s.output)
    y = _count_y(pr)
    link = "probit" if o.get("probit") else "logit"
    main = block_from_pr(pr, pr.smp.depname)
    _set_offset(main, pr, o.get("offset"), o.get("exposure"))
    infl = block_from_text(pr, "inflate", infl_vars, constant=True)
    _set_offset(infl, pr, infl_opts.get("offset"))
    for nm in infl.notes:
        s.output.write(f"note: {nm} omitted because of collinearity\n", "text")
    w = pr.w
    n_theta = 1 if nb else 0

    def lik(etas, th):
        if nb:
            ll, d1, d2, d3 = L.zinb(y, etas[0], etas[1], th[0], link)
            return ll, [d1, d2], d3.reshape(-1, 1)
        ll, d1, d2 = L.zip_(y, etas[0], etas[1], link)
        return ll, [d1, d2], None

    eng = Engine(lik, [main, infl], w, n_theta)
    off = main.offset
    pres, ybar = _poisson_start(y, main.X, w, off)
    mu0 = np.exp(np.clip(main.X @ pres.b + (off if off is not None else 0.0), -700, 700))
    # inflation start: the excess of zeros over what the Poisson fit predicts
    excess = float(np.clip(np.sum(w * (y == 0)) / np.sum(w) - np.sum(w * np.exp(-mu0)) / np.sum(w), 0.02, 0.8))
    g0 = np.zeros(infl.k)
    g0[-1] = float(np.log(excess / (1 - excess))) if link == "logit" else float(st.norm.ppf(excess))
    b0 = np.concatenate([pres.b, g0] + ([[np.log(_alpha_start(y, mu0, w))]] if nb else []))
    res = pr.maximize(eng.fun, b0)
    V, n_clust = pr.vce(res.H, eng.scores(res.b) if pr.robust else None)
    null = _null_fit(pr, lik, [main, infl], w, n_theta, res.b)
    if null is not None:
        fit_log(s, pr, "Fitting constant-only model:", null)
    fit_log(s, pr, "Fitting full model:", res)
    ll0 = null.ll if null else SYS
    k_main = main.k
    # comparison models: Poisson (zip) or zip (zinb) for the Vuong and alpha tests
    extra_scalars, chibar = [], None
    if nb:
        lna = float(res.b[-1])
        se_lna = float(np.sqrt(V[-1, -1]))
        # zip with the same equations: alpha = 0 is the limit lnalpha -> -inf
        zeng = Engine(lambda e, th: (lambda ll, a, b: (ll, [a, b], None))(*L.zip_(y, e[0], e[1], link)),
                      [main, infl], w, 0)
        rz = maximize(zeng.fun, np.concatenate([pres.b, g0]), maxiter=300)
        chibar = max(2 * (res.ll - rz.ll), 0.0)
        p_chibar = 0.5 * _chi2p(chibar, 1) if chibar > 0 else 0.5
        extra_scalars = [("alpha", float(np.exp(lna)))] + ([] if pr.robust else [("chi2_c", chibar), ("p_c", p_chibar)])
    vuong = None
    if o.get("vuong"):
        if nb:
            # standard negative binomial
            def nbl(etas, th):
                ll, d1, d2 = L.gnbreg(y, etas[0], np.full(len(y), th[0]))
                return ll, [d1], d2.reshape(-1, 1)
            e1 = Engine(nbl, [main], w, 1)
            r1 = maximize(e1.fun, np.concatenate([pres.b, [float(np.log(_alpha_start(y, mu0, w)))]]), maxiter=300)
            f1 = e1.ll_obs(r1.b)
            name = "zinb vs. standard negative binomial"
        else:
            ep = Engine(lambda e, th: (lambda ll, d: (ll, [d], None))(*_pois_ll(y, e[0])), [main], w, 0)
            r1 = maximize(ep.fun, pres.b, maxiter=300)
            f1 = ep.ll_obs(r1.b)
            name = "zip vs. standard Poisson"
        m = eng.ll_obs(res.b) - f1
        n = len(m)
        zv = float(np.sqrt(n) * m.mean() / m.std()) if m.std() > 0 else SYS
        vuong = (name, zv, float(st.norm.sf(zv)) if zv < SYS else SYS)

    def aux_rows(level):
        if not nb:
            return []
        return [CoefRow("sep"), aux_row("/lnalpha", lna, se_lna, level), CoefRow("sep"),
                transformed_row("alpha", np.exp, np.exp, lna, se_lna, level)]

    n_zero = _wsum(w, y == 0)
    n_nonzero = _wsum(w, y > 0)
    title = "Zero-inflated negative binomial regression" if nb else "Zero-inflated Poisson regression"

    def hdr(est):
        x = est.extra
        ll4 = g(x["ll"], "%9.0g")
        left = [title, "", "", "", f"Inflation model = {link}", f"{lab(est):<16}= {ll4}"]
        right = [("Number of obs", comma(x["N"])), ("Nonzero obs", comma(n_nonzero)), ("Zero obs", comma(n_zero)),
                 ("", ""), (f"{x['chi2type']} chi2({int(x['df_m'])})", g(x["chi2"], "%10.2f")),
                 ("Prob > chi2", g(x["p"], "%10.4f"))]   # VERIFY layout
        return left, right

    def footer(s_, est):
        out = s_.output
        if vuong:
            nm, zv, pz = vuong
            out.write(f"{'Vuong test of ' + nm + ':':<51}z = {g(zv, '%8.2f'):>8}  Pr>z = {g(pz, '%6.4f')}\n", "text")
        if nb and not pr.robust and lna > -10:     # VERIFY: Stata omits the test at the boundary
            out.write(f"Likelihood-ratio test of alpha=0:  chibar2(01) = {chibar:7.2f} Prob>=chibar2 = {p_chibar:5.3f}\n",
                      "text")   # VERIFY

    est = finish_model(s, pr, cmd, args, title=title, blocks=[main, infl], eng=eng, res=res, V=V,
                       n_clust=n_clust, ll0=ll0, aux_names=[("lnalpha", "_cons")] if nb else [],
                       aux_rows_fn=aux_rows, hdr_fn=hdr, footer=footer, k_eq=3 if nb else 2,
                       extra_scalars=extra_scalars,
                       extra_macros=[("inflate", link)], test_idx=list(range(main.nslopes)),
                       wald=pr.robust, predict_opts={"xbinflate": 3})
    est.n_zero = (n_zero, n_nonzero)

    def pred(what, etas, th, e):
        mu = np.exp(etas[0])
        pi = special.expit(etas[1]) if link == "logit" else st.norm.cdf(etas[1])
        if what == "n":
            return (1 - pi) * mu
        if what == "xbinflate":
            return etas[1]
        if what == "pr":
            if nb:
                a = np.exp(th[0])
                p0 = (1 + a * mu) ** (-1 / a)
            else:
                p0 = np.exp(-mu)
            return pi + (1 - pi) * p0
        return None
    est.predict = make_predict(pred, "n")
    est.pred_label = "number of events"
    if o.get("irr"):
        est.coef_title, est.eform_on = "IRR", True
    display_eq(s, est)
    est.coef_title, est.eform_on = "Coef.", False


def _pois_ll(y, eta):
    mu = np.exp(np.clip(eta, -700, 700))
    return -mu + y * eta - special.gammaln(y + 1), y - mu


def cmd_zip(s: "Session", args: str) -> None:
    _zero_inflated(s, args, "zip")


def cmd_zinb(s: "Session", args: str) -> None:
    _zero_inflated(s, args, "zinb")


# ---------------------------------------------------------------------------
# tpoisson and tnbreg
# ---------------------------------------------------------------------------

def _truncation(o, pr) -> float:
    if o.get("ul") not in (None, False):
        raise StataError(198, "option ul() is not available yet in OpenDTA")
    if o.get("ll") in (None, False, True):
        return 0.0
    return float(np.floor(float(pr.s.eval(str(o["ll"])))))


def _truncated(s: "Session", args: str, cmd: str) -> None:
    nb = cmd == "tnbreg"
    if _replay(s, cmd, args):
        return
    spec = {"ll": 2, "ul": 2, "irr": 3, "offset": 3, "exposure": 3, "dispersion": 4}
    o = peek(args, spec)
    pr = _Prep(s, args, spec, extra_vars=names_in(o.get("offset"), o.get("exposure")))
    pr.write_notes(s.output)
    lower = _truncation(o, pr)
    disp = str(o.get("dispersion", "mean")).strip().lower() or "mean"
    if nb and not ("mean".startswith(disp) and disp):
        raise StataError(198, "option dispersion(constant) is not available yet in OpenDTA")
    y = _count_y(pr)
    if np.any(y <= lower):
        raise StataError(459, f"{pr.smp.depname} <= {lower:g} in sample")   # VERIFY
    main = block_from_pr(pr, pr.smp.depname)
    _set_offset(main, pr, o.get("offset"), o.get("exposure"))
    w = pr.w
    n_theta = 1 if nb else 0

    def lik(etas, th):
        if nb:
            ll, d1, d2 = L.tnbreg(y, etas[0], th[0], lower)
            return ll, [d1], d2.reshape(-1, 1)
        ll, d1 = L.tpoisson(y, etas[0], lower)
        return ll, [d1], None

    eng = Engine(lik, [main], w, n_theta)
    pres, _ = _poisson_start(y, main.X, w, main.offset)
    mu0 = np.exp(np.clip(main.X @ pres.b + (main.offset if main.offset is not None else 0.0), -700, 700))
    b0 = np.concatenate([pres.b] + ([[np.log(_alpha_start(y, mu0, w))]] if nb else []))
    res = pr.maximize(eng.fun, b0)
    V, n_clust = pr.vce(res.H, eng.scores(res.b) if pr.robust else None)
    # constant-only model: start from the log of the mean
    null = None
    if main.constant:
        eng0 = Engine(lik, [constant_only(main)], w, n_theta)
        c0 = float(np.log(max(np.sum(w * y) / np.sum(w) - lower, 0.5)))
        null = maximize(eng0.fun, np.array([c0] + ([float(res.b[-1])] if nb else [])), maxiter=300)
    ll0 = null.ll if null else SYS
    ptr = None
    if nb:
        ep = Engine(lambda e, th: (lambda ll, d: (ll, [d], None))(*L.tpoisson(y, e[0], lower)), [main], w, 0)
        ptr = maximize(ep.fun, pres.b, maxiter=300)
        fit_log(s, pr, "Fitting truncated Poisson model:", ptr)
        if null is not None:
            fit_log(s, pr, "Fitting constant-only model:", null)
        fit_log(s, pr, "Fitting full model:", res)
    else:
        _log(s, pr, res)
    if nb:
        lna, se_lna = float(res.b[-1]), float(np.sqrt(V[-1, -1]))

    def aux_rows(level):
        if not nb:
            return []
        return [CoefRow("sep"), aux_row("/lnalpha", lna, se_lna, level), CoefRow("sep"),
                transformed_row("alpha", np.exp, np.exp, lna, se_lna, level)]

    title = "Truncated negative binomial regression" if nb else "Truncated Poisson regression"

    def hdr(est):
        x = est.extra
        left = [title, f"Truncation point: {lower:g}", "Dispersion     = mean" if nb else "",
                f"{lab(est)} = {g(x['ll'], '%10.0g'):>10}"]
        return left, std_right(x)

    est = finish_model(s, pr, cmd, args, title=title, blocks=[main], eng=eng, res=res, V=V, n_clust=n_clust,
                       ll0=ll0, aux_names=[("lnalpha", "_cons")] if nb else [], aux_rows_fn=aux_rows,
                       hdr_fn=hdr, k_eq=2 if nb else 1, wald=pr.robust,
                       extra_scalars=[("llopt", lower)] + ([("alpha", float(np.exp(lna)))] if nb else []))
    if nb and not pr.robust:
        chibar = max(2 * (res.ll - ptr.ll), 0.0)
        p_chibar = 0.5 * _chi2p(chibar, 1) if chibar > 0 else 0.5
        est.footer = lambda s_, e: s_.output.write(
            f"Likelihood-ratio test of alpha=0:  chibar2(01) = {chibar:7.2f} Prob>=chibar2 = {p_chibar:5.3f}\n", "text")
    mean_fn = _trunc_mean(nb, lower)
    est.predict = make_predict(lambda what, etas, th, e: mean_fn(np.exp(etas[0]), th) if what == "n" else None, "n")
    est.pred_label = "number of events"
    if o.get("irr"):
        est.coef_title, est.eform_on = "IRR", True
    display_eq(s, est)
    est.coef_title, est.eform_on = "Coef.", False


def _trunc_mean(nb: bool, lower: float):
    def mean(mu, th):
        # E[y | y > lower] = (mu - sum_{j <= lower} j p_j) / (1 - sum_{j <= lower} p_j)
        num = mu.copy()
        den = np.ones_like(mu)
        for j in range(int(lower) + 1):
            if nb:
                m = np.exp(-th[0])
                pj = np.exp(L._nb_terms(np.full_like(mu, j), mu, m)[0])
            else:
                pj = st.poisson.pmf(j, mu)
            num -= j * pj
            den -= pj
        return num / den
    return mean


def cmd_tpoisson(s: "Session", args: str) -> None:
    _truncated(s, args, "tpoisson")


def cmd_tnbreg(s: "Session", args: str) -> None:
    _truncated(s, args, "tnbreg")


# ---------------------------------------------------------------------------
# gnbreg
# ---------------------------------------------------------------------------

def cmd_gnbreg(s: "Session", args: str) -> None:
    if _replay(s, "gnbreg", args):
        return
    spec = {"lnalpha": 3, "irr": 3, "offset": 3, "exposure": 3}
    o = peek(args, spec)
    la_vars, la_opts = sub_options(str(o["lnalpha"]), {"offset": 3}) if o.get("lnalpha") not in (None, True) \
        else ("", {})
    pr = _Prep(s, args, spec, extra_vars=names_in(la_vars, o.get("offset"), o.get("exposure"),
                                                  la_opts.get("offset")))
    pr.write_notes(s.output)
    y = _count_y(pr)
    main = block_from_pr(pr, pr.smp.depname)
    _set_offset(main, pr, o.get("offset"), o.get("exposure"))
    la = block_from_text(pr, "lnalpha", la_vars, constant=True)
    _set_offset(la, pr, la_opts.get("offset"))
    for nm in la.notes:
        s.output.write(f"note: {nm} omitted because of collinearity\n", "text")
    w = pr.w

    def lik(etas, th):
        ll, d1, d2 = L.gnbreg(y, etas[0], etas[1])
        return ll, [d1, d2], None

    eng = Engine(lik, [main, la], w, 0)
    pres, _ = _poisson_start(y, main.X, w, main.offset)
    mu0 = np.exp(np.clip(main.X @ pres.b + (main.offset if main.offset is not None else 0.0), -700, 700))
    g0 = np.zeros(la.k)
    g0[-1] = np.log(_alpha_start(y, mu0, w))
    res = pr.maximize(eng.fun, np.concatenate([pres.b, g0]))
    V, n_clust = pr.vce(res.H, eng.scores(res.b) if pr.robust else None)
    # constant-only model: both equations with the constant only
    null = None
    if main.constant:
        e0 = Engine(lik, [constant_only(main), constant_only(la)], w, 0)
        null = maximize(e0.fun, np.array([res.b[main.k - 1], res.b[-1]]), maxiter=300)
    if null is not None:
        fit_log(s, pr, "Fitting constant-only model:", null)
    fit_log(s, pr, "Fitting full model:", res)
    ll0 = null.ll if null else SYS
    est = finish_model(s, pr, "gnbreg", args, title="Generalized negative binomial regression", blocks=[main, la],
                       eng=eng, res=res, V=V, n_clust=n_clust, ll0=ll0, aux_names=[],
                       aux_rows_fn=lambda lv: [], hdr_fn=_hdr_logit_style(), k_eq=2, wald=pr.robust)
    est.predict = make_predict(lambda what, etas, th, e: np.exp(etas[0]) if what == "n" else None, "n")
    est.pred_label = "number of events"
    if o.get("irr"):
        est.coef_title, est.eform_on = "IRR", True
    display_eq(s, est)
    est.coef_title, est.eform_on = "Coef.", False


# ---------------------------------------------------------------------------
# binreg
# ---------------------------------------------------------------------------

_BINREG_LINKS = {"or": ("logit", "Odds Ratio", True), "rr": ("log", "Risk Ratio", True),
                 "hr": ("cloglog", "Haz. Ratio", True), "rd": ("identity", "Risk Diff.", False)}


def cmd_binreg(s: "Session", args: str) -> None:
    """binreg is glm with the binomial family: the option picks the link (logit, log, cloglog,
    identity) and the title of the coefficients."""
    from ..lang.syntax import find_top, parse_options
    from .glm import cmd_glm
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, "glm", args)
        return
    t = args.strip()
    c = find_top(t, ",")
    head, tail = (t[:c], t[c + 1:]) if c != -1 else (t, "")
    link, title, nbin = None, ("Coef.", False), None
    keep = []
    for name, arg in parse_options(tail):
        low = name.lower()
        hit = [k for k in _BINREG_LINKS if k.startswith(low) and len(low) >= 2 or k == low]
        if hit:
            if link is not None:
                raise StataError(198, "only one of or, rr, hr and rd may be specified")   # VERIFY
            link = _BINREG_LINKS[hit[0]][0]
            title = (_BINREG_LINKS[hit[0]][1], _BINREG_LINKS[hit[0]][2])
        elif low == "n" and arg is not None:
            nbin = arg.strip()
        elif low in ("ml", "irls") or (low.startswith("mu") and arg is not None):
            if low.startswith("mu"):
                raise StataError(198, "option mu() is not available yet in OpenDTA")
            keep.append(low)
        else:
            keep.append(f"{name}({arg})" if arg is not None else name)
    link = link or "logit"
    if not {"ml", "irls"} & set(keep):
        keep.append("irls")        # binreg uses IRLS unless ml is given
    opts = ["family(binomial" + (f" {nbin}" if nbin else "") + ")", f"link({link})"] + keep
    cmd_glm(s, f"{head} , " + " ".join(opts), titles=title)
    if "cmdline" in s.e:
        s.e["cmdline"] = f"binreg {args.strip()}"


# ---------------------------------------------------------------------------
# glogit
# ---------------------------------------------------------------------------

def cmd_glogit(s: "Session", args: str) -> None:
    """Weighted least squares logit on grouped data: regress ln((pos + .5) / (pop - pos + .5)) on
    the covariates with weights (pos + .5)(pop - pos + .5)/(pop + 1) ([R] glogit)."""
    from ..commands._util import touse
    from ..core.dataset import Variable
    from ..lang.syntax import find_top, match_options, parse_standard
    p = parse_standard(args)
    o = match_options(p.options, {"or": 2, "level": 1}) if p.options.strip() else {}
    if o.get("or"):
        raise StataError(198, "option or is not available yet in OpenDTA")
    toks = p.varlist.split()
    if len(toks) < 2:
        raise StataError(102, "too few variables specified")
    pos, pop = toks[0], toks[1]
    ds = s.data
    n = ds.nobs
    a, b = ds.get(pos).data, ds.get(pop).data
    ok = touse(s, p) & (a < SYS) & (b < SYS)
    if np.any(ok & ((a < 0) | (a > b))):
        raise StataError(459, f"{pos} must be between 0 and {pop}")   # VERIFY
    yv = np.full(n, SYS)
    wv = np.full(n, SYS)
    yv[ok] = np.log((a[ok] + .5) / (b[ok] - a[ok] + .5))
    wv[ok] = (a[ok] + .5) * (b[ok] - a[ok] + .5) / (b[ok] + 1)
    t_pos, t_pop, t_w = "_odta_gl_pos", "_odta_gl_pop", "_odta_gl_w"
    for nm in (t_pos, t_pop, t_w):
        if ds.has(nm):
            raise StataError(110, f"{nm} already defined")
    ds.rename(pos, t_pos)
    ds.rename(pop, t_pop)
    try:
        for nm, vals in ((pos, yv), (t_w, wv)):
            v = Variable(nm, "double", np.full(n, SYS))
            v.data = vals
            ds.add(v)
        c = find_top(args, ",")
        head = args[:c] if c != -1 else args
        rest = " ".join(toks[2:])
        lvl = f", level({o['level']})" if o.get("level") else ""
        # the sample is the observations where the logit is defined; if/in were applied above
        s.run(f"quietly regress {pos} {rest} [aw={t_w}]{lvl}")
        if "cmdline" in s.e:
            s.e["cmdline"] = f"glogit {args.strip()}"
        s.output.write("\nWeighted LS logistic regression for grouped data\n", "text")
        est = s._estimates
        est.depvar_cell = ""
        s.run("regress")
        del est.depvar_cell
    finally:
        ds.drop_vars([pos, t_w])
        ds.rename(t_pos, pos)
        ds.rename(t_pop, pop)
        s.notify_state()
