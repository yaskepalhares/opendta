"""Treatment-effects models: etregress ([R] etregress). The log likelihood is in likelihoods.py."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..lazy import lazy
from . import likelihoods as L
from .censmodels import _chi2p, _replay, _rho_rows, _weighted_probit, _lambda_row
from .mleq import (Engine, aux_row, block_from_pr, block_from_text, display_eq, finish_model, lab,
                   make_predict, names_in, peek, transformed_row)
from .models2 import _Prep, _log, g
from .regress import comma
from .results import CoefRow

st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def _treat_spec(args: str):
    """(rewritten args with the treatment variable at the end of the outcome varlist, treatment
    variable, its varlist, noconstant flag)."""
    o = peek(args, {"treat": 3, "twostep": 3})
    txt = o.get("treat")
    if txt in (None, True) or "=" not in str(txt):
        raise StataError(198, "option treat() required")
    tname, rest = (t.strip() for t in str(txt).split("=", 1))
    zvars, noc = rest, False
    if "," in rest:
        zvars, tail = (t.strip() for t in rest.split(",", 1))
        noc = "noconstant".startswith(tail.split()[0]) if tail.split() else False
    head, sep, opts = args.partition(",")
    return f"{head.strip()} {tname}{sep}{opts}", tname, zvars, noc, o


def _setup(s, args, spec):
    args2, tname, zvars, noc, o = _treat_spec(args)
    pr = _Prep(s, args2, spec, extra_vars=names_in(zvars, tname))
    pr.write_notes(s.output)
    t = (s.data.get(tname).data[pr.smp.mask] != 0).astype(float)
    main = block_from_pr(pr, pr.smp.depname)
    trt = block_from_text(pr, tname, zvars, constant=not noc)
    for nm in trt.notes:
        s.output.write(f"note: {nm} omitted because of collinearity\n", "text")
    return pr, pr.y, t, main, trt, tname, o


def _label_trt(est, tname: str) -> None:
    """The maximum-likelihood tables name the treatment indicator in the outcome equation 1.t."""
    for r in est.rows:
        if r.kind == "coef" and r.label == tname:
            r.label = f"1.{tname}"
            break


def cmd_etregress(s: "Session", args: str) -> None:
    if _replay(s, "etregress", args):
        return
    pr, y, t, main, trt, tname, o = _setup(s, args, {"treat": 3, "twostep": 3})
    w = pr.w
    n = len(y)
    if o.get("twostep"):
        return _twostep(s, args, pr, y, t, main, trt, tname)

    def lik(etas, th):
        ll, d_eta, d_zeta, d_ath, d_ls = L.etreg(t, y, etas[0], etas[1], th[0], th[1])
        return ll, [d_eta, d_zeta], np.column_stack([d_ath, d_ls])

    eng = Engine(lik, [main, trt], w, 2)
    rp = _weighted_probit(t, trt.X, w)
    sw = np.sqrt(w)
    bo = np.linalg.lstsq(main.X * sw[:, None], y * sw, rcond=None)[0]
    rss = float(np.sum(w * (y - main.X @ bo) ** 2))
    n1 = float(w.sum())
    b0 = np.concatenate([bo, rp.b, [0.0, 0.5 * np.log(rss / n1)]])
    res = pr.maximize(eng.fun, b0)
    _log(s, pr, res)
    V, n_clust = pr.vce(res.H, eng.scores(res.b) if pr.robust else None)
    ath, lns = float(res.b[-2]), float(res.b[-1])
    V2 = V[-2:, -2:]
    ll_ols = -0.5 * n1 * (np.log(2 * np.pi) + 1 + np.log(rss / n1))
    lr = max(2 * (res.ll - (rp.ll + ll_ols)), 0.0)
    p_lr = _chi2p(lr, 1)
    w_rho = ath ** 2 / float(V[-2, -2])
    p_w = _chi2p(w_rho, 1)

    def aux_rows(level):
        return [CoefRow("sep"), aux_row("/athrho", ath, float(np.sqrt(V[-2, -2])), level, name="_cons"),
                aux_row("/lnsigma", lns, float(np.sqrt(V[-1, -1])), level, name="_cons"), CoefRow("sep"),
                transformed_row("rho", np.tanh, lambda x: 1 - np.tanh(x) ** 2, ath, float(np.sqrt(V[-2, -2])), level),
                transformed_row("sigma", np.exp, np.exp, lns, float(np.sqrt(V[-1, -1])), level),
                _lambda_row(level, ath, lns, V2)]

    def hdr(est):
        x = est.extra
        left_ = [est.title, "Estimator: maximum likelihood", f"{lab(est)} = {g(x['ll'], '%10.0g').strip():>9}"]
        right = [("Number of obs", comma(x["N"])),
                 (f"{x['chi2type']} chi2({int(x['df_m'])})", g(x["chi2"], "%10.2f")),
                 ("Prob > chi2", g(x["p"], "%10.4f"))]
        return left_, right

    def footer(s_, est):
        if pr.robust:
            s_.output.write(f"Wald test of indep. eqns. (rho = 0): chi2(1) = {w_rho:8.2f}   Prob > chi2 = {p_w:6.4f}\n", "text")
        else:
            s_.output.write(f"LR test of indep. eqns. (rho = 0):   chi2(1) = {lr:8.2f}   Prob > chi2 = {p_lr:6.4f}\n", "text")

    rho, sg = float(np.tanh(ath)), float(np.exp(lns))
    est = finish_model(s, pr, "etregress", args, title="Linear regression with endogenous treatment",
                       blocks=[main, trt], eng=eng, res=res, V=V, n_clust=n_clust, ll0=SYS,
                       aux_names=[("athrho", "_cons"), ("lnsigma", "_cons")], aux_rows_fn=aux_rows, hdr_fn=hdr,
                       depvar_cell="", footer=footer, wald=True, k_eq=4, test_idx=list(range(main.nslopes)),
                       extra_scalars=[("rho", rho), ("sigma", sg), ("lambda", rho * sg), ("chi2_c", lr),
                                      ("p_c", p_lr)],
                       predict_opts={"xbtrt": 4, "ptrt": 4})
    est.predict = make_predict(_predict, "xb")
    est.pred_label = "fitted values"
    _label_trt(est, tname)
    display_eq(s, est)


def _predict(what, etas, th, est):
    if what == "xb":
        return etas[0]
    if what == "xbtrt":
        return etas[1]
    if what in ("ptrt", "pr"):
        return st.norm.cdf(etas[1])
    return None


def _twostep(s, args, pr, y, t, main, trt, tname):
    """Two-step estimates: probit for the treatment, then OLS of the outcome on the covariates and
    the hazard (Maddala 1983; Greene). VERIFY: variance and sigma formulas."""
    from .ml import MLResult
    if pr.smp.wtype:
        raise StataError(101, "weights not allowed")   # VERIFY
    n = len(y)
    ps = _weighted_probit(t, trt.X, np.ones(n))
    Vg = np.linalg.inv(-ps.H)
    zg = trt.X @ ps.b
    q = 2 * t - 1
    lam = q * np.exp(st.norm.logpdf(zg) - st.norm.logcdf(q * zg))
    delta = lam * (lam + zg)
    Xs = np.column_stack([main.X, lam])
    coef = np.linalg.lstsq(Xs, y, rcond=None)[0]
    beta, theta = coef[:-1], coef[-1]
    e = y - Xs @ coef
    sig2 = float(e @ e / n + theta ** 2 * delta.mean())
    sig = float(np.sqrt(sig2))
    rho = float(np.clip(theta / sig, -1.0, 1.0))
    XtXi = np.linalg.inv(Xs.T @ Xs)
    Fm = Xs.T @ (trt.X * delta[:, None])
    mid = Xs.T @ (Xs * (1 - rho ** 2 * delta)[:, None]) + rho ** 2 * Fm @ Vg @ Fm.T
    Vs = sig2 * XtXi @ mid @ XtXi
    k1, k2 = main.k, trt.k
    b_all = np.concatenate([beta, ps.b, [theta]])
    V_all = np.zeros((k1 + k2 + 1, k1 + k2 + 1))
    pos = list(range(k1)) + [k1 + k2]
    V_all[np.ix_(pos, pos)] = Vs
    V_all[k1:k1 + k2, k1:k1 + k2] = Vg
    res = MLResult(b_all, SYS, np.zeros(len(b_all)), np.zeros((len(b_all), len(b_all))), 0, True, [])
    se_l = float(np.sqrt(Vs[-1, -1]))

    def aux_rows(level):
        return [CoefRow("sep"), CoefRow("eq", "hazard"), CoefRow("coef", "lambda", theta, se_l, name="lambda"),
                CoefRow("sep"), CoefRow("value", "rho", rho, eq=f"{rho:.5f}"),
                CoefRow("value", "sigma", sig, eq=g(sig, "%10.8g"))]

    def hdr(est):
        x = est.extra
        def write(s_):
            o_ = s_.output
            o_.write("\n", "text")
            for lft, lb, val in (("Linear regression with endogenous treatment", "Number of obs", comma(x["N"])),
                                 ("Estimator: two-step", f"Wald chi2({int(x['df_m'])})", g(x["chi2"], "%10.2f")),
                                 ("", "Prob > chi2", g(x["p"], "%10.4f"))):
                o_.write(f"{lft:<48}{lb:<19}=", "text")
                o_.write(f"{val:>10}\n", "result")
            o_.write("\n", "text")
        return write

    eng = Engine(lambda e_, th: (np.zeros(n), [np.zeros(n)], None), [main], np.ones(n), 0)
    order = (["N", "df_m", "chi2", "p", "rho", "sigma", "lambda"],
             ["cmdline", "cmd", "title", "chi2type", "vce", "depvar", "properties"], [])
    est = finish_model(s, pr, "etregress", args, title="Linear regression with endogenous treatment",
                       blocks=[main, trt], eng=eng, res=res, V=V_all, n_clust=None, ll0=SYS,
                       aux_names=[("hazard", "lambda")], aux_rows_fn=aux_rows, hdr_fn=hdr, depvar_cell="",
                       wald=True, k_eq=3, test_idx=list(range(main.nslopes)) + [main.k + trt.k], order=order,
                       extra_scalars=[("rho", rho), ("sigma", sig), ("lambda", theta)],
                       predict_opts={"xbtrt": 4, "ptrt": 4})
    est.predict = make_predict(_predict, "xb")
    est.pred_label = "fitted values"
    display_eq(s, est)


def cmd_etpoisson(s: "Session", args: str) -> None:
    if _replay(s, "etpoisson", args):
        return
    pr, y, t, main, trt, tname, o = _setup(s, args, {"treat": 3, "intpoints": 3})
    w = pr.w
    K = int(o["intpoints"]) if o.get("intpoints") else 24
    if np.any(y < 0) or np.any(y != np.round(y)):
        raise StataError(459, f"{pr.smp.depname} must be a nonnegative integer")   # VERIFY wording

    def lik(etas, th):
        ll, d_eta, d_zeta, d_ath, d_ls = L.etpois(t, y, etas[0], etas[1], th[0], th[1], K)
        return ll, [d_eta, d_zeta], np.column_stack([d_ath, d_ls])

    eng = Engine(lik, [main, trt], w, 2)
    from ..core.treatfx import glm_fit
    rp = _weighted_probit(t, trt.X, w)
    bo = glm_fit("poisson", main.X, y, w)
    b0 = np.concatenate([bo, rp.b, [0.0, np.log(0.5)]])
    res = pr.maximize(eng.fun, b0)
    _log(s, pr, res)
    V, n_clust = pr.vce(res.H, eng.scores(res.b) if pr.robust else None)
    ath, lns = float(res.b[-2]), float(res.b[-1])
    se_a, se_l = float(np.sqrt(V[-2, -2])), float(np.sqrt(V[-1, -1]))
    w_rho = ath ** 2 / V[-2, -2]
    p_w = _chi2p(w_rho, 1)

    def aux_rows(level):
        return [CoefRow("sep"), aux_row("/athrho", ath, se_a, level, name="_cons"),
                aux_row("/lnsigma", lns, se_l, level, name="_cons"), CoefRow("sep"),
                transformed_row("rho", np.tanh, lambda x: 1 - np.tanh(x) ** 2, ath, se_a, level),
                transformed_row("sigma", np.exp, np.exp, lns, se_l, level)]

    def hdr(est):
        x = est.extra
        left_ = [est.title, f"({K} quadrature points)", f"{lab(est)} = {g(x['ll'], '%10.0g').strip():>9}"]
        right = [("Number of obs", comma(x["N"])),
                 (f"{x['chi2type']} chi2({int(x['df_m'])})", g(x["chi2"], "%10.2f")),
                 ("Prob > chi2", g(x["p"], "%10.4f"))]
        return left_, right

    def footer(s_, est):
        s_.output.write(f"Wald test of indep. eqns. (rho = 0): chi2(1) = {w_rho:8.2f}   Prob > chi2 = {p_w:6.4f}\n", "text")

    est = finish_model(s, pr, "etpoisson", args, title="Poisson regression with endogenous treatment",
                       blocks=[main, trt], eng=eng, res=res, V=V, n_clust=n_clust, ll0=SYS,
                       aux_names=[("athrho", "_cons"), ("lnsigma", "_cons")], aux_rows_fn=aux_rows, hdr_fn=hdr,
                       depvar_cell="", footer=footer, wald=True, k_eq=4, test_idx=list(range(main.nslopes)),
                       extra_scalars=[("rho", float(np.tanh(ath))), ("sigma", float(np.exp(lns))),
                                      ("chi2_c", w_rho), ("p_c", p_w)],
                       predict_opts={"xbtrt": 4, "ptrt": 4})
    est.predict = make_predict(_predict, "xb")
    est.pred_label = "linear prediction"
    _label_trt(est, tname)
    display_eq(s, est)
