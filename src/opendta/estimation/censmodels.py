"""Interval, truncated and selection models: intreg, truncreg, heckman, heckprobit,
heckoprobit and biprobit ([R] intreg, [R] truncreg, [R] heckman, [R] heckprobit,
[R] heckoprobit, [R] biprobit). The log likelihoods are in likelihoods.py."""

from __future__ import annotations

import dataclasses
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.stats import chi2_tail
from ..lazy import lazy
from . import fvars
from . import likelihoods as L
from .ml import maximize
from .mleq import (Block, Engine, aux_row, block_from_pr, block_from_text, constant_only, display_eq,
                   finish_model, fit_log, lab, make_predict, names_in, peek, std_right, sub_options,
                   transformed_row)
from .models import _probit
from .models2 import _Prep, _log, g
from .regress import comma
from .results import CoefRow, crit

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


def _chi2p(x: float, df: int) -> float:
    return float(chi2_tail(x, df)) if df > 0 else SYS


def _split_vars(args: str, k: int) -> tuple[list[str], str]:
    """The first k plain variable names of a command line and the rest of the line."""
    text = args.strip()
    names = []
    for _ in range(k):
        parts = text.split(None, 1)
        if not parts or parts[0].startswith(",") or "(" in parts[0]:
            raise StataError(102, "too few variables specified")
        names.append(parts[0])
        text = parts[1] if len(parts) > 1 else ""
    return names, text


def _values(ds, name: str, mask) -> np.ndarray:
    """Variable values on the sample with missing as nan."""
    v = ds.get(name).data[mask].astype(float)
    return np.where(v >= SYS, np.nan, v)


def _weighted_probit(y, X, w, off=None):
    fun, _ = _probit(y, X, w)
    return maximize(fun, np.zeros(X.shape[1]), maxiter=300)


# ---------------------------------------------------------------------------
# intreg
# ---------------------------------------------------------------------------

def cmd_intreg(s: "Session", args: str) -> None:
    if _replay(s, "intreg", args):
        return
    (d1, d2), rest = _split_vars(args, 2)
    ds = s.data
    spec = {"offset": 3, "het": 3}
    o = peek(rest, spec)
    het_vars = str(o["het"]) if o.get("het") not in (None, True) else None
    both_missing = (ds.get(d1).data >= SYS) & (ds.get(d2).data >= SYS)
    pr = _Prep(s, rest, spec, need_dep=False, restrict=~both_missing, extra_vars=names_in(het_vars))
    pr.smp.depname = d1
    pr.write_notes(s.output)
    lo = _values(ds, d1, pr.smp.mask)
    hi = _values(ds, d2, pr.smp.mask)
    if np.any(lo > hi):
        raise StataError(498, f"{d1} > {d2} in {int(np.sum(lo > hi))} observations")   # VERIFY
    main = block_from_pr(pr, "model" if het_vars else d1)
    hetb = block_from_text(pr, "lnsigma", het_vars, constant=True) if het_vars else None
    w = pr.w
    exact = lo == hi
    left = np.isnan(lo)
    right = np.isnan(hi)
    inter = ~exact & ~left & ~right

    if hetb is None:
        def lik(etas, th):
            ll, d1_, d2_ = L.intreg(lo, hi, etas[0], th[0])
            return ll, [d1_], d2_.reshape(-1, 1)
        eng = Engine(lik, [main], w, 1)
    else:
        def lik(etas, th):
            ll, d1_, d2_ = L.intreg(lo, hi, etas[0], etas[1])
            return ll, [d1_, d2_], None
        eng = Engine(lik, [main, hetb], w, 0)
    # start: OLS on the exact value, the interval midpoint or the censoring point
    yy = np.where(exact, lo, np.where(inter, (lo + hi) / 2, np.where(left, hi, lo)))
    sw = np.sqrt(w)
    beta = np.linalg.lstsq(main.X * sw[:, None], yy * sw, rcond=None)[0]
    sig = float(np.sqrt(np.sum(w * (yy - main.X @ beta) ** 2) / np.sum(w)))
    if hetb is None:
        b0 = np.concatenate([beta, [np.log(max(sig, 1e-6))]])
    else:
        h0 = np.zeros(hetb.k)
        h0[-1] = np.log(max(sig, 1e-6))
        b0 = np.concatenate([beta, h0])
    res = pr.maximize(eng.fun, b0)
    V, n_clust = pr.vce(res.H, eng.scores(res.b) if pr.robust else None)
    null = None
    if main.constant and hetb is None:
        e0 = Engine(lik, [constant_only(main)], w, 1)
        mu0 = float(np.sum(w * yy) / np.sum(w))
        null = maximize(e0.fun, np.array([mu0, np.log(max(float(np.std(yy)), 1e-6))]), maxiter=300)
        fit_log(s, pr, "Fitting constant-only model:", null)
    fit_log(s, pr, "Fitting full model:", res)
    # observation counts: weighted only with fweights (case 1909: svy: intreg)
    counts = [float(w[m].sum()) if pr.smp.wtype == "fweight" else float(np.sum(m)) for m in (left, exact, right, inter)]
    if hetb is None:
        lns, se_lns = float(res.b[-1]), float(np.sqrt(V[-1, -1]))

        def aux_rows(level):
            return [CoefRow("sep"), aux_row("/lnsigma", lns, se_lns, level), CoefRow("sep"),
                    transformed_row("sigma", np.exp, np.exp, lns, se_lns, level)]
    else:
        def aux_rows(level):
            return []

    def hdr(est):
        x = est.extra
        left_ = [est.title, "", f"{lab(est)} = {g(x['ll'], '%10.0g'):>10}"]
        return left_, std_right({**x, "r2_p": None})

    def footer(s_, est):
        out = s_.output
        for n_, lbl in zip(counts, ("left-censored", "uncensored", "right-censored", "interval")):
            out.write(f"{int(n_):>14}{lbl:>15} observations\n", "text")

    est = finish_model(s, pr, "intreg", args, title="Interval regression", blocks=[main] + ([hetb] if hetb else []),
                       eng=eng, res=res, V=V, n_clust=n_clust, ll0=null.ll if null else SYS,
                       aux_names=[] if hetb else [("lnsigma", "_cons")],
                       aux_rows_fn=aux_rows, hdr_fn=hdr, depvar_cell="", footer=footer,
                       wald=pr.robust or hetb is not None, k_eq=2, r2p=None,
                       extra_scalars=([] if hetb else [("sigma", float(np.exp(lns)))])
                       + [("N_lc", counts[0]), ("N_unc", counts[1]), ("N_rc", counts[2]), ("N_int", counts[3])],
                       extra_macros=[("depvar", f"{d1} {d2}")])
    est.depvar = d1
    est.predict = make_predict(lambda what, etas, th, e: etas[0] if what == "xb" else None, "xb")
    est.pred_label = "fitted values"
    display_eq(s, est)


# ---------------------------------------------------------------------------
# truncreg
# ---------------------------------------------------------------------------

def _limit(s, text) -> float:
    if text in (None, False, True):
        return None
    return float(s.eval(str(text)))


def cmd_truncreg(s: "Session", args: str) -> None:
    if _replay(s, "truncreg", args):
        return
    spec = {"ll": 2, "ul": 2, "offset": 3}
    o = peek(args, spec)
    lower = _limit(s, o.get("ll"))
    upper = _limit(s, o.get("ul"))
    if lower is None and upper is None:
        raise StataError(198, "must specify at least one of ll() and ul()")   # VERIFY
    lo_ = -np.inf if lower is None else lower
    hi_ = np.inf if upper is None else upper
    pr = _Prep(s, args, spec, extra_vars=names_in(o.get("offset")))
    out_ = (pr.y <= lo_) | (pr.y >= hi_)
    n_trunc = int(out_.sum())
    if n_trunc:        # observations outside the limits are dropped
        yfull = s.data.get(args.split()[0]).data
        pr = _Prep(s, args, spec, extra_vars=names_in(o.get("offset")),
                   restrict=~((yfull < SYS) & ((yfull <= lo_) | (yfull >= hi_))))
    pr.write_notes(s.output)
    s.output.write(f"(note: {n_trunc} obs. truncated)\n", "text")
    y = pr.y
    main = block_from_pr(pr, pr.smp.depname)
    w = pr.w

    def lik(etas, th):
        ll, d1, d2 = L.truncreg(y, etas[0], th[0], lo_, hi_)
        return ll, [d1], d2.reshape(-1, 1)

    eng = Engine(lik, [main], w, 1)
    sw = np.sqrt(w)
    beta = np.linalg.lstsq(main.X * sw[:, None], y * sw, rcond=None)[0]
    sig = float(np.sqrt(np.sum(w * (y - main.X @ beta) ** 2) / np.sum(w)))
    res = pr.maximize(eng.fun, np.concatenate([beta, [np.log(max(sig, 1e-6))]]))
    fit_log(s, pr, "Fitting full model:", res)
    Vl, n_clust = pr.vce(res.H, eng.scores(res.b) if pr.robust else None)
    # Stata's e(b) holds sigma, not ln sigma: delta method
    sg = float(np.exp(res.b[-1]))
    J = np.eye(len(res.b))
    J[-1, -1] = sg
    V = J @ Vl @ J.T
    b_rep = np.concatenate([res.b[:-1], [sg]])
    res2 = dataclasses.replace(res, b=b_rep)
    null = None
    if main.constant:
        e0 = Engine(lik, [constant_only(main)], w, 1)
        null = maximize(e0.fun, np.array([float(np.average(y, weights=w)), np.log(max(float(np.std(y)), 1e-6))]),
                        maxiter=300)
    se_sg = float(np.sqrt(V[-1, -1]))

    def aux_rows(level):
        return [CoefRow("sep"), aux_row("/sigma", sg, se_sg, level)]

    def hdr(est):
        x = est.extra
        lo_t = "-inf" if lower is None else f"{lower:g}"
        hi_t = "+inf" if upper is None else f"{upper:g}"
        left_ = [est.title, f"Limit:   lower = {lo_t:>10}", f"         upper = {hi_t:>10}",
                 f"{lab(est)} = {g(x['ll'], '%10.0g'):>10}"]
        right = [("", ""), ("Number of obs", comma(x["N"])),
                 (f"Wald chi2({int(x['df_m'])})", g(x["chi2"], "%10.2f")), ("Prob > chi2", g(x["p"], "%10.4f"))]
        return left_, right   # VERIFY layout

    est = finish_model(s, pr, "truncreg", args, title="Truncated regression", blocks=[main], eng=eng,
                       res=res2, V=V, n_clust=n_clust, ll0=null.ll if null else SYS,
                       aux_names=[("sigma", "_cons")], aux_rows_fn=aux_rows, hdr_fn=hdr, wald=True, k_eq=2,
                       heads=False,
                       extra_scalars=[("sigma", sg)] + ([("llopt", lower)] if lower is not None else [])
                       + ([("ulopt", upper)] if upper is not None else []))
    est.predict = make_predict(lambda what, etas, th, e: etas[0] if what == "xb" else None, "xb")
    est.pred_label = "fitted values"
    display_eq(s, est)


# ---------------------------------------------------------------------------
# selection models: heckman, heckprobit, heckoprobit
# ---------------------------------------------------------------------------

def _selection_setup(s, args: str, spec: dict, outcome_constant: bool | None = None):
    """(dep name, selection variable name or None, selection varlist, rest of the options,
    _Prep) of a selection command."""
    (dep,), rest = _split_vars(args, 1)
    o = peek(rest, {**spec, "select": 3})
    if o.get("select") in (None, True):
        raise StataError(198, "option select() required")
    sel_text = str(o["select"])
    if "=" in sel_text:
        sel_name, sel_vars = (t.strip() for t in sel_text.split("=", 1))
    else:
        sel_name, sel_vars = None, sel_text.strip()
    ds = s.data
    yv = ds.get(dep).data
    n = ds.nobs
    if sel_name:
        sv = ds.get(sel_name).data
        sel_miss = sv >= SYS
        selected = (sv != 0) & ~sel_miss
        restrict = ~sel_miss & (~selected | (yv < SYS))
    else:
        selected = yv < SYS
        restrict = np.ones(n, bool)
    pr = _Prep(s, rest, spec_with_select(spec), need_dep=False, extra_vars=names_in(sel_vars, sel_name),
               restrict=restrict, constant=outcome_constant)
    pr.smp.depname = dep
    return dep, sel_name, sel_vars, pr, selected[pr.smp.mask], o


def spec_with_select(spec: dict) -> dict:
    return {**spec, "select": 3}


def _selection_common(s, pr, dep, sel_vars, selected, ncat=None, sel_name=None):
    ds = s.data
    yfull = ds.get(dep).data[pr.smp.mask].astype(float)
    y = np.where(selected & (yfull < SYS), yfull, 0.0)
    main = block_from_pr(pr, dep)
    sel_blk = block_from_text(pr, sel_name or "select", sel_vars, constant=True)
    for nm in sel_blk.notes:
        s.output.write(f"note: {nm} omitted because of collinearity\n", "text")
    return y, main, sel_blk


def _selected_subset(blk: Block, m) -> Block:
    return Block(blk.name, blk.cols, blk.terms, blk.constant, blk.X[m])


def _lambda_row(level, ath, lns, V2):
    """lambda = rho sigma with the delta-method standard error."""
    rho, sg = np.tanh(ath), np.exp(lns)
    gr = np.array([(1 - rho ** 2) * sg, rho * sg])
    se = float(np.sqrt(gr @ V2 @ gr))
    c = crit("z", level, None)
    val = float(rho * sg)
    return aux_row("lambda", val, se, level, ci=(val - c * se, val + c * se))


def cmd_heckman(s: "Session", args: str) -> None:
    if _replay(s, "heckman", args):
        return
    spec = {"twostep": 3, "offset": 3, "mills": 3}
    dep, sel_name, sel_vars, pr, selected, o = _selection_setup(s, args, spec)
    if o.get("twostep"):
        return _heckman_twostep(s, args, dep, sel_vars, pr, selected, o, sel_name)
    pr.write_notes(s.output)
    y, main, sel_blk = _selection_common(s, pr, dep, sel_vars, selected, sel_name=sel_name)
    sel = selected.astype(float)
    w = pr.w
    n = len(y)
    sw_ = selected

    def lik(etas, th):
        ll, d_eta, d_zeta, d_ath, d_ls = L.heckman(sel, y, etas[0], etas[1], th[0], th[1])
        return ll, [d_eta, d_zeta], np.column_stack([d_ath, d_ls])

    eng = Engine(lik, [main, sel_blk], w, 2)
    b0 = _heckman_start(y, main, sel_blk, sel, w)
    res = pr.maximize(eng.fun, b0)
    _log(s, pr, res)
    V, n_clust = pr.vce(res.H, eng.scores(res.b) if pr.robust else None)
    ath, lns = float(res.b[-2]), float(res.b[-1])
    V2 = V[-2:, -2:]
    # independent equations: probit on the selection + OLS on the selected observations
    rp = _weighted_probit(sel, sel_blk.X, w)
    Xs, ys, ws = main.X[sw_], y[sw_], w[sw_]
    swr = np.sqrt(ws)
    bo = np.linalg.lstsq(Xs * swr[:, None], ys * swr, rcond=None)[0]
    rss = float(np.sum(ws * (ys - Xs @ bo) ** 2))
    n1 = float(ws.sum())
    ll_ols = -0.5 * n1 * (np.log(2 * np.pi) + 1 + np.log(rss / n1))
    lr = max(2 * (res.ll - (rp.ll + ll_ols)), 0.0)
    p_lr = _chi2p(lr, 1)
    n_cens, n_unc = float(w[~sw_].sum()), float(w[sw_].sum())

    def aux_rows(level):
        return [CoefRow("sep"), aux_row("/athrho", ath, float(np.sqrt(V[-2, -2])), level, name="_cons"),
                aux_row("/lnsigma", lns, float(np.sqrt(V[-1, -1])), level, name="_cons"), CoefRow("sep"),
                transformed_row("rho", np.tanh, lambda x: 1 - np.tanh(x) ** 2, ath, float(np.sqrt(V[-2, -2])), level),
                transformed_row("sigma", np.exp, np.exp, lns, float(np.sqrt(V[-1, -1])), level),
                _lambda_row(level, ath, lns, V2)]

    def hdr(est):
        x = est.extra
        left_ = [est.title, "(regression model with sample selection)", "", "", "",
                 f"{lab(est)} = {g(x['ll'], '%9.0g'):>9}"]
        right = [("Number of obs", comma(x["N"])), ("Censored obs", comma(n_cens)),
                 ("Uncensored obs", comma(n_unc)), ("", ""),
                 (f"{x['chi2type']} chi2({int(x['df_m'])})", g(x["chi2"], "%10.2f")),
                 ("Prob > chi2", g(x["p"], "%10.4f"))]
        return left_, right

    w_rho = (ath / float(np.sqrt(V[-2, -2]))) ** 2
    p_w = _chi2p(w_rho, 1)

    def footer(s_, est):
        if pr.robust:    # with robust or cluster, the LR test is not valid: Wald test of /athrho (compat 1907)
            s_.output.write(f"Wald test of indep. eqns. (rho = 0): chi2(1) = {w_rho:8.2f}   Prob > chi2 = {p_w:6.4f}\n",
                            "text")
        else:
            s_.output.write(f"LR test of indep. eqns. (rho = 0):   chi2(1) = {lr:8.2f}   Prob > chi2 = {p_lr:6.4f}\n",
                            "text")   # VERIFY

    rho, sg = float(np.tanh(ath)), float(np.exp(lns))
    est = finish_model(s, pr, "heckman", args, title="Heckman selection model", blocks=[main, sel_blk],
                       eng=eng, res=res, V=V, n_clust=n_clust, ll0=SYS,
                       aux_names=[("athrho", "_cons"), ("lnsigma", "_cons")], aux_rows_fn=aux_rows, hdr_fn=hdr,
                       depvar_cell="", footer=footer, wald=True, k_eq=4, test_idx=list(range(main.nslopes)),
                       extra_scalars=[("N_selected", n_unc), ("N_nonselected", n_cens), ("rho", rho),
                                      ("sigma", sg), ("lambda", rho * sg), ("chi2_c", lr), ("p_c", p_lr)],
                       predict_opts={"xbsel": 4, "ycond": 2, "yexpected": 2, "psel": 4})
    est.predict = make_predict(_predict_heckman, "xb")
    est.pred_label = "fitted values"
    display_eq(s, est)


def _predict_heckman(what, etas, th, est):
    ath, lns = th[0], th[1]
    rho, sg = np.tanh(ath), np.exp(lns)
    if what == "xb":
        return etas[0]
    if what == "xbsel":
        return etas[1]
    if what == "pr":
        return st.norm.cdf(etas[1])
    if what == "ycond":
        return etas[0] + rho * sg * np.exp(st.norm.logpdf(etas[1]) - st.norm.logcdf(etas[1]))
    if what == "yexpected":
        return st.norm.cdf(etas[1]) * (etas[0] + rho * sg * np.exp(st.norm.logpdf(etas[1]) - st.norm.logcdf(etas[1])))
    return None


def _heckman_start(y, main, sel_blk, sel, w):
    """Two-step estimates as starting values."""
    rp = _weighted_probit(sel, sel_blk.X, w)
    zg = sel_blk.X @ rp.b
    m = sel > 0
    lam = np.exp(st.norm.logpdf(zg) - st.norm.logcdf(zg))
    Xs = np.column_stack([main.X[m], lam[m]])
    sw = np.sqrt(w[m])
    coef = np.linalg.lstsq(Xs * sw[:, None], y[m] * sw, rcond=None)[0]
    beta, theta = coef[:-1], coef[-1]
    resid = y[m] - main.X[m] @ beta
    sig2 = float(np.sum(w[m] * resid ** 2) / np.sum(w[m]))
    sig = np.sqrt(max(sig2, 1e-8))
    rho = float(np.clip(theta / sig, -0.9, 0.9))
    return np.concatenate([beta, rp.b, [np.arctanh(rho), np.log(sig)]])


def _heckman_twostep(s, args, dep, sel_vars, pr, selected, o, sel_name=None):
    """Heckman's two-step estimates: probit for the selection, then OLS of the outcome on the
    covariates and the inverse Mills ratio, with the variance corrected for the estimated
    lambda (Heckman 1979; Greene, Econometric Analysis)."""
    from .ml import MLResult
    from .models2 import _finish  # noqa: F401  (finish_model posts through it)
    if pr.smp.wtype:
        raise StataError(101, "weights not allowed")   # VERIFY
    pr.write_notes(s.output)
    y, main, sel_blk = _selection_common(s, pr, dep, sel_vars, selected, sel_name=sel_name)
    sel = selected.astype(float)
    n = len(y)
    w = np.ones(n)
    ps = _weighted_probit(sel, sel_blk.X, w)
    Vg = np.linalg.inv(-ps.H)
    zg = sel_blk.X @ ps.b
    lam = np.exp(st.norm.logpdf(zg) - st.norm.logcdf(zg))
    delta = lam * (lam + zg)
    m = selected
    n1 = int(m.sum())
    Xs = np.column_stack([main.X[m], lam[m]])
    coef = np.linalg.lstsq(Xs, y[m], rcond=None)[0]
    beta, theta = coef[:-1], coef[-1]
    e = y[m] - Xs @ coef
    sig2 = float(e @ e / n1 + theta ** 2 * delta[m].mean())
    sig = float(np.sqrt(sig2))
    rho = float(np.clip(theta / sig, -1.0, 1.0))
    XtXi = np.linalg.inv(Xs.T @ Xs)
    D = delta[m]
    Fm = Xs.T @ (sel_blk.X[m] * D[:, None])
    mid = Xs.T @ (Xs * (1 - rho ** 2 * D)[:, None]) + rho ** 2 * Fm @ Vg @ Fm.T
    Vs = sig2 * XtXi @ mid @ XtXi
    k1, k2 = main.k, sel_blk.k
    b_all = np.concatenate([beta, ps.b, [theta]])
    V_all = np.zeros((k1 + k2 + 1, k1 + k2 + 1))
    pos_a = list(range(k1)) + [k1 + k2]
    V_all[np.ix_(pos_a, pos_a)] = Vs
    V_all[k1:k1 + k2, k1:k1 + k2] = Vg
    res = MLResult(b_all, SYS, np.zeros(len(b_all)), np.zeros((len(b_all), len(b_all))), 0, True, [])
    n_cens, n_unc = float((~m).sum()), float(n1)
    se_l = float(np.sqrt(Vs[-1, -1]))

    def aux_rows(level):
        return [CoefRow("sep"), CoefRow("eq", "mills"), CoefRow("coef", "lambda", theta, se_l, name="lambda"),
                CoefRow("sep"), CoefRow("value", "rho", rho, eq=f"{rho:.5f}"),
                CoefRow("value", "sigma", sig, eq=g(sig, "%10.8g"))]

    def hdr(est):
        x = est.extra
        left_ = ["Heckman selection model -- two-step estimates", "(regression model with sample selection)",
                 "", "", ""]
        right = [("Number of obs", comma(x["N"])), ("Censored obs", comma(n_cens)),
                 ("Uncensored obs", comma(n_unc)), ("", ""),
                 (f"Wald chi2({int(x['df_m'])})", g(x["chi2"], "%10.2f")), ("Prob > chi2", g(x["p"], "%10.4f"))]
        return left_, right   # VERIFY layout

    eng = Engine(lambda e_, th: (np.zeros(n), [np.zeros(n)], None), [main], w, 0)
    order = (["N", "N_selected", "N_nonselected", "df_m", "chi2", "p", "rho", "sigma", "lambda"],
             ["cmdline", "cmd", "title", "chi2type", "vce", "depvar", "properties"], [])
    idx = list(range(main.nslopes))
    est = finish_model(s, pr, "heckman", args, title="Heckman selection model -- two-step estimates",
                       blocks=[main, sel_blk], eng=eng, res=res, V=V_all, n_clust=None, ll0=SYS,
                       aux_names=[("mills", "lambda")], aux_rows_fn=aux_rows, hdr_fn=hdr, depvar_cell="",
                       wald=True, k_eq=3, test_idx=idx, order=order,
                       extra_scalars=[("N_selected", n_unc), ("N_nonselected", n_cens), ("rho", rho),
                                      ("sigma", sig), ("lambda", theta)],
                       predict_opts={"xbsel": 4, "psel": 4})
    est.predict = make_predict(_predict_heckman_two, "xb")
    est.pred_label = "fitted values"
    display_eq(s, est)


def _predict_heckman_two(what, etas, th, est):
    if what == "xb":
        return etas[0]
    if what == "xbsel":
        return etas[1]
    if what == "pr":
        return st.norm.cdf(etas[1])
    return None


def _sel_header(title, n_cens, n_unc, line2=""):
    def hdr(est):
        x = est.extra
        left_ = [title] + ([line2] if line2 else [""]) + ["", "", "", f"{lab(est)} = {g(x['ll'], '%9.0g'):>9}"]
        right = [("Number of obs", comma(x["N"])), ("Censored obs", comma(n_cens)),
                 ("Uncensored obs", comma(n_unc)), ("", ""),
                 (f"{x['chi2type']} chi2({int(x['df_m'])})", g(x["chi2"], "%10.2f")),
                 ("Prob > chi2", g(x["p"], "%10.4f"))]
        return left_, right   # VERIFY layout
    return hdr


def _rho_rows(ath, se_ath, level):
    return [aux_row("/athrho", ath, se_ath, level, name="_cons"), CoefRow("sep"),
            transformed_row("rho", np.tanh, lambda x: 1 - np.tanh(x) ** 2, ath, se_ath, level)]


def _lr_footer_rho(lr, p_lr, label="LR test of indep. eqns. (rho = 0):  "):
    def footer(s_, est):
        s_.output.write(f"{label} chi2(1) = {lr:8.2f}   Prob > chi2 = {p_lr:6.4f}\n", "text")   # VERIFY
    return footer


# ---------------------------------------------------------------------------
# heckprobit
# ---------------------------------------------------------------------------

def cmd_heckprobit(s: "Session", args: str) -> None:
    if _replay(s, "heckprobit", args):
        return
    dep, sel_name, sel_vars, pr, selected, o = _selection_setup(s, args, {"offset": 3})
    pr.write_notes(s.output)
    y, main, sel_blk = _selection_common(s, pr, dep, sel_vars, selected, sel_name=sel_name)
    yb = np.where(selected, (y != 0).astype(float), 0.0)
    if np.all(yb[selected] == yb[selected][0]):
        raise StataError(2000, "outcome does not vary; remember:\n"
                         "                                  0 = negative outcome,\n"
                         "        all other nonmissing values = positive outcome")
    sel = selected.astype(float)
    w = pr.w

    def lik(etas, th):
        ll, d_eta, d_zeta, d_ath = L.heckprobit(sel, yb, etas[0], etas[1], th[0])
        return ll, [d_eta, d_zeta], d_ath.reshape(-1, 1)

    eng = Engine(lik, [main, sel_blk], w, 1)
    ps = _weighted_probit(sel, sel_blk.X, w)
    po = _weighted_probit(yb[selected], main.X[selected], w[selected])
    res = pr.maximize(eng.fun, np.concatenate([po.b, ps.b, [0.0]]))
    fit_log(s, pr, "Fitting probit model:", po)
    fit_log(s, pr, "Fitting selection model:", ps)
    s.output.write(f"\nComparison:    log likelihood = {g(ps.ll + po.ll, '%10.0g'):>10}\n", "text")
    # starting values: the outcome probit with the inverse Mills ratio of the selection equation
    zg = sel_blk.X @ ps.b
    lam_ = np.exp(st.norm.logpdf(zg) - st.norm.logcdf(zg))
    ysel_ = yb[selected]
    pf_, _ = _probit(ysel_, np.column_stack([main.X[selected], lam_[selected]]), w[selected])
    pst = maximize(pf_, np.concatenate([po.b, [0.0]]), maxiter=300)
    fit_log(s, pr, "Fitting starting values:", pst)
    fit_log(s, pr, "Fitting full model:", res)
    V, n_clust = pr.vce(res.H, eng.scores(res.b) if pr.robust else None)
    ath, se_ath = float(res.b[-1]), float(np.sqrt(V[-1, -1]))
    lr = max(2 * (res.ll - (ps.ll + po.ll)), 0.0)
    p_lr = _chi2p(lr, 1)
    n_cens, n_unc = float(w[~selected].sum()), float(w[selected].sum())
    est = finish_model(s, pr, "heckprobit", args, title="Probit model with sample selection",
                       blocks=[main, sel_blk], eng=eng, res=res, V=V, n_clust=n_clust, ll0=SYS,
                       aux_names=[("athrho", "_cons")], aux_rows_fn=lambda lv: [CoefRow("sep")] + _rho_rows(ath, se_ath, lv),
                       hdr_fn=_sel_header("Probit model with sample selection", n_cens, n_unc), depvar_cell="",
                       footer=_lr_footer_rho(lr, p_lr), wald=True, k_eq=3,
                       test_idx=list(range(main.nslopes)),
                       extra_scalars=[("N_selected", n_unc), ("N_nonselected", n_cens), ("rho", float(np.tanh(ath))),
                                      ("chi2_c", lr), ("p_c", p_lr)],
                       predict_opts={"pmargin": 4, "psel": 4, "xbsel": 4})
    est.predict = make_predict(_predict_heckprobit, "pmargin")
    est.pred_label = f"Pr({est.depvar})"
    display_eq(s, est)


def _predict_heckprobit(what, etas, th, est):
    if what == "xb":
        return etas[0]
    if what == "xbsel":
        return etas[1]
    if what == "pmargin":
        return st.norm.cdf(etas[0])
    if what == "psel":
        return st.norm.cdf(etas[1])
    return None


# ---------------------------------------------------------------------------
# heckoprobit
# ---------------------------------------------------------------------------

def cmd_heckoprobit(s: "Session", args: str) -> None:
    if _replay(s, "heckoprobit", args):
        return
    dep, sel_name, sel_vars, pr, selected, o = _selection_setup(s, args, {"offset": 3}, outcome_constant=False)
    pr.write_notes(s.output)
    y, main, sel_blk = _selection_common(s, pr, dep, sel_vars, selected, sel_name=sel_name)
    levels = sorted(set(y[selected].tolist()))
    J = len(levels)
    if J < 3:
        raise StataError(198, "dependent variable must have at least 3 categories")   # VERIFY
    yi = np.where(selected, np.searchsorted(levels, y), 0)
    sel = selected.astype(float)
    w = pr.w
    n = len(y)

    def lik(etas, th):
        ll, d_eta, d_zeta, d_ath, d_cut = L.heckoprobit(sel, yi, J, etas[0], etas[1], th[:J - 1], th[J - 1])
        return ll, [d_eta, d_zeta], np.column_stack([d_cut, d_ath])

    class _Nan(Exception):
        pass

    eng = Engine(lik, [main, sel_blk], w, J)
    ps = _weighted_probit(sel, sel_blk.X, w)
    cum = np.array([np.sum(w[selected & (yi <= j)]) for j in range(J - 1)]) / np.sum(w[selected])
    cuts0 = st.norm.ppf(np.clip(cum, 1e-4, 1 - 1e-4))
    # ordered probit among the selected: the model of the LR test (the selection index at +inf)
    big = np.full(n, 1e10)
    sub = selected

    def lik_o(etas, th):
        ll, d_eta, _, _, d_cut = L.heckoprobit(np.ones(sub.sum()), yi[sub], J, etas[0], big[sub], th[:J - 1], 0.0)
        return ll, [d_eta], d_cut

    mo = Block(main.name, main.cols, main.terms, main.constant, main.X[sub])
    eo = Engine(lik_o, [mo], w[sub], J - 1)
    ro = maximize(eo.fun, np.concatenate([np.zeros(main.k), cuts0]), maxiter=300)
    res = pr.maximize(eng.fun, np.concatenate([ro.b[:main.k], ps.b, ro.b[main.k:], [0.0]]))
    _log(s, pr, res)
    V, n_clust = pr.vce(res.H, eng.scores(res.b) if pr.robust else None)
    ath, se_ath = float(res.b[-1]), float(np.sqrt(V[-1, -1]))
    lr = max(2 * (res.ll - (ps.ll + ro.ll)), 0.0)
    p_lr = _chi2p(lr, 1)
    n_cens, n_unc = float(w[~selected].sum()), float(w[selected].sum())
    cut_names = [(f"cut{j + 1}", "_cons") for j in range(J - 1)]

    def aux_rows(level):
        rows = [CoefRow("sep")]
        for j in range(J - 1):
            rows.append(aux_row(f"/cut{j + 1}", float(res.b[-J + j]), float(np.sqrt(V[-J + j, -J + j])), level,
                                name="_cons"))
        return rows + _rho_rows(ath, se_ath, level)

    est = finish_model(s, pr, "heckoprobit", args, title="Ordered probit model with sample selection",
                       blocks=[main, sel_blk], eng=eng, res=res, V=V, n_clust=n_clust, ll0=SYS,
                       aux_names=cut_names + [("athrho", "_cons")], aux_rows_fn=aux_rows,
                       hdr_fn=_sel_header("Ordered probit model with sample selection", n_cens, n_unc),
                       depvar_cell="", footer=_lr_footer_rho(lr, p_lr), wald=True, k_eq=3 + J - 1,
                       test_idx=list(range(main.nslopes)),
                       extra_scalars=[("N_selected", n_unc), ("N_nonselected", n_cens), ("k_cat", float(J)),
                                      ("rho", float(np.tanh(ath))), ("chi2_c", lr), ("p_c", p_lr)])
    est.predict = make_predict(lambda what, etas, th, e: etas[0] if what == "xb" else
                               (etas[1] if what == "xbsel" else None), "xb")
    est.pred_label = "linear prediction"
    display_eq(s, est)


# ---------------------------------------------------------------------------
# biprobit
# ---------------------------------------------------------------------------

def _biprobit_equations(args: str):
    """((y1, x1), (y2, x2), rest of the line) for `biprobit (y1 = x1) (y2 = x2) ...`, or None."""
    t = args.strip()
    if not t.startswith("("):
        return None
    eqs = []
    for _ in range(2):
        if not t.startswith("("):
            raise StataError(198, "invalid syntax")
        depth, i = 0, 0
        for i, ch in enumerate(t):
            depth += ch == "("
            depth -= ch == ")"
            if depth == 0:
                break
        body, t = t[1:i], t[i + 1:].lstrip()
        if "=" not in body:
            raise StataError(198, "invalid syntax")
        yv, xv = body.split("=", 1)
        eqs.append((yv.strip(), xv.strip()))
    return eqs[0], eqs[1], t


def cmd_biprobit(s: "Session", args: str) -> None:
    if _replay(s, "biprobit", args):
        return
    eqs = _biprobit_equations(args)
    ds = s.data
    if eqs is None:
        (d1, d2), rest = _split_vars(args, 2)
        x1 = x2 = None
        first_vars = None
        rest_args = rest
        # shared regressors: the varlist after the two dependent variables
        from ..lang.syntax import parse_standard
        shared = parse_standard(rest).varlist
        x1 = x2 = shared
        pr_args = rest
    else:
        (d1, x1), (d2, x2), rest = eqs
        pr_args = f"{x1} {rest}" if not rest.startswith(",") else f"{x1} {rest}"
    bad = (ds.get(d1).data >= SYS) | (ds.get(d2).data >= SYS)
    spec = {"partial": 3, "offset1": 7, "offset2": 7}
    pr = _Prep(s, pr_args, spec, need_dep=False, restrict=~bad, extra_vars=names_in(x2))
    pr.smp.depname = d1
    pr.write_notes(s.output)
    mask = pr.smp.mask
    y1 = (ds.get(d1).data[mask] != 0).astype(float)
    y2 = (ds.get(d2).data[mask] != 0).astype(float)
    for yy in (y1, y2):
        if np.all(yy == yy[0]):
            raise StataError(2000, "outcome does not vary; remember:\n"
                             "                                  0 = negative outcome,\n"
                             "        all other nonmissing values = positive outcome")
    b1 = block_from_pr(pr, d1)
    b2 = block_from_text(pr, d2, x2, constant=pr.constant) if eqs is not None or x2 != x1 else None
    if b2 is None:
        b2 = Block(d2, b1.cols, b1.terms, b1.constant, b1.X.copy())
        # an independent copy of the columns: omitted terms are marked per equation
    for nm in b2.notes:
        s.output.write(f"note: {nm} omitted because of collinearity\n", "text")
    w = pr.w

    def lik(etas, th):
        ll, d1_, d2_, d_ath = L.biprobit(y1, y2, etas[0], etas[1], th[0])
        return ll, [d1_, d2_], d_ath.reshape(-1, 1)

    eng = Engine(lik, [b1, b2], w, 1)
    p1 = _weighted_probit(y1, b1.X, w)
    p2 = _weighted_probit(y2, b2.X, w)
    res = pr.maximize(eng.fun, np.concatenate([p1.b, p2.b, [0.0]]))
    fit_log(s, pr, "Fitting comparison equation 1:", p1)
    fit_log(s, pr, "Fitting comparison equation 2:", p2)
    s.output.write(f"\nComparison:    log likelihood = {g(p1.ll + p2.ll, '%10.0g'):>10}\n", "text")
    fit_log(s, pr, "Fitting full model:", res)
    V, n_clust = pr.vce(res.H, eng.scores(res.b) if pr.robust else None)
    ath, se_ath = float(res.b[-1]), float(np.sqrt(V[-1, -1]))
    lr = max(2 * (res.ll - (p1.ll + p2.ll)), 0.0)
    p_lr = _chi2p(lr, 1)
    title_bp = "Seemingly unrelated bivariate probit" if eqs is not None else "Bivariate probit regression"

    def bp_footer(s_, est):
        s_.output.write(f"Likelihood-ratio test of rho=0:     chi2(1) ={g(lr, '%8.0g'):>9}    Prob > chi2 = {p_lr:6.4f}\n",
                        "text")
    slopes = list(range(b1.nslopes)) + [b1.k + i for i in range(b2.nslopes)]

    def hdr(est):
        x = est.extra
        left_ = [est.title, "", f"{lab(est)} = {g(x['ll'], '%10.0g'):>10}"]
        right = [("Number of obs", comma(x["N"])), (f"{x['chi2type']} chi2({int(x['df_m'])})", g(x["chi2"], "%10.2f")),
                 ("Prob > chi2", g(x["p"], "%10.4f"))]
        return left_, right   # VERIFY layout

    est = finish_model(s, pr, "biprobit", args, title=title_bp, blocks=[b1, b2],
                       eng=eng, res=res, V=V, n_clust=n_clust, ll0=SYS, aux_names=[("athrho", "_cons")],
                       aux_rows_fn=lambda lv: [CoefRow("sep")] + _rho_rows(ath, se_ath, lv), hdr_fn=hdr,
                       depvar_cell="", wald=True, k_eq=3, test_idx=slopes,
                       footer=bp_footer,
                       extra_scalars=[("rho", float(np.tanh(ath))), ("chi2_c", lr), ("p_c", p_lr)],
                       extra_macros=[("depvar", f"{d1} {d2}")],
                       predict_opts={"xb1": 3, "xb2": 3, "pmarg1": 6, "pmarg2": 6, "p11": 3, "p10": 3, "p01": 3, "p00": 3})
    est.depvar = d1
    est.predict = make_predict(_predict_biprobit, "pmarg1")
    est.pred_label = f"Pr({d1})"
    display_eq(s, est)


def _predict_biprobit(what, etas, th, est):
    from ..core.bvnorm import bvn_cdf
    rho = np.tanh(th[0])
    if what == "xb1":
        return etas[0]
    if what == "xb2":
        return etas[1]
    if what == "pmarg1":
        return st.norm.cdf(etas[0])
    if what == "pmarg2":
        return st.norm.cdf(etas[1])
    if what == "p11":
        return bvn_cdf(etas[0], etas[1], rho)
    if what == "p10":
        return bvn_cdf(etas[0], -etas[1], -rho)
    if what == "p01":
        return bvn_cdf(-etas[0], etas[1], -rho)
    if what == "p00":
        return bvn_cdf(-etas[0], -etas[1], rho)
    return None
