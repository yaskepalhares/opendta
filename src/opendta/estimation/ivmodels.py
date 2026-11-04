"""Models with an endogenous regressor: ivprobit, ivtobit and ivpoisson ([R] ivprobit,
[R] ivtobit, [R] ivpoisson). ivprobit and ivtobit are fitted by maximum likelihood (one
endogenous regressor); ivpoisson by GMM or by the control-function approach. The log
likelihoods are in likelihoods.py."""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.stats import chi2_tail
from ..lazy import lazy
from . import likelihoods as L
from .censmodels import _replay, _weighted_probit
from .mleq import (Block, Engine, aux_row, block_from_pr, block_from_text, display_eq, finish_model, fit_log,
                   lab, names_in, peek, transformed_row)
from .models2 import _Prep, _tobit_footer, g
from .regress import comma
from .results import CoefRow

st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def _split_iv(args: str):
    """'y x1 (y2 = z1 z2), opts' -> (dependent, exogenous text, endogenous, instruments, options)."""
    m = re.search(r"\(([^=()]+)=([^()]*)\)", args)
    if not m:
        raise StataError(198, "equation not specified")   # VERIFY
    endog, inst = m.group(1).strip(), m.group(2).strip()
    rest = (args[:m.start()] + " " + args[m.end():]).strip()
    head, _, opts = rest.partition(",")
    parts = head.split(None, 1)
    if not parts:
        raise StataError(102, "too few variables specified")
    return parts[0], (parts[1] if len(parts) > 1 else ""), endog, inst, opts


def _ols(y, X, w):
    sw = np.sqrt(w)
    b = np.linalg.lstsq(X * sw[:, None], y * sw, rcond=None)[0]
    e = y - X @ b
    return b, float(np.sqrt(np.sum(w * e * e) / np.sum(w)))


def _iv_ml(s: "Session", args: str, cmd: str) -> None:
    if _replay(s, cmd, args):
        return
    dep, exog, endog, inst, opts = _split_iv(args)
    tobit = cmd == "ivtobit"
    spec = {"twostep": 4, "first": 3, "mfx": 3}
    if tobit:
        spec.update({"ll": 2, "ul": 2})
    o = peek(f"{dep} , {opts}", spec)
    if o.get("twostep"):
        return _iv_twostep(s, args, cmd, dep, exog, endog, inst, opts, spec)
    if len(endog.split()) != 1:
        # the manual gives the likelihood but not the parameterization of Sigma that Stata reports
        raise StataError(198, "only one endogenous regressor is allowed with maximum likelihood in OpenDTA "
                              "for now; use twostep")   # VERIFY
    pr = _Prep(s, f"{dep} {endog} {exog} , {opts}".replace(", ,", ","), spec, extra_vars=inst)
    pr.write_notes(s.output)
    ds = s.data
    mask = pr.smp.mask
    main = block_from_pr(pr, dep)
    aux = block_from_text(pr, endog, f"{inst} {exog}", constant=pr.constant)
    for nm in aux.notes:
        s.output.write(f"note: {nm} omitted because of collinearity\n", "text")
    y2 = ds.get(endog).data[mask].astype(float)
    w = pr.w
    ylab = pr.y
    if tobit:
        lo = _limit(o.get("ll"))
        hi = _limit(o.get("ul"))
        if lo is None and hi is None:
            raise StataError(198, "must specify at least one of ll() and ul()")   # VERIFY
        lo_ = -np.inf if lo is None else lo
        hi_ = np.inf if hi is None else hi
        y = ylab.astype(float)
        b1, sg = _ols(y, main.X, w)
        b2, sv = _ols(y2, aux.X, w)

        def lik(etas, th):
            ll, d1, d2, da, dls, dlv = L.ivtobit(y, y2, etas[0], etas[1], th[0], th[1], th[2], lo_, hi_)
            return ll, [d1, d2], np.column_stack([da, dls, dlv])

        eng = Engine(lik, [main, aux], w, 3)
        start = np.concatenate([b1, b2, [0.0, np.log(max(sg, 1e-6)), np.log(max(sv, 1e-6))]])
        aux_names = [("alpha", "_cons"), ("lns", "_cons"), ("lnv", "_cons")]
        title = "Tobit model with endogenous regressors"
    else:
        y = (ylab != 0).astype(float)
        if np.all(y == y[0]):
            raise StataError(2000, "outcome does not vary; remember:\n"
                             "                                  0 = negative outcome,\n"
                             "        all other nonmissing values = positive outcome")
        p1 = _weighted_probit(y, main.X, w)
        b2, sv = _ols(y2, aux.X, w)

        def lik(etas, th):
            ll, d1, d2, da, dls = L.ivprobit(y, y2, etas[0], etas[1], th[0], th[1])
            return ll, [d1, d2], np.column_stack([da, dls])

        eng = Engine(lik, [main, aux], w, 2)
        start = np.concatenate([p1.b, b2, [0.0, np.log(max(sv, 1e-6))]])
        aux_names = [("athrho", "_cons"), ("lnsigma", "_cons")]
        title = "Probit model with endogenous regressors"
    res = pr.maximize(eng.fun, start)
    if tobit:
        if not pr.o.get("nolog"):
            s.output.write("\nFitting exogenous tobit model\n", "text")
    else:
        fit_log(s, pr, "Fitting exogenous probit model", p1)
    fit_log(s, pr, "Fitting full model", res)
    V, n_clust = pr.vce(res.H, eng.scores(res.b) if pr.robust else None)
    K = len(res.b)
    n_th = 3 if tobit else 2
    ath, se_ath = float(res.b[K - n_th]), float(np.sqrt(V[K - n_th, K - n_th]))
    w_exog = (ath / se_ath) ** 2 if se_ath > 0 else SYS
    p_exog = float(chi2_tail(w_exog, 1)) if w_exog < SYS else SYS
    lns, se_lns = float(res.b[K - n_th + 1]), float(np.sqrt(V[K - n_th + 1, K - n_th + 1]))
    if tobit:
        lsv, se_lsv = float(res.b[-1]), float(np.sqrt(V[-1, -1]))

        def aux_rows(level):
            return [CoefRow("sep"), aux_row("/alpha", ath, se_ath, level), aux_row("/lns", lns, se_lns, level),
                    aux_row("/lnv", lsv, se_lsv, level), CoefRow("sep"),
                    transformed_row("s", np.exp, np.exp, lns, se_lns, level),
                    transformed_row("v", np.exp, np.exp, lsv, se_lsv, level)]
    else:
        def aux_rows(level):
            return [CoefRow("sep"), aux_row("/athrho", ath, se_ath, level), aux_row("/lnsigma", lns, se_lns, level),
                    CoefRow("sep"),
                    transformed_row("rho", np.tanh, lambda x: 1 - np.tanh(x) ** 2, ath, se_ath, level),
                    transformed_row("sigma", np.exp, np.exp, lns, se_lns, level)]
    slopes = list(range(main.nslopes))
    exog_names = [c.name for c in main.cols if c.values is not None and not c.base and c.name != endog]
    all_inst = " ".join(exog_names + [c.name for c in aux.cols if c.values is not None and not c.base
                                      and c.name not in exog_names])

    def hdr(est):
        x = est.extra
        left_ = [est.title, "", f"{lab(est)} = {g(x['ll'], '%10.0g'):>10}"]
        right = [("Number of obs", comma(x["N"])), (f"Wald chi2({int(x['df_m'])})", g(x["chi2"], "%10.2f")),
                 ("Prob > chi2", g(x["p"], "%10.4f"))]
        return left_, right

    def footer(s_, est):
        s_.output.write(f"Instrumented:  {endog}\nInstruments:   {all_inst}\n", "text")
        s_.output.write("-" * 78 + "\n", "text")
        if tobit:
            s_.output.write(f"Wald test of exogeneity (/alpha = 0): chi2(1) = {w_exog:8.2f}  "
                            f"Prob > chi2 = {p_exog:6.4f}\n\n", "text")
            _tobit_footer(s_, est)
        else:
            s_.output.write(f"Wald test of exogeneity (/athrho = 0): chi2(1) = {w_exog:8.2f} "
                            f"Prob > chi2 = {p_exog:6.4f}\n", "text")

    est = finish_model(s, pr, cmd, args, title=title, blocks=[main, aux], eng=eng, res=res, V=V,
                       n_clust=n_clust, ll0=SYS, aux_names=aux_names, aux_rows_fn=aux_rows, hdr_fn=hdr,
                       depvar_cell="", wald=True, k_eq=len(aux_names) + 2, test_idx=slopes, footer=footer,
                       extra_scalars=([] if tobit else [("rho", float(np.tanh(ath)))]) + [("chi2_exog", w_exog), ("p_exog", p_exog)],
                       extra_macros=[("instd", endog), ("insts", all_inst), ("method", "ml"),
                                     ("depvar", dep)],
                       predict_opts={"xb": 3})
    est.depvar = dep
    if tobit:
        est.cens = (float(w[y <= lo_].sum()), float(w[(y > lo_) & (y < hi_)].sum()), float(w[y >= hi_].sum()),
                    lo, hi)
    full_rows = est.rows_fn
    est.rows_fn = lambda lv: [r for r in full_rows(lv)[:next(i for i, r in enumerate(full_rows(lv)) if r.kind == "sep")]
                              if r.kind != "eq"] + aux_rows(lv)
    from .mleq import make_predict
    est.predict = make_predict(lambda what, etas, th, e: etas[0] if what == "xb" else None, "xb")
    est.pred_label = "linear prediction"
    display_eq(s, est)


def _tobit_fit(y, X, w, lo, hi):
    """Tobit by maximum likelihood: (b, V of b) with V the inverse information of (b, ln sigma)
    restricted to b."""
    from .ml import maximize
    from .models2 import num_hess
    left, right = y <= lo, y >= hi
    mid = ~(left | right)

    def grad_ll(par):
        b, ls = par[:-1], par[-1]
        sg = np.exp(ls)
        eta = X @ b
        z = (y - eta) / sg
        ll = np.where(mid, st.norm.logpdf(z) - ls, 0.0)
        d_eta = np.where(mid, z / sg, 0.0)
        d_ls = np.where(mid, z * z - 1, 0.0)
        if np.isfinite(lo):
            zl = (lo - eta) / sg
            lc = st.norm.logcdf(zl)
            r = np.exp(st.norm.logpdf(zl) - lc)
            ll = np.where(left, lc, ll)
            d_eta = np.where(left, -r / sg, d_eta)
            d_ls = np.where(left, -r * zl, d_ls)
        if np.isfinite(hi):
            zh = (eta - hi) / sg
            lc = st.norm.logcdf(zh)
            r = np.exp(st.norm.logpdf(zh) - lc)
            ll = np.where(right, lc, ll)
            d_eta = np.where(right, r / sg, d_eta)
            d_ls = np.where(right, -r * zh, d_ls)
        return float(w @ ll), np.concatenate([X.T @ (w * d_eta), [w @ d_ls]])

    def fun(par, todo):
        with np.errstate(all="ignore"):
            ll, gr = grad_ll(par)
        if not np.isfinite(ll):
            return -np.inf, None, None
        if todo == 0:
            return ll, None, None
        return ll, gr, num_hess(lambda q: grad_ll(q)[1], par)
    b0, s0 = _ols(y, X, w)
    r = maximize(fun, np.concatenate([b0, [np.log(max(s0, 1e-6))]]), maxiter=300, tol=1e-12, ltol=1e-14, nrtol=1e-12)
    V = np.linalg.inv(-r.H)
    k = X.shape[1]
    return r.b[:k], V[:k, :k]


def _iv_twostep(s: "Session", args: str, cmd: str, dep: str, exog: str, endog: str, inst: str, opts: str,
                spec: dict) -> None:
    """Newey's (1987) minimum chi-squared estimator ([R] ivprobit and [R] ivtobit, Methods and
    formulas): (1) each endogenous regressor on all the exogenous variables X; (2) probit (tobit)
    of y1 on X and the first-stage residuals, giving alpha~ and lambda; (3) the same on z = (y2, x1)
    and the residuals (2SIV), giving beta; (4) Omega = the covariance of alpha~ plus that of the
    regression of y2 (lambda - beta) on X; (5) delta = (D' Omega^-1 D)^-1 D' Omega^-1 alpha~ with
    D = (Pi, I1); (6) the Wald test of lambda = 0 in the 2SIV fit."""
    from .models2 import _header
    tobit = cmd == "ivtobit"
    pr = _Prep(s, f"{dep} {endog} {exog} , {opts}".replace(", ,", ","), spec, extra_vars=inst,
               weights=("fweight", "aweight"))
    if pr.robust:
        raise StataError(198, "vce(robust) not allowed with twostep")   # VERIFY
    pr.write_notes(s.output)
    ds = s.data
    mask = pr.smp.mask
    o = pr.o
    w = pr.w
    endogs = endog.split()
    p = len(endogs)
    main = block_from_pr(pr, dep)                    # z = (y2, x1, _cons)
    xb = block_from_text(pr, "x", f"{exog} {inst}", constant=pr.constant)    # X = (x1, x2, _cons)
    for nm in xb.notes:
        s.output.write(f"note: {nm} omitted because of collinearity\n", "text")
    X = xb.X
    Y2 = np.column_stack([ds.get(v).data[mask].astype(float) for v in endogs])
    Pi = np.column_stack([np.linalg.lstsq(X * np.sqrt(w)[:, None], Y2[:, j] * np.sqrt(w), rcond=None)[0]
                          for j in range(p)])
    Vhat = Y2 - X @ Pi
    # I1: the columns of X that are the included exogenous variables of z (and the constant)
    zcols = [c.name for c in main.cols if c.values is not None]
    xcols = [c.name for c in xb.cols if c.values is not None]
    if pr.constant:
        zcols, xcols = zcols + ["_cons"], xcols + ["_cons"]
    incl = [nm for nm in zcols if nm not in endogs]
    if len(xcols) - len(incl) < p:
        raise StataError(481, "equation not identified; must have at least as many instruments not in\n"
                              "the regression as there are instrumented variables")
    I1 = np.zeros((len(xcols), len(incl)))
    for j, nm in enumerate(incl):
        I1[xcols.index(nm), j] = 1.0
    ylab = pr.y
    if tobit:
        lo = _limit(o.get("ll"))
        hi = _limit(o.get("ul"))
        if lo is None and hi is None:
            raise StataError(198, "must specify at least one of ll() and ul()")   # VERIFY
        lo_ = -np.inf if lo is None else lo
        hi_ = np.inf if hi is None else hi
        y = ylab.astype(float)

        def fit(A):
            return _tobit_fit(y, A, w, lo_, hi_)
    else:
        y = (ylab != 0).astype(float)
        if np.all(y == y[0]):
            raise StataError(2000, "outcome does not vary")
        s.output.write("Checking reduced-form model...\n", "text")

        def fit(A):
            from .ml import maximize
            from .models import _probit
            fun, _ = _probit(y, A, w)
            r = maximize(fun, np.zeros(A.shape[1]), maxiter=300, tol=1e-12, ltol=1e-14, nrtol=1e-12)
            return r.b, np.linalg.inv(-r.H)
    if o.get("first"):
        # the first-stage regressions, on the excluded instruments and then the included exogenous
        # variables (compat 1926)
        from ..commands.syntaxcmd import _temp_name
        from ..core.dataset import Variable
        tmp = _temp_name(s)
        ds.add(Variable(tmp, "byte", mask.astype(float)))
        try:
            for v in endogs:
                s.output.write("First-stage regression\n", "text")
                s.interp._execute_expanded(f"regress {v} {inst} {exog} if {tmp}")
        finally:
            ds.drop_vars([tmp])
    kx = X.shape[1]
    # (2) alpha~ and lambda
    a_full, Va_full = fit(np.column_stack([X, Vhat]))
    alpha, lam = a_full[:kx], a_full[kx:]
    J_inv = Va_full[:kx, :kx]
    # (3) 2SIV on z = (y2, x1) with the residuals; z's columns are ordered as in main (endogenous first)
    Zm = main.X
    zpos = {nm: j for j, nm in enumerate(zcols)}
    b2, V2 = fit(np.column_stack([Zm, Vhat]))
    kz = Zm.shape[1]
    beta = np.array([b2[zpos[v]] for v in endogs])
    lam2, Vlam2 = b2[kz:], V2[kz:, kz:]
    # (4) Omega
    u = Y2 @ (lam - beta)
    sw = np.sqrt(w)
    XtX = (X * w[:, None]).T @ X
    cu = np.linalg.solve(XtX, X.T @ (w * u))
    e = u - X @ cu
    N = float(w.sum())
    s2 = float(w @ (e * e)) / (N - kx)                  # VERIFY: the residual degrees of freedom
    Omega = J_inv + s2 * np.linalg.inv(XtX)
    # (5) delta, ordered as z (endogenous, included exogenous, _cons)
    D = np.column_stack([Pi, I1])
    Oi = np.linalg.inv(Omega)
    Vd = np.linalg.inv(D.T @ Oi @ D)
    delta = Vd @ D.T @ Oi @ alpha
    order = [p + incl.index(nm) if nm not in endogs else endogs.index(nm) for nm in zcols]
    b = delta[order]
    V = Vd[np.ix_(order, order)]
    V = (V + V.T) / 2
    # (6) exogeneity
    w_exog = float(lam2 @ np.linalg.solve(Vlam2, lam2))
    p_exog = float(chi2_tail(w_exog, p))
    slopes = [j for j, nm in enumerate(zcols) if nm != "_cons"]
    Vs = V[np.ix_(slopes, slopes)]
    chi2 = float(b[slopes] @ np.linalg.solve(Vs, b[slopes]))
    df_m = len(slopes)
    pv = float(chi2_tail(chi2, df_m))
    exog_names = [nm for nm in incl if nm != "_cons"]
    all_inst = " ".join(exog_names + [nm for nm in xcols if nm != "_cons" and nm not in exog_names])
    from .results import Estimates, coef_table, post, rows_from_columns
    from .regress import cmdline
    bmap = {nm: float(b[j]) for j, nm in enumerate(zcols)}
    semap = {nm: float(np.sqrt(max(V[j, j], 0))) for j, nm in enumerate(zcols)}
    rows = rows_from_columns(main.cols, bmap, semap, ds, constant=pr.constant)
    names = [c.name for c in main.cols if c.values is not None] + (["_cons"] if pr.constant else [])
    est = Estimates(cmd, dep, names, b, V, int(N), "z", None)
    est.level = pr.smp.level
    est.vcetype, est.clustvar, est.n_clust = "", "", None
    est.rows = rows
    title = "Two-step tobit with endogenous regressors" if tobit else "Two-step probit with endogenous regressors"
    est.title = title
    cens = None
    if tobit:
        cens = (float(w[y <= lo_].sum()), float(w[(y > lo_) & (y < hi_)].sum()), float(w[y >= hi_].sum()), lo, hi)
        est.cens = cens

    def display(s_, e_, *, header=True, table=True, level=None, **_kw):
        out = s_.output
        if header:
            # compat 1926: the labels start in column 51
            out.write("\n", "text")
            for left, lab_, val in ((title, "Number of obs", comma(N)), ("", f"Wald chi2({df_m})", g(chi2, "%10.2f")),
                                    ("", "Prob > chi2", g(pv, "%10.4f"))):
                out.write(f"{left:<50}{lab_:<16}= ", "text")
                out.write(f"{val:>10}\n", "result")
            out.write("\n", "text")
        if not table:
            return
        coef_table(s_, "", e_.rows, stat="z", df=None, level=level or e_.level, vcetype="")
        out.write(f"Instrumented:  {endog}\nInstruments:   {all_inst}\n", "text")
        out.write("-" * 78 + "\n", "text")
        out.write(f"Wald test of exogeneity:{'':5}chi2({p}) = {w_exog:8.2f}{'':11}Prob > chi2 = {p_exog:6.4f}\n",
                  "text")
        if tobit:
            nl_, nu_, nr_, lo_v, hi_v = e_.cens
            out.write("\n", "text")
            out.write(f"{comma(nl_):>14}  left-censored observations"
                      + (f" at {dep} <= {g(lo_v, '%9.0g')}" if lo_v is not None and nl_ else "") + "\n", "text")
            out.write(f"{comma(nu_):>14}     uncensored observations\n", "text")
            out.write(f"{comma(nr_):>14} right-censored observations"
                      + (f" at {dep} >= {g(hi_v, '%9.0g')}" if hi_v is not None and nr_ else "") + "\n", "text")
    est.display = display
    # Stata's lists and order (compat 1926)
    scalars = [("N", N), ("chi2_exog", w_exog), ("df_exog", float(p)), ("p_exog", p_exog), ("chi2", chi2),
               ("df_m", float(df_m)), ("p", pv)]
    rank = float(np.linalg.matrix_rank(V))
    if tobit:
        nl_, nu_, nr_ = cens[:3]
        scalars += [("endog_ct", float(p)), ("rank", rank), ("N_unc", nu_), ("N_lc", nl_), ("N_rc", nr_)]
        if lo is not None:
            scalars.append(("llopt", float(lo)))
        if hi is not None:
            scalars.append(("ulopt", float(hi)))   # VERIFY
        macros = [("cmdline", cmdline(cmd, args)), ("cmd", cmd),
                  ("marginsok", "XB E(passthru) Pr(passthru) YStar(passthru) default"), ("predict", "ivtobit_p"),
                  ("footnote", "ivtobit_footnote")]
    else:
        scalars += [("rank", rank), ("N_cds", 0.0), ("N_cdf", 0.0)]
        macros = [("cmdline", cmdline(cmd, args)), ("cmd", cmd), ("marginsok", "Pr XB default"),
                  ("predict", "ivprobit_p"), ("estat_cmd", "ivprobit_estat"), ("footnote", "ivprobit_footnote"),
                  ("vce", "twostep")]
    macros += [("depvar", f"{dep} {endog}"), ("insts", all_inst), ("instd", endog), ("method", "twostep"),
               ("properties", "b V")]
    extra = []
    if not tobit:
        # perfect-predictor rules of probit: none dropped here (compat 1926: 2 x 4 with two
        # regressors; VERIFY the contents)
        from ..commands.matrix import Matrix
        extra = [("rules", Matrix(np.zeros((len(slopes), 4))))]
    post(s, est, scalars, macros, mask, extra_mats=extra)
    display(s, est)


def _limit(v):
    if v in (None, True, ""):
        return None
    try:
        return float(str(v).strip())
    except ValueError:
        raise StataError(198, "ll() and ul() must be numbers in OpenDTA")   # VERIFY


def cmd_ivprobit(s: "Session", args: str) -> None:
    _iv_ml(s, args, "ivprobit")


def cmd_ivtobit(s: "Session", args: str) -> None:
    _iv_ml(s, args, "ivtobit")


# ---------------------------------------------------------------------------
# ivpoisson
# ---------------------------------------------------------------------------

def _cluster_meat(G: np.ndarray, cl: np.ndarray | None) -> np.ndarray:
    """Sum of the outer products of the observation moments (summed within clusters)."""
    if cl is None:
        return G.T @ G
    u, inv = np.unique(cl, return_inverse=True)
    T = np.zeros((len(u), G.shape[1]))
    np.add.at(T, inv, G)
    return T.T @ T


def _poisson_start(y, X, w=None):
    b = np.zeros(X.shape[1])
    b[-1] = np.log(max(float(np.mean(y)), 1e-8))
    for _ in range(100):
        mu = np.exp(np.clip(X @ b, -30, 30))
        step = np.linalg.solve((X * mu[:, None]).T @ X + 1e-10 * np.eye(len(b)), X.T @ (y - mu))
        b = b + step
        if np.max(np.abs(step)) < 1e-10:
            break
    return b


def cmd_ivpoisson(s: "Session", args: str) -> None:
    from . import fvars
    from .postest import replay
    from .regress import Sample, cmdline
    from .results import Estimates, coef_table, post, rows_from_columns
    from .models2 import _header
    t = args.strip()
    if not t or t.startswith(","):
        replay(s, "ivpoisson", args)
        return
    m = re.match(r"^(\w+)\s+(.*)$", t, re.S)
    if not m or m.group(1) not in ("gmm", "cfunction"):
        raise StataError(198, "estimator must be gmm or cfunction")   # VERIFY
    method, rest = m.group(1), m.group(2)
    dep, exog_txt, endog_txt, inst_txt, opts = _split_iv(rest)
    ds = s.data
    spec = {"multiplicative": 4, "additive": 3, "onestep": 4, "twostep": 4, "wmatrix": 2, "first": 3,
            "irr": 3, "exposure": 3, "offset": 3}
    smp_args = f"{dep} {exog_txt}" + (f", {opts}" if opts.strip() else "")
    smp = Sample(s, smp_args, opts_spec=spec, weights=("fweight", "pweight", "aweight"),
                 extra_vars=endog_txt + " " + inst_txt)
    mult = bool(smp.o.get("multiplicative"))
    if mult and method == "cfunction":
        raise StataError(198, "option multiplicative is allowed only with gmm")   # VERIFY
    mask = smp.mask
    y = smp.y.astype(float)
    if np.any(y < 0):
        raise StataError(459, "depvar must be nonnegative")   # VERIFY
    N = len(y)
    endog_terms = fvars.expand_fv(ds, endog_txt)
    inst_terms = fvars.expand_fv(ds, inst_txt)
    endog_cols = fvars.build_design(ds, endog_terms, mask)
    inst_cols = fvars.build_design(ds, inst_terms, mask)
    ex_cols = [c for c in smp.cols if c.values is not None and not c.base]
    Xn = [c for c in endog_cols if c.values is not None and not c.base]
    Zc = [c for c in inst_cols if c.values is not None and not c.base]
    one = np.ones(N)
    X = np.column_stack([c.values for c in Xn] + [c.values for c in ex_cols] + [one])
    Z = np.column_stack([c.values for c in Zc] + [c.values for c in ex_cols] + [one])
    k = X.shape[1]
    if Z.shape[1] < k:
        raise StataError(481, "equation not identified; must have at least as many instruments not in\n"
                         "the regression as there are instrumented variables")
    cl = smp.cluster if smp.vce == "cluster" else None
    names = [c.name for c in Xn] + [c.name for c in ex_cols] + ["_cons"]
    out = s.output
    gcrit = None
    wm_label = ""
    if method == "gmm":
        def resid(b):
            xb = np.exp(np.clip(X @ b, -30, 30))
            return (y / xb - 1.0, -(y / xb)) if mult else (y - xb, -xb)

        def moments(b):
            r, _ = resid(b)
            return Z * r[:, None]

        def step(b, W, label):
            out.write(f"\n{label}\n", "text")
            q_old = None
            for it in range(100):
                r, dr = resid(b)
                gbar = (Z * r[:, None]).mean(0)
                G = (Z * dr[:, None]).T @ X / N
                q = float(gbar @ W @ gbar)
                out.write(f"Iteration {it}:   GMM criterion Q(b) = {g(q, '%10.0g'):>10}\n", "text")
                if q_old is not None and abs(q_old - q) < 1e-9 * (1 + abs(q)) or q < 1e-30:
                    break
                d = np.linalg.solve(G.T @ W @ G + 1e-12 * np.eye(k), G.T @ W @ gbar)
                lam, q_new = 1.0, None
                while lam > 1e-6:                      # step halving
                    bn = b - lam * d
                    rn, _ = resid(bn)
                    gn = (Z * rn[:, None]).mean(0)
                    q_new = float(gn @ W @ gn)
                    if np.isfinite(q_new) and q_new <= q + 1e-12:
                        break
                    lam /= 2
                q_old = q
                b = bn
            return b

        W1 = np.linalg.inv(Z.T @ Z / N)
        b0 = _poisson_start(y, X)
        b = step(b0, W1, "Step 1")
        if smp.o.get("onestep"):
            W = W1
            wm_label = "Unadjusted"
        else:
            S = _cluster_meat(moments(b), cl) / N
            W = np.linalg.inv(S)
            b = step(b, W, "Step 2")
            wm_label = f"Cluster ({smp.clustvar})" if cl is not None else "Robust"
        r, dr = resid(b)
        gbar = (Z * r[:, None]).mean(0)
        G = (Z * dr[:, None]).T @ X / N
        Q = float(gbar @ W @ gbar)
        S2 = _cluster_meat(moments(b), cl) / N
        Ainv = np.linalg.inv(G.T @ W @ G)
        V = Ainv @ (G.T @ W @ S2 @ W @ G) @ Ainv / N
        gcrit = Q
        title = "Exponential mean model with endogenous regressors"
        rows_names = names
        bfull = b
        Vfull = V
        n_par, n_mom = k, Z.shape[1]
        eqs_c = None
    else:
        # control function: first-stage residual enters the exponential mean; fitted as an
        # exactly identified GMM problem (first-stage and mean moments together)
        E = len(Xn)
        y2s = [c.values for c in Xn]
        kz = Z.shape[1]
        pis = [np.linalg.lstsq(Z, y2_, rcond=None)[0] for y2_ in y2s]

        def moms(theta):
            vs = [y2s[e] - Z @ theta[e * kz:(e + 1) * kz] for e in range(E)]
            be = theta[E * kz:]
            X2_ = np.column_stack([X] + vs)
            mu = np.exp(np.clip(X2_ @ be, -30, 30))
            return np.column_stack([Z * v_[:, None] for v_ in vs]
                                   + [X2_ * (y / mu - 1.0)[:, None]])   # multiplicative residual

        def jac(theta):
            h = 1e-6
            A_ = np.zeros((len(theta), len(theta)))
            for j in range(len(theta)):
                tp, tm = theta.copy(), theta.copy()
                tp[j] += h
                tm[j] -= h
                A_[:, j] = (moms(tp).mean(0) - moms(tm).mean(0)) / (2 * h)
            return A_
        theta = np.concatenate(pis + [_poisson_start(y, X), np.zeros(E)])
        M0 = moms(theta)
        W1 = np.linalg.inv(M0.T @ M0 / N + 1e-12 * np.eye(M0.shape[1]))
        out.write("\nStep 1\n", "text")
        q_old = None
        for it in range(50):
            gb = moms(theta).mean(0)
            q = float(gb @ W1 @ gb)
            out.write(f"Iteration {it}:   GMM criterion Q(b) = {g(q, '%10.0g'):>10}\n", "text")
            if q < 1e-14 or (q_old is not None and abs(q_old - q) < 1e-12 * (1 + abs(q))):
                break
            d = np.linalg.lstsq(jac(theta), gb, rcond=None)[0]
            lam, tn = 1.0, theta
            while lam > 1e-8:
                tn = theta - lam * d
                gn = moms(tn).mean(0)
                if np.isfinite(gn).all() and float(gn @ W1 @ gn) <= q:
                    break
                lam /= 2
            q_old, theta = q, tn
        out.write("\nnote: model is exactly identified\n", "text")
        theta_init = np.concatenate(pis + [_poisson_start(y, X), np.zeros(E)])
        gbf = moms(theta).mean(0)
        Q = float(gbf @ W1 @ gbf)
        Gm = moms(theta)
        A = jac(theta)
        B = _cluster_meat(Gm, cl) / N
        Ainv = np.linalg.inv(A)
        Vall = Ainv @ B @ Ainv.T / N
        wm_label = f"Cluster ({smp.clustvar})" if cl is not None else "Robust"
        n_ex = len(ex_cols)
        b0_ = E * kz                                  # start of the mean-equation coefficients
        main_idx = [b0_ + E + i for i in range(n_ex)] + [b0_ + j for j in range(E)] + [b0_ + E + n_ex]
        first_idx = [kz - 1 - n_ex + i for i in range(n_ex)] + list(range(len(Zc))) + [kz - 1]
        c_idx = [b0_ + k + j for j in range(E)]
        order = main_idx + [e * kz + i for e in range(E) for i in first_idx] + c_idx
        eqs_c = ([smp.depname] * len(main_idx) + [Xn[e].name for e in range(E) for _ in first_idx]
                 + ["c_" + Xn[e].name for e in range(E)])
        nm_main = [c.name for c in ex_cols] + [c.name for c in Xn] + ["_cons"]
        nm_first = [c.name for c in ex_cols] + [c.name for c in Zc] + ["_cons"]
        rows_names = nm_main + nm_first * E + ["_cons"] * E
        bfull = theta[order]
        Vfull = Vall[np.ix_(order, order)]
        title = "Exponential mean model with endogenous regressors"
        n_par = n_mom = len(theta)
    Vfull = (Vfull + Vfull.T) / 2
    kk = len(bfull)
    bmap = dict(zip(rows_names, bfull))
    semap = {nm: float(np.sqrt(max(Vfull[i, i], 0))) for i, nm in enumerate(rows_names)}
    allcols = Xn + ex_cols
    if method == "gmm":
        rows = rows_from_columns(allcols, {nm: bmap[nm] for nm in rows_names if nm in bmap},
                                 semap, ds, constant=True)
    else:
        se_all = np.sqrt(np.maximum(np.diag(Vfull), 0))
        nmain, nfirst = len(nm_main), len(nm_first)

        def eq_rows(cols, lo, hi):
            bm = {nm: float(bfull[lo + i]) for i, nm in enumerate(rows_names[lo:hi])}
            sm = {nm: float(se_all[lo + i]) for i, nm in enumerate(rows_names[lo:hi])}
            return rows_from_columns(cols, bm, sm, ds, constant=True)
        rows = [CoefRow("eq", smp.depname)] + eq_rows(ex_cols + Xn, 0, nmain)
        for e_ in range(len(Xn)):
            lo_ = nmain + e_ * nfirst
            rows += [CoefRow("sep"), CoefRow("eq", Xn[e_].name)] + eq_rows(ex_cols + Zc, lo_, lo_ + nfirst)
        rows.append(CoefRow("sep"))
        for e_ in range(len(Xn)):
            j_ = len(bfull) - len(Xn) + e_
            rows.append(aux_row("/c_" + Xn[e_].name, float(bfull[j_]), float(se_all[j_]), smp.level))
    est = Estimates("ivpoisson", smp.depname, rows_names, bfull, Vfull, N, "z", None)
    est.rows = rows
    if eqs_c is not None:
        est.eqnames = eqs_c
    est.cols = Xn + ex_cols
    est.terms = endog_terms + smp.terms
    est.constant, est.depcomp = True, smp.depcomp
    est.level = smp.level
    est.vcetype = "Robust"
    est.clustvar, est.n_clust = smp.clustvar, (len(np.unique(cl)) if cl is not None else None)
    inst_all = " ".join([c.name for c in ex_cols] + [c.name for c in Zc])
    est.extra = dict(N=N, df_m=kk, instd=" ".join(c.name for c in Xn), insts=inst_all,
                     title=title, wm_label=wm_label, Q=gcrit, n_par=n_par, n_mom=n_mom)
    from .panel import _retarget
    _retarget(est, Xn, ex_cols, endog_terms, smp.terms)

    def display(s_, e, *, header=True, table=True, level=None, **_kw):
        o_ = s_.output
        x = e.extra
        if header:
            o_.write(f"\n{x['title']}\n\n", "text")
            o_.write(f"Number of parameters = {x['n_par']:>3}{' ' * 25}Number of obs  = {comma(x['N']):>10}\n", "text")
            o_.write(f"Number of moments    = {x['n_mom']:>3}\n", "text")
            o_.write(f"Initial weight matrix: Unadjusted\n", "text")
            o_.write(f"GMM weight matrix:     {x['wm_label']}\n\n", "text")   # VERIFY onestep
        if not table:
            return
        if e.clustvar:
            # compat 1920
            o_.write(f"{'(Std. Err. adjusted for ' + format(int(e.n_clust), ',') + ' clusters in ' + e.clustvar + ')':>78}\n",
                     "text")
        coef_table(s_, e.depvar, e.rows, stat="z", df=None, level=level or e.level, vcetype=e.vcetype)
        o_.write(f"Instrumented:  {x['instd']}\n", "text")
        o_.write(f"Instruments:   {x['insts']}\n", "text")
    est.display = display
    est.predict_opts = {"n": 3, "xb": 3}
    scalars = [("N", float(N)), ("k", float(kk)), ("rank", float(kk))]
    if gcrit is not None:
        scalars += [("Q", gcrit), ("J", N * gcrit)]
    macros = [("cmdline", cmdline("ivpoisson", args)), ("cmd", "ivpoisson"), ("method", method),
              ("depvar", smp.depname), ("instd", est.extra["instd"]), ("insts", inst_all),
              ("title", title), ("properties", "b V")]   # VERIFY e() lists
    extra = []
    if method == "gmm":
        macros.append(("wmatrix", wm_label.lower()))
    else:
        # Stata's lists and order for cfunction (compat 1920)
        E_ = len(Xn)
        vce_ = "cluster" if cl is not None else "robust"
        scalars = [("N", float(N)), ("k", float(kk)), ("k_eq", float(1 + 2 * E_)), ("k_aux", float(E_)),
                   ("k_dv", float(1 + E_)), ("Q", Q), ("J", N * Q), ("J_df", 0.0), ("rank", float(kk))]
        if cl is not None:
            scalars.append(("N_clust", float(len(np.unique(cl)))))   # VERIFY place
        scalars.append(("converged", 1.0))
        macros = [("cmd", "ivpoisson"), ("cmdline", cmdline("ivpoisson", args)), ("depvar", smp.depname),
                  ("instd", est.extra["instd"]), ("insts", inst_all), ("estimator", "cfunction"),
                  ("multiplicative", "multiplicative"), ("title", title), ("winit", "Unadjusted"),
                  ("gmmestimator", "onestep"), ("wmatrix", vce_), ("vce", vce_), ("vcetype", "Robust")]
        if cl is not None:
            macros.append(("clustvar", smp.clustvar))   # VERIFY place
        macros += [("technique", "gn"), ("properties", "b V"), ("estat_cmd", "ivpoisson_estat"),
                   ("predict", "ivpoisson_p"), ("footnote", "ivpoisson_footnote"),
                   ("marginsok", "xb default n pr XBTotal")]
        from ..commands.matrix import Matrix
        extra = [("init", Matrix(theta_init[order].reshape(1, -1))), ("W", Matrix(W1)), ("S", Matrix(B)),
                 ("V_modelbased", Matrix(Vfull))]   # VERIFY contents (only the sizes are in the log)
    post(s, est, scalars, macros, mask, extra_mats=extra)
    display(s, est)
