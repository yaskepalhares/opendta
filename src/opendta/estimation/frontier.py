"""frontier ([R] frontier): stochastic frontier models by maximum likelihood, with half-normal,
exponential or truncated-normal inefficiency. Layout from the manual; lines marked VERIFY wait for
the Stata 14 log of case 1908."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.stats import chi2_tail
from ..lazy import lazy
from .ml import maximize
from .mleq import (Engine, aux_row, block_from_pr, constant_only, display_eq, finish_model, lab, peek, std_right)
from .models2 import _Prep, g
from .results import CoefRow
from .regress import comma

special = lazy("scipy.special")
st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS

_TITLES = {"hnormal": "Stoc. frontier normal/half-normal model",
           "exponential": "Stoc. frontier normal/exponential model",
           "tnormal": "Stoc. frontier normal/truncated-normal model"}


def _ll(e: np.ndarray, lsv: float, lsu: float, mu: float, dist: str) -> np.ndarray:
    """Log density of the composed error e = v - u (u subtracted after the sign for costs)."""
    sv2, su2 = np.exp(lsv), np.exp(lsu)
    sv, su = np.sqrt(sv2), np.sqrt(su2)
    s2 = sv2 + su2
    if dist == "hnormal":
        lam = su / sv
        return 0.5 * np.log(2 / np.pi) - 0.5 * np.log(s2) + special.log_ndtr(-e * lam / np.sqrt(s2)) - e * e / (2 * s2)
    if dist == "exponential":
        return -np.log(su) + 0.5 * sv2 / su2 + special.log_ndtr(-(e + sv2 / su) / sv) + e / su
    mt = (-e * su2 + mu * sv2) / s2
    st_ = su * sv / np.sqrt(s2)
    return (-0.5 * np.log(2 * np.pi) - 0.5 * np.log(s2) - special.log_ndtr(mu / su)
            + special.log_ndtr(mt / st_) - 0.5 * (e + mu) ** 2 / s2)


def cmd_frontier(s: "Session", args: str) -> None:
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, "frontier", args)
        return
    spec = {"distribution": 3, "cost": 4, "uhet": 4, "vhet": 4, "from": 4, "cm": 2, "constraints": 3}
    pr = _Prep(s, args, spec)
    pr.write_notes(s.output)
    o = pr.o
    for bad in ("uhet", "vhet", "cm", "constraints"):
        if o.get(bad):
            raise StataError(198, f"option {bad}() is not available yet in frontier")   # VERIFY
    dist = str(o.get("distribution", "hnormal")).strip().lower()
    for full in ("hnormal", "exponential", "tnormal"):
        if dist and full.startswith(dist):
            dist = full
            break
    else:
        raise StataError(198, "option distribution() invalid")   # VERIFY
    sign = -1.0 if o.get("cost") else 1.0
    y = pr.y
    main = block_from_pr(pr, pr.smp.depname)
    w = pr.w
    n = len(y)
    n_theta = 3 if dist == "tnormal" else 2

    def lik(etas, th):
        th = np.asarray(th, dtype=float)
        e0 = sign * (y - etas[0])
        mu = th[2] if dist == "tnormal" else 0.0
        ll = _ll(e0, th[0], th[1], mu, dist)
        h = 1e-6
        d_e = (_ll(e0 + h, th[0], th[1], mu, dist) - _ll(e0 - h, th[0], th[1], mu, dist)) / (2 * h)
        d_eta = -sign * d_e
        dth = np.zeros((n, n_theta))
        for a in range(n_theta):
            tu, td = th.copy(), th.copy()
            tu[a] += h
            td[a] -= h
            mu_u = tu[2] if dist == "tnormal" else 0.0
            mu_d = td[2] if dist == "tnormal" else 0.0
            dth[:, a] = (_ll(e0, tu[0], tu[1], mu_u, dist) - _ll(e0, td[0], td[1], mu_d, dist)) / (2 * h)
        return ll, [d_eta], dth

    eng = Engine(lik, [main], w, n_theta)
    # starting values: OLS, with the residual variance split between v and u
    b_ols = np.linalg.lstsq(main.X, y, rcond=None)[0]
    r = y - main.X @ b_ols
    s2 = float(r @ r / max(n - main.k, 1))
    th0 = np.array([np.log(s2 / 2), np.log(s2 / 2)] + ([0.0] if dist == "tnormal" else []))
    b0 = np.concatenate([b_ols, th0])
    from scipy import optimize as opt
    try:
        rr = opt.minimize(lambda b: tuple(-v for v in eng.fun(b, 1)[:2]), b0, jac=True, method="BFGS",
                          options={"gtol": 1e-7, "maxiter": 500})
        if np.isfinite(rr.fun):
            b0 = rr.x
    except Exception:
        pass
    res = pr.maximize(eng.fun, b0)
    V, n_clust = pr.vce(res.H, eng.scores(res.b) if pr.robust else None)
    # the model with sigma_u = 0 is least squares
    ll_ols = float(-0.5 * n * (np.log(2 * np.pi) + np.log(float(r @ r) / n) + 1))
    lsv, lsu = float(res.b[main.k]), float(res.b[main.k + 1])
    lr = max(2 * (res.ll - ll_ols), 0.0)
    p_lr = 0.5 * float(chi2_tail(lr, 1)) if lr > 0 else 0.5
    off = main.k
    sv2, su2 = np.exp(lsv), np.exp(lsu)
    Vt = V[off:off + 2, off:off + 2]

    def delta(f_val, grad):
        g_ = np.asarray(grad, dtype=float)
        return float(np.sqrt(max(g_ @ Vt @ g_, 0)))

    def aux_rows(level):
        from .results import crit
        c = crit("z", level, None)
        rows = [CoefRow("sep"),
                aux_row("/lnsig2v", lsv, float(np.sqrt(V[off, off])), level),
                aux_row("/lnsig2u", lsu, float(np.sqrt(V[off + 1, off + 1])), level)]
        if dist == "tnormal":
            rows.append(aux_row("/mu", float(res.b[off + 2]), float(np.sqrt(V[off + 2, off + 2])), level))
        rows.append(CoefRow("sep"))
        stats_ = [("sigma_v", np.exp(lsv / 2), [0.5 * np.exp(lsv / 2), 0.0]),
                  ("sigma_u", np.exp(lsu / 2), [0.0, 0.5 * np.exp(lsu / 2)]),
                  ("sigma2", sv2 + su2, [sv2, su2]),
                  ("lambda", np.exp((lsu - lsv) / 2), [-0.5 * np.exp((lsu - lsv) / 2), 0.5 * np.exp((lsu - lsv) / 2)])]
        for nm, val, gr in stats_:
            se_ = delta(val, gr)
            if nm.startswith("sigma_"):
                # the interval is the transformed interval of the log variance
                i = 0 if nm == "sigma_v" else 1
                l_, se_l = (lsv, np.sqrt(V[off, off])) if i == 0 else (lsu, np.sqrt(V[off + 1, off + 1]))
                rows.append(aux_row(nm, float(val), float(se_), level,
                                    ci=(float(np.exp((l_ - c * se_l) / 2)), float(np.exp((l_ + c * se_l) / 2)))))
            else:
                rows.append(aux_row(nm, float(val), float(se_), level, ci=(float(val - c * se_), float(val + c * se_))))
        return rows

    def hdr(est):
        x = est.extra
        ll4 = g(x["ll"], "%9.0g")
        left = [est.title, "", f"{lab(est)} = {ll4:>9}"]
        right = [("Number of obs", comma(x["N"])), (f"Wald chi2({int(x['df_m'])})", g(x["chi2"], "%10.2f")),
                 ("Prob > chi2", g(x["p"], "%10.4f"))]
        return left, right

    def footer(s_, est):
        s_.output.write(f"LR test of sigma_u=0: chibar2(01) = {lr:.2f}"
                        f"{' ' * 16}Prob>=chibar2 = {p_lr:.3f}\n", "text")   # VERIFY layout

    names = [("lnsig2v", "_cons"), ("lnsig2u", "_cons")] + ([("mu", "_cons")] if dist == "tnormal" else [])
    est = finish_model(s, pr, "frontier", args, title=_TITLES[dist], blocks=[main], eng=eng, res=res, V=V,
                       n_clust=n_clust, ll0=SYS, aux_names=names, aux_rows_fn=aux_rows, hdr_fn=hdr, wald=True,
                       footer=footer, k_eq=1 + len(names),
                       extra_scalars=[("sigma_u", float(np.exp(lsu / 2))), ("sigma_v", float(np.exp(lsv / 2))),
                                      ("lambda", float(np.exp((lsu - lsv) / 2))), ("sigma2", float(sv2 + su2)),
                                      ("chi2_c", lr), ("p_c", p_lr)],
                       extra_macros=[("function", "cost" if o.get("cost") else "production"), ("dist", dist)])
    est.coef_title, est.eform_on = "Coef.", False
    display_eq(s, est)
