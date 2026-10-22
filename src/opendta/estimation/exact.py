"""exlogistic and expoisson ([R] exlogistic, [R] expoisson): exact conditional inference on the
distribution of the sufficient statistics (core/exact.py).

Layouts come from the manuals; the lines marked VERIFY wait for the Stata 14 log of case 1908."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ..core import exact as X
from ..core import missing as M
from ..core.errors import StataError
from . import fvars
from .models2 import _header, g
from .regress import Sample, cmdline, comma
from .results import CoefRow, Estimates, post

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def cmd_exlogistic(s: "Session", args: str) -> None:
    _exact(s, args, "exlogistic")


def cmd_expoisson(s: "Session", args: str) -> None:
    _exact(s, args, "expoisson")


def _exact(s: "Session", args: str, cmd: str) -> None:
    from .postest import replay
    if not args.strip() or args.strip().startswith(","):
        replay(s, cmd, args)
        return
    logistic = cmd == "exlogistic"
    smp = Sample(s, args, opts_spec={"test": 2, "condvars": 4, "coef": 4, "midp": 4, "nolog": 4, "memory": 3,
                                     "saving": 3, "estconstant": 4, "terms": 5, "group": 3, "binomial": 3,
                                     "exposure": 3, "offset": 3, "irr": 3, "or": 2},
                 weights=("fweight",))
    o = smp.o
    for bad in ("terms", "group", "binomial", "exposure", "offset", "condvars", "saving", "estconstant"):
        if o.get(bad):
            raise StataError(198, f"option {bad}() is not available yet in {cmd}")   # VERIFY
    test = str(o.get("test", "sufficient")).strip().lower()   # VERIFY: probability keeps the sufficient layout
    if test not in ("score", "sufficient", "probability", "s", "su", "p", "pr", "sc"):
        raise StataError(198, "option test() invalid")   # VERIFY
    est_cols = [c for c in smp.cols if c.values is not None]
    if not est_cols:
        raise StataError(102, "too few variables specified")
    y = smp.y.astype(float)
    n = smp.n
    if logistic:
        if not np.all((y == 0) | (y == 1)):
            raise StataError(450, f"{smp.depname} must be 0 or 1")   # VERIFY
        if np.all(y == y[0]):
            raise StataError(2000, "outcome does not vary")
    else:
        if np.any(y < 0) or np.any(y != np.floor(y)):
            raise StataError(459, f"{smp.depname} must be a nonnegative integer")   # VERIFY
    Xm = np.column_stack([np.ones(n)] + [c.values for c in est_cols])
    obs = tuple(round(float(v), 9) for v in (Xm * y[:, None]).sum(axis=0))
    if logistic:
        dist = X.suff_dist_logistic(Xm)
        pbar = float(y.mean())
    else:
        total = int(y.sum())
        try:
            dist = X.suff_dist_poisson(Xm, total)
        except ValueError:
            raise StataError(459, "the covariates of expoisson must be integers")   # VERIFY
        except MemoryError:
            raise StataError(950, "insufficient memory")
        pbar = float(y.mean())
    sc, p_sc = X.score_test(dist, obs, pbar, Xm, "logistic" if logistic else "poisson")
    if logistic:
        sc *= (n - 1) / n
    mid = bool(o.get("midp"))
    alpha = 1 - smp.level / 100
    rows = []
    for j, c in enumerate(est_cols, start=1):
        sup, wts = X.conditional(dist, obs, j)
        t = obs[j]
        b = X.cmle(sup, wts, t)
        is_mue = b is None
        if is_mue:
            b = X.mue(sup, wts, t, mid)
        lo, hi = X.exact_ci(sup, wts, t, alpha, mid)
        if test in ("score", "sc"):
            tstat, pv = X.score_conditional(sup, wts, t)
        else:
            tstat, pv = t, X.exact_p(sup, wts, t, mid)
        rows.append((c.name, b, is_mue, tstat, pv, lo, hi))
    names = [r[0] for r in rows]
    b_vec = np.array([r[1] for r in rows])
    est = Estimates(cmd, smp.depname, names, b_vec, np.zeros((len(names), len(names))), n, "z", None)
    est.level = smp.level
    est.vcetype, est.clustvar, est.n_clust = "", "", None
    est.extra = {"N": n, "score": sc, "p_score": p_sc, "rows": rows, "coef": bool(o.get("coef")) if logistic else not o.get("irr"), "score_cols": test in ("score", "sc"),
                 "mue": any(r[2] for r in rows)}
    est.display = _display_exact
    est.rows = []
    title = "Exact logistic regression" if logistic else "Exact Poisson regression"
    est.title = title
    scalars = [("N", n), ("sum_y", float(y.sum())), ("n_possible", float(len(dist))), ("p_score", p_sc),
               ("score", sc)]
    macros = [("cmdline", cmdline(cmd, args)), ("title", title), ("depvar", smp.depname), ("cmd", cmd),
              ("properties", "b"), ("level", f"{smp.level:g}")]
    post(s, est, scalars, macros, smp.mask)
    _display_exact(s, est)


def _display_exact(s, est, *, header=True, table=True, level=None, **_kw) -> None:
    x = est.extra
    out = s.output
    lv = f"{est.level:g}"
    ratio = "Odds Ratio" if est.cmd == "exlogistic" else "IRR"
    first = "Coef." if x["coef"] else ratio
    if header:
        out.write("\n", "text")
        if est.cmd == "exlogistic":
            rows_ = [(est.title, "Number of obs", comma(x["N"])), ("", "Model score", g(x["score"], "%9.0g")),
                     ("", "Pr >= score", g(x["p_score"], "%10.4f"))]
            labw = 14
        else:
            rows_ = [(est.title, None, None), ("", "Number of obs", comma(x["N"]))]
            labw = 14
        for lft, lb, val in rows_:
            if lb is None:
                out.write(lft + "\n", "text")
                continue
            out.write(f"{lft:<49}{lb:<{labw}}= ", "text")
            out.write(f"{val:>10}\n", "result")
    t1, t2 = ("Score", "Pr>=Score") if x["score_cols"] else ("Suff.", "2*Pr(Suff.)")
    out.write("-" * 75 + "\n", "text")
    out.write(f"{est.depvar:>12} |{first:>11}{t1:>12}{t2:>13}{f'[{lv}% Conf. Interval]':>25}\n", "text")
    out.write("-" * 13 + "+" + "-" * 61 + "\n", "text")
    conv = (lambda v: v) if x["coef"] else (lambda v: float(np.exp(v)) if np.isfinite(v) else v)
    for nm, b, is_mue, t, p, lo, hi in x["rows"]:
        def fmt(v):
            if v == np.inf:
                return "+Inf"
            if v == -np.inf:
                return "-Inf"
            return g(conv(v), "%9.0g")
        star = "*" if is_mue else " "
        bs = g(conv(b), '%9.0g') + ("*" if star == "*" else "")
        out.write(f"{nm:>12} |{bs:>11}{g(t, '%9.0g'):>12}{g(p, '%10.4f'):>12}"
                  f"{fmt(lo):>14}{fmt(hi):>12}\n", "text")
    out.write("-" * 75 + "\n", "text")
    if x["mue"]:
        out.write("(*) median unbiased estimates (MUE)\n", "text")
