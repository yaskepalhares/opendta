"""Receiver operating characteristic analysis ([R] roc, [R] rocfit, [R] lroc,
[R] lsens): roctab (nonparametric ROC curve of one classifier), roccomp
(comparison of the areas of several classifiers), rocfit (binormal model fit by
maximum likelihood), and lroc/lsens after logit, logistic or probit.

The area is the Mann-Whitney statistic; its standard error is DeLong, DeLong and
Clarke-Pearson's (1988) by default, or Hanley and McNeil's (1982). roccomp tests
the equality of areas with DeLong's covariance (correlated classifiers) or with
independent variances (by()). Layout details are marked VERIFY and checked in
compat/do/1402_kappa_roc.do."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core import roc as R
from ..core.errors import StataError
from ..lang.syntax import match_options, parse_standard
from ..lazy import lazy
from ._util import touse
from .epitab import _set_r, _strata
from .registry import command
from .summarize import _level, name12, weights

st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

NAN = float("nan")
SYS = M.SYSMISS


# ---------------------------------------------------------------------------
# area and its variance
# ---------------------------------------------------------------------------

def placements(pos: np.ndarray, neg: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """V10 (one per positive) and V01 (one per negative): the share of the other group a
    score beats, ties counting one half."""
    ns = np.sort(neg)
    ps = np.sort(pos)
    v10 = (np.searchsorted(ns, pos, "left") + np.searchsorted(ns, pos, "right")) / (2 * len(neg))
    v01 = 1 - (np.searchsorted(ps, neg, "left") + np.searchsorted(ps, neg, "right")) / (2 * len(pos))
    return v10, v01


def auc_delong(pos: np.ndarray, neg: np.ndarray) -> tuple[float, float]:
    v10, v01 = placements(pos, neg)
    a = float(v10.mean())
    m, n = len(pos), len(neg)
    var = (v10.var(ddof=1) if m > 1 else 0.0) / m + (v01.var(ddof=1) if n > 1 else 0.0) / n
    return a, math.sqrt(var)


def se_hanley(a: float, m: int, n: int) -> float:
    q1, q2 = a / (2 - a), 2 * a * a / (1 + a)
    var = (a * (1 - a) + (m - 1) * (q1 - a * a) + (n - 1) * (q2 - a * a)) / (m * n)
    return math.sqrt(max(var, 0.0))


def delong_cov(pos_list, neg_list) -> tuple[np.ndarray, np.ndarray]:
    """Areas and DeLong covariance matrix of classifiers scored on the same subjects."""
    V10 = []
    V01 = []
    for pos, neg in zip(pos_list, neg_list):
        v10, v01 = placements(pos, neg)
        V10.append(v10)
        V01.append(v01)
    V10, V01 = np.array(V10), np.array(V01)
    m, n = V10.shape[1], V01.shape[1]
    areas = V10.mean(1)
    S = np.cov(V10) / m + np.cov(V01) / n
    return areas, np.atleast_2d(S)


def roc_points(pos: np.ndarray, neg: np.ndarray):
    """Cutpoints (each distinct score, classified + if score >= cutpoint) with sensitivity
    and specificity, plus the final 'greater than the maximum' row."""
    cuts = np.unique(np.concatenate([pos, neg]))
    rows = []
    for c in cuts:
        rows.append((c, float((pos >= c).mean()), float((neg < c).mean()), ">="))
    rows.append((cuts[-1], 0.0, 1.0, ">"))
    return rows


def _expand(x: np.ndarray, w: np.ndarray | None) -> np.ndarray:
    return x if w is None else np.repeat(x, w.astype(int))


def _ci(a: float, se: float, level: float) -> tuple[float, float]:
    z = float(st.norm.ppf(1 - (1 - level / 100) / 2))
    return a - z * se, a + z * se


def _f(x, spec) -> str:
    return "." if x is None or (isinstance(x, float) and math.isnan(x)) else format(x, spec)


# ---------------------------------------------------------------------------
# roctab
# ---------------------------------------------------------------------------

def _pos_neg(s: "Session", p, ref: str, cls: list[str], extra_mask=None):
    ds = s.data
    rv = ds.get(s.expand_varlist(ref)[0])
    cv = [ds.get(s.expand_varlist(c)[0]) for c in cls]
    mask = touse(s, p) & (rv.data < SYS)
    for c in cv:
        if c.is_string:
            raise StataError(109, "type mismatch")
        mask &= c.data < SYS
    if extra_mask is not None:
        mask &= extra_mask
    w, _, mask = weights(s, p, mask, ("fweight",))
    vals = set(rv.data[mask].tolist())
    if not vals <= {0.0, 1.0}:
        # Stata: the reference variable must be 0/1 (VERIFY: rc and message)
        raise StataError(450, f"{rv.name} not a 0/1 variable")
    return rv, cv, mask, w


@command("roctab", "roctab")
def cmd_roctab(s: "Session", args: str) -> None:
    p = parse_standard(args)
    names = p.varlist.split()
    if len(names) != 2:
        raise StataError(198, "invalid syntax")
    o = match_options(p.options, {"detail": 1, "lorenz": 3, "binomial": 3, "nolabel": 5, "table": 3,
                                  "bamber": 3, "hanley": 3, "graph": 2, "norefline": 5, "summary": 3,
                                  "specificity": 4, "level": 1})
    level = _level(s, o)
    rv, (cv,), mask, w = _pos_neg(s, p, names[0], [names[1]])
    wm = None if w is None else w[mask]
    y = _expand(rv.data[mask], wm)
    x = _expand(cv.data[mask], wm)
    pos, neg = x[y == 1], x[y == 0]
    if not len(pos) or not len(neg):
        raise StataError(2000, "no observations")   # VERIFY
    _a, se = auc_delong(pos, neg)
    a = R.area_trapezoid(pos, neg)
    if o.get("hanley") and o.get("bamber"):
        raise StataError(198, "options hanley and bamber may not be combined")   # VERIFY
    if o.get("hanley"):
        se = se_hanley(a, len(pos), len(neg))
    elif o.get("bamber"):
        se = R.se_bamber(pos, neg, a)
    N = len(pos) + len(neg)
    lo, hi = R.binomial_ci(a, N, level) if o.get("binomial") else _ci(a, se, level)
    out = s.output
    rows = roc_points(pos, neg)
    if o.get("table"):
        # the two-way table of tabulate (compat 1402)
        from .tabulate import _twoway
        _twoway(s, rv, cv, mask, w, "fweight" if w is not None else "", {})
    if o.get("detail"):
        out.write("\nDetailed report of sensitivity and specificity\n", "text")
        out.write("-" * 78 + "\n", "text")
        out.write(" " * 43 + "Correctly\n", "text")
        out.write("Cutpoint      Sensitivity   Specificity   Classified          LR+          LR-\n", "text")
        out.write("-" * 78 + "\n", "text")
        N = len(pos) + len(neg)
        for c, se_, sp, op in rows:
            cc = (se_ * len(pos) + sp * len(neg)) / N
            lrp = se_ / (1 - sp) if sp < 1 else NAN
            lrn = (1 - se_) / sp if sp > 0 else NAN
            cut = f"( {op:<2} {c:g} )"
            out.write(f"{cut:<14}", "text")
            out.write(f"{100 * se_:>10.2f}%{100 * sp:>13.2f}%{100 * cc:>12.2f}%"
                      f"{_f(lrp, '>13.4f') if not math.isnan(lrp) else ' ' * 13}"
                      f"{_f(lrn, '>13.4f') if not math.isnan(lrn) else ''}\n", "result")
        out.write("-" * 78 + "\n", "text")
    show_area = not o.get("graph") or o.get("summary") or o.get("detail") or o.get("table")
    if o.get("lorenz") and not o.get("summary"):
        show_area = False   # lorenz reports the indices instead (compat 1601)
    if show_area:
        # the standard error's method above its column, the interval's above the limits
        # (VERIFY: exact placement and dashes)
        se_label = "Hanley" if o.get("hanley") else "Bamber" if o.get("bamber") else ""
        title = "-- Binomial Exact --" if o.get("binomial") else "-Asymptotic Normal--"
        # ROC, 6 spaces, the standard error's method in 6, 8 spaces, the interval (compat 1402)
        out.write("\n" + " " * 22 + "ROC" + " " * 6 + f"{se_label:<6}" + " " * 8 + title + "\n", "text")
        out.write(f"{'Obs':>14}{'Area':>11}     Std. Err.      [{level:g}% Conf. Interval]\n", "text")
        out.write(" " * 5 + "-" * 60 + "\n", "text")
        out.write(f"{len(pos) + len(neg):>14,}{a:>11.4f}{se:>13.4f}{lo:>15.5f}{hi:>12.5f}\n", "result")
    if o.get("lorenz"):
        lx, ly, gini, pietra = R.lorenz(pos, neg)
        out.write("\n Lorenz curve \n" + "-" * 27 + "\n", "text")
        out.write("  Pietra index =  ", "text")
        out.write(f"{pietra:>7.4f}\n", "result")
        out.write("  Gini index   =  ", "text")
        out.write(f"{gini:>7.4f}\n", "result")
    if o.get("graph"):
        if o.get("lorenz"):
            _lorenz_graph(s, rv.name, lx, ly, refline=not o.get("norefline"))
        else:
            _roc_graph(s, [(cv.name, rows, a)], refline=not o.get("norefline"),
                       title=f"Area under ROC curve = {a:.4f}")
    if o.get("lorenz") and not show_area:
        _set_r(s, [("gini", gini), ("pietra", pietra)])   # only these (compat 1601)
    else:
        _set_r(s, [("ub", hi), ("lb", lo), ("se", se), ("area", a), ("N", float(N))])


def _lorenz_graph(s: "Session", refname: str, x, y, refline: bool) -> None:
    from ..graphics.model import Graph, Layer, Panel
    from .graphcmd import finish
    layers = [Layer("connected", list(map(float, x)), list(map(float, y)), {}, {}, "Lorenz curve")]
    if refline:
        layers.append(Layer("line", [0.0, 1.0], [0.0, 1.0], {}, {"lcolor": "gs8"}, "Reference"))
    opts = {"_xtitle": f"cumulative % of {refname}=0", "_ytitle": f"cumulative % of {refname}=1",
            "aspectratio": 1.0, "legend": {"off": True}}
    finish(s, Graph(family="twoway", opts=opts, panels=[Panel("", layers)]), {})


def _roc_graph(s: "Session", curves, refline: bool, title: str = "", legend: bool = False) -> None:
    from ..graphics.model import Graph, Layer, Panel
    from .graphcmd import finish
    layers = []
    for name, rows, a in curves:
        pts = sorted({(1 - sp, se_) for _, se_, sp, _ in rows} | {(0.0, 0.0), (1.0, 1.0)})
        label = f"{name} ROC area: {a:.4f}" if legend else name
        layers.append(Layer("connected", [q[0] for q in pts], [q[1] for q in pts], {}, {}, label))
    if refline:
        layers.append(Layer("line", [0.0, 1.0], [0.0, 1.0], {}, {"lcolor": "gs8"}, "Reference"))
    opts = {"_xtitle": "1 - Specificity", "_ytitle": "Sensitivity", "aspectratio": 1.0}
    if title:
        opts["note"] = {"lines": [title]}
    if not legend:
        opts["legend"] = {"off": True}
    finish(s, Graph(family="twoway", opts=opts, panels=[Panel("", layers)]), {})


# ---------------------------------------------------------------------------
# roccomp
# ---------------------------------------------------------------------------

@command("roccomp", "roccomp")
def cmd_roccomp(s: "Session", args: str) -> None:
    p = parse_standard(args)
    names = p.varlist.split()
    if len(names) < 2:
        raise StataError(198, "invalid syntax")
    o = match_options(p.options, {"by": 2, "test": 1, "graph": 2, "norefline": 5, "separate": 3,
                                  "summary": 3, "binormal": 3, "line": 4, "level": 1})
    level = _level(s, o)
    binorm = bool(o.get("binormal"))
    out = s.output
    if isinstance(o.get("by"), str):
        if len(names) != 2:
            raise StataError(198, "only one classifier allowed with by()")   # VERIFY
        rv, (cv,), mask, w = _pos_neg(s, p, names[0], [names[1]])
        byname, groups = _strata(s, str(o["by"]), mask)
        items = []
        for lab, sel in groups:
            if binorm:
                out.write(f"Fitting binormal model for: {byname}=={lab}\n", "text")
            ws = None if w is None else w[sel]
            y, x = _expand(rv.data[sel], ws), _expand(cv.data[sel], ws)
            pos, neg = x[y == 1], x[y == 0]
            a, se = _binormal_area(pos, neg) if binorm else auc_delong(pos, neg)
            if not binorm:
                a = R.area_trapezoid(pos, neg)
            items.append((lab, len(pos) + len(neg), a, se, roc_points(pos, neg)))
        areas = np.array([i[2] for i in items])
        S = np.diag([i[3] ** 2 for i in items])
    else:
        rv, cvs, mask, w = _pos_neg(s, p, names[0], names[1:])
        wm = None if w is None else w[mask]
        y = _expand(rv.data[mask], wm)
        xs = [_expand(c.data[mask], wm) for c in cvs]
        pos_l = [x[y == 1] for x in xs]
        neg_l = [x[y == 0] for x in xs]
        areas, S = delong_cov(pos_l, neg_l)
        areas = np.array([R.area_trapezoid(pp, nn) for pp, nn in zip(pos_l, neg_l)])
        ses = np.sqrt(np.diag(S))
        if binorm:
            # binormal areas and standard errors; the test uses DeLong's covariance of the
            # correlated classifiers (VERIFY)
            fits = []
            for c, pp, nn in zip(cvs, pos_l, neg_l):
                out.write(f"Fitting binormal model for: {c.name}\n", "text")
                fits.append(_binormal_area(pp, nn))
            areas = np.array([f[0] for f in fits])
            ses = np.array([f[1] for f in fits])
        items = [(c.name, len(y), areas[k], ses[k], roc_points(pos_l[k], neg_l[k]))
                 for k, c in enumerate(cvs)]
    # binormal fits show no interval title (compat 1601)
    out.write("\n" + " " * 30 + "ROC" + ("" if binorm else " " * 20 + "-Asymptotic Normal--") + "\n", "text")
    byname = name12(str(o["by"]).split()[0]) if isinstance(o.get("by"), str) else ""
    out.write(f"{byname:<17}{'Obs':>5}{'Area':>11}     Std. Err.      [{level:g}% Conf. Interval]\n", "text")
    out.write("-" * 73 + "\n", "text")
    for lab, n, a, se, _ in items:
        lo, hi = _ci(a, se, level)
        out.write(f"{name12(str(lab)):<17}", "text")
        out.write(f"{n:>5}{a:>11.4f}{se:>13.4f}{lo:>15.5f}{hi:>12.5f}\n", "result")
    out.write("-" * 73 + "\n", "text")
    k = len(items)
    L = np.zeros((k - 1, k))
    for i in range(k - 1):
        L[i, i], L[i, i + 1] = 1, -1
    if isinstance(o.get("test"), str):
        m = getattr(s, "matrices", {}).get(str(o["test"]))
        if m is None:
            raise StataError(111, f"matrix {o['test']} not found")
        L = np.asarray(m.data, float)
    d = L @ areas
    try:
        chi2 = float(d @ np.linalg.solve(L @ S @ L.T, d))
    except np.linalg.LinAlgError:
        chi2 = NAN
    df = L.shape[0]
    pv = float(st.chi2.sf(chi2, df)) if not math.isnan(chi2) else NAN
    labs = [str(i[0]) for i in items]
    ho = " = ".join(f"area({lab})" for lab in labs) if not isinstance(o.get("test"), str) else "Comparison as defined by contrast matrix: " + str(o["test"])
    out.write(f"Ho: {ho}\n", "text")
    out.write(f"    chi2({df}) = ", "text")
    out.write(f"{_f(chi2, '>8.2f')}", "result")
    out.write("       Prob>chi2 = ", "text")
    out.write(f"{_f(pv, '>8.4f')}\n", "result")
    if o.get("graph"):
        _roc_graph(s, [(str(i[0]), i[4], i[2]) for i in items], refline=not o.get("norefline"), legend=True)
    from .matrix import Matrix
    # return list order of Stata 14 (compat 1402)
    s.r = {"p": pv if not math.isnan(pv) else SYS, "df": float(df), "chi2": chi2 if not math.isnan(chi2) else SYS,
           "N_g": float(k), "V": Matrix(S, labs, labs)}


def _binormal_area(pos: np.ndarray, neg: np.ndarray) -> tuple[float, float]:
    """Area under the binormal ROC curve fit to the categories of the scores, and its standard
    error ([R] rocfit)."""
    from ..estimation import binormal
    cats = np.unique(np.concatenate([pos, neg]))
    if len(cats) < 3:
        raise StataError(459, "classvar must have at least 3 categories for the binormal fit")   # VERIFY
    n0 = np.bincount(np.searchsorted(cats, neg), minlength=len(cats)).astype(float)
    n1 = np.bincount(np.searchsorted(cats, pos), minlength=len(cats)).astype(float)
    bf = binormal.fit(n0, n1)
    est, se = binormal.indices(bf.b[0], bf.b[1], bf.V[:2, :2])
    return float(est[0]), float(se[0])


# ---------------------------------------------------------------------------
# lroc, lsens (after logit, logistic, probit)
# ---------------------------------------------------------------------------

def _fitted(s: "Session", varlist: str):
    from ..estimation.estat import _binary
    from ..estimation.results import current
    est = current(s)
    y, pr, mask = _binary(s, est)
    if varlist.strip():
        dv = s.data.get(s.expand_varlist(varlist)[0])
        y = dv.data[mask] != 0
    return est, y, pr, mask


@command("lroc", "lroc")
def cmd_lroc(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"nograph": 3, "all": 3, "beta": 4})
    for k in ("all", "beta"):
        if o.get(k):
            raise StataError(198, f"option {k} of lroc is not available yet in OpenDTA")
    est, y, pr, mask = _fitted(s, p.varlist)
    pos, neg = pr[y], pr[~y]
    a = R.area_trapezoid(pos, neg)   # VERIFY for lroc
    title = ("Probit" if est.cmd == "probit" else "Logistic") + f" model for {est.depvar}"
    out = s.output
    out.write(f"\n{title}\n\n", "text")
    out.write("number of observations = ", "text")
    out.write(f"{len(y):>8}\n", "result")
    out.write("area under ROC curve   =", "text")
    out.write(f"{a:>9.4f}\n", "result")
    if not o.get("nograph"):
        _roc_graph(s, [("", roc_points(pos, neg), a)], refline=True, title=f"Area under ROC curve = {a:.4f}")
    _set_r(s, [("N", float(len(y))), ("area", a)])


@command("lsens", "lsens")
def cmd_lsens(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"all": 3, "genprob": 4, "gensens": 4, "genspec": 4, "replace": 3,
                                  "nograph": 3})
    if o.get("all"):
        raise StataError(198, "option all of lsens is not available yet in OpenDTA")
    est, y, pr, mask = _fitted(s, p.varlist)
    pos, neg = pr[y], pr[~y]
    cuts = np.unique(pr)
    sens = np.array([(pos >= c).mean() for c in cuts])
    spec = np.array([(neg < c).mean() for c in cuts])
    ds = s.data
    # the variables hold one row per observed cutoff; the graph adds the ends 0 and 1
    new = [(o.get("genprob"), cuts), (o.get("gensens"), sens), (o.get("genspec"), spec)]
    cuts = np.concatenate([[0.0], cuts, [1.0]])
    sens = np.concatenate([[1.0], sens, [0.0]])
    spec = np.concatenate([[0.0], spec, [1.0]])
    for name, vals in new:
        if isinstance(name, str):
            # one row per cutoff, in the first observations (VERIFY: Stata's layout)
            if ds.has(name):
                if not o.get("replace"):
                    raise StataError(110, f"{name} already defined")
                ds.drop_vars([name])
            if len(vals) > ds.nobs:
                raise StataError(198, "too few observations to store the cutoffs")   # VERIFY
            col = np.full(ds.nobs, SYS)
            col[:len(vals)] = vals
            from ..core.dataset import Variable
            ds.add(Variable(name, "float", col))
    if any(isinstance(n, str) for n, _ in new):
        s.notify_state()
    if not o.get("nograph"):
        from ..graphics.model import Graph, Layer, Panel
        from .graphcmd import finish
        layers = [Layer("connected", cuts.tolist(), sens.tolist(), {}, {}, "Sensitivity"),
                  Layer("connected", cuts.tolist(), spec.tolist(), {}, {}, "Specificity")]
        opts = {"_xtitle": "Probability cutoff", "_ytitle": "Sensitivity/Specificity"}
        finish(s, Graph(family="twoway", opts=opts, panels=[Panel("", layers)]), {})


# ---------------------------------------------------------------------------
# rocfit (binormal model)
# ---------------------------------------------------------------------------

@command("rocfit", "rocfit")
def cmd_rocfit(s: "Session", args: str) -> None:
    from ..estimation import binormal, ml
    from ..estimation.results import Estimates, post
    p = parse_standard(args)
    names = p.varlist.split()
    if len(names) != 2:
        raise StataError(198, "invalid syntax")
    o = match_options(p.options, {"continuous": 4, "generate": 3, "level": 1, "iterate": 4,
                                  "tolerance": 3, "ltolerance": 4, "nolog": 3, "log": 3})
    level = _level(s, o)
    rv, (cv,), mask, w = _pos_neg(s, p, names[0], [names[1]])
    x = cv.data[mask]
    if isinstance(o.get("continuous"), str):
        G = int(o["continuous"])
        qs = np.quantile(x, np.linspace(0, 1, G + 1)[1:-1])
        x = np.searchsorted(qs, x, side="right").astype(float)
    cats = np.unique(x)
    K = len(cats)
    if K < 3:
        raise StataError(459, "classvar must have at least 3 categories for the binormal fit")   # VERIFY
    wv = np.ones(mask.sum()) if w is None else w[mask]
    k_idx = np.searchsorted(cats, x)
    y = rv.data[mask] == 1
    n0 = np.bincount(k_idx[~y], weights=wv[~y], minlength=K)
    n1 = np.bincount(k_idx[y], weights=wv[y], minlength=K)
    out = s.output
    bf = binormal.fit(n0, n1, maxiter=int(o["iterate"]) if isinstance(o.get("iterate"), str) else 16000)
    res = bf.res
    if not o.get("nolog"):
        out.write("\nFitting binormal model:\n\n", "text")
        ml.print_log(out, res)
    b, V, p0, p1 = bf.b, bf.V, bf.p0, bf.p1
    # Pearson goodness of fit over the 2 x K table (df = K - 3; VERIFY)
    e0, e1 = n0.sum() * p0, n1.sum() * p1
    gof = float(((n0 - e0) ** 2 / e0).sum() + ((n1 - e1) ** 2 / e1).sum())
    df_gof = K - 3
    p_gof = float(st.chi2.sf(gof, df_gof)) if df_gof > 0 else NAN
    N = int(n0.sum() + n1.sum())
    est_i, se_i = binormal.indices(b[0], b[1], V[:2, :2])
    names_b = ["intercept", "slope"] + [f"cut{i + 1}" for i in range(K - 1)]
    est = Estimates("rocfit", rv.name, ["_cons"] * len(names_b), b.copy(), V, N, stat="z", eqnames=names_b)
    est.level = level
    est.extra = {"title": f"Binormal model of {rv.name} on {cv.name}", "gof": gof, "df_gof": df_gof,
                 "p_gof": p_gof, "ll": res.ll, "idx": (est_i, se_i)}
    est.display = _rocfit_display
    # e() in Stata 14's ereturn list order (compat 1402)
    scalars = [("rank", float(len(b))), ("N", N), ("ic", float(res.iterations)), ("k", float(len(b))),
               ("k_eq", float(len(b))), ("k_dv", 2.0), ("converged", 1.0 if res.converged else 0.0),
               ("rc", 0.0), ("ll", res.ll), ("k_eq_model", 1.0), ("df_m", 0.0), ("chi2_gf", gof),
               ("df_gf", df_gof), ("p_gf", p_gof), ("area", est_i[0]), ("se_area", se_i[0]),
               ("deltam", est_i[1]), ("se_delm", se_i[1]), ("de", est_i[2]), ("se_da", se_i[3]),
               ("da", est_i[3]), ("se_de", se_i[2])]
    wt = [("wexp", f"= {p.weight[1]}")] if p.weight else []
    wtt = [("wtype", "fweight")] if p.weight else []
    macros = ([("cmd", "rocfit"), ("cmdline", f"rocfit {args.strip()}"), ("depvar", f"{rv.name} {cv.name}"),
               ("title", est.extra["title"])] + wt + [("chi2type", "GOF")] + wtt
              + [("opt", "moptimize"), ("vce", "oim"), ("user", "rocf_lf"), ("ml_method", "lf"),
                 ("technique", "nr"), ("which", "max"), ("properties", "b V")])
    from .matrix import Matrix
    grad = Matrix(res.g.reshape(1, -1), ["r1"], ["_cons"] * len(b))
    grad.coleq = names_b
    post(s, est, scalars, macros, mask, extra_mats=[("ilog", Matrix(ml.ilog(res), ["r1"],
                                                                     [f"c{i}" for i in range(1, 21)])),
                                                   ("gradient", grad)])
    _rocfit_display(s, est)


def _rocfit_display(s: "Session", est, *, header: bool = True, table: bool = True, level=None, **_kw) -> None:
    """rocfit's report: intercept and slope (with the z test that the slope is 1), the
    cutpoints, then the indices of the fitted curve; numbers as %9.6f (compat 1402)."""
    out = s.output
    x = est.extra
    lv = level or est.level
    z = float(st.norm.ppf(1 - (1 - lv / 100) / 2))
    if header:
        out.write(f"\n{x['title']:<50}Number of obs   = {est.N:>10}\n", "text")
        out.write(f"Goodness-of-fit chi2({x['df_gof']}) = ", "text")
        out.write(f"{x['gof']:>11.2f}\n", "result")
        out.write("Prob > chi2             = ", "text")
        out.write(f"{_f(x['p_gof'], '>11.4f')}\n", "result")
        out.write("Log likelihood          = ", "text")
        out.write(f"{_g10(x['ll']):>11}\n", "result")
    if not table:
        return
    from ..estimation.results import pvalue
    b, V = est.b, est.V
    se = np.sqrt(np.diag(V))

    def row(label, bb, ss, zz):
        out.write(f"{label:>12} |", "text")
        out.write(f"{bb:>11.6f}{ss:>11.6f}{zz:>9.2f}{pvalue('z', zz, None):>8.3f}{bb - z * ss:>13.6f}"
                  f"{bb + z * ss:>12.6f}\n", "result")
    out.write("\n" + "-" * 78 + "\n", "text")
    out.write(f"{'':>12} |      Coef.   Std. Err.      z    P>|z|     [{lv:g}% Conf. Interval]\n", "text")
    out.write("-" * 13 + "+" + "-" * 64 + "\n", "text")
    row("intercept", b[0], se[0], b[0] / se[0])
    row("slope (*)", b[1], se[1], (b[1] - 1) / se[1])
    out.write("-" * 13 + "+" + "-" * 64 + "\n", "text")
    for i in range(2, len(b)):
        row(f"/cut{i - 1}", b[i], se[i], b[i] / se[i])
    out.write("-" * 78 + "\n", "text")
    est_i, se_i = x["idx"]
    out.write("-" * 78 + "\n", "text")
    out.write(f"{'':>12} |                    Indices from binormal fit\n", "text")
    out.write(f"{'Index':>12} |   Estimate   Std. Err.                     [{lv:g}% Conf. Interval]\n", "text")
    out.write("-" * 13 + "+" + "-" * 64 + "\n", "text")
    for lab, e_, s_ in zip(("ROC area", "delta(m)", "d(e)", "d(a)"), est_i, se_i):
        out.write(f"{lab:>12} |", "text")
        out.write(f"{e_:>11.6f}{s_:>11.6f}{'':>17}{e_ - z * s_:>13.6f}{e_ + z * s_:>12.6f}\n", "result")
    out.write("-" * 78 + "\n", "text")
    out.write("(*) z test for slope==1\n", "text")


def _g10(x: float) -> str:
    from ..core.formats import format_value
    return format_value(float(x), "%10.0g", pad=False).strip()
