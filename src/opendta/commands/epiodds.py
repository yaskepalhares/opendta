"""Odds by exposure level and Mantel-Haenszel odds ratios ([ST] epitab: tabodds,
mhodds) and tests of table symmetry ([ST] symmetry). Layout details not shown
in the manual are marked VERIFY and checked in compat/do/1401_epitab.do."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from ..core import epi as E
from ..core import missing as M
from ..core.errors import StataError
from ..lang.syntax import match_options, parse_standard
from ..lazy import lazy
from ._util import touse
from .epitab import _p4, _r, _set_r, _strata
from .registry import command
from .summarize import _level, name12, weights

st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

NAN = float("nan")
SYS = M.SYSMISS


def _prep(s: "Session", p, case: str, extra: list[str]):
    ds = s.data
    names = [s.expand_varlist(case)[0]] + [s.expand_varlist(n)[0] for n in extra]
    mask = touse(s, p)
    for n in names:
        v = ds.get(n)
        mask &= (np.array([x != "" for x in v.raw]) if v.is_string else v.data < SYS)
    w, _, mask = weights(s, p, mask, ("fweight",))
    if w is None:
        w = np.ones(ds.nobs)
    return mask, w


def _levels(s: "Session", var, mask):
    from .inspect import cell_text
    vals = sorted(set(var.data[mask].tolist()))
    labels = []
    for v in vals:
        i = int(np.flatnonzero(mask & (var.data == v))[0])
        labels.append(cell_text(s.data, var, i))
    return vals, labels


def _trend(strata_tables, scores) -> tuple[float, float]:
    """Score (Cochran-Armitage / Mantel extension) statistic for a linear trend in the log
    odds, summed over strata: U = sum x_k (a_k - E_k), V its hypergeometric variance."""
    U = V = 0.0
    x = np.asarray(scores, dtype=float)
    for cases, ctrls in strata_tables:
        a = np.asarray(cases, float)
        n_k = a + np.asarray(ctrls, float)
        N, A = n_k.sum(), a.sum()
        if N < 2:
            continue
        U += float((x * (a - n_k * A / N)).sum())
        V += A * (N - A) / (N * N * (N - 1)) * (N * float((n_k * x * x).sum()) - float((n_k * x).sum()) ** 2)
    return U, V


def _homogeneity(strata_tables) -> tuple[float, int]:
    """Generalized Cochran-Mantel-Haenszel test that the odds are equal across the K levels
    (K-1 df), summed over strata; with one stratum it is N-1/N times Pearson's chi2
    (VERIFY: Stata may report Pearson's chi2 without adjust())."""
    K = len(strata_tables[0][0])
    d = np.zeros(K - 1)
    Vm = np.zeros((K - 1, K - 1))
    for cases, ctrls in strata_tables:
        a = np.asarray(cases, float)
        n_k = a + np.asarray(ctrls, float)
        N, A = n_k.sum(), a.sum()
        if N < 2:
            continue
        d += (a - n_k * A / N)[:-1]
        f = A * (N - A) / (N * N * (N - 1))
        nk = n_k[:-1]
        Vm += f * (N * np.diag(nk) - np.outer(nk, nk))
    keep = np.diag(Vm) > 0
    d, Vm = d[keep], Vm[np.ix_(keep, keep)]
    if not len(d):
        return NAN, 0
    chi2 = float(d @ np.linalg.solve(Vm, d))
    return chi2, len(d)


# ---------------------------------------------------------------------------
# tabodds
# ---------------------------------------------------------------------------

@command("tabodds", "tabodds")
def cmd_tabodds(s: "Session", args: str) -> None:
    p = parse_standard(args)
    names = p.varlist.split()
    if not 1 <= len(names) <= 2:
        raise StataError(198, "invalid syntax")
    o = match_options(p.options, {"or": 2, "adjust": 3, "base": 2, "cornfield": 4, "tb": 2,
                                  "woolf": 1, "graph": 2, "ciplot": 2, "level": 1})
    level = _level(s, o)
    adj = s.expand_varlist(str(o["adjust"])) if isinstance(o.get("adjust"), str) else []
    expn = names[1] if len(names) == 2 else None
    mask, w = _prep(s, p, names[0], ([expn] if expn else []) + adj)
    ds = s.data
    case = ds.get(s.expand_varlist(names[0])[0]).data != 0
    if expn is None:
        raise StataError(198, "expvar required")   # VERIFY: tabodds casevar alone (odds of the sample)
    ev = ds.get(s.expand_varlist(expn)[0])
    vals, labels = _levels(s, ev, mask)
    K = len(vals)
    strata = [("", mask)] if not adj else _strata(s, " ".join(adj), mask)[1]

    def counts(sel):
        return ([float(w[sel & case & (ev.data == v)].sum()) for v in vals],
                [float(w[sel & ~case & (ev.data == v)].sum()) for v in vals])

    tabs = [counts(sel) for _, sel in strata]
    a, b = counts(mask)
    out = s.output
    z = E.z_of(level)
    rows = []
    if not o.get("or"):
        out.write("\n" + "-" * 74 + "\n", "text")
        out.write(f"{name12(ev.name):>12} |     cases    controls       odds      [{level:g}% Conf. Interval]\n", "text")
        out.write("-" * 13 + "+" + "-" * 60 + "\n", "text")
        for lab, ai, bi in zip(labels, a, b):
            odds = E._div(ai, bi)
            if ai > 0 and bi > 0:
                se = math.sqrt(1 / ai + 1 / bi)
                lo, hi = odds * math.exp(-z * se), odds * math.exp(z * se)
            else:
                lo = hi = NAN
            rows.append((odds, lo, hi))
            out.write(f"{name12(lab):>12} |", "text")
            out.write(f"{ai:>10.0f}{bi:>12.0f}   {_f5(odds):>9}   {_f5(lo):>9}   {_f5(hi):>9}\n", "result")
        out.write("-" * 74 + "\n", "text")
    else:
        bi = 0
        if isinstance(o.get("base"), str):
            bi = int(o["base"]) - 1
            if not 0 <= bi < K:
                raise StataError(198, "base() out of range")   # VERIFY
        out.write("\n" + "-" * 79 + "\n", "text")
        out.write(f"{name12(ev.name):>12} |  Odds Ratio       chi2       P>chi2     [{level:g}% Conf. Interval]\n", "text")
        out.write("-" * 13 + "+" + "-" * 65 + "\n", "text")
        for k, lab in enumerate(labels):
            out.write(f"{name12(lab):>12} |", "text")
            if k == bi:
                out.write(f"{1.0:>12.6f}{'.':>11}{'.':>13}{'.':>15}{'.':>12}\n", "result")
                rows.append((1.0, NAN, NAN))
                continue
            t2 = [(ca[k], ca[bi], co[k], co[bi]) for ca, co in tabs]
            est, lo, hi = E.mh_or(t2, level)
            chi2, pv = E.mh_chi2(t2)
            rows.append((est, lo, hi))
            out.write(f"{_f6(est):>12}{_f2(chi2):>11}{_p4(pv):>13}{_f6(lo):>15}{_f6(hi):>12}\n", "result")
        out.write("-" * 79 + "\n", "text")
    chi_h, df_h = _homogeneity(tabs) if adj else _pearson_2xk(a, b)
    ph = float(st.chi2.sf(chi_h, df_h)) if df_h else NAN
    scores = [float(v) for v in vals]
    U, V = _trend(tabs, scores)
    chi_t = U * U / V if V else NAN
    pt = float(st.chi2.sf(chi_t, 1)) if V else NAN
    if adj:
        out.write(f"\nMantel-Haenszel odds ratios adjusted for {' '.join(adj)}\n", "text") if o.get("or") else None
    out.write(f"\nTest of homogeneity (equal odds): chi2({df_h})  = ", "text")
    out.write(f"{_f2(chi_h):>8}\n", "result")
    out.write(f"{'Pr>chi2  = ':>45}", "text")
    out.write(f"{_p4(ph):>8}\n", "result")
    out.write("\nScore test for trend of odds:     chi2(1)  = ", "text")
    out.write(f"{_f2(chi_t):>8}\n", "result")
    out.write(f"{'Pr>chi2  = ':>45}", "text")
    out.write(f"{_p4(pt):>8}\n", "result")
    if o.get("graph") or o.get("ciplot"):
        _tabodds_graph(s, ev, vals, rows, bool(o.get("ciplot")), bool(o.get("or")))
    _set_r(s, [("p_trend", pt), ("chi2_tr", chi_t), ("p_hom", ph), ("df_hom", float(df_h)),
               ("chi2_hom", chi_h)])


def _pearson_2xk(a, b) -> tuple[float, int]:
    a, b = np.asarray(a, float), np.asarray(b, float)
    n_k = a + b
    keep = n_k > 0
    a, b, n_k = a[keep], b[keep], n_k[keep]
    N, A = n_k.sum(), a.sum()
    if N == 0 or A == 0 or A == N:
        return NAN, 0
    ea, eb = n_k * A / N, n_k * (N - A) / N
    chi2 = float(((a - ea) ** 2 / ea).sum() + ((b - eb) ** 2 / eb).sum())
    return chi2, len(a) - 1


def _f5(x) -> str:
    return "." if x is None or math.isnan(x) or math.isinf(x) else f"{x:.5f}"


def _f6(x) -> str:
    return "." if x is None or math.isnan(x) or math.isinf(x) else f"{x:.6f}"


def _f2(x) -> str:
    return "." if x is None or math.isnan(x) else f"{x:.2f}"


def _tabodds_graph(s: "Session", ev, vals, rows, ci: bool, is_or: bool) -> None:
    from ..graphics.model import Graph, Layer, Panel
    from .graphcmd import finish
    x = [float(v) for v in vals]
    y = [r[0] for r in rows]
    layers = [Layer("connected", x, y, {}, {}, "Odds Ratio" if is_or else "Odds")]
    if ci:
        lo = [r[1] for r in rows]
        hi = [r[2] for r in rows]
        layers.append(Layer("rcap", x, lo, {"y2": hi}, {}, f"{s.settings.get('level', '95')}% CI"))
    opts = {"_xtitle": ev.name, "_ytitle": "Odds Ratio" if is_or else "Odds"}
    finish(s, Graph(family="twoway", opts=opts, panels=[Panel("", layers)]), {})


# ---------------------------------------------------------------------------
# mhodds
# ---------------------------------------------------------------------------

@command("mhodds", "mhodds")
def cmd_mhodds(s: "Session", args: str) -> None:
    p = parse_standard(args)
    names = s.expand_varlist(p.varlist) if p.varlist.strip() else []
    if len(names) < 2:
        raise StataError(198, "invalid syntax")
    o = match_options(p.options, {"by": 2, "binary": 3, "compare": 3, "level": 1})
    level = _level(s, o)
    case_n, exp_n, adj = names[0], names[1], names[2:]
    byv = s.expand_varlist(str(o["by"])) if isinstance(o.get("by"), str) else []
    mask, w = _prep(s, p, case_n, [exp_n] + adj + byv)
    ds = s.data
    case = ds.get(case_n).data != 0
    ev = ds.get(exp_n)
    vals, _ = _levels(s, ev, mask)
    out = s.output
    z = E.z_of(level)
    compare = None
    if isinstance(o.get("compare"), str):
        parts = [x.strip() for x in str(o["compare"]).split(",")]
        if len(parts) != 2:
            raise StataError(198, "compare() requires two values")   # VERIFY
        compare = (float(parts[0]), float(parts[1]))
    binary = bool(compare or len(vals) == 2 or o.get("binary"))
    hi_v, lo_v = compare if compare else (vals[-1], vals[0])

    def estimate(sel_all, control):
        """(or, lb, ub, chi2, p, b, var(b)) for the observations sel_all, stratified by the
        variables control; b is the log odds ratio (or the score slope)."""
        strata = [("", sel_all)] if not control else _strata(s, " ".join(control), sel_all)[1]
        if binary:
            tabs = []
            for _, sel in strata:
                e1, e0 = sel & (ev.data == hi_v), sel & (ev.data == lo_v)
                if o.get("binary") and not compare:
                    e1, e0 = sel & (ev.data != 0), sel & (ev.data == 0)
                tabs.append((float(w[e1 & case].sum()), float(w[e0 & case].sum()),
                             float(w[e1 & ~case].sum()), float(w[e0 & ~case].sum())))
            est, lo, hi = E.mh_or(tabs, level)
            chi2, pv = E.mh_chi2(tabs)
            b = math.log(est) if est > 0 and not math.isinf(est) else NAN
            vb = ((math.log(hi) - b) / z) ** 2 if hi > 0 and not math.isinf(hi) and not math.isnan(b) else NAN
            return est, lo, hi, chi2, pv, b, vb
        tabs = [([float(w[sel & case & (ev.data == v)].sum()) for v in vals],
                 [float(w[sel & ~case & (ev.data == v)].sum()) for v in vals]) for _, sel in strata]
        U, V = _trend(tabs, [float(v) for v in vals])
        chi2 = U * U / V if V else NAN
        pv = float(st.chi2.sf(chi2, 1)) if V else NAN
        # one-step (score) estimate of the odds ratio per unit of the exposure
        b1 = U / V if V else NAN
        return (math.exp(b1), math.exp(b1 - z / math.sqrt(V)), math.exp(b1 + z / math.sqrt(V)),
                chi2, pv, b1, 1 / V if V else NAN)

    if binary:
        title = "Mantel-Haenszel estimate of the odds ratio" if adj or byv else "Maximum likelihood estimate of the odds ratio"   # VERIFY
        out.write(f"\n{title}\n", "text")
        out.write(f"Comparing {ev.name}=={_num(hi_v)} vs. {ev.name}=={_num(lo_v)}\n", "text")
    else:
        out.write(f"\nScore test for trend of odds with {ev.name}\n", "text")
    if adj:
        out.write(f"controlling for {' '.join(adj)}\n", "text")
    if byv:
        _mhodds_by(s, o, estimate, mask, adj, byv, level)
        return
    if not binary:
        out.write("(The Odds Ratio estimate is an approximation to the odds ratio for a one unit increase in "
                  f"{ev.name})\n", "text")
    est, lo, hi, chi2, pv, _, _ = estimate(mask, adj)
    _mhodds_row(s, level, est, chi2, pv, lo, hi)
    _set_r(s, [("ub_or", hi), ("lb_or", lo), ("p", pv), ("chi2", chi2), ("or", est)])


def _mhodds_row(s: "Session", level: float, est, chi2, pv, lo, hi, indent: str = "", lead: str = "\n") -> None:
    out = s.output
    out.write(lead + indent + "-" * 64 + "\n", "text")
    out.write(indent + f" Odds Ratio    chi2(1)        P>chi2        [{level:g}% Conf. Interval]\n", "text")
    out.write(indent + "-" * 64 + "\n", "text")
    out.write(indent + f"{_f6(est):>11}{_f2(chi2):>11}{_p4(pv):>14}{_f5(lo):>16}{_f5(hi):>10}\n", "result")
    out.write(indent + "-" * 64 + "\n", "text")


def _mhodds_by(s: "Session", o, estimate, mask, adj, byv, level) -> None:
    """by(): the estimate within each group of the by() variables, the estimate controlling
    for them, and a test of homogeneity of the group estimates,
    sum (b_j - b)^2 / var(b_j) on the log odds ratios (VERIFY: layout and Stata's test)."""
    out = s.output
    byname, groups = _strata(s, " ".join(byv), mask)
    out.write(f"by {byname}\n", "text")
    out.write("\n" + "-" * 79 + "\n", "text")
    out.write(f"{name12(byname):>12} |  Odds Ratio    chi2(1)       P>chi2       [{level:g}% Conf. Interval]\n", "text")
    out.write("-" * 13 + "+" + "-" * 65 + "\n", "text")
    parts = []
    for lab, sel in groups:
        est, lo, hi, chi2, pv, b, vb = estimate(sel, adj)
        parts.append((b, vb))
        out.write(f"{name12(lab):>12} |", "text")
        out.write(f"{_f6(est):>12}{_f2(chi2):>11}{_p4(pv):>13}{_f6(lo):>15}{_f6(hi):>12}\n", "result")
    out.write("-" * 79 + "\n", "text")
    est, lo, hi, chi2, pv, b, _ = estimate(mask, adj + byv)
    out.write(f"\n    Mantel-Haenszel estimate controlling for {' '.join(adj + byv)}\n", "text")
    _mhodds_row(s, level, est, chi2, pv, lo, hi, indent="    ", lead="")
    ok = [(bj, vj) for bj, vj in parts if not math.isnan(bj) and vj > 0 and not math.isinf(vj)]
    chi_h = sum((bj - b) ** 2 / vj for bj, vj in ok) if not math.isnan(b) else NAN
    df_h = len(ok) - 1
    ph = float(st.chi2.sf(chi_h, df_h)) if df_h > 0 and not math.isnan(chi_h) else NAN
    out.write(f"\nTest of homogeneity of ORs (approx): chi2({df_h})   = ", "text")
    out.write(f"{_f2(chi_h):>8}\n", "result")
    out.write(f"{'Pr>chi2   = ':>49}", "text")
    out.write(f"{_p4(ph):>8}\n", "result")
    _set_r(s, [("ub_or", hi), ("lb_or", lo), ("p", pv), ("chi2", chi2), ("or", est),
               ("df_hom", float(df_h)), ("chi2_hom", chi_h), ("p_hom", ph)])


def _num(v: float) -> str:
    return f"{v:g}"


# ---------------------------------------------------------------------------
# symmetry
# ---------------------------------------------------------------------------

@command("symmetry", "symmetry")
def cmd_symmetry(s: "Session", args: str) -> None:
    p = parse_standard(args)
    names = s.expand_varlist(p.varlist) if p.varlist.strip() else []
    if len(names) != 2:
        raise StataError(198, "invalid syntax")
    o = match_options(p.options, {"notable": 4, "contrib": 4, "exact": 2, "mh": 2, "trend": 2, "cc": 2})
    mask, w = _prep(s, p, names[0], [names[1]])
    ds = s.data
    r, c = ds.get(names[0]), ds.get(names[1])
    levels = sorted(set(r.data[mask].tolist()) | set(c.data[mask].tolist()))
    K = len(levels)
    pos = {v: i for i, v in enumerate(levels)}
    T = np.zeros((K, K))
    for i in np.flatnonzero(mask):
        T[pos[float(r.data[i])], pos[float(c.data[i])]] += w[i]
    out = s.output
    if not o.get("notable"):
        width = 9
        out.write("\n" + "-" * (14 + width * K + 11) + "\n", "text")
        out.write(f"{name12(r.name):>12}  |{name12(c.name):^{width * K}}|\n", "text")
        out.write(" " * 14 + "|" + "".join(f"{_num(v):>{width}}" for v in levels) + f"|{'Total':>10}\n", "text")
        out.write("-" * 14 + "+" + "-" * (width * K) + "+" + "-" * 10 + "\n", "text")
        for i, v in enumerate(levels):
            out.write(f"{_num(v):>12}  |", "text")
            out.write("".join(f"{T[i, j]:>{width}.0f}" for j in range(K)) + f"|{T[i].sum():>10.0f}\n", "result")
        out.write("-" * 14 + "+" + "-" * (width * K) + "+" + "-" * 10 + "\n", "text")
        out.write(f"{'Total':>12}  |", "text")
        out.write("".join(f"{T[:, j].sum():>{width}.0f}" for j in range(K)) + f"|{T.sum():>10.0f}\n", "result")
        out.write("-" * (14 + width * K + 11) + "\n", "text")
    # Bowker's test of symmetry
    chi_s = 0.0
    df_s = 0
    contrib = []
    for i in range(K):
        for j in range(i + 1, K):
            nij, nji = T[i, j], T[j, i]
            if nij + nji > 0:
                if o.get("cc"):
                    term = (abs(nij - nji) - 1) ** 2 / (nij + nji) if abs(nij - nji) >= 1 else 0.0
                else:
                    term = (nij - nji) ** 2 / (nij + nji)
                chi_s += term
                df_s += 1
                contrib.append((levels[i], levels[j], nij, nji, term))
    p_s = float(st.chi2.sf(chi_s, df_s)) if df_s else NAN
    # Stuart-Maxwell test of marginal homogeneity
    d = (T.sum(1) - T.sum(0))[:-1]
    Vm = -(T + T.T)[:-1, :-1]
    np.fill_diagonal(Vm, (T.sum(1) + T.sum(0) - 2 * np.diag(T))[:-1])
    try:
        chi_m = float(d @ np.linalg.solve(Vm, d))
        df_m = K - 1
    except np.linalg.LinAlgError:
        chi_m, df_m = NAN, K - 1
    p_m = float(st.chi2.sf(chi_m, df_m)) if df_m and not math.isnan(chi_m) else NAN
    if o.get("contrib"):
        out.write("\n" + " " * 4 + "Contribution to symmetry chi-squared statistic\n", "text")
        out.write(f"{'cell':>14}{'n(i,j)':>9}{'n(j,i)':>9}{'chi2':>10}\n", "text")
        for vi, vj, nij, nji, term in contrib:
            out.write(f"{'(' + _num(vi) + ',' + _num(vj) + ')':>14}", "text")
            out.write(f"{nij:>9.0f}{nji:>9.0f}{term:>10.3f}\n", "result")
    out.write("\n" + " " * 42 + "chi2      df   Prob>chi2\n", "text")
    out.write("-" * 66 + "\n", "text")
    rows = [("Symmetry (asymptotic)", chi_s, df_s, p_s),
            ("Marginal homogeneity (Stuart-Maxwell)", chi_m, df_m, p_m)]
    res = [("p_sm", p_m), ("df_sm", float(df_m)), ("chi2_sm", chi_m), ("p", p_s), ("df", float(df_s)),
           ("chi2", chi_s), ("N_pair", T.sum())]
    if o.get("mh"):
        # Bickenboller's statistic: Stuart-Maxwell without the diagonal covariance terms (VERIFY)
        Vb = np.diag((T.sum(1) + T.sum(0) - 2 * np.diag(T))[:-1])
        try:
            chi_b = float(d @ np.linalg.solve(Vb, d))
        except np.linalg.LinAlgError:
            chi_b = NAN
        p_b = float(st.chi2.sf(chi_b, K - 1)) if not math.isnan(chi_b) else NAN
        rows.append(("Marginal homogeneity (Bickenboller)", chi_b, K - 1, p_b))
        res = [("p_b", p_b), ("chi2_b", chi_b)] + res
    for label, chi, df, pv in rows:
        out.write(f"{label:<38}", "text")
        out.write(f"{_f2(chi):>10}{df:>8}{_p4(pv):>12}\n", "result")
    out.write("-" * 66 + "\n", "text")
    if o.get("exact"):
        # exact significance of the symmetry test: product of conditional binomials over the
        # off-diagonal pairs, enumerated by the statistic (VERIFY: Stata's exact test)
        p_e = _symmetry_exact([(int(t[2]), int(t[3])) for t in contrib], chi_s)
        out.write(f"{'Symmetry (exact significance probability)':<54}", "text")
        out.write(f"{_p4(p_e):>12}\n", "result")
        res = [("p_exact", p_e)] + res
    _set_r(s, res)


def _symmetry_exact(pairs, observed: float) -> float:
    dist = {0.0: 1.0}
    for nij, nji in pairs:
        n = nij + nji
        new = {}
        for k in range(n + 1):
            pk = float(st.binom.pmf(k, n, 0.5))
            term = round((2 * k - n) ** 2 / n, 10)
            for val, pv in dist.items():
                key = round(val + term, 10)
                new[key] = new.get(key, 0.0) + pv * pk
        dist = new
    return min(1.0, sum(pv for val, pv in dist.items() if val >= observed - 1e-9))
