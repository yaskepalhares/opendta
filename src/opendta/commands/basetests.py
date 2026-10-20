"""Small hypothesis tests of the base-statistics block: bitest/bitesti ([R] bitest), ksmirnov
([R] ksmirnov), median ([R] median), nptrend ([R] nptrend), runtest ([R] runtest), tetrachoric
([R] tetrachoric), ktau ([R] spearman) and pcorr ([R] pcorr).

The computations are in opendta.core.basestats. The layouts follow the manual's examples where
they show one; the rest is marked VERIFY and checked in compat/do/1903_base_stats.do."""

from __future__ import annotations

import math
import re
from typing import TYPE_CHECKING

import numpy as np

from ..core import basestats as B
from ..core import missing as M
from ..core import stats as S
from ..core.errors import StataError
from ..lang.syntax import find_top, match_options, parse_standard
from ..lazy import lazy
from ._util import eval_vector, touse
from .nptests import _groups, _set_r
from .registry import command
from .summarize import _corr_table, _level, g9, name12, numeric_vars, weights

st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def _p(x: float) -> str:
    return "." if math.isnan(x) else f"{x:.4f}"


def _kv(s: "Session", label: str, value: str, width: int = 14) -> None:
    s.output.write(f"{label:>{width}} = ", "text")
    s.output.write(f"{value}\n", "result")


def _split_eq(args: str):
    """'varname == #p [if] [in] [, options]' -> (varname, rest-for-parse, number text, options)."""
    c = find_top(args, ",")
    body, opts = (args, "") if c == -1 else (args[:c], args[c + 1:])
    m = re.match(r"\s*(\S+?)\s*==?\s*(\S+)(.*)$", body, re.S)
    if not m:
        raise StataError(198, "invalid syntax")
    return m.group(1), m.group(3), m.group(2), opts


# ---------------------------------------------------------------------------
# bitest, bitesti
# ---------------------------------------------------------------------------

def _bit_report(s: "Session", name: str, n: int, k: int, p0: float, detail: bool) -> None:
    upper, lower, two = B.binomial_test(n, k, p0)
    out = s.output
    out.write("\n" + f"{'Variable':>12} |{'N':>9}{'Observed k':>14}{'Expected k':>13}"
              f"{'Assumed p':>12}{'Observed p':>13}\n", "text")
    out.write("-" * 13 + "+" + "-" * 61 + "\n", "text")
    out.write(f"{name12(name):>12} |", "text")
    out.write(f"{n:>9}{k:>14}{n * p0:>13.5g}{p0:>12.5f}{k / n:>13.5f}\n", "result")
    out.write("\n", "text")
    pm = st.binom.pmf(np.arange(n + 1), n, p0)
    pk = st.binom.pmf(k, n, p0)
    other = np.flatnonzero(pm <= pk * (1 + 1e-7))
    if k <= n * p0:
        lo, hi = k, int(other[other > n * p0].min()) if (other > n * p0).any() else n + 1
    else:
        hi, lo = k, int(other[other < n * p0].max()) if (other < n * p0).any() else -1
    out.write(f"  Pr(k >= {k})".ljust(26) + "= ", "text")
    out.write(f"{upper:.6f}  (one-sided test)\n", "result")
    out.write(f"  Pr(k <= {k})".ljust(26) + "= ", "text")
    out.write(f"{lower:.6f}  (one-sided test)\n", "result")
    out.write(f"  Pr(k <= {lo} or k >= {hi})".ljust(26) + "= ", "text")
    out.write(f"{two:.6f}  (two-sided test)\n", "result")
    _set_r(s, [("N", float(n)), ("k", float(k)), ("p0", p0), ("p_l", lower), ("p_u", upper),
               ("p", two)])


@command("bitest", "bitest", byable=True)
def cmd_bitest(s: "Session", args: str) -> None:
    name, rest, pt, opts = _split_eq(args)
    o = match_options(opts, {"detail": 1}) if opts.strip() else {}
    p = parse_standard(name + " " + rest)
    v = numeric_vars(s, p.varlist, allow_string=False)[0]
    mask = touse(s, p) & S.valid(v.data)
    x = v.data[mask]
    if mask.sum() == 0:
        raise StataError(2000, "no observations")
    if not np.isin(x, (0, 1)).all():
        raise StataError(450, f"{v.name} not binary")      # VERIFY: error text
    p0 = float(eval_vector(s, pt)) if not re.fullmatch(r"[0-9.eE+-]+", pt) else float(pt)
    if not 0 < p0 < 1:
        raise StataError(198, "p0 must be between 0 and 1")
    _bit_report(s, v.name, int(mask.sum()), int(x.sum()), p0, bool(o.get("detail")))


@command("bitesti", "bitesti")
def cmd_bitesti(s: "Session", args: str) -> None:
    c = find_top(args, ",")
    body, opts = (args, "") if c == -1 else (args[:c], args[c + 1:])
    a = body.split()
    if len(a) != 3:
        raise StataError(198, "invalid syntax")
    n, k, p0 = int(float(a[0])), int(float(a[1])), float(a[2])
    if n < 1 or not 0 <= k <= n or not 0 < p0 < 1:
        raise StataError(198, "invalid arguments")
    o = match_options(opts, {"detail": 1}) if opts.strip() else {}
    _bit_report(s, "", n, k, p0, bool(o.get("detail")))


# ---------------------------------------------------------------------------
# ksmirnov
# ---------------------------------------------------------------------------

@command("ksmirnov", "ksmirnov")
def cmd_ksmirnov(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"by": 2, "exact": 1}) if p.options.strip() else {}
    v = numeric_vars(s, p.varlist, allow_string=False)[0]
    mask = touse(s, p) & S.valid(v.data)
    out = s.output
    if p.exp is not None:
        if o.get("by"):
            raise StataError(198, "option by() not allowed")
        th = eval_vector(s, p.exp)
        th = np.broadcast_to(np.asarray(th, dtype=float), (s.data.nobs,))
        # The expression is evaluated at the sorted data; it is a cdf of the variable itself
        x = v.data[mask]
        order = np.argsort(x, kind="stable")
        xs = x[order]
        fvals = np.asarray(th)[mask][order]
        if (np.diff(fvals) < -1e-12).any():
            raise StataError(198, "the expression is not a nondecreasing function of the variable")  # VERIFY
        res = B.ks_one_sample(xs, lambda _: fvals)
        out.write("\nOne-sample Kolmogorov-Smirnov test against theoretical distribution\n", "text")
        out.write(f"{p.exp:>{len(p.exp) + 11}}\n\n", "text")
        out.write(" Smaller group       D       P-value  Corrected\n", "text")
        out.write(" " + "-" * 46 + "\n", "text")
        corr = res["pexact"] if o.get("exact") or res["n"] < 100 else math.nan
        for lab, d, pv in ((f"{v.name}:", res["dplus"], res["p1"]),
                           ("Cumulative:", res["dminus"], res["p2"])):
            out.write(f" {lab:<18}", "text")
            out.write(f"{d:>8.4f}{pv:>11.3f}\n", "result")
        out.write(f" {'Combined K-S:':<18}", "text")
        out.write(f"{res['d']:>8.4f}{res['p']:>11.3f}", "result")
        out.write(f"{corr:>11.3f}\n" if not math.isnan(corr) else "\n", "result")
        out.write("\n", "text")
        s.r = dict(reversed([("N", float(res["n"])), ("D", res["d"]), ("p", res["p"]),
                             ("D_1", res["dplus"]), ("p_1", res["p1"]), ("D_2", res["dminus"]),
                             ("p_2", res["p2"]), ("p_exact", res["pexact"])]))
        return
    if not o.get("by"):
        raise StataError(198, "option by() required")
    _, grp = _groups(s, o["by"], mask)
    if len(grp) != 2:
        raise StataError(198, "group variable must have exactly 2 distinct values")   # VERIFY
    (g1, t1, m1), (g2, t2, m2) = grp
    res = B.ks_two_sample(v.data[m1], v.data[m2])
    out.write("\nTwo-sample Kolmogorov-Smirnov test for equality of distribution functions\n\n",
              "text")
    out.write(" Smaller group       D       P-value  Corrected\n", "text")
    out.write(" " + "-" * 46 + "\n", "text")
    for lab, d, pv in ((t1, res["d1"], res["p1"]), (t2, res["d2"], res["p2"])):
        out.write(f" {lab + ':':<18}", "text")
        out.write(f"{d:>8.4f}{pv:>11.3f}\n", "result")
    corr = res["pexact"] if o.get("exact") or min(res["n1"], res["n2"]) < 50 else math.nan
    out.write(f" {'Combined K-S:':<18}", "text")
    out.write(f"{res['d']:>8.4f}{res['p']:>11.3f}", "result")
    out.write(f"{corr:>11.3f}\n" if not math.isnan(corr) else "\n", "result")
    out.write("\n", "text")
    s.r = dict(reversed([("D_1", res["d1"]), ("p_1", res["p1"]), ("D_2", res["d2"]),
                         ("p_2", res["p2"]), ("D", res["d"]), ("p", res["p"]),
                         ("p_exact", res["pexact"])]))


# ---------------------------------------------------------------------------
# median
# ---------------------------------------------------------------------------

@command("median", "median")
def cmd_median(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"by": 2, "exact": 1, "medianties": 7, "split": 2}) \
        if p.options.strip() else {}
    if not o.get("by"):
        raise StataError(198, "option by() required")
    v = numeric_vars(s, p.varlist, allow_string=False)[0]
    mask = touse(s, p) & S.valid(v.data)
    gv, grp = _groups(s, o["by"], mask)
    if len(grp) < 2:
        raise StataError(198, "need at least 2 groups")
    allm = np.zeros(len(mask), dtype=bool)
    for *_, m in grp:
        allm |= m
    x = v.data[allm]
    med = float(np.median(x))
    ties = str(o.get("medianties", "below"))
    above_n = np.array([float((v.data[m] > med).sum()) for *_, m in grp])
    below_n = np.array([float((v.data[m] < med).sum()) for *_, m in grp])
    eq_n = np.array([float((v.data[m] == med).sum()) for *_, m in grp])
    if ties == "above":
        above_n = above_n + eq_n
    elif ties == "below":
        below_n = below_n + eq_n
    elif ties == "split":
        above_n = above_n + eq_n / 2
        below_n = below_n + eq_n / 2
    elif ties != "drop":
        raise StataError(198, "medianties() must be above, below, drop or split")
    tab = np.vstack([below_n, above_n])
    out = s.output
    k = len(grp)
    out.write("\nMedian test\n\n", "text")
    out.write(f"{'Greater than':>14} |\n", "text")
    out.write(f"{'the median':>14} | " + "".join(f"{lab:>10}" for _, lab, _m in grp)
              + f" |{'Total':>10}\n", "text")
    sep = "-" * 15 + "+" + "-" * (10 * k + 1) + "+" + "-" * 11 + "\n"
    out.write(sep, "text")
    for lab, row in (("no", tab[0]), ("yes", tab[1])):
        out.write(f"{lab:>14} | ", "text")
        out.write("".join(f"{c:>10.0f}" for c in row) + f" |{row.sum():>10.0f}\n", "result")
    out.write(sep, "text")
    out.write(f"{'Total':>14} | ", "text")
    cs = tab.sum(0)
    out.write("".join(f"{c:>10.0f}" for c in cs) + f" |{cs.sum():>10.0f}\n\n", "result")
    chi2, pv, df, _ = st.chi2_contingency(tab, correction=False)
    out.write(f"{'Pearson chi2(' + str(df) + ')':>26} = ", "text")
    out.write(f"{chi2:>8.4f}   Pr = {pv:.3f}\n", "result")
    res = [("N", float(tab.sum())), ("chi2", float(chi2)), ("p", float(pv)), ("median", med)]
    if k == 2:
        chi2c, pvc, *_ = st.chi2_contingency(tab, correction=True)
        out.write(f"{'Continuity corrected:':>26}\n", "text")
        out.write(f"{'Pearson chi2(1)':>26} = ", "text")
        out.write(f"{chi2c:>8.4f}   Pr = {pvc:.3f}\n", "result")
        res += [("chi2_cc", float(chi2c)), ("p_cc", float(pvc))]
        if o.get("exact"):
            odds, pf = st.fisher_exact(tab)
            out.write(f"{'Fisher exact':>26} = ", "text")
            out.write(f"{pf:>8.3f}\n", "result")
            pf1 = st.fisher_exact(tab, alternative="greater" if odds >= 1 else "less")[1]
            out.write(f"{'1-sided Fisher exact':>26} = ", "text")
            out.write(f"{pf1:>8.3f}\n", "result")
            res += [("p_exact", float(pf)), ("p1_exact", float(pf1))]
    out.write("\n", "text")
    _set_r(s, res)


# ---------------------------------------------------------------------------
# nptrend
# ---------------------------------------------------------------------------

@command("nptrend", "nptrend")
def cmd_nptrend(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"by": 2, "nodetail": 3, "score": 2, "carmitage": 1}) \
        if p.options.strip() else {}
    if not o.get("by"):
        raise StataError(198, "option by() required")
    v = numeric_vars(s, p.varlist, allow_string=False)[0]
    mask = touse(s, p) & S.valid(v.data)
    gv, grp = _groups(s, o["by"], mask)
    gnum = s.data.get(s.expand_varlist(o["by"])[0])
    if gnum.is_string:
        raise StataError(109, "type mismatch")
    if o.get("score"):
        sv = numeric_vars(s, o["score"], allow_string=False)[0]
        mask &= S.valid(sv.data)
        gv, grp = _groups(s, o["by"], mask)
    y = v.data[mask]
    g = gnum.data[mask]
    if o.get("score"):
        score = {}
        for val, _lab, m in grp:
            sc = np.unique(sv.data[m])
            if len(sc) != 1:
                raise StataError(198, "score must be constant within group")    # VERIFY
            score[val] = float(sc[0])
    else:
        score = {val: float(val) for val, *_ in grp}
    res = B.cuzick(y, g, score)
    out = s.output
    if not o.get("nodetail"):
        out.write(f"\n{'group':>12} |{'score':>9}{'obs':>9}{'sum of ranks':>15}\n", "text")
        out.write("-" * 13 + "+" + "-" * 33 + "\n", "text")
        for i, (val, lab, _m) in enumerate(grp):
            out.write(f"{lab[:12]:>12} |", "text")
            out.write(f"{res['l'][i]:>9g}{int(res['n'][i]):>9}{res['R'][i]:>15.1f}\n", "result")
    out.write("\n", "text")
    _kv(s, "z", f"{res['z']:.2f}", 14)
    _kv(s, "Prob > |z|", f"{res['p']:.4f}", 14)
    out.write("\n", "text")
    _set_r(s, [("p", res["p"]), ("z", res["z"]), ("T", res["T"]), ("N", float(len(y)))])


# ---------------------------------------------------------------------------
# runtest
# ---------------------------------------------------------------------------

@command("runtest", "runtest")
def cmd_runtest(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"threshold": 1, "continuity": 1, "drop": 1, "split": 1}) \
        if p.options.strip() else {}
    v = numeric_vars(s, p.varlist, allow_string=False)[0]
    mask = touse(s, p) & S.valid(v.data)
    x = v.data[mask]
    if len(x) < 2:
        raise StataError(2001, "insufficient observations")
    thr = float(o["threshold"]) if isinstance(o.get("threshold"), str) else float(np.median(x))
    if o.get("drop"):
        keep = x != thr
        x = x[keep]
    above = x > thr
    if o.get("split") and not o.get("drop"):
        ties = np.flatnonzero(x == thr)
        for i, j in enumerate(ties):
            above[j] = i % 2 == 0
    res = B.runs_test(above, bool(o.get("continuity")))
    out = s.output
    tx = f"{thr:g}"
    out.write("\n", "text")
    _kv(s, f"N({v.name} <= {tx})", f"{res['n0']}", 22)
    _kv(s, f"N({v.name} > {tx})", f"{res['n1']}", 22)
    _kv(s, "obs", f"{res['n0'] + res['n1']}", 22)
    _kv(s, "N(runs)", f"{res['runs']}", 22)
    _kv(s, "z", f"{res['z']:.2f}", 22)
    _kv(s, "Prob>|z|", f"{res['p']:.4f}", 22)
    out.write("\n", "text")
    _set_r(s, [("N", float(res["n0"] + res["n1"])), ("N_below", float(res["n0"])),
               ("N_above", float(res["n1"])), ("n_runs", float(res["runs"])),
               ("mu", res["mu"]), ("sigma2", res["var"]), ("z", res["z"]), ("p", res["p"]),
               ("threshold", thr)])


# ---------------------------------------------------------------------------
# tetrachoric
# ---------------------------------------------------------------------------

def _binary(v, m):
    vals = np.unique(v.data[m])
    if len(vals) > 2:
        raise StataError(198, f"{v.name} has more than two values")      # VERIFY
    return (v.data[m] == vals[-1]).astype(int) if len(vals) == 2 else np.zeros(int(m.sum()), int)


@command("tetrachoric", "tetrachoric")
def cmd_tetrachoric(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"stats": 2, "print": 1, "star": 2, "bonferroni": 1, "sidak": 2,
                                  "pw": 2, "zeroadjust": 1, "matrix": 3}) \
        if p.options.strip() else {}
    vars_ = numeric_vars(s, p.varlist, allow_string=False)
    if len(vars_) < 2:
        raise StataError(102, "too few variables specified")
    mask = touse(s, p)
    for v in vars_:
        mask &= S.valid(v.data)
    n = int(mask.sum())
    if n < 2:
        raise StataError(2001, "insufficient observations")
    cols = [_binary(v, mask) for v in vars_]
    k = len(vars_)
    R = np.eye(k)
    SE = np.zeros((k, k))
    P = np.zeros((k, k))
    for i in range(k):
        for j in range(i):
            t = np.zeros((2, 2))
            for a, b in zip(cols[i], cols[j]):
                t[1 - a, 1 - b] += 1       # row/col 0 = the higher value, so a positive rho is positive
            r = B.tetrachoric(t[::-1, ::-1])
            R[i, j] = R[j, i] = r["rho"]
            SE[i, j] = SE[j, i] = r["se"]
            P[i, j] = P[j, i] = r["p"]
    out = s.output
    names = [v.name for v in vars_]
    from .matrix import Matrix
    if k == 2 and not o.get("matrix"):
        out.write("\nTetrachoric correlation\n\n", "text")
        _kv(s, "Number of obs", f"{n}", 16)
        _kv(s, "Tetrachoric rho", f"{R[0, 1]:.4f}", 16)
        _kv(s, "Std. Err.", "." if math.isnan(SE[0, 1]) else f"{SE[0, 1]:.4f}", 16)
        out.write(f"\nTest of Ho: {names[0]} and {names[1]} are independent\n", "text")
        _kv(s, "2-sided exact P", f"{P[0, 1]:.4f}", 16)
        out.write("\n", "text")
        _set_r(s, [("N", float(n)), ("rho", float(R[0, 1])), ("se_rho", float(SE[0, 1])),
                   ("p", float(P[0, 1]))])
        return
    out.write(f"(obs={n:,})\n", "text")
    cells = [[f"{R[i, j]:.4f}" for j in range(k)] for i in range(k)]
    extra = None
    if o.get("stats") or o.get("print"):
        extra = [[[_p(P[i, j]) if i != j else "" for j in range(k)] for i in range(k)]]
    _corr_table(s, names, cells, extra=extra)
    s.r = {"N": float(n), "Rho": Matrix(R.copy(), names, names)}


# ---------------------------------------------------------------------------
# ktau
# ---------------------------------------------------------------------------

@command("ktau", "ktau")
def cmd_ktau(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"stats": 2, "print": 2, "star": 2, "bonferroni": 1, "sidak": 2,
                                  "pw": 2, "matrix": 3}) if p.options.strip() else {}
    vars_ = numeric_vars(s, p.varlist, allow_string=False)
    if len(vars_) < 2:
        raise StataError(102, "too few variables specified")
    mask = touse(s, p)
    for v in vars_:
        mask &= S.valid(v.data)
    n = int(mask.sum())
    if n < 2:
        raise StataError(2001, "insufficient observations")
    out = s.output
    if len(vars_) == 2 and not o.get("matrix"):
        r = B.kendall(vars_[0].data[mask], vars_[1].data[mask])
        out.write("\n", "text")
        _kv(s, "Number of obs", f"{n:>8}", 16)
        _kv(s, "Kendall's tau-a", f"{r['taua']:>8.4f}", 16)
        _kv(s, "Kendall's tau-b", f"{r['taub']:>8.4f}", 16)
        _kv(s, "Kendall's score", f"{r['S']:>8.0f}", 16)
        _kv(s, "SE of score", f"{r['se']:>8.3f}", 16)
        out.write(f"\nTest of Ho: {vars_[0].name} and {vars_[1].name} are independent\n", "text")
        _kv(s, "Prob > |z|", f"{r['p']:>8.4f}", 16)
        out.write("(continuity corrected)\n", "text")
        out.write("\n", "text")
        _set_r(s, [("N", float(n)), ("tau_a", r["taua"]), ("tau_b", r["taub"]),
                   ("score", r["S"]), ("se_score", r["se"]), ("p", r["p"])])
        return
    k = len(vars_)
    names = [v.name for v in vars_]
    T = np.eye(k)
    Pm = np.zeros((k, k))
    for i in range(k):
        for j in range(i):
            if o.get("pw"):
                m = S.valid(vars_[i].data) & S.valid(vars_[j].data) & touse(s, p)
            else:
                m = mask
            r = B.kendall(vars_[i].data[m], vars_[j].data[m])
            T[i, j] = T[j, i] = r["taub"]
            Pm[i, j] = Pm[j, i] = r["p"]
    out.write(f"(obs={n:,})\n", "text")
    cells = [[f"{T[i, j]:.4f}" for j in range(k)] for i in range(k)]
    extra = None
    if o.get("stats") or o.get("print"):
        extra = [[[_p(Pm[i, j]) if i != j else "" for j in range(k)] for i in range(k)]]
    _corr_table(s, names, cells, extra=extra)
    from .matrix import Matrix
    s.r = {"N": float(n), "Rho": Matrix(T.copy(), names, names)}


# ---------------------------------------------------------------------------
# pcorr
# ---------------------------------------------------------------------------

@command("pcorr", "pcorr")
def cmd_pcorr(s: "Session", args: str) -> None:
    p = parse_standard(args)
    vars_ = numeric_vars(s, p.varlist, allow_string=False)
    if len(vars_) < 2:
        raise StataError(102, "too few variables specified")
    mask = touse(s, p)
    w, wtype, mask = weights(s, p, mask, ("aweight", "fweight", "pweight"))
    for v in vars_:
        mask &= S.valid(v.data)
    n = int(mask.sum())
    if n <= len(vars_):
        raise StataError(2001, "insufficient observations")
    y = vars_[0].data[mask]
    X = np.column_stack([v.data[mask] for v in vars_[1:]])
    ww = w[mask] if w is not None else None
    res = B.partial_corr(y, X, ww)
    out = s.output
    nobs = float(ww.sum()) if ww is not None and wtype == "fweight" else n
    out.write(f"(obs={nobs:,.0f})\n\n", "text")
    out.write(f"Partial and semipartial correlations of {vars_[0].name} with\n\n", "text")
    out.write(f"{'':>13}|{'Partial':>11}{'Semipartial':>14}{'Partial':>12}{'Semipartial':>14}"
              f"{'Significance':>15}\n", "text")
    out.write(f"{'Variable':>12} |{'Corr.':>11}{'Corr.':>14}{'Corr.^2':>12}{'Corr.^2':>14}"
              f"{'Value':>15}\n", "text")
    out.write("-" * 13 + "+" + "-" * 66 + "\n", "text")
    for i, v in enumerate(vars_[1:]):
        pc, sc = res["partial"][i], res["semipartial"][i]
        out.write(f"{name12(v.name):>12} |", "text")
        out.write(f"{pc:>11.4f}{sc:>14.4f}{pc * pc:>12.4f}{sc * sc:>14.4f}{res['p'][i]:>15.4f}\n",
                  "result")
    out.write("\n", "text")
    from .matrix import Matrix
    names = [v.name for v in vars_[1:]]
    col = lambda a: Matrix(np.asarray(a, dtype=float).reshape(-1, 1).copy(), names, ["c1"])  # noqa: E731
    s.r = {"N": nobs, "df": float(res["df"]), "p_corr": col(res["partial"]),
           "sp_corr": col(res["semipartial"]), "p_corr2": col(res["partial"] ** 2),
           "sp_corr2": col(res["semipartial"] ** 2), "sig": col(res["p"])}
