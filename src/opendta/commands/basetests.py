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


def _g0(x: float) -> str:
    """%g without the leading zero, as display shows numbers."""
    r = f"{x:g}"
    return r.replace("0.", ".", 1) if r.startswith(("0.", "-0.")) else r


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
    heads = f"{'N':>9}{'Observed k':>13}{'Expected k':>13}{'Assumed p':>12}{'Observed p':>13}"
    vals = f"{n:>9}{k:>11}{n * p0:>13.5g}{p0:>14.5f}{k / n:>13.5f}"
    if name:
        out.write("\n" + f"{'Variable':>12} |{heads}\n", "text")
        out.write("-" * 13 + "+" + "-" * 60 + "\n", "text")
        out.write(f"{name12(name):>12} |", "text")
    else:
        out.write("\n" + heads + "\n", "text")
        out.write("-" * 60 + "\n", "text")      # VERIFY: bitesti without the variable column
    out.write(vals + "\n", "result")
    out.write("\n", "text")
    pm = st.binom.pmf(np.arange(n + 1), n, p0)
    pk = st.binom.pmf(k, n, p0)
    other = np.flatnonzero(pm <= pk * (1 + 1e-7))
    if k <= n * p0:
        lo, hi = k, int(other[other > n * p0].min()) if (other > n * p0).any() else n + 1
    else:
        hi, lo = k, int(other[other < n * p0].max()) if (other < n * p0).any() else -1
    rows = ((f"  Pr(k >= {k})", upper, "one-sided"), (f"  Pr(k <= {k})", lower, "one-sided"),
            (f"  Pr(k <= {lo} or k >= {hi})", two, "two-sided"))
    w = max(len(r[0]) for r in rows)
    for lab, val, kind in rows:
        out.write(lab.ljust(w) + " = ", "text")
        out.write(f"{val:.6f}  ({kind} test)\n", "result")
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

def _ties_note(s: "Session", x: np.ndarray, what: str) -> None:
    u = len(np.unique(x))
    if u < len(x):
        s.output.write(f"\nNote: Ties exist in {what};\n      there are {u} unique values out of "
                       f"{len(x)} observations.\n", "text")      # VERIFY: shown only with ties


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
        out.write(" Smaller group       D       P-value  \n", "text")
        out.write(" " + "-" * 35 + "\n", "text")
        for lab, d, pv in ((f"{v.name}:", res["dplus"], res["p1"]),
                           ("Cumulative:", -res["dminus"], res["p2"])):     # VERIFY: sign of D-
            out.write(f" {lab:<18}", "text")
            out.write(f"{d:>8.4f}{pv:>9.3f}\n", "result")
        out.write(f" {'Combined K-S:':<18}", "text")
        out.write(f"{res['d']:>8.4f}{res['p']:>9.3f}\n", "result")
        _ties_note(s, xs, "dataset")
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
    exact = bool(o.get("exact"))
    out.write(" Smaller group       D       P-value" + ("      Exact" if exact else "  ") + "\n", "text")
    out.write(" " + "-" * (46 if exact else 35) + "\n", "text")
    for lab, d, pv in ((t1, res["d1"], res["p1"]), (t2, res["d2"], res["p2"])):
        out.write(f" {lab + ':':<18}", "text")
        out.write(f"{d:>8.4f}{pv:>9.3f}\n", "result")
    out.write(f" {'Combined K-S:':<18}", "text")
    out.write(f"{res['d']:>8.4f}{res['p']:>9.3f}", "result")
    out.write(f"{res['pexact']:>11.3f}\n" if exact else "\n", "result")
    _ties_note(s, np.concatenate([v.data[m1], v.data[m2]]), "combined dataset")
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
    gname = s.expand_varlist(o["by"])[0]
    pad = -(-(11 * k - len(gname)) // 2)
    out.write(f"{'Greater':>10} |\n", "text")
    out.write(f"{'than the':>10} |{'':>{pad}}{gname}\n", "text")

    def _row(lab: str, cells, tot: str, kind: str) -> None:
        out.write(f"{lab:>10} |", "text")
        out.write("".join(f"{c:>{10 if i == 0 else 11}}" for i, c in enumerate(cells)), kind)
        out.write(" |", "text")
        out.write(f"{tot:>10} \n" if kind == "result" else f"{tot:>10}\n", kind)
    _row("median", [lab for _, lab, _m in grp], "Total", "text")
    sep = "-" * 11 + "+" + "-" * (11 * k) + "+" + "-" * 10 + "\n"
    out.write(sep, "text")
    for lab, row in (("no", tab[0]), ("yes", tab[1])):
        _row(lab, [f"{c:.0f}" for c in row], f"{row.sum():.0f}", "result")
    out.write(sep, "text")
    cs = tab.sum(0)
    _row("Total", [f"{c:.0f}" for c in cs], f"{cs.sum():.0f}", "result")
    out.write("\n", "text")
    chi2, pv, df, _ = st.chi2_contingency(tab, correction=False)
    out.write(f"{'Pearson chi2(' + str(df) + ')':>25} = ", "text")
    out.write(f"{chi2:>8.4f}   Pr = {pv:.3f}\n", "result")
    res = [("N", float(tab.sum())), ("chi2", float(chi2)), ("p", float(pv)), ("median", med)]
    if k == 2:
        if o.get("exact"):
            odds, pf = st.fisher_exact(tab)
            out.write(f"{'Fisher' + chr(39) + 's exact':>25} =", "text")
            out.write(f"{pf:>22.3f}\n", "result")
            pf1 = st.fisher_exact(tab, alternative="greater" if odds >= 1 else "less")[1]
            out.write(f"{'1-sided Fisher' + chr(39) + 's exact':>25} =", "text")
            out.write(f"{pf1:>22.3f}\n", "result")
            res += [("p_exact", float(pf)), ("p1_exact", float(pf1))]
        chi2c, pvc, *_ = st.chi2_contingency(tab, correction=True)
        out.write(f"\n{'Continuity corrected:':>24}\n", "text")
        out.write(f"{'Pearson chi2(1)':>25} = ", "text")
        out.write(f"{chi2c:>8.4f}   Pr = {pvc:.3f}\n", "result")
        res += [("chi2_cc", float(chi2c)), ("p_cc", float(pvc))]
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
    gnum = s.data.get(s.expand_varlist(o["by"])[0])
    if gnum.is_string:
        raise StataError(109, "type mismatch")
    mask &= S.valid(gnum.data)
    gv, grp = _groups(s, o["by"], mask)
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
        out.write(f"\n{gnum.name:>14}{'score':>13}{'obs':>13}{'sum of ranks':>15}  \n", "text")
        for i, (val, lab, _m) in enumerate(grp):
            out.write(f"{lab[:14]:>14}", "result")
            out.write(f"{res['l'][i]:>13g}{int(res['n'][i]):>13}{res['R'][i]:>15.9g}  \n", "result")
    out.write("\n", "text")
    out.write(f"{'z':>11}  = ", "text")
    out.write(f"{res['z']:>5.2f}\n", "result")
    out.write("  Prob > |z| = ", "text")
    out.write(f"{res['p']:>5.3f}\n", "result")
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
    zr = round(res["z"], 2)
    pr = round(res["p"], 2)         # VERIFY: Stata shows z and p rounded to two decimals
    for lab, val in ((f" N({v.name} <= {tx})", f"{res['n0']}"), (f" N({v.name} >  {tx})", f"{res['n1']}"),
                     ("obs", f"{res['n0'] + res['n1']}"), ("N(runs)", f"{res['runs']}"),
                     ("z ", _g0(zr)), ("Prob>|z|", _g0(pr))):
        out.write(f"{lab:>14} = ", "text")
        out.write(f"{val}\n", "result")
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
        out.write("\n", "text")
        _kv(s, "Number of obs", f"{n:>12}", 16)
        _kv(s, "Tetrachoric rho", f"{R[0, 1]:>12.4f}", 16)
        _kv(s, "Std error", f"{'.' if math.isnan(SE[0, 1]) else format(SE[0, 1], '.4f'):>12}", 16)
        out.write(f"\nTest of Ho: {names[0]} and {names[1]} are independent\n", "text")
        _kv(s, "2-sided exact P", f"{P[0, 1]:>12.4f}", 16)
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
        _kv(s, "Number of obs", f"{n:>7}", 15)
        _kv(s, "Kendall's tau-a", f"{r['taua']:>12.4f}", 15)
        _kv(s, "Kendall's tau-b", f"{r['taub']:>12.4f}", 15)
        _kv(s, "Kendall's score", f"{r['S']:>7.0f}", 15)
        ties = len(np.unique(vars_[0].data[mask])) < n or len(np.unique(vars_[1].data[mask])) < n
        s.output.write(f"{'SE of score':>15} = ", "text")
        s.output.write(f"{r['se']:>11.3f}" + ("   (corrected for ties)" if ties else "") + "\n", "result")
        out.write(f"\nTest of Ho: {vars_[0].name} and {vars_[1].name} are independent\n", "text")
        s.output.write(f"{'Prob > |z|':>15} = ", "text")
        s.output.write(f"{r['p']:>12.4f}  (continuity corrected)\n", "result")
        out.write("\n", "text")
        _set_r(s, [("N", float(n)), ("tau_a", r["taua"]), ("tau_b", r["taub"]),
                   ("score", r["S"]), ("se_score", r["se"]), ("p", r["p"])])
        return
    k = len(vars_)
    names = [v.name for v in vars_]
    stats = str(o.get("stats", "taua")).split() if isinstance(o.get("stats"), str) else ["taua"]
    labels = {"taua": "tau_a", "taub": "tau_b", "score": "score", "obs": "Number of obs",
              "p": "Sig. level"}
    for sname in stats:
        if sname not in labels:
            raise StataError(198, f"stats({sname}) invalid")          # VERIFY: error text
    pmax = float(o["print"]) if isinstance(o.get("print"), str) else None
    T = np.eye(k)
    Pm = np.zeros((k, k))
    cell: dict = {}
    for i in range(k):
        for j in range(i + 1):
            if o.get("pw"):
                m = S.valid(vars_[i].data) & S.valid(vars_[j].data) & touse(s, p)
            else:
                m = mask
            r = B.kendall(vars_[i].data[m], vars_[j].data[m])
            if i != j:
                T[i, j] = T[j, i] = r["taub"]
                Pm[i, j] = Pm[j, i] = r["p"]
            vals = {"taua": f"{r['taua']:.4f}", "taub": "1.0000" if i == j else f"{r['taub']:.4f}",
                    "score": f"{r['S']:.0f}", "obs": f"{int(m.sum())}",
                    "p": "" if i == j else f"{r['p']:.4f}"}
            if pmax is not None and i != j and r["p"] > pmax:
                vals = {key: "" for key in vals}
            cell[i, j] = vals
    out.write(f"(obs={n:,})\n\n", "text")
    if len(stats) > 1:
        w = 17
        out.write("+" + "-" * w + "+\n", "text")
        out.write("|" + "  Key".ljust(w) + "|\n", "text")
        out.write("|" + "-" * w + "|\n", "text")
        for sname in stats:
            out.write("|" + ("   " + labels[sname]).ljust(w) + "|\n", "text")
        out.write("+" + "-" * w + "+\n\n", "text")
    out.write(f"{'':>13}|" + "".join(f"{nm:>11}" for nm in names) + "\n", "text")
    out.write("-" * 13 + "+" + "-" * (11 * k) + "\n", "text")
    for i in range(k):
        for si, sname in enumerate(stats):
            lab = name12(names[i]) if si == 0 else ""
            out.write(f"{lab:>12} |", "text")
            row = "".join(f"{cell[i, j][sname]:>11}" for j in range(i + (1 if sname != "p" else 0)))
            out.write(row + " \n", "result")
        if len(stats) > 1:
            out.write(f"{'':>13}|\n", "text")
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
    out.write(f"{'':>15}Partial   Semipartial      Partial   Semipartial   Significance\n", "text")
    out.write(f"{'Variable':>11} |    Corr.         Corr.      Corr.^2       Corr.^2          Value\n",
              "text")
    out.write("-" * 12 + "+" + "-" * 65 + "\n", "text")
    for i, v in enumerate(vars_[1:]):
        pc, sc = res["partial"][i], res["semipartial"][i]
        out.write(f"{name12(v.name):>11} |", "text")
        out.write(f"{pc:>9.4f}{sc:>14.4f}{pc * pc:>13.4f}{sc * sc:>14.4f}{res['p'][i]:>15.4f}\n",
                  "result")
    out.write("\n", "text")
    from .matrix import Matrix
    names = [v.name for v in vars_[1:]]
    col = lambda a: Matrix(np.asarray(a, dtype=float).reshape(-1, 1).copy(), names, ["c1"])  # noqa: E731
    s.r = {"N": nobs, "df": float(res["df"]), "p_corr": col(res["partial"]),
           "sp_corr": col(res["semipartial"]), "p_corr2": col(res["partial"] ** 2),
           "sp_corr2": col(res["semipartial"] ** 2), "sig": col(res["p"])}
