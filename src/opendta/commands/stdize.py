"""Direct and indirect standardization of rates ([R] dstdize: dstdize, istdize).

dstdize applies each group's stratum-specific rates to a standard population
(the populations of all groups combined, a group chosen with base(), or one read
from a dataset with using) and reports the adjusted rate with a normal
confidence interval. istdize applies standard stratum rates (from a dataset) to
each group's population and reports the standardized mortality ratio (SMR)
with an exact Poisson interval. Layout details are marked VERIFY and checked
in compat/do/1401_epitab.do."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..lang.syntax import match_options, parse_standard
from ..lazy import lazy
from ._util import touse
from .combine import load_using
from .epitab import _strata
from .matrix import Matrix
from .registry import command
from .summarize import _level, name12

st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

NAN = float("nan")
SYS = M.SYSMISS


def _key_rows(ds, names, mask):
    """{stratum key: [rows]} and the label of each key, for the strata variables."""
    from .inspect import cell_text
    cols = [ds.get(n) for n in names]
    keys: dict[tuple, list[int]] = {}
    for i in np.flatnonzero(mask):
        key = tuple((c.raw[i] if c.is_string else float(c.data[i])) for c in cols)
        keys.setdefault(key, []).append(i)
    labels = {k: " ".join(cell_text(ds, c, rows[0]) for c in cols) for k, rows in keys.items()}
    return keys, labels


def _nonmissing(ds, names, mask):
    for n in names:
        v = ds.get(n)
        mask = mask & (np.array([x != "" for x in v.raw]) if v.is_string else v.data < SYS)
    return mask


def _fmt(x: float, fmt: str = "%10.0g") -> str:
    from .summarize import g9
    return g9(SYS if x is None or math.isnan(x) else x, fmt)


# ---------------------------------------------------------------------------
# dstdize
# ---------------------------------------------------------------------------

@command("dstdize", "dstdize")
def cmd_dstdize(s: "Session", args: str) -> None:
    p = parse_standard(args)
    names = s.expand_varlist(p.varlist) if p.varlist.strip() else []
    if len(names) < 3:
        raise StataError(198, "invalid syntax")
    o = match_options(p.options, {"by": 2, "using": 5, "base": 4, "level": 3, "saving": 4,
                                  "format": 3, "print": 2, "nores": 5})
    if not isinstance(o.get("by"), str):
        raise StataError(198, "option by() required")
    level = _level(s, o)
    z = float(st.norm.ppf(1 - (1 - level / 100) / 2))
    fmt = str(o["format"]) if isinstance(o.get("format"), str) else "%10.0g"
    charn, popn, strat = names[0], names[1], names[2:]
    ds = s.data
    byvars = s.expand_varlist(str(o["by"]))
    mask = _nonmissing(ds, [charn, popn] + strat + byvars, touse(s, p))
    keys, labels = _key_rows(ds, strat, mask)
    order = sorted(keys)
    cases = ds.get(charn).data
    pop = ds.get(popn).data
    byname, groups = _strata(s, " ".join(byvars), mask)
    # standard population
    if p.using is not None or isinstance(o.get("using"), str):
        U = load_using(p.using or str(o["using"]))
        for n in strat + [popn]:
            if not U.has(n):
                raise StataError(111, f"variable {n} not found in using data")   # VERIFY
        umask = _nonmissing(U, strat + [popn], np.ones(U.nobs, dtype=bool))
        ukeys, _ = _key_rows(U, strat, umask)
        upop = U.get(popn).data
        std = {k: float(upop[ukeys[k]].sum()) if k in ukeys else 0.0 for k in order}
    elif o.get("base") not in (None, True):
        b = str(o["base"]).strip().strip('"')
        hit = [sel for lab, sel in groups if lab == b or _same_num(lab, b)]
        if not hit:
            raise StataError(198, f"base({b}) not found")   # VERIFY
        std = {k: float(pop[[i for i in keys[k] if hit[0][i]]].sum()) for k in order}
    else:
        std = {k: float(pop[keys[k]].sum()) for k in order}
    tot_std = sum(std.values())
    P = {k: (std[k] / tot_std if tot_std else NAN) for k in order}
    out = s.output
    summary = []
    for glab, gsel in groups:
        rows = []
        for k in order:
            idx = [i for i in keys[k] if gsel[i]]
            n_k, c_k = float(pop[idx].sum()), float(cases[idx].sum())
            rows.append((k, n_k, c_k))
        N = sum(r[1] for r in rows)
        C = sum(r[2] for r in rows)
        adj = 0.0
        var = 0.0
        # layout of Stata 14 (compat 1401)
        if not o.get("nores"):
            out.write(f"\n{'-' * 58}\n-> {byname}= {glab} \n", "text")
            out.write(" " * 25 + "-----Unadjusted-----  Std.\n", "text")
            out.write(" " * 32 + "Pop.  Stratum  Pop.  \n", "text")
            out.write("  Stratum       Pop.     Cases  Dist. Rate[s] Dst[P]  s*P\n", "text")
            out.write("-" * 58 + "\n", "text")
        for k, n_k, c_k in rows:
            rate = c_k / n_k if n_k else 0.0
            adj += rate * P[k]
            if n_k:
                var += P[k] ** 2 * rate * (1 - rate) / n_k
            if not o.get("nores"):
                out.write(f"{name12(labels[k]):>9}", "text")
                out.write(f"{n_k:>11.0f}{c_k:>10.0f}{(n_k / N if N else 0):>7.3f}{rate:>7.4f}{P[k]:>7.3f}"
                          f"{rate * P[k]:>7.4f}\n", "result")
        crude = C / N if N else NAN
        se = math.sqrt(var)
        lo, hi = adj - z * se, adj + z * se
        if not o.get("nores"):
            out.write("-" * 58 + "\n", "text")
            out.write(f"{'Totals:':<9}", "text")
            out.write(f"{N:>11.0f}{C:>10.0f}", "result")
            out.write("    Adjusted Cases:", "text")
            out.write(f"{adj * N:>9.1f}\n", "result")
            for lab, val in (("Crude Rate:", crude), ("Adjusted Rate:", adj)):
                out.write(f"{lab:>49}", "text")
                out.write(f"{val:>9.4f}\n", "result")
            out.write(f"{f'{level:g}% Conf. Interval:':>42} [", "text")
            out.write(f"{lo:.4f}, {hi:.4f}", "result")
            out.write("]\n", "text")
        summary.append((glab, N, crude, adj, lo, hi, se))
    out.write("\nSummary of Study Populations:\n", "text")
    out.write(f"{name12(byname):>9}{'N':>14}{'Crude':>11}{'Adj_Rate':>13}       Confidence Interval\n", "text")
    out.write(" " + "-" * 74 + "\n", "text")
    for glab, N, crude, adj, lo, hi, se in summary:
        out.write(f"{name12(glab):>9}", "text")
        out.write(f"{N:>14.0f}{crude:>11.6f}{adj:>13.6f}    [{lo:>10.6f},{hi:>12.6f}]\n", "result")
    G = len(summary)
    cols = [f"c{j + 1}" for j in range(G)]

    def mat(j, row="r1"):
        return Matrix(np.array([[g[j] for g in summary]], dtype=float), [row], cols)
    # return list order of Stata 14; the group names go in macros c1, c2, ... (compat 1401;
    # VERIFY: the row names of the matrices other than r(adj))
    s.r = {"Nobs": mat(1), "crude": mat(2), "adj": mat(3, "Adjusted"), "lb_adj": mat(4), "ub_adj": mat(5),
           "se": mat(6)}
    for j, g in enumerate(summary):
        s.r[f"c{j + 1}"] = f'`"{g[0]}"\''
    s.r["by"] = " ".join(byvars)
    s.r["k"] = float(G)


def _same_num(a: str, b: str) -> bool:
    try:
        return float(a) == float(b)
    except ValueError:
        return False


# ---------------------------------------------------------------------------
# istdize
# ---------------------------------------------------------------------------

@command("istdize", "istdize")
def cmd_istdize(s: "Session", args: str) -> None:
    p = parse_standard(args)
    names = s.expand_varlist(p.varlist) if p.varlist.strip() else []
    if len(names) < 3 or p.using is None:
        raise StataError(198, "invalid syntax" if len(names) < 3 else "using required")
    o = match_options(p.options, {"popvars": 4, "rate": 4, "by": 2, "level": 3, "print": 2, "format": 3})
    level = _level(s, o)
    a = (1 - level / 100) / 2
    z = float(st.norm.ppf(1 - a))
    casen, popn, strat = names[0], names[1], names[2:]
    ds = s.data
    byvars = s.expand_varlist(str(o["by"])) if isinstance(o.get("by"), str) else []
    mask = _nonmissing(ds, [casen, popn] + strat + byvars, touse(s, p))
    U = load_using(p.using)
    umask = _nonmissing(U, strat, np.ones(U.nobs, dtype=bool))
    ukeys, _ = _key_rows(U, strat, umask)
    # standard rates by stratum, and the crude standard rate
    if isinstance(o.get("popvars"), str):
        pv = str(o["popvars"]).split()
        if len(pv) != 2:
            raise StataError(198, "popvars() requires two variables")   # VERIFY
        uc, up = U.get(pv[0]).data, U.get(pv[1]).data
        R = {k: float(uc[r].sum()) / float(up[r].sum()) for k, r in ukeys.items() if float(up[r].sum())}
        crude_std = float(uc[umask].sum()) / float(up[umask].sum())
    elif isinstance(o.get("rate"), str):
        rv = str(o["rate"]).split()
        if len(rv) != 2:
            raise StataError(198, "rate() requires a variable and a crude rate")   # VERIFY
        ur = U.get(rv[0]).data
        R = {k: float(ur[r][0]) for k, r in ukeys.items()}
        try:
            crude_std = float(rv[1])
        except ValueError:
            crude_std = float(U.get(rv[1]).data[umask][0])
    else:
        raise StataError(198, "popvars() or rate() required")
    keys, labels = _key_rows(ds, strat, mask)
    order = sorted(keys)
    missing = [k for k in order if k not in R]
    if missing:
        raise StataError(459, "standard rates missing for some strata")   # VERIFY rc/message
    cases, pop = ds.get(casen).data, ds.get(popn).data
    groups = [("", mask)] if not byvars else _strata(s, " ".join(byvars), mask)[1]
    out = s.output
    summary = []
    for glab, gsel in groups:
        # the cases variable holds the group's total and must be constant in it (compat 1401)
        if len(set(cases[gsel].tolist())) > 1:
            raise StataError(499, f"values of {casen} must be the same for all observations in a group")
    for glab, gsel in groups:
        if byvars:
            out.write(f"\n{'-' * 79}\n-> {' '.join(byvars)}= {glab}\n", "text")   # VERIFY
        out.write("\n" + " " * 26 + "-----Study Population-----\n", "text")
        out.write(f"{'':26}{'Cases':>8}{'Stratum':>13}\n", "text")
        out.write(f"{'Stratum':>11}{'Rate':>12}{'Population':>13}{'Observed':>10}{'Rate':>10}"
                  f"{'Expected':>11}\n", "text")
        out.write("-" * 67 + "\n", "text")
        O = float(cases[gsel][0]) if gsel.any() else 0.0
        Ex = N = 0.0
        for k in order:
            idx = [i for i in keys[k] if gsel[i]]
            n_k = float(pop[idx].sum())
            e_k = n_k * R[k]
            Ex += e_k
            N += n_k
            # VERIFY: the stratum observed cases and rate are unknown from a group total
            out.write(f"{name12(labels[k]):>11}", "text")
            out.write(f"{R[k]:>12.4f}{n_k:>13,.0f}{'-':>10}{'-':>10}{e_k:>11.2f}\n", "result")
        out.write("-" * 67 + "\n", "text")
        smr = O / Ex if Ex else NAN
        lo = float(st.chi2.ppf(a, 2 * O)) / 2 / Ex if O > 0 else 0.0
        hi = float(st.chi2.ppf(1 - a, 2 * (O + 1))) / 2 / Ex
        crude = O / N if N else NAN
        adj = crude_std * smr
        alo, ahi = crude_std * lo, crude_std * hi
        out.write(f"{'Totals:':>11}", "text")
        out.write(f"{N:>25,.0f}{O:>10,.0f}", "result")
        out.write("\n", "text")
        for lab, val, f in (("Observed Cases:", O, "{:>11,.0f}"), ("Expected Cases:", Ex, "{:>11.1f}"),
                            ("SMR (Obs/Exp):", smr, "{:>11.2f}")):
            out.write(f"{lab:>52}", "text")
            out.write(f.format(val) + "\n", "result")
        out.write(f"{f'SMR exact {level:g}% Conf. Interval: [':>52}", "text")
        out.write(f"{lo:.4f}, {hi:.4f}", "result")
        out.write("]\n", "text")
        for lab, val in (("Crude Rate:", crude), ("Adjusted Rate:", adj)):
            out.write(f"{lab:>52}", "text")
            out.write(f"{val:>11.4f}\n", "result")
        out.write(f"{f'{level:g}% Conf. Interval: [':>52}", "text")
        out.write(f"{alo:.4f}, {ahi:.4f}", "result")
        out.write("]\n", "text")
        summary.append((glab, O, Ex, smr, lo, hi, crude, adj, alo, ahi))
    if byvars:
        out.write("\nSummary of Study Populations (Rates):\n\n", "text")
        out.write(f"{'Cases':>27}\n", "text")
        out.write(f" {' '.join(byvars):<12}{'Observed':>12}{'Crude':>10}{'Adj_Rate':>10}        Confidence Interval\n",
                  "text")
        out.write(" " + "-" * 74 + "\n", "text")
        for glab, O, Ex, smr, lo, hi, crude, adj, alo, ahi in summary:
            out.write(f" {name12(glab):<12}", "text")
            out.write(f"{O:>12,.0f}{crude:>10.6f}{adj:>10.6f}    [{alo:>10.6f}, {ahi:>10.6f}]\n", "result")
    cols = [g[0] or "r1" for g in summary]

    def mat(j):
        return Matrix(np.array([[g[j] for g in summary]], dtype=float), ["r1"], cols)
    s.r = {"k": float(len(summary)), "cases_obs": mat(1), "cases_exp": mat(2), "smr": mat(3),
           "lb_smr": mat(4), "ub_smr": mat(5), "crude": mat(6), "adj": mat(7), "lb_adj": mat(8), "ub_adj": mat(9)}
