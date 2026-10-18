"""Tables for epidemiologists ([ST] epitab): cs/csi (cohort studies), cc/cci
(case-control), ir/iri (incidence rates) and mcc/mcci (matched case-control).
The formulas are in opendta.core.epi; this module reads the data and draws the
tables. The column widths follow the manual's examples; details the manual
does not show are marked VERIFY and checked in compat/do/1401_epitab.do."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from ..core import epi as E
from ..core import missing as M
from ..core.errors import StataError
from ..lang.syntax import find_top, match_options, parse_standard
from ._util import touse
from .registry import command
from .summarize import _level, g9, name12, weights

if TYPE_CHECKING:
    from ..session import Session

NAN = float("nan")
SYS = M.SYSMISS


def _r(x: float) -> float:
    """Result for r(): missing instead of nan/inf."""
    return SYS if x is None or math.isnan(x) or math.isinf(x) else float(x)


def _g(x: float, fmt: str = "%9.0g") -> str:
    return g9(_r(x), fmt)


def _p4(p: float) -> str:
    return "." if math.isnan(p) else f"{p:.4f}"


_NOT_YET = ("estandard", "istandard", "standard", "binomial", "pool", "rd", "noshow", "bd", "tarone", "ird")


def _not_yet(o) -> None:
    """Options of epitab that OpenDTA does not implement yet: an error rather than a silently
    different result."""
    for k in _NOT_YET:
        if o.get(k):
            raise StataError(198, f"option {k} is not available yet in OpenDTA")


# ---------------------------------------------------------------------------
# reading the data
# ---------------------------------------------------------------------------

def _var01(s: "Session", name: str, what: str):
    v = s.data.get(s.expand_varlist(name)[0])
    if v.is_string:
        raise StataError(109, "type mismatch")
    return v


def _strata(s: "Session", byvars: str | None, mask: np.ndarray):
    """[(label, mask)] for the groups of by(), in sort order (VERIFY: several by variables)."""
    if not byvars:
        return None, None
    from .inspect import cell_text
    names = s.expand_varlist(byvars)
    ds = s.data
    cols = [ds.get(n) for n in names]
    for c in cols:
        ok = np.array([x != "" for x in c.raw]) if c.is_string else c.data < SYS
        mask = mask & ok
    keys = {}
    idx = np.flatnonzero(mask)
    for i in idx:
        key = tuple((c.raw[i] if c.is_string else float(c.data[i])) for c in cols)
        keys.setdefault(key, []).append(i)
    out = []
    for key in sorted(keys):
        rows = keys[key]
        sel = np.zeros(ds.nobs, dtype=bool)
        sel[rows] = True
        label = ", ".join(cell_text(ds, c, rows[0]) for c in cols)
        out.append((label, sel))
    return " ".join(names), out


def _counts(s: "Session", p, case: str, exp: str, extra_ok=None):
    """Mask of usable observations and the weight vector (fweights)."""
    cv, ev = _var01(s, case, "case"), _var01(s, exp, "exposure")
    mask = touse(s, p) & (cv.data < SYS) & (ev.data < SYS)
    if extra_ok is not None:
        mask &= extra_ok
    w, wtype, mask = weights(s, p, mask, ("fweight",))
    if w is None:
        w = np.ones(s.data.nobs)
    return cv, ev, mask, w


def _cells(cv, ev, w, sel) -> tuple[float, float, float, float]:
    # nonzero = case / exposed (VERIFY: Stata may require 0/1)
    c1, e1 = cv.data != 0, ev.data != 0
    a = float(w[sel & c1 & e1].sum())
    b = float(w[sel & c1 & ~e1].sum())
    c = float(w[sel & ~c1 & e1].sum())
    d = float(w[sel & ~c1 & ~e1].sum())
    return a, b, c, d


# ---------------------------------------------------------------------------
# drawing
# ---------------------------------------------------------------------------

LBL = 16          # width of the label column


def _head(s: "Session", title: str, cols=("Exposed", "Unexposed"), total="Total", extra: str = "") -> None:
    out = s.output
    out.write("\n", "text")
    if extra:
        out.write(" " * 53 + extra.split("\n")[0] + "\n", "text")
    if title:
        out.write(" " * (LBL + 1) + "| " + f"{title:<23}" + "|\n", "text")
    tail = f"{total:>11}" + (f"{extra.split(chr(10))[-1]:>12}" if extra else "")
    out.write(" " * (LBL + 1) + f"|{cols[0]:>10}{cols[1]:>12}  |" + tail + "\n", "text")


def _sep(s: "Session", wide: bool = False) -> None:
    s.output.write("-" * (LBL + 1) + "+" + "-" * 24 + "+" + "-" * (24 if wide else 12) + "\n", "text")


def _row(s: "Session", label: str, x, y, t, fmt="count", extra=None) -> None:
    out = s.output
    out.write(f"{label:>{LBL}} |", "text")
    if fmt == "count":
        out.write(f"{_count(x):>10}{_count(y):>12}  |{_count(t):>11}", "result")
    else:
        out.write(f"{_g(x):>10}{_g(y):>12}  |{_g(t):>11}", "result")
    if extra is not None:
        out.write(f"{_g(extra):>12}", "result")
    out.write("\n", "text")


def _count(x: float) -> str:
    return f"{x:,.0f}" if float(x).is_integer() else g9(x).strip()


def _blank(s: "Session") -> None:
    s.output.write(" " * (LBL + 1) + "|" + " " * 24 + "|\n", "text")


def _est_head(s: "Session", level: float) -> None:
    out = s.output
    out.write(" " * (LBL + 1) + f"|      Point estimate    |    [{level:g}% Conf. Interval]\n", "text")
    out.write(" " * (LBL + 1) + "|" + "-" * 24 + "+" + "-" * 24 + "\n", "text")


def _est(s: "Session", label: str, est, lo=None, hi=None, note: str = "") -> None:
    out = s.output
    out.write(f"{label:>{LBL}} |", "text")
    out.write(f"{_g(est):>17}       |", "result")
    if lo is not None:
        out.write(f"  {_g(lo):>10}  {_g(hi):>10}", "result")
        if note:
            out.write(f" {note}", "text")
    out.write("\n", "text")


def _est_end(s: "Session") -> None:
    s.output.write(" " * (LBL + 1) + "+" + "-" * 49 + "\n", "text")


def _chi2_line(s: "Session", chi2: float, p: float, df: int = 1) -> None:
    out = s.output
    out.write(f"{'chi2(' + str(df) + ') =':>39}", "text")
    out.write(f"{_r(chi2):>9.2f}" if not math.isnan(chi2) else f"{'.':>9}", "result")
    out.write("  Pr>chi2 = ", "text")
    out.write(f"{_p4(p)}\n", "result")


def _fisher_lines(s: "Session", a, b, c, d) -> tuple[float, float]:
    p1, p2 = E.fisher_2x2(int(a), int(b), int(c), int(d))
    out = s.output
    out.write(f"{'1-sided Fisher' + chr(39) + 's exact P =':>47}", "text")
    out.write(f"{p1:>8.4f}\n", "result")
    out.write(f"{'2-sided Fisher' + chr(39) + 's exact P =':>47}", "text")
    out.write(f"{p2:>8.4f}\n", "result")
    return p1, p2


def _strat_table(s: "Session", byname: str, measure: str, level: float, rows, crude, pooled,
                 pooled_label: str = "M-H combined", note: str = "") -> None:
    """Stratified summary: one line per stratum, then Crude and the pooled estimate."""
    out = s.output
    out.write("\n", "text")
    out.write(f"{name12(byname):>{LBL}} |{measure:>9}       [{level:g}% Conf. Interval]   M-H Weight\n", "text")
    out.write("-" * (LBL + 1) + "+" + "-" * 49 + "\n", "text")
    for label, est, lo, hi, w in rows:
        out.write(f"{name12(label):>{LBL}} |", "text")
        out.write(f" {_g(est):>9}    {_g(lo):>9}   {_g(hi):>9}    {_g(w):>9}", "result")
        out.write(f" {note}\n" if note else "\n", "text")
    out.write("-" * (LBL + 1) + "+" + "-" * 49 + "\n", "text")
    if crude is not None:
        out.write(f"{'Crude':>{LBL}} |", "text")
        out.write(f" {_g(crude[0]):>9}    {_g(crude[1]):>9}   {_g(crude[2]):>9}", "result")
        out.write(f"{'':13} {note}\n" if note else "\n", "text")
    out.write(f"{pooled_label:>{LBL}} |", "text")
    out.write(f" {_g(pooled[0]):>9}    {_g(pooled[1]):>9}   {_g(pooled[2]):>9}\n", "result")
    out.write("-" * (LBL + 1) + "+" + "-" * 49 + "\n", "text")


def _homog(s: "Session", label: str, chi2: float, p: float, df: int) -> None:
    out = s.output
    out.write(f"{label:<31}chi2({df}) =", "text")
    out.write(f"{_r(chi2):>9.2f}" if not math.isnan(chi2) else f"{'.':>9}", "result")
    out.write("  Pr>chi2 = ", "text")
    out.write(f"{_p4(p)}\n", "result")


def _set_r(s: "Session", pairs) -> None:
    s.r = {k: _r(v) for k, v in pairs}


# ---------------------------------------------------------------------------
# cs, csi
# ---------------------------------------------------------------------------

_CS_OPTS = {"by": 2, "level": 1, "or": 2, "tb": 2, "woolf": 1, "exact": 2, "pool": 1,
            "nocrude": 5, "nohom": 5, "rd": 2, "binomial": 3, "estandard": 4, "istandard": 4,
            "standard": 4, "noshow": 6}


@command("cs", "cs")
def cmd_cs(s: "Session", args: str) -> None:
    p = parse_standard(args)
    names = p.varlist.split()
    if len(names) != 2:
        raise StataError(198, "invalid syntax" if len(names) < 2 else "too many variables specified")
    o = match_options(p.options, _CS_OPTS)
    _not_yet(o)
    cv, ev, mask, w = _counts(s, p, names[0], names[1])
    if o.get("by"):
        byname, strata = _strata(s, str(o["by"]), mask)
        _cs_by(s, byname, [(lab, _cells(cv, ev, w, sel)) for lab, sel in strata], o)
        return
    _cs_table(s, _cells(cv, ev, w, mask), o, ev.name, cv.name)


def _immediate_counts(args: str, k: int):
    c = find_top(args, ",")
    body, opt = (args, "") if c == -1 else (args[:c], args[c + 1:])
    try:
        vals = [float(x) for x in body.split()]
    except ValueError:
        raise StataError(198, "invalid syntax")
    if len(vals) != k:
        raise StataError(198, "invalid syntax")
    if any(v < 0 for v in vals):
        raise StataError(411, "negative counts not allowed")   # VERIFY rc and message
    return vals, opt


@command("csi", "csi")
def cmd_csi(s: "Session", args: str) -> None:
    (a, b, c, d), opt = _immediate_counts(args, 4)
    o = match_options(opt, _CS_OPTS)
    _not_yet(o)
    _cs_table(s, (a, b, c, d), o, "", "")


def _cs_table(s: "Session", cells, o, expname: str, casename: str) -> None:
    a, b, c, d = cells
    level = _level(s, o)
    out = s.output
    tb = bool(o.get("tb"))
    _head(s, expname)
    _sep(s)
    _row(s, "Cases", a, b, a + b)
    _row(s, "Noncases", c, d, c + d)
    _sep(s)
    _row(s, "Total", a + c, b + d, a + b + c + d)
    _blank(s)
    r1, r0 = E._div(a, a + c), E._div(b, b + d)
    rt = E._div(a + b, a + b + c + d)
    _row(s, "Risk", r1, r0, rt, fmt="g")
    _blank(s)
    _est_head(s, level)
    chi2, pchi = E.pearson_2x2(a, b, c, d)
    rd = r1 - r0
    rd_lo, rd_hi = E.rd_ci(a, b, c, d, level, tb, chi2)
    rr = E.risk_ratio(a, b, c, d)
    rr_lo, rr_hi = E.rr_ci(a, b, c, d, level, tb, chi2)
    note = "(tb)" if tb else ""
    _est(s, "Risk difference", rd, rd_lo, rd_hi, note)
    _est(s, "Risk ratio", rr, rr_lo, rr_hi, note)
    pairs = []
    if o.get("or"):
        method = "woolf" if o.get("woolf") else "tb" if tb else "cornfield"   # VERIFY default
        orr = E.odds_ratio(a, b, c, d)
        or_lo, or_hi = E.or_ci(a, b, c, d, level, method, chi2)
        _est(s, "Odds ratio", orr, or_lo, or_hi, {"woolf": "(Woolf)", "tb": "(tb)", "cornfield": "(Cornfield)"}[method])
        pairs += [("or", orr), ("lb_or", or_lo), ("ub_or", or_hi)]
    if rr >= 1 or math.isnan(rr):
        af, af_lo, af_hi = E.af_from_ratio(rr, rr_lo, rr_hi)
        afp = E._div(rt - r0, rt)
        _est(s, "Attr. frac. ex.", af, af_lo, af_hi, note)
        _est(s, "Attr. frac. pop", afp)
        pairs += [("afe", af), ("lb_afe", af_lo), ("ub_afe", af_hi), ("afp", afp)]
    else:
        # protective exposure: prevented fractions
        pf, pf_lo, pf_hi = 1 - rr, 1 - rr_hi, 1 - rr_lo
        pfp = E._div(r0 - rt, r0)
        _est(s, "Prev. frac. ex.", pf, pf_lo, pf_hi, note)
        _est(s, "Prev. frac. pop", pfp)
        pairs += [("pe", pf), ("lb_pe", pf_lo), ("ub_pe", pf_hi), ("pp", pfp)]
    _est_end(s)
    _chi2_line(s, chi2, pchi)
    exact = []
    if o.get("exact"):
        p1, p2 = _fisher_lines(s, a, b, c, d)
        exact = [("p_exact", p2), ("p1_exact", p1)]
    _set_r(s, [("p", pchi), ("chi2", chi2)] + exact
           + [("rd", rd), ("lb_rd", rd_lo), ("ub_rd", rd_hi), ("rr", rr), ("lb_rr", rr_lo),
              ("ub_rr", rr_hi)] + pairs)


def _cs_by(s: "Session", byname: str, strata, o) -> None:
    level = _level(s, o)
    tables = [cells for _, cells in strata]
    rows, logv = [], []
    for label, (a, b, c, d) in strata:
        rr = E.risk_ratio(a, b, c, d)
        lo, hi = E.rr_ci(a, b, c, d, level)
        n = a + b + c + d
        rows.append((label, rr, lo, hi, E._div(b * (a + c), n)))
        logv.append(1 / a - 1 / (a + c) + 1 / b - 1 / (b + d) if a > 0 and b > 0 else NAN)
    tot = tuple(sum(t[i] for t in tables) for i in range(4))
    crude = (E.risk_ratio(*tot), *E.rr_ci(*tot, level))
    pooled = E.mh_rr(tables, level)
    _strat_table(s, byname, "RR", level, rows, None if o.get("nocrude") else crude, pooled)
    chi2, ph, df = E.homogeneity_log([r[1] for r in rows], logv, pooled[0])
    if not o.get("nohom"):
        _homog(s, "Test of homogeneity (M-H)", chi2, ph, df)
    _set_r(s, [("df", df), ("chi2_mh", chi2), ("rr", pooled[0]), ("lb_rr", pooled[1]),
               ("ub_rr", pooled[2]), ("crude", crude[0]), ("lb_crude", crude[1]), ("ub_crude", crude[2])])


# ---------------------------------------------------------------------------
# cc, cci
# ---------------------------------------------------------------------------

_CC_OPTS = {"by": 2, "level": 1, "exact": 2, "cornfield": 4, "woolf": 1, "tb": 2, "bd": 2,
            "tarone": 3, "nocrude": 5, "nohom": 5, "pool": 1, "binomial": 3, "estandard": 4,
            "istandard": 4, "standard": 4}


def _cc_method(o) -> str:
    return "cornfield" if o.get("cornfield") else "woolf" if o.get("woolf") else "tb" if o.get("tb") else "exact"


@command("cc", "cc")
def cmd_cc(s: "Session", args: str) -> None:
    p = parse_standard(args)
    names = p.varlist.split()
    if len(names) != 2:
        raise StataError(198, "invalid syntax" if len(names) < 2 else "too many variables specified")
    o = match_options(p.options, _CC_OPTS)
    _not_yet(o)
    cv, ev, mask, w = _counts(s, p, names[0], names[1])
    if o.get("by"):
        byname, strata = _strata(s, str(o["by"]), mask)
        _cc_by(s, byname, [(lab, _cells(cv, ev, w, sel)) for lab, sel in strata], o)
        return
    _cc_table(s, _cells(cv, ev, w, mask), o, ev.name)


@command("cci", "cci")
def cmd_cci(s: "Session", args: str) -> None:
    (a, b, c, d), opt = _immediate_counts(args, 4)
    o = match_options(opt, _CC_OPTS)
    _not_yet(o)
    _cc_table(s, (a, b, c, d), o, "")


def _cc_table(s: "Session", cells, o, expname: str) -> None:
    a, b, c, d = cells
    level = _level(s, o)
    method = _cc_method(o)
    _head(s, expname, extra="Proportion\nExposed")
    _sep(s, wide=True)
    _row(s, "Cases", a, b, a + b, extra=E._div(a, a + b))
    _row(s, "Controls", c, d, c + d, extra=E._div(c, c + d))
    _sep(s, wide=True)
    _row(s, "Total", a + c, b + d, a + b + c + d, extra=E._div(a + c, a + b + c + d))
    _blank(s)
    _est_head(s, level)
    chi2, pchi = E.pearson_2x2(a, b, c, d)
    orr = E.odds_ratio(a, b, c, d)
    lo, hi = E.or_ci(a, b, c, d, level, method, chi2)
    note = {"exact": "(exact)", "cornfield": "(Cornfield)", "woolf": "(Woolf)", "tb": "(tb)"}[method]
    _est(s, "Odds ratio", orr, lo, hi, note)
    pa = E._div(a, a + b)
    if orr >= 1 or math.isnan(orr):
        af, af_lo, af_hi = E.af_from_ratio(orr, lo, hi)
        afp = af * pa if not math.isnan(af) else NAN
        _est(s, "Attr. frac. ex.", af, af_lo, af_hi, note)
        _est(s, "Attr. frac. pop", afp)
        frac = [("afe", af), ("lb_afe", af_lo), ("ub_afe", af_hi), ("afp", afp)]
    else:
        pf, pf_lo, pf_hi = 1 - orr, 1 - hi, 1 - lo
        pc = E._div(c, c + d)
        pfp = pc * pf
        _est(s, "Prev. frac. ex.", pf, pf_lo, pf_hi, note)
        _est(s, "Prev. frac. pop", pfp)
        frac = [("pe", pf), ("lb_pe", pf_lo), ("ub_pe", pf_hi), ("pp", pfp)]
    _est_end(s)
    _chi2_line(s, chi2, pchi)
    exact = []
    if o.get("exact"):
        p1, p2 = _fisher_lines(s, a, b, c, d)
        exact = [("p_exact", p2), ("p1_exact", p1)]
    _set_r(s, [("p", pchi), ("chi2", chi2)] + exact + [("or", orr), ("lb_or", lo), ("ub_or", hi)] + frac)


def _cc_by(s: "Session", byname: str, strata, o) -> None:
    level = _level(s, o)
    method = "woolf" if o.get("woolf") else "tb" if o.get("tb") else "cornfield" if o.get("cornfield") else "exact"
    tables = [cells for _, cells in strata]
    rows = []
    for label, (a, b, c, d) in strata:
        chi2, _ = E.pearson_2x2(a, b, c, d)
        lo, hi = E.or_ci(a, b, c, d, level, method, chi2)
        rows.append((label, E.odds_ratio(a, b, c, d), lo, hi, E._div(b * c, a + b + c + d)))
    tot = tuple(sum(t[i] for t in tables) for i in range(4))
    chi2t, _ = E.pearson_2x2(*tot)
    crude = (E.odds_ratio(*tot), *E.or_ci(*tot, level, method, chi2t))
    pooled = E.mh_or(tables, level)
    note = {"exact": "(exact)", "cornfield": "(Cornfield)", "woolf": "(Woolf)", "tb": "(tb)"}[method]
    _strat_table(s, byname, "OR", level, rows, None if o.get("nocrude") else crude, pooled, note=note)
    bd, pbd, df = E.breslow_day(tables, pooled[0])
    if not o.get("nohom"):
        _homog(s, "Test of homogeneity (M-H)", bd, pbd, df)
    mh, pmh = E.mh_chi2(tables)
    out = s.output
    out.write("\n" + " " * 19 + "Test that combined OR = 1:\n", "text")
    out.write(f"{'Mantel-Haenszel chi2(1) =':>56}", "text")
    out.write(f"{_r(mh):>10.2f}\n", "result")
    out.write(f"{'Pr>chi2 =':>56}", "text")
    out.write(f"{pmh:>10.4f}\n", "result")
    _set_r(s, [("df", df), ("chi2_bd" if o.get("bd") else "chi2_mh", bd), ("chi2_p", mh), ("p", pmh),
               ("or", pooled[0]), ("lb_or", pooled[1]), ("ub_or", pooled[2]),
               ("crude", crude[0]), ("lb_crude", crude[1]), ("ub_crude", crude[2])])


# ---------------------------------------------------------------------------
# ir, iri
# ---------------------------------------------------------------------------

_IR_OPTS = {"by": 2, "level": 1, "tb": 2, "exact": 2, "nocrude": 5, "nohom": 5, "pool": 1,
            "ird": 3, "estandard": 4, "istandard": 4, "standard": 4}


@command("ir", "ir")
def cmd_ir(s: "Session", args: str) -> None:
    p = parse_standard(args)
    names = p.varlist.split()
    if len(names) != 3:
        raise StataError(198, "invalid syntax" if len(names) < 3 else "too many variables specified")
    o = match_options(p.options, _IR_OPTS)
    _not_yet(o)
    tv = _var01(s, names[2], "time")
    cv, ev, mask, w = _counts(s, p, names[0], names[1], extra_ok=tv.data < SYS)
    t = tv.data

    def cells(sel):
        e1 = ev.data != 0
        # cases are counts: the case variable may hold the number of cases (VERIFY)
        cases = np.where(cv.data < SYS, cv.data, 0)
        a = float((w * cases)[sel & e1].sum())
        b = float((w * cases)[sel & ~e1].sum())
        t1 = float((w * t)[sel & e1].sum())
        t0 = float((w * t)[sel & ~e1].sum())
        return a, b, t1, t0

    if o.get("by"):
        byname, strata = _strata(s, str(o["by"]), mask)
        _ir_by(s, byname, [(lab, cells(sel)) for lab, sel in strata], o)
        return
    _ir_table(s, cells(mask), o, ev.name)


@command("iri", "iri")
def cmd_iri(s: "Session", args: str) -> None:
    (a, b, t1, t0), opt = _immediate_counts(args, 4)
    o = match_options(opt, _IR_OPTS)
    _not_yet(o)
    _ir_table(s, (a, b, t1, t0), o, "")


def _ir_table(s: "Session", cells, o, expname: str) -> None:
    a, b, t1, t0 = cells
    level = _level(s, o)
    out = s.output
    _head(s, expname)
    _sep(s)
    _row(s, "Cases", a, b, a + b)
    _row(s, "Person-time", t1, t0, t1 + t0)
    _sep(s)
    _blank(s)
    i1, i0, it = E._div(a, t1), E._div(b, t0), E._div(a + b, t1 + t0)
    _row(s, "Incidence rate", i1, i0, it, fmt="g")
    _blank(s)
    _est_head(s, level)
    z = E.z_of(level)
    ird = i1 - i0
    se = math.sqrt(a / t1 ** 2 + b / t0 ** 2) if t1 and t0 else NAN
    ird_lo, ird_hi = ird - z * se, ird + z * se
    irr = E._div(i1, i0)
    if o.get("tb"):
        # test-based limits from the score chi2 of the rate comparison
        n = a + b
        e = n * t1 / (t1 + t0)
        v = n * t1 * t0 / (t1 + t0) ** 2
        chi2 = (a - e) ** 2 / v if v else NAN
        lo, hi = E.ci_test_based(irr, chi2, level)
        if chi2 > 0:
            ird_lo, ird_hi = ird - z * abs(ird) / math.sqrt(chi2), ird + z * abs(ird) / math.sqrt(chi2)
        note = "(tb)"
    else:
        lo, hi = E.irr_exact_ci(int(a), int(b), t1, t0, level)
        note = "(exact)"
    _est(s, "Inc. rate diff.", ird, ird_lo, ird_hi, "(tb)" if o.get("tb") else "")
    _est(s, "Inc. rate ratio", irr, lo, hi, note)
    if irr >= 1 or math.isnan(irr):
        af, af_lo, af_hi = E.af_from_ratio(irr, lo, hi)
        afp = af * E._div(a, a + b) if not math.isnan(af) else NAN
        _est(s, "Attr. frac. ex.", af, af_lo, af_hi, note)
        _est(s, "Attr. frac. pop", afp)
        frac = [("afe", af), ("lb_afe", af_lo), ("ub_afe", af_hi), ("afp", afp)]
    else:
        pf, pf_lo, pf_hi = 1 - irr, 1 - hi, 1 - lo
        pfp = pf * E._div(t1 * i0, t1 * i0 + t0 * i0) if i0 else NAN   # VERIFY: prevented fraction in population
        _est(s, "Prev. frac. ex.", pf, pf_lo, pf_hi, note)
        _est(s, "Prev. frac. pop", pfp)
        frac = [("pe", pf), ("lb_pe", pf_lo), ("ub_pe", pf_hi), ("pp", pfp)]
    _est_end(s)
    if o.get("exact"):
        p1, p2 = E.irr_exact_p(int(a), int(b), t1, t0)
        tag = "(exact)"
    else:
        p1, p2 = E.irr_midp(int(a), int(b), t1, t0)
        tag = "(midp)"
    expected = (a + b) * t1 / (t1 + t0) if t1 + t0 else 0
    rel = ">=" if a >= expected else "<="
    k = f"{a:g}"
    out.write(f"{tag:>24}   Pr(k{rel}{k}) = ", "text")
    out.write(f"{p1:>18.4f}", "result")
    out.write(f" {tag}\n", "text")
    out.write(f"{tag:>24} 2*Pr(k{rel}{k}) = ", "text")
    out.write(f"{p2:>18.4f}", "result")
    out.write(f" {tag}\n", "text")
    _set_r(s, [("p", p2), ("ird", ird), ("lb_ird", ird_lo), ("ub_ird", ird_hi), ("irr", irr),
               ("lb_irr", lo), ("ub_irr", hi)] + frac)


def _ir_by(s: "Session", byname: str, strata, o) -> None:
    level = _level(s, o)
    rows, logv = [], []
    for label, (a, b, t1, t0) in strata:
        irr = E._div(E._div(a, t1), E._div(b, t0))
        lo, hi = E.irr_exact_ci(int(a), int(b), t1, t0, level) if a + b > 0 else (NAN, NAN)
        rows.append((label, irr, lo, hi, E._div(b * t1, t1 + t0)))
        logv.append(1 / a + 1 / b if a > 0 and b > 0 else NAN)
    tables = [c for _, c in strata]
    tot = tuple(sum(t[i] for t in tables) for i in range(4))
    a, b, t1, t0 = tot
    crude = (E._div(E._div(a, t1), E._div(b, t0)), *E.irr_exact_ci(int(a), int(b), t1, t0, level))
    pooled = E.mh_irr(tables, level)
    _strat_table(s, byname, "IRR", level, rows, None if o.get("nocrude") else crude, pooled, note="(exact)")
    chi2, ph, df = E.homogeneity_log([r[1] for r in rows], logv, pooled[0])
    if not o.get("nohom"):
        _homog(s, "Test of homogeneity (M-H)", chi2, ph, df)
    _set_r(s, [("df", df), ("chi2_mh", chi2), ("irr", pooled[0]), ("lb_irr", pooled[1]),
               ("ub_irr", pooled[2]), ("crude", crude[0]), ("lb_crude", crude[1]), ("ub_crude", crude[2])])


# ---------------------------------------------------------------------------
# mcc, mcci
# ---------------------------------------------------------------------------

@command("mcc", "mcc")
def cmd_mcc(s: "Session", args: str) -> None:
    p = parse_standard(args)
    names = p.varlist.split()
    if len(names) != 2:
        raise StataError(198, "invalid syntax" if len(names) < 2 else "too many variables specified")
    o = match_options(p.options, {"level": 1, "tb": 2})
    case_e, ctrl_e, mask, w = _counts(s, p, names[0], names[1])
    x, y = case_e.data != 0, ctrl_e.data != 0
    a = float(w[mask & x & y].sum())
    b = float(w[mask & x & ~y].sum())
    c = float(w[mask & ~x & y].sum())
    d = float(w[mask & ~x & ~y].sum())
    _mcc_table(s, (a, b, c, d), o)


@command("mcci", "mcci")
def cmd_mcci(s: "Session", args: str) -> None:
    (a, b, c, d), opt = _immediate_counts(args, 4)
    o = match_options(opt, {"level": 1, "tb": 2})
    _mcc_table(s, (a, b, c, d), o)


def _mcc_table(s: "Session", cells, o) -> None:
    """a: both exposed; b: case exposed, control not; c: control exposed, case not; d: neither."""
    a, b, c, d = cells
    n = a + b + c + d
    level = _level(s, o)
    z = E.z_of(level)
    out = s.output
    out.write("\n" + " " * (LBL + 1) + "| Controls               |\n", "text")
    out.write(f"{'Cases':<{LBL}} |   Exposed   Unexposed  |      Total\n", "text")
    _sep(s)
    _row(s, "Exposed", a, b, a + b)
    _row(s, "Unexposed", c, d, c + d)
    _sep(s)
    _row(s, "Total", a + c, b + d, n)
    chi2 = (b - c) ** 2 / (b + c) if b + c else NAN
    from ..lazy import lazy
    st = lazy("scipy.stats")
    pchi = float(st.chi2.sf(chi2, 1)) if not math.isnan(chi2) else NAN
    k = int(min(b, c))
    pex = min(1.0, 2 * float(st.binom.cdf(k, int(b + c), 0.5))) if b + c else NAN
    out.write("\nMcNemar's chi2(1) =", "text")
    out.write(f"{_r(chi2):>10.2f}" if not math.isnan(chi2) else f"{'.':>10}", "result")
    out.write("    Prob > chi2 = ", "text")
    out.write(f"{_p4(pchi)}\n", "result")
    out.write("Exact McNemar significance probability       = ", "text")
    out.write(f"{_p4(pex)}\n", "result")
    # proportions with the factor and their comparisons (paired data)
    p1, p2 = E._div(a + b, n), E._div(a + c, n)
    diff = p1 - p2
    # Wald limits with a continuity correction of 1/n (Fleiss)
    se_d = math.sqrt((b + c) - (b - c) ** 2 / n) / n if n else NAN
    cc_d = 1 / n if n else NAN
    ratio = E._div(p1, p2)
    se_lr = math.sqrt((b + c) / ((a + b) * (a + c))) if (a + b) and (a + c) else NAN
    rd_ = E._div(b - c, b + d)
    # relative difference (b - c)/(b + d); variance from Fleiss (VERIFY against Stata)
    se_rel = math.sqrt(((b + c + d) * (b * c + b * d + c * d) - b * c * d) / ((b + d) ** 3 * (b + c + d))) \
        if b + d and b + c + d else NAN
    orr = E._div(b, c)
    if o.get("tb"):
        or_lo, or_hi = E.ci_test_based(orr, chi2, level)
        note = "(tb)"
    else:
        lo_p = 0.0 if b == 0 else float(st.beta.ppf((1 - level / 100) / 2, b, c + 1))
        hi_p = 1.0 if c == 0 else float(st.beta.ppf(1 - (1 - level / 100) / 2, b + 1, c))
        or_lo = lo_p / (1 - lo_p) if lo_p < 1 else math.inf
        or_hi = hi_p / (1 - hi_p) if hi_p < 1 else math.inf
        note = "(exact)"
    out.write("\nProportion with factor\n", "text")
    out.write(f"{'Cases':>14}       ", "text")
    out.write(f"{_g(p1)}\n", "result")
    out.write(f"{'Controls':>17}    ", "text")
    out.write(f"{_g(p2)}", "result")
    out.write(f"     [{level:g}% Conf. Interval]\n", "text")
    out.write(" " * 21 + "---------    --------------------\n", "text")
    # ratio: log limits with var = (b + c)/((a + b)(a + c)) (VERIFY: Stata's limits are narrower)
    rows = [("difference", diff, diff - z * se_d - cc_d, diff + z * se_d + cc_d),
            ("ratio", ratio, ratio * math.exp(-z * se_lr), ratio * math.exp(z * se_lr)),
            ("rel. diff.", rd_, rd_ - z * se_rel, rd_ + z * se_rel)]
    for label, est, lo, hi in rows:
        out.write(f"{label:>19}  ", "text")
        out.write(f"{_g(est)}    {_g(lo)}   {_g(hi)}\n", "result")
    out.write(f"\n{'odds ratio':>19}  ", "text")
    out.write(f"{_g(orr)}    {_g(or_lo)}   {_g(or_hi)}", "result")
    out.write(f"   {note}\n", "text")
    _set_r(s, [("p_exact", pex), ("or", orr), ("lb_or", or_lo), ("ub_or", or_hi),
               ("D_f", diff), ("lb_D_f", rows[0][2]), ("ub_D_f", rows[0][3]),
               ("R_f", ratio), ("lb_R_f", rows[1][2]), ("ub_R_f", rows[1][3]),
               ("RD_f", rd_), ("lb_RD_f", rows[2][2]), ("ub_RD_f", rows[2][3]), ("chi2", chi2)])
