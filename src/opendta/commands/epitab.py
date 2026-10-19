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


_NOT_YET = ("noshow",)
# options that only make sense with stratified data
_BY_ONLY = ("estandard", "istandard", "standard", "pool", "rd", "ird", "bd", "tarone")


def _not_yet(o) -> None:
    """Options of epitab that OpenDTA does not implement yet: an error rather than a silently
    different result."""
    for k in _NOT_YET:
        if o.get(k):
            raise StataError(198, f"option {k} is not available yet in OpenDTA")


def _check_by(o, immediate: bool = False) -> None:
    """The immediate forms take neither the stratified options nor binomial(); without by(),
    the stratified options are ignored, as Stata does with pool (compat 1601; VERIFY the
    others); with by(), rd and ird need a standardization ("invalid syntax", compat 1601)."""
    if immediate:
        for k in _BY_ONLY + ("binomial",):
            if o.get(k):
                raise StataError(198, f"option {k} not allowed")
        return
    if not o.get("by"):
        return
    if sum(bool(o.get(k)) for k in ("estandard", "istandard", "standard")) > 1:
        raise StataError(198, "only one of estandard, istandard, or standard() is allowed")   # VERIFY
    for k in ("rd", "ird"):
        if o.get(k) and not any(o.get(x) for x in ("estandard", "istandard", "standard")):
            raise StataError(198, "invalid syntax")


def _std_spec(s: "Session", o, strata):
    """The weights of the standardized estimate: ("I. Standardized", None) or
    ("E. Standardized", None) for the internal/external weights the caller computes from the
    cells, ("Standardized", [w per stratum]) for standard(varname), or None."""
    if o.get("istandard"):
        return "I. Standardized", None
    if o.get("estandard"):
        return "E. Standardized", None
    if not isinstance(o.get("standard"), str):
        return None
    v = s.data.get(s.expand_varlist(str(o["standard"]))[0])
    if v.is_string:
        raise StataError(109, "type mismatch")
    ws = []
    for _, sel in strata:
        vals = set(v.data[sel].tolist())
        if len(vals) != 1 or next(iter(vals)) >= SYS:
            raise StataError(198, f"{v.name} not constant within strata")   # VERIFY
        w = next(iter(vals))
        if w < 0:
            raise StataError(198, f"{v.name} has negative values")   # VERIFY
        ws.append(float(w))
    return "Standardized", ws


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


def _binomial(s: "Session", o, mask):
    """binomial(varname): the number of subjects of each binomial frequency record, whose
    case variable holds the number of cases. Returns (n or None, mask)."""
    if not isinstance(o.get("binomial"), str):
        return None, mask
    nv = _var01(s, str(o["binomial"]), "binomial")
    return nv.data, mask & (nv.data < SYS)


def _cells(cv, ev, w, sel, n=None) -> tuple[float, float, float, float]:
    e1 = ev.data != 0
    if n is not None:
        # binomial records: cases and noncases (n - cases) of each record (VERIFY: checks
        # that cases <= n)
        cases, rest = cv.data, n - cv.data
        a = float((w * cases)[sel & e1].sum())
        b = float((w * cases)[sel & ~e1].sum())
        c = float((w * rest)[sel & e1].sum())
        d = float((w * rest)[sel & ~e1].sum())
        return a, b, c, d
    # nonzero = case / exposed (VERIFY: Stata may require 0/1)
    c1 = cv.data != 0
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
        out.write(" " * 57 + extra.split("\n")[0] + "\n", "text")
    if title and not extra:      # cc shows no variable name (compat 1401)
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
        out.write(f"{extra:>13.4f}" if not math.isnan(extra) else f"{'.':>13}", "result")
    out.write("\n", "text")


def _count(x: float) -> str:
    # no thousands separator (compat 1401: 1250, 28010)
    return f"{x:.0f}" if float(x).is_integer() else g9(x).strip()


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
    out.write(f"{'chi2(' + str(df) + ') =':>40}", "text")
    out.write(f"{_r(chi2):>9.2f}" if not math.isnan(chi2) else f"{'.':>9}", "result")
    out.write("  Pr>chi2 = ", "text")
    out.write(f"{_p4(p)}\n", "result")


def _fisher_lines(s: "Session", a, b, c, d) -> tuple[float, float]:
    p1, p2 = E.fisher_2x2(int(a), int(b), int(c), int(d))
    out = s.output
    out.write(f"{'1-sided Fisher' + chr(39) + 's exact P =':>60}", "text")
    out.write(f" {p1:.4f}\n", "result")
    out.write(f"{'2-sided Fisher' + chr(39) + 's exact P =':>60}", "text")
    out.write(f" {p2:.4f}\n", "result")
    return p1, p2


# widths of the estimate and lower-limit columns of the stratified tables (compat 1601)
_COLW = {"cs": (12, 13), "cc": (11, 14), "ir": (11, 14)}


def _strat_table(s: "Session", cmd: str, byname: str, measure: str, level: float, rows, crude, combined,
                 wlabel: str = "M-H Weight", note: str = "") -> None:
    """Stratified summary: one line per stratum, then Crude and the combined estimates
    ([(label, (est, lb, ub), note)]: M-H, directly pooled or standardized)."""
    ew, lw = _COLW[cmd]
    out = s.output
    out.write("\n", "text")
    out.write(f"{name12(byname):>{LBL}} |{measure:>9}       [{level:g}% Conf. Interval]{wlabel:>13}\n", "text")
    out.write("-" * (LBL + 1) + "+" + "-" * 49 + "\n", "text")
    for label, est, lo, hi, w in rows:
        out.write(f"{name12(label):>{LBL}} |", "text")
        out.write(f"{_g(est):>{ew}}{_g(lo):>{lw}}{_g(hi):>11}{_g(w):>13}", "result")
        out.write(f" {note}\n" if note else "\n", "text")
    out.write("-" * (LBL + 1) + "+" + "-" * 49 + "\n", "text")
    if crude is not None:
        out.write(f"{'Crude':>{LBL}} |", "text")
        out.write(f"{_g(crude[0]):>{ew}}{_g(crude[1]):>{lw}}{_g(crude[2]):>11}", "result")
        out.write(f"{'':13} {note}\n" if note else "\n", "text")
    for label, (est, lo, hi), rnote in combined:
        out.write(f"{label:>{LBL}} |", "text")
        out.write(f"{_g(est):>{ew}}{_g(lo):>{lw}}{_g(hi):>11}", "result")
        out.write(f"{'':13} {rnote}\n" if rnote else "\n", "text")


def _homog_lines(s: "Session", cmd: str, tests) -> None:
    """The tests of homogeneity under a closing rule, each command with its own layout
    (compat 1601): tests = [(name, chi2, p, df)]."""
    if not tests:
        return
    out = s.output
    out.write("-" * 67 + "\n", "text")
    for name, chi2, p, df in tests:
        label = f"Test of homogeneity ({name})"
        if cmd == "ir":
            out.write(f" {label:<29}chi2({df}) =", "text")
            out.write(f"{_r(chi2):>10.2f}" if not math.isnan(chi2) else f"{'.':>10}", "result")
        else:
            out.write(f"{label:<31}chi2({df}) =", "text")
            fmt = ">9.3f" if cmd == "cs" else ">9.2f"
            out.write(format(_r(chi2), fmt) if not math.isnan(chi2) else f"{'.':>9}", "result")
        out.write("  Pr>chi2 = ", "text")
        out.write(f"{_p4(p)}\n", "result")


class _R(dict):
    """r() under construction: setting a result again moves it to the top of return list,
    as Stata does (compat 1601)."""

    def put(self, key, value) -> None:
        self.pop(key, None)
        self[key] = value


def _f32(x: float) -> float:
    """Stata stores these test statistics in single precision (compat 1601)."""
    return float(np.float32(x)) if x is not None and not math.isnan(x) else x


def _tb_chi2(a, b, c, d) -> float:
    """The chi2 behind test-based limits: Mantel-Haenszel's, (n - 1)/n times Pearson's
    (compat 1401)."""
    n = a + b + c + d
    chi2, _ = E.pearson_2x2(a, b, c, d)
    return chi2 * (n - 1) / n if n > 1 else NAN


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
    _check_by(o)
    cv, ev, mask, w = _counts(s, p, names[0], names[1])
    n, mask = _binomial(s, o, mask)
    if o.get("by"):
        byname, strata = _strata(s, str(o["by"]), mask)
        _cs_by(s, byname, [(lab, _cells(cv, ev, w, sel, n)) for lab, sel in strata], o,
               _std_spec(s, o, strata))
        return
    _cs_table(s, _cells(cv, ev, w, mask, n), o, ev.name, cv.name)


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
    _check_by(o, immediate=True)
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
    chi_tb = _tb_chi2(a, b, c, d)
    rd = r1 - r0
    rd_lo, rd_hi = E.rd_ci(a, b, c, d, level, tb, chi_tb)
    rr = E.risk_ratio(a, b, c, d)
    rr_lo, rr_hi = E.rr_ci(a, b, c, d, level, tb, chi_tb)
    note = "(tb)" if tb else ""
    _est(s, "Risk difference", rd, rd_lo, rd_hi, note)
    _est(s, "Risk ratio", rr, rr_lo, rr_hi, note)
    # r() keeps the names afe/afp for prevented fractions too (compat 1401)
    if rr >= 1 or math.isnan(rr):
        af, af_lo, af_hi = E.af_from_ratio(rr, rr_lo, rr_hi)
        afp = E._div(rt - r0, rt)
        _est(s, "Attr. frac. ex.", af, af_lo, af_hi, note)
        _est(s, "Attr. frac. pop", afp)
    else:
        af, af_lo, af_hi = 1 - rr, 1 - rr_hi, 1 - rr_lo
        afp = E._div(r0 - rt, r0)
        _est(s, "Prev. frac. ex.", af, af_lo, af_hi, note)
        _est(s, "Prev. frac. pop", afp)
    pairs = []
    if o.get("or"):
        # after the fractions (compat 1401)
        method = "woolf" if o.get("woolf") else "tb" if tb else "cornfield"
        orr = E.odds_ratio(a, b, c, d)
        or_lo, or_hi = E.or_ci(a, b, c, d, level, method, chi_tb)
        _est(s, "Odds ratio", orr, or_lo, or_hi, {"woolf": "(Woolf)", "tb": "(tb)", "cornfield": "(Cornfield)"}[method])
        pairs += [("or", orr), ("lb_or", or_lo), ("ub_or", or_hi)]   # VERIFY: position in r()
    _est_end(s)
    exact = []
    if o.get("exact"):
        # Fisher's exact test replaces the chi2 line (compat 1401)
        p1, p2 = _fisher_lines(s, a, b, c, d)
        exact = [("p_exact", p2), ("p1_exact", p1)]   # VERIFY: names and position in r()
    else:
        _chi2_line(s, chi2, pchi)
    _set_r(s, [("chi2", chi2), ("p", pchi)] + exact
           + [("rd", rd), ("lb_rd", rd_lo), ("ub_rd", rd_hi), ("rr", rr), ("lb_rr", rr_lo),
              ("ub_rr", rr_hi)] + pairs + [("afe", af), ("lb_afe", af_lo), ("ub_afe", af_hi), ("afp", afp)])


def _summary(s: "Session", cmd: str, o, byname: str, measure: str, level: float, rows, crude, mh, direct,
             std, note: str = "", mh_note: str = "") -> list:
    """Draws the stratified table with its combined rows: "Pooled (direct)" with pool, then the
    standardized estimate with estandard/istandard/standard() or else "M-H combined" (compat
    1601). Returns the combined rows shown."""
    combined = []
    if direct is not None:
        combined.append(("Pooled (direct)", direct[:3], ""))
    if std is not None:
        combined.append((std[0], std[1], ""))
    else:
        combined.append(("M-H combined", mh, mh_note))
    _strat_table(s, cmd, byname, measure, level, rows, None if o.get("nocrude") else crude, combined,
                 wlabel="Weight" if std is not None else "M-H Weight", note=note)
    return combined


def _std_weights(spec, inside, outside):
    """Weights per stratum for _std_spec: the internal (exposed) or external (unexposed) group
    sizes, or the user's."""
    label, ws = spec
    return label, ws if ws is not None else inside if label.startswith("I.") else outside


def _cs_by(s: "Session", byname: str, strata, o, spec=None) -> None:
    level = _level(s, o)
    tables = [cells for _, cells in strata]
    rd = bool(o.get("rd"))
    rows, var, m1, v1, m0, v0 = [], [], [], [], [], []
    for label, (a, b, c, d) in strata:
        n1, n0 = a + c, b + d
        r1, r0 = E._div(a, n1), E._div(b, n0)
        m1.append(r1)
        m0.append(r0)
        v1.append(E._div(r1 * (1 - r1), n1))
        v0.append(E._div(r0 * (1 - r0), n0))
        if rd:
            est, (lo, hi) = r1 - r0, E.rd_ci(a, b, c, d, level)
            var.append(v1[-1] + v0[-1])
        else:
            est, (lo, hi) = E.risk_ratio(a, b, c, d), E.rr_ci(a, b, c, d, level)
            var.append(1 / a - 1 / (a + c) + 1 / b - 1 / (b + d) if a > 0 and b > 0 else NAN)
        rows.append([label, est, lo, hi, E._div(b * (a + c), a + b + c + d)])
    tot = tuple(sum(t[i] for t in tables) for i in range(4))
    crude = (E.risk_difference(*tot), *E.rd_ci(*tot, level)) if rd else (E.risk_ratio(*tot), *E.rr_ci(*tot, level))
    mh = None if rd else E.mh_rr(tables, level)
    direct = E.pooled([r[1] for r in rows], var, level, log=not rd) if o.get("pool") else None
    std = None
    if spec is not None:
        label, w = _std_weights(spec, [t[0] + t[2] for t in tables], [t[1] + t[3] for t in tables])
        for r, wi in zip(rows, w):
            r[4] = wi
        std = (label, E.standardized(m1, v1, m0, v0, w, level, difference=rd))
    combined = _summary(s, "cs", o, byname, "RD" if rd else "RR", level, rows, crude, mh, direct, std)
    main = combined[-1][1]
    tests = []
    if direct is not None:
        tests.append(("direct", *direct[3:6]))
    if std is None:
        chi2, ph, df = E.homogeneity_log([r[1] for r in rows], var, mh[0])
        tests.append(("M-H", chi2, ph, df))
    if not o.get("nohom"):
        _homog_lines(s, "cs", tests)
    # r(): the population attributable fraction of the crude table, the combined estimate and
    # the tests (compat 1601; VERIFY the M-H test's names)
    a_, b_, c_, d_ = tot
    rt, r0 = E._div(a_ + b_, a_ + b_ + c_ + d_), E._div(b_, b_ + d_)
    r = _R()
    r.put("afp", E._div(rt - r0, rt))
    key = "rd" if rd else "rr"
    r.put(key, main[0])
    r.put("lb_" + key, main[1])
    r.put("ub_" + key, main[2])
    for name, chi2, _p, df in tests:
        if name == "direct":
            r.put("df", df)
            r.put("chi2_p", _f32(chi2))
        else:
            # cs lists chi2_mh before df (compat 1401)
            r.put("df", df)
            r.put("chi2_mh", _f32(chi2))
    _set_r(s, list(r.items()))


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
    _check_by(o)
    cv, ev, mask, w = _counts(s, p, names[0], names[1])
    n, mask = _binomial(s, o, mask)
    if o.get("by"):
        byname, strata = _strata(s, str(o["by"]), mask)
        _cc_by(s, byname, [(lab, _cells(cv, ev, w, sel, n)) for lab, sel in strata], o,
               _std_spec(s, o, strata))
        return
    _cc_table(s, _cells(cv, ev, w, mask, n), o, ev.name)


@command("cci", "cci")
def cmd_cci(s: "Session", args: str) -> None:
    (a, b, c, d), opt = _immediate_counts(args, 4)
    o = match_options(opt, _CC_OPTS)
    _not_yet(o)
    _check_by(o, immediate=True)
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
    lo, hi = E.or_ci(a, b, c, d, level, method, _tb_chi2(a, b, c, d))
    note = {"exact": "(exact)", "cornfield": "(Cornfield)", "woolf": "(Woolf)", "tb": "(tb)"}[method]
    _est(s, "Odds ratio", orr, lo, hi, note)
    pa = E._div(a, a + b)
    # r() keeps the names afe/afp for prevented fractions too (VERIFY for cc)
    if orr >= 1 or math.isnan(orr):
        af, af_lo, af_hi = E.af_from_ratio(orr, lo, hi)
        afp = af * pa if not math.isnan(af) else NAN
        _est(s, "Attr. frac. ex.", af, af_lo, af_hi, note)
        _est(s, "Attr. frac. pop", afp)
    else:
        af, af_lo, af_hi = 1 - orr, 1 - hi, 1 - lo
        afp = E._div(c, c + d) * af
        _est(s, "Prev. frac. ex.", af, af_lo, af_hi, note)
        _est(s, "Prev. frac. pop", afp)
    _est_end(s)
    exact = []
    if o.get("exact"):
        p1, p2 = _fisher_lines(s, a, b, c, d)
        exact = [("p_exact", p2), ("p1_exact", p1)]   # VERIFY: names and position in r()
    else:
        _chi2_line(s, chi2, pchi)
    # Stata 14's order: afe before its limits (compat 1401)
    _set_r(s, [("p", pchi), ("chi2", chi2)] + exact + [("or", orr), ("lb_or", lo), ("ub_or", hi),
                                                       ("lb_afe", af_lo), ("ub_afe", af_hi),
                                                       ("afe", af), ("afp", afp)])


def _cc_by(s: "Session", byname: str, strata, o, spec=None) -> None:
    level = _level(s, o)
    method = "woolf" if o.get("woolf") else "tb" if o.get("tb") else "cornfield" if o.get("cornfield") else "exact"
    tables = [cells for _, cells in strata]
    rows, var, m1, v1, m0, v0 = [], [], [], [], [], []
    for label, (a, b, c, d) in strata:
        lo, hi = E.or_ci(a, b, c, d, level, method, _tb_chi2(a, b, c, d))
        rows.append([label, E.odds_ratio(a, b, c, d), lo, hi, E._div(b * c, a + b + c + d)])
        var.append(1 / a + 1 / b + 1 / c + 1 / d if min(a, b, c, d) > 0 else NAN)
        # cases per control among the exposed and the unexposed, with the variance of a ratio
        # of counts (VERIFY: Stata's standardized odds ratio)
        o1, o0 = E._div(a, c), E._div(b, d)
        m1.append(o1)
        m0.append(o0)
        v1.append(o1 * o1 * (1 / a + 1 / c) if a > 0 and c > 0 else 0.0)
        v0.append(o0 * o0 * (1 / b + 1 / d) if b > 0 and d > 0 else 0.0)
    tot = tuple(sum(t[i] for t in tables) for i in range(4))
    crude = (E.odds_ratio(*tot), *E.or_ci(*tot, level, method, _tb_chi2(*tot)))
    mh = E.mh_or(tables, level)
    mhc, pmh = E.mh_chi2(tables)
    mh_note = ""
    if method == "tb":
        # test-based limits from the Mantel-Haenszel chi2 of the combined odds ratio (compat 1601)
        mh = (mh[0], *E.ci_test_based(mh[0], mhc, level))
        mh_note = "(tb)"
    direct = E.pooled([r[1] for r in rows], var, level) if o.get("pool") else None
    std = None
    if spec is not None:
        # internal weights: the exposed controls; external: the unexposed controls
        label, w = _std_weights(spec, [t[2] for t in tables], [t[3] for t in tables])
        for r, wi in zip(rows, w):
            r[4] = wi
        std = (label, E.standardized(m1, v1, m0, v0, w, level))   # VERIFY: Stata's limits are wider
    note = {"exact": "(exact)", "cornfield": "(Cornfield)", "woolf": "(Woolf)", "tb": "(tb)"}[method]
    combined = _summary(s, "cc", o, byname, "OR", level, rows, crude, mh, direct, std, note=note,
                        mh_note=mh_note)
    main = combined[-1][1]
    tests = []
    if direct is not None:
        tests.append(("direct", *direct[3:6]))
    if std is None:
        # the M-H test uses the logs of the stratum odds ratios around the M-H estimate
        tests.append(("M-H", *E.homogeneity_log([r[1] for r in rows], var, mh[0])))
    extra = []
    if o.get("bd"):
        extra.append(("B-D", *E.breslow_day(tables, mh[0])))
    if o.get("tarone"):
        extra.append(("Tarone", *E.tarone(tables, mh[0])))
    if not o.get("nohom"):
        _homog_lines(s, "cc", tests + extra)
    if std is None:
        out = s.output
        out.write("\n" + " " * 19 + "Test that combined OR = 1:\n", "text")
        out.write(f"{'Mantel-Haenszel chi2(1) =':>57}", "text")
        out.write(f"{_r(mhc):>10.2f}\n", "result")
        out.write(f"{'Pr>chi2 =':>57}", "text")
        out.write(f"{pmh:>10.4f}\n", "result")
    # r() in Stata 14's order (compat 1601)
    a_, b_, c_, d_ = tot
    oc = E.odds_ratio(*tot)
    af = (oc - 1) / oc if oc > 0 and not math.isinf(oc) else NAN
    r = _R()
    r.put("afp", af * E._div(a_, a_ + b_))   # VERIFY: protective exposure
    r.put("chi2", mhc if std is None else NAN)
    for name, chi2, _p, df in tests:
        r.put("df", df)
        r.put("chi2_p", _f32(chi2))
    for name, chi2, _p, df in extra:
        k = "bd" if name == "B-D" else "t"
        r.put(f"chi2_{k}", chi2)
        r.put(f"df_{k}", df)
    r.put("or", main[0])
    r.put("lb_or", main[1])
    r.put("ub_or", main[2])
    _set_r(s, list(r.items()))


# ---------------------------------------------------------------------------
# ir, iri
# ---------------------------------------------------------------------------

# no exact option: Stata 14 answers "option exact not allowed" (compat 1401)
_IR_OPTS = {"by": 2, "level": 1, "tb": 2, "nocrude": 5, "nohom": 5, "pool": 1,
            "ird": 3, "estandard": 4, "istandard": 4, "standard": 4}


@command("ir", "ir")
def cmd_ir(s: "Session", args: str) -> None:
    p = parse_standard(args)
    names = p.varlist.split()
    if len(names) != 3:
        raise StataError(198, "invalid syntax" if len(names) < 3 else "too many variables specified")
    o = match_options(p.options, _IR_OPTS)
    _not_yet(o)
    _check_by(o)
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
        _ir_by(s, byname, [(lab, cells(sel)) for lab, sel in strata], o, _std_spec(s, o, strata))
        return
    # with variables, the rows are named after the case and time variables (compat 1401)
    _ir_table(s, cells(mask), o, ev.name, (cv.name, tv.name))


@command("iri", "iri")
def cmd_iri(s: "Session", args: str) -> None:
    (a, b, t1, t0), opt = _immediate_counts(args, 4)
    o = match_options(opt, _IR_OPTS)
    _not_yet(o)
    _check_by(o, immediate=True)
    _ir_table(s, (a, b, t1, t0), o, "")


def _ir_table(s: "Session", cells, o, expname: str, labels=("Cases", "Person-time")) -> None:
    a, b, t1, t0 = cells
    level = _level(s, o)
    out = s.output
    _head(s, expname)
    _sep(s)
    _row(s, labels[0], a, b, a + b)
    _row(s, labels[1], t1, t0, t1 + t0)
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
        # r() keeps the names afe/afp (as cs, compat 1401; VERIFY for ir)
        frac = [("afe", pf), ("lb_afe", pf_lo), ("ub_afe", pf_hi), ("afp", pfp)]
    _est_end(s)
    p1, p2 = E.irr_midp(int(a), int(b), t1, t0)
    expected = (a + b) * t1 / (t1 + t0) if t1 + t0 else 0
    rel = ">=" if a >= expected else "<="
    k = f"{a:g}"
    # the probabilities end at column 67, followed by "(exact)" (compat 1401)
    for lab in (f"{'(midp)':>27}   Pr(k{rel}{k}) =", f"{'(midp)':>27} 2*Pr(k{rel}{k}) ="):
        out.write(lab, "text")
        out.write(f"{(p1 if lab.endswith(chr(41) + ' =') and '2*' not in lab else p2):.4f}".rjust(67 - len(lab)), "result")
        out.write(" (exact)\n", "text")
    # r(p) is the one-sided mid-p (compat 1401)
    _set_r(s, [("ird", ird), ("lb_ird", ird_lo), ("ub_ird", ird_hi), ("irr", irr),
               ("lb_irr", lo), ("ub_irr", hi), ("p", p1)] + frac)


def _ir_by(s: "Session", byname: str, strata, o, spec=None) -> None:
    level = _level(s, o)
    ird = bool(o.get("ird"))
    z = E.z_of(level)
    rows, var, m1, v1, m0, v0 = [], [], [], [], [], []
    for label, (a, b, t1, t0) in strata:
        i1, i0 = E._div(a, t1), E._div(b, t0)
        m1.append(i1)
        m0.append(i0)
        v1.append(E._div(a, t1 * t1))
        v0.append(E._div(b, t0 * t0))
        if ird:
            est = i1 - i0
            se = math.sqrt(v1[-1] + v0[-1])
            lo, hi = est - z * se, est + z * se
            var.append(v1[-1] + v0[-1])
        else:
            est = E._div(i1, i0)
            lo, hi = E.irr_exact_ci(int(a), int(b), t1, t0, level) if a + b > 0 else (NAN, NAN)
            var.append(1 / a + 1 / b if a > 0 and b > 0 else NAN)
        rows.append([label, est, lo, hi, E._div(b * t1, t1 + t0)])
    tables = [c for _, c in strata]
    tot = tuple(sum(t[i] for t in tables) for i in range(4))
    a, b, t1, t0 = tot
    if ird:
        est = E._div(a, t1) - E._div(b, t0)
        se = math.sqrt(E._div(a, t1 * t1) + E._div(b, t0 * t0))
        crude = (est, est - z * se, est + z * se)
    else:
        crude = (E._div(E._div(a, t1), E._div(b, t0)), *E.irr_exact_ci(int(a), int(b), t1, t0, level))
    mh = None if ird else E.mh_irr(tables, level)
    direct = E.pooled([r[1] for r in rows], var, level, log=not ird) if o.get("pool") else None
    std = None
    if spec is not None:
        label, w = _std_weights(spec, [t[2] for t in tables], [t[3] for t in tables])
        for r, wi in zip(rows, w):
            r[4] = wi
        std = (label, E.standardized(m1, v1, m0, v0, w, level, difference=ird))
    combined = _summary(s, "ir", o, byname, "IRD" if ird else "IRR", level, rows, crude, mh, direct, std,
                        note="" if ird else "(exact)")
    main = combined[-1][1]
    tests = []
    if direct is not None:
        tests.append(("direct", *direct[3:6]))
    if std is None:
        tests.append(("M-H", *E.homogeneity_log([r[1] for r in rows], var, mh[0])))
    if not o.get("nohom"):
        _homog_lines(s, "ir", tests)
    # r(): the one-sided mid-p and the attributable fractions of the crude table, the combined
    # estimate (the other measure missing) and the tests (compat 1601)
    a_, b_, t1_, t0_ = tot
    irr_c = E._div(E._div(a_, t1_), E._div(b_, t0_))
    lo_c, hi_c = E.irr_exact_ci(int(a_), int(b_), t1_, t0_, level)
    af, af_lo, af_hi = E.af_from_ratio(irr_c, lo_c, hi_c)   # VERIFY: protective exposure
    r = _R()
    r.put("p", E.irr_midp(int(a_), int(b_), t1_, t0_)[0])
    r.put("afe", af)
    r.put("lb_afe", af_lo)
    r.put("ub_afe", af_hi)
    r.put("afp", af * E._div(a_, a_ + b_))
    # the requested measure is put last, so it is listed first (compat 1601)
    for key in (("irr", "ird") if ird else ("ird", "irr")):
        on = (key == "ird") == ird
        r.put(key, main[0] if on else NAN)
        r.put("lb_" + key, main[1] if on else NAN)
        r.put("ub_" + key, main[2] if on else NAN)
    for name, chi2, _p, df in tests:
        if name == "direct":
            r.put("df", df)
            r.put("chi2_p", _f32(chi2))
        else:
            r.put("chi2_mh", _f32(chi2))
            r.put("df", df)
    _set_r(s, list(r.items()))


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
    # relative difference (b - c)/(b + d), standard error (b + c)/(b + d)^(3/2) (compat 1401)
    se_rel = (b + c) / (b + d) ** 1.5 if b + d else NAN
    rows = [("difference", diff, diff - z * se_d - cc_d, diff + z * se_d + cc_d),
            ("ratio", ratio, ratio * math.exp(-z * se_lr), ratio * math.exp(z * se_lr)),
            ("rel. diff.", rd_, rd_ - z * se_rel, rd_ + z * se_rel)]

    def line(label, est, lo=None, hi=None, tail=""):
        out.write(f"        {label:<11}", "text")
        out.write(_g(est), "result")
        if lo is not None:
            out.write(f"{'':5}{_g(lo)}  {_g(hi)}", "result")
        out.write(tail + "\n", "text")
    # layout of compat 1401
    out.write("\nProportion with factor\n", "text")
    line("Cases", p1)
    out.write(f"        {'Controls':<11}", "text")
    out.write(_g(p2), "result")
    out.write(f"     [{level:g}% Conf. Interval]\n", "text")
    out.write(" " * 19 + "---------     --------------------\n", "text")
    for label, est, lo, hi in rows:
        line(label, est, lo, hi)
    out.write("\n", "text")
    line("odds ratio", orr, or_lo, or_hi, f"   {note}")
    # return list order of Stata 14 (compat 1401)
    _set_r(s, [("chi2", chi2), ("or", orr), ("lb_or", or_lo), ("ub_or", or_hi),
               ("D_f", diff), ("lb_D_f", rows[0][2]), ("ub_D_f", rows[0][3]),
               ("R_f", ratio), ("lb_R_f", rows[1][2]), ("ub_R_f", rows[1][3]),
               ("RD_f", rd_), ("lb_RD_f", rows[2][2]), ("ub_RD_f", rows[2][3]), ("p_exact", pex)])
