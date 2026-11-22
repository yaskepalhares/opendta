"""The pk suite ([R] pk): pkexamine, pksumm, pkcollapse, pkshape, pkequiv and pkcross. The
measures and tests are in core/pk.py; the layouts follow compat 2008."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core import pk as P
from ..core.dataset import Dataset, Variable
from ..core.errors import StataError
from ..core.varlist import expand
from ..lang.syntax import match_options, parse_standard
from ..lang.words import strip_outer_quotes
from ._util import touse
from .registry import command

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS
_MEASURES = ["auc", "aucline", "aucexp", "auclog", "half", "ke", "cmax", "tomc", "tmax"]
_KEY = {"auc": "auc", "aucline": "auc_line", "aucexp": "auc_exp", "auclog": "auc_ln", "half": "half",
        "ke": "ke", "cmax": "cmax", "tomc": "tomc", "tmax": "tmax"}


def _num(x: float) -> float:
    return SYS if x is None or not math.isfinite(x) else float(x)


def _vars(s: "Session", p, n: int) -> list[str]:
    if not p.varlist.strip():
        raise StataError(100, "varlist required")   # compat 2005
    names = expand(s.data, p.varlist)
    if len(names) < n:
        raise StataError(102, "too few variables specified")
    if len(names) > n:
        raise StataError(103, "too many variables specified")
    return names


def _columns(s: "Session", names: list[str], mask: np.ndarray) -> list[np.ndarray]:
    cols = [s.data.get(n).data.astype(float) for n in names]
    ok = mask.copy()
    for c in cols:
        ok &= c < SYS
    return [c[ok] for c in cols]


def _fmt(x: float, f: str = "%9.0g") -> str:
    from ..core.formats import format_value
    return "." if not math.isfinite(x) else format_value(x, f, pad=False).strip()


# ---------------------------------------------------------------------------
# pkexamine
# ---------------------------------------------------------------------------

@command("pkexamine")
def cmd_pkexamine(s: "Session", args: str) -> None:
    """pkexamine time concentration [if] [in] [, fit(#) trapezoid]."""
    p = parse_standard(args)
    tname, cname = _vars(s, p, 2)
    o = match_options(p.options, {"fit": 1, "trapezoid": 4, "graph": 2, "line": 1, "log": 3,
                                  "exp": 3})
    for opt in ("graph", "line", "log", "exp"):
        if o.get(opt):
            raise StataError(198, f"option {opt} is not available yet in OpenDTA")
    t, c = _columns(s, [tname, cname], touse(s, p))
    if len(t) < 3:
        raise StataError(2001, "insufficient observations")
    fit = int(o["fit"]) if o.get("fit") else 3
    m = P.measures(t, c, fit, bool(o.get("trapezoid")))
    out = s.output
    # layout of compat 2008
    if not (np.any((t == 0) & (c == 0))):
        out.write("\n    Warning: the point (0, 0) is not in your data.  It will be added.\n", "text")
    out.write("\n", "text")
    for label, k, f in (("Maximum concentration", "cmax", "%9.0g"),
                        ("Time of maximum concentration", "tomc", "%9.0g"),
                        ("Time of last observation (Tmax)", "tmax", "%9.0g"),
                        ("Elimination rate", "ke", "%9.4f"), ("Half life", "half", "%9.4f")):
        out.write(f"{label:>63} = {_fmt(m[k], f):>9}\n", "result")
    rule = "    " + "-" * 71 + "\n"
    out.write("\n    Area under the curve\n" + rule, "text")
    out.write("                    |    AUC [0, inf.)    | AUC [0, inf.) |  AUC [0, inf.)\n"
              "     AUC [0, Tmax]  | Linear of log conc. |   Linear fit  | Exponential fit\n"
              "    ----------------+---------------------+---------------+----------------\n",
              "text")
    out.write(f"    {_fmt(m['auc'], '%9.2f'):>10}      |{_fmt(m['auc_ln'], '%9.3f'):>15}      |"
              f"{_fmt(m['auc_line'], '%9.3f'):>11}    |{_fmt(m['auc_exp'], '%9.3f'):>12}\n", "result")
    out.write(rule + f"    Fit based on last {min(fit, len(t))} points.\n", "text")
    # return list shows auc_line auc_exp auc_ln auc half ke cmax tomc tmax
    s.r = {k: _num(m[v]) for k, v in (("tmax", "tmax"), ("tomc", "tomc"), ("cmax", "cmax"),
                                      ("ke", "ke"), ("half", "half"), ("auc", "auc"),
                                      ("auc_ln", "auc_ln"), ("auc_exp", "auc_exp"),
                                      ("auc_line", "auc_line"))}


def _per_subject(s: "Session", idname: str, tname: str, cname: str, mask, fit: int,
                 trapezoid: bool) -> tuple[list, list[dict]]:
    ids, t, c = _columns(s, [idname, tname, cname], mask)
    keys, rows = [], []
    for k in np.unique(ids):
        sel = ids == k
        if sel.sum() < 3:
            continue
        keys.append(float(k))
        rows.append(P.measures(t[sel], c[sel], fit, trapezoid))
    if not rows:
        raise StataError(2001, "insufficient observations")
    return keys, rows


# ---------------------------------------------------------------------------
# pksumm
# ---------------------------------------------------------------------------

@command("pksumm")
def cmd_pksumm(s: "Session", args: str) -> None:
    """pksumm id time concentration [, fit(#) trapezoid stat(measure) nodots notimechk].
    The p-value is Shapiro-Wilk's, missing below 4 subjects (compat 2008 shows . with 3)."""
    p = parse_standard(args)
    idn, tn, cn = _vars(s, p, 3)
    o = match_options(p.options, {"fit": 1, "trapezoid": 4, "stat": 4, "nodots": 5, "notimechk": 6,
                                  "graph": 2})
    if o.get("graph"):
        raise StataError(198, "option graph is not available yet in OpenDTA")
    keys, rows = _per_subject(s, idn, tn, cn, touse(s, p), int(o["fit"]) if o.get("fit") else 3,
                              bool(o.get("trapezoid")))
    if o.get("stat"):
        st_ = strip_outer_quotes(str(o["stat"]))
        if st_ not in _KEY:
            raise StataError(198, f"stat({st_}) invalid")   # VERIFY
    out = s.output
    if not o.get("nodots"):
        out.write("." * len(keys) + "\n", "text")
    # compat 2008: the table has every measure even with stat(), which picks the graphed one
    out.write("\n Summary statistics for the pharmacokinetic measures\n\n", "text")
    out.write(f"{'':49}Number of observations = {len(keys):>5}\n\n", "text")
    out.write("     Measure |      Mean     Median     Variance  Skewness   Kurtosis   p-value\n"
              + "-" * 13 + "+" + "-" * 65 + "\n", "text")
    for m in _MEASURES:
        x = np.array([r[_KEY[m]] for r in rows], dtype=float)
        z = P.summary(x)
        cells = [_fmt(z[k], "%9.2f") for k in ("mean", "median", "var", "skew", "kurt", "p")]
        out.write(f"{m:>12} |" + "".join(f"{v:>{w}}" for v, w in zip(cells, (10, 11, 13, 10, 11, 10)))
                  + "\n", "result")


# ---------------------------------------------------------------------------
# pkcollapse
# ---------------------------------------------------------------------------

@command("pkcollapse")
def cmd_pkcollapse(s: "Session", args: str) -> None:
    """pkcollapse time concentration [if], id(idvar) [stat(list) keep(varlist) fit(#)
    trapezoid force nodots]: one observation per subject with the measures, named
    measure_concentration, in the order of compat 2008."""
    p = parse_standard(args)
    tn, cn = _vars(s, p, 2)
    o = match_options(p.options, {"id": 2, "stat": 4, "keep": 4, "fit": 1, "trapezoid": 4,
                                  "force": 5, "nodots": 5})
    if not o.get("id"):
        raise StataError(198, "option id() required")
    idn = expand(s.data, str(o["id"]))[0]
    keys, rows = _per_subject(s, idn, tn, cn, touse(s, p), int(o["fit"]) if o.get("fit") else 3,
                              bool(o.get("trapezoid")))
    # tmax comes before tomc here (compat 2008)
    stats = strip_outer_quotes(str(o["stat"])).split() if o.get("stat") else \
        ["auc", "aucline", "aucexp", "auclog", "half", "ke", "cmax", "tmax", "tomc"]
    if not o.get("nodots"):
        s.output.write("." * len(keys) + "\n", "text")
    keep = [k for k in expand(s.data, str(o["keep"])) if k != idn] if o.get("keep") else []
    old = s.data
    ids = old.get(idn).data
    new = Dataset()
    new.set_obs(len(keys))
    new.add(Variable(idn, old.get(idn).vtype, np.array(keys)))
    for k in keep:
        v = old.get(k)
        first = [v.data[np.flatnonzero(ids == key)[0]] for key in keys]
        new.add(Variable(k, v.vtype, np.array(first, dtype=object if v.is_string else float)))
    for m in stats:
        if m not in _KEY:
            raise StataError(198, f"stat({m}) invalid")   # VERIFY
        new.add(Variable(f"{m}_{cn}"[:32], "float", np.array([_num(r[_KEY[m]]) for r in rows])))
    new.changed = True
    s.data = new
    s.notify_state()


# ---------------------------------------------------------------------------
# pkshape
# ---------------------------------------------------------------------------

@command("pkshape")
def cmd_pkshape(s: "Session", args: str) -> None:
    """pkshape id sequence period1 period2 [...], order(string) [outcome() treatment()
    carryover() sequence() period()]: one observation per subject and period, with the
    treatment of each period from order() (the letters A, B... coded 1, 2...) and the carryover
    of the previous period's treatment (0 in the first). VERIFY: the coding and the order of the
    variables."""
    p = parse_standard(args)
    if not p.varlist.strip():
        raise StataError(100, "varlist required")
    names = expand(s.data, p.varlist)
    if len(names) < 4:
        raise StataError(102, "too few variables specified")
    o = match_options(p.options, {"order": 1, "outcome": 3, "treatment": 3, "carryover": 5,
                                  "sequence": 3, "period": 3})
    if not o.get("order"):
        raise StataError(198, "option order() required")
    orders = strip_outer_quotes(str(o["order"])).upper().split()
    idn, seqn, pers = names[0], names[1], names[2:]
    if any(len(x) != len(pers) for x in orders):
        raise StataError(198, "order() does not match the number of periods")   # VERIFY
    letters = sorted({ch for x in orders for ch in x})
    code = {ch: i + 1 for i, ch in enumerate(letters)}
    ds = s.data
    ids, seq = ds.get(idn).data, ds.get(seqn).data.astype(float)
    rows = []
    # one block per period, the subjects in their order inside it (compat 2008)
    for i in range(ds.nobs):
        k = int(seq[i]) - 1
        if not 0 <= k < len(orders):
            raise StataError(198, f"sequence {seq[i]:g} has no order()")   # VERIFY
        prev = 0
        for j, pv in enumerate(pers):
            tr = code[orders[k][j]]
            rows.append((ids[i], k + 1, ds.get(pv).data[i], tr, prev, j + 1))
            prev = tr
    rows.sort(key=lambda r: r[5])
    new = Dataset()
    new.set_obs(len(rows))
    cols = list(zip(*rows)) if rows else [[]] * 6
    out_names = [idn, str(o.get("sequence") or "sequence"), str(o.get("outcome") or "outcome"),
                 str(o.get("treatment") or "treat"), str(o.get("carryover") or "carry"),
                 str(o.get("period") or "period")]
    order = [0, 1, 2, 3, 4, 5]
    for j in order:
        new.add(Variable(out_names[j], "float", np.array(cols[j], dtype=float)))
    new.changed = True
    s.data = new
    s.notify_state()


# ---------------------------------------------------------------------------
# pkequiv and pkcross
# ---------------------------------------------------------------------------

@command("pkequiv")
def cmd_pkequiv(s: "Session", args: str) -> None:
    """pkequiv outcome treatment period sequence id [, compare(# #) limit(#) level(#)]: the
    classic confidence interval for bioequivalence (layout of compat 2008). fieller, symmetric,
    anderson and tost are not available yet."""
    p = parse_standard(args)
    yn, tn, pn, sn, idn = _vars(s, p, 5)
    o = match_options(p.options, {"compare": 4, "limit": 3, "level": 1, "noboot": 4, "fieller": 4,
                                  "symmetric": 3, "anderson": 3, "tost": 4})
    for opt in ("fieller", "symmetric", "anderson", "tost"):
        if o.get(opt):
            raise StataError(198, f"option {opt} is not available yet in OpenDTA")
    y, t, per, seq, ident = _columns(s, [yn, tn, pn, sn, idn], touse(s, p))
    levels = sorted(set(t.tolist()))
    if o.get("compare"):
        test, ref = (float(x) for x in str(o["compare"]).split())
    else:
        if len(levels) != 2:
            raise StataError(198, "compare() required with more than two treatments")   # VERIFY
        test, ref = levels[1], levels[0]   # the first treatment is the reference (compat 2008)
    limit = float(o["limit"]) if o.get("limit") else 0.2
    level = float(o["level"]) if o.get("level") else 90.0
    r = P.equivalence(y, t, per, seq, ident, test, ref, level, limit)
    out = s.output
    rule = "      " + "-" * 70 + "\n"
    out.write("\n      Classic confidence interval for bioequivalence\n" + rule
              + "                         [equivalence limits]        [    test limits   ]\n"
              + rule, "text")
    out.write(f"{'difference:':>22}{_fmt(r['eq_lo'], '%10.3f'):>10}{_fmt(r['eq_hi'], '%10.3f'):>13}"
              f"{_fmt(r['lo'], '%10.3f'):>15}{_fmt(r['hi'], '%10.3f'):>13}\n", "result")
    out.write(f"{'ratio:':>22}{_fmt(100 * (1 - limit), '%9.0g') + '%':>10}"
              f"{_fmt(100 * (1 + limit), '%9.0g') + '%':>13}{_fmt(100 * r['rlo'], '%9.3f') + '%':>16}"
              f"{_fmt(100 * r['rhi'], '%9.3f') + '%':>13}\n", "result")
    # VERIFY: the probability (compat 2008 shows .6550 and .9700 where this gives others)
    out.write(rule + f"      Probability test limits are within equivalence limits = {r['p']:>9.4f}\n"
              f"      Note: Reference treatment = {_fmt(ref)}\n", "result")
    # return list shows stddev uci lci (compat 2008); stddev is sqrt(MSE/2)
    s.r = {"lci": r["lo"], "uci": r["hi"], "stddev": math.sqrt(r["mse"] / 2)}


@command("pkcross")
def cmd_pkcross(s: "Session", args: str) -> None:
    """pkcross outcome [if] [in] [, sequence() treatment() period() id() carryover(varname|none)]:
    the ANOVA of a crossover design. VERIFY: the layout, model() and the default variable
    names (those pkshape creates)."""
    p = parse_standard(args)
    yn = _vars(s, p, 1)[0]
    o = match_options(p.options, {"sequence": 3, "treatment": 3, "period": 3, "id": 2,
                                  "carryover": 5, "model": 3, "sequential": 3, "param": 3})
    for opt in ("model", "sequential", "param"):
        if o.get(opt):
            raise StataError(198, f"option {opt} is not available yet in OpenDTA")
    ds = s.data
    seqn = str(o.get("sequence") or "sequence")
    tn = str(o.get("treatment") or "treat")
    pn = str(o.get("period") or "period")
    idn = str(o.get("id") or "id")
    carry = str(o.get("carryover") or "carry")
    names = [yn, seqn, idn, pn, tn] + ([] if carry == "none" else [carry])
    for n in names:
        if not ds.has(n):
            raise StataError(111, f"variable {n} not found")
    cols = _columns(s, names, touse(s, p))
    a = P.anova_crossover(*cols[:5], cols[5] if len(cols) > 5 else None)
    out = s.output
    # layout of compat 2008 (a 2x2 design)
    shown = "" if carry == "none" else carry
    for label, v in (("sequence variable", seqn), ("period variable", pn),
                     ("treatment variable", tn), ("carryover variable", shown), ("id variable", idn)):
        out.write(f"{label:>66} = {v}\n", "text")
    nseq, nper = len(np.unique(cols[1])), len(np.unique(cols[3]))
    design = f"a {nseq}x{nper} crossover study" if nseq == nper == 2 else "a crossover study"   # VERIFY
    rule = "    " + "-" * 22 + "+" + "-" * 50 + "\n"
    out.write(f"\n               Analysis of variance (ANOVA) for {design}\n"
              "     Source of Variation  | Partial SS   df        MS          F     Prob > F\n"
              + rule, "text")

    def line(label, ss, df, ms=None, F=None, pv=None):
        txt = f"{label:>25} |{_fmt(ss, '%10.2f'):>11}{df:>5}"
        if ms is not None:
            txt += f"{_fmt(ms, '%10.2f'):>12}"
        if F is not None:
            txt += f"{_fmt(F, '%9.2f'):>10}{_fmt(pv, '%9.4f'):>12}"
        out.write(txt + "\n", "result")
    rows = {r["term"]: r for r in a["rows"] if r["df"]}   # aliased terms are left out
    out.write("     Intersubjects        |\n", "text")
    for term, label in (("sequence", "Sequence effect"), ("id|sequence", "Residuals")):
        if term in rows:
            r = rows[term]
            # the subjects within sequence are tested against the residual
            line(label, r["ss"], r["df"], r["ms"], r["F"], r["p"])
    out.write(rule + "     Intrasubjects        |\n", "text")
    for term, label in (("treatment", "Treatment effect"), ("carryover", "Carryover effect"),
                        ("period", "Period effect")):
        if term in rows:
            r = rows[term]
            line(label, r["ss"], r["df"], r["ms"], r["F"], r["p"])
    line("Residuals", a["rss"], a["df_r"], a["mse"])
    out.write(rule, "text")
    line("Total", a["tss"], a["n"] - 1)
    if carry != "none":
        sep = P.separability(cols[4], cols[5])
        out.write(f"\n    Omnibus measure of separability of treatment and carryover = "
                  f"{_fmt(100 * sep, '%9.4f'):>9}%\n", "result")
