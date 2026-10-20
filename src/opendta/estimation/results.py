"""Estimation results: e(b), e(V), e(sample) and the coefficient table
([R] estimation commands, [P] ereturn, [R] estimates).

The table follows the Stata 14 layout (78 columns):

    ------------------------------------------------------------------------------
           price |      Coef.   Std. Err.      t    P>|t|     [95% Conf. Interval]
    -------------+----------------------------------------------------------------
             mpg |  -49.51222   86.15604    -0.57   0.567    -221.3025     122.278
    ...
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.formats import format_value
from ..lazy import lazy

# imported on first use (start-up time; see opendta/lazy.py)
st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session
    from .fvars import Column

SYS = M.SYSMISS
LINE = "-" * 78


def g9(x: float, fmt: str = "%9.0g") -> str:
    if x is None or not np.isfinite(x) or abs(x) >= SYS:
        return "."
    return format_value(float(x), fmt, pad=False).strip()


def abbrev(name: str, width: int = 12) -> str:
    """Abbreviates as Stata does in tables: start~end."""
    if len(name) <= width:
        return name
    return name[: width - 1] + "~" if "#" not in name else name[: width - 1] + "~"


@dataclass
class CoefRow:
    """One row of the table."""
    kind: str                     # "coef", "header", "blank", "base", "omitted", "empty"
    label: str = ""
    b: float = 0.0
    se: float = 0.0
    eq: str = ""
    name: str = ""                # name in e(b), for estimates table
    noeform: bool = False         # never exponentiated (e.g. /ln_p in streg)
    tight: bool = False           # label over 13 columns, touching the bar (var(_cons)| in me)
    zero_lo: bool = False         # variance: the lower confidence limit is truncated at zero


@dataclass
class Estimates:
    """What an estimation leaves for postestimation."""
    cmd: str
    depvar: str
    names: list[str]                       # colnames of e(b), with b./o.
    b: np.ndarray
    V: np.ndarray
    N: int
    stat: str = "t"                        # "t" ou "z"
    df_r: float | None = None
    rows: list[CoefRow] = field(default_factory=list)
    eqnames: list[str] = field(default_factory=list)


def post(s: "Session", est: Estimates, scalars: list[tuple[str, float]],
         macros: list[tuple[str, str]], sample: np.ndarray, *, extra_mats=()) -> None:
    """Stores e(): Stata lists in the reverse order of storage (see progcmd)."""
    from ..commands.matrix import Matrix
    e: dict = {}
    for k, v in reversed(macros):
        if v != "":
            e[k] = v
    for k, v in reversed(scalars):
        e[k] = float(v) if v is not None else SYS
    k = len(est.names)
    eqs = est.eqnames or ["_"] * k
    mats = [("b", Matrix(est.b.reshape(1, -1), ["y1"], list(est.names))),
            ("V", Matrix(est.V, list(est.names), list(est.names)))] + list(extra_mats)
    for name, m in reversed(mats):
        m.roweq = eqs if name == "V" else ["_"]
        m.coleq = eqs
        e[name] = m
    # matrices first in insertion = listed last; scalars and macros keep the order above
    s.e = e
    s.data.esample = np.asarray(sample, dtype=bool).copy()
    s._estimates = est


def datasignature(s: "Session", varnames: list[str], sample: np.ndarray) -> str:
    """e(datasignature): "N:k:checksum1:checksum2". Stata's checksum algorithm is undocumented;
    OpenDTA uses two CRC-32 sums of the values (tools/compare.py ignores the checksums)."""
    import zlib
    m = np.asarray(sample, dtype=bool)
    c1 = c2 = 0
    for name in varnames:
        v = s.data.get(name)
        vals = v.data[m]
        raw = "\x00".join(vals.tolist()).encode() if v.is_string else \
            np.ascontiguousarray(vals, dtype=np.float64).tobytes()
        c1 = zlib.crc32(name.encode() + raw, c1)
        c2 = zlib.adler32(raw, c2)
    return f"{int(m.sum())}:{len(varnames)}:{c1}:{c2}"


def current(s: "Session") -> Estimates:
    from ..core.errors import StataError
    est = getattr(s, "_estimates", None)
    if est is None or "b" not in s.e:
        raise StataError(301, "last estimates not found")
    return est


# ---------------------------------------------------------------------------
# tabela
# ---------------------------------------------------------------------------

def pvalue(stat: str, x: float, df: float | None) -> float:
    if not np.isfinite(x):
        return SYS
    if stat == "t" and df is not None and df > 0:
        return float(2 * st.t.sf(abs(x), df))
    return float(2 * st.norm.sf(abs(x)))


def crit(stat: str, level: float, df: float | None) -> float:
    a = 1 - level / 100
    if stat == "t" and df is not None and df > 0:
        return float(st.t.ppf(1 - a / 2, df))
    return float(st.norm.ppf(1 - a / 2))


def level_text(level: float) -> str:
    return f"{level:g}"


def coef_table(s: "Session", depvar: str, rows: list[CoefRow], *, stat: str, df: float | None,
               level: float, vcetype: str = "", coef_title: str = "Coef.",
               footer: bool = True, eq_headers: bool = False, eform: bool = False,
               top_line: str | None = None, depvar_left: bool = False) -> None:
    out = s.output
    lv = level_text(level)
    ci = f"[{lv}% Conf. Interval]"
    # a label that touches the bar and does not fit in 13 columns widens the first column
    W = max([13] + [len(r.label) + 1 for r in rows if r.tight])
    X = W - 1
    LINE_ = "-" * (len(LINE) + W - 13)
    out.write(LINE_ + "\n", "text")
    if top_line is not None:
        # bootstrap: "Observed   Bootstrap ... Normal-based" (VERIFY alinhamento)
        out.write(top_line + "\n", "text")
    elif vcetype:
        # VERIFY alignment of the label (Robust, Bootstrap...)
        lab = vcetype
        start = 28 + (9 - len(lab)) // 2      # centrado sobre "Std. Err."
        out.write(f"{'':>{X}} |" + " " * max(start - 14, 1) + lab + "\n", "text")
    stat_h = "t    P>|t|" if stat == "t" else "z    P>|z|"
    # depvar_left: dfuller, regress writes "D.w          |" (compat 1101)
    dep = f"{abbrev(depvar):<{W}}|" if depvar_left else f"{abbrev(depvar):>{X}} |"
    out.write(f"{dep} {coef_title:>10}   Std. Err.      {stat_h}     {ci:>20}\n", "text")
    out.write("-" * W + "+" + "-" * 64 + "\n", "text")
    c = crit(stat, level, df)
    for r in rows:
        if r.kind == "eq":
            out.write(f"{abbrev(r.label):<{X}} |", "result")
            out.write(("  " + r.eq if r.eq else "") + "\n", "text")   # ex.: (base outcome)
            continue
        if r.kind == "header":
            out.write(f"{abbrev(r.label):>{X}} |\n", "text")
            continue
        if r.kind == "blank":
            out.write(f"{'':>{X}} |\n", "text")
            continue
        if r.kind == "sep":
            out.write("-" * W + "+" + "-" * 64 + "\n", "text")
            continue
        if r.tight:
            out.write(f"{abbrev(r.label, W):>{W}}|", "text")
        else:
            out.write(f"{abbrev(r.label):>{X}} |", "text")
        if r.kind == "base":
            out.write("          0  (base)\n", "result")
            continue
        if r.kind == "omitted":
            out.write("          0  (omitted)\n", "result")
            continue
        if r.kind == "empty":
            out.write("          0  (empty)\n", "result")
            continue
        if r.kind == "cons":
            # a parameter fixed by the model (svar: entries of A and B)
            out.write(f"  {g9(r.b):>9}  (constrained)\n", "result")
            continue
        if r.kind == "value":
            # a single number without standard error (rho in prais)
            out.write(f"  {g9(r.b):>9}\n", "result")
            continue
        if r.kind == "fixed":
            # exposure()/offset(): coefficient fixed at 1 (compat 0507)
            out.write(f"          1  {r.eq}\n", "result")
            continue
        b, se = r.b, r.se
        if r.kind == "aux":
            # ancillary parameters (/cut1, /sigma, alpha): no z and P>|z|; the interval may come
            # ready (transformed) in r.eq = "lo hi"
            if not (np.isfinite(se) and abs(se) < SYS):
                # missing standard error (alpha at the boundary, compat 0602)
                out.write(f"  {g9(b):>9}  {'.':>9}" + " " * 17 + f"    {'.':>9}   {'.':>9}\n", "result")
                continue
            if r.eq:
                lo, hi = (float(v) for v in r.eq.split())
            else:
                lo, hi = b - c * se, b + c * se
            out.write(f"  {g9(b):>9}  {g9(se):>9}" + " " * 17 + f"    {g9(lo):>9}   {g9(hi):>9}\n",
                      "result")
            continue
        if se > 0 and np.isfinite(se) and abs(se) < SYS:
            tval = b / se
            p = pvalue(stat, tval, df)
            lo, hi = b - c * se, b + c * se
            if r.zero_lo:       # variance: one-sided test against zero, interval truncated at zero
                p = float(st.norm.sf(tval))
                lo = max(lo, 0.0)
            if eform and not r.noeform:
                b, se, lo, hi = np.exp(b), np.exp(b) * se, np.exp(lo), np.exp(hi)
            out.write(f"  {g9(b):>9}  {g9(se):>9}  {tval:>7.2f}  {p:>6.3f}    {g9(lo):>9}   {g9(hi):>9}\n",
                      "result")
        else:
            # VERIFY: zero or missing standard error (constrained)
            out.write(f"  {g9(b):>9}  {g9(se) if se == 0 else '.':>9}  {'.':>7}  {'.':>6}    "
                      f"{'.':>9}   {'.':>9}\n", "result")
    if footer:
        out.write(LINE_ + "\n", "text")


def rows_from_columns(cols: "list[Column]", b: dict[str, float], se: dict[str, float],
                      ds, *, constant: bool, show_base: bool = False) -> list[CoefRow]:
    """Table rows: factor terms and variables with time-series operators
    get a header and a blank line before and after."""
    rows: list[CoefRow] = []
    i = 0
    n = len(cols)

    def blank_before():
        if rows and rows[-1].kind != "blank":
            rows.append(CoefRow("blank"))

    while i < n:
        c = cols[i]
        j = i + 1
        if c.group:
            while j < n and cols[j].group == c.group:
                j += 1
        group = cols[i:j]
        shown = [g for g in group if not g.base or show_base]
        interaction = "#" in c.name
        grouped = bool(c.group) and (len(shown) > 1 or (c.group[0] == "t" and interaction)
                                     or (c.group[0] == "v" and c.rowlab != "--."))
        if (not grouped and c.group and c.group[0] == "t" and not interaction and len(shown) == 1
                and _has_value_label(ds, shown[0])):
            # factor with a single level and a value label: header and the label, without a blank
            # line after (compat 0601: "sexo |" / "fem  |")
            blank_before()
            g = shown[0]
            rows.append(CoefRow("header", c.labels[0]))
            lab = _cell_label(ds, g) + " "
            if g.omitted:
                rows.append(CoefRow("omitted", lab, name=g.name))
            else:
                rows.append(CoefRow("coef", lab, b[g.name], se[g.name], name=g.name))
            i = j
            continue
        if not grouped:
            # a single line: "x", "1.h" (factor with one level), "c.x#c.x"; interactions go between
            # blank lines (compat 0502)
            if interaction:
                blank_before()
            for g in shown:
                lab = g.labels[0] if (g.labels and not g.group) or (g.group and g.group[0] == "v") \
                    else _strip_b(g.name)
                if g.base:
                    rows.append(CoefRow("base", lab, name=g.name))
                elif g.omitted:
                    rows.append(CoefRow("omitted", lab, name=g.name))
                else:
                    rows.append(CoefRow("coef", lab, b[g.name], se[g.name], name=g.name))
            if interaction:
                rows.append(CoefRow("blank"))
            i = j
            continue
        blank_before()
        rows.append(CoefRow("header", c.labels[0] if c.group[0] == "t" else c.group[1]))
        for g in group:
            lab = (_cell_label(ds, g) + " ") if g.group[0] == "t" else g.rowlab
            if g.base:
                if show_base:
                    rows.append(CoefRow("base", lab, name=g.name))
                continue
            if g.omitted:
                rows.append(CoefRow("omitted", lab, name=g.name))
                continue
            rows.append(CoefRow("coef", lab, b[g.name], se[g.name], name=g.name))
        rows.append(CoefRow("blank"))
        i = j
    if constant:
        rows.append(CoefRow("coef", "_cons", b["_cons"], se["_cons"], name="_cons"))
    elif rows and rows[-1].kind == "blank":
        rows.pop()
    return rows


def _has_value_label(ds, col) -> bool:
    """Does the factor level have a value label?"""
    try:
        var = col.name.split(".", 1)[1]
        v = ds.get(var)
        lv = col.cell[0] if col.cell else None
        return bool(v.value_label and lv is not None
                    and ds.value_labels.get(v.value_label, {}).get(int(lv)))
    except Exception:   # noqa: BLE001
        return False


def _head_label(head: str) -> str:
    return head


def _strip_b(name: str) -> str:
    import re
    return "#".join(re.sub(r"^(\d+)(?:b|bn|o)+\.", r"\1.", p) for p in name.split("#"))


def _cell_label(ds, col) -> str:
    """Level of a factor: its value label, if any (VERIFY width)."""
    from .fvars import _level_text
    parts = col.name.split("#")
    texts = []
    k = 0
    for p in parts:
        if p.startswith("c.") or "." not in p or not p[0].isdigit():
            continue
        var = p.split(".", 1)[1]
        lv = col.cell[k] if k < len(col.cell) else None
        k += 1
        txt = _level_text(lv) if lv is not None else ""
        try:
            v = ds.get(var)
            if v.value_label and lv is not None:
                lab = ds.value_labels.get(v.value_label, {}).get(int(lv))
                if lab:
                    txt = lab
        except Exception:   # noqa: BLE001
            pass
        texts.append(txt)
    return "#".join(texts) if texts else col.labels[1]
