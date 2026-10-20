"""Displays and variable builders of the base-statistics block: cumul ([R] cumul), stem
([R] stem), lv ([R] lv) and smooth ([R] smooth). The smoothers are in opendta.core.basestats.
Layouts the manual does not fully show are marked VERIFY (compat/do/1903_base_stats.do)."""

from __future__ import annotations

import math
import re
from typing import TYPE_CHECKING

import numpy as np

from ..core import basestats as B
from ..core import missing as M
from ..core import stats as S
from ..core.dataset import Variable
from ..core.errors import StataError
from ..lang.syntax import match_options, parse_standard
from ..lazy import lazy
from ._util import touse
from .nptests import _set_r
from .registry import command
from .summarize import numeric_vars, weights

st = lazy("scipy.stats")

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


# ---------------------------------------------------------------------------
# cumul
# ---------------------------------------------------------------------------

@command("cumul", "cumul")
def cmd_cumul(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"generate": 1, "freq": 1, "equal": 1}) if p.options.strip() else {}
    if not isinstance(o.get("generate"), str):
        raise StataError(198, "option generate() required")
    new = o["generate"]
    if s.data.has(new):
        raise StataError(110, f"{new} already defined")
    v = numeric_vars(s, p.varlist, allow_string=False)[0]
    mask = touse(s, p) & S.valid(v.data)
    w, wtype, mask = weights(s, p, mask, ("fweight", "aweight", "pweight"))
    x = v.data[mask]
    ww = w[mask] if w is not None else np.ones(len(x))
    if len(x) == 0:
        raise StataError(2000, "no observations")
    order = np.argsort(x, kind="stable")
    xs, wsrt = x[order], ww[order]
    cum = np.cumsum(wsrt)
    if not o.get("equal"):
        # tied values get the highest cumulative value
        last = np.searchsorted(xs, xs, side="right") - 1
        cum = cum[last]
    total = float(ww.sum())
    out = cum if o.get("freq") else cum / total
    col = np.full(s.data.nobs, SYS)
    col[np.flatnonzero(mask)[order]] = out
    s.data.add(Variable(new, "float", col))


# ---------------------------------------------------------------------------
# stem
# ---------------------------------------------------------------------------

_SYMS = {1: [""], 2: ["*", "."], 5: ["*", "t", "f", "s", "."]}


def _rows(q: np.ndarray, L: int):
    """Row key (negative?, magnitude row) of each integer q for L lines per stem."""
    w = 10 // L
    mag = np.abs(q) // w
    return [(bool(a < 0), int(b)) for a, b in zip(q, mag)]


def _order(keys):
    """Keys from the most negative to the largest positive."""
    neg = sorted({k for k in keys if k[0]}, key=lambda k: -k[1])
    pos = sorted({k for k in keys if not k[0]}, key=lambda k: k[1])
    return neg, pos


def _span(keys, L: int) -> int:
    neg, pos = _order(keys)
    if neg and pos:
        return neg[0][1] + pos[-1][1] + 2
    if neg:
        return neg[0][1] - neg[-1][1] + 1
    return pos[-1][1] - pos[0][1] + 1


@command("stem", "stem", byable=True)
def cmd_stem(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"digits": 1, "lines": 1, "round": 1, "prune": 1}) \
        if p.options.strip() else {}
    v = numeric_vars(s, p.varlist, allow_string=False)[0]
    mask = touse(s, p) & S.valid(v.data)
    x = v.data[mask]
    n = len(x)
    if n == 0:
        raise StataError(2000, "no observations")
    if isinstance(o.get("round"), str):
        r = float(o["round"])
        x = np.round(x / r) * r
    target = 10 * math.log10(n)
    best = None
    cands = []
    if isinstance(o.get("digits"), str):
        es = [-int(o["digits"])]
    else:
        spread = float(x.max() - x.min()) or 1.0
        e0 = int(math.floor(math.log10(spread)))
        es = range(e0 - 3, e0 + 2)
    Ls = [int(o["lines"])] if isinstance(o.get("lines"), str) else [1, 2, 5]
    for e in es:
        q = np.round(x / 10.0 ** e).astype(np.int64)
        for L in Ls:
            if L not in _SYMS:
                raise StataError(198, "lines() must be 1, 2 or 5")
            rows = _span(_rows(q, L), L)
            cands.append((abs(rows - target) + (0.001 if rows > target else 0), e, L, rows))
    cands.sort()
    _, e, L, rows = cands[0]
    q = np.round(x / 10.0 ** e).astype(np.int64)
    q.sort()
    keys = _rows(q, L)
    neg, pos = _order(keys)
    if neg and pos:
        full = [(True, m) for m in range(neg[0][1], -1, -1)] + [(False, m) for m in range(0, pos[-1][1] + 1)]
    elif neg:
        full = [(True, m) for m in range(neg[0][1], neg[-1][1] - 1, -1)]
    else:
        full = [(False, m) for m in range(pos[0][1], pos[-1][1] + 1)]
    w = 10 // L
    leaves: dict = {k: [] for k in full}
    for qq, k in zip(q, keys):
        leaves[k].append(str(abs(int(qq)) % 10))
    labels = []
    for neg_, m in full:
        stem = (m * w) // 10
        sym = _SYMS[L][((m * w) % 10) // w]
        labels.append((("-" if neg_ else "") + str(stem) + sym, "".join(leaves[(neg_, m)])))
    if o.get("prune"):
        labels = [t for t in labels if t[1]]
    out = s.output
    head = f"Stem-and-leaf plot for {v.name}"
    if v.label:
        head += f" ({v.label})"
    out.write("\n" + head + "\n\n", "text")
    if e > 0:
        unit = f"{10 ** e:g}"
        out.write(f"{v.name} rounded to nearest multiple of {unit}\nplot in units of {unit}\n\n", "text")
    elif e < 0:
        out.write(f"{v.name} rounded to {-e} decimal place{'s' if -e > 1 else ''}\n"
                  f"plot in units of {10.0 ** e:g}\n\n", "text")     # VERIFY
    wd = max(len(a) for a, _b in labels)
    for lab, lv in labels:
        out.write(f"{lab:>{wd + 2}} | ", "text")
        out.write(lv + "\n", "result")
    out.write("\n", "text")


# ---------------------------------------------------------------------------
# lv
# ---------------------------------------------------------------------------

_LETTERS = "MHEDCBAZYXWVU"


@command("lv", "lv", byable=True)
def cmd_lv(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"generate": 1, "tail": 1}) if p.options.strip() else {}
    if o:
        raise StataError(198, f"option {next(iter(o))} is not available yet in OpenDTA")
    vars_ = numeric_vars(s, p.varlist, allow_string=False)
    mask = touse(s, p)
    for v in vars_:
        x = v.data[mask & S.valid(v.data)]
        n = len(x)
        if n == 0:
            raise StataError(2000, "no observations")
        xs = np.sort(x)
        depths = [(n + 1) / 2]
        while depths[-1] > 1:
            depths.append((math.floor(depths[-1]) + 1) / 2)
        rows = []
        for i, d in enumerate(depths):
            lo = 0.5 * (xs[int(math.floor(d)) - 1] + xs[int(math.ceil(d)) - 1])
            hi = 0.5 * (xs[n - int(math.floor(d))] + xs[n - int(math.ceil(d))])
            rows.append((_LETTERS[i], d, lo, hi))
        out = s.output
        out.write("\n" + f"{'#':>33} {n:<9}" + (f"{v.label}" if v.label else "") + "\n", "text")
        hspread = None
        for i, (let, d, lo, hi) in enumerate(rows):
            mid = 0.5 * (lo + hi)
            ds = f"{d:g}"
            out.write(f"{let:>4}{ds:>9}   |", "text")
            if i == 0:
                out.write(f"{'':>14}{mid:>10.7g}{'':>14} |\n", "result")
            else:
                spread = hi - lo
                tail = (d - 1 / 3) / (n + 1 / 3)
                pseudo = spread / (2 * float(st.norm.ppf(1 - tail))) if tail < 0.5 else math.nan
                out.write(f"{lo:>14.7g}{mid:>10.7g}{hi:>14.7g} |{spread:>10.7g}{pseudo:>10.4g}\n", "result")
                if let == "H":
                    hspread = spread
        if hspread is not None:
            hl, hu = rows[1][2], rows[1][3]
            out.write(f"{'':>13}{'inner fence':>14}|{hl - 1.5 * hspread:>14.7g}{'':>10}{hu + 1.5 * hspread:>14.7g}\n", "text")
            out.write(f"{'':>13}{'outer fence':>14}|{hl - 3 * hspread:>14.7g}{'':>10}{hu + 3 * hspread:>14.7g}\n", "text")
        out.write("\n", "text")


# ---------------------------------------------------------------------------
# smooth
# ---------------------------------------------------------------------------

@command("smooth", "sm")
def cmd_smooth(s: "Session", args: str) -> None:
    m = re.match(r"\s*(\w+)\s*(,\s*twice\b)?\s*(.*)$", args, re.S)
    if not m:
        raise StataError(198, "invalid syntax")
    spec, twice, rest = m.group(1), bool(m.group(2)), m.group(3)
    p = parse_standard(rest)
    o = match_options(p.options, {"generate": 1, "twice": 1}) if p.options.strip() else {}
    twice = twice or bool(o.get("twice"))
    if not isinstance(o.get("generate"), str):
        raise StataError(198, "option generate() required")
    new = o["generate"]
    if s.data.has(new):
        raise StataError(110, f"{new} already defined")
    v = numeric_vars(s, p.varlist, allow_string=False)[0]
    mask = touse(s, p) & S.valid(v.data)
    x = v.data[mask]
    if len(x) == 0:
        raise StataError(2000, "no observations")
    try:
        y = B.tukey_smooth(x, spec, twice)
    except ValueError as exc:
        raise StataError(198, f"invalid smoother: {exc}")
    col = np.full(s.data.nobs, SYS)
    col[mask] = y
    s.data.add(Variable(new, "float", col))
