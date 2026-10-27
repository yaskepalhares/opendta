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
from .summarize import g9, numeric_vars, weights

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
    if o.get("equal"):
        # equal: tied values get the same (highest) cumulative value; otherwise each observation
        # of a tie has its own, in the order of the data
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

_SYMS = {1: ["*"], 2: ["*", "."], 5: ["*", "t", "f", "s", "."]}


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


def _natural_unit(x: np.ndarray) -> int:
    """Largest 10**u that divides every value (u <= 0 for values with decimals)."""
    u = 0
    # float variables hold 2.01 as 2.0099999: compare with a tolerance relative to the values,
    # or the unit falls to 1e-8 and the plot has hundreds of millions of rows
    while u > -8 and not np.all(np.abs(x - np.round(x / 10.0 ** u) * 10.0 ** u)
                                <= 1e-6 * np.maximum(np.abs(x), 10.0 ** u)):
        u -= 1
    return u


def _stem_rows(q: np.ndarray, d: int, L: int):
    """(rows, leaves-by-row) of a stem-and-leaf plot of integers q with d digits per leaf and L
    lines per stem."""
    B = 10 ** d
    w = B // L
    mag = np.abs(q) // w
    keys = [(bool(a < 0), int(b)) for a, b in zip(q, mag)]
    return keys, w


@command("stem", "stem", byable=True)
def cmd_stem(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"digits": 1, "lines": 1, "round": 1, "prune": 1, "width": 1}) \
        if p.options.strip() else {}
    if isinstance(o.get("width"), str):
        # width(#): the width of a stem, another way of giving lines() = 10 / width
        # (VERIFY: widths that do not divide 10)
        if isinstance(o.get("lines"), str):
            raise StataError(198, "options lines() and width() may not be combined")
        try:
            wv = int(o["width"])
        except ValueError:
            raise StataError(198, "option width() incorrectly specified")
        if wv not in (2, 5, 10):
            raise StataError(198, "width() must be 2, 5 or 10")
        o["lines"] = str(10 // wv)
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
    if not isinstance(o.get("digits"), str) and not isinstance(o.get("round"), str):
        return _stem_natural(s, v, x, o, target)
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


def _stem_natural(s: "Session", v, x: np.ndarray, o: dict, target: float) -> None:
    """stem without digits() or round(): the leaves keep every digit of the data. lines() given:
    one digit per leaf; otherwise the stem unit and the number of lines come from the number of
    rows nearest to 10*log10(n). VERIFY: the choice of the unit (compat/do/1903_base_stats.do)."""
    u = _natural_unit(x)
    q = np.round(x / 10.0 ** u).astype(np.int64)
    if isinstance(o.get("lines"), str):
        L = int(o["lines"])
        if L not in _SYMS:
            raise StataError(198, "lines() must be 1, 2 or 5")
        d = 1
    else:
        top = int(np.abs(q).max())
        cands = []
        for d_ in range(1, len(str(top)) + 2):
            for L_ in ((1, 2, 5) if d_ == 1 else (1, 2)):
                keys, _w = _stem_rows(q, d_, L_)
                rows = _span(keys, L_)
                cands.append((abs(rows - target) + (0.001 if rows > target else 0), d_, L_))
        cands.sort()
        _, d, L = cands[0]
    keys, w = _stem_rows(q, d, L)
    B = 10 ** d
    neg, pos = _order(keys)
    if neg and pos:
        full = [(True, m) for m in range(neg[0][1], -1, -1)] + [(False, m) for m in range(0, pos[-1][1] + 1)]
    elif neg:
        full = [(True, m) for m in range(neg[0][1], neg[-1][1] - 1, -1)]
    else:
        full = [(False, m) for m in range(pos[0][1], pos[-1][1] + 1)]
    leaves: dict = {k: [] for k in full}
    qs = np.sort(q)
    for qq, k in zip(qs, _stem_rows(qs, d, L)[0]):
        leaves[k].append(f"{abs(int(qq)) % B:0{d}d}")
    labels = []
    for neg_, m in full:
        stem = (m * w) // B
        if d == 1:
            sym = _SYMS[L][((m * w) % B) // w]
        else:
            sym = "*" * d                      # VERIFY: label of stems with multi-digit leaves
        labels.append((("-" if neg_ else "") + str(stem) + sym, (",".join if d > 1 else "".join)(leaves[(neg_, m)])))
    if o.get("prune"):
        labels = [t for t in labels if t[1]]
    out = s.output
    head = f"Stem-and-leaf plot for {v.name}"
    if v.label:
        head += f" ({v.label})"
    out.write("\n" + head + "\n\n", "text")
    if u:
        # also for decimal units: "rounded to nearest multiple of .01" (compat 1913)
        unit = f"{10.0 ** u:g}".removeprefix("0")
        out.write(f"{v.name} rounded to nearest multiple of {unit}\nplot in units of {unit}\n\n", "text")
    wd = max(len(a) for a, _b in labels)
    for lab, lv in labels:
        out.write(f"{lab:>{wd + 2}} | ", "text")
        out.write(lv + "\n", "result")
    out.write("\n", "text")


# ---------------------------------------------------------------------------
# lv
# ---------------------------------------------------------------------------

_LETTERS = "MFEDCBAZYXWVU"


@command("lv", "lv", byable=True)
def cmd_lv(s: "Session", args: str) -> None:
    p = parse_standard(args)
    o = match_options(p.options, {"generate": 1, "tail": 1}) if p.options.strip() else {}
    tail_max = None
    if isinstance(o.get("tail"), str):
        try:
            tail_max = int(o["tail"])
        except ValueError:
            tail_max = 0
        if tail_max < 1:
            raise StataError(198, "option tail() incorrectly specified")         # VERIFY
    vars_ = numeric_vars(s, p.varlist, allow_string=False)
    if o.get("generate"):
        if len(vars_) > 1:
            raise StataError(103, "too many variables specified")             # VERIFY
        for nm in ("_mid", "_psigma"):
            if s.data.has(nm):
                s.data.drop_vars([nm])          # generate replaces them (compat 1913)
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
        if tail_max is not None:
            # the letter of 1/2^(i+1) is shown while 2^i < tail(), then the extremes (compat 1913:
            # tail(3) and tail(4) show M and F, tail(8) down to E)
            keep = [i for i in range(len(depths) - 1) if 2 ** i < tail_max]
            depths = [depths[i] for i in keep] + [1.0]
        rows = []
        for i, d in enumerate(depths):
            lo = 0.5 * (xs[int(math.floor(d)) - 1] + xs[int(math.ceil(d)) - 1])
            hi = 0.5 * (xs[n - int(math.floor(d))] + xs[n - int(math.ceil(d))])
            rows.append((_LETTERS[i] if i < len(depths) - 1 else " ", d, lo, hi))
        out = s.output
        ttl = v.label or v.name
        head = [" "] * 70
        head[1] = "#"
        txt = f"{n}"
        head[7:7 + len(txt)] = txt
        pos = (59 - len(ttl) + 1) // 2
        head[pos:pos + len(ttl)] = ttl
        out.write("\n" + "".join(head).rstrip() + "\n", "text")
        out.write(" " * 13 + "-" * 33 + "\n", "text")
        hspread = None
        for i, (let, d, lo, hi) in enumerate(rows):
            ds = f"{d:g}"
            ip, dp = ds.split(".") if "." in ds else (ds, "")
            left = f" {let}".ljust(9 - len(ip)) + ip + ("." + dp if dp else "")
            out.write(f"{left:<12}|", "text")
            mid = 0.5 * (lo + hi)
            if i == 0:
                out.write(f"{'':>10}{g9(mid):>11}{'':>11} |", "result")
                out.write(f"{'spread':>10}{'pseudosigma':>13}\n", "text")
            else:
                spread = hi - lo
                # VERIFY: the pseudosigma of the last row (depth 1) is not given by the same tail
                tail = (d - 1 / 3) / (n + 1 / 3)
                pseudo = spread / (2 * float(st.norm.ppf(1 - tail))) if tail < 0.5 else math.nan
                out.write(f"{g9(lo):>10}{g9(mid):>11}{g9(hi):>11} |", "result")
                out.write(f"{g9(spread):>10}{g9(pseudo):>12}\n", "result")
                if let == "F":
                    hspread = spread
        if o.get("generate"):
            # _mid and _psigma: one observation per letter value, the extremes left out (compat 1913)
            rows_g = rows[:-1]
            k = len(rows_g)
            if s.data.nobs < k:
                raise StataError(198, "too few observations for generate")      # VERIFY
            mids, psig = np.full(s.data.nobs, SYS), np.full(s.data.nobs, SYS)
            for i, (_let, d, lo, hi) in enumerate(rows_g):
                mids[i] = 0.5 * (lo + hi)
                if i:
                    tl = (d - 1 / 3) / (n + 1 / 3)
                    if tl < 0.5:
                        psig[i] = (hi - lo) / (2 * float(st.norm.ppf(1 - tl)))
            s.data.add(Variable("_mid", "float", mids))
            s.data.add(Variable("_psigma", "float", psig))
        out.write(" " * 12 + "|" + " " * 33 + "|\n", "text")
        out.write(" " * 12 + "|" + " " * 33 + "|" + f"{'# below':>10}{'# above':>12}\n", "text")
        if hspread is not None:
            hl, hu = rows[1][2], rows[1][3]
            for lab, k in (("inner fence", 1.5), ("outer fence", 3.0)):
                flo, fhi = hl - k * hspread, hu + k * hspread
                nb, na = int((xs < flo).sum()), int((xs > fhi).sum())
                out.write(f"{lab:<12}|{g9(flo):>10}{'':>11}{g9(fhi):>11} |", "text")
                out.write(f"{nb:>10}{na:>12}\n", "result")


# ---------------------------------------------------------------------------
# smooth
# ---------------------------------------------------------------------------

@command("smooth", "sm")
def cmd_smooth(s: "Session", args: str) -> None:
    # smoother[,twice] varname: the comma and twice are glued to the smoother, as in Stata
    m = re.match(r"\s*(\w+)(,twice\b)?\s*,?\s*(.*)$", args, re.S)
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
