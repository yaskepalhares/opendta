"""Immediate commands: cii, ttesti, prtesti and tabi ([R] ci, ttest, prtest,
tabulate twoway). They take the statistics (or the table) on the command
line itself, without using the data in memory."""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np
from scipy import stats as st

from ..core import missing as M
from ..core import stats as S
from ..core.errors import StataError
from ..lang.syntax import find_top, match_options
from .registry import command
from .summarize import _clopper, _level, _prtest_one, _ttest_one, _ttest_two, g9

if TYPE_CHECKING:
    from ..session import Session


def _split(text: str) -> tuple[list[str], str]:
    c = find_top(text, ",")
    body, opt = (text, "") if c == -1 else (text[:c], text[c + 1:])
    return body.split(), opt


def _nums(words: list[str]) -> list[float]:
    try:
        return [float(w) for w in words]
    except ValueError:
        raise StataError(198, "invalid syntax")


def _synthetic(n: int, mean: float, sd: float) -> np.ndarray:
    """Sample with exactly the given n, mean and standard deviation (to reuse
    the code of the tests that take the data)."""
    z = np.arange(n, dtype=float)
    z = (z - z.mean()) / z.std(ddof=1)
    return mean + sd * z


# ---------------------------------------------------------------------------
# cii
# ---------------------------------------------------------------------------

@command("cii", "cii")
def cmd_cii(s: "Session", args: str) -> None:
    words, opt = _split(args)
    kind = None
    if words and not _isnum(words[0]):
        head = words.pop(0)
        for k in ("means", "proportions", "variances"):
            if k.startswith(head) and len(head) >= 4 or head == k[:-1]:
                kind = k
                break
        if kind is None:
            raise StataError(198, "invalid syntax")
    vals = _nums(words)
    if kind is None:   # old syntax: 3 numbers = mean; 2 = proportion
        kind = "means" if len(vals) == 3 else "proportions"
    o = match_options(opt, {"level": 1, "exact": 2, "wald": 2, "wilson": 2, "agresti": 2,
                            "jeffreys": 2, "poisson": 1, "exposure": 3})
    level = _level(s, o)
    out = s.output
    if kind == "means":
        if len(vals) != 3:
            raise StataError(198, "invalid syntax")
        n, mean, sd = vals
        if n < 1 or n != int(n):
            raise StataError(198, "n must be a positive integer")   # VERIFY
        se = sd / math.sqrt(n)
        lo, hi = S.t_ci(mean, se, n - 1, level) if n > 1 else (M.SYSMISS, M.SYSMISS)
        out.write("\n", "text")
        out.write(f"    Variable |        Obs        Mean    Std. Err.       [{level:g}% Conf. Interval]\n", "text")
        out.write("-" * 13 + "+" + "-" * 63 + "\n", "text")
        out.write(" " * 12 + " |", "text")
        out.write(f" {n:>10,.0f}   {g9(mean)}   {g9(se)}    {g9(lo).rjust(12)}   {g9(hi).rjust(9)}\n", "result")
        s.r = {"ub": hi, "lb": lo, "se": se, "mean": mean, "N": n, "level": level}
        return
    if kind == "proportions" and o.get("poisson"):
        # cii #exposure #events, poisson (Stata 14.0 syntax): exact Poisson CI from the chi-squared
        # ([R] ci, Methods and formulas)
        if len(vals) != 2:
            raise StataError(198, "invalid syntax")
        E, k = vals
        rate = k / E
        se = math.sqrt(k) / E
        a = (100 - level) / 200
        lo = float(st.chi2.ppf(a, 2 * k)) / 2 / E if k > 0 else 0.0
        hi = float(st.chi2.ppf(1 - a, 2 * (k + 1))) / 2 / E
        out.write("\n" + " " * 48 + "-- Poisson  Exact --".rjust(30) + "\n", "text")   # VERIFY title
        out.write(f"    Variable |   Exposure        Mean    Std. Err.       [{level:g}% Conf. Interval]\n", "text")
        out.write("-" * 13 + "+" + "-" * 63 + "\n", "text")
        out.write(" " * 12 + " |", "text")
        out.write(f" {E:>10,.0f}   {g9(rate)}   {g9(se)}    {g9(lo).rjust(12)}   {g9(hi).rjust(9)}\n", "result")
        s.r = {"ub": hi, "lb": lo, "se": se, "mean": rate, "N": E, "level": level}
        return
    if kind == "proportions":
        if len(vals) != 2:
            raise StataError(198, "invalid syntax")
        n, k = vals
        if 0 < k < 1:          # proportion instead of count (VERIFY)
            k = round(k * n)
        ph = k / n
        se = math.sqrt(ph * (1 - ph) / n)
        method = next((m for m in ("wald", "wilson", "agresti", "jeffreys") if o.get(m)), "exact")
        lo, hi = _prop_ci(k, n, level, method)
        title = {"exact": "-- Binomial Exact --", "wald": "-- Binomial Wald ---",
                 "wilson": "------ Wilson ------", "agresti": "-- Agresti-Coull ---",
                 "jeffreys": "----- Jeffreys -----"}[method]
        out.write("\n" + " " * 48 + title.rjust(30) + "\n", "text")
        out.write(f"    Variable |        Obs        Mean    Std. Err.       [{level:g}% Conf. Interval]\n", "text")
        out.write("-" * 13 + "+" + "-" * 63 + "\n", "text")
        out.write(" " * 12 + " |", "text")
        out.write(f" {n:>10,.0f}   {g9(ph)}   {g9(se)}    {g9(lo).rjust(12)}   {g9(hi).rjust(9)}\n", "result")
        s.r = {"ub": hi, "lb": lo, "se": se, "proportion": ph, "N": n, "level": level}
        return
    raise StataError(198, "cii variances not implemented")   # VERIFY


def _isnum(t: str) -> bool:
    try:
        float(t)
        return True
    except ValueError:
        return False


def _prop_ci(k: float, n: float, level: float, method: str) -> tuple[float, float]:
    """[R] ci intervals for proportions."""
    z = S.z_crit(level)
    ph = k / n
    if method == "exact":
        return _clopper(k, n, level)
    if method == "wald":
        se = math.sqrt(ph * (1 - ph) / n)
        return ph - z * se, ph + z * se
    if method == "wilson":
        c = (k + z * z / 2) / (n + z * z)
        h = z * math.sqrt(n) / (n + z * z) * math.sqrt(ph * (1 - ph) + z * z / (4 * n))
        return c - h, c + h
    if method == "agresti":
        nt = n + z * z
        pt = (k + z * z / 2) / nt
        h = z * math.sqrt(pt * (1 - pt) / nt)
        return max(0.0, pt - h), min(1.0, pt + h)
    from scipy import stats as st
    a = 1 - level / 100
    lo = 0.0 if k == 0 else float(st.beta.ppf(a / 2, k + 0.5, n - k + 0.5))
    hi = 1.0 if k == n else float(st.beta.ppf(1 - a / 2, k + 0.5, n - k + 0.5))
    return lo, hi


# ---------------------------------------------------------------------------
# ttesti e prtesti
# ---------------------------------------------------------------------------

@command("ttesti", "ttesti")
def cmd_ttesti(s: "Session", args: str) -> None:
    words, opt = _split(args)
    vals = _nums(words)
    o = match_options(opt, {"level": 1, "unequal": 3, "welch": 1})
    level = _level(s, o)
    if len(vals) == 4:
        n, mean, sd, mu = vals
        _ttest_one(s, "x", int(n), mean, sd, mu, level)
        return
    if len(vals) == 6:
        n1, m1, sd1, n2, m2, sd2 = vals
        if min(n1, n2) < 2:
            raise StataError(2001, "insufficient observations")
        x1, x2 = _synthetic(int(n1), m1, sd1), _synthetic(int(n2), m2, sd2)
        _ttest_two(s, "x", "", x1, x2, level, bool(o.get("unequal")), bool(o.get("welch")), names=("x", "y"))
        return
    raise StataError(198, "invalid syntax")


@command("prtesti", "prtesti")
def cmd_prtesti(s: "Session", args: str) -> None:
    words, opt = _split(args)
    vals = _nums(words)
    o = match_options(opt, {"level": 1, "count": 1})
    level = _level(s, o)
    if len(vals) == 3:
        n, p1, p0 = vals
        if o.get("count"):
            p1 = p1 / n
        _prtest_one(s, "x", int(n), p1, p0, level)
        return
    raise StataError(198, "two-sample prtesti not implemented")   # VERIFY


# ---------------------------------------------------------------------------
# tabi
# ---------------------------------------------------------------------------

@command("tabi", "tabi")
def cmd_tabi(s: "Session", args: str) -> None:
    from ..core.dataset import Dataset, Variable
    c = find_top(args, ",")
    body, opt = (args, "") if c == -1 else (args[:c], args[c + 1:])
    rows = [r.split() for r in body.split("\\")]
    if not rows or any(not r for r in rows):
        raise StataError(198, "invalid syntax")
    ncol = len(rows[0])
    if any(len(r) != ncol for r in rows):
        raise StataError(198, "rows have unequal number of columns")   # VERIFY
    vals = [_nums(r) for r in rows]
    if any(v < 0 or v != int(v) for r in vals for v in r):
        raise StataError(198, "may not use noninteger frequency weights")   # VERIFY
    opts = [w for w in opt.split()]
    keep = any(w.startswith("rep") for w in opts)
    opts = [w for w in opts if not w.startswith("rep")]
    ds = Dataset()
    cells = [(i + 1, j + 1, v) for i, r in enumerate(vals) for j, v in enumerate(r)]
    ds.set_obs(len(cells))
    ds.add(Variable("row", "byte", np.array([c[0] for c in cells], dtype=float)))
    ds.add(Variable("col", "byte", np.array([c[1] for c in cells], dtype=float)))
    ds.add(Variable("pop", "long", np.array([c[2] for c in cells], dtype=float)))
    saved = s.data
    s.data = ds
    try:
        tail = (", " + " ".join(opts)) if opts else ""
        s.interp._execute_expanded(f"tabulate row col [fweight=pop]{tail}")
    finally:
        if not keep:
            s.data = saved
        else:
            s.notify_state()
