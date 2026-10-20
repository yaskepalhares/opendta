"""Helpers shared by the time-series commands of phase 18 (tsfilter, tssmooth, tsappend, tsfill,
tsreport, rolling, pergram, xcorr, the unit-root tests...): the time grid of each panel, time
units, creating variables and adding observations."""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.dataset import Variable
from ..core.errors import StataError
from ..estimation.tsops import require_ts

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS

# periods in a year for the %t formats (Ravn-Uhlig scaling and the defaults of the filters)
PERIODS_PER_YEAR = {"d": 365.0, "w": 52.0, "m": 12.0, "q": 4.0, "h": 2.0, "y": 1.0}
UNIT_NAME = {"d": "daily", "w": "weekly", "m": "monthly", "q": "quarterly", "h": "halfyearly",
             "y": "yearly"}
TYPES = ("byte", "int", "long", "float", "double")


def unit_code(s: "Session") -> str:
    """Letter of the %t format of the time variable ('' when it has none)."""
    info = require_ts(s.data)
    fmt = s.data.get(info.tvar).fmt
    c = fmt[2:3] if fmt.startswith("%t") else ""
    return c if c in PERIODS_PER_YEAR else ""


def unit_name(s: "Session") -> str:
    return UNIT_NAME.get(unit_code(s), ".")


def to_nan(x) -> np.ndarray:
    """Float vector with NaN in place of Stata's missing values."""
    v = np.array(x, dtype=np.float64)
    v[v >= SYS] = np.nan
    return v


def from_nan(x: np.ndarray) -> np.ndarray:
    v = np.array(x, dtype=np.float64)
    v[np.isnan(v)] = SYS
    return v


def panels(s: "Session", mask: np.ndarray | None = None):
    """[(panel value or None, rows ordered by time, times)] for the observations of `mask` that
    have a time value (and a panel value)."""
    ds = s.data
    info = require_ts(ds)
    t = ds.get(info.tvar).data.astype(np.float64)
    ok = t < SYS
    if mask is not None:
        ok &= mask
    if info.panel:
        pv = ds.get(info.panel).data.astype(np.float64)
        ok &= pv < SYS
        keys = np.unique(pv[ok])
    else:
        pv, keys = None, [None]
    out = []
    for k in keys:
        sel = ok & (pv == k) if pv is not None else ok
        rows = np.flatnonzero(sel)
        order = np.argsort(t[rows], kind="stable")
        rows = rows[order]
        out.append((k, rows, t[rows]))
    return out


def new_variable(s: "Session", name: str, values: np.ndarray, *, vtype: str | None = None,
                 replace: bool = False, label: str = "") -> None:
    ds = s.data
    if ds.has(name):
        if not replace:
            raise StataError(110, f"{name} already defined")
        ds.drop_vars([name])
    ds.add(Variable(name, vtype or s.settings.get("type", "float"), from_nan(values), label=label))


def split_type(words: list[str]) -> tuple[str | None, list[str]]:
    if words and words[0] in TYPES:
        return words[0], words[1:]
    return None, words


def append_periods(s: "Session", entries: list[tuple[float | None, float]]) -> np.ndarray:
    """Adds one observation for each (panel value, time) and re-sorts by panel and time. Returns
    the positions of the new observations before the sort as a boolean mask over the new data
    (so callers can refer to them)."""
    ds = s.data
    info = require_ts(ds)
    n0 = ds.nobs
    k = len(entries)
    if k == 0:
        return np.zeros(n0, dtype=bool)
    ds.set_obs(n0 + k)
    rows = np.arange(n0, n0 + k)
    ds.set_numeric(ds.get(info.tvar), np.array([e[1] for e in entries], dtype=float), rows)
    if info.panel:
        ds.set_numeric(ds.get(info.panel), np.array([e[0] for e in entries], dtype=float), rows)
    flag = np.zeros(n0 + k, dtype=bool)
    flag[rows] = True
    keys = [x for x in (info.panel, info.tvar) if x]
    from ..core.sorting import sort_order
    order = sort_order(ds, keys)
    ds.reorder_obs(order)
    ds.sortlist = list(keys)
    return flag[order]


def kv(s: "Session", pairs: list) -> None:
    """r() from (name, value) pairs in the order of storage."""
    s.r = dict(reversed(pairs))
