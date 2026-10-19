"""Business calendars ([D] bcal, [D] datetime business calendars): the .stbcal
files, the %tb display format and the bofd()/dofb() conversions.

A calendar lists the days of a range that are business days; business date 0 is
the center date, and the others are numbered consecutively over the business
days only. Supported directives: version, purpose, dateformat, range,
centerdate, omit dayofweek, omit date (a year may be *), omit dowinmonth and
include date. The "and"/"if" extensions of omit are not supported yet.
Calendars are found in the current directory and the folders registered with
add_path (the adopath)."""

from __future__ import annotations

import datetime as _dt
import re
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

from .errors import StataError

EPOCH = _dt.date(1960, 1, 1)
_MONTHS = ["jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"]
_DOW = {"su": 6, "mo": 0, "tu": 1, "we": 2, "th": 3, "fr": 4, "sa": 5}

_PATHS: list[Path] = []
_CACHE: dict[str, "Calendar"] = {}


@dataclass
class Calendar:
    name: str
    purpose: str = ""
    lo: int = 0                 # range, in %td units
    hi: int = 0
    center: int = 0
    days: np.ndarray = field(default_factory=lambda: np.zeros(0, dtype=np.int64))   # business days (%td)
    path: str = ""

    def bofd(self, td: float) -> float:
        i = int(np.searchsorted(self.days, td))
        if i < len(self.days) and self.days[i] == td:
            return float(i - self.c_index)
        return float("nan")

    def dofb(self, b: float) -> float:
        i = int(b) + self.c_index
        if b == int(b) and 0 <= i < len(self.days):
            return float(self.days[i])
        return float("nan")

    @property
    def c_index(self) -> int:
        return int(np.searchsorted(self.days, self.center))


def add_path(p: str | Path) -> None:
    p = Path(p)
    if p not in _PATHS:
        _PATHS.append(p)


def clear_cache() -> None:
    _CACHE.clear()


def find(name: str) -> Path | None:
    for d in [Path.cwd()] + _PATHS:
        f = d / f"{name}.stbcal"
        if f.exists():
            return f
    return None


def get(name: str) -> Calendar:
    cal = _CACHE.get(name)
    if cal is None:
        f = find(name)
        if f is None:
            raise StataError(601, f"business calendar {name} (file {name}.stbcal) not found")
        cal = parse(name, f.read_text(encoding="utf-8"), str(f))
        _CACHE[name] = cal
    return cal


def _date(tok: str, order: str, year: int | None = None) -> _dt.date:
    """A date written as in the calendar's dateformat: 2011nov01, 01nov2011, 2011-11-01...
    With year given, '*' in the token stands for it."""
    t = tok.strip().lower()
    parts = re.findall(r"[a-z]+|\d+|\*", t)
    if len(parts) != 3:
        raise StataError(198, f"{tok}: invalid date")
    vals = dict(zip(order, parts))
    y = vals["y"]
    y = year if y == "*" else int(y)
    m = vals["m"]
    m = _MONTHS.index(m[:3]) + 1 if m.isalpha() else int(m)
    d = int(vals["d"])
    return _dt.date(y, m, d)


def _td(d: _dt.date) -> int:
    return (d - EPOCH).days


def parse(name: str, text: str, path: str = "") -> Calendar:
    cal = Calendar(name, path=path)
    order = "ymd"
    omit_dow: set[int] = set()
    omit_dates: set[int] = set()
    yearly: list[str] = []
    dowinmonth: list[tuple[int, int, int]] = []     # (k, weekday, month)
    include: set[int] = set()
    center_tok = None
    rng = None
    for raw in text.splitlines():
        line = raw.split("//")[0].strip()
        if not line or line.startswith("*"):
            continue
        words = line.split()
        key = words[0].lower()
        if key == "version":
            continue
        if key == "purpose":
            cal.purpose = line[len(words[0]):].strip().strip('"')
        elif key == "dateformat":
            order = words[1].lower()
            if sorted(order) != ["d", "m", "y"]:
                raise StataError(198, f"dateformat {words[1]} invalid")
        elif key == "range":
            rng = (words[1], words[2])
        elif key == "centerdate":
            center_tok = words[1]
        elif key in ("omit", "include"):
            what = words[1].lower()
            rest = words[2:]
            if any(w.lower() in ("and", "if") for w in rest):
                raise StataError(198, "omit ... and/if is not supported yet in OpenDTA")
            if key == "omit" and what == "dayofweek":
                for w in " ".join(rest).strip("()").split():
                    omit_dow.add(_DOW[w.lower()[:2]])
            elif what == "date":
                tok = rest[0]
                if "*" in tok:
                    if key == "include":
                        raise StataError(198, "include date with * is not supported")   # VERIFY
                    yearly.append(tok)
                else:
                    target = include if key == "include" else omit_dates
                    target.add(_td(_date(tok, order)))
            elif key == "omit" and what == "dowinmonth":
                # omit dowinmonth +4 Th of Nov
                k = int(rest[0])
                wd = _DOW[rest[1].lower()[:2]]
                mon = rest[3].lower()
                m = _MONTHS.index(mon[:3]) + 1 if mon.isalpha() else int(mon)
                dowinmonth.append((k, wd, m))
            else:
                raise StataError(198, f"{line}: invalid directive")
        else:
            raise StataError(198, f"{line}: invalid directive")
    if rng is None:
        raise StataError(198, "business calendar has no range")
    lo, hi = _date(rng[0], order), _date(rng[1], order)
    cal.lo, cal.hi = _td(lo), _td(hi)
    days = []
    d = lo
    while d <= hi:
        t = _td(d)
        off = d.weekday() in omit_dow or t in omit_dates
        if not off:
            for tok in yearly:
                try:
                    if _date(tok, order, d.year) == d:
                        off = True
                        break
                except ValueError:
                    pass
        if not off:
            for k, wd, m in dowinmonth:
                if d.month == m and d.weekday() == wd:
                    if k > 0 and (d.day - 1) // 7 + 1 == k:
                        off = True
                    elif k < 0:
                        last = (_dt.date(d.year + (m == 12), m % 12 + 1, 1) - _dt.timedelta(days=1)).day
                        if (last - d.day) // 7 + 1 == -k:
                            off = True
        if t in include:
            off = False
        if not off:
            days.append(t)
        d += _dt.timedelta(days=1)
    cal.days = np.array(days, dtype=np.int64)
    cal.center = _td(_date(center_tok, order)) if center_tok else cal.lo
    if cal.center not in set(days):
        raise StataError(198, "the center date is not a business day")   # VERIFY
    return cal
