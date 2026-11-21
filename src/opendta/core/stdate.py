"""Stata dates and times as numbers (manual [D] datetime): elapsed dates (%td), months, quarters,
half-years and weeks since 1960, milliseconds since 01jan1960 00:00:00 without (%tc) and with
(%tC) leap seconds, and the reading of date strings with a mask (date(), clock(), monthly()...).

Plain functions of floats and strings; missing values are handled by the caller."""

from __future__ import annotations

import datetime as _dt
import math
import re

EPOCH = _dt.date(1960, 1, 1)
MS_DAY = 86_400_000.0

# Days after which a leap second was added (23:59:60 UTC), from the IERS bulletins.
_LEAP_DAYS = [
    (1972, 6, 30), (1972, 12, 31), (1973, 12, 31), (1974, 12, 31), (1975, 12, 31),
    (1976, 12, 31), (1977, 12, 31), (1978, 12, 31), (1979, 12, 31), (1981, 6, 30),
    (1982, 6, 30), (1983, 6, 30), (1985, 6, 30), (1987, 12, 31), (1989, 12, 31),
    (1990, 12, 31), (1992, 6, 30), (1993, 6, 30), (1994, 6, 30), (1995, 12, 31),
    (1997, 6, 30), (1998, 12, 31), (2005, 12, 31), (2008, 12, 31), (2012, 6, 30),
    (2015, 6, 30), (2016, 12, 31),
]
# %tc value of the midnight that follows each leap second
_LEAP_TC = [((_dt.date(*d) - EPOCH).days + 1) * MS_DAY for d in _LEAP_DAYS]

_MONTHS = ["jan", "feb", "mar", "apr", "may", "jun", "jul", "aug", "sep", "oct", "nov", "dec"]


def mdy(m: float, d: float, y: float) -> float | None:
    try:
        return float((_dt.date(int(y), int(m), int(d)) - EPOCH).days)
    except (ValueError, OverflowError):
        return None


def to_date(td: float) -> _dt.date | None:
    try:
        return EPOCH + _dt.timedelta(days=math.floor(td))
    except OverflowError:
        return None


def week(d: _dt.date) -> int:
    """Stata week of the year: 1 to 52, the 52nd week has 8 or 9 days."""
    return min((d.timetuple().tm_yday - 1) // 7 + 1, 52)


def dofw(w: float) -> float | None:
    y, k = 1960 + math.floor(w / 52), int(w) % 52
    start = mdy(1, 1, y)
    return None if start is None else start + 7 * k


def dofperiod(p: float, per_year: int) -> float | None:
    """First day of period p of a %tm (12), %tq (4) or %th (2) calendar."""
    y, k = 1960 + math.floor(p / per_year), int(p) % per_year
    return mdy(1 + k * (12 // per_year), 1, y)


def hms(h: float, m: float, s: float) -> float:
    return 1000.0 * (3600.0 * h + 60.0 * m + s)


def leap_before(c: float) -> int:
    """Leap seconds added before the %tc time c."""
    return sum(1 for b in _LEAP_TC if c >= b)


def Cofc(c: float) -> float:
    return c + 1000.0 * leap_before(c)


def cofC(C: float) -> float:
    """%tc time of the %tC time C. VERIFY: a time inside a leap second (23:59:60) maps to
    23:59:59 plus the fraction."""
    n = 0
    for i, b in enumerate(_LEAP_TC):
        if C >= b + 1000.0 * i:
            n = i + 1
    return C - 1000.0 * n


# -- strings with a mask -----------------------------------------------------------------------

_DIGITS = re.compile(r"\d+")


def parse(text: str, mask: str, topyear: float | None = None) -> dict[str, float] | None:
    """Reads text following mask (letters D M Y h m s, plus Q H W for periods, # to skip one
    element, digits before Y as the century of two-digit years). Returns the fields found or
    None when text does not follow the mask."""
    fields = [c for c in mask if c in "DMYhmsQHW#" or c.isdigit()]
    century = "".join(c for c in fields if c.isdigit())
    fields = [c for c in fields if not c.isdigit()]
    t = text.strip().lower()
    pos, out = 0, {}
    for k, f in enumerate(fields):
        while pos < len(t) and not t[pos].isalnum():
            pos += 1
        if f == "M" and pos < len(t) and t[pos].isalpha():
            word = re.match(r"[a-z]+", t[pos:]).group(0)
            if word[:3] in _MONTHS and len(word) >= 3:
                out["M"] = float(_MONTHS.index(word[:3]) + 1)
                pos += len(word)
                continue
        # letters between numbers (2020m1, 2020q1, 10h30) are separators
        while pos < len(t) and t[pos].isalpha():
            pos += 1
        m = _DIGITS.match(t, pos)
        if m is None:
            return None
        run = m.group(0)
        nxt = fields[k + 1] if k + 1 < len(fields) else None
        width = 4 if f == "Y" else 2
        if nxt is not None and len(run) > width:
            run = run[:width]
        pos += len(run)
        val = float(run)
        if f == "s" and t[pos:pos + 1] == "." and t[pos + 1:pos + 2].isdigit():
            frac = _DIGITS.match(t, pos + 1).group(0)
            val += float("0." + frac)
            pos += 1 + len(frac)
        if f == "Y" and len(run) <= 2:
            if century:
                val += 100 * int(century)
            elif topyear is not None:
                val = topyear - ((topyear - val) % 100)
            else:
                return None
        if f != "#":
            out[f] = val
    rest = t[pos:].strip(" .,:/-")
    if rest in ("am", "pm", "a.m", "p.m") and "h" in out:
        out["h"] = out["h"] % 12 + (12 if rest.startswith("p") else 0)
    elif rest:
        return None
    return out


def daily(text: str, mask: str, topyear: float | None = None) -> float | None:
    f = parse(text, mask, topyear)
    if f is None or not {"D", "M", "Y"} <= f.keys():
        return None
    return mdy(f["M"], f["D"], f["Y"])


def clock(text: str, mask: str, topyear: float | None = None) -> float | None:
    """%tc time; the mask needs the date (DMY) and the parts of the time it lists."""
    f = parse(text, mask, topyear)
    if f is None:
        return None
    if {"D", "M", "Y"} <= f.keys():
        d = mdy(f["M"], f["D"], f["Y"])
        if d is None:
            return None
    elif any(c in f for c in "DMY"):
        return None
    else:
        d = 0.0
    h, mi, s = f.get("h", 0.0), f.get("m", 0.0), f.get("s", 0.0)
    if h > 23 or mi > 59 or s >= 61:
        return None
    return d * MS_DAY + hms(h, mi, s)


def period(text: str, mask: str, letter: str, per_year: int,
           topyear: float | None = None) -> float | None:
    """monthly() (M, 12), quarterly() (Q, 4), halfyearly() (H, 2), weekly() (W, 52)."""
    f = parse(text, mask, topyear)
    if f is None or "Y" not in f or letter not in f:
        return None
    k = f[letter]
    if not 1 <= k <= per_year:
        return None
    return per_year * (f["Y"] - 1960) + k - 1
