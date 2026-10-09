"""Séries temporais e painéis: tsset, xtset e os operadores L. F. D. S.
([TS] tsset, [XT] xtset, [U] 11.4.4 time-series varlists).

A declaração fica nas características de _dta, como no Stata (_TStvar,
_TSpanel, _TSdelta), e por isso acompanha o arquivo salvo.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError

if TYPE_CHECKING:
    from ..core.dataset import Dataset
    from ..session import Session

SYS = M.SYSMISS


@dataclass
class TSInfo:
    tvar: str
    panel: str = ""
    delta: float = 1.0


def ts_info(ds: "Dataset") -> TSInfo | None:
    ch = ds.chars.get("_dta", {})
    tvar = ch.get("_TStvar", "")
    if not tvar or not ds.has(tvar):
        return None
    panel = ch.get("_TSpanel", "")
    if panel and not ds.has(panel):
        panel = ""
    try:
        delta = float(ch.get("_TSdelta", "1") or 1)
    except ValueError:
        delta = 1.0
    return TSInfo(tvar, panel, delta)


def panel_info(ds: "Dataset") -> str:
    """Variável de painel declarada por xtset (com ou sem tempo)."""
    ch = ds.chars.get("_dta", {})
    p = ch.get("_TSpanel", "") or ch.get("iis", "")
    return p if p and ds.has(p) else ""


def require_ts(ds: "Dataset") -> TSInfo:
    info = ts_info(ds)
    if info is None:
        # VERIFICAR texto exato
        raise StataError(111, "time variable not set, use tsset varname ...")
    return info


def shift(ds: "Dataset", values: np.ndarray, k: int) -> np.ndarray:
    """Valor de `values` k períodos antes (k > 0, lag) ou depois (k < 0, lead),
    no mesmo painel; missing quando o período não existe."""
    info = require_ts(ds)
    t = ds.get(info.tvar).data
    n = ds.nobs
    if info.panel:
        pv = ds.get(info.panel).data
    else:
        pv = np.zeros(n)
    index: dict[tuple[float, float], int] = {}
    for i in range(n):
        if t[i] < SYS and pv[i] < SYS:
            index[(pv[i], t[i])] = i
    out = np.full(n, SYS)
    step = k * info.delta
    for i in range(n):
        if t[i] >= SYS or pv[i] >= SYS:
            continue
        j = index.get((pv[i], t[i] - step))
        if j is not None:
            out[i] = values[j]
    return out


def apply_ops(ds: "Dataset", ops: list[tuple[str, int]], values: np.ndarray) -> np.ndarray:
    """Aplica uma sequência de operadores (letra, ordem) a um vetor."""
    x = np.asarray(values, dtype=np.float64)
    for op, k in ops:
        for _ in range(k if op in "DS" else 1):
            if op == "L":
                y = shift(ds, x, k)
            elif op == "F":
                y = shift(ds, x, -k)
            elif op == "D":
                prev = shift(ds, x, 1)
                y = np.where((x < SYS) & (prev < SYS), x - prev, SYS)
            elif op == "S":
                # S.x sazonal: x - L.x; Sk.x = x - Lk.x (aplicado uma vez)
                prev = shift(ds, x, k)
                y = np.where((x < SYS) & (prev < SYS), x - prev, SYS)
                x = y
                break
            else:
                raise StataError(198, f"{op}: unknown operator")
            x = y
    return x


def ops_name(ops: list[tuple[str, int]]) -> str:
    """Prefixo canônico: [('L', 2), ('D', 1)] -> 'L2D'. VERIFICAR a ordem."""
    out = []
    for op, k in ops:
        out.append(op + (str(k) if k != 1 else ""))
    return "".join(out)


# ---------------------------------------------------------------------------
# tsset / xtset
# ---------------------------------------------------------------------------

_UNITS = {"d": "day", "w": "week", "m": "month", "q": "quarter", "h": "halfyear", "y": "year",
          "c": "ms", "C": "ms"}


def _fmt_time(ds: "Dataset", tvar: str, x: float) -> str:
    from ..core.formats import format_value
    v = ds.get(tvar)
    if v.fmt.startswith("%t"):
        return format_value(x, v.fmt, pad=False).strip()
    return format_value(x, "%9.0g", pad=False).strip()


def _declare(s: "Session", args: str, *, xt: bool) -> None:
    from ..commands._util import plural  # noqa: F401
    from ..core.varlist import resolve_name
    from ..lang.syntax import match_options
    ds = s.data
    head, comma, opts = args.partition(",")
    words = head.split()
    o = match_options(opts, {"delta": 1, "clear": 5, "daily": 5, "weekly": 6, "monthly": 7,
                             "quarterly": 9, "yearly": 6, "generic": 7, "format": 3}) \
        if comma and opts.strip() else {}
    ch = ds.chars.setdefault("_dta", {})
    if o.get("clear"):
        for k in ("_TStvar", "_TSpanel", "_TSdelta", "_TSitrvl", "tis", "iis"):
            ch.pop(k, None)
        return
    if not words:
        info = ts_info(ds)
        if xt:
            p = panel_info(ds)
            if not p:
                raise StataError(459, "panel variable not set; use xtset varname ...")   # VERIFICAR
            _report(s, p, info.tvar if info else "", info.delta if info else 1.0, xt=True)
        else:
            if info is None:
                raise StataError(111, "time variable not set, use tsset varname ...")   # VERIFICAR
            _report(s, info.panel, info.tvar, info.delta, xt=False)
        return
    names = [resolve_name(ds, w) for w in words]
    if len(names) > 2:
        raise StataError(103, "too many variables specified")
    if xt:
        panel = names[0]
        tvar = names[1] if len(names) > 1 else ""
    else:
        panel, tvar = (names[0], names[1]) if len(names) == 2 else ("", names[0])
    for nm in [panel, tvar]:
        if nm and ds.get(nm).is_string:
            raise StataError(109, f"variable {nm} must be numeric")   # VERIFICAR
    delta = 1.0
    if o.get("delta"):
        delta = float(s.eval(str(o["delta"]).split()[0]))
    if tvar:
        t = ds.get(tvar).data
        if np.any((t < SYS) & (t != np.trunc(t))):
            raise StataError(451, "time variable must contain only integer values")   # VERIFICAR
        pv = ds.get(panel).data if panel else np.zeros(ds.nobs)
        seen: set = set()
        for i in range(ds.nobs):
            if t[i] >= SYS or pv[i] >= SYS:
                continue
            key = (pv[i], t[i])
            if key in seen:
                raise StataError(451, "repeated time values within panel" if panel
                                 else "repeated time values in sample")
            seen.add(key)
    ch.pop("_TStvar", None)
    ch.pop("_TSpanel", None)
    ch.pop("_TSdelta", None)
    if panel:
        ch["_TSpanel"] = panel
    if tvar:
        ch["_TStvar"] = tvar
        ch["_TSdelta"] = f"{delta:g}"
    # tsset/xtset ordena os dados pelo painel e pelo tempo
    from ..core.sorting import sort
    keys = [k for k in (panel, tvar) if k]
    sort(ds, keys)
    _report(s, panel, tvar, delta, xt=xt)
    s.notify_state()


def _report(s: "Session", panel: str, tvar: str, delta: float, *, xt: bool) -> None:
    ds = s.data
    out = s.output
    s.r = {}
    pv = ds.get(panel).data if panel else np.zeros(ds.nobs)
    t = ds.get(tvar).data if tvar else None
    if panel:
        balance = "strongly balanced"
        if t is not None:
            ok = (pv < SYS) & (t < SYS)
            sets: dict[float, set] = {}
            for p_, tt in zip(pv[ok], t[ok]):
                sets.setdefault(p_, set()).add(tt)
            vals = list(sets.values())
            if vals and any(v != vals[0] for v in vals):
                lens = {len(v) for v in vals}
                balance = "weakly balanced" if len(lens) == 1 else "unbalanced"
        out.write(f"{'panel variable:':>21}  ", "text")
        out.write(f"{panel} ({balance})\n", "result")
    if t is not None:
        ok = t < SYS
        if panel:
            ok &= pv < SYS
        tmin, tmax = (float(t[ok].min()), float(t[ok].max())) if ok.any() else (SYS, SYS)
        gaps = False
        if ok.any():
            groups: dict[float, list[float]] = {}
            for p_, tt in zip(pv[ok], t[ok]):
                groups.setdefault(p_, []).append(tt)
            for ts in groups.values():
                ts = sorted(ts)
                if any(abs((b - a) - delta) > 1e-9 for a, b in zip(ts, ts[1:])):
                    gaps = True
                    break
        line = f"{tvar}, {_fmt_time(ds, tvar, tmin)} to {_fmt_time(ds, tvar, tmax)}"
        if gaps:
            line += ", but with gaps"   # VERIFICAR
        out.write(f"{'time variable:':>21}  ", "text")
        out.write(line + "\n", "result")
        fmt = ds.get(tvar).fmt
        unit = "unit" if not fmt.startswith("%t") else _UNITS.get(fmt[2:3], "unit")
        dtext = f"{delta:g} {unit}" + ("s" if delta != 1 else "")
        out.write(f"{'delta:':>21}  ", "text")
        out.write(dtext + "\n", "result")
        s.r.update({"tmax": tmax, "tmin": tmin, "tdelta": float(delta), "gaps": 1.0 if gaps else 0.0})
        s.r["timevar"] = tvar
        s.r["unit"] = unit if unit != "unit" else "generic"
        s.r["tsfmt"] = ds.get(tvar).fmt
    if panel:
        s.r["panelvar"] = panel


def cmd_tsset(s: "Session", args: str) -> None:
    _declare(s, args, xt=False)


def cmd_xtset(s: "Session", args: str) -> None:
    _declare(s, args, xt=True)
