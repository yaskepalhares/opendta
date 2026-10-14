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
        raise StataError(111, "time variable not set")
    return info


class _TSIndex:
    """Índice (painel, tempo) -> observação, vetorizado: as chaves viram
    inteiros painel*(amplitude+1) + (t - tmin) e a busca é binária."""

    def __init__(self, ds: "Dataset", info: TSInfo):
        n = ds.nobs
        t = ds.get(info.tvar).data.astype(np.float64)
        pv = ds.get(info.panel).data.astype(np.float64) if info.panel else np.zeros(n)
        self.ok = (t < SYS) & (pv < SYS)
        self.t, self.pv = t, pv
        self.n = n
        tt = t[self.ok]
        if not len(tt):
            self.keys = np.zeros(0, dtype=np.int64)
            self.order = np.zeros(0, dtype=np.int64)
            return
        self.tmin = float(tt.min())
        span = float(tt.max()) - self.tmin
        codes = np.zeros(n, dtype=np.int64)
        codes[self.ok] = np.unique(pv[self.ok], return_inverse=True)[1]
        self.codes = codes
        npan = int(codes.max()) + 1
        # margem para t - passo sair do intervalo sem colidir com outro painel
        self.width = int(3 * span + 3)
        self.vector = (npan + 1) * self.width < 2 ** 62
        if self.vector:
            key = codes * self.width + (t - self.tmin + span + 1)
            key = np.where(self.ok, key, -1).astype(np.int64)
            idx = np.nonzero(self.ok)[0]
            self.order = idx[np.argsort(key[idx], kind="stable")]
            self.keys = key[self.order]
            self.span = span
        else:
            self.lookup = {(pv[i], t[i]): i for i in np.nonzero(self.ok)[0]}

    def source(self, step: float) -> np.ndarray:
        """Para cada observação, a linha com tempo t - step no mesmo painel (-1 se não há)."""
        out = np.full(self.n, -1, dtype=np.int64)
        if not self.ok.any():
            return out
        if not self.vector:
            for i in np.nonzero(self.ok)[0]:
                out[i] = self.lookup.get((self.pv[i], self.t[i] - step), -1)
            return out
        tq = self.t - step
        inside = self.ok & (np.abs(tq - self.tmin) <= 2 * self.span + 1)
        q = (self.codes * self.width + (tq - self.tmin + self.span + 1)).astype(np.int64)
        pos = np.searchsorted(self.keys, q)
        pos = np.clip(pos, 0, max(len(self.keys) - 1, 0))
        hit = inside & (self.keys[pos] == q) & (np.round(tq) == tq)
        out[hit] = self.order[pos[hit]]
        return out


def shift(ds: "Dataset", values: np.ndarray, k: int, index: _TSIndex | None = None) -> np.ndarray:
    """Valor de `values` k períodos antes (k > 0, lag) ou depois (k < 0, lead),
    no mesmo painel; missing quando o período não existe."""
    info = require_ts(ds)
    index = index or _TSIndex(ds, info)
    src = index.source(k * info.delta)
    vals = np.asarray(values, dtype=np.float64)
    return np.where(src >= 0, vals[np.clip(src, 0, None)], SYS) if len(vals) else vals.copy()


def apply_ops(ds: "Dataset", ops: list[tuple[str, int]], values: np.ndarray) -> np.ndarray:
    """Aplica os operadores (letra, ordem) a um vetor. L e F se somam num
    deslocamento líquido (LF.y = y, como no Stata, compat 0503); D e S são
    aplicados antes do deslocamento (os operadores comutam)."""
    x = np.asarray(values, dtype=np.float64)
    index = _TSIndex(ds, require_ts(ds))
    net = sum(k for op, k in ops if op == "L") - sum(k for op, k in ops if op == "F")
    rest = [(op, k) for op, k in ops if op not in "LF"]
    if net:
        rest.append(("L", net) if net > 0 else ("F", -net))
    for op, k in rest:
        for _ in range(k if op in "DS" else 1):
            if op == "L":
                y = shift(ds, x, k, index)
            elif op == "F":
                y = shift(ds, x, -k, index)
            elif op == "D":
                prev = shift(ds, x, 1, index)
                y = np.where((x < SYS) & (prev < SYS), x - prev, SYS)
            elif op == "S":
                # S.x sazonal: x - L.x; Sk.x = x - Lk.x (aplicado uma vez)
                prev = shift(ds, x, k, index)
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
                raise StataError(111, "time variable not set")
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
        out.write(f"{'panel variable:':>22}  ", "text")
        out.write(f"{panel} ({balance})\n", "result")
    if t is not None:
        ok = t < SYS
        if panel:
            ok &= pv < SYS
        tmin, tmax = (float(t[ok].min()), float(t[ok].max())) if ok.any() else (SYS, SYS)
        ngaps = 0
        if ok.any():
            groups: dict[float, list[float]] = {}
            for p_, tt in zip(pv[ok], t[ok]):
                groups.setdefault(p_, []).append(tt)
            for ts in groups.values():
                ts = sorted(ts)
                ngaps += sum(1 for a, b in zip(ts, ts[1:]) if abs((b - a) - delta) > 1e-9)
        gaps = ngaps > 0
        line = f"{tvar}, {_fmt_time(ds, tvar, tmin)} to {_fmt_time(ds, tvar, tmax)}"
        if gaps:
            line += ", but with a gap" if ngaps == 1 else ", but with gaps"
        out.write(f"{'time variable:':>22}  ", "text")
        out.write(line + "\n", "result")
        fmt = ds.get(tvar).fmt
        unit = "unit" if not fmt.startswith("%t") else _UNITS.get(fmt[2:3], "unit")
        dtext = f"{delta:g} {unit}" + ("s" if delta != 1 else "")
        out.write(f"{'delta:':>22}  ", "text")
        out.write(dtext + "\n", "result")
        tfmt = ds.get(tvar).fmt
        code = tfmt[2:3] if tfmt.startswith("%t") else ""
        tmins, tmaxs = _fmt_time(ds, tvar, tmin), _fmt_time(ds, tvar, tmax)
        res = {}
        if panel:
            res["balanced"] = balance
        res["tmins"], res["tmaxs"], res["tdeltas"] = tmins, tmaxs, dtext
        res["tsfmt"] = tfmt
        res["unit1"] = code or "."     # VERIFICAR para formatos %t
        if code:
            res["unit"] = {"d": "daily", "w": "weekly", "m": "monthly", "q": "quarterly",
                           "h": "halfyearly", "y": "yearly", "c": "clocktime", "C": "clocktime"}[code]
        res["timevar"] = tvar
        if panel:
            res["panelvar"] = panel
            pvals = pv[pv < SYS]
            res["imin"] = float(pvals.min()) if len(pvals) else SYS
            res["imax"] = float(pvals.max()) if len(pvals) else SYS
        res["tmin"], res["tmax"], res["tdelta"] = tmin, tmax, float(delta)
        s.r = res          # o Stata lista na ordem inversa da gravação
        return
    if panel:
        s.r["panelvar"] = panel


def cmd_tsset(s: "Session", args: str) -> None:
    _declare(s, args, xt=False)


def cmd_xtset(s: "Session", args: str) -> None:
    _declare(s, args, xt=True)
