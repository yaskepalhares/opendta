"""estimates table ([R] estimates table): coeficientes de vários modelos
lado a lado, com erros-padrão, estatísticas e estrelas opcionais."""

from __future__ import annotations

import re
from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.formats import format_value
from ..lang.syntax import match_options

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS
COLW = 13


def _fmt(x: float, fmt: str) -> str:
    if x is None or not np.isfinite(x) or abs(x) >= SYS:
        return "."
    return format_value(float(x), fmt, pad=False).strip()


def _center(text: str, w: int = COLW) -> str:
    left = (w - len(text)) // 2
    return " " * left + text + " " * (w - len(text) - left)


def estimates_table(s: "Session", args: str) -> None:
    from .postest import _restore, _snapshot, _stored, _strip_marks
    from .results import current, pvalue
    head, comma, opts = args.partition(",")
    o = match_options(opts, {"b": 1, "se": 2, "t": 1, "p": 1, "stats": 2, "star": 4, "keep": 4,
                             "drop": 4, "title": 2, "varlabel": 4, "stfmt": 4, "newpanel": 3,
                             "style": 3, "varwidth": 4, "modelwidth": 6}) if comma and opts.strip() else {}
    names = head.split() or ["."]
    st_ = _stored(s)
    saved = _snapshot(s)
    models = []
    try:
        for n in names:
            if n in (".", "_active"):
                _restore(s, saved)
                est = current(s)
                label = "active" if names != ["."] else "active"
            else:
                if n not in st_:
                    raise StataError(111, f"estimation result {n} not found")
                _restore(s, st_[n])
                est = current(s)
                label = n
            e = dict(s.e)
            models.append((label, est, e))
    finally:
        _restore(s, saved)
    bfmt = _optfmt(o.get("b"), "%9.0g")
    sefmt = _optfmt(o.get("se"), bfmt)
    tfmt = _optfmt(o.get("t"), "%7.2f")
    pfmt = _optfmt(o.get("p"), "%7.4f")
    stfmt = _optfmt(o.get("stfmt"), "%9.0g")
    # linhas: união dos coeficientes, na ordem do primeiro aparecimento
    rows: list[str] = []
    for _, est, _ in models:
        for nm in est.names:
            if re.search(r"(^|#)\d+(b|o)\w*\.", nm) or nm.startswith("o."):
                continue
            k = _strip_marks(nm)
            if k not in rows:
                rows.append(k)
    keep = str(o.get("keep") or "").split()
    drop = str(o.get("drop") or "").split()
    if keep:
        rows = [r_ for r_ in rows if r_ in keep]
    if drop:
        rows = [r_ for r_ in rows if r_ not in drop]
    out = s.output
    width = 13 + 1 + COLW * len(models)
    if o.get("title"):
        out.write("\n" + str(o["title"]).strip('"') + "\n", "text")
    out.write("\n" + "-" * width + "\n", "text")
    out.write(f"{'Variable':>12} |" + "".join(_center(m[0][:12]) for m in models) + "\n", "text")
    out.write("-" * 13 + "+" + "-" * (COLW * len(models)) + "\n", "text")
    stars = bool(o.get("star"))
    for rname in rows:
        cells, ses, ts, ps = [], [], [], []
        for _, est, _ in models:
            idx = next((i for i, n in enumerate(est.names) if _strip_marks(n) == rname), None)
            if idx is None:
                cells.append("")
                ses.append("")
                ts.append("")
                ps.append("")
                continue
            b = est.b[idx]
            se = np.sqrt(est.V[idx, idx]) if est.V[idx, idx] > 0 else 0.0
            tv = b / se if se > 0 else SYS
            pv = pvalue(est.stat, tv, est.df_r) if se > 0 else SYS
            txt = _fmt(b, bfmt)
            if stars:
                mark = "***" if pv < .001 else "**" if pv < .01 else "*" if pv < .05 else ""
                txt += mark.ljust(3)
            cells.append(txt)
            ses.append(_fmt(se, sefmt))
            ts.append(_fmt(tv, tfmt))
            ps.append(_fmt(pv, pfmt))
        lab = rname if len(rname) <= 12 else rname[:11] + "~"
        out.write(f"{lab:>12} |", "text")
        out.write("".join(_cell(c, stars) for c in cells) + "\n", "result")
        for flag, vals in (("se", ses), ("t", ts), ("p", ps)):
            if o.get(flag) is not None and o.get(flag) is not False:
                out.write(f"{'':>12} |", "text")
                out.write("".join(_cell(c, stars) for c in vals) + "\n", "result")
    stats = str(o.get("stats") or "").split()
    if stats:
        out.write("-" * 13 + "+" + "-" * (COLW * len(models)) + "\n", "text")
        for stn in stats:
            out.write(f"{stn:>12} |", "text")
            vals = []
            for _, _, e in models:
                v = e.get(stn)
                vals.append(_fmt(v, stfmt) if isinstance(v, float) else "")
            out.write("".join(_cell(c, stars) for c in vals) + "\n", "result")
    out.write("-" * width + "\n", "text")
    if stars:
        out.write(f"{'legend: * p<0.05; ** p<0.01; *** p<0.001':>{width}}\n", "text")


def _cell(text: str, stars: bool) -> str:
    if stars:
        return f" {text:>12}" if text else " " * COLW
    return f" {text:>10}  "


def _optfmt(v, default: str) -> str:
    if isinstance(v, str) and v.strip().startswith("%"):
        return v.strip()
    return default
