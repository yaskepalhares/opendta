"""estimates table ([R] estimates table): coeficientes de vários modelos
lado a lado, com erros-padrão, estatísticas e estrelas opcionais.

Layout observado no Stata 14 (compat 0504): cada modelo ocupa a largura do
formato + 3 colunas (+3 com estrelas); o nome do modelo vem centrado à
esquerda; os fatores aparecem agrupados como na tabela de coeficientes.
"""

from __future__ import annotations

from typing import TYPE_CHECKING

import numpy as np

from ..core import missing as M
from ..core.errors import StataError
from ..core.formats import format_value, parse_format
from ..lang.syntax import match_options

if TYPE_CHECKING:
    from ..session import Session

SYS = M.SYSMISS


def _fmt(x, fmt: str) -> str:
    if x is None or not isinstance(x, float) or not np.isfinite(x) or abs(x) >= SYS:
        return "."
    return format_value(float(x), fmt, pad=False).strip()


def _optfmt(v, default: str) -> str:
    if isinstance(v, str) and v.strip().startswith("%"):
        return v.strip()
    return default


def _merge_rows(models) -> list:
    """União das linhas (CoefRow) dos modelos, na ordem do primeiro em que
    cada coeficiente aparece; cabeçalhos e brancos vêm junto."""
    rows = []
    seen = set()
    for _, est, _ in models:
        pending = []
        insert_at = None
        for r in est.rows:
            if r.kind in ("coef", "omitted", "base"):
                if r.name in seen:
                    # já existe: os novos que vieram antes entram antes dele
                    if pending:
                        pos = next(i for i, x in enumerate(rows) if x.name == r.name)
                        start = pos
                        while start > 0 and rows[start - 1].kind in ("header", "blank") and \
                                any(p.kind == "header" for p in pending):
                            start -= 1
                        rows[start:start] = pending
                        pending = []
                    continue
                seen.add(r.name)
                pending.append(r)
            elif r.kind in ("header", "blank"):
                pending.append(r)
        if pending:
            # novos no fim: antes de _cons, se houver
            pos = next((i for i, x in enumerate(rows) if x.name == "_cons"), len(rows))
            rows[pos:pos] = pending
    # limpa brancos repetidos e cabeçalhos sem linhas
    out = []
    for r in rows:
        if r.kind == "blank" and (not out or out[-1].kind == "blank"):
            continue
        out.append(r)
    if out and out[-1].kind == "blank":
        out.pop()
    return out


def estimates_table(s: "Session", args: str) -> None:
    from .postest import _restore, _snapshot, _stored
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
                label = "active"
            else:
                if n not in st_:
                    raise StataError(111, f"estimation result {n} not found")
                _restore(s, st_[n])
                label = n
            est = current(s)
            models.append((label, est, dict(s.e)))
    finally:
        _restore(s, saved)
    bfmt = _optfmt(o.get("b"), "%10.0g")     # VERIFICAR: padrão observado com 8 algarismos
    sefmt = _optfmt(o.get("se"), bfmt)
    tfmt = _optfmt(o.get("t"), bfmt)
    pfmt = _optfmt(o.get("p"), bfmt)
    stfmt = _optfmt(o.get("stfmt"), bfmt)
    fw = parse_format(bfmt).width
    stars = bool(o.get("star"))
    W = fw + 3 + (3 if stars else 0)
    rows = _merge_rows(models)
    keep = str(o.get("keep") or "").split()
    drop = str(o.get("drop") or "").split()
    if keep or drop:
        rows = [r for r in rows if r.kind not in ("coef", "omitted", "base")
                or ((not keep or _shortname(r.name) in keep) and _shortname(r.name) not in drop)]
    out = s.output
    width = 13 + 1 + W * len(models)
    if o.get("title"):
        out.write("\n" + str(o["title"]).strip('"') + "\n", "text")
    out.write("\n" + "-" * width + "\n", "text")
    hdr = ""
    for m in models:
        nm = m[0][: W - 1]
        left = (W - 1 - len(nm)) // 2
        hdr += " " * left + nm + " " * (W - left - len(nm))
    out.write(f"{'Variable':>12} |" + hdr + "\n", "text")
    out.write("-" * 13 + "+" + "-" * (W * len(models)) + "\n", "text")

    def cell(text: str, star: str = "") -> str:
        if not text:
            return " " * W
        return f"{text:>{fw + 1}}" + (star.ljust(3) if stars else "") + "  "

    show = [k for k in ("se", "t", "p") if o.get(k) is not None and o.get(k) is not False]
    for r in rows:
        if r.kind == "blank":
            out.write(f"{'':>12} |\n", "text")
            continue
        if r.kind == "header":
            out.write(f"{r.label:>12} |\n", "text")
            continue
        if r.kind != "coef":
            continue
        line, extra = "", {k: "" for k in show}
        for _, est, _ in models:
            idx = next((i for i, nn in enumerate(est.names) if nn == r.name), None)
            if idx is None:
                line += cell("")
                for k in show:
                    extra[k] += cell("")
                continue
            b = float(est.b[idx])
            se = float(np.sqrt(est.V[idx, idx])) if est.V[idx, idx] > 0 else 0.0
            tv = b / se if se > 0 else SYS
            pv = pvalue(est.stat, tv, est.df_r) if se > 0 else SYS
            mark = ""
            if stars and pv < SYS:
                mark = "***" if pv < .001 else "**" if pv < .01 else "*" if pv < .05 else ""
            line += cell(_fmt(b, bfmt), mark)
            vals = {"se": _fmt(se, sefmt), "t": _fmt(tv, tfmt), "p": _fmt(pv, pfmt)}
            for k in show:
                extra[k] += cell(vals[k])
        lab = r.label if len(r.label) <= 12 else r.label[:11] + "~"
        out.write(f"{lab:>12} |", "text")
        out.write(line + "\n", "result")
        for k in show:
            out.write(f"{'':>12} |", "text")
            out.write(extra[k] + "\n", "result")
    stats = str(o.get("stats") or "").split()
    if stats:
        out.write("-" * 13 + "+" + "-" * (W * len(models)) + "\n", "text")
        for stn in stats:
            out.write(f"{stn:>12} |", "text")
            vals = "".join(cell(_stat_text(m[2].get(stn), stfmt)) for m in models)
            out.write(vals + "\n", "result")
    out.write("-" * width + "\n", "text")
    legend = []
    if show:
        legend.append("legend: " + "/".join(["b"] + show))
    if stars:
        legend.append("legend: * p<0.05; ** p<0.01; *** p<0.001")
    for lg in legend:
        out.write(f"{lg:>{width}}\n", "text")


def _stat_text(v, fmt: str) -> str:
    """Estatísticas inteiras (N) saem sem casas decimais (compat 0504)."""
    if isinstance(v, float) and abs(v) < SYS and v == int(v):
        return f"{int(v)}"
    return _fmt(v, fmt)


def _shortname(name: str) -> str:
    from .postest import _strip_marks
    return _strip_marks(name)
