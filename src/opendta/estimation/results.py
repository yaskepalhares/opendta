"""Resultados de estimação: e(b), e(V), e(sample) e a tabela de coeficientes
([R] estimation commands, [P] ereturn, [R] estimates).

A tabela segue o layout do Stata 14 (78 colunas):

    ------------------------------------------------------------------------------
           price |      Coef.   Std. Err.      t    P>|t|     [95% Conf. Interval]
    -------------+----------------------------------------------------------------
             mpg |  -49.51222   86.15604    -0.57   0.567    -221.3025     122.278
    ...
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import TYPE_CHECKING

import numpy as np
from scipy import stats as st

from ..core import missing as M
from ..core.formats import format_value

if TYPE_CHECKING:
    from ..session import Session
    from .fvars import Column

SYS = M.SYSMISS
LINE = "-" * 78


def g9(x: float, fmt: str = "%9.0g") -> str:
    if x is None or not np.isfinite(x) or abs(x) >= SYS:
        return "."
    return format_value(float(x), fmt, pad=False).strip()


def abbrev(name: str, width: int = 12) -> str:
    """Abrevia como o Stata nas tabelas: começo~fim."""
    if len(name) <= width:
        return name
    return name[: width - 1] + "~" if "#" not in name else name[: width - 1] + "~"


@dataclass
class CoefRow:
    """Uma linha da tabela."""
    kind: str                     # "coef", "header", "blank", "base", "omitted", "empty"
    label: str = ""
    b: float = 0.0
    se: float = 0.0
    eq: str = ""
    name: str = ""                # nome em e(b), para estimates table


@dataclass
class Estimates:
    """O que uma estimação deixa para a pós-estimação."""
    cmd: str
    depvar: str
    names: list[str]                       # colnames de e(b), com b./o.
    b: np.ndarray
    V: np.ndarray
    N: int
    stat: str = "t"                        # "t" ou "z"
    df_r: float | None = None
    rows: list[CoefRow] = field(default_factory=list)
    eqnames: list[str] = field(default_factory=list)


def post(s: "Session", est: Estimates, scalars: list[tuple[str, float]],
         macros: list[tuple[str, str]], sample: np.ndarray, *, extra_mats=()) -> None:
    """Grava e(): o Stata lista na ordem inversa da gravação (ver progcmd)."""
    from ..commands.matrix import Matrix
    e: dict = {}
    for k, v in reversed(macros):
        if v != "":
            e[k] = v
    for k, v in reversed(scalars):
        e[k] = float(v) if v is not None else SYS
    k = len(est.names)
    eqs = est.eqnames or ["_"] * k
    mats = [("b", Matrix(est.b.reshape(1, -1), ["y1"], list(est.names))),
            ("V", Matrix(est.V, list(est.names), list(est.names)))] + list(extra_mats)
    for name, m in reversed(mats):
        m.roweq = eqs if name == "V" else ["_"]
        m.coleq = eqs
        e[name] = m
    # matrizes primeiro na inserção = listadas por último; escalares e macros
    # mantêm a ordem acima
    s.e = e
    s.data.esample = np.asarray(sample, dtype=bool).copy()
    s._estimates = est


def current(s: "Session") -> Estimates:
    from ..core.errors import StataError
    est = getattr(s, "_estimates", None)
    if est is None or "b" not in s.e:
        raise StataError(301, "last estimates not found")
    return est


# ---------------------------------------------------------------------------
# tabela
# ---------------------------------------------------------------------------

def pvalue(stat: str, x: float, df: float | None) -> float:
    if not np.isfinite(x):
        return SYS
    if stat == "t" and df is not None and df > 0:
        return float(2 * st.t.sf(abs(x), df))
    return float(2 * st.norm.sf(abs(x)))


def crit(stat: str, level: float, df: float | None) -> float:
    a = 1 - level / 100
    if stat == "t" and df is not None and df > 0:
        return float(st.t.ppf(1 - a / 2, df))
    return float(st.norm.ppf(1 - a / 2))


def level_text(level: float) -> str:
    return f"{level:g}"


def coef_table(s: "Session", depvar: str, rows: list[CoefRow], *, stat: str, df: float | None,
               level: float, vcetype: str = "", coef_title: str = "Coef.",
               footer: bool = True, eq_headers: bool = False, eform: bool = False) -> None:
    out = s.output
    lv = level_text(level)
    ci = f"[{lv}% Conf. Interval]"
    out.write(LINE + "\n", "text")
    if vcetype:
        # VERIFICAR alinhamento do rótulo (Robust, Bootstrap...)
        lab = vcetype
        start = 28 + (9 - len(lab)) // 2      # centrado sobre "Std. Err."
        out.write(f"{'':>12} |" + " " * max(start - 14, 1) + lab + "\n", "text")
    stat_h = "t    P>|t|" if stat == "t" else "z    P>|z|"
    out.write(f"{abbrev(depvar):>12} | {coef_title:>10}   Std. Err.      {stat_h}     {ci:>20}\n", "text")
    out.write("-" * 13 + "+" + "-" * 64 + "\n", "text")
    c = crit(stat, level, df)
    for r in rows:
        if r.kind == "eq":
            out.write(f"{abbrev(r.label):<12} |", "result")
            out.write(("  " + r.eq if r.eq else "") + "\n", "text")   # ex.: (base outcome)
            continue
        if r.kind == "header":
            out.write(f"{abbrev(r.label):>12} |\n", "text")
            continue
        if r.kind == "blank":
            out.write(f"{'':>12} |\n", "text")
            continue
        if r.kind == "sep":
            out.write("-" * 13 + "+" + "-" * 64 + "\n", "text")
            continue
        out.write(f"{abbrev(r.label):>12} |", "text")
        if r.kind == "base":
            out.write("          0  (base)\n", "result")
            continue
        if r.kind == "omitted":
            out.write("          0  (omitted)\n", "result")
            continue
        if r.kind == "empty":
            out.write("          0  (empty)\n", "result")
            continue
        b, se = r.b, r.se
        if r.kind == "aux":
            # parâmetros auxiliares (/cut1, /sigma, alpha): sem z e P>|z|; o
            # intervalo pode vir pronto (transformado) em r.eq = "lo hi"
            if r.eq:
                lo, hi = (float(v) for v in r.eq.split())
            else:
                lo, hi = b - c * se, b + c * se
            out.write(f"  {g9(b):>9}  {g9(se):>9}" + " " * 17 + f"    {g9(lo):>9}   {g9(hi):>9}\n",
                      "result")
            continue
        if se > 0 and np.isfinite(se) and abs(se) < SYS:
            tval = b / se
            p = pvalue(stat, tval, df)
            lo, hi = b - c * se, b + c * se
            if eform:
                b, se, lo, hi = np.exp(b), np.exp(b) * se, np.exp(lo), np.exp(hi)
            out.write(f"  {g9(b):>9}  {g9(se):>9}  {tval:>7.2f}  {p:>6.3f}    {g9(lo):>9}   {g9(hi):>9}\n",
                      "result")
        else:
            # VERIFICAR: erro-padrão zero ou missing (constrained)
            out.write(f"  {g9(b):>9}  {g9(se) if se == 0 else '.':>9}  {'.':>7}  {'.':>6}    "
                      f"{'.':>9}   {'.':>9}\n", "result")
    if footer:
        out.write(LINE + "\n", "text")


def rows_from_columns(cols: "list[Column]", b: dict[str, float], se: dict[str, float],
                      ds, *, constant: bool, show_base: bool = False) -> list[CoefRow]:
    """Linhas da tabela: termos fatoriais e variáveis com operadores de séries
    temporais ganham cabeçalho e linha em branco antes e depois."""
    rows: list[CoefRow] = []
    i = 0
    n = len(cols)

    def blank_before():
        if rows and rows[-1].kind != "blank":
            rows.append(CoefRow("blank"))

    while i < n:
        c = cols[i]
        j = i + 1
        if c.group:
            while j < n and cols[j].group == c.group:
                j += 1
        group = cols[i:j]
        shown = [g for g in group if not g.base or show_base]
        interaction = "#" in c.name
        grouped = bool(c.group) and (len(shown) > 1 or (c.group[0] == "t" and interaction)
                                     or (c.group[0] == "v" and c.rowlab != "--."))
        if not grouped:
            # uma linha só: "x", "1.h" (fator com um nível), "c.x#c.x"; as
            # interações ficam entre linhas em branco (compat 0502)
            if interaction:
                blank_before()
            for g in shown:
                lab = g.labels[0] if (g.labels and not g.group) or (g.group and g.group[0] == "v") \
                    else _strip_b(g.name)
                if g.base:
                    rows.append(CoefRow("base", lab, name=g.name))
                elif g.omitted:
                    rows.append(CoefRow("omitted", lab, name=g.name))
                else:
                    rows.append(CoefRow("coef", lab, b[g.name], se[g.name], name=g.name))
            if interaction:
                rows.append(CoefRow("blank"))
            i = j
            continue
        blank_before()
        rows.append(CoefRow("header", c.labels[0] if c.group[0] == "t" else c.group[1]))
        for g in group:
            lab = (_cell_label(ds, g) + " ") if g.group[0] == "t" else g.rowlab
            if g.base:
                if show_base:
                    rows.append(CoefRow("base", lab, name=g.name))
                continue
            if g.omitted:
                rows.append(CoefRow("omitted", lab, name=g.name))
                continue
            rows.append(CoefRow("coef", lab, b[g.name], se[g.name], name=g.name))
        rows.append(CoefRow("blank"))
        i = j
    if constant:
        rows.append(CoefRow("coef", "_cons", b["_cons"], se["_cons"], name="_cons"))
    elif rows and rows[-1].kind == "blank":
        rows.pop()
    return rows


def _head_label(head: str) -> str:
    return head


def _strip_b(name: str) -> str:
    import re
    return "#".join(re.sub(r"^(\d+)(?:b|bn|o)+\.", r"\1.", p) for p in name.split("#"))


def _cell_label(ds, col) -> str:
    """Nível de um fator: o rótulo de valor, se houver (VERIFICAR largura)."""
    from .fvars import _level_text
    parts = col.name.split("#")
    texts = []
    k = 0
    for p in parts:
        if p.startswith("c.") or "." not in p or not p[0].isdigit():
            continue
        var = p.split(".", 1)[1]
        lv = col.cell[k] if k < len(col.cell) else None
        k += 1
        txt = _level_text(lv) if lv is not None else ""
        try:
            v = ds.get(var)
            if v.value_label and lv is not None:
                lab = ds.value_labels.get(v.value_label, {}).get(int(lv))
                if lab:
                    txt = lab
        except Exception:   # noqa: BLE001
            pass
        texts.append(txt)
    return "#".join(texts) if texts else col.labels[1]
