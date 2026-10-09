"""Exibição de resultados do Mata (o que aparece ao digitar uma expressão).

Escalar: dois espaços e o valor. Matriz: números das colunas centrados
sobre cada coluna, caixa com +---+ e o número de cada linha:

           1   2
        +---------+
      1 |  1   2  |
      2 |  3   4  |
        +---------+

Cada coluna tem a largura do seu maior elemento; os números usam %9.0g.
Matrizes simétricas aparecem como [symmetric], só com o triângulo de
baixo. Detalhes marcados VERIFICAR têm casos em compat/do/04*.do.
"""

from __future__ import annotations

import numpy as np

from ..core import missing as M
from ..core.formats import format_value
from .values import MV, Pointer

REAL_FORMAT = "%9.0g"      # VERIFICAR: formato dos números na exibição do Mata


def cell(x, t: str) -> str:
    if t == "real":
        x = float(x)
        if x >= M.SYSMISS:
            return M.missing_name(x)
        return format_value(x, REAL_FORMAT, pad=False).strip()
    if t == "string":
        return str(x)
    if t == "complex":
        z = complex(x)
        re_ = format_value(z.real, REAL_FORMAT, pad=False).strip()
        im = format_value(abs(z.imag), REAL_FORMAT, pad=False).strip()
        return f"{re_} {'-' if z.imag < 0 else '+'} {im}i"
    if t == "pointer":
        if x is None:
            return "NULL"
        return f"0x{id(x) & 0xffffffff:x}"
    if t == "struct":
        return f"0x{id(x) & 0xffffffff:x}"
    return str(x)


def _symmetric(v: MV) -> bool:
    if v.t != "real" or v.rows != v.cols or v.rows < 2:
        return False
    return bool(np.array_equal(v.a, v.a.T))


def render(v: MV) -> list[str]:
    r, c = v.a.shape
    if r == 0 or c == 0:
        return []
    if r == 1 and c == 1:
        return ["  " + cell(v.a[0, 0], v.t)]
    sym = _symmetric(v)
    texts = [[cell(v.a[i, j], v.t) if (not sym or j <= i) else "" for j in range(c)] for i in range(r)]
    widths = [max(len(texts[i][j]) for i in range(r)) for j in range(c)]
    widths = [max(w, 1) for w in widths]
    lab_w = max(3, len(str(r)))
    lead = lab_w + 4                 # "  1 |" + 2 espaços antes da 1ª coluna
    lines: list[str] = []
    if sym:
        lines.append("[symmetric]")
    head = [" "] * (lead + sum(widths) + 3 * (c - 1) + 2)
    pos = lead
    for j, w in enumerate(widths):
        lab = str(j + 1)
        start = pos + max(0, (w - len(lab)) // 2)
        if w < len(lab):
            start = pos + w - len(lab)
        for k, ch in enumerate(lab):
            if 0 <= start + k < len(head):
                head[start + k] = ch
        pos += w + 3
    lines.append("".join(head).rstrip())
    inner = 2 + sum(widths) + 3 * (c - 1) + 2
    bar = " " * (lab_w + 1) + "+" + "-" * inner + "+"
    lines.append(bar)
    for i in range(r):
        # VERIFICAR: strings alinhadas à direita como os números
        body = "   ".join(texts[i][j].rjust(widths[j]) for j in range(c))
        lines.append(f"{i + 1:>{lab_w}} |  {body}  |")
    lines.append(bar)
    return lines


def is_displayable(v) -> bool:
    return isinstance(v, MV) and not isinstance(v.a, Pointer)
